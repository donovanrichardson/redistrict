"""
Redistricting CLI.

Usage:
    redistrict                      (interactive prompts)
    redistrict --continue <run_id>  (re-run with curated water links)
    redistrict --help
"""

import sys

from InquirerPy import inquirer
from InquirerPy.base.control import Choice
from tqdm import tqdm

from redistrict import db, graph, partition


_FIPS_TO_NAME = {
    "01": "Alabama",       "02": "Alaska",        "04": "Arizona",
    "05": "Arkansas",      "06": "California",    "08": "Colorado",
    "09": "Connecticut",   "10": "Delaware",      "11": "District of Columbia",
    "12": "Florida",       "13": "Georgia",       "15": "Hawaii",
    "16": "Idaho",         "17": "Illinois",      "18": "Indiana",
    "19": "Iowa",          "20": "Kansas",        "21": "Kentucky",
    "22": "Louisiana",     "23": "Maine",         "24": "Maryland",
    "25": "Massachusetts", "26": "Michigan",      "27": "Minnesota",
    "28": "Mississippi",   "29": "Missouri",      "30": "Montana",
    "31": "Nebraska",      "32": "Nevada",        "33": "New Hampshire",
    "34": "New Jersey",    "35": "New Mexico",    "36": "New York",
    "37": "North Carolina","38": "North Dakota",  "39": "Ohio",
    "40": "Oklahoma",      "41": "Oregon",        "42": "Pennsylvania",
    "44": "Rhode Island",  "45": "South Carolina","46": "South Dakota",
    "47": "Tennessee",     "48": "Texas",         "49": "Utah",
    "50": "Vermont",       "51": "Virginia",      "53": "Washington",
    "54": "West Virginia", "55": "Wisconsin",     "56": "Wyoming",
    "72": "Puerto Rico",   "78": "U.S. Virgin Islands",
}


_GEOGRAPHY_LABELS = {
    "tracts":       "Census tracts",
    "block_groups": "Census block groups",
    "counties":     "Counties",
}


def _state_choices(conn, geography: str) -> list[Choice]:
    available = db.fetch_available_states(conn, geography)
    choices = []
    for fips in available:
        name = _FIPS_TO_NAME.get(fips, fips)
        choices.append(Choice(value=fips, name=f"{name} ({fips})"))
    return choices


def _ensure_adjacency(conn, geography: str, nodes: list[dict]) -> set[tuple[str, str]]:
    geoids = [n["geoid"] for n in nodes]
    missing = db.get_missing_adjacency_geoids(conn, geography, geoids)
    label = _GEOGRAPHY_LABELS.get(geography, geography)

    if missing:
        print(f"\nCalculating adjacency for {len(missing)} {label} "
              f"(first time for this state)...")
        bar = tqdm(total=len(missing), unit=geography, desc="Adjacency")

        def _cb(done, total, geoid):
            bar.update(1)

        inserted = db.compute_and_store_adjacency(
            conn, geography, missing, progress_callback=_cb
        )
        bar.close()
        print(f"  Inserted {inserted:,} new adjacency pairs.")
    else:
        print(f"Adjacency data already complete for all {label} in this state.")

    return db.fetch_adjacency(conn, geography, geoids)


def run(
    statefp: str,
    geography: str,
    n_districts: int,
    ncuts: int = 10,
    niter: int = 20,
    formula: str = "original",
    recursive: bool = False,
    water_penalty: float = graph.WATER_PENALTY,
) -> int:
    """
    Execute one redistricting run and return its run_id.

    Parameters
    ----------
    statefp:       2-digit FIPS code (e.g. '44' for Rhode Island).
    geography:     'tracts', 'block_groups', or 'counties'.
    n_districts:   Number of districts to produce.
    ncuts:         METIS independent attempts; best result kept (default 10).
    niter:         METIS refinement iterations per attempt (default 20).
    formula:       Edge weight formula: "original", "uniform",
                   "original_clamped", or "blend".
    recursive:     Use recursive bisection instead of k-way (disables contig).
    water_penalty: Divisor applied to non-rook-contiguous edge weights.
    """
    conn = db.connect()
    try:
        db.ensure_tables(conn, geography)

        print(f"\nFetching {geography} for state {statefp}...")
        nodes = db.fetch_nodes(conn, geography, statefp)
        if not nodes:
            print(f"No {geography} found for state {statefp}. Is the data loaded?")
            return
        print(f"  {len(nodes):,} nodes loaded.")

        active_nodes = [n for n in nodes if n["pop"] > 0]
        zero_pop_nodes = [n for n in nodes if n["pop"] == 0]
        if zero_pop_nodes:
            print(f"  {len(zero_pop_nodes):,} zero-pop nodes excluded from graph "
                  f"(will be assigned post-METIS).")

        if n_districts >= len(active_nodes):
            print(f"Error: n_districts ({n_districts}) must be less than "
                  f"active node count ({len(active_nodes)}).")
            return

        adjacent_pairs = _ensure_adjacency(conn, geography, nodes)
        print(f"  {len(adjacent_pairs):,} adjacency pairs loaded.")

        print("\nBuilding spherical Delaunay triangulation...")
        triangles = graph.spherical_delaunay_triangles(active_nodes)
        print(f"  {len(triangles):,} triangles.")

        print("Deriving Urquhart graph...")
        edges = graph.urquhart_edges(active_nodes, triangles)
        print(f"  {len(edges):,} edges.")

        n_water = sum(
            1 for i, j in edges
            if (min(active_nodes[i]["geoid"], active_nodes[j]["geoid"]),
                max(active_nodes[i]["geoid"], active_nodes[j]["geoid"])) not in adjacent_pairs
        )
        print(f"  {n_water:,} non-adjacent (water) edges.")

        params = {
            "status":           "pending",
            "water_penalty":    water_penalty,
            "edge_weight_scale": graph.EDGE_WEIGHT_SCALE,
            "formula":          formula,
            "recursive":        recursive,
            "ncuts":            ncuts,
            "niter":            niter,
            "n_nodes":          len(active_nodes),
            "n_zero_pop_nodes": len(zero_pop_nodes),
            "n_edges":          len(edges),
            "n_water_edges":    n_water,
            "n_triangles":      len(triangles),
        }

        print("\nSaving edges to database...")
        run_id = db.write_run(conn, geography, statefp, n_districts, params)
        db.write_edges(conn, run_id, active_nodes, edges, adjacent_pairs)

        state_name = _FIPS_TO_NAME.get(statefp, statefp)
        label = _GEOGRAPHY_LABELS.get(geography, geography)
        print(f"\n{state_name}: {len(nodes):,} {label}, {len(edges):,} edges saved.")
        print(f"\nNext steps:")
        print(f"  1. Open QGIS and filter redistrict_edges:")
        print(f'       "run_id" = {run_id} AND NOT "is_adjacent"')
        print(f"  2. Delete the non-adjacent edges you don't want.")
        print(f"  3. Run:  redistrict --continue {run_id}")

        return run_id

    finally:
        conn.close()


def continue_run(
    parent_run_id: int,
    formula: str | None = None,
    water_penalty: float | None = None,
    ncuts: int | None = None,
    niter: int | None = None,
) -> int:
    """
    Re-run using edges stored in redistrict_edges for parent_run_id.

    The user has opened the run in QGIS, filtered to is_adjacent=false,
    and deleted the non-adjacent edges they don't want. This function
    reads whatever edges remain, applies the water penalty to non-adjacent
    ones, checks connectivity, then runs METIS with the same parameters.
    Saves as a new run referencing the parent. No Urquhart recomputation.

    Optional keyword arguments override the stored params from the parent run.

    Returns the new run_id.
    """
    conn = db.connect()
    try:
        parent = db.fetch_run(conn, parent_run_id)
        geography   = parent["geography"]
        statefp     = parent["statefp"]
        n_districts = parent["n_districts"]
        params_orig = parent["params"]

        db.ensure_tables(conn, geography)

        print(f"\nLoading edges from run {parent_run_id}...")
        adj_geoid_pairs, non_adj_geoid_pairs = db.fetch_edges(conn, parent_run_id)
        n_orig_nonadj = params_orig.get("n_water_edges", len(non_adj_geoid_pairs))
        n_removed = n_orig_nonadj - len(non_adj_geoid_pairs)
        print(f"  {len(adj_geoid_pairs):,} adjacent edges.")
        print(f"  {len(non_adj_geoid_pairs):,} non-adjacent edges kept "
              f"({n_removed} removed by user).")

        nodes = db.fetch_nodes(conn, geography, statefp)
        active_nodes = [n for n in nodes if n["pop"] > 0]
        zero_pop_nodes = [n for n in nodes if n["pop"] == 0]
        if zero_pop_nodes:
            print(f"  {len(zero_pop_nodes):,} zero-pop nodes will be assigned post-METIS.")

        geoid_to_idx = {n["geoid"]: i for i, n in enumerate(active_nodes)}

        # Rebuild index-pair edge set from geoid pairs.
        all_geoid_pairs = adj_geoid_pairs | non_adj_geoid_pairs
        edges: set[tuple[int, int]] = set()
        for ga, gb in all_geoid_pairs:
            if ga in geoid_to_idx and gb in geoid_to_idx:
                i, j = geoid_to_idx[ga], geoid_to_idx[gb]
                edges.add((min(i, j), max(i, j)))

        # Connectivity check — auto-bridge any disconnected components.
        components = graph.check_connectivity(active_nodes, edges)
        if len(components) > 1:
            sizes = sorted((len(c) for c in components), reverse=True)
            print(f"\n  Graph has {len(components)} disconnected components "
                  f"(sizes: {sizes}). Auto-connecting...")
            bridges = graph.reconnect_components(active_nodes, components)
            for i, j in bridges:
                ga = min(active_nodes[i]["geoid"], active_nodes[j]["geoid"])
                gb = max(active_nodes[i]["geoid"], active_nodes[j]["geoid"])
                d = graph.haversine_km(
                    active_nodes[i]["lat"], active_nodes[i]["lon"],
                    active_nodes[j]["lat"], active_nodes[j]["lon"],
                )
                print(f"    Bridge: {ga} — {gb}  ({d:.1f} km, non-adjacent)")
                non_adj_geoid_pairs.add((ga, gb))
            edges |= bridges

        formula       = formula       if formula       is not None else params_orig.get("formula", "original")
        recursive     = params_orig.get("recursive", False)
        water_penalty = water_penalty if water_penalty is not None else params_orig.get("water_penalty", graph.WATER_PENALTY)
        ncuts         = ncuts         if ncuts         is not None else params_orig.get("ncuts", 10)
        niter         = niter         if niter         is not None else params_orig.get("niter", 20)

        print("Building METIS graph...")
        adj_lists, eweights, nweights = graph.build_metis_graph(
            active_nodes, edges, adj_geoid_pairs,
            water_penalty=water_penalty, formula=formula,
        )

        print(f"\nRunning PyMETIS: {len(active_nodes)} nodes -> {n_districts} districts "
              f"(ncuts={ncuts}, niter={niter}, recursive={recursive})...")
        edge_cut, membership = partition.partition(
            adj_lists, eweights, nweights, n_districts,
            ncuts=ncuts, niter=niter, recursive=recursive,
        )
        print(f"  Edge cut weight: {edge_cut:,}")

        geoid_to_district = {active_nodes[i]["geoid"]: membership[i]
                             for i in range(len(active_nodes))}

        # Post-assign zero-pop nodes to the district of their nearest active node.
        if zero_pop_nodes:
            print(f"  Assigning {len(zero_pop_nodes):,} zero-pop nodes to nearest district...")
            for zn in zero_pop_nodes:
                nearest_idx = graph.nearest_node(zn, active_nodes)
                nearest_geoid = active_nodes[nearest_idx]["geoid"]
                geoid_to_district[zn["geoid"]] = geoid_to_district[nearest_geoid]

        pop_per_district: dict[int, int] = {}
        for i, d in enumerate(membership):
            pop_per_district[d] = pop_per_district.get(d, 0) + active_nodes[i]["pop"]
        total_pop = sum(nweights)
        ideal = total_pop / n_districts
        print("\nDistrict populations:")
        for d in sorted(pop_per_district):
            pct_dev = 100 * (pop_per_district[d] - ideal) / ideal
            print(f"  District {d}: {pop_per_district[d]:,}  ({pct_dev:+.1f}% from ideal)")

        params = {
            **params_orig,
            "parent_run_id":        parent_run_id,
            "non_adj_edges_kept":   len(non_adj_geoid_pairs),
            "non_adj_edges_removed": n_removed,
            "edge_cut":             edge_cut,
            "n_edges":              len(edges),
        }

        print("\nWriting results to database...")
        run_id = db.write_run(conn, geography, statefp, n_districts, params)
        db.write_assignments(conn, run_id, geoid_to_district)
        db.write_district_geoms(conn, run_id, geography, geoid_to_district)
        db.write_edges(conn, run_id, active_nodes, edges, adj_geoid_pairs)
        db.update_run_params(conn, parent_run_id, {"continued_as": run_id})
        print(f"  Run ID: {run_id}  (parent: {parent_run_id})")

        state_name = _FIPS_TO_NAME.get(statefp, statefp)
        label = _GEOGRAPHY_LABELS.get(geography, geography)
        print(f"\nDone. {state_name}: {n_districts} districts from {len(nodes):,} {label}.")
        return run_id

    finally:
        conn.close()


def _parse_flag(flag: str, cast):
    """Return cast(value) for --flag value in sys.argv, or None if absent."""
    if flag in sys.argv:
        i = sys.argv.index(flag)
        try:
            return cast(sys.argv[i + 1])
        except (IndexError, ValueError):
            print(f"Usage: redistrict --continue <run_id> [{flag} <value>]")
            sys.exit(1)
    return None


def main() -> None:
    # Handle --continue <run_id> before entering interactive mode.
    if "--continue" in sys.argv:
        idx = sys.argv.index("--continue")
        try:
            run_id = int(sys.argv[idx + 1])
        except (IndexError, ValueError):
            print("Usage: redistrict --continue <run_id> [--formula <f>] "
                  "[--water-penalty <w>] [--ncuts <n>] [--niter <n>]")
            sys.exit(1)
        continue_run(
            run_id,
            formula=_parse_flag("--formula", str),
            water_penalty=_parse_flag("--water-penalty", float),
            ncuts=_parse_flag("--ncuts", int),
            niter=_parse_flag("--niter", int),
        )
        return

    geography: str = inquirer.select(
        message="Select geography level:",
        choices=[
            Choice(value="tracts",       name="Census tracts"),
            Choice(value="block_groups", name="Census block groups"),
            Choice(value="counties",     name="Counties"),
        ],
    ).execute()

    conn = db.connect()
    try:
        db.ensure_tables(conn, geography)
        state_choices = _state_choices(conn, geography)
    finally:
        conn.close()

    if not state_choices:
        label = _GEOGRAPHY_LABELS.get(geography, geography)
        print(f"No states found for {label}. Run the aggregation scripts first.")
        sys.exit(1)

    statefp: str = inquirer.select(
        message="Select a state:",
        choices=state_choices,
    ).execute()

    n_districts_str: str = inquirer.text(
        message="Number of districts:",
        validate=lambda x: x.isdigit() and int(x) >= 2,
        invalid_message="Enter an integer >= 2.",
    ).execute()

    run(statefp, geography, int(n_districts_str))
