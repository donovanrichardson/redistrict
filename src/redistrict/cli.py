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

        if n_districts >= len(nodes):
            print(f"Error: n_districts ({n_districts}) must be less than "
                  f"node count ({len(nodes)}).")
            return

        adjacent_pairs = _ensure_adjacency(conn, geography, nodes)
        print(f"  {len(adjacent_pairs):,} adjacency pairs loaded.")

        print("\nBuilding spherical Delaunay triangulation...")
        triangles = graph.spherical_delaunay_triangles(nodes)
        print(f"  {len(triangles):,} triangles.")

        print("Deriving Urquhart graph...")
        edges = graph.urquhart_edges(nodes, triangles)
        print(f"  {len(edges):,} edges.")

        n_water = sum(
            1 for i, j in edges
            if (min(nodes[i]["geoid"], nodes[j]["geoid"]),
                max(nodes[i]["geoid"], nodes[j]["geoid"])) not in adjacent_pairs
        )
        print(f"  {n_water:,} non-adjacent (water) edges.")

        print("Building METIS graph...")
        adj_lists, eweights, nweights = graph.build_metis_graph(
            nodes, edges, adjacent_pairs,
            water_penalty=water_penalty, formula=formula,
        )

        print(f"\nRunning PyMETIS: {len(nodes)} nodes -> {n_districts} districts "
              f"(ncuts={ncuts}, niter={niter}, recursive={recursive})...")
        edge_cut, membership = partition.partition(
            adj_lists, eweights, nweights, n_districts,
            ncuts=ncuts, niter=niter, recursive=recursive,
        )
        print(f"  Edge cut weight: {edge_cut:,}")

        geoid_to_district = {nodes[i]["geoid"]: membership[i] for i in range(len(nodes))}

        pop_per_district: dict[int, int] = {}
        for i, d in enumerate(membership):
            pop_per_district[d] = pop_per_district.get(d, 0) + nodes[i]["pop"]
        total_pop = sum(nweights)
        ideal = total_pop / n_districts
        print("\nDistrict populations:")
        for d in sorted(pop_per_district):
            pct_dev = 100 * (pop_per_district[d] - ideal) / ideal
            print(f"  District {d}: {pop_per_district[d]:,}  ({pct_dev:+.1f}% from ideal)")

        params = {
            "water_penalty":    water_penalty,
            "edge_weight_scale": graph.EDGE_WEIGHT_SCALE,
            "formula":          formula,
            "recursive":        recursive,
            "ncuts":            ncuts,
            "niter":            niter,
            "edge_cut":         edge_cut,
            "n_nodes":          len(nodes),
            "n_edges":          len(edges),
            "n_water_edges":    n_water,
            "n_triangles":      len(triangles),
        }

        print("\nWriting results to database...")
        run_id = db.write_run(conn, geography, statefp, n_districts, params)
        db.write_assignments(conn, run_id, geoid_to_district)
        db.write_district_geoms(conn, run_id, geography, geoid_to_district)
        db.write_edges(conn, run_id, nodes, edges, adjacent_pairs)
        print(f"  Run ID: {run_id}  (redistrict_runs table)")
        print(f"  {n_water:,} non-adjacent edges saved (filter in QGIS, then --continue {run_id})")

        state_name = _FIPS_TO_NAME.get(statefp, statefp)
        label = _GEOGRAPHY_LABELS.get(geography, geography)
        print(f"\nDone. {state_name}: {n_districts} districts from {len(nodes):,} {label}.")

        return run_id

    finally:
        conn.close()


def continue_run(parent_run_id: int) -> int:
    """
    Re-run using edges stored in redistrict_edges for parent_run_id.

    The user has opened the run in QGIS, filtered to is_adjacent=false,
    and deleted the non-adjacent edges they don't want. This function
    reads whatever edges remain, applies the water penalty to non-adjacent
    ones, checks connectivity, then runs METIS with the same parameters.
    Saves as a new run referencing the parent. No Urquhart recomputation.

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
        geoid_to_idx = {n["geoid"]: i for i, n in enumerate(nodes)}

        # Rebuild index-pair edge set from geoid pairs.
        all_geoid_pairs = adj_geoid_pairs | non_adj_geoid_pairs
        edges: set[tuple[int, int]] = set()
        for ga, gb in all_geoid_pairs:
            if ga in geoid_to_idx and gb in geoid_to_idx:
                i, j = geoid_to_idx[ga], geoid_to_idx[gb]
                edges.add((min(i, j), max(i, j)))

        # Connectivity check.
        components = graph.check_connectivity(nodes, edges)
        if len(components) > 1:
            sizes = sorted((len(c) for c in components), reverse=True)
            print(f"\nWARNING: graph has {len(components)} disconnected components "
                  f"(sizes: {sizes}).")
            print("  METIS contig enforcement may produce non-contiguous districts.")
            answer = input("  Proceed anyway? [y/N] ").strip().lower()
            if answer != "y":
                print("Aborted.")
                return -1

        formula       = params_orig.get("formula", "original")
        recursive     = params_orig.get("recursive", False)
        water_penalty = params_orig.get("water_penalty", graph.WATER_PENALTY)
        ncuts         = params_orig.get("ncuts", 10)
        niter         = params_orig.get("niter", 20)

        print("Building METIS graph...")
        adj_lists, eweights, nweights = graph.build_metis_graph(
            nodes, edges, adj_geoid_pairs,
            water_penalty=water_penalty, formula=formula,
        )

        print(f"\nRunning PyMETIS: {len(nodes)} nodes -> {n_districts} districts "
              f"(ncuts={ncuts}, niter={niter}, recursive={recursive})...")
        edge_cut, membership = partition.partition(
            adj_lists, eweights, nweights, n_districts,
            ncuts=ncuts, niter=niter, recursive=recursive,
        )
        print(f"  Edge cut weight: {edge_cut:,}")

        geoid_to_district = {nodes[i]["geoid"]: membership[i] for i in range(len(nodes))}

        pop_per_district: dict[int, int] = {}
        for i, d in enumerate(membership):
            pop_per_district[d] = pop_per_district.get(d, 0) + nodes[i]["pop"]
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
        db.write_edges(conn, run_id, nodes, edges, adj_geoid_pairs)
        print(f"  Run ID: {run_id}  (parent: {parent_run_id})")

        state_name = _FIPS_TO_NAME.get(statefp, statefp)
        label = _GEOGRAPHY_LABELS.get(geography, geography)
        print(f"\nDone. {state_name}: {n_districts} districts from {len(nodes):,} {label}.")
        return run_id

    finally:
        conn.close()


def main() -> None:
    # Handle --continue <run_id> before entering interactive mode.
    if "--continue" in sys.argv:
        idx = sys.argv.index("--continue")
        try:
            run_id = int(sys.argv[idx + 1])
        except (IndexError, ValueError):
            print("Usage: redistrict --continue <run_id>")
            sys.exit(1)
        continue_run(run_id)
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
