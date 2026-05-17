#!/usr/bin/env python3
"""
Run the redistricting pipeline for every state in the DB (except NY/36),
using run 37's settings: block_groups, uniform formula, water_penalty=40.0,
ncuts=10, niter=20, k-way, zero-pop nodes excluded and post-assigned.

House seat counts from 2020 Census apportionment.
States with 1 seat are skipped (nothing to partition).
One DB connection per state; results and GeoJSONs written to output/all_states/.
"""
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))

from tqdm import tqdm

from redistrict import db, graph, partition

FORMULA       = "uniform"
WATER_PENALTY = 40.0
NCUTS         = 10
NITER         = 20
RECURSIVE     = False
GEOGRAPHY     = "block_groups"

# 2020 Census apportionment — FIPS -> House seats.
# States with 1 seat are skipped (no partition possible).
SEATS: dict[str, int] = {
    "01": 7,   # Alabama
    "04": 9,   # Arizona
    "05": 4,   # Arkansas
    "06": 52,  # California
    "08": 8,   # Colorado
    "09": 5,   # Connecticut
    "10": 1,   # Delaware        — skip
    "11": 1,   # DC              — skip
    "12": 28,  # Florida
    "13": 14,  # Georgia
    "16": 2,   # Idaho
    "17": 17,  # Illinois
    "18": 9,   # Indiana
    "19": 4,   # Iowa
    "20": 4,   # Kansas
    "21": 6,   # Kentucky
    "22": 6,   # Louisiana
    "23": 2,   # Maine
    "24": 8,   # Maryland
    "25": 9,   # Massachusetts
    "26": 13,  # Michigan
    "27": 8,   # Minnesota
    "28": 4,   # Mississippi
    "29": 8,   # Missouri
    "30": 1,   # Montana         — skip
    "31": 3,   # Nebraska
    "32": 4,   # Nevada
    "33": 2,   # New Hampshire
    "34": 12,  # New Jersey
    "35": 3,   # New Mexico
    # "36" NY skipped per user request
    "37": 14,  # North Carolina
    "38": 1,   # North Dakota    — skip
    "39": 15,  # Ohio
    "40": 5,   # Oklahoma
    "41": 6,   # Oregon
    "42": 17,  # Pennsylvania
    "44": 2,   # Rhode Island
    "45": 7,   # South Carolina
    "46": 1,   # South Dakota    — skip
    "47": 9,   # Tennessee
    "48": 38,  # Texas
    "49": 4,   # Utah
    "50": 1,   # Vermont         — skip
    "51": 11,  # Virginia
    "53": 10,  # Washington
    "54": 2,   # West Virginia
    "55": 8,   # Wisconsin
    "56": 1,   # Wyoming         — skip
}

_FIPS_TO_NAME = {
    "01": "Alabama",       "04": "Arizona",        "05": "Arkansas",
    "06": "California",    "08": "Colorado",        "09": "Connecticut",
    "10": "Delaware",      "11": "DC",              "12": "Florida",
    "13": "Georgia",       "16": "Idaho",           "17": "Illinois",
    "18": "Indiana",       "19": "Iowa",            "20": "Kansas",
    "21": "Kentucky",      "22": "Louisiana",       "23": "Maine",
    "24": "Maryland",      "25": "Massachusetts",   "26": "Michigan",
    "27": "Minnesota",     "28": "Mississippi",     "29": "Missouri",
    "30": "Montana",       "31": "Nebraska",        "32": "Nevada",
    "33": "New Hampshire", "34": "New Jersey",      "35": "New Mexico",
    "37": "North Carolina","38": "North Dakota",    "39": "Ohio",
    "40": "Oklahoma",      "41": "Oregon",          "42": "Pennsylvania",
    "44": "Rhode Island",  "45": "South Carolina",  "46": "South Dakota",
    "47": "Tennessee",     "48": "Texas",           "49": "Utah",
    "50": "Vermont",       "51": "Virginia",        "53": "Washington",
    "54": "West Virginia", "55": "Wisconsin",       "56": "Wyoming",
}

OUTPUT_DIR = os.path.join(os.path.dirname(__file__), "..", "output", "all_states")


def run_state(statefp: str, n_districts: int) -> dict:
    """
    Run the full pipeline for one state. Returns a summary dict.
    Raises on error so the caller can log and continue.
    """
    state_name = _FIPS_TO_NAME.get(statefp, statefp)
    conn = db.connect()
    try:
        db.ensure_tables(conn, GEOGRAPHY)

        nodes = db.fetch_nodes(conn, GEOGRAPHY, statefp)
        if not nodes:
            raise RuntimeError(f"No nodes found for {statefp}")

        active_nodes = [n for n in nodes if n["pop"] > 0]
        zero_pop_nodes = [n for n in nodes if n["pop"] == 0]

        if n_districts >= len(active_nodes):
            raise RuntimeError(
                f"n_districts ({n_districts}) >= active nodes ({len(active_nodes)})"
            )

        # Adjacency (usually already cached after first run)
        geoids = [n["geoid"] for n in nodes]
        missing = db.get_missing_adjacency_geoids(conn, GEOGRAPHY, geoids)
        if missing:
            print(f"  Computing adjacency for {len(missing):,} geoids...")
            bar = tqdm(total=len(missing), unit="bg", desc="  Adjacency", leave=False)
            db.compute_and_store_adjacency(
                conn, GEOGRAPHY, missing,
                progress_callback=lambda d, t, g: bar.update(1),
            )
            bar.close()
        adjacent_pairs = db.fetch_adjacency(conn, GEOGRAPHY, geoids)

        triangles = graph.spherical_delaunay_triangles(active_nodes)
        edges = graph.urquhart_edges(active_nodes, triangles)

        n_water = sum(
            1 for i, j in edges
            if (min(active_nodes[i]["geoid"], active_nodes[j]["geoid"]),
                max(active_nodes[i]["geoid"], active_nodes[j]["geoid"])) not in adjacent_pairs
        )

        adj_lists, eweights, nweights = graph.build_metis_graph(
            active_nodes, edges, adjacent_pairs,
            water_penalty=WATER_PENALTY, formula=FORMULA,
        )

        edge_cut, membership = partition.partition(
            adj_lists, eweights, nweights, n_districts,
            ncuts=NCUTS, niter=NITER, recursive=RECURSIVE,
        )

        geoid_to_district = {active_nodes[i]["geoid"]: membership[i]
                             for i in range(len(active_nodes))}

        for zn in zero_pop_nodes:
            nearest_idx = graph.nearest_node(zn, active_nodes)
            geoid_to_district[zn["geoid"]] = geoid_to_district[active_nodes[nearest_idx]["geoid"]]

        pop_per_district: dict[int, int] = {}
        for i, d in enumerate(membership):
            pop_per_district[d] = pop_per_district.get(d, 0) + active_nodes[i]["pop"]

        total_pop = sum(nweights)
        ideal = total_pop / n_districts
        worst = max(abs(100 * (p - ideal) / ideal) for p in pop_per_district.values())

        params = {
            "status":            "complete",
            "formula":           FORMULA,
            "water_penalty":     WATER_PENALTY,
            "ncuts":             NCUTS,
            "niter":             NITER,
            "recursive":         RECURSIVE,
            "n_nodes":           len(active_nodes),
            "n_zero_pop_nodes":  len(zero_pop_nodes),
            "n_edges":           len(edges),
            "n_water_edges":     n_water,
            "n_triangles":       len(triangles),
            "edge_cut":          edge_cut,
            "edge_weight_scale": graph.EDGE_WEIGHT_SCALE,
        }

        run_id = db.write_run(conn, GEOGRAPHY, statefp, n_districts, params)
        db.write_assignments(conn, run_id, geoid_to_district)
        db.write_district_geoms(conn, run_id, GEOGRAPHY, geoid_to_district)
        db.write_edges(conn, run_id, active_nodes, edges, adjacent_pairs)

        slug = state_name.lower().replace(" ", "_")
        geojson_path = os.path.join(OUTPUT_DIR, f"{slug}_run{run_id}.geojson")
        db.export_geojson(conn, run_id, geojson_path)

        log_lines = [
            f"# {state_name} — block_groups — {n_districts} districts — Run {run_id}",
            f"",
            f"water_penalty={WATER_PENALTY}  ncuts={NCUTS}  niter={NITER}  formula={FORMULA}",
            f"nodes={len(active_nodes):,}  zero_pop={len(zero_pop_nodes):,}  "
            f"edges={len(edges):,}  water_edges={n_water:,}",
            f"districts={n_districts}  ideal={ideal:,.0f}  worst_deviation={worst:.1f}%",
            f"",
            f"| District | Population | Deviation |",
            f"|----------|-----------|-----------|",
        ]
        for d in sorted(pop_per_district):
            pop = pop_per_district[d]
            pct = 100 * (pop - ideal) / ideal
            log_lines.append(f"| {d} | {pop:,} | {pct:+.1f}% |")
        log_path = os.path.join(OUTPUT_DIR, f"{slug}_run{run_id}_deviations.md")
        with open(log_path, "w") as fh:
            fh.write("\n".join(log_lines) + "\n")

        return {
            "statefp":      statefp,
            "state_name":   state_name,
            "run_id":       run_id,
            "n_nodes":      len(active_nodes),
            "n_zero_pop":   len(zero_pop_nodes),
            "n_districts":  n_districts,
            "worst_dev":    worst,
            "geojson":      geojson_path,
            "log":          log_path,
        }
    finally:
        conn.close()


def main() -> None:
    os.makedirs(OUTPUT_DIR, exist_ok=True)

    to_run = [(fips, seats) for fips, seats in sorted(SEATS.items()) if seats >= 2]
    skipped_single = [fips for fips, seats in SEATS.items() if seats < 2]

    print(f"States to run: {len(to_run)}")
    if skipped_single:
        names = [_FIPS_TO_NAME.get(f, f) for f in sorted(skipped_single)]
        print(f"Skipped (1 seat, nothing to partition): {', '.join(names)}")
    print()

    results: list[dict] = []
    failures: list[tuple[str, str]] = []

    for statefp, n_districts in to_run:
        state_name = _FIPS_TO_NAME.get(statefp, statefp)
        print(f"[{statefp}] {state_name} — {n_districts} districts")
        try:
            r = run_state(statefp, n_districts)
            results.append(r)
            print(f"  run_id={r['run_id']}  nodes={r['n_nodes']:,}  "
                  f"zero_pop={r['n_zero_pop']:,}  worst_dev={r['worst_dev']:.1f}%")
        except Exception as exc:
            print(f"  FAILED: {exc}")
            failures.append((statefp, str(exc)))

    # Summary log
    summary_path = os.path.join(OUTPUT_DIR, "summary.md")
    lines = [
        "# All-States Redistricting — Summary",
        f"",
        f"formula={FORMULA}  water_penalty={WATER_PENALTY}  ncuts={NCUTS}  niter={NITER}",
        f"geography=block_groups  zero-pop nodes excluded and post-assigned",
        f"",
        f"| State | FIPS | Districts | Nodes | Zero-pop | Worst Dev | Run ID |",
        f"|-------|------|-----------|-------|----------|-----------|--------|",
    ]
    for r in results:
        lines.append(
            f"| {r['state_name']} | {r['statefp']} | {r['n_districts']} | "
            f"{r['n_nodes']:,} | {r['n_zero_pop']:,} | {r['worst_dev']:.1f}% | {r['run_id']} |"
        )
    if failures:
        lines += ["", "## Failures", ""]
        for fips, msg in failures:
            lines.append(f"- {_FIPS_TO_NAME.get(fips, fips)} ({fips}): {msg}")
    with open(summary_path, "w") as fh:
        fh.write("\n".join(lines) + "\n")

    print(f"\nDone. {len(results)} states completed, {len(failures)} failed.")
    print(f"Summary -> {summary_path}")


if __name__ == "__main__":
    main()
