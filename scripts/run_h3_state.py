#!/usr/bin/env python3
"""
H3-based redistricting from census blocks for a single state.

Usage:
    uv run --env-file .env python scripts/run_h3_state.py <statefp> <n_districts>

Example:
    uv run --env-file .env python scripts/run_h3_state.py 44 2   # Rhode Island, 2 districts
    uv run --env-file .env python scripts/run_h3_state.py 36 26  # New York, 26 districts

Pipeline:
  1. Fetch census blocks from DB (blocks_2020 table)
  2. Compute the 99th-percentile-by-pop-mass threshold
  3. Assign each block centroid to H3 resolution 15
  4. Bottom-up aggregate: find coarsest hex ≤ threshold per group
  5. Compute population-weighted centroid per aggregated hex
  6. Build adjacency graph via H3 res-15 neighbour lookup
  7. Run PyMETIS (uniform edge weights, pop node weights)
  8. Disaggregate: assign every block to its hex's district
  9. Write run record, assignments, district geometries to DB
  10. Export GeoJSON to output/h3_runs/
"""
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))

from redistrict import db, graph, partition
from redistrict import h3_graph

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
    "36": "New York",      "37": "North Carolina",  "38": "North Dakota",
    "39": "Ohio",          "40": "Oklahoma",        "41": "Oregon",
    "42": "Pennsylvania",  "44": "Rhode Island",    "45": "South Carolina",
    "46": "South Dakota",  "47": "Tennessee",       "48": "Texas",
    "49": "Utah",          "50": "Vermont",         "51": "Virginia",
    "53": "Washington",    "54": "West Virginia",   "55": "Wisconsin",
    "56": "Wyoming",
}

OUTPUT_DIR = os.path.join(os.path.dirname(__file__), "..", "output", "h3_runs")
NCUTS = 10
NITER = 20
RECURSIVE = False


def main(statefp: str, n_districts: int) -> int:
    os.makedirs(OUTPUT_DIR, exist_ok=True)
    state_name = _FIPS_TO_NAME.get(statefp, statefp)
    conn = db.connect()
    try:
        db.ensure_tables(conn, "blocks")

        print(f"\nFetching census blocks for {state_name} (statefp={statefp})...")
        blocks = db.fetch_nodes(conn, "blocks", statefp)
        if not blocks:
            print("No blocks found. Is blocks_2020 loaded?")
            return
        print(f"  {len(blocks):,} blocks loaded.")

        active_blocks = [b for b in blocks if int(b["pop"]) > 0]
        zero_pop_blocks = [b for b in blocks if int(b["pop"]) == 0]
        total_pop = sum(int(b["pop"]) for b in active_blocks)
        print(f"  {len(active_blocks):,} active blocks, "
              f"{len(zero_pop_blocks):,} zero-pop (post-assigned after METIS).")
        print(f"  Total population: {total_pop:,}")

        # --- Threshold ---
        threshold = h3_graph.compute_threshold(active_blocks)
        above = sum(1 for b in active_blocks if int(b["pop"]) > threshold)
        print(f"\nThreshold (99th pct by pop mass): {threshold:,}")
        print(f"  {above:,} blocks individually above threshold (stay at res-15).")

        # --- H3 res-15 assignment (active only) ---
        print("\nAssigning blocks to H3 resolution-15...")
        geoid_to_res15 = h3_graph.assign_h3_res15(active_blocks)
        print(f"  {len(set(geoid_to_res15.values())):,} unique res-15 cells.")

        # --- Bottom-up aggregation ---
        print("Aggregating H3 cells bottom-up...")
        block_pops = {b["geoid"]: int(b["pop"]) for b in active_blocks}
        geoid_to_cell = h3_graph.aggregate_h3_cells(
            geoid_to_res15, block_pops, threshold,
        )
        n_cells = len(set(geoid_to_cell.values()))
        print(f"  {n_cells:,} aggregated cells (from {len(active_blocks):,} active blocks).")

        if n_districts >= n_cells:
            print(f"Error: n_districts ({n_districts}) >= cell count ({n_cells}).")
            return

        # --- Population-weighted centroids ---
        print("Computing population-weighted centroids...")
        blocks_by_geoid = {b["geoid"]: b for b in blocks}
        cell_nodes = h3_graph.weighted_centroids(geoid_to_cell, blocks_by_geoid)
        print(f"  {len(cell_nodes):,} cell nodes.")

        # --- Adjacency ---
        print("Building H3 adjacency graph...")
        adjacency = h3_graph.build_h3_adjacency(geoid_to_cell)
        print(f"  {len(adjacency):,} edges.")

        # --- Connectivity check ---
        components = h3_graph.check_connectivity(cell_nodes, adjacency)
        if len(components) > 1:
            print(f"  WARNING: graph has {len(components)} disconnected components. "
                  f"Largest: {max(len(c) for c in components):,} cells.")

        # --- METIS ---
        adj_lists, eweights, nweights = h3_graph.build_metis_graph(cell_nodes, adjacency)
        print(f"\nRunning PyMETIS: {len(cell_nodes):,} nodes -> {n_districts} districts "
              f"(ncuts={NCUTS}, niter={NITER})...")
        edge_cut, membership = partition.partition(
            adj_lists, eweights, nweights, n_districts,
            ncuts=NCUTS, niter=NITER, recursive=RECURSIVE,
        )
        print(f"  Edge cut: {edge_cut:,}")

        # --- Disaggregate: cell -> district, then block -> district ---
        cell_to_district = {cell_nodes[i]["cell"]: membership[i]
                            for i in range(len(cell_nodes))}
        geoid_to_district = {g: cell_to_district[c] for g, c in geoid_to_cell.items()}

        # --- Post-assign zero-pop blocks to nearest active block's district ---
        if zero_pop_blocks:
            print(f"  Assigning {len(zero_pop_blocks):,} zero-pop blocks to nearest district...")
            for zb in zero_pop_blocks:
                nearest_idx = graph.nearest_node(zb, active_blocks)
                geoid_to_district[zb["geoid"]] = geoid_to_district[active_blocks[nearest_idx]["geoid"]]

        # --- Population stats ---
        pop_per_district: dict[int, int] = {}
        for b in active_blocks:
            d = geoid_to_district[b["geoid"]]
            pop_per_district[d] = pop_per_district.get(d, 0) + int(b["pop"])
        ideal = total_pop / n_districts
        print("\nDistrict populations:")
        for d in sorted(pop_per_district):
            pct = 100 * (pop_per_district[d] - ideal) / ideal
            print(f"  District {d}: {pop_per_district[d]:,}  ({pct:+.1f}%)")
        worst = max(abs(100 * (p - ideal) / ideal) for p in pop_per_district.values())

        # --- Write to DB ---
        params = {
            "status":          "complete",
            "method":          "h3_blocks",
            "threshold":       threshold,
            "n_blocks":        len(blocks),
            "n_h3_cells":      n_cells,
            "n_edges":         len(adjacency),
            "n_components":    len(components),
            "edge_cut":        edge_cut,
            "ncuts":           NCUTS,
            "niter":           NITER,
            "recursive":       RECURSIVE,
        }
        print("\nWriting results to database...")
        run_id = db.write_run(conn, "blocks", statefp, n_districts, params)
        db.write_assignments(conn, run_id, geoid_to_district)
        db.write_district_geoms(conn, run_id, "blocks", geoid_to_district)
        print(f"  Run ID: {run_id}")

        # --- GeoJSON export ---
        slug = state_name.lower().replace(" ", "_")
        geojson_path = os.path.join(OUTPUT_DIR, f"{slug}_h3_run{run_id}.geojson")
        db.export_geojson(conn, run_id, geojson_path)
        print(f"  GeoJSON -> {geojson_path}")

        # --- Deviation log ---
        log_path = os.path.join(OUTPUT_DIR, f"{slug}_h3_run{run_id}_deviations.md")
        lines = [
            f"# {state_name} — H3 blocks — {n_districts} districts — Run {run_id}",
            f"",
            f"method=h3_blocks  threshold={threshold}  ncuts={NCUTS}  niter={NITER}",
            f"blocks={len(blocks):,}  h3_cells={n_cells:,}  edges={len(adjacency):,}",
            f"districts={n_districts}  ideal={ideal:,.0f}  worst_deviation={worst:.1f}%",
            f"",
            f"| District | Population | Deviation |",
            f"|----------|-----------|-----------|",
        ]
        for d in sorted(pop_per_district):
            pop = pop_per_district[d]
            pct = 100 * (pop - ideal) / ideal
            lines.append(f"| {d} | {pop:,} | {pct:+.1f}% |")
        with open(log_path, "w") as fh:
            fh.write("\n".join(lines) + "\n")
        print(f"  Deviations -> {log_path}")

        return run_id

    finally:
        conn.close()


if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: run_h3_state.py <statefp> <n_districts>")
        sys.exit(1)
    statefp = sys.argv[1]
    try:
        n_districts = int(sys.argv[2])
    except ValueError:
        print("n_districts must be an integer")
        sys.exit(1)
    main(statefp, n_districts)
