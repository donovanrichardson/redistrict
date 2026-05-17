#!/usr/bin/env python3
"""
Reproduce run 37 settings (NY block_groups, 26 districts, uniform, water_penalty=40.0)
with zero-pop node exclusion and post-METIS assignment.

Full one-shot pipeline — no QGIS curation step.
"""
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))

from tqdm import tqdm

from redistrict import db, graph, partition

NY_FIPS      = "36"
GEOGRAPHY    = "block_groups"
N_DISTRICTS  = 26
FORMULA      = "uniform"
WATER_PENALTY = 40.0
NCUTS        = 10
NITER        = 20
RECURSIVE    = False

OUTPUT_DIR = os.path.join(os.path.dirname(__file__), "..", "output", "ny_experiments")


def main() -> None:
    os.makedirs(OUTPUT_DIR, exist_ok=True)
    conn = db.connect()
    try:
        db.ensure_tables(conn, GEOGRAPHY)

        print(f"Fetching {GEOGRAPHY} for NY (statefp={NY_FIPS})...")
        nodes = db.fetch_nodes(conn, GEOGRAPHY, NY_FIPS)
        print(f"  {len(nodes):,} nodes loaded.")

        active_nodes = [n for n in nodes if n["pop"] > 0]
        zero_pop_nodes = [n for n in nodes if n["pop"] == 0]
        print(f"  {len(active_nodes):,} active nodes, {len(zero_pop_nodes):,} zero-pop (will be assigned post-METIS).")

        # Adjacency
        geoids = [n["geoid"] for n in nodes]
        missing = db.get_missing_adjacency_geoids(conn, GEOGRAPHY, geoids)
        if missing:
            print(f"Computing adjacency for {len(missing):,} geoids...")
            bar = tqdm(total=len(missing), unit=GEOGRAPHY, desc="Adjacency")
            db.compute_and_store_adjacency(conn, GEOGRAPHY, missing,
                                           progress_callback=lambda d, t, g: bar.update(1))
            bar.close()
        adjacent_pairs = db.fetch_adjacency(conn, GEOGRAPHY, geoids)
        print(f"  {len(adjacent_pairs):,} adjacency pairs loaded.")

        print("Building spherical Delaunay triangulation...")
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

        print("Building METIS graph...")
        adj_lists, eweights, nweights = graph.build_metis_graph(
            active_nodes, edges, adjacent_pairs,
            water_penalty=WATER_PENALTY, formula=FORMULA,
        )

        print(f"Running PyMETIS: {len(active_nodes):,} nodes -> {N_DISTRICTS} districts "
              f"(ncuts={NCUTS}, niter={NITER})...")
        edge_cut, membership = partition.partition(
            adj_lists, eweights, nweights, N_DISTRICTS,
            ncuts=NCUTS, niter=NITER, recursive=RECURSIVE,
        )
        print(f"  Edge cut weight: {edge_cut:,}")

        geoid_to_district = {active_nodes[i]["geoid"]: membership[i]
                             for i in range(len(active_nodes))}

        if zero_pop_nodes:
            print(f"  Assigning {len(zero_pop_nodes):,} zero-pop nodes to nearest district...")
            for zn in zero_pop_nodes:
                nearest_idx = graph.nearest_node(zn, active_nodes)
                geoid_to_district[zn["geoid"]] = geoid_to_district[active_nodes[nearest_idx]["geoid"]]

        pop_per_district: dict[int, int] = {}
        for i, d in enumerate(membership):
            pop_per_district[d] = pop_per_district.get(d, 0) + active_nodes[i]["pop"]

        total_pop = sum(nweights)
        ideal = total_pop / N_DISTRICTS
        print("\nDistrict populations:")
        for d in sorted(pop_per_district):
            pct = 100 * (pop_per_district[d] - ideal) / ideal
            print(f"  District {d}: {pop_per_district[d]:,}  ({pct:+.1f}%)")

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
            "repro_of":          37,
        }

        print("\nWriting results to database...")
        run_id = db.write_run(conn, GEOGRAPHY, NY_FIPS, N_DISTRICTS, params)
        db.write_assignments(conn, run_id, geoid_to_district)
        db.write_district_geoms(conn, run_id, GEOGRAPHY, geoid_to_district)
        db.write_edges(conn, run_id, active_nodes, edges, adjacent_pairs)
        print(f"  Run ID: {run_id}")

        geojson_path = os.path.join(OUTPUT_DIR, f"run_{run_id}_run37_repro_zero_pop.geojson")
        db.export_geojson(conn, run_id, geojson_path)
        print(f"  GeoJSON -> {geojson_path}")

        log_path = os.path.join(OUTPUT_DIR, f"run_{run_id}_run37_repro_zero_pop_deviations.md")
        worst = max(abs(100 * (p - ideal) / ideal) for p in pop_per_district.values())
        lines = [
            f"# NY block_groups k-way uniform — water8x — zero-pop excluded — Run {run_id}",
            f"",
            f"Reproduces run 37 settings with zero-pop node exclusion.",
            f"zero_pop_nodes={len(zero_pop_nodes)} post-assigned to nearest active node.",
            f"",
            f"water_penalty={WATER_PENALTY}  ncuts={NCUTS}  niter={NITER}  formula={FORMULA}",
            f"districts={N_DISTRICTS}  ideal={ideal:,.0f}",
            f"",
            f"| District | Population | Deviation |",
            f"|----------|-----------|-----------|",
        ]
        for d in sorted(pop_per_district):
            pop = pop_per_district[d]
            pct = 100 * (pop - ideal) / ideal
            lines.append(f"| {d} | {pop:,} | {pct:+.1f}% |")
        lines += ["", f"Max deviation: {worst:.1f}%",
                  f"GeoJSON: run_{run_id}_run37_repro_zero_pop.geojson"]
        with open(log_path, "w") as fh:
            fh.write("\n".join(lines) + "\n")
        print(f"  Deviations -> {log_path}")

    finally:
        conn.close()


if __name__ == "__main__":
    main()
