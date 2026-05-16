#!/usr/bin/env python3
"""
High-quality run of the uniform edge-weight formula for New York census tracts.
State: New York (FIPS 36), 26 districts, ncuts=30, niter=50.

Exports GeoJSON and writes a deviation log to output/ny_experiments/.
"""
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))

from redistrict import cli, db

NY_FIPS = "36"
GEOGRAPHY = "tracts"
N_DISTRICTS = 26
NCUTS = 30
NITER = 50
FORMULA = "uniform"

OUTPUT_DIR = os.path.join(
    os.path.dirname(__file__), "..", "output", "ny_experiments"
)


def main() -> None:
    os.makedirs(OUTPUT_DIR, exist_ok=True)

    run_id = cli.run(
        NY_FIPS, GEOGRAPHY, N_DISTRICTS,
        ncuts=NCUTS, niter=NITER, formula=FORMULA,
    )

    conn = db.connect()
    try:
        # GeoJSON
        geojson_path = os.path.join(OUTPUT_DIR, f"run_{run_id}_{FORMULA}_hq.geojson")
        db.export_geojson(conn, run_id, geojson_path)
        print(f"\nGeoJSON -> {geojson_path}")

        # Deviation log
        pops = db.fetch_district_populations(conn, run_id)
    finally:
        conn.close()

    total_pop = sum(pops.values())
    ideal = total_pop / N_DISTRICTS

    log_lines = [
        f"# NY Tracts — Uniform HQ — Run {run_id}",
        f"",
        f"ncuts={NCUTS}  niter={NITER}  formula={FORMULA}",
        f"nodes=5411  districts={N_DISTRICTS}  ideal={ideal:,.0f}",
        f"",
        f"| District | Population | Deviation |",
        f"|----------|-----------|-----------|",
    ]
    for d in sorted(pops):
        pop = pops[d]
        pct = 100 * (pop - ideal) / ideal
        sign = "+" if pct >= 0 else ""
        log_lines.append(f"| {d} | {pop:,} | {sign}{pct:.1f}% |")

    worst = max(abs(100 * (p - ideal) / ideal) for p in pops.values())
    log_lines += [f"", f"Worst deviation: {worst:.1f}%"]

    log_path = os.path.join(OUTPUT_DIR, f"run_{run_id}_{FORMULA}_hq_deviations.md")
    with open(log_path, "w") as fh:
        fh.write("\n".join(log_lines) + "\n")
    print(f"Deviations -> {log_path}")


if __name__ == "__main__":
    main()
