#!/usr/bin/env python3
"""
Run the uniform edge-weight formula across four geography/algorithm variants
for New York (FIPS 36, 26 districts):

  1. tracts    + k-way      (recursive=False)
  2. tracts    + recursive  (recursive=True)
  3. block_groups + k-way
  4. block_groups + recursive

Exports GeoJSONs and deviation logs to output/ny_experiments/.
"""
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))

from redistrict import cli, db

NY_FIPS = "36"
N_DISTRICTS = 26
NCUTS = 10
NITER = 20
FORMULA = "uniform"

VARIANTS = [
    ("tracts",       False, "tracts_kway"),
    ("tracts",       True,  "tracts_recursive"),
    ("block_groups", False, "block_groups_kway"),
    ("block_groups", True,  "block_groups_recursive"),
]

OUTPUT_DIR = os.path.join(
    os.path.dirname(__file__), "..", "output", "ny_experiments"
)


def write_deviation_log(path: str, run_id: int, label: str, pops: dict[int, int]) -> None:
    total = sum(pops.values())
    ideal = total / N_DISTRICTS
    lines = [
        f"# NY {label} — Uniform — Run {run_id}",
        f"",
        f"ncuts={NCUTS}  niter={NITER}  formula={FORMULA}  label={label}",
        f"districts={N_DISTRICTS}  ideal={ideal:,.0f}",
        f"",
        f"| District | Population | Deviation |",
        f"|----------|-----------|-----------|",
    ]
    for d in sorted(pops):
        pop = pops[d]
        pct = 100 * (pop - ideal) / ideal
        sign = "+" if pct >= 0 else ""
        lines.append(f"| {d} | {pop:,} | {sign}{pct:.1f}% |")
    worst = max(abs(100 * (p - ideal) / ideal) for p in pops.values())
    lines += ["", f"Worst deviation: {worst:.1f}%"]
    with open(path, "w") as fh:
        fh.write("\n".join(lines) + "\n")


def main() -> None:
    os.makedirs(OUTPUT_DIR, exist_ok=True)
    summary: list[tuple[str, int]] = []

    for geography, recursive, label in VARIANTS:
        print(f"\n{'=' * 60}")
        print(f"Variant: {label}")
        print("=" * 60)

        run_id = cli.run(
            NY_FIPS, geography, N_DISTRICTS,
            ncuts=NCUTS, niter=NITER,
            formula=FORMULA, recursive=recursive,
        )
        summary.append((label, run_id))

        conn = db.connect()
        try:
            geojson_path = os.path.join(OUTPUT_DIR, f"run_{run_id}_{label}.geojson")
            db.export_geojson(conn, run_id, geojson_path)
            print(f"  GeoJSON -> {geojson_path}")

            pops = db.fetch_district_populations(conn, run_id)
        finally:
            conn.close()

        log_path = os.path.join(OUTPUT_DIR, f"run_{run_id}_{label}_deviations.md")
        write_deviation_log(log_path, run_id, label, pops)
        print(f"  Deviations -> {log_path}")

    print("\n\nSummary — Run IDs:")
    for label, run_id in summary:
        print(f"  {label:30s}  run_id={run_id}")


if __name__ == "__main__":
    main()
