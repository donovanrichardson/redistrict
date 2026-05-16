#!/usr/bin/env python3
"""
Repeat run 33 (block_groups, k-way, uniform) with escalating water penalties:
  2x (10.0), 4x (20.0), 8x (40.0) the baseline of 5.0.

Exports GeoJSONs and deviation logs to output/ny_experiments/.
"""
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))

from redistrict import cli, db

NY_FIPS = "36"
GEOGRAPHY = "block_groups"
N_DISTRICTS = 26
NCUTS = 10
NITER = 20
FORMULA = "uniform"
RECURSIVE = False
BASE_PENALTY = 5.0

PENALTIES = [
    (BASE_PENALTY * 2, "water2x"),
    (BASE_PENALTY * 4, "water4x"),
    (BASE_PENALTY * 8, "water8x"),
]

OUTPUT_DIR = os.path.join(
    os.path.dirname(__file__), "..", "output", "ny_experiments"
)


def write_deviation_log(path: str, run_id: int, label: str, penalty: float, pops: dict[int, int]) -> None:
    total = sum(pops.values())
    ideal = total / N_DISTRICTS
    lines = [
        f"# NY block_groups k-way uniform — {label} — Run {run_id}",
        f"",
        f"water_penalty={penalty}  ncuts={NCUTS}  niter={NITER}  formula={FORMULA}",
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
    summary: list[tuple[str, float, int]] = []

    for penalty, label in PENALTIES:
        print(f"\n{'=' * 60}")
        print(f"water_penalty={penalty}  ({label})")
        print("=" * 60)

        run_id = cli.run(
            NY_FIPS, GEOGRAPHY, N_DISTRICTS,
            ncuts=NCUTS, niter=NITER,
            formula=FORMULA, recursive=RECURSIVE,
            water_penalty=penalty,
        )
        summary.append((label, penalty, run_id))

        conn = db.connect()
        try:
            geojson_path = os.path.join(OUTPUT_DIR, f"run_{run_id}_{label}.geojson")
            db.export_geojson(conn, run_id, geojson_path)
            print(f"  GeoJSON -> {geojson_path}")

            pops = db.fetch_district_populations(conn, run_id)
        finally:
            conn.close()

        log_path = os.path.join(OUTPUT_DIR, f"run_{run_id}_{label}_deviations.md")
        write_deviation_log(log_path, run_id, label, penalty, pops)
        print(f"  Deviations -> {log_path}")

    print("\n\nSummary — Run IDs:")
    for label, penalty, run_id in summary:
        print(f"  {label} (penalty={penalty:.1f})  run_id={run_id}")


if __name__ == "__main__":
    main()
