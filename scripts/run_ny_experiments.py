#!/usr/bin/env python3
"""
Run four edge-weight formula experiments on New York census tracts.
State: New York (FIPS 36), 26 districts.

Formulas:
  1. original          - k-medoids cost formula
  2. uniform           - all land edges equal weight
  3. original_clamped  - original + constant so w_max/w_min == 4
  4. blend             - 50/50 normalised sum of original and 1/haversine

Prints run IDs for each formula and exports district GeoJSONs to
output/ny_experiments/.
"""
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))

from redistrict import cli, db

NY_FIPS = "36"
GEOGRAPHY = "tracts"
N_DISTRICTS = 26

FORMULAS = [
    "original",
    "uniform",
    "original_clamped",
    "blend",
]

OUTPUT_DIR = os.path.join(
    os.path.dirname(__file__), "..", "output", "ny_experiments"
)


def main() -> None:
    os.makedirs(OUTPUT_DIR, exist_ok=True)
    run_ids: dict[str, int] = {}

    for formula in FORMULAS:
        print(f"\n{'=' * 60}")
        print(f"Formula: {formula}")
        print("=" * 60)
        run_id = cli.run(NY_FIPS, GEOGRAPHY, N_DISTRICTS, formula=formula)
        run_ids[formula] = run_id

        conn = db.connect()
        try:
            path = os.path.join(OUTPUT_DIR, f"run_{run_id}_{formula}.geojson")
            db.export_geojson(conn, run_id, path)
            print(f"  GeoJSON -> {path}")
        finally:
            conn.close()

    print("\n\nSummary — Run IDs:")
    for formula, run_id in run_ids.items():
        print(f"  {formula:20s}  run_id={run_id}")


if __name__ == "__main__":
    main()
