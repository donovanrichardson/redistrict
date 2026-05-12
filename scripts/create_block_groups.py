#!/usr/bin/env python3
"""
Create (or extend) block_groups_2020 by aggregating blocks_2020.

Each block group row contains:
  - geoid20      TEXT PK  (12-char: LEFT(geoid20, 12))
  - statefp20    TEXT
  - countyfp20   TEXT
  - tractce20    TEXT
  - geom         MultiPolygon(4269) via ST_Union of block geometries
  - pop20        INTEGER  sum of block pop20
  - aland20      NUMERIC  sum of block aland20
  - awater20     NUMERIC  sum of block awater20
  - intptlat20   DOUBLE PRECISION  population-weighted centroid lat
  - intptlon20   DOUBLE PRECISION  population-weighted centroid lon
  - centroid_geom  Point(4269)  built from intptlat20/intptlon20

Population-weighted centroid:
  lat = SUM(pop20 * intptlat20::float) / SUM(pop20)
  lon = SUM(pop20 * intptlon20::float) / SUM(pop20)
  Falls back to geometric centroid of ST_Union when total pop = 0.

The table is created if it does not exist. Processing is county-by-county
within each state (tqdm progress). Each state is committed as one transaction.
"""

import logging
import os
import sys
import time
from datetime import datetime
from pathlib import Path

import psycopg2
from tqdm import tqdm

# States to process: Rhode Island (44), Vermont (50), New Hampshire (33).
# Override with STATE_FIPS env var (comma-separated) or set to empty for all states.
_env_fips = os.getenv("STATE_FIPS", "50,33")
STATE_FIPS_LIST: list[str] = [s.strip() for s in _env_fips.split(",") if s.strip()] if _env_fips else []

DB_CRED = {
    "host":     os.getenv("POSTGRES_HOST", "localhost"),
    "port":     int(os.getenv("POSTGRES_HOST_PORT", "5433")),
    "dbname":   os.getenv("POSTGRES_DB",   "block-county"),
    "user":     os.getenv("POSTGRES_USER", "block-county"),
    "password": os.getenv("POSTGRES_PASSWORD"),
}

CREATE_TABLE_SQL = """
CREATE TABLE IF NOT EXISTS public.block_groups_2020 (
    geoid20       text                         NOT NULL,
    statefp20     text,
    countyfp20    text,
    tractce20     text,
    geom          geometry(MultiPolygon, 4269),
    pop20         integer,
    aland20       numeric,
    awater20      numeric,
    intptlat20    double precision,
    intptlon20    double precision,
    centroid_geom geometry(Point, 4269),
    PRIMARY KEY (geoid20)
);

CREATE INDEX IF NOT EXISTS block_groups_2020_geom_idx
    ON public.block_groups_2020 USING GIST (geom);
CREATE INDEX IF NOT EXISTS block_groups_2020_centroid_idx
    ON public.block_groups_2020 USING GIST (centroid_geom);
CREATE INDEX IF NOT EXISTS block_groups_2020_statefp_idx
    ON public.block_groups_2020 (statefp20);
"""

COUNTIES_SQL = """
SELECT DISTINCT countyfp20
FROM public.blocks_2020
WHERE statefp20 = %(statefp)s
ORDER BY countyfp20;
"""

UPSERT_COUNTY_SQL = """
INSERT INTO public.block_groups_2020
    (geoid20, statefp20, countyfp20, tractce20,
     geom, pop20, aland20, awater20,
     intptlat20, intptlon20, centroid_geom)
WITH agg AS (
    SELECT
        LEFT(geoid20, 12)                          AS geoid20,
        statefp20,
        countyfp20,
        tractce20,
        ST_Union(geom)                             AS geom,
        SUM(pop20)::integer                        AS pop20,
        SUM(aland20)                               AS aland20,
        SUM(awater20)                              AS awater20,
        SUM(pop20::float * intptlat20::float)      AS weighted_lat_sum,
        SUM(pop20::float * intptlon20::float)      AS weighted_lon_sum,
        SUM(pop20::float)                          AS pop_float,
        ST_Y(ST_Centroid(ST_Union(geom)))          AS geom_centroid_lat,
        ST_X(ST_Centroid(ST_Union(geom)))          AS geom_centroid_lon
    FROM public.blocks_2020
    WHERE geom IS NOT NULL
      AND statefp20  = %(statefp)s
      AND countyfp20 = %(countyfp)s
    GROUP BY LEFT(geoid20, 12), statefp20, countyfp20, tractce20
)
SELECT
    geoid20,
    statefp20,
    countyfp20,
    tractce20,
    ST_Multi(geom)::geometry(MultiPolygon, 4269),
    pop20,
    aland20,
    awater20,
    CASE WHEN pop_float > 0 THEN weighted_lat_sum / pop_float ELSE geom_centroid_lat END,
    CASE WHEN pop_float > 0 THEN weighted_lon_sum / pop_float ELSE geom_centroid_lon END,
    ST_SetSRID(ST_MakePoint(
        CASE WHEN pop_float > 0 THEN weighted_lon_sum / pop_float ELSE geom_centroid_lon END,
        CASE WHEN pop_float > 0 THEN weighted_lat_sum / pop_float ELSE geom_centroid_lat END
    ), 4269)::geometry(Point, 4269)
FROM agg
ON CONFLICT (geoid20) DO UPDATE SET
    statefp20     = EXCLUDED.statefp20,
    countyfp20    = EXCLUDED.countyfp20,
    tractce20     = EXCLUDED.tractce20,
    geom          = EXCLUDED.geom,
    pop20         = EXCLUDED.pop20,
    aland20       = EXCLUDED.aland20,
    awater20      = EXCLUDED.awater20,
    intptlat20    = EXCLUDED.intptlat20,
    intptlon20    = EXCLUDED.intptlon20,
    centroid_geom = EXCLUDED.centroid_geom;
"""


def setup_logging() -> logging.Logger:
    logs_dir = Path(__file__).parent.parent / "logs"
    logs_dir.mkdir(exist_ok=True)
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    log_path = logs_dir / f"create_block_groups_{timestamp}.log"

    logger = logging.getLogger("create_block_groups")
    logger.setLevel(logging.DEBUG)
    fmt = logging.Formatter("%(asctime)s %(levelname)s %(message)s")

    fh = logging.FileHandler(log_path)
    fh.setFormatter(fmt)
    logger.addHandler(fh)

    sh = logging.StreamHandler(sys.stdout)
    sh.setFormatter(fmt)
    logger.addHandler(sh)

    logger.info(f"Log file: {log_path}")
    return logger


def process_state(cur, statefp: str, log: logging.Logger) -> int:
    cur.execute(COUNTIES_SQL, {"statefp": statefp})
    counties = [row[0] for row in cur.fetchall()]
    if not counties:
        log.warning(f"State {statefp}: no counties found in blocks_2020, skipping.")
        return 0

    total_rows = 0
    for countyfp in tqdm(counties, desc=f"State {statefp}", unit="county", leave=True):
        cur.execute(UPSERT_COUNTY_SQL, {"statefp": statefp, "countyfp": countyfp})
        total_rows += cur.rowcount

    return total_rows


def main():
    log = setup_logging()
    scope = ", ".join(STATE_FIPS_LIST) if STATE_FIPS_LIST else "all states"
    log.info(f"Processing states: {scope}")

    conn = psycopg2.connect(**DB_CRED)
    conn.autocommit = False

    try:
        with conn.cursor() as cur:
            log.info("Ensuring block_groups_2020 table and indexes exist...")
            cur.execute(CREATE_TABLE_SQL)
            conn.commit()

        states_to_run = STATE_FIPS_LIST or _get_all_states(conn)

        for statefp in states_to_run:
            log.info(f"Starting state {statefp}...")
            t0 = time.time()
            with conn.cursor() as cur:
                rows = process_state(cur, statefp, log)
            conn.commit()
            elapsed = time.time() - t0
            log.info(f"State {statefp} done in {elapsed:.1f}s — {rows:,} rows upserted.")

        with conn.cursor() as cur:
            cur.execute("SELECT COUNT(*) FROM public.block_groups_2020")
            total = cur.fetchone()[0]
            cur.execute("SELECT COUNT(*) FROM public.block_groups_2020 WHERE pop20 = 0")
            zero_pop = cur.fetchone()[0]
            cur.execute("SELECT COUNT(DISTINCT statefp20) FROM public.block_groups_2020")
            states = cur.fetchone()[0]

        log.info(f"Table totals: rows={total:,} zero_pop={zero_pop:,} states={states}")

    except Exception as e:
        conn.rollback()
        log.exception(f"Fatal error: {e}")
        sys.exit(1)
    finally:
        conn.close()


def _get_all_states(conn) -> list[str]:
    with conn.cursor() as cur:
        cur.execute("SELECT DISTINCT statefp20 FROM public.blocks_2020 ORDER BY statefp20")
        return [row[0] for row in cur.fetchall()]


if __name__ == "__main__":
    main()
