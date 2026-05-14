#!/usr/bin/env python3
"""
Create (or extend) tracts_2020 by aggregating block_groups_2020.

Each tract row contains:
  - geoid20      TEXT PK  (11-char: LEFT(geoid20, 11))
  - statefp20    TEXT
  - countyfp20   TEXT
  - tractce20    TEXT
  - geom         MultiPolygon(4269) via ST_Union of block group geometries
  - pop20        INTEGER  sum of block group pop20
  - aland20      NUMERIC  sum of block group aland20
  - awater20     NUMERIC  sum of block group awater20
  - intptlat20   DOUBLE PRECISION  population-weighted centroid lat
  - intptlon20   DOUBLE PRECISION  population-weighted centroid lon
  - centroid_geom  Point(4269)  built from intptlat20/intptlon20

Population-weighted centroid:
  lat = SUM(pop20 * intptlat20) / SUM(pop20)
  lon = SUM(pop20 * intptlon20) / SUM(pop20)
  Falls back to geometric centroid of ST_Union when total pop = 0.

Goal
----
Tracts are the primary geography for the redistricting algorithm. They are
derived by aggregating block_groups_2020 (which were themselves derived from
blocks_2020), preserving population-weighted centroids up the hierarchy so
that each tract centroid reflects where people actually live.

Processing
----------
By default, processes all states that are 'done' in block_groups_state_log
but not yet recorded in tracts_state_log, in FIPS order. Override with
STATE_FIPS env var (comma-separated) to target specific states. Each state
is one transaction. On failure the data is rolled back, a 'failed' row is
written to tracts_state_log, and the script moves on — that state will not
be retried on subsequent runs unless its log row is deleted or STATE_FIPS
override is used.
"""

import logging
import os
import sys
import time
from datetime import datetime
from pathlib import Path

import psycopg2
from tqdm import tqdm

_env_fips = os.getenv("STATE_FIPS", "")
STATE_FIPS_OVERRIDE: list[str] = [s.strip() for s in _env_fips.split(",") if s.strip()]

DB_CRED = {
    "host":     os.getenv("POSTGRES_HOST", "localhost"),
    "port":     int(os.getenv("POSTGRES_HOST_PORT", "5433")),
    "dbname":   os.getenv("POSTGRES_DB",   "block-county"),
    "user":     os.getenv("POSTGRES_USER", "block-county"),
    "password": os.getenv("POSTGRES_PASSWORD"),
}

CREATE_TABLE_SQL = """
CREATE TABLE IF NOT EXISTS public.tracts_state_log (
    statefp20   text        NOT NULL PRIMARY KEY,
    status      text        NOT NULL CHECK (status IN ('done', 'failed')),
    updated_at  timestamptz NOT NULL DEFAULT now(),
    error_msg   text
);

CREATE TABLE IF NOT EXISTS public.tracts_2020 (
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

CREATE INDEX IF NOT EXISTS tracts_2020_geom_idx
    ON public.tracts_2020 USING GIST (geom);
CREATE INDEX IF NOT EXISTS tracts_2020_centroid_idx
    ON public.tracts_2020 USING GIST (centroid_geom);
CREATE INDEX IF NOT EXISTS tracts_2020_statefp_idx
    ON public.tracts_2020 (statefp20);
"""

# States that are done in block_groups_state_log but not yet in tracts_state_log.
PENDING_STATES_SQL = """
SELECT l.statefp20
FROM public.block_groups_state_log l
WHERE l.status = 'done'
  AND NOT EXISTS (
      SELECT 1 FROM public.tracts_state_log tl
      WHERE tl.statefp20 = l.statefp20
  )
ORDER BY l.statefp20;
"""

LOG_STATE_SQL = """
INSERT INTO public.tracts_state_log (statefp20, status, updated_at, error_msg)
VALUES (%(statefp)s, %(status)s, now(), %(error_msg)s)
ON CONFLICT (statefp20) DO UPDATE SET
    status     = EXCLUDED.status,
    updated_at = EXCLUDED.updated_at,
    error_msg  = EXCLUDED.error_msg;
"""

COUNTIES_SQL = """
SELECT DISTINCT countyfp20
FROM public.block_groups_2020
WHERE statefp20 = %(statefp)s
ORDER BY countyfp20;
"""

UPSERT_COUNTY_SQL = """
INSERT INTO public.tracts_2020
    (geoid20, statefp20, countyfp20, tractce20,
     geom, pop20, aland20, awater20,
     intptlat20, intptlon20, centroid_geom)
WITH agg AS (
    SELECT
        LEFT(geoid20, 11)                          AS geoid20,
        statefp20,
        countyfp20,
        tractce20,
        ST_Union(geom)                             AS geom,
        SUM(pop20)::integer                        AS pop20,
        SUM(aland20)                               AS aland20,
        SUM(awater20)                              AS awater20,
        SUM(pop20::float * intptlat20)             AS weighted_lat_sum,
        SUM(pop20::float * intptlon20)             AS weighted_lon_sum,
        SUM(pop20::float)                          AS pop_float,
        ST_Y(ST_Centroid(ST_Union(geom)))          AS geom_centroid_lat,
        ST_X(ST_Centroid(ST_Union(geom)))          AS geom_centroid_lon
    FROM public.block_groups_2020
    WHERE geom IS NOT NULL
      AND statefp20  = %(statefp)s
      AND countyfp20 = %(countyfp)s
    GROUP BY LEFT(geoid20, 11), statefp20, countyfp20, tractce20
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
    log_path = logs_dir / f"create_tracts_{timestamp}.log"

    logger = logging.getLogger("create_tracts")
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


def get_pending_states(conn) -> list[str]:
    """Return states done in block_groups_state_log but not yet in tracts_state_log."""
    with conn.cursor() as cur:
        cur.execute(PENDING_STATES_SQL)
        return [row[0] for row in cur.fetchall()]


def process_state(cur, statefp: str, log: logging.Logger) -> int:
    """Upsert all tracts for one state, county by county. Returns row count."""
    cur.execute(COUNTIES_SQL, {"statefp": statefp})
    counties = [row[0] for row in cur.fetchall()]
    if not counties:
        log.warning(f"State {statefp}: no counties found in block_groups_2020, skipping.")
        return 0

    total_rows = 0
    for countyfp in tqdm(counties, desc=f"State {statefp}", unit="county", leave=True):
        cur.execute(UPSERT_COUNTY_SQL, {"statefp": statefp, "countyfp": countyfp})
        total_rows += cur.rowcount
    return total_rows


def main():
    log = setup_logging()

    conn = psycopg2.connect(**DB_CRED)
    conn.autocommit = False

    try:
        with conn.cursor() as cur:
            log.info("Ensuring tracts_2020 table and indexes exist...")
            cur.execute(CREATE_TABLE_SQL)
            conn.commit()

        if STATE_FIPS_OVERRIDE:
            states_to_run = STATE_FIPS_OVERRIDE
            log.info(f"Processing specified states: {', '.join(states_to_run)}")
        else:
            states_to_run = get_pending_states(conn)
            log.info(f"Found {len(states_to_run)} pending state(s): {', '.join(states_to_run)}")

        if not states_to_run:
            log.info("Nothing to do — all eligible states already present in tracts_2020.")
            return

        failed: list[str] = []

        for statefp in states_to_run:
            log.info(f"Starting state {statefp}...")
            t0 = time.time()
            try:
                with conn.cursor() as cur:
                    rows = process_state(cur, statefp, log)
                    cur.execute(LOG_STATE_SQL, {"statefp": statefp, "status": "done", "error_msg": None})
                # ── state-level transaction commit ──
                conn.commit()
                elapsed = time.time() - t0
                log.info(f"State {statefp} done in {elapsed:.1f}s — {rows:,} rows upserted.")
            except Exception as e:
                conn.rollback()
                elapsed = time.time() - t0
                log.error(f"State {statefp} failed after {elapsed:.1f}s — {e}. Skipping.")
                failed.append(statefp)
                with conn.cursor() as cur:
                    cur.execute(LOG_STATE_SQL, {"statefp": statefp, "status": "failed", "error_msg": str(e)})
                conn.commit()

        with conn.cursor() as cur:
            cur.execute("SELECT COUNT(*) FROM public.tracts_2020")
            total = cur.fetchone()[0]
            cur.execute("SELECT COUNT(*) FROM public.tracts_2020 WHERE pop20 = 0")
            zero_pop = cur.fetchone()[0]
            cur.execute("SELECT COUNT(DISTINCT statefp20) FROM public.tracts_2020")
            states_done = cur.fetchone()[0]

        log.info(f"Table totals: rows={total:,} zero_pop={zero_pop:,} states={states_done}")
        if failed:
            log.warning(f"Failed states (not retried): {', '.join(failed)}")

    finally:
        conn.close()


if __name__ == "__main__":
    main()
