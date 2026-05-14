"""Database access: connection, adjacency management, results persistence."""

import json
import os

import psycopg2
import psycopg2.extensions
from psycopg2.extras import execute_values

DB_CRED = {
    "host":     os.getenv("POSTGRES_HOST", "localhost"),
    "port":     int(os.getenv("POSTGRES_HOST_PORT", "5433")),
    "dbname":   os.getenv("POSTGRES_DB",   "block-county"),
    "user":     os.getenv("POSTGRES_USER", "block-county"),
    "password": os.getenv("POSTGRES_PASSWORD"),
}

_GEOGRAPHY_TABLES = {"tracts", "counties", "block_groups", "blocks"}

_ADJACENCY_DDL = """
CREATE TABLE IF NOT EXISTS public.tract_adjacency (
    geoid_a  text NOT NULL,
    geoid_b  text NOT NULL,
    PRIMARY KEY (geoid_a, geoid_b),
    CHECK (geoid_a < geoid_b)
);
CREATE INDEX IF NOT EXISTS tract_adjacency_a_idx ON public.tract_adjacency (geoid_a);
CREATE INDEX IF NOT EXISTS tract_adjacency_b_idx ON public.tract_adjacency (geoid_b);

CREATE TABLE IF NOT EXISTS public.tract_adjacency_log (
    geoid20       text        NOT NULL PRIMARY KEY,
    calculated_at timestamptz NOT NULL DEFAULT now()
);
"""

_RESULTS_DDL = """
CREATE TABLE IF NOT EXISTS public.redistrict_runs (
    id          serial      PRIMARY KEY,
    geography   text        NOT NULL,
    statefp20   text        NOT NULL,
    n_districts integer     NOT NULL,
    created_at  timestamptz NOT NULL DEFAULT now(),
    params      jsonb
);

CREATE TABLE IF NOT EXISTS public.redistrict_districts (
    run_id      integer NOT NULL REFERENCES public.redistrict_runs(id),
    district_id integer NOT NULL,
    geom        geometry(MultiPolygon, 4269),
    pop20       integer,
    PRIMARY KEY (run_id, district_id)
);

CREATE TABLE IF NOT EXISTS public.redistrict_assignments (
    run_id      integer NOT NULL REFERENCES public.redistrict_runs(id),
    geoid20     text    NOT NULL,
    district_id integer NOT NULL,
    PRIMARY KEY (run_id, geoid20)
);
"""


def connect() -> psycopg2.extensions.connection:
    return psycopg2.connect(**DB_CRED)


def ensure_tables(conn: psycopg2.extensions.connection) -> None:
    """Create adjacency and results tables if they do not exist."""
    with conn.cursor() as cur:
        cur.execute(_ADJACENCY_DDL)
        cur.execute(_RESULTS_DDL)
    conn.commit()


def fetch_tracts(conn: psycopg2.extensions.connection, statefp: str) -> list[dict]:
    """Return all tracts for a state as dicts with geoid, pop, lat, lon."""
    with conn.cursor() as cur:
        cur.execute(
            """
            SELECT geoid20, pop20, intptlat20, intptlon20
            FROM public.tracts_2020
            WHERE statefp20 = %s
            ORDER BY geoid20
            """,
            (statefp,),
        )
        rows = cur.fetchall()
    return [{"geoid": r[0], "pop": r[1] or 0, "lat": r[2], "lon": r[3]} for r in rows]


def fetch_available_states(conn: psycopg2.extensions.connection) -> list[tuple[str, str]]:
    """Return (statefp20, state_name_placeholder) for all states in tracts_2020."""
    with conn.cursor() as cur:
        cur.execute(
            "SELECT DISTINCT statefp20 FROM public.tracts_2020 ORDER BY statefp20"
        )
        return [(r[0], r[0]) for r in cur.fetchall()]


def get_missing_adjacency_geoids(
    conn: psycopg2.extensions.connection, geoids: list[str]
) -> list[str]:
    """Return the subset of geoids not yet in tract_adjacency_log."""
    if not geoids:
        return []
    with conn.cursor() as cur:
        cur.execute(
            """
            SELECT g FROM unnest(%s::text[]) AS g
            WHERE NOT EXISTS (
                SELECT 1 FROM public.tract_adjacency_log WHERE geoid20 = g
            )
            ORDER BY g
            """,
            (geoids,),
        )
        return [row[0] for row in cur.fetchall()]


def compute_and_store_adjacency(
    conn: psycopg2.extensions.connection,
    geoids: list[str],
    progress_callback=None,
) -> int:
    """
    For each geoid, find all rook-contiguous (shared-edge) tracts in the same
    state, insert new pairs into tract_adjacency, and record the geoid in
    tract_adjacency_log. Commits after each geoid so progress survives
    interruption.

    Rook contiguity: ST_Touches AND shared boundary has dimension >= 1 (line).

    Returns total new adjacency pairs inserted.
    """
    if not geoids:
        return 0

    total = 0
    for i, geoid in enumerate(geoids):
        with conn.cursor() as cur:
            cur.execute(
                """
                INSERT INTO public.tract_adjacency (geoid_a, geoid_b)
                SELECT
                    LEAST(a.geoid20, b.geoid20),
                    GREATEST(a.geoid20, b.geoid20)
                FROM public.tracts_2020 a
                JOIN public.tracts_2020 b
                  ON a.statefp20 = b.statefp20
                 AND a.geoid20 <> b.geoid20
                 AND ST_Touches(a.geom, b.geom)
                 AND ST_Dimension(ST_Intersection(a.geom, b.geom)) >= 1
                WHERE a.geoid20 = %s
                ON CONFLICT (geoid_a, geoid_b) DO NOTHING
                """,
                (geoid,),
            )
            total += cur.rowcount

            cur.execute(
                """
                INSERT INTO public.tract_adjacency_log (geoid20)
                VALUES (%s)
                ON CONFLICT (geoid20) DO NOTHING
                """,
                (geoid,),
            )
        conn.commit()
        if progress_callback:
            progress_callback(i + 1, len(geoids), geoid)

    return total


def fetch_adjacency(
    conn: psycopg2.extensions.connection, geoids: list[str]
) -> set[tuple[str, str]]:
    """Return all adjacent pairs where both geoids are in the provided list."""
    if not geoids:
        return set()
    with conn.cursor() as cur:
        cur.execute(
            """
            SELECT geoid_a, geoid_b
            FROM public.tract_adjacency
            WHERE geoid_a = ANY(%s) AND geoid_b = ANY(%s)
            """,
            (geoids, geoids),
        )
        return {(r[0], r[1]) for r in cur.fetchall()}


def write_run(
    conn: psycopg2.extensions.connection,
    geography: str,
    statefp: str,
    n_districts: int,
    params: dict,
) -> int:
    """Insert a run record and return its auto-generated id."""
    with conn.cursor() as cur:
        cur.execute(
            """
            INSERT INTO public.redistrict_runs (geography, statefp20, n_districts, params)
            VALUES (%s, %s, %s, %s)
            RETURNING id
            """,
            (geography, statefp, n_districts, json.dumps(params)),
        )
        run_id = cur.fetchone()[0]
    conn.commit()
    return run_id


def write_assignments(
    conn: psycopg2.extensions.connection,
    run_id: int,
    geoid_to_district: dict[str, int],
) -> None:
    """Bulk-insert geoid -> district_id assignments for a run."""
    rows = [(run_id, geoid, dist) for geoid, dist in geoid_to_district.items()]
    with conn.cursor() as cur:
        execute_values(
            cur,
            """
            INSERT INTO public.redistrict_assignments (run_id, geoid20, district_id)
            VALUES %s
            """,
            rows,
        )
    conn.commit()


def write_district_geoms(
    conn: psycopg2.extensions.connection,
    run_id: int,
    geography: str,
    geoid_to_district: dict[str, int],
) -> None:
    """
    Aggregate member geographies into district polygons via ST_Union and write
    them to redistrict_districts.
    """
    if geography not in _GEOGRAPHY_TABLES:
        raise ValueError(f"Unknown geography '{geography}'. Must be one of {_GEOGRAPHY_TABLES}.")

    table = f"public.{geography}_2020"
    districts: dict[int, list[str]] = {}
    for geoid, dist_id in geoid_to_district.items():
        districts.setdefault(dist_id, []).append(geoid)

    with conn.cursor() as cur:
        for dist_id, geoids in sorted(districts.items()):
            cur.execute(
                f"""
                INSERT INTO public.redistrict_districts (run_id, district_id, geom, pop20)
                SELECT %s, %s,
                    ST_Multi(ST_Union(geom))::geometry(MultiPolygon, 4269),
                    SUM(pop20)
                FROM {table}
                WHERE geoid20 = ANY(%s)
                """,
                (run_id, dist_id, geoids),
            )
    conn.commit()
