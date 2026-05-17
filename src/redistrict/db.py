"""Database access: connection, adjacency management, results persistence."""

import json
import math
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

CREATE TABLE IF NOT EXISTS public.redistrict_edges (
    id          serial  PRIMARY KEY,
    run_id      integer NOT NULL REFERENCES public.redistrict_runs(id),
    geoid_a     text    NOT NULL,
    geoid_b     text    NOT NULL,
    is_adjacent boolean NOT NULL,
    dist_km     float,
    geom        geometry(LineString, 4326),
    UNIQUE (run_id, geoid_a, geoid_b),
    CHECK (geoid_a < geoid_b)
);
CREATE INDEX IF NOT EXISTS redistrict_edges_run_idx
    ON public.redistrict_edges (run_id);
CREATE INDEX IF NOT EXISTS redistrict_edges_nonadj_idx
    ON public.redistrict_edges (run_id) WHERE NOT is_adjacent;
"""


def _validate_geography(geography: str) -> None:
    if geography not in _GEOGRAPHY_TABLES:
        raise ValueError(
            f"Unknown geography '{geography}'. Must be one of {_GEOGRAPHY_TABLES}."
        )


def _adjacency_ddl(geography: str) -> str:
    """Generate DDL for the adjacency and adjacency log tables for a geography."""
    return f"""
CREATE TABLE IF NOT EXISTS public.{geography}_adjacency (
    geoid_a  text NOT NULL,
    geoid_b  text NOT NULL,
    PRIMARY KEY (geoid_a, geoid_b),
    CHECK (geoid_a < geoid_b)
);
CREATE INDEX IF NOT EXISTS {geography}_adjacency_a_idx
    ON public.{geography}_adjacency (geoid_a);
CREATE INDEX IF NOT EXISTS {geography}_adjacency_b_idx
    ON public.{geography}_adjacency (geoid_b);

CREATE TABLE IF NOT EXISTS public.{geography}_adjacency_log (
    geoid20       text        NOT NULL PRIMARY KEY,
    calculated_at timestamptz NOT NULL DEFAULT now()
);
"""


def connect() -> psycopg2.extensions.connection:
    return psycopg2.connect(**DB_CRED)


def ensure_tables(
    conn: psycopg2.extensions.connection,
    geography: str = "tracts",
) -> None:
    """Create adjacency tables for the given geography and results tables."""
    _validate_geography(geography)
    with conn.cursor() as cur:
        cur.execute(_adjacency_ddl(geography))
        cur.execute(_RESULTS_DDL)
    conn.commit()


def fetch_nodes(
    conn: psycopg2.extensions.connection,
    geography: str,
    statefp: str,
) -> list[dict]:
    """Return all geographies for a state as dicts with geoid, pop, lat, lon."""
    _validate_geography(geography)
    table = f"public.{geography}_2020"
    with conn.cursor() as cur:
        cur.execute(
            f"""
            SELECT geoid20, pop20, intptlat20, intptlon20
            FROM {table}
            WHERE statefp20 = %s
            ORDER BY geoid20
            """,
            (statefp,),
        )
        rows = cur.fetchall()
    return [{"geoid": r[0], "pop": r[1] or 0, "lat": r[2], "lon": r[3]} for r in rows]


def fetch_available_states(
    conn: psycopg2.extensions.connection,
    geography: str = "tracts",
) -> list[str]:
    """Return sorted list of statefp20 values present in the geography table."""
    _validate_geography(geography)
    table = f"public.{geography}_2020"
    with conn.cursor() as cur:
        cur.execute(f"SELECT DISTINCT statefp20 FROM {table} ORDER BY statefp20")
        return [r[0] for r in cur.fetchall()]


def get_missing_adjacency_geoids(
    conn: psycopg2.extensions.connection,
    geography: str,
    geoids: list[str],
) -> list[str]:
    """Return the subset of geoids not yet in the geography's adjacency log."""
    if not geoids:
        return []
    _validate_geography(geography)
    log_table = f"public.{geography}_adjacency_log"
    with conn.cursor() as cur:
        cur.execute(
            f"""
            SELECT g FROM unnest(%s::text[]) AS g
            WHERE NOT EXISTS (
                SELECT 1 FROM {log_table} WHERE geoid20 = g
            )
            ORDER BY g
            """,
            (geoids,),
        )
        return [row[0] for row in cur.fetchall()]


def compute_and_store_adjacency(
    conn: psycopg2.extensions.connection,
    geography: str,
    geoids: list[str],
    progress_callback=None,
) -> int:
    """
    For each geoid, find all rook-contiguous (shared-edge) neighbors in the
    same state, insert new pairs into {geography}_adjacency, and record the
    geoid in {geography}_adjacency_log. Commits after each geoid so progress
    survives interruption.

    Rook contiguity: ST_Touches AND shared boundary has dimension >= 1 (line).

    Returns total new adjacency pairs inserted.
    """
    if not geoids:
        return 0
    _validate_geography(geography)
    geo_table = f"public.{geography}_2020"
    adj_table = f"public.{geography}_adjacency"
    log_table = f"public.{geography}_adjacency_log"

    total = 0
    for i, geoid in enumerate(geoids):
        with conn.cursor() as cur:
            cur.execute(
                f"""
                INSERT INTO {adj_table} (geoid_a, geoid_b)
                SELECT
                    LEAST(a.geoid20, b.geoid20),
                    GREATEST(a.geoid20, b.geoid20)
                FROM {geo_table} a
                JOIN {geo_table} b
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
                f"""
                INSERT INTO {log_table} (geoid20)
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
    conn: psycopg2.extensions.connection,
    geography: str,
    geoids: list[str],
) -> set[tuple[str, str]]:
    """Return all adjacent pairs where both geoids are in the provided list."""
    if not geoids:
        return set()
    _validate_geography(geography)
    adj_table = f"public.{geography}_adjacency"
    with conn.cursor() as cur:
        cur.execute(
            f"""
            SELECT geoid_a, geoid_b
            FROM {adj_table}
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


def fetch_run(
    conn: psycopg2.extensions.connection,
    run_id: int,
) -> dict:
    """Return the run record for run_id as a dict, or raise ValueError."""
    with conn.cursor() as cur:
        cur.execute(
            """
            SELECT geography, statefp20, n_districts, params
            FROM public.redistrict_runs
            WHERE id = %s
            """,
            (run_id,),
        )
        row = cur.fetchone()
    if row is None:
        raise ValueError(f"No run found with id={run_id}")
    return {
        "geography":   row[0],
        "statefp":     row[1],
        "n_districts": row[2],
        "params":      row[3] or {},
    }


def write_edges(
    conn: psycopg2.extensions.connection,
    run_id: int,
    nodes: list[dict],
    edges: set[tuple[int, int]],
    adjacent_geoid_pairs: set[tuple[str, str]],
) -> int:
    """
    Insert all Urquhart edges for a run with an is_adjacent flag.

    Each edge becomes a LineString between the two node centroids (EPSG:4326).
    Returns the number of rows inserted.
    """
    with conn.cursor() as cur:
        for i, j in edges:
            ga = min(nodes[i]["geoid"], nodes[j]["geoid"])
            gb = max(nodes[i]["geoid"], nodes[j]["geoid"])
            is_adj = (ga, gb) in adjacent_geoid_pairs
            loni, lati = nodes[i]["lon"], nodes[i]["lat"]
            lonj, latj = nodes[j]["lon"], nodes[j]["lat"]
            dlat = math.radians(latj - lati)
            dlon = math.radians(lonj - loni)
            a = (math.sin(dlat / 2) ** 2
                 + math.cos(math.radians(lati)) * math.cos(math.radians(latj))
                 * math.sin(dlon / 2) ** 2)
            dist_km = 2 * 6371.0 * math.asin(math.sqrt(a))
            cur.execute(
                """
                INSERT INTO public.redistrict_edges
                    (run_id, geoid_a, geoid_b, is_adjacent, dist_km, geom)
                VALUES (%s, %s, %s, %s, %s,
                    ST_SetSRID(ST_MakeLine(
                        ST_Point(%s, %s),
                        ST_Point(%s, %s)
                    ), 4326))
                ON CONFLICT (run_id, geoid_a, geoid_b) DO NOTHING
                """,
                (run_id, ga, gb, is_adj, dist_km, loni, lati, lonj, latj),
            )
    conn.commit()
    return len(edges)


def fetch_edges(
    conn: psycopg2.extensions.connection,
    run_id: int,
) -> tuple[set[tuple[str, str]], set[tuple[str, str]]]:
    """
    Return all remaining edges for run_id as two sets of geoid pairs:
      (adjacent_pairs, non_adjacent_pairs)

    Non-adjacent edges may have been deleted by the user in QGIS.
    """
    with conn.cursor() as cur:
        cur.execute(
            """
            SELECT geoid_a, geoid_b, is_adjacent
            FROM public.redistrict_edges
            WHERE run_id = %s
            """,
            (run_id,),
        )
        adjacent: set[tuple[str, str]] = set()
        non_adjacent: set[tuple[str, str]] = set()
        for geoid_a, geoid_b, is_adj in cur.fetchall():
            (adjacent if is_adj else non_adjacent).add((geoid_a, geoid_b))
    return adjacent, non_adjacent


def fetch_district_populations(
    conn: psycopg2.extensions.connection,
    run_id: int,
) -> dict[int, int]:
    """Return {district_id: pop20} for a completed run."""
    with conn.cursor() as cur:
        cur.execute(
            """
            SELECT district_id, pop20
            FROM public.redistrict_districts
            WHERE run_id = %s
            ORDER BY district_id
            """,
            (run_id,),
        )
        return {row[0]: row[1] for row in cur.fetchall()}


def export_geojson(
    conn: psycopg2.extensions.connection,
    run_id: int,
    output_path: str,
) -> None:
    """Write districts for run_id to a GeoJSON FeatureCollection file."""
    with conn.cursor() as cur:
        cur.execute(
            """
            SELECT json_build_object(
                'type', 'FeatureCollection',
                'features', json_agg(
                    json_build_object(
                        'type', 'Feature',
                        'geometry', ST_AsGeoJSON(ST_Transform(geom, 4326))::json,
                        'properties', json_build_object(
                            'run_id',      run_id,
                            'district_id', district_id,
                            'pop20',       pop20
                        )
                    )
                )
            )
            FROM public.redistrict_districts
            WHERE run_id = %s
            """,
            (run_id,),
        )
        result = cur.fetchone()[0]
    with open(output_path, "w") as fh:
        json.dump(result, fh)


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
    _validate_geography(geography)
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
