# redistrict

Population-balanced redistricting using PyMETIS graph partitioning on US Census geographies.

Given a state and a geographic unit (counties, tracts, block groups, or blocks), the app builds a graph from census centroids, weights edges by population, and uses METIS to partition the graph into districts with roughly equal population and cuts that fall in low-population areas.

## Requirements

- Python 3.12+
- [uv](https://docs.astral.sh/uv/)
- Docker + Docker Compose
- cmake (for building PyMETIS): `brew install cmake`

## Setup

```bash
# Install dependencies
uv sync

# Copy and configure environment
cp .env.example .env
# Edit .env as needed (see Database section below)

# Start the database
docker compose up -d
```

## Database

By default, `docker compose up` starts a fresh local PostgreSQL instance on port 5432.

To use an existing database (e.g. a local hex-maps instance), set `DB_COMPOSE_PATH` in your `.env`:

```
DB_COMPOSE_PATH=../hex-maps/docker-compose.yml
POSTGRES_HOST_PORT=5433
POSTGRES_DB=block-county
POSTGRES_USER=block-county
POSTGRES_PASSWORD=your_db_password
```

## Loading geography tables

`run_tract_metis.py` reads from `blocks_2020`, which must already be loaded into the database. Block groups, tracts, and counties are derived tables built by aggregating blocks upward through the hierarchy. Run the scripts in this order after blocks are loaded:

```bash
# 1. Aggregate blocks → block groups
uv run --env-file .env python scripts/create_block_groups.py

# 2. Aggregate block groups → tracts
uv run --env-file .env python scripts/create_tracts.py

# 3. Aggregate tracts → counties
uv run --env-file .env python scripts/create_counties.py
```

Each script processes all states incrementally (skipping states already marked done in its state log table) and can be safely re-run. No TIGER imports are needed beyond blocks.

## Running

The below command will create 2 districts from the Census blocks of Rhode Island, which is identified by the FIPS code 44. Choose other states' FIPS and other district number values to customize your run.

```bash
uv run --env-file .env python scripts/run_tract_metis.py 44 2
```

## Tests

```bash
uv run pytest
```

## Dependencies

| Package | Purpose |
|---|---|
| pymetis | Graph partitioning (METIS wrapper, builds from source — requires cmake) |
| scipy | Spherical convex hull construction for adjacency and bridging |
| numpy | Coordinate projection and array operations |
| shapely | Geometry union, WKB loading, spatial indexing |
| psycopg2-binary | PostgreSQL connection |
| tqdm | Progress bars for bisection and bridging loops |
