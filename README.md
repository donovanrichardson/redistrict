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

## Running

```bash
uv run redistrict
```

## Tests

```bash
uv run pytest
```

## Dependencies

| Package | Purpose |
|---|---|
| pymetis | Graph partitioning (METIS wrapper, builds from source — requires cmake) |
| scipy | Delaunay triangulation |
| networkx | Graph construction and contiguity checking |
| shapely / pyproj | Geometry and projection |
| textual | TUI framework (extensible to full GUI) |
| inquirerpy | Interactive CLI prompts (maintained replacement for abandoned PyInquirer) |
| psycopg2-binary | PostgreSQL connection |
