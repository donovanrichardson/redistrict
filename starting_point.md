# Redistricting Project: Starting Point

## Context

- Working directory: `codebases2026/redistrict`
- Related prior work: `codebases2025` contains several districting projects using k-medoids algorithms
- One of those projects uses a Python package called **PyMETIS** (METIS graph partitioning library)

---

## Database

- PostgreSQL instance running in a Docker container named **HexMaps**
- Currently not running due to a port conflict with another PostgreSQL Docker container (same port mapping)
- The database is believed to contain:
  - US Census geographic units: counties, tracts, block groups, blocks
  - Population data for all geographic units
  - Possibly pre-calculated adjacency data (uncertain)

### Immediate Action Needed

Fix the Docker port conflict so HexMaps can run, then introspect the database to confirm:
- What tables and schemas are present
- Whether population data is loaded
- Whether adjacency or graph data was previously calculated

---

## Algorithm Design

### Graph Construction

1. Take all census geographies within a given boundary (e.g., a state)
2. Compute a **Delaunay triangulation** (or convex hull) over the centroids of those geographies
3. This produces a graph where nodes are census units and edges connect nearby centroids

### Edge Weights (for k-medoids / exploration)

Edges get a weight based on the populations of the two connected nodes:

```
w(a, b) = 0.5 * (1 / sqrt(pop_a)) + 0.5 * (1 / sqrt(pop_b))
```

Effect: edges connecting two highly populated areas get low weight (easy to traverse); edges connecting sparsely populated or distant areas get high weight.

### Edge Weights (for PyMETIS partitioning)

PyMETIS performs a **minimum edge-weight cut**. To get district boundaries that fall in low-population areas, the weights must be **inverted** relative to the k-medoids formula above:

```
w_metis(a, b) = sqrt(pop_a) * sqrt(pop_b)   (approximately)
```

Or more precisely, the inverse of the k-medoids weight, so that:
- Edges between dense urban areas have **high** weight (METIS avoids cutting them)
- Edges between sparse rural areas have **low** weight (METIS prefers to cut here)

### Node Weights

Each node carries a weight equal to its population. PyMETIS will try to keep node-weight totals roughly equal across all partitions, producing **population-balanced districts**.

---

## Application Plan

### Stack

- **CLI interface**: Textual (extensible to full GUI) + PyInquirer
- **Graph construction**: Delaunay triangulation via SciPy or similar
- **Projection**: Either a planar projection or globe-based (globe-based was used in HexMaps previously)
- **Partitioning**: PyMETIS

### User Flow

1. User is prompted to select a **state** (or other boundary)
2. User selects a **geographic unit**: counties, tracts, block groups, or blocks
3. User specifies **number of districts** (partitions)
4. App fetches geometries and populations from the database
5. App constructs the Delaunay triangulation graph
6. App computes edge and node weights
7. App runs PyMETIS partitioning
8. Results are stored to the database (see below)

### First Test Run

- **State**: Rhode Island (smallest state, lowest complexity)
- **Geography**: Census tracts
- **Districts**: 2

---

## Database Output Schema

Each algorithm run should persist:

1. **Run record**: metadata (state, geography type, number of districts, timestamp, parameters)
2. **District geometries**: one geographic feature (polygon/multipolygon) per district, associated with the run
3. **Assignment table**: join table mapping each census unit to its district within a given run

---

## Development Process

1. Fix Docker port conflict, bring HexMaps online
2. Introspect the database schema and confirm available data
3. Build the CLI app using **test-driven development**
4. Write user-facing documentation (setup, usage, how to run for a new state)
5. Execute a real run on Rhode Island (2 districts)
6. Draft a blog post about the process and the Rhode Island result

---

## Blog Post (later)

- Audience: people interested in computational redistricting, graph algorithms, or census data
- Content: the full process, the algorithm choices, and the Rhode Island result
- Style: the user will rewrite AI-drafted content; the AI's role is to provide an outline and a first draft to react to
