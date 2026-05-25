# Changelog

## [Unreleased] — feature/two-pass-metis

### Added

**Tract-based METIS pipeline** (`scripts/run_tract_metis.py`)
- New redistricting pipeline that uses census tracts as the natural subcluster unit instead of arbitrary graph partitions.
- Populated blocks are grouped by tract (first 11 chars of GEOID20). Each tract's block geometries are unioned and exploded into singlepart polygons, so geographically disconnected tract fragments become separate subclusters without any rook-adjacency graph traversal.
- Threshold for bisection is 2× the median singlepart population (robust to large outlier tracts).
- Iterative METIS bisection (k=2, haversine-inverse edge weights) splits oversized subclusters. Subclusters that cannot be reduced below threshold are tracked in a `give_up` set and skipped in subsequent iterations to prevent infinite loops. tqdm progress bar on each bisection iteration.
- Zero-pop blocks are assigned after subclustering via a BFS wave-front: each zero-pop block is assigned to the nearest (by centroid distance) adjacent subcluster. Fully isolated zero-pop groups form their own subclusters.
- Subcluster adjacency uses sphere convex hull filtering: block-level rook pairs are promoted to subcluster edges only if the subcluster centroids are also connected in the 3D convex hull of all centroids projected onto the unit sphere (equivalent to spherical Delaunay triangulation). Subclusters pruned out entirely fall back to their closest rook neighbor.
- Bridge edges added via Kruskal's MST on exterior-ring-filtered nodes to connect disconnected subcluster components before the final METIS pass.
- Final METIS k-way pass (uniform edge weights, k = n\_districts) on subclusters produces district assignments.
- Adjacency visualisation plot after bridging: subcluster polygons colored by ID, rook edges in blue, bridge edges in red.
- Outputs: district GeoJSON, subcluster GeoJSON with component IDs, deviation log, and full run/assignment/geometry rows written to DB.

**Two-pass block METIS pipeline** (`scripts/run_two_pass_metis.py`)
- Pass 1: METIS k-way on all census blocks → clusters. Pass 2: METIS k-way on clusters → districts.
- `rook_adjacency` snapshot taken before synthetic bridge edges are added; repair step uses rook-only adjacency to avoid treating bridged-but-non-touching blocks as contiguous (fixes MultiPolygon cluster geometry bug on states with island blocks such as Illinois/Kaskaskia).
- `_add_bridge_edges` accepts `node_geoms` parameter: unions component geometries, extracts exterior rings via `_iter_exteriors` helper (handles Polygon/MultiPolygon/GeometryCollection), and restricts bridge candidates to nodes whose geometry intersects the outer boundary. Prevents bay-crossing bridge edges (e.g. Chesapeake Bay).
- Zero-pop block assignment replaced with BFS wave-front from populated clusters outward through rook adjacency; isolated zero-pop groups form their own clusters.
- Cluster-level bridge call passes `node_geoms` for exterior-ring-aware candidate selection.

### Fixed

- Two-pass pipeline: synthetic bridge edges no longer contaminate the contiguity repair check, eliminating MultiPolygon district geometries caused by geographically separate blocks being placed in the same cluster by METIS.
- Tract pipeline bisection loop: subclusters that METIS cannot split below threshold no longer cause an infinite retry loop.

### Stashed (not yet merged)

- `git stash "took too long to load"`: initial subcluster construction using `_rook_components` per tract to detect non-contiguous populated blocks. Replaced by geometry explosion approach which is faster and handles zero-pop bridge blocks correctly.

---

## Earlier work (on `main` / prior branches)

### H3-based pipeline (`scripts/run_h3_metis.py` and related)
- Adaptive H3-block aggregation pipeline: aggregate census blocks to H3 cells at a target resolution, with adjacency computed via `ST_Touches` with shared-edge (rook) requirement.
- Post-pass: zero-pop bridge leaves connected via BFS within top-level H3 scope.
- Boundary cleanup, partition visualisation, block-union geometries.
- Fixed adjacency for mixed-resolution H3 cells; excluded zero-pop blocks from graph; fixed Decimal cast from psycopg2.

### Core redistricting app (`src/redistrict/`)
- PyMETIS wrapper (`partition.py`) with fixed seed=42, configurable `ncuts`/`niter`, contiguity enforcement for k-way.
- `db.py`: block/tract/county fetch, rook adjacency cache (`compute_and_store_adjacency_bulk` using `ST_Touches` + `ST_Dimension(ST_Intersection) >= 1`), run/assignment/geometry write, GeoJSON export.
- `--continue` workflow: store all edges in DB, auto-bridge disconnected components after manual QGIS edits.
- Zero-pop node exclusion from graph; post-assign to nearest active node.
- Water penalty sweep experiments; multi-formula experiments across NY tracts/block groups/counties.

### Data preparation scripts
- `create_block_groups.py`, `create_tracts.py`, `create_counties.py`: aggregate Census geometries and populations into DB tables.
- All-states batch runner with per-state logging and skip-on-failure.
