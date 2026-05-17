# Blog Post Outline: Automated Congressional Redistricting with Graph Partitioning

## What we built and why

- Redistricting is usually done by humans, subject to partisan manipulation. We wanted to see how far a purely algorithmic, geography-driven approach could get.
- The tool takes Census block group (or tract, or county) data and produces congressional district maps for any US state, optimising for geographic compactness and population equality.
- Built entirely in Python on top of open Census data in PostGIS.

## Data foundation

- Census 2020 TIGER/Line geometries at three resolutions: counties, tracts, block groups (~15,000–200,000 units per state).
- Population from the 2020 Decennial Census (PL 94-171 redistricting file).
- Rook contiguity (shared-edge adjacency) computed in PostGIS using `ST_Touches` + `ST_Dimension(ST_Intersection(...)) >= 1` — distinguishes true land borders from point-only touches.
- Zero-population block groups (parks, water bodies, industrial zones) are excluded from the graph and assigned to their nearest populated neighbour after partitioning.

## Graph construction: Spherical Delaunay + Urquhart

- Naive k-nearest-neighbour graphs produce too many long-range edges across water bodies. We needed a sparser, geographically coherent structure.
- Project each unit's centroid onto the unit sphere and compute the 3D convex hull — this gives a spherical Delaunay triangulation via scipy.
- Apply the Urquhart reduction: remove the longest edge (by haversine distance) from each triangle. The result is a planar-ish graph that respects geographic proximity without hard-coding any distance cutoff.
- Non-rook-contiguous edges (crossing water) are kept but penalised in edge weight — METIS can still cut them but prefers not to.

## Edge weighting and the formula search

- Explored four weighting strategies to understand the tradeoff between population awareness and geographic compactness:
  - **original**: `SCALE / (dist/(2√pop_a) + dist/(2√pop_b))` — population-informed; encourages METIS to keep dense areas together.
  - **original_clamped**: same formula but with a constant offset so the weight range is exactly 4:1 — reduces extreme edges without losing the population signal.
  - **blend**: 50% normalised inverse-cost + 50% normalised inverse-distance — a middle ground.
  - **uniform**: all land edges equal weight (`SCALE`), non-adjacent edges divided by a water penalty — METIS cuts on topology alone.
- Ran all four on New York block groups (26 districts). Uniform produced the most compact, geographically coherent districts with the least urban fragmentation.

## Water penalty and QGIS curation

- A configurable `water_penalty` divisor makes cross-water edges cheaper to cut. Tested 1×, 2×, 4×, 8× the base penalty on New York.
- 8× produced the best island/peninsula handling (Long Island, Manhattan, Staten Island assigned sensibly) while still allowing the algorithm to bridge where needed.
- Built a QGIS curation workflow: all Urquhart edges are stored in PostGIS with an `is_adjacent` flag. The user filters to non-adjacent edges in QGIS and deletes any physically impossible links (e.g., open-ocean crossings with no road or ferry). A `--continue <run_id>` CLI command re-runs METIS with the curated edge set.
- If the user deletes too many edges and disconnects the graph, the tool auto-bridges the smallest component back to the main graph via a nearest-neighbour search (bounding-box pre-filtered to ~1% of candidates for performance).

## Partitioning with PyMETIS

- PyMETIS (Python bindings for the METIS multilevel k-way graph partitioner) does the heavy lifting.
- Node weights = population; edge weights = the formula output. METIS minimises total cut weight subject to balanced node weights.
- `contig=1` enforces geographic contiguity (k-way only; recursive bisection silently ignores it).
- `ufactor=8` allows ±0.8% population deviation per district — tight but achievable at block-group granularity.
- Population deviations across all tested states are consistently within ±1.5%, comparable to hand-drawn maps.

## Results

- New York (26 districts, block groups): uniform + 8× water penalty + QGIS-curated water links gives clean, compact districts respecting borough and county boundaries where population allows.
- All 41 continental states with ≥2 congressional seats ran successfully with the same settings (runs 45–85), maps produced in a single batch. 7 states skipped — all have only 1 House seat (Delaware, DC, Montana, North Dakota, South Dakota, Vermont, Wyoming).
- 38 of 41 states achieve worst-case deviation ≤0.8%. Ohio reaches 1.7%; California 8.9% and Texas 14.3% are outliers — both have very large-population urban block groups that create a hard floor on achievable balance at this geography level.
- Zero-population units (parks, water, industrial) — ranging from 0 (Nebraska, West Virginia) to 170 (Michigan) per state — are correctly excluded from the partitioning graph and post-assigned to their nearest populated neighbour.

## What this doesn't do (yet)

- Partisan fairness metrics (efficiency gap, mean-median, etc.) — the algorithm is neutral by construction but the output hasn't been audited against these.
- Preserving existing political subdivisions (counties, municipalities) as soft constraints.
- Handling multi-member or ranked-choice districts.
- Alaska and Hawaii (non-contiguous; the spherical Delaunay would connect them to the mainland).
