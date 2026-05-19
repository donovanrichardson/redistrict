# H3-Based Census Block Redistricting

## Motivation

The block-group / tract pipeline uses pre-aggregated Census units as graph nodes. These units have widely varying populations — a Manhattan block group can have 5× the population of a rural one — which limits how finely the METIS partitioner can balance districts.

Census blocks are the atomic unit: the smallest geography the Census publishes population for (~8M nationwide, ~288K for NY). Using them directly as METIS nodes would produce extremely large graphs (slow) with many zero-pop nodes (parks, water, industrial land).

The H3 approach solves both problems: adaptively aggregate blocks into hexagonal cells whose population sits at a uniform threshold, then run METIS on those cells. Dense urban areas get many small cells; rural areas get a single large coarse cell.

## Threshold Definition

Sort census blocks descending by population. Walk the sorted list accumulating population until the running total reaches 1% of the state total. The population of the last block included at that point is the **threshold**.

This is the block population at the 99th percentile of population mass — meaning 1% of the state's people live in blocks at least this dense.

Example: 1 block pop=100, 1500 other blocks total pop=9900, state total=10000.
Target = 1% × 10000 = 100. The first block (pop=100) reaches the target. Threshold = 100.

## H3 Assignment and Aggregation

1. Map each block's centroid (intptlat20, intptlon20) to H3 resolution 15 using the Python `h3` library. Resolution 15 has ~0.9 m² average cell area — effectively one cell per point.

2. Bottom-up aggregation from res 14 down to min_res (default 4):
   - Group current cells by their parent at the current resolution.
   - If the parent group's total population ≤ threshold, replace the children with the parent (merge).
   - If not, leave the children as separate nodes.
   - This naturally handles blocks above threshold: a group containing a high-pop block will have total > threshold and won't merge.

3. The result is a set of H3 cells at mixed resolutions:
   - Dense urban cores: many res-15 cells (each covering one or a few blocks)
   - Suburban areas: res-12 to res-10 cells (covering several blocks each)
   - Rural areas: res-8 to res-5 cells (covering entire townships)

## Population-Weighted Centroids

For each aggregated H3 cell, compute a population-weighted centroid from the block centroids it contains:

    lat = Σ(pop_i × lat_i) / Σ(pop_i)
    lon = Σ(pop_i × lon_i) / Σ(pop_i)

Cells with zero total population fall back to the H3 geometric centroid.

## Adjacency Graph

Two aggregated cells are adjacent if any of their constituent res-15 cells are H3 neighbours. The implementation:

1. Build a mapping: res-15 cell → final aggregated ancestor cell.
2. For each res-15 cell, look up its 6 H3 neighbours (also at res-15).
3. If a neighbour maps to a different final cell, add an edge between the two final cells.

This handles mixed resolutions correctly: a coarse res-8 cell is adjacent to any finer cell whose res-15 descendants border a res-15 descendant of the coarse cell.

## METIS Partition

- **Edge weights**: uniform (1). METIS cuts based on topology and population balance only, not geographic distance.
- **Node weights**: cell population (floor: 1).
- `contig=1` enforces geographic contiguity (k-way partitioning).
- `ufactor=8` allows ±0.8% population deviation.

## Disaggregation

After METIS assigns each cell to a district:
1. Map each census block to its aggregated H3 cell.
2. Assign the block to the district of its cell.

District geometries are produced by `ST_Union` over all blocks in each district (using the existing `db.write_district_geoms` with `geography="blocks"`).

## Comparison with Block-Group Pipeline

| Aspect | Block-group pipeline | H3-blocks pipeline |
|--------|---------------------|-------------------|
| Node count (NY) | ~15,700 | ~varies by threshold |
| Geography | Fixed Census units | Adaptive hexagons |
| Urban granularity | Block-group (~1500 pop) | Sub-block-group (threshold) |
| Rural granularity | Block-group | Coarse hexagons |
| Edge weights | Uniform / formula-based | Uniform (1) |
| Water handling | QGIS curation workflow | H3 topology only |
| Adjacency | Urquhart graph | H3 native neighbours |

## Files

- `src/redistrict/h3_graph.py` — core algorithm module
- `tests/test_h3_graph.py` — unit tests (no DB required)
- `scripts/run_h3_state.py` — CLI: `python run_h3_state.py <statefp> <n_districts>`
- `output/h3_runs/` — GeoJSON and deviation logs per run

## Running

```bash
# Rhode Island (2 districts, fast smoke test)
uv run --env-file .env python scripts/run_h3_state.py 44 2

# New York (26 districts)
uv run --env-file .env python scripts/run_h3_state.py 36 26
```
