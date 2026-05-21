# H3 Leaf-Adjacency Pipeline

Branch: `feature/h3-leaf-adjacency`

## Motivation

The original H3 pipeline (`feature/h3-blocks`) used a bottom-up merge approach: group sibling cells by parent, merge if parent pop ≤ threshold. Adjacency was then computed as a separate post-processing step with three cases (same-resolution, coarser ancestor, finer descendants).

The leaf-adjacency pipeline unifies these two phases. Leaves are determined by a single traversal; provisional edges are emitted during that same traversal and resolved in one postprocessing pass. The result is structurally cleaner and guarantees that adjacency is derived from the same data that produced the leaf nodes.

## Algorithm

### Threshold

Same as before: sort blocks descending by pop, walk until cumulative sum ≥ 1% of state total, return that block's population. This is the block population at the 99th percentile of population mass.

### Phase 0 — Build cell hierarchy

Walk from res-15 upward to res-0, accumulating population and parent-child links. Only cells with population > 0 are included at any level.

```
cell_pop[res15_cell] = sum of block pops mapped to that cell
for res 14 down to 0:
    for each populated cell at res+1:
        parent = cell_to_parent(child, res)
        cell_pop[parent] += cell_pop[child]
        cell_direct_children[parent].add(child)
        populated_at_res[res].add(parent)
```

### Phase 1 — Leaf determination (fine to coarse)

```
has_leaf_descendant = set()   # ancestors of confirmed leaves
leaves = set()

for res R from 14 down to 0:
    for each populated cell P at res R (not in has_leaf_descendant):
        effective_children = [c for c in direct_children[P]
                              if c not in has_leaf_descendant]
        if cell_pop[P] > threshold and len(effective_children) > 1:
            for each child in effective_children:
                mark_leaf(child)

# Cleanup: any res-0 cell still with no leaf descendants becomes a leaf
for P in populated_at_res[0]:
    if P not in has_leaf_descendant and P not in leaves:
        mark_leaf(P)
```

**mark_leaf(cell)**:
1. Add cell to `leaves`
2. Walk ancestor chain (res-1 down to 0), adding each to `has_leaf_descendant` (stop early if already marked)
3. Emit provisional edges to same-resolution populated neighbors that have no leaf descendants

### Phase 2 — Provisional edges (triggered inside mark_leaf)

When a leaf cell is marked, check its H3 ring-1 neighbors at the same resolution:

```
for neighbor in grid_disk(cell, 1) - {cell}:
    if neighbor in cell_pop and neighbor not in has_leaf_descendant:
        provisional_edges.add((cell, neighbor))   # directed: cell is the leaf
```

Neighbors with `has_leaf_descendant` are skipped because their own fine-grained leaves will emit edges when they are marked.

### Phase 3 — Resolve provisional edges

Each provisional edge `(A, N)` is directed: A is the originating leaf, N is a raw H3 neighbor at the same resolution.

```
for each (A, N) in provisional_edges:
    if N in leaves:
        confirm edge(A, N)
    elif any ancestor of N is in leaves:
        confirm edge(A, leaf_ancestor(N))
    else:
        discard
```

The "leaf ancestor" case handles mixed resolutions: N may be a fine-grained cell inside a coarser leaf, meaning N's area was aggregated to a coarser level.

## Leaf semantics

A cell is a leaf when:
- Its parent's population exceeds the threshold AND the parent has more than one populated effective child (cells with no leaf descendants).

Dense urban areas produce many fine-grained leaves (res-15 or res-14). Rural areas collapse to coarse leaves (res-8 to res-4). Isolated single-child chains percolate to res-0.

## District geometries in Python

District polygons are computed entirely in Python using `h3.cell_to_boundary` + Shapely `unary_union`. No PostGIS geometry computation is required.

```python
for each district:
    polys = [Polygon(cell_to_boundary(cell)) for cell in district_cells]
    geom  = MultiPolygon(unary_union(polys))
```

Progress is shown with tqdm. The result is inserted as WKT via `ST_GeomFromText`.

## Files

| File | Role |
|------|------|
| `src/redistrict/h3_graph.py` | `build_cell_hierarchy`, `compute_leaves_and_provisional_edges`, `resolve_edges`, `assign_geoids_to_leaves`, `build_district_geoms` |
| `src/redistrict/db.py` | `write_district_geoms_wkt` — insert pre-computed WKT district geometries |
| `tests/test_h3_graph.py` | `TestBuildCellHierarchy`, `TestComputeLeavesAndProvisionalEdges`, `TestResolveEdges`, `TestAssignGeooidsToLeaves` (48 tests total) |
| `scripts/run_h3_state.py` | Updated to use new pipeline; tqdm on zero-pop assignment and district geometry |

## Key design decisions

**Directed provisional edges.** Provisional edges are stored as `(leaf, neighbor)` ordered pairs, not sorted tuples. Phase 3 needs to know which endpoint is the confirmed leaf (A) to correctly redirect the other endpoint (N) to its leaf ancestor. Sorting would lose this distinction.

**`has_leaf_descendant` early-break.** When marking ancestors of a new leaf, iteration stops as soon as an ancestor is already in `has_leaf_descendant`. By induction, all further ancestors are already marked from a prior `mark_leaf` call.

**Only populated cells evaluated.** `build_cell_hierarchy` and all subsequent phases operate exclusively on cells with population > 0. Empty H3 cells are never instantiated.

**District geometry uses H3 cell boundaries, not census block boundaries.** The district polygons are hexagonal approximations of the true block-union boundaries. For redistricting visualization this is sufficient; the authoritative geography is in the block assignment table.

## Running

```bash
uv run --env-file .env python scripts/run_h3_state.py 44 2   # Rhode Island, 2 districts
uv run --env-file .env python scripts/run_h3_state.py 36 26  # New York, 26 districts
```
