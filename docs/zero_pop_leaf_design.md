# Zero-Pop Leaf Design

Branch: `feature/h3-blocks`

## Problem

The H3 leaf-adjacency pipeline produced 60 disconnected components for Rhode Island, crashing
PyMETIS (which requires a connected input graph). The root cause: zero-pop geographic areas
(water, parks, uninhabited land) have no populated H3 cells, so no leaves and no edges cross
them. Populated clusters on different sides of Narragansett Bay were completely isolated in the
leaf graph.

## Why a Post-Pass Fails

An early design called `find_zero_pop_leaves` to run after `compute_leaves_and_provisional_edges`,
adding zero-pop siblings of leaf ancestors as a separate step. This fails because:

1. Provisional edges are already emitted during leaf marking, filtered to `cell_pop` neighbors.
   Zero-pop cells never appear as provisional edge targets, so even if we add zero-pop leaves
   after the fact, Phase 3 (`resolve_edges`) has no provisional edges to confirm between them
   and their populated neighbors.

2. It adds global complexity without being anchored to the level-by-level traversal, making it
   harder to reason about scope and convergence.

## Integrated Per-Level Algorithm

Zero-pop bridge leaves are determined **inside** `compute_leaves_and_provisional_edges`, inline
with the regular traversal.

### Two scenarios at each level R

**Scenario 1 — populated neighbor.** When a leaf is marked at res R, its parent (at res R-1)
may have populated neighbors at res R-1. Those neighbors are already in `populated_at_res[R-1]`
and will be evaluated naturally in the next loop iteration. No special handling required.

**Scenario 2 — zero-pop neighbor.** When a leaf is marked at res R, its parent's zero-pop
neighbors at res R-1 are queued as candidates:

```python
if res > 0:
    parent = h3.cell_to_parent(cell, res - 1)
    for nbr_parent in h3.grid_disk(parent, 1) - {parent}:
        if (nbr_parent not in cell_pop
                and nbr_parent not in has_leaf_descendant
                and nbr_parent not in leaves):
            zero_pop_candidates.setdefault(res - 1, set()).add(nbr_parent)
```

At each level R, after processing populated cells, these zero-pop candidates are evaluated:

```python
for P in zero_pop_candidates.get(res, set()):
    if P not in has_leaf_descendant and P not in leaves:
        _mark_leaf(P)
```

Zero-pop leaves have no population, so `cell_pop[P] > threshold` is never true. They cannot
trigger their children to split. `_mark_leaf(P)` just adds P to leaves, marks its ancestors,
emits provisional edges, and queues P's parent's zero-pop neighbors at the next coarser level.
This cascades upward.

### Provisional edge condition broadened

The original condition `if neighbor in cell_pop and neighbor not in has_leaf_descendant` was
too restrictive: zero-pop neighbors could never receive edges. The new condition:

```
if neighbor not in has_leaf_descendant:
```

This allows provisional edges to reach zero-pop cells. Phase 3 (`resolve_edges`) then:
- Confirms if the neighbor is a leaf (zero-pop bridge)
- Redirects to the neighbor's leaf ancestor (Case 3: finer cell inside coarser zero-pop leaf)
- Discards if no leaf ancestor exists

### Convergence and top-level scope

**Convergence is defined by populated leaves only.** Zero-pop bridge cells do not influence when
the expansion stops, because they are connective tissue, not substantive units. The stopping
condition at each level R:

```python
pop_leaves = {L for L in leaves if L in cell_pop}
ancestors = {cell_to_parent(L, R) if res(L) > R else L for L in pop_leaves}
if len(ancestors) == 1:
    top_level = ancestors; break
```

When all populated leaves share a single ancestor at resolution R, that ancestor becomes the
`top_level`. Zero-pop bridge leaves accumulate naturally within the scope of that ancestor;
out-of-scope strays (from neighboring states or unrelated geography) are created but then
discarded by the scope filter.

**Scope filter:** After the main loop, leaves and provisional edges are filtered:

```python
def _in_scope(cell):
    tl_res = get_resolution(top_level_cell)
    cell_res = get_resolution(cell)
    if cell_res > tl_res:
        return cell_to_parent(cell, tl_res) == top_level_cell
    return cell == top_level_cell

leaves = {L for L in leaves if _in_scope(L)}
provisional_edges = {(A, N) for A, N in provisional_edges if _in_scope(A) and _in_scope(N)}
```

This removes zero-pop bridge cells that cascaded into neighboring states before convergence
was reached.

## Example: Rhode Island

RI has populated blocks spread across many res-14 cells. Narragansett Bay separates the eastern
and western parts. At res=14, there are many populated parents with different ancestors → no
convergence at fine resolutions. The zero-pop expansion bridges the Bay: water cells at various
resolutions become bridge leaves, creating edges that span the water. Convergence eventually
occurs at res=2 or res=3 (all of RI within a single ancestor cell). The scope filter discards
any MA/CT cells that crept in during the expansion.

## Effect on `run_h3_state.py`

- After `weighted_centroids`, zero-pop bridge leaves are appended to `cell_nodes` with pop=0
  and centroid from `h3.cell_to_latlng`.
- `show_leaves` is called **before** PyMETIS, displaying all leaves (populated + bridges)
  colored by resolution.
- Zero-pop bridge nodes in the METIS graph have pop=0. After partitioning, those nodes are
  assigned to some district but no geoids map to them, so they don't appear in the final
  `geoid_to_district` mapping. District geometries exclude them (only geoid-assigned cells
  contribute to the polygon union).

## Files Changed

| File | Change |
|------|--------|
| `src/redistrict/h3_graph.py` | `compute_leaves_and_provisional_edges`: inline zero-pop candidates, broadened edge condition, convergence check, scope filter |
| `scripts/run_h3_state.py` | Append zero-pop nodes; move `show_leaves` before METIS |
| `tests/test_h3_graph.py` | `TestIntegratedZeroPopLeaves`, `TestProvisionalEdgeCaseThree` |
