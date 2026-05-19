"""
H3-based graph construction for redistricting directly from census blocks.

Pipeline:
  1. compute_threshold    - block pop at the 99th-percentile point by pop mass
  2. assign_h3_res15      - map each block centroid to H3 resolution 15
  3. aggregate_h3_cells   - bottom-up: find coarsest res where cell pop ≤ threshold
  4. weighted_centroids   - pop-weighted centroid per aggregated cell
  5. build_h3_adjacency   - adjacency via res-15 neighbour lookup
  6. build_metis_graph    - uniform edge weight=1, node weight=pop

Threshold definition:
  Sort blocks descending by population. Walk the list accumulating population
  until the running total reaches 1% of the state total. The population of the
  last block added at that point is the threshold — i.e. the block population
  at the 99th percentile of population mass.

  Example: 1 block pop=100, 1500 blocks each pop≈6.6, total=10000.
  1% target = 100. The first block (pop=100) hits the target → threshold = 100.

Aggregation:
  Blocks above threshold stay as individual res-15 nodes. All others are merged
  bottom-up: at each resolution from 14 down to min_res, candidate cells are
  grouped by their parent. If the parent's total pop ≤ threshold the children
  are replaced by the parent. This continues until no further merging is
  possible at the current resolution.
"""

from __future__ import annotations

import math
from typing import Sequence

import h3


MIN_CELL_POP = 1     # floor used as node weight in METIS
DEFAULT_MIN_RES = 4  # coarsest H3 resolution permitted (avg ~1730 km²)
THRESHOLD_PCT = 0.01 # fraction of state pop defining the threshold


# ---------------------------------------------------------------------------
# Threshold
# ---------------------------------------------------------------------------

def compute_threshold(blocks: Sequence[dict], pct: float = THRESHOLD_PCT) -> int:
    """
    Return the population of the block at the (1-pct) quantile of pop mass.

    Sort blocks descending by pop, accumulate until sum >= pct * total.
    The population of the last block included is returned as the threshold.
    """
    total = sum(b["pop"] for b in blocks)
    if total == 0:
        return 1
    target = pct * total
    sorted_pops = sorted((b["pop"] for b in blocks if b["pop"] > 0), reverse=True)
    accumulated = 0
    for pop in sorted_pops:
        accumulated += pop
        if accumulated >= target:
            return pop
    return sorted_pops[-1] if sorted_pops else 1


# ---------------------------------------------------------------------------
# H3 assignment
# ---------------------------------------------------------------------------

def assign_h3_res15(blocks: Sequence[dict]) -> dict[str, str]:
    """
    Map geoid -> H3 resolution-15 cell for each block's centroid.

    Each block dict must have 'geoid', 'lat', 'lon'.
    """
    return {
        b["geoid"]: h3.latlng_to_cell(float(b["lat"]), float(b["lon"]), 15)
        for b in blocks
    }


# ---------------------------------------------------------------------------
# Bottom-up aggregation
# ---------------------------------------------------------------------------

def aggregate_h3_cells(
    geoid_to_res15: dict[str, str],
    block_pops: dict[str, int],
    threshold: int,
    min_res: int = DEFAULT_MIN_RES,
) -> dict[str, str]:
    """
    Return geoid -> final H3 cell after bottom-up aggregation.

    Starting at res 15, iteratively replace groups of sibling cells with
    their common parent whenever the parent's total population ≤ threshold.
    Stops at min_res.

    Cells whose individual population already exceeds threshold remain at
    res 15 (they prevent their parent from being merged too).
    """
    # Build: current_cell -> list[geoid]
    current: dict[str, list[str]] = {}
    for geoid, cell in geoid_to_res15.items():
        current.setdefault(cell, []).append(geoid)

    # cell -> total pop (sum of blocks inside)
    cell_pop: dict[str, int] = {
        cell: sum(block_pops.get(g, 0) for g in geoids)
        for cell, geoids in current.items()
    }

    for res in range(14, min_res - 1, -1):
        # Group current cells by parent at this resolution
        parent_to_children: dict[str, list[str]] = {}
        for cell in current:
            parent = h3.cell_to_parent(cell, res)
            parent_to_children.setdefault(parent, []).append(cell)

        for parent, children in parent_to_children.items():
            total = sum(cell_pop[c] for c in children)
            if total <= threshold:
                merged: list[str] = []
                for c in children:
                    merged.extend(current.pop(c))
                    del cell_pop[c]
                current[parent] = merged
                cell_pop[parent] = total

    # Build final geoid -> cell mapping
    result: dict[str, str] = {}
    for cell, geoids in current.items():
        for g in geoids:
            result[g] = cell
    return result


# ---------------------------------------------------------------------------
# Population-weighted centroids
# ---------------------------------------------------------------------------

def weighted_centroids(
    geoid_to_cell: dict[str, str],
    blocks_by_geoid: dict[str, dict],
) -> list[dict]:
    """
    Compute a population-weighted centroid for each aggregated H3 cell.

    Returns a list of node dicts: {"cell", "pop", "lat", "lon"}.
    Falls back to the H3 geometric centroid when a cell contains only
    zero-population blocks.
    """
    cell_to_geoids: dict[str, list[str]] = {}
    for geoid, cell in geoid_to_cell.items():
        cell_to_geoids.setdefault(cell, []).append(geoid)

    nodes: list[dict] = []
    for cell, geoids in cell_to_geoids.items():
        total_pop = sum(blocks_by_geoid[g]["pop"] for g in geoids)
        if total_pop > 0:
            lat = (
                sum(blocks_by_geoid[g]["pop"] * float(blocks_by_geoid[g]["lat"])
                    for g in geoids) / total_pop
            )
            lon = (
                sum(blocks_by_geoid[g]["pop"] * float(blocks_by_geoid[g]["lon"])
                    for g in geoids) / total_pop
            )
        else:
            lat, lon = h3.cell_to_latlng(cell)
        nodes.append({"cell": cell, "pop": total_pop, "lat": lat, "lon": lon})
    return nodes


# ---------------------------------------------------------------------------
# Adjacency
# ---------------------------------------------------------------------------

def build_h3_adjacency(
    geoid_to_cell: dict[str, str],
    geoid_to_res15: dict[str, str],
) -> set[tuple[str, str]]:
    """
    Build adjacency between final (possibly mixed-resolution) H3 cells.

    Strategy: map each res-15 cell to its final aggregated ancestor, then
    for each res-15 cell inspect its 6 H3 neighbours.  If a neighbour maps
    to a different final cell, the two final cells are adjacent.

    This correctly handles mixed resolutions: a coarse cell is adjacent to
    any fine cell whose res-15 representative neighbours a res-15 descendant
    of the coarse cell.
    """
    # res-15 cell -> final cell
    res15_to_final: dict[str, str] = {}
    for geoid, res15 in geoid_to_res15.items():
        res15_to_final[res15] = geoid_to_cell[geoid]

    edges: set[tuple[str, str]] = set()
    for res15_cell, final_cell in res15_to_final.items():
        for neighbour in set(h3.grid_disk(res15_cell, 1)):
            if neighbour == res15_cell:
                continue
            neighbour_final = res15_to_final.get(neighbour)
            if neighbour_final and neighbour_final != final_cell:
                a, b = min(final_cell, neighbour_final), max(final_cell, neighbour_final)
                edges.add((a, b))
    return edges


# ---------------------------------------------------------------------------
# METIS graph
# ---------------------------------------------------------------------------

def build_metis_graph(
    cell_nodes: list[dict],
    adjacency: set[tuple[str, str]],
) -> tuple[list[list[int]], list[int], list[int]]:
    """
    Build PyMETIS-ready structures from H3 cell nodes and their adjacency.

    All edges have weight 1 (uniform — METIS cuts on topology/population only).
    Node weights are cell populations.

    Returns (adjacency_lists, eweights, nweights).
    """
    cell_to_idx = {n["cell"]: i for i, n in enumerate(cell_nodes)}
    n = len(cell_nodes)

    adj_map: list[dict[int, int]] = [{} for _ in range(n)]
    for ca, cb in adjacency:
        ia = cell_to_idx.get(ca)
        ib = cell_to_idx.get(cb)
        if ia is None or ib is None:
            continue
        adj_map[ia][ib] = 1
        adj_map[ib][ia] = 1

    adjacency_lists: list[list[int]] = []
    eweights: list[int] = []
    for i in range(n):
        neighbours = sorted(adj_map[i].keys())
        adjacency_lists.append(neighbours)
        for _ in neighbours:
            eweights.append(1)

    nweights = [max(MIN_CELL_POP, node["pop"]) for node in cell_nodes]
    return adjacency_lists, eweights, nweights


# ---------------------------------------------------------------------------
# Connectivity helpers (reuse from graph module logic)
# ---------------------------------------------------------------------------

def check_connectivity(
    cell_nodes: list[dict],
    adjacency: set[tuple[str, str]],
) -> list[list[str]]:
    """
    Return connected components as lists of cell ids.
    A single-element list means the graph is fully connected.
    """
    cell_ids = [n["cell"] for n in cell_nodes]
    adj: dict[str, list[str]] = {c: [] for c in cell_ids}
    for a, b in adjacency:
        if a in adj:
            adj[a].append(b)
        if b in adj:
            adj[b].append(a)

    visited: set[str] = set()
    components: list[list[str]] = []
    for start in cell_ids:
        if start in visited:
            continue
        component: list[str] = []
        stack = [start]
        while stack:
            v = stack.pop()
            if v in visited:
                continue
            visited.add(v)
            component.append(v)
            stack.extend(adj[v])
        components.append(component)
    return components
