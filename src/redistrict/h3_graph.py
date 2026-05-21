"""
H3-based graph construction for redistricting directly from census blocks.

Leaf-adjacency pipeline:
  1. compute_threshold             - block pop at the 99th-percentile point by pop mass
  2. assign_h3_res15               - map each block centroid to H3 resolution 15
  3. build_cell_hierarchy          - accumulate pop and parent-child links (populated cells only)
  4. compute_leaves_and_provisional_edges - bottom-up leaf determination + provisional edges
  5. resolve_edges                 - confirm/redirect/discard provisional edges
  6. assign_geoids_to_leaves       - map each geoid to its leaf cell
  7. weighted_centroids            - pop-weighted centroid per leaf cell
  8. build_metis_graph             - uniform edge weight=1, node weight=pop

Threshold definition:
  Sort blocks descending by population. Walk the sorted list accumulating population
  until the running total reaches 1% of the state total. The population of the
  last block included is returned as the threshold — i.e. the block population
  at the 99th percentile of population mass.

Leaf definition:
  A cell C at resolution R is a LEAF when its parent P satisfies:
    cell_pop[P] > threshold  AND  P has >1 populated effective children (no leaf descendants).
  Remaining res-0 cells with no leaf descendants are also marked as leaves.
"""

from __future__ import annotations

import math
from typing import Sequence

import h3
from shapely.geometry import MultiPolygon, Polygon
from shapely.ops import unary_union
from tqdm import tqdm


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
    total = sum(int(b["pop"]) for b in blocks)
    if total == 0:
        return 1
    target = pct * total
    sorted_pops = sorted((int(b["pop"]) for b in blocks if b["pop"] > 0), reverse=True)
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
    # Build: current_cell -> list[geoid], excluding zero-pop blocks
    current: dict[str, list[str]] = {}
    for geoid, cell in geoid_to_res15.items():
        if int(block_pops.get(geoid, 0)) > 0:
            current.setdefault(cell, []).append(geoid)

    # cell -> total pop (sum of blocks inside)
    cell_pop: dict[str, int] = {
        cell: sum(int(block_pops.get(g, 0)) for g in geoids)
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
# Leaf-adjacency pipeline (replaces aggregate_h3_cells + build_h3_adjacency)
# ---------------------------------------------------------------------------

def build_cell_hierarchy(
    geoid_to_res15: dict[str, str],
    block_pops: dict[str, int],
) -> tuple[dict[str, int], dict[str, set[str]], dict[int, set[str]]]:
    """
    Accumulate population and parent-child relationships upward from res-15.
    Only cells with population > 0 are included at any level.

    Returns:
        cell_pop: cell -> total population
        cell_direct_children: cell -> set of direct populated children (at res+1)
        populated_at_res: resolution -> set of populated cells at that resolution
    """
    cell_pop: dict[str, int] = {}
    for geoid, cell in geoid_to_res15.items():
        pop = block_pops.get(geoid, 0)
        if pop > 0:
            cell_pop[cell] = cell_pop.get(cell, 0) + pop

    populated_at_res: dict[int, set[str]] = {15: set(cell_pop.keys())}
    cell_direct_children: dict[str, set[str]] = {}

    for res in range(14, -1, -1):
        level: set[str] = set()
        for child in populated_at_res.get(res + 1, set()):
            parent = h3.cell_to_parent(child, res)
            cell_pop[parent] = cell_pop.get(parent, 0) + cell_pop[child]
            cell_direct_children.setdefault(parent, set()).add(child)
            level.add(parent)
        populated_at_res[res] = level

    return cell_pop, cell_direct_children, populated_at_res


def compute_leaves_and_provisional_edges(
    cell_pop: dict[str, int],
    cell_direct_children: dict[str, set[str]],
    populated_at_res: dict[int, set[str]],
    threshold: int,
) -> tuple[set[str], set[tuple[str, str]]]:
    """
    Phase 1 (fine to coarse, res 14 -> 0):
      For each populated cell P with no leaf descendants:
        If cell_pop[P] > threshold and P has >1 effective children → mark children as leaves.
    Phase 2 (triggered on mark_leaf):
      Add provisional edges to same-resolution populated neighbors not yet a leaf ancestor.
    Cleanup:
      Remaining res-0 cells with no leaf descendants become leaves.

    Returns (leaves, provisional_edges). Only populated cells become leaves here.
    Zero-pop bridge leaves are added separately via find_zero_pop_leaves.
    """
    leaves: set[str] = set()
    has_leaf_descendant: set[str] = set()
    provisional_edges: set[tuple[str, str]] = set()

    def _mark_leaf(cell: str) -> None:
        leaves.add(cell)
        res = h3.get_resolution(cell)
        for r in range(res - 1, -1, -1):
            anc = h3.cell_to_parent(cell, r)
            if anc in has_leaf_descendant:
                break
            has_leaf_descendant.add(anc)
        for neighbor in set(h3.grid_disk(cell, 1)) - {cell}:
            if neighbor in cell_pop and neighbor not in has_leaf_descendant:
                provisional_edges.add((cell, neighbor))

    for res in range(14, -1, -1):
        for P in populated_at_res.get(res, set()):
            if P in has_leaf_descendant:
                continue
            effective_children = [
                c for c in cell_direct_children.get(P, set())
                if c not in has_leaf_descendant
            ]
            if cell_pop[P] > threshold and len(effective_children) > 1:
                for child in effective_children:
                    _mark_leaf(child)

    for P in populated_at_res.get(0, set()):
        if P not in has_leaf_descendant and P not in leaves:
            _mark_leaf(P)

    return leaves, provisional_edges


def compute_top_level(leaves: set[str]) -> set[str]:
    """
    Find the tightest resolution at which all leaves share a single common ancestor.
    Returns that single-element ancestor set, or all res-0 ancestors if none found.
    """
    if not leaves:
        return set()
    prev: set[str] = set()
    for res in range(0, 16):
        cur: set[str] = set()
        for L in leaves:
            L_res = h3.get_resolution(L)
            cur.add(h3.cell_to_parent(L, res) if L_res > res else L)
        if len(cur) > 1:
            return prev if prev else cur
        prev = cur
    return prev


def find_zero_pop_leaves(
    populated_leaves: set[str],
    top_level: set[str],
) -> set[str]:
    """
    Post-pass: within the scope of top_level, every cell that is not a populated leaf
    and has no populated-leaf descendant becomes a zero-pop leaf at the coarsest
    resolution that is not an ancestor of a populated leaf.

    BFS from each top_level cell downward. A cell is either:
      - a populated leaf (terminal, no children needed)
      - an ancestor of a populated leaf (expand its children)
      - otherwise a zero-pop leaf (terminal)
    """
    has_leaf_descendant: set[str] = set()
    for L in populated_leaves:
        res = h3.get_resolution(L)
        for r in range(res - 1, -1, -1):
            anc = h3.cell_to_parent(L, r)
            if anc in has_leaf_descendant:
                break
            has_leaf_descendant.add(anc)

    zero_pop_leaves: set[str] = set()

    for top_cell in top_level:
        top_res = h3.get_resolution(top_cell)
        if top_cell in populated_leaves:
            continue
        if top_cell not in has_leaf_descendant:
            zero_pop_leaves.add(top_cell)
            continue
        frontier: list[str] = list(h3.cell_to_children(top_cell, top_res + 1))
        for res in range(top_res + 1, 16):
            next_frontier: list[str] = []
            for cell in frontier:
                if cell in populated_leaves:
                    pass
                elif cell in has_leaf_descendant:
                    if res < 15:
                        next_frontier.extend(h3.cell_to_children(cell, res + 1))
                else:
                    zero_pop_leaves.add(cell)
            frontier = next_frontier

    return zero_pop_leaves


def resolve_edges(
    provisional_edges: set[tuple[str, str]],
    leaves: set[str],
) -> set[tuple[str, str]]:
    """
    Phase 3: For each directed provisional_edge(A, N) where A is the originating leaf:
      - N is a leaf              → confirm edge(A, N)
      - N has a leaf ancestor    → redirect to edge(A, leaf_ancestor)
      - else                     → discard

    Returns undirected edges as sorted (min, max) tuples.
    """
    def _leaf_ancestor(cell: str) -> str | None:
        res = h3.get_resolution(cell)
        for r in range(res - 1, -1, -1):
            anc = h3.cell_to_parent(cell, r)
            if anc in leaves:
                return anc
        return None

    confirmed: set[tuple[str, str]] = set()
    for A, N in provisional_edges:
        if N in leaves:
            if A != N:
                confirmed.add((min(A, N), max(A, N)))
        else:
            anc = _leaf_ancestor(N)
            if anc is not None and anc != A:
                confirmed.add((min(A, anc), max(A, anc)))
    return confirmed


def assign_geoids_to_leaves(
    geoid_to_res15: dict[str, str],
    leaves: set[str],
) -> dict[str, str]:
    """Map each geoid to its leaf cell by walking the ancestor chain."""
    result: dict[str, str] = {}
    for geoid, res15_cell in geoid_to_res15.items():
        if res15_cell in leaves:
            result[geoid] = res15_cell
            continue
        res = h3.get_resolution(res15_cell)
        for r in range(res - 1, -1, -1):
            anc = h3.cell_to_parent(res15_cell, r)
            if anc in leaves:
                result[geoid] = anc
                break
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
        total_pop = sum(int(blocks_by_geoid[g]["pop"]) for g in geoids)
        if total_pop > 0:
            lat = (
                sum(int(blocks_by_geoid[g]["pop"]) * float(blocks_by_geoid[g]["lat"])
                    for g in geoids) / total_pop
            )
            lon = (
                sum(int(blocks_by_geoid[g]["pop"]) * float(blocks_by_geoid[g]["lon"])
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
    min_res: int = DEFAULT_MIN_RES,
) -> set[tuple[str, str]]:
    """
    Build adjacency between final (possibly mixed-resolution) H3 cells.

    For each final cell C, H3 neighbours are fetched at C's own resolution.
    Three cases handle mixed resolutions:

      1. Neighbour N is directly in the final cell set — same-resolution edge.
      2. N is not in the set but one of its coarser ancestors is — the blocks
         in N were merged upward into that ancestor.
      3. N is not in the set but it is a coarser parent of finer final cells —
         the blocks inside N stayed at finer resolutions (dense area).
    """
    cell_set: set[str] = set(geoid_to_cell.values())

    # ancestor_cell -> set of final cells that are its descendants (case 3)
    parent_to_finals: dict[str, set[str]] = {}
    for cell in cell_set:
        res = h3.get_resolution(cell)
        for r in range(min_res, res):
            anc = h3.cell_to_parent(cell, r)
            parent_to_finals.setdefault(anc, set()).add(cell)

    edges: set[tuple[str, str]] = set()

    def _add(a: str, b: str) -> None:
        if a != b:
            edges.add((min(a, b), max(a, b)))

    for final_cell in cell_set:
        for neighbour in set(h3.grid_disk(final_cell, 1)) - {final_cell}:
            # Case 1: direct neighbour at same resolution
            if neighbour in cell_set:
                _add(final_cell, neighbour)
                continue

            # Case 2: neighbour was merged into a coarser ancestor
            n_res = h3.get_resolution(neighbour)
            for r in range(n_res - 1, min_res - 1, -1):
                anc = h3.cell_to_parent(neighbour, r)
                if anc in cell_set:
                    _add(final_cell, anc)
                    break

            # Case 3: neighbour is a coarser placeholder containing finer finals
            if neighbour in parent_to_finals:
                for desc in parent_to_finals[neighbour]:
                    _add(final_cell, desc)

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


# ---------------------------------------------------------------------------
# District geometry (Python/Shapely, no PostGIS required)
# ---------------------------------------------------------------------------

def build_district_geoms(
    cell_to_district: dict[str, int],
    pop_per_district: dict[int, int],
) -> dict[int, tuple[str, int]]:
    """
    Compute district MultiPolygon geometries in Python using H3 cell boundaries.

    For each district, union the Shapely polygons of all its leaf cells.
    Returns district_id -> (wkt_geometry, population).
    """
    district_to_cells: dict[int, list[str]] = {}
    for cell, dist in cell_to_district.items():
        district_to_cells.setdefault(dist, []).append(cell)

    result: dict[int, tuple[str, int]] = {}
    for dist_id, cells in tqdm(sorted(district_to_cells.items()),
                               desc="Building district geometries", unit="district"):
        polys: list[Polygon] = []
        for cell in cells:
            boundary = h3.cell_to_boundary(cell)
            # h3 returns (lat, lon) pairs; Shapely expects (lon, lat)
            coords = [(lon, lat) for lat, lon in boundary]
            polys.append(Polygon(coords))
        geom = unary_union(polys)
        if geom.geom_type == "Polygon":
            geom = MultiPolygon([geom])
        result[dist_id] = (geom.wkt, pop_per_district.get(dist_id, 0))
    return result


# ---------------------------------------------------------------------------
# Leaf visualization
# ---------------------------------------------------------------------------

def show_leaves(
    leaves: set[str],
    title: str = "H3 Leaf Cells",
) -> None:
    """
    Open a matplotlib window showing every leaf cell as an H3 hexagon,
    coloured by resolution (fine = dark, coarse = light).
    """
    import geopandas as gpd
    import matplotlib.pyplot as plt
    import matplotlib.colors as mcolors
    import pandas as pd

    rows = []
    for cell in leaves:
        boundary = h3.cell_to_boundary(cell)
        coords = [(lon, lat) for lat, lon in boundary]
        res = h3.get_resolution(cell)
        rows.append({"geometry": Polygon(coords), "resolution": res})

    gdf = gpd.GeoDataFrame(rows, crs="EPSG:4326")

    min_res = gdf["resolution"].min()
    max_res = gdf["resolution"].max()

    # Normalize resolution to [0, 1]: finer res → darker colour
    cmap = plt.cm.get_cmap("YlOrRd_r")
    norm = mcolors.Normalize(vmin=min_res, vmax=max(max_res, min_res + 1))
    gdf["colour"] = gdf["resolution"].apply(lambda r: mcolors.to_hex(cmap(norm(r))))

    fig, ax = plt.subplots(figsize=(12, 10))
    gdf.plot(ax=ax, color=gdf["colour"], edgecolor="white", linewidth=0.3, alpha=0.85)

    # Legend: one patch per resolution present
    from matplotlib.patches import Patch
    unique_res = sorted(gdf["resolution"].unique())
    legend_elements = [
        Patch(facecolor=mcolors.to_hex(cmap(norm(r))), edgecolor="grey",
              label=f"res {r}  ({(gdf['resolution'] == r).sum():,} cells)")
        for r in unique_res
    ]
    ax.legend(handles=legend_elements, loc="lower left", fontsize=8, title="Resolution")

    ax.set_title(f"{title}\n{len(leaves):,} leaf cells  "
                 f"(res {min_res}–{max_res})", fontsize=12)
    ax.set_xlabel("Longitude")
    ax.set_ylabel("Latitude")
    ax.set_aspect("equal")
    plt.tight_layout()
    plt.show()
