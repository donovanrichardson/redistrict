"""
Graph construction for redistricting.

Pipeline:
  1. spherical_delaunay_triangles  - 3D convex hull over unit-sphere points
  2. urquhart_edges                - remove longest edge per Delaunay triangle
  3. build_metis_graph             - produce adjacency lists + integer weights for PyMETIS

Four edge-weight formulas are available via the `formula` parameter:

  "original" (default)
    cost(a,b) = dist/(2*sqrt(pop_a)) + dist/(2*sqrt(pop_b))
    w = SCALE / cost

  "uniform"
    All land edges receive weight SCALE. Only the water penalty differentiates.

  "original_clamped"
    Adds a constant C to each cost before inversion so that the resulting
    weight range is exactly 4:1 (w_max = 4 * w_min):
      C = (cost_max - 4 * cost_min) / 3
    w = SCALE / (cost + C)

  "blend"
    Normalises two independent signals each to [0, 0.5] and sums them:
      component 1: 1/cost_orig  (original formula inverted)
      component 2: 1/dist_km    (pure inverse distance)
    w = SCALE * (norm(1/cost_orig) + norm(1/dist_km))

For all formulas, non-rook-contiguous (water) edges are divided by
WATER_PENALTY after the formula weight is computed.
"""

import math
from typing import Sequence

import numpy as np
from scipy.spatial import ConvexHull

# Scale factor: raw costs are small floats; this keeps METIS weights in a
# comfortable integer range (~hundreds to ~millions).
EDGE_WEIGHT_SCALE = 10_000

# Non-adjacent (water-boundary) edges get this multiplier applied to their
# cost before inversion, making them cheaper for METIS to cut.
WATER_PENALTY = 5.0

# Minimum METIS weight so no edge has weight zero.
MIN_EDGE_WEIGHT = 1

# Minimum population used in the sqrt denominator to avoid division by zero
# for zero-pop nodes.
MIN_POP = 1


def haversine_km(lat1: float, lon1: float, lat2: float, lon2: float) -> float:
    """Return great-circle distance in kilometres between two WGS-84 points."""
    r = 6_371.0
    phi1, phi2 = math.radians(lat1), math.radians(lat2)
    d_phi = math.radians(lat2 - lat1)
    d_lam = math.radians(lon2 - lon1)
    a = math.sin(d_phi / 2) ** 2 + math.cos(phi1) * math.cos(phi2) * math.sin(d_lam / 2) ** 2
    return 2 * r * math.asin(math.sqrt(a))


def _to_unit_sphere(lat_deg: float, lon_deg: float) -> tuple[float, float, float]:
    lat = math.radians(lat_deg)
    lon = math.radians(lon_deg)
    return (
        math.cos(lat) * math.cos(lon),
        math.cos(lat) * math.sin(lon),
        math.sin(lat),
    )


def spherical_delaunay_triangles(
    nodes: Sequence[dict],
) -> list[tuple[int, int, int]]:
    """
    Compute a spherical Delaunay triangulation via 3D convex hull.

    Each node dict must have 'lat' and 'lon' keys (degrees).
    Returns a list of (i, j, k) index triples referencing nodes.
    Requires at least 4 non-coplanar nodes.
    """
    pts = np.array([_to_unit_sphere(n["lat"], n["lon"]) for n in nodes])
    hull = ConvexHull(pts)
    return [tuple(sorted(simplex)) for simplex in hull.simplices]


def urquhart_edges(
    nodes: Sequence[dict],
    triangles: list[tuple[int, int, int]],
) -> set[tuple[int, int]]:
    """
    Derive the Urquhart graph from a Delaunay triangulation.

    Removes the longest edge (by haversine distance) from each triangle.
    Returns a set of (i, j) pairs with i < j.
    """
    all_edges: set[tuple[int, int]] = set()
    longest_per_triangle: set[tuple[int, int]] = set()

    for tri in triangles:
        i, j, k = tri
        edge_distances = [
            ((min(i, j), max(i, j)),
             haversine_km(nodes[i]["lat"], nodes[i]["lon"],
                          nodes[j]["lat"], nodes[j]["lon"])),
            ((min(j, k), max(j, k)),
             haversine_km(nodes[j]["lat"], nodes[j]["lon"],
                          nodes[k]["lat"], nodes[k]["lon"])),
            ((min(i, k), max(i, k)),
             haversine_km(nodes[i]["lat"], nodes[i]["lon"],
                          nodes[k]["lat"], nodes[k]["lon"])),
        ]
        for edge, _ in edge_distances:
            all_edges.add(edge)
        longest_edge = max(edge_distances, key=lambda x: x[1])[0]
        longest_per_triangle.add(longest_edge)

    return all_edges - longest_per_triangle


def _edge_cost(
    nodes: Sequence[dict],
    i: int,
    j: int,
    dist_km: float | None = None,
) -> float:
    """
    k-medoids cost for an edge between nodes i and j.

    cost = dist_km / (2 * sqrt(pop_i)) + dist_km / (2 * sqrt(pop_j))
    """
    if dist_km is None:
        dist_km = haversine_km(
            nodes[i]["lat"], nodes[i]["lon"],
            nodes[j]["lat"], nodes[j]["lon"],
        )
    pop_i = max(nodes[i]["pop"], MIN_POP)
    pop_j = max(nodes[j]["pop"], MIN_POP)
    return dist_km / (2 * math.sqrt(pop_i)) + dist_km / (2 * math.sqrt(pop_j))


def _nearest_in_set(
    u: int,
    candidates: list[int],
    nodes: Sequence[dict],
) -> tuple[int, float]:
    """
    Return (v, dist_km) for the candidate node closest to u by haversine,
    pre-filtered to ~1% of candidates via a 2-pass bounding-box sort.

    Pass 1: keep the closest sqrt(0.01) fraction by |lat| difference.
    Pass 2: keep the closest sqrt(0.01) fraction of those by |lon| difference.
    Haversine is computed only on the remaining ~1% subset.
    """
    frac = math.sqrt(0.01)
    lat_u = nodes[u]["lat"]
    lon_u = nodes[u]["lon"]

    k_lat = max(1, math.ceil(frac * len(candidates)))
    by_lat = sorted(candidates, key=lambda v: abs(nodes[v]["lat"] - lat_u))[:k_lat]

    k_lon = max(1, math.ceil(frac * len(by_lat)))
    by_lon = sorted(by_lat, key=lambda v: abs(nodes[v]["lon"] - lon_u))[:k_lon]

    best_v, best_d = -1, float("inf")
    for v in by_lon:
        d = haversine_km(lat_u, lon_u, nodes[v]["lat"], nodes[v]["lon"])
        if d < best_d:
            best_d, best_v = d, v
    return best_v, best_d


def reconnect_components(
    nodes: Sequence[dict],
    components: list[list[int]],
) -> set[tuple[int, int]]:
    """
    Return the minimum set of bridging edges needed to make the graph connected.

    For each disconnected component, finds the closest node pair between that
    component and the already-connected nodes and adds one bridge edge. All
    returned edges should be treated as non-adjacent (water) by the caller.
    """
    if len(components) <= 1:
        return set()

    components = sorted(components, key=len, reverse=True)
    connected: list[int] = list(components[0])
    remaining = [list(c) for c in components[1:]]
    new_edges: set[tuple[int, int]] = set()

    while remaining:
        best_dist = float("inf")
        best_edge: tuple[int, int] = (-1, -1)
        best_idx = 0

        for idx, comp in enumerate(remaining):
            for u in comp:
                v, d = _nearest_in_set(u, connected, nodes)
                if d < best_dist:
                    best_dist = d
                    best_edge = (min(u, v), max(u, v))
                    best_idx = idx

        new_edges.add(best_edge)
        connected.extend(remaining[best_idx])
        remaining.pop(best_idx)

    return new_edges


def check_connectivity(
    nodes: Sequence[dict],
    edges: set[tuple[int, int]],
) -> list[list[int]]:
    """
    Return connected components as lists of node indices.

    A single-element list means the graph is fully connected.
    Multiple elements means the graph is disconnected — METIS contig will fail.
    """
    n = len(nodes)
    adj: list[list[int]] = [[] for _ in range(n)]
    for i, j in edges:
        adj[i].append(j)
        adj[j].append(i)

    visited = [False] * n
    components: list[list[int]] = []
    for start in range(n):
        if visited[start]:
            continue
        component: list[int] = []
        stack = [start]
        while stack:
            v = stack.pop()
            if visited[v]:
                continue
            visited[v] = True
            component.append(v)
            stack.extend(adj[v])
        components.append(component)
    return components


def build_metis_graph(
    nodes: Sequence[dict],
    edges: set[tuple[int, int]],
    adjacent_geoid_pairs: set[tuple[str, str]],
    water_penalty: float = WATER_PENALTY,
    scale: int = EDGE_WEIGHT_SCALE,
    formula: str = "original",
) -> tuple[list[list[int]], list[int], list[int]]:
    """
    Build PyMETIS-ready graph structures from nodes and Urquhart edges.

    Parameters
    ----------
    nodes:
        Ordered list of node dicts (geoid, pop, lat, lon).
    edges:
        Set of (i, j) pairs with i < j (Urquhart graph).
    adjacent_geoid_pairs:
        Set of (geoid_a, geoid_b) pairs (geoid_a < geoid_b) from rook
        contiguity. Edges NOT in this set are divided by water_penalty.
    water_penalty:
        Divisor applied to non-rook-contiguous edge weights.
    scale:
        Integer scaling factor.
    formula:
        One of "original", "uniform", "original_clamped", "blend".
        See module docstring for details.

    Returns
    -------
    adjacency_lists : list[list[int]]
    eweights : list[int]
    nweights : list[int]
    """
    if formula not in ("original", "uniform", "original_clamped", "blend"):
        raise ValueError(f"Unknown formula {formula!r}")

    n = len(nodes)
    edge_list = sorted(edges)

    # Compute distances once per edge.
    dists = {
        (i, j): haversine_km(
            nodes[i]["lat"], nodes[i]["lon"],
            nodes[j]["lat"], nodes[j]["lon"],
        )
        for i, j in edge_list
    }

    # Compute raw weights (before water penalty).
    if formula == "uniform":
        raw: dict[tuple[int, int], float] = {e: float(scale) for e in edge_list}

    elif formula == "original":
        raw = {
            e: scale / _edge_cost(nodes, e[0], e[1], dist_km=dists[e])
            for e in edge_list
        }

    elif formula == "original_clamped":
        costs = {e: _edge_cost(nodes, e[0], e[1], dist_km=dists[e]) for e in edge_list}
        if costs:
            cmin = min(costs.values())
            cmax = max(costs.values())
            # C chosen so w_max / w_min == 4, i.e. (cost_max - 4*cost_min) / 3.
            c_offset = max(0.0, (cmax - 4.0 * cmin) / 3.0)
        else:
            c_offset = 0.0
        raw = {e: scale / (costs[e] + c_offset) for e in edge_list}

    else:  # "blend"
        costs = {e: _edge_cost(nodes, e[0], e[1], dist_km=dists[e]) for e in edge_list}
        w_orig = {e: 1.0 / c for e, c in costs.items()}
        w_inv  = {e: 1.0 / dists[e] for e in edge_list}

        wo_min, wo_max = min(w_orig.values()), max(w_orig.values())
        wi_min, wi_max = min(w_inv.values()),  max(w_inv.values())

        raw = {}
        for e in edge_list:
            n1 = 0.5 * (w_orig[e] - wo_min) / (wo_max - wo_min) if wo_max > wo_min else 0.25
            n2 = 0.5 * (w_inv[e]  - wi_min) / (wi_max - wi_min) if wi_max > wi_min else 0.25
            raw[e] = scale * (n1 + n2)

    # Apply water penalty and build adjacency map.
    adjacency_map: list[dict[int, int]] = [{} for _ in range(n)]
    for i, j in edge_list:
        geoid_pair = (
            min(nodes[i]["geoid"], nodes[j]["geoid"]),
            max(nodes[i]["geoid"], nodes[j]["geoid"]),
        )
        w = raw[(i, j)]
        if geoid_pair not in adjacent_geoid_pairs:
            w /= water_penalty
        w = max(MIN_EDGE_WEIGHT, round(w))
        adjacency_map[i][j] = w
        adjacency_map[j][i] = w

    adjacency_lists: list[list[int]] = []
    eweights: list[int] = []
    for idx in range(n):
        neighbors = sorted(adjacency_map[idx].keys())
        adjacency_lists.append(neighbors)
        for nb in neighbors:
            eweights.append(adjacency_map[idx][nb])

    nweights = [max(MIN_POP, node["pop"]) for node in nodes]
    return adjacency_lists, eweights, nweights
