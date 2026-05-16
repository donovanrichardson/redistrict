"""
Graph construction for redistricting.

Pipeline:
  1. spherical_delaunay_triangles  - 3D convex hull over unit-sphere points
  2. urquhart_edges                - remove longest edge per Delaunay triangle
  3. build_metis_graph             - produce adjacency lists + integer weights for PyMETIS

Edge weight formula (k-medoids cost):
    cost(a, b) = dist_km / (2 * sqrt(pop_a)) + dist_km / (2 * sqrt(pop_b))

PyMETIS weight (higher = harder to cut):
    w_metis(a, b) = round(EDGE_WEIGHT_SCALE / cost(a, b))

For edges that are NOT rook-contiguous (cross water / have no shared boundary),
the cost is multiplied by WATER_PENALTY before inversion, making those edges
1/WATER_PENALTY as heavy and therefore easier for METIS to cut.

Example with EDGE_WEIGHT_SCALE=10000, WATER_PENALTY=3:
  pop_a=4000, pop_b=3600, dist=5 km, land border
    cost  = 5/(2*63.2) + 5/(2*60.0) = 0.0813
    w     = 10000 / 0.0813 = 123,000  (hard to cut)
  Same pair, water border:
    cost  = 0.0813 * 3 = 0.244
    w     = 10000 / 0.244 = 41,000   (3x easier to cut)
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


def build_metis_graph(
    nodes: Sequence[dict],
    edges: set[tuple[int, int]],
    adjacent_geoid_pairs: set[tuple[str, str]],
    water_penalty: float = WATER_PENALTY,
    scale: int = EDGE_WEIGHT_SCALE,
) -> tuple[list[list[int]], list[int], list[int]]:
    """
    Build PyMETIS-ready graph structures from nodes and Urquhart edges.

    Parameters
    ----------
    nodes:
        Ordered list of node dicts (geoid, pop, lat, lon). Index in this list
        is the node index used throughout.
    edges:
        Set of (i, j) pairs with i < j (Urquhart graph).
    adjacent_geoid_pairs:
        Set of (geoid_a, geoid_b) pairs (geoid_a < geoid_b) from rook
        contiguity. Edges NOT in this set receive the water penalty.
    water_penalty:
        Cost multiplier for non-rook-contiguous edges (default 3).
    scale:
        Integer scaling factor for converting float costs to METIS weights.

    Returns
    -------
    adjacency_lists : list[list[int]]
        adjacency_lists[i] = list of neighbor indices of node i.
    eweights : list[int]
        Flat edge-weight list aligned with adjacency_lists.
    nweights : list[int]
        Node weight = population (min 1).
    """
    geoid_index = {n["geoid"]: idx for idx, n in enumerate(nodes)}
    n = len(nodes)

    # Build symmetric adjacency with weights.
    # adjacency_map[i][j] = METIS edge weight
    adjacency_map: list[dict[int, int]] = [{} for _ in range(n)]

    for i, j in edges:
        dist = haversine_km(
            nodes[i]["lat"], nodes[i]["lon"],
            nodes[j]["lat"], nodes[j]["lon"],
        )
        cost = _edge_cost(nodes, i, j, dist_km=dist)

        geoid_pair = (
            min(nodes[i]["geoid"], nodes[j]["geoid"]),
            max(nodes[i]["geoid"], nodes[j]["geoid"]),
        )
        if geoid_pair not in adjacent_geoid_pairs:
            cost *= water_penalty

        w = max(MIN_EDGE_WEIGHT, round(scale / cost))
        adjacency_map[i][j] = w
        adjacency_map[j][i] = w

    adjacency_lists: list[list[int]] = []
    eweights: list[int] = []

    for i in range(n):
        neighbors = sorted(adjacency_map[i].keys())
        adjacency_lists.append(neighbors)
        for nb in neighbors:
            eweights.append(adjacency_map[i][nb])

    nweights = [max(MIN_POP, n["pop"]) for n in nodes]

    return adjacency_lists, eweights, nweights
