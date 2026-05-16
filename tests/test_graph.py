"""Tests for redistrict.graph (pure, no DB)."""

import math

import pytest

from redistrict.graph import (
    EDGE_WEIGHT_SCALE,
    MIN_EDGE_WEIGHT,
    WATER_PENALTY,
    _edge_cost,
    build_metis_graph,
    haversine_km,
    spherical_delaunay_triangles,
    urquhart_edges,
)


# ---------------------------------------------------------------------------
# haversine_km
# ---------------------------------------------------------------------------

class TestHaversine:
    def test_same_point_is_zero(self):
        assert haversine_km(41.0, -71.0, 41.0, -71.0) == pytest.approx(0.0)

    def test_known_distance(self):
        # Providence, RI -> Boston, MA: ~66 km great-circle
        dist = haversine_km(41.824, -71.412, 42.360, -71.058)
        assert 63 < dist < 70

    def test_symmetry(self):
        a = haversine_km(40.0, -75.0, 41.0, -74.0)
        b = haversine_km(41.0, -74.0, 40.0, -75.0)
        assert a == pytest.approx(b)

    def test_equatorial_degree(self):
        # One degree of longitude on the equator ~ 111.3 km
        dist = haversine_km(0.0, 0.0, 0.0, 1.0)
        assert 111.0 < dist < 112.0


# ---------------------------------------------------------------------------
# spherical_delaunay_triangles
# ---------------------------------------------------------------------------

class TestSphericalDelaunay:
    def _rhode_island_nodes(self):
        return [
            {"geoid": "A", "pop": 1000, "lat": 41.70, "lon": -71.55},
            {"geoid": "B", "pop": 2000, "lat": 41.72, "lon": -71.48},
            {"geoid": "C", "pop": 1500, "lat": 41.76, "lon": -71.45},
            {"geoid": "D", "pop": 1200, "lat": 41.80, "lon": -71.50},
            {"geoid": "E", "pop": 900,  "lat": 41.83, "lon": -71.57},
        ]

    def test_returns_triangles(self):
        nodes = self._rhode_island_nodes()
        triangles = spherical_delaunay_triangles(nodes)
        assert len(triangles) > 0
        for tri in triangles:
            assert len(tri) == 3

    def test_indices_in_range(self):
        nodes = self._rhode_island_nodes()
        triangles = spherical_delaunay_triangles(nodes)
        n = len(nodes)
        for tri in triangles:
            for idx in tri:
                assert 0 <= idx < n

    def test_sorted_indices(self):
        nodes = self._rhode_island_nodes()
        triangles = spherical_delaunay_triangles(nodes)
        for tri in triangles:
            assert tri == tuple(sorted(tri))


# ---------------------------------------------------------------------------
# urquhart_edges
# ---------------------------------------------------------------------------

class TestUrquhartEdges:
    def _five_nodes(self):
        return [
            {"geoid": "A", "pop": 1000, "lat": 41.70, "lon": -71.55},
            {"geoid": "B", "pop": 2000, "lat": 41.72, "lon": -71.48},
            {"geoid": "C", "pop": 1500, "lat": 41.76, "lon": -71.45},
            {"geoid": "D", "pop": 1200, "lat": 41.80, "lon": -71.50},
            {"geoid": "E", "pop": 900,  "lat": 41.83, "lon": -71.57},
        ]

    def test_returns_set_of_pairs(self):
        nodes = self._five_nodes()
        triangles = spherical_delaunay_triangles(nodes)
        edges = urquhart_edges(nodes, triangles)
        assert isinstance(edges, set)
        for edge in edges:
            i, j = edge
            assert i < j

    def test_urquhart_subset_of_delaunay(self):
        nodes = self._five_nodes()
        triangles = spherical_delaunay_triangles(nodes)
        urq = urquhart_edges(nodes, triangles)

        all_delaunay = set()
        for tri in triangles:
            i, j, k = tri
            all_delaunay.add((min(i, j), max(i, j)))
            all_delaunay.add((min(j, k), max(j, k)))
            all_delaunay.add((min(i, k), max(i, k)))

        assert urq.issubset(all_delaunay)

    def test_fewer_edges_than_delaunay(self):
        nodes = self._five_nodes()
        triangles = spherical_delaunay_triangles(nodes)
        urq = urquhart_edges(nodes, triangles)

        all_delaunay = set()
        for tri in triangles:
            i, j, k = tri
            all_delaunay.add((min(i, j), max(i, j)))
            all_delaunay.add((min(j, k), max(j, k)))
            all_delaunay.add((min(i, k), max(i, k)))

        assert len(urq) <= len(all_delaunay)


# ---------------------------------------------------------------------------
# _edge_cost
# ---------------------------------------------------------------------------

class TestEdgeCost:
    def _nodes(self):
        return [
            {"geoid": "A", "pop": 4000, "lat": 41.70, "lon": -71.55},
            {"geoid": "B", "pop": 3600, "lat": 41.72, "lon": -71.48},
        ]

    def test_cost_is_positive(self):
        nodes = self._nodes()
        cost = _edge_cost(nodes, 0, 1)
        assert cost > 0

    def test_cost_equals_haversine(self):
        nodes = self._nodes()
        dist = haversine_km(nodes[0]["lat"], nodes[0]["lon"],
                            nodes[1]["lat"], nodes[1]["lon"])
        assert _edge_cost(nodes, 0, 1) == pytest.approx(dist)

    def test_longer_distance_higher_cost(self):
        near = [
            {"geoid": "A", "pop": 1000, "lat": 41.70, "lon": -71.55},
            {"geoid": "B", "pop": 1000, "lat": 41.71, "lon": -71.54},
        ]
        far = [
            {"geoid": "A", "pop": 1000, "lat": 41.70, "lon": -71.55},
            {"geoid": "B", "pop": 1000, "lat": 41.90, "lon": -71.35},
        ]
        assert _edge_cost(near, 0, 1) < _edge_cost(far, 0, 1)


# ---------------------------------------------------------------------------
# build_metis_graph
# ---------------------------------------------------------------------------

class TestBuildMetisGraph:
    def _five_nodes(self):
        return [
            {"geoid": "A", "pop": 4000, "lat": 41.70, "lon": -71.55},
            {"geoid": "B", "pop": 3600, "lat": 41.72, "lon": -71.48},
            {"geoid": "C", "pop": 2000, "lat": 41.76, "lon": -71.45},
            {"geoid": "D", "pop": 1200, "lat": 41.80, "lon": -71.50},
            {"geoid": "E", "pop": 900,  "lat": 41.83, "lon": -71.57},
        ]

    def _edges(self):
        return {(0, 1), (1, 2), (2, 3), (3, 4), (0, 4)}

    def test_output_lengths(self):
        nodes = self._five_nodes()
        edges = self._edges()
        adj, ew, nw = build_metis_graph(nodes, edges, set())
        assert len(adj) == len(nodes)
        assert len(nw) == len(nodes)
        assert len(ew) == sum(len(nb) for nb in adj)

    def test_symmetric_adjacency(self):
        nodes = self._five_nodes()
        edges = self._edges()
        adj, ew, nw = build_metis_graph(nodes, edges, set())
        for i, neighbors in enumerate(adj):
            for j in neighbors:
                assert i in adj[j]

    def test_water_penalty_lowers_weight(self):
        nodes = self._five_nodes()
        edges = {(0, 1)}

        # A-B is adjacent (land border)
        adj_land, ew_land, _ = build_metis_graph(
            nodes, edges, {("A", "B")}
        )
        # A-B is NOT adjacent (water)
        adj_water, ew_water, _ = build_metis_graph(
            nodes, edges, set()
        )
        # Land weight should be ~3x the water weight
        assert ew_land[0] > ew_water[0]
        ratio = ew_land[0] / ew_water[0]
        assert pytest.approx(ratio, rel=0.01) == WATER_PENALTY

    def test_node_weights_are_populations(self):
        nodes = self._five_nodes()
        _, _, nw = build_metis_graph(nodes, self._edges(), set())
        for i, node in enumerate(nodes):
            assert nw[i] == node["pop"]

    def test_minimum_edge_weight_enforced(self):
        nodes = self._five_nodes()
        edges = self._edges()
        adj, ew, _ = build_metis_graph(nodes, edges, set())
        assert all(w >= MIN_EDGE_WEIGHT for w in ew)
