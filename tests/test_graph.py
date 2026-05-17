"""Tests for redistrict.graph (pure, no DB)."""

import math

import pytest

from redistrict.graph import (
    EDGE_WEIGHT_SCALE,
    MIN_EDGE_WEIGHT,
    WATER_PENALTY,
    _edge_cost,
    build_metis_graph,
    check_connectivity,
    haversine_km,
    reconnect_components,
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

    def test_higher_pop_lower_cost(self):
        low_pop = [
            {"geoid": "A", "pop": 100, "lat": 41.70, "lon": -71.55},
            {"geoid": "B", "pop": 100, "lat": 41.72, "lon": -71.48},
        ]
        high_pop = [
            {"geoid": "A", "pop": 10000, "lat": 41.70, "lon": -71.55},
            {"geoid": "B", "pop": 10000, "lat": 41.72, "lon": -71.48},
        ]
        assert _edge_cost(low_pop, 0, 1) > _edge_cost(high_pop, 0, 1)

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


# ---------------------------------------------------------------------------
# check_connectivity
# ---------------------------------------------------------------------------

class TestCheckConnectivity:
    def _nodes(self, n: int):
        return [{"geoid": str(i), "pop": 100, "lat": 41.0 + i * 0.01, "lon": -71.0}
                for i in range(n)]

    def test_fully_connected_single_component(self):
        nodes = self._nodes(4)
        edges = {(0, 1), (1, 2), (2, 3)}
        components = check_connectivity(nodes, edges)
        assert len(components) == 1
        assert len(components[0]) == 4

    def test_two_components(self):
        nodes = self._nodes(4)
        edges = {(0, 1), (2, 3)}  # two disconnected pairs
        components = check_connectivity(nodes, edges)
        assert len(components) == 2
        sizes = sorted(len(c) for c in components)
        assert sizes == [2, 2]

    def test_isolated_node_is_its_own_component(self):
        nodes = self._nodes(3)
        edges = {(0, 1)}  # node 2 is isolated
        components = check_connectivity(nodes, edges)
        assert len(components) == 2

    def test_complete_graph_single_component(self):
        nodes = self._nodes(4)
        edges = {(0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)}
        components = check_connectivity(nodes, edges)
        assert len(components) == 1

    def test_empty_edges_all_isolated(self):
        nodes = self._nodes(3)
        components = check_connectivity(nodes, set())
        assert len(components) == 3


# ---------------------------------------------------------------------------
# reconnect_components
# ---------------------------------------------------------------------------

class TestReconnectComponents:
    def _geo_nodes(self):
        # Spread across RI for realistic haversine distances
        return [
            {"geoid": "A", "pop": 100, "lat": 41.70, "lon": -71.55},
            {"geoid": "B", "pop": 100, "lat": 41.72, "lon": -71.48},
            {"geoid": "C", "pop": 100, "lat": 41.76, "lon": -71.45},
            {"geoid": "D", "pop": 100, "lat": 41.80, "lon": -71.50},
            {"geoid": "E", "pop": 100, "lat": 41.83, "lon": -71.57},
        ]

    def test_already_connected_returns_empty(self):
        nodes = self._geo_nodes()
        components = [[0, 1, 2, 3, 4]]
        bridges = reconnect_components(nodes, components)
        assert bridges == set()

    def test_two_components_returns_one_bridge(self):
        nodes = self._geo_nodes()
        components = [[0, 1, 2], [3, 4]]
        bridges = reconnect_components(nodes, components)
        assert len(bridges) == 1

    def test_bridge_connects_components(self):
        nodes = self._geo_nodes()
        components = [[0, 1, 2], [3, 4]]
        bridges = reconnect_components(nodes, components)
        edge = next(iter(bridges))
        # One node from each component
        i, j = edge
        assert (i in {0, 1, 2}) != (j in {0, 1, 2})

    def test_three_components_returns_two_bridges(self):
        nodes = self._geo_nodes()
        components = [[0, 1], [2, 3], [4]]
        bridges = reconnect_components(nodes, components)
        assert len(bridges) == 2

    def test_result_makes_graph_connected(self):
        nodes = self._geo_nodes()
        components = [[0, 1], [2], [3, 4]]
        bridges = reconnect_components(nodes, components)
        all_edges = bridges.copy()
        # Add original intra-component edges
        all_edges |= {(0, 1), (3, 4)}
        final_components = check_connectivity(nodes, all_edges)
        assert len(final_components) == 1
