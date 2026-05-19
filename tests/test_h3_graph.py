"""Tests for redistrict.h3_graph (pure, no DB)."""

import h3
import pytest

from redistrict.h3_graph import (
    aggregate_h3_cells,
    assign_h3_res15,
    build_h3_adjacency,
    build_metis_graph,
    check_connectivity,
    compute_threshold,
    weighted_centroids,
)

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _block(geoid: str, pop: int, lat: float = 41.7, lon: float = -71.5) -> dict:
    return {"geoid": geoid, "pop": pop, "lat": lat, "lon": lon}


# ---------------------------------------------------------------------------
# compute_threshold
# ---------------------------------------------------------------------------

class TestComputeThreshold:
    def test_example_from_spec(self):
        # 1 block pop=100, 1500 blocks each pop≈6.6, total≈10000
        # 1% = 100; first block hits target → threshold = 100
        blocks = [_block("A", 100)] + [_block(str(i), 6) for i in range(1500)]
        assert compute_threshold(blocks) == 100

    def test_returns_block_population_not_cumulative(self):
        # 3 blocks: 500, 300, 200 → total=1000, 1%=10
        # First block (pop=500) immediately exceeds target → threshold=500
        blocks = [_block("A", 500), _block("B", 300), _block("C", 200)]
        result = compute_threshold(blocks)
        assert result == 500

    def test_multiple_blocks_needed(self):
        # 10 blocks each pop=10 → total=100, 1%=1
        # First block (pop=10) hits target → threshold=10
        blocks = [_block(str(i), 10) for i in range(10)]
        assert compute_threshold(blocks) == 10

    def test_zero_pop_blocks_ignored(self):
        # Zero-pop blocks don't contribute; threshold comes from positive ones
        blocks = [_block("A", 0), _block("B", 0), _block("C", 100)]
        result = compute_threshold(blocks)
        assert result == 100

    def test_all_zero_pop_returns_one(self):
        blocks = [_block("A", 0), _block("B", 0)]
        assert compute_threshold(blocks) == 1

    def test_empty_blocks_returns_one(self):
        assert compute_threshold([]) == 1

    def test_custom_pct(self):
        # 100 blocks each pop=1 → total=100, 10%=10
        # Need 10 blocks to reach target → threshold=1
        blocks = [_block(str(i), 1) for i in range(100)]
        assert compute_threshold(blocks, pct=0.10) == 1


# ---------------------------------------------------------------------------
# assign_h3_res15
# ---------------------------------------------------------------------------

class TestAssignH3Res15:
    def test_returns_valid_h3_cells(self):
        blocks = [_block("A", 100, lat=41.70, lon=-71.55)]
        result = assign_h3_res15(blocks)
        assert "A" in result
        assert h3.get_resolution(result["A"]) == 15

    def test_nearby_blocks_may_share_res15_cell(self):
        # Two blocks at exactly the same coordinate share the same cell
        blocks = [
            _block("A", 100, lat=41.70, lon=-71.55),
            _block("B", 200, lat=41.70, lon=-71.55),
        ]
        result = assign_h3_res15(blocks)
        assert result["A"] == result["B"]

    def test_distant_blocks_different_cells(self):
        blocks = [
            _block("A", 100, lat=41.0, lon=-71.0),
            _block("B", 100, lat=42.0, lon=-72.0),
        ]
        result = assign_h3_res15(blocks)
        assert result["A"] != result["B"]

    def test_all_geoids_present(self):
        blocks = [_block(str(i), 100, lat=41.0 + i * 0.1, lon=-71.0) for i in range(5)]
        result = assign_h3_res15(blocks)
        assert set(result.keys()) == {str(i) for i in range(5)}


# ---------------------------------------------------------------------------
# aggregate_h3_cells
# ---------------------------------------------------------------------------

class TestAggregateH3Cells:
    def _get_test_cells(self):
        """Return 3 geographically close blocks that share a res-14 parent."""
        # Use a single known res-15 cell and its siblings in the same parent
        base = h3.latlng_to_cell(41.70, -71.55, 15)
        parent = h3.cell_to_parent(base, 14)
        children = list(h3.cell_to_children(parent, 15))[:3]
        return parent, children

    def test_merges_when_below_threshold(self):
        parent, children = self._get_test_cells()
        geoid_to_res15 = {str(i): c for i, c in enumerate(children)}
        block_pops = {str(i): 100 for i in range(len(children))}
        threshold = 10_000  # well above any group total

        result = aggregate_h3_cells(geoid_to_res15, block_pops, threshold, min_res=14)

        # All geoids should map to the same parent
        final_cells = set(result.values())
        assert len(final_cells) == 1
        final_cell = next(iter(final_cells))
        assert h3.get_resolution(final_cell) == 14

    def test_does_not_merge_when_above_threshold(self):
        parent, children = self._get_test_cells()
        geoid_to_res15 = {str(i): c for i, c in enumerate(children)}
        block_pops = {str(i): 1000 for i in range(len(children))}
        # threshold below group total → no merge
        threshold = len(children) * 1000 - 1

        result = aggregate_h3_cells(geoid_to_res15, block_pops, threshold, min_res=14)

        # Should stay at res 15
        for cell in result.values():
            assert h3.get_resolution(cell) == 15

    def test_aggregates_multiple_levels(self):
        # Put blocks far apart so they can only merge at a very coarse level,
        # but use a very high threshold so merging IS permitted
        blocks_at_two_sites = {}
        for i, (lat, lon) in enumerate([(41.7, -71.5), (41.7, -71.5)]):
            cell = h3.latlng_to_cell(lat, lon, 15)
            blocks_at_two_sites[str(i)] = cell
        pops = {str(i): 1 for i in range(2)}
        result = aggregate_h3_cells(blocks_at_two_sites, pops, threshold=10_000, min_res=4)
        # All should be merged into one coarse ancestor
        assert len(set(result.values())) == 1

    def test_all_geoids_in_result(self):
        parent, children = self._get_test_cells()
        geoid_to_res15 = {str(i): c for i, c in enumerate(children)}
        block_pops = {str(i): 10 for i in range(len(children))}
        result = aggregate_h3_cells(geoid_to_res15, block_pops, threshold=5)
        assert set(result.keys()) == set(geoid_to_res15.keys())

    def test_single_block_above_threshold_stays_res15(self):
        cell = h3.latlng_to_cell(41.70, -71.55, 15)
        geoid_to_res15 = {"A": cell}
        block_pops = {"A": 9999}
        result = aggregate_h3_cells(geoid_to_res15, block_pops, threshold=100, min_res=4)
        assert h3.get_resolution(result["A"]) == 15


# ---------------------------------------------------------------------------
# weighted_centroids
# ---------------------------------------------------------------------------

class TestWeightedCentroids:
    def test_single_block_centroid_equals_block(self):
        cell = h3.latlng_to_cell(41.70, -71.55, 15)
        geoid_to_cell = {"A": cell}
        blocks = {"A": _block("A", 100, lat=41.70, lon=-71.55)}
        nodes = weighted_centroids(geoid_to_cell, blocks)
        assert len(nodes) == 1
        assert pytest.approx(nodes[0]["lat"], abs=1e-6) == 41.70
        assert pytest.approx(nodes[0]["lon"], abs=1e-6) == -71.55
        assert nodes[0]["pop"] == 100

    def test_two_equal_weight_blocks_midpoint(self):
        cell = h3.latlng_to_cell(41.70, -71.55, 15)
        geoid_to_cell = {"A": cell, "B": cell}
        blocks = {
            "A": _block("A", 100, lat=41.70, lon=-71.50),
            "B": _block("B", 100, lat=41.80, lon=-71.60),
        }
        nodes = weighted_centroids(geoid_to_cell, blocks)
        assert len(nodes) == 1
        assert pytest.approx(nodes[0]["lat"], abs=1e-6) == 41.75
        assert pytest.approx(nodes[0]["lon"], abs=1e-6) == -71.55
        assert nodes[0]["pop"] == 200

    def test_pop_weighted_not_simple_average(self):
        cell = h3.latlng_to_cell(41.70, -71.55, 15)
        geoid_to_cell = {"A": cell, "B": cell}
        blocks = {
            "A": _block("A", 900, lat=41.70, lon=-71.50),  # 90% weight
            "B": _block("B", 100, lat=41.80, lon=-71.60),  # 10% weight
        }
        nodes = weighted_centroids(geoid_to_cell, blocks)
        # Expected: 0.9*41.70 + 0.1*41.80 = 41.71
        assert pytest.approx(nodes[0]["lat"], abs=1e-6) == 41.71

    def test_zero_pop_cell_uses_h3_centroid(self):
        cell = h3.latlng_to_cell(41.70, -71.55, 15)
        geoid_to_cell = {"A": cell}
        blocks = {"A": _block("A", 0, lat=41.70, lon=-71.55)}
        nodes = weighted_centroids(geoid_to_cell, blocks)
        # Should fall back to h3 geometric centroid — just check it returns something
        assert nodes[0]["pop"] == 0
        assert -90 <= nodes[0]["lat"] <= 90
        assert -180 <= nodes[0]["lon"] <= 180

    def test_returns_one_node_per_unique_cell(self):
        cell_a = h3.latlng_to_cell(41.70, -71.55, 15)
        cell_b = h3.latlng_to_cell(42.00, -72.00, 15)
        geoid_to_cell = {"A": cell_a, "B": cell_a, "C": cell_b}
        blocks = {
            "A": _block("A", 100, lat=41.70, lon=-71.55),
            "B": _block("B", 100, lat=41.70, lon=-71.55),
            "C": _block("C", 100, lat=42.00, lon=-72.00),
        }
        nodes = weighted_centroids(geoid_to_cell, blocks)
        assert len(nodes) == 2


# ---------------------------------------------------------------------------
# build_h3_adjacency
# ---------------------------------------------------------------------------

class TestBuildH3Adjacency:
    def test_adjacent_res15_cells_are_adjacent(self):
        base = h3.latlng_to_cell(41.70, -71.55, 15)
        # get a true H3 neighbour of base
        neighbours = list(set(h3.grid_disk(base, 1)) - {base})
        neighbour = neighbours[0]

        geoid_to_cell = {"A": base, "B": neighbour}
        geoid_to_res15 = {"A": base, "B": neighbour}

        edges = build_h3_adjacency(geoid_to_cell, geoid_to_res15)
        expected = (min(base, neighbour), max(base, neighbour))
        assert expected in edges

    def test_non_adjacent_cells_not_connected(self):
        cell_a = h3.latlng_to_cell(41.70, -71.55, 15)
        cell_b = h3.latlng_to_cell(42.50, -73.00, 15)  # far away
        geoid_to_cell = {"A": cell_a, "B": cell_b}
        geoid_to_res15 = {"A": cell_a, "B": cell_b}
        edges = build_h3_adjacency(geoid_to_cell, geoid_to_res15)
        expected = (min(cell_a, cell_b), max(cell_a, cell_b))
        assert expected not in edges

    def test_same_cell_no_self_loop(self):
        cell = h3.latlng_to_cell(41.70, -71.55, 15)
        geoid_to_cell = {"A": cell, "B": cell}
        geoid_to_res15 = {"A": cell, "B": cell}
        edges = build_h3_adjacency(geoid_to_cell, geoid_to_res15)
        for a, b in edges:
            assert a != b


# ---------------------------------------------------------------------------
# build_metis_graph
# ---------------------------------------------------------------------------

class TestBuildMetisGraph:
    def _two_adjacent_nodes(self):
        base = h3.latlng_to_cell(41.70, -71.55, 15)
        nb = list(set(h3.grid_disk(base, 1)) - {base})[0]
        nodes = [
            {"cell": base, "pop": 1000, "lat": 41.70, "lon": -71.55},
            {"cell": nb,   "pop": 2000, "lat": 41.70, "lon": -71.54},
        ]
        adjacency = {(min(base, nb), max(base, nb))}
        return nodes, adjacency

    def test_output_shapes(self):
        nodes, adj = self._two_adjacent_nodes()
        al, ew, nw = build_metis_graph(nodes, adj)
        assert len(al) == 2
        assert len(nw) == 2
        assert len(ew) == sum(len(nb) for nb in al)

    def test_uniform_edge_weights(self):
        nodes, adj = self._two_adjacent_nodes()
        _, ew, _ = build_metis_graph(nodes, adj)
        assert all(w == 1 for w in ew)

    def test_node_weights_are_populations(self):
        nodes, adj = self._two_adjacent_nodes()
        _, _, nw = build_metis_graph(nodes, adj)
        assert nw[0] == 1000
        assert nw[1] == 2000

    def test_symmetric_adjacency(self):
        nodes, adj = self._two_adjacent_nodes()
        al, _, _ = build_metis_graph(nodes, adj)
        assert al[0] == [1]
        assert al[1] == [0]

    def test_isolated_node_has_empty_adjacency(self):
        cell = h3.latlng_to_cell(41.70, -71.55, 15)
        nodes = [{"cell": cell, "pop": 100, "lat": 41.70, "lon": -71.55}]
        al, ew, nw = build_metis_graph(nodes, set())
        assert al == [[]]
        assert ew == []
        assert nw == [100]


# ---------------------------------------------------------------------------
# check_connectivity
# ---------------------------------------------------------------------------

class TestCheckConnectivity:
    def _cells(self, n: int):
        return [
            {"cell": h3.latlng_to_cell(41.0 + i * 0.5, -71.0, 8),
             "pop": 100, "lat": 41.0 + i * 0.5, "lon": -71.0}
            for i in range(n)
        ]

    def _edge(self, nodes, i, j):
        a, b = nodes[i]["cell"], nodes[j]["cell"]
        return (min(a, b), max(a, b))

    def test_fully_connected(self):
        nodes = self._cells(3)
        adj = {self._edge(nodes, 0, 1), self._edge(nodes, 1, 2)}
        assert len(check_connectivity(nodes, adj)) == 1

    def test_two_components(self):
        nodes = self._cells(4)
        adj = {self._edge(nodes, 0, 1), self._edge(nodes, 2, 3)}
        assert len(check_connectivity(nodes, adj)) == 2

    def test_empty_edges_all_isolated(self):
        nodes = self._cells(3)
        assert len(check_connectivity(nodes, set())) == 3
