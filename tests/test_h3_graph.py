"""Tests for redistrict.h3_graph (pure, no DB)."""

import h3
import pytest

from redistrict.h3_graph import (
    aggregate_h3_cells,
    assign_geoids_to_leaves,
    assign_h3_res15,
    build_cell_hierarchy,
    build_h3_adjacency,
    build_metis_graph,
    check_connectivity,
    compute_leaves_and_provisional_edges,
    compute_threshold,
    compute_top_level,
    find_zero_pop_leaves,
    resolve_edges,
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
        # 1 block pop=100, 1500 blocks each pop≈6, total≈9100
        # 1% target ≈ 91; first block (pop=100) hits target → threshold=100
        blocks = [_block("A", 100)] + [_block(str(i), 6) for i in range(1500)]
        assert compute_threshold(blocks) == 100

    def test_returns_block_population_not_cumulative(self):
        # 3 blocks: 5000, 3000, 2000 → total=10000, 1%=100
        # First block (pop=5000) immediately exceeds target → threshold=5000
        blocks = [_block("A", 5000), _block("B", 3000), _block("C", 2000)]
        result = compute_threshold(blocks)
        assert result == 5000

    def test_multiple_blocks_needed(self):
        # 200 blocks each pop=50 → total=10000, 1%=100
        # Need 2 blocks to accumulate 100 → threshold=50
        blocks = [_block(str(i), 50) for i in range(200)]
        assert compute_threshold(blocks) == 50

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
    def test_adjacent_cells_are_adjacent(self):
        # Use a coarser resolution so neighbours are more meaningful
        base = h3.latlng_to_cell(41.70, -71.55, 8)
        neighbour = list(set(h3.grid_disk(base, 1)) - {base})[0]

        geoid_to_cell = {"A": base, "B": neighbour}
        edges = build_h3_adjacency(geoid_to_cell)
        expected = (min(base, neighbour), max(base, neighbour))
        assert expected in edges

    def test_non_adjacent_cells_not_connected(self):
        cell_a = h3.latlng_to_cell(41.70, -71.55, 8)
        cell_b = h3.latlng_to_cell(42.50, -73.00, 8)  # far away
        geoid_to_cell = {"A": cell_a, "B": cell_b}
        edges = build_h3_adjacency(geoid_to_cell)
        expected = (min(cell_a, cell_b), max(cell_a, cell_b))
        assert expected not in edges

    def test_same_cell_no_self_loop(self):
        cell = h3.latlng_to_cell(41.70, -71.55, 8)
        geoid_to_cell = {"A": cell, "B": cell}
        edges = build_h3_adjacency(geoid_to_cell)
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


# ---------------------------------------------------------------------------
# build_cell_hierarchy
# ---------------------------------------------------------------------------

class TestBuildCellHierarchy:
    def _res15_siblings(self):
        """Return three res-15 cells sharing the same res-14 parent."""
        base = h3.latlng_to_cell(41.70, -71.55, 15)
        parent14 = h3.cell_to_parent(base, 14)
        children = list(h3.cell_to_children(parent14, 15))[:3]
        return parent14, children

    def test_cell_pop_sums_blocks(self):
        parent14, children = self._res15_siblings()
        geoid_to_res15 = {str(i): c for i, c in enumerate(children)}
        block_pops = {str(i): (i + 1) * 10 for i in range(len(children))}
        cell_pop, _, _ = build_cell_hierarchy(geoid_to_res15, block_pops)
        for i, c in enumerate(children):
            assert cell_pop[c] == (i + 1) * 10
        assert cell_pop[parent14] == sum((i + 1) * 10 for i in range(len(children)))

    def test_zero_pop_blocks_excluded(self):
        parent14, children = self._res15_siblings()
        geoid_to_res15 = {"A": children[0], "B": children[1]}
        block_pops = {"A": 0, "B": 50}
        cell_pop, _, populated_at_res = build_cell_hierarchy(geoid_to_res15, block_pops)
        assert children[0] not in cell_pop
        assert children[1] in cell_pop
        assert children[0] not in populated_at_res[15]

    def test_cell_direct_children(self):
        parent14, children = self._res15_siblings()
        geoid_to_res15 = {str(i): c for i, c in enumerate(children)}
        block_pops = {str(i): 10 for i in range(len(children))}
        _, cell_direct_children, _ = build_cell_hierarchy(geoid_to_res15, block_pops)
        assert set(children) == cell_direct_children[parent14]

    def test_populated_at_res_contains_only_populated(self):
        _, children = self._res15_siblings()
        geoid_to_res15 = {"A": children[0]}
        block_pops = {"A": 100}
        _, _, populated_at_res = build_cell_hierarchy(geoid_to_res15, block_pops)
        assert children[0] in populated_at_res[15]
        # res-0 should contain exactly one cell (the ancestor chain)
        assert len(populated_at_res[0]) == 1

    def test_pop_accumulates_through_all_levels(self):
        _, children = self._res15_siblings()
        geoid_to_res15 = {"A": children[0]}
        block_pops = {"A": 777}
        cell_pop, _, populated_at_res = build_cell_hierarchy(geoid_to_res15, block_pops)
        for res in range(0, 16):
            for cell in populated_at_res.get(res, set()):
                assert cell_pop[cell] == 777


# ---------------------------------------------------------------------------
# compute_leaves_and_provisional_edges
# ---------------------------------------------------------------------------

class TestComputeLeavesAndProvisionalEdges:
    def _hierarchy_two_siblings(self, pop_each: int = 100):
        """Two res-15 cells sharing a res-14 parent."""
        base = h3.latlng_to_cell(41.70, -71.55, 15)
        parent14 = h3.cell_to_parent(base, 14)
        children = list(h3.cell_to_children(parent14, 15))[:2]
        geoid_to_res15 = {str(i): c for i, c in enumerate(children)}
        block_pops = {str(i): pop_each for i in range(len(children))}
        return build_cell_hierarchy(geoid_to_res15, block_pops)

    def test_children_become_leaves_when_parent_exceeds_threshold(self):
        parent14 = h3.cell_to_parent(h3.latlng_to_cell(41.70, -71.55, 15), 14)
        children = list(h3.cell_to_children(parent14, 15))[:2]
        geoid_to_res15 = {str(i): c for i, c in enumerate(children)}
        block_pops = {str(i): 1000 for i in range(2)}
        cell_pop, cell_direct_children, populated_at_res = build_cell_hierarchy(
            geoid_to_res15, block_pops
        )
        # threshold = 1 so parent pop (2000) > threshold → children become leaves
        leaves, _ = compute_leaves_and_provisional_edges(
            cell_pop, cell_direct_children, populated_at_res, threshold=1
        )
        assert set(children).issubset(leaves)

    def test_single_child_does_not_split(self):
        cell15 = h3.latlng_to_cell(41.70, -71.55, 15)
        geoid_to_res15 = {"A": cell15}
        block_pops = {"A": 9999}
        cell_pop, cell_direct_children, populated_at_res = build_cell_hierarchy(
            geoid_to_res15, block_pops
        )
        leaves, _ = compute_leaves_and_provisional_edges(
            cell_pop, cell_direct_children, populated_at_res, threshold=1
        )
        # Only one child at every level → no splits, cleanup marks the res-0 ancestor
        assert len(leaves) == 1
        leaf = next(iter(leaves))
        assert h3.get_resolution(leaf) == 0

    def test_below_threshold_collapses_to_coarse_leaf(self):
        parent14 = h3.cell_to_parent(h3.latlng_to_cell(41.70, -71.55, 15), 14)
        children = list(h3.cell_to_children(parent14, 15))[:3]
        geoid_to_res15 = {str(i): c for i, c in enumerate(children)}
        block_pops = {str(i): 1 for i in range(3)}
        cell_pop, cell_direct_children, populated_at_res = build_cell_hierarchy(
            geoid_to_res15, block_pops
        )
        # threshold >> total pop → never splits → collapses to single res-0 leaf
        leaves, _ = compute_leaves_and_provisional_edges(
            cell_pop, cell_direct_children, populated_at_res, threshold=10_000
        )
        assert len(leaves) == 1
        assert h3.get_resolution(next(iter(leaves))) == 0

    def test_orphaned_branch_gets_coarsest_leaf(self):
        # res-13 parent P13 has two res-14 children: A14 (dense, splits) and B14 (sparse).
        # A14's blocks create fine leaves → P13 enters has_leaf_descendant.
        # B14 has no leaf descendants → would be orphaned by old res-0-only cleanup.
        # New boundary pass should mark B14 (or its res-13 ancestor) as a leaf.
        p13 = h3.latlng_to_cell(41.70, -71.55, 13)
        children14 = sorted(h3.cell_to_children(p13, 14))
        a14, b14 = children14[0], children14[-1]
        # Dense cluster under a14: 3 res-15 children with high pop
        a15s = list(h3.cell_to_children(a14, 15))[:3]
        # Sparse cluster under b14: 2 res-15 children with very low pop
        b15s = list(h3.cell_to_children(b14, 15))[:2]
        geoid_to_res15 = {str(i): c for i, c in enumerate(a15s + b15s)}
        block_pops = {str(i): 5000 for i in range(len(a15s))}  # dense
        block_pops.update({str(len(a15s) + i): 1 for i in range(len(b15s))})  # sparse
        cell_pop, cdc, par = build_cell_hierarchy(geoid_to_res15, block_pops)
        # threshold=100 → a14 (pop≈15000) splits; b14 (pop=2) doesn't
        leaves, _ = compute_leaves_and_provisional_edges(cell_pop, cdc, par, threshold=100)
        # All b15s blocks must have a leaf ancestor
        assigned = assign_geoids_to_leaves(geoid_to_res15, leaves)
        for i in range(len(a15s), len(a15s) + len(b15s)):
            assert str(i) in assigned, f"Block {i} (sparse branch) has no leaf"

    def test_provisional_edges_added_between_neighboring_leaves(self):
        # Build two adjacent res-14 groups (each with 2 res-15 children).
        # threshold=1 so both parents exceed threshold → all children become leaves.
        # Neighboring leaves should have provisional edges between them.
        parent14 = h3.cell_to_parent(h3.latlng_to_cell(41.70, -71.55, 15), 14)
        children14 = list(h3.cell_to_children(parent14, 15))[:2]
        nb14 = list(set(h3.grid_disk(parent14, 1)) - {parent14})[0]
        children_nb14 = list(h3.cell_to_children(nb14, 15))[:2]

        all_cells = children14 + children_nb14
        geoid_to_res15 = {str(i): c for i, c in enumerate(all_cells)}
        block_pops = {str(i): 1000 for i in range(len(all_cells))}
        cell_pop, cell_direct_children, populated_at_res = build_cell_hierarchy(
            geoid_to_res15, block_pops
        )
        leaves, provisional = compute_leaves_and_provisional_edges(
            cell_pop, cell_direct_children, populated_at_res, threshold=1
        )
        assert len(provisional) > 0


# ---------------------------------------------------------------------------
# resolve_edges
# ---------------------------------------------------------------------------

class TestResolveEdges:
    def _two_leaf_cells(self):
        base = h3.latlng_to_cell(41.70, -71.55, 8)
        nb = list(set(h3.grid_disk(base, 1)) - {base})[0]
        return base, nb

    def test_confirms_edge_when_both_are_leaves(self):
        A, N = self._two_leaf_cells()
        leaves = {A, N}
        provisional = {(A, N)}  # directed: A is the originating leaf
        confirmed = resolve_edges(provisional, leaves)
        assert (min(A, N), max(A, N)) in confirmed

    def test_redirects_to_leaf_ancestor(self):
        # N_fine is a res-15 neighbor; its res-8 ancestor N_coarse is a leaf.
        # Edge should be redirected from (A, N_fine) to (A, N_coarse).
        N_fine = h3.latlng_to_cell(41.70, -71.55, 15)
        N_coarse = h3.cell_to_parent(N_fine, 8)
        A = h3.latlng_to_cell(42.00, -72.00, 15)
        leaves = {A, N_coarse}
        provisional = {(A, N_fine)}  # directed: A is leaf, N_fine is the raw neighbor
        confirmed = resolve_edges(provisional, leaves)
        expected = (min(A, N_coarse), max(A, N_coarse))
        assert expected in confirmed

    def test_discards_when_no_leaf_ancestor(self):
        A = h3.latlng_to_cell(41.70, -71.55, 15)
        N = h3.latlng_to_cell(42.00, -72.00, 15)
        leaves = {A}  # N has no leaf ancestor
        provisional = {(A, N)}
        confirmed = resolve_edges(provisional, leaves)
        assert len(confirmed) == 0

    def test_no_self_loops(self):
        A = h3.latlng_to_cell(41.70, -71.55, 8)
        leaves = {A}
        provisional = {(A, A)}
        confirmed = resolve_edges(provisional, leaves)
        assert (A, A) not in confirmed


# ---------------------------------------------------------------------------
# assign_geoids_to_leaves
# ---------------------------------------------------------------------------

class TestAssignGeooidsToLeaves:
    def test_direct_assignment_when_res15_is_leaf(self):
        cell = h3.latlng_to_cell(41.70, -71.55, 15)
        leaves = {cell}
        result = assign_geoids_to_leaves({"A": cell}, leaves)
        assert result["A"] == cell

    def test_ancestor_lookup(self):
        cell15 = h3.latlng_to_cell(41.70, -71.55, 15)
        leaf = h3.cell_to_parent(cell15, 8)
        leaves = {leaf}
        result = assign_geoids_to_leaves({"A": cell15}, leaves)
        assert result["A"] == leaf

    def test_all_geoids_assigned(self):
        parent14 = h3.cell_to_parent(h3.latlng_to_cell(41.70, -71.55, 15), 14)
        cells = list(h3.cell_to_children(parent14, 15))[:4]
        leaf = h3.cell_to_parent(cells[0], 8)
        leaves = {leaf}
        geoid_to_res15 = {str(i): c for i, c in enumerate(cells)}
        result = assign_geoids_to_leaves(geoid_to_res15, leaves)
        assert set(result.keys()) == set(geoid_to_res15.keys())
        assert all(v == leaf for v in result.values())


# ---------------------------------------------------------------------------
# TestIntegratedZeroPopLeaves
# ---------------------------------------------------------------------------

# ---------------------------------------------------------------------------
# TestComputeTopLevel
# ---------------------------------------------------------------------------

class TestComputeTopLevel:
    def test_single_leaf_returns_its_res0_ancestor(self):
        cell = h3.latlng_to_cell(41.70, -71.55, 5)
        result = compute_top_level({cell})
        assert len(result) == 1
        (tl,) = result
        assert h3.get_resolution(tl) <= 5

    def test_two_leaves_same_parent_returns_parent(self):
        parent = h3.latlng_to_cell(41.70, -71.55, 3)
        children = list(h3.cell_to_children(parent, 4))[:2]
        result = compute_top_level(set(children))
        assert result == {parent}

    def test_two_leaves_different_parents_returns_common_ancestor(self):
        parent = h3.latlng_to_cell(41.70, -71.55, 2)
        children3 = list(h3.cell_to_children(parent, 3))[:2]
        result = compute_top_level(set(children3))
        (tl,) = result
        assert h3.get_resolution(tl) <= 2
        for L in children3:
            assert h3.cell_to_parent(L, h3.get_resolution(tl)) == tl

    def test_empty_leaves_returns_empty(self):
        assert compute_top_level(set()) == set()

    def test_leaves_at_mixed_resolutions(self):
        parent3 = h3.latlng_to_cell(41.70, -71.55, 3)
        child4 = list(h3.cell_to_children(parent3, 4))[0]
        result = compute_top_level({parent3, child4})
        assert result == {parent3}


# ---------------------------------------------------------------------------
# TestFindZeroPopLeaves
# ---------------------------------------------------------------------------

class TestFindZeroPopLeaves:
    def _setup(self):
        """
        One res-13 parent with two populated clusters (res-14 A and C).
        Remaining res-14 siblings of parent13 are zero-pop.
        """
        parent13 = h3.latlng_to_cell(41.70, -71.55, 13)
        children14 = sorted(h3.cell_to_children(parent13, 14))
        parent14_A = children14[0]
        parent14_C = children14[-1]

        children15_A = list(h3.cell_to_children(parent14_A, 15))[:2]
        children15_C = list(h3.cell_to_children(parent14_C, 15))[:2]
        all_res15 = children15_A + children15_C
        geoid_to_res15 = {str(i): c for i, c in enumerate(all_res15)}
        block_pops = {str(i): 1000 for i in range(len(all_res15))}
        cell_pop, cdc, par = build_cell_hierarchy(geoid_to_res15, block_pops)
        pop_leaves, _ = compute_leaves_and_provisional_edges(cell_pop, cdc, par, threshold=1)
        top_level = compute_top_level(pop_leaves)
        return parent13, parent14_A, parent14_C, pop_leaves, top_level

    def test_populated_leaves_not_in_zero_pop_result(self):
        _, _, _, pop_leaves, top_level = self._setup()
        zero_pop = find_zero_pop_leaves(pop_leaves, top_level)
        assert pop_leaves.isdisjoint(zero_pop)

    def test_zero_pop_res14_siblings_become_leaves(self):
        parent13, parent14_A, parent14_C, pop_leaves, top_level = self._setup()
        zero_pop = find_zero_pop_leaves(pop_leaves, top_level)
        all_children14 = set(h3.cell_to_children(parent13, 14))
        expected_bridges = all_children14 - {parent14_A, parent14_C}
        assert expected_bridges.issubset(zero_pop)

    def test_all_leaves_cover_full_scope(self):
        # Every res-14 cell under parent13 must be either a populated leaf (at res-15)
        # or a zero-pop leaf (at res-14).
        parent13, _, _, pop_leaves, top_level = self._setup()
        zero_pop = find_zero_pop_leaves(pop_leaves, top_level)
        all_leaves = pop_leaves | zero_pop
        for cell14 in h3.cell_to_children(parent13, 14):
            is_pop_ancestor = any(
                h3.cell_to_parent(L, 14) == cell14
                for L in pop_leaves
                if h3.get_resolution(L) >= 14
            )
            if is_pop_ancestor:
                pass  # children of this cell are in pop_leaves
            else:
                assert cell14 in all_leaves, f"{cell14} not covered"

    def test_empty_top_level_returns_empty(self):
        _, _, _, pop_leaves, _ = self._setup()
        assert find_zero_pop_leaves(pop_leaves, set()) == set()

    def test_top_level_is_populated_leaf_returns_empty(self):
        # If the top_level cell is itself a populated leaf, no zero-pop leaves needed.
        cell = h3.latlng_to_cell(41.70, -71.55, 5)
        pop_leaves = {cell}
        top_level = {cell}
        result = find_zero_pop_leaves(pop_leaves, top_level)
        assert result == set()
