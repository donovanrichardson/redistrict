"""Tests for pure functions in scripts/run_tract_metis.py (no DB, no METIS required)."""

import importlib.util
import os
import sys

import pytest

# ---------------------------------------------------------------------------
# Load the script as a module without executing main()
# ---------------------------------------------------------------------------

_SCRIPT_PATH = os.path.join(
    os.path.dirname(__file__), "..", "scripts", "run_tract_metis.py"
)

spec = importlib.util.spec_from_file_location("run_tract_metis", _SCRIPT_PATH)
tract_metis = importlib.util.module_from_spec(spec)
# Insert src/ so the script's own sys.path manipulation works
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))
spec.loader.exec_module(tract_metis)

_tract_id            = tract_metis._tract_id
_haversine_km        = tract_metis._haversine_km
_build_rook_nbrs     = tract_metis._build_rook_nbrs
_rook_components     = tract_metis._rook_components
_subcluster_pop      = tract_metis._subcluster_pop
_subcluster_centroid = tract_metis._subcluster_centroid
_build_subcluster_nodes = tract_metis._build_subcluster_nodes
_build_subcluster_adj   = tract_metis._build_subcluster_adj
_bisect              = tract_metis._bisect


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _make_block(geoid, pop, lat, lon):
    return {"geoid": geoid, "pop": pop, "lat": lat, "lon": lon}


# ---------------------------------------------------------------------------
# _tract_id
# ---------------------------------------------------------------------------

class TestTractId:
    def test_returns_first_11_chars(self):
        assert _tract_id("012345678901234") == "01234567890"

    def test_with_real_format(self):
        # State(2) + county(3) + tract(6) + block(4)
        geoid = "25017341100001"  # Massachusetts, Bristol county
        assert _tract_id(geoid) == "25017341100"

    def test_length_is_11(self):
        assert len(_tract_id("A" * 15)) == 11


# ---------------------------------------------------------------------------
# _haversine_km
# ---------------------------------------------------------------------------

class TestHaversineKm:
    def test_same_point_is_zero(self):
        assert _haversine_km(42.0, -71.0, 42.0, -71.0) == pytest.approx(0.0)

    def test_known_distance_boston_providence(self):
        distance = _haversine_km(42.360, -71.058, 41.824, -71.412)
        assert 63 < distance < 70

    def test_symmetry(self):
        a = _haversine_km(40.0, -75.0, 41.0, -74.0)
        b = _haversine_km(41.0, -74.0, 40.0, -75.0)
        assert a == pytest.approx(b)

    def test_longer_distance_is_larger(self):
        near = _haversine_km(42.0, -71.0, 42.01, -71.0)
        far  = _haversine_km(42.0, -71.0, 43.0,  -71.0)
        assert far > near


# ---------------------------------------------------------------------------
# _build_rook_nbrs
# ---------------------------------------------------------------------------

class TestBuildRookNbrs:
    def test_empty_adjacency(self):
        result = _build_rook_nbrs(set())
        assert result == {}

    def test_single_pair(self):
        result = _build_rook_nbrs({("A", "B")})
        assert "B" in result["A"]
        assert "A" in result["B"]

    def test_symmetric(self):
        pairs = {("A", "B"), ("B", "C"), ("A", "C")}
        result = _build_rook_nbrs(pairs)
        for geoid, neighbours in result.items():
            for neighbour in neighbours:
                assert geoid in result[neighbour]

    def test_all_geoids_present_as_keys(self):
        pairs = {("A", "B"), ("C", "D")}
        result = _build_rook_nbrs(pairs)
        assert set(result.keys()) == {"A", "B", "C", "D"}

    def test_no_self_loops(self):
        result = _build_rook_nbrs({("A", "B"), ("B", "C")})
        for geoid, neighbours in result.items():
            assert geoid not in neighbours


# ---------------------------------------------------------------------------
# _rook_components
# ---------------------------------------------------------------------------

class TestRookComponents:
    def _chain_nbrs(self):
        # A-B-C-D chain
        return _build_rook_nbrs({("A", "B"), ("B", "C"), ("C", "D")})

    def test_fully_connected_is_one_component(self):
        neighbours = self._chain_nbrs()
        components = _rook_components(["A", "B", "C", "D"], neighbours)
        assert len(components) == 1
        assert set(components[0]) == {"A", "B", "C", "D"}

    def test_two_disconnected_components(self):
        neighbours = _build_rook_nbrs({("A", "B"), ("C", "D")})
        components = _rook_components(["A", "B", "C", "D"], neighbours)
        assert len(components) == 2
        component_sets = [set(c) for c in components]
        assert {"A", "B"} in component_sets
        assert {"C", "D"} in component_sets

    def test_isolated_node_is_own_component(self):
        neighbours = _build_rook_nbrs({("A", "B")})
        components = _rook_components(["A", "B", "C"], neighbours)
        assert len(components) == 2

    def test_single_node(self):
        components = _rook_components(["A"], {})
        assert len(components) == 1
        assert components[0] == ["A"]

    def test_subset_of_adjacency_used(self):
        # Full adjacency has A-B-C-D, but we only query A and C
        neighbours = self._chain_nbrs()
        components = _rook_components(["A", "C"], neighbours)
        # A and C are not directly adjacent, so two components
        assert len(components) == 2

    def test_all_geoids_covered(self):
        neighbours = self._chain_nbrs()
        components = _rook_components(["A", "B", "C", "D"], neighbours)
        all_geoids = {geoid for component in components for geoid in component}
        assert all_geoids == {"A", "B", "C", "D"}


# ---------------------------------------------------------------------------
# _subcluster_pop
# ---------------------------------------------------------------------------

class TestSubclusterPop:
    def _blocks(self):
        return {
            "A": _make_block("A", 100, 42.0, -71.0),
            "B": _make_block("B", 200, 42.1, -71.1),
            "C": _make_block("C",   0, 42.2, -71.2),
        }

    def test_sum_of_populations(self):
        blocks = self._blocks()
        assert _subcluster_pop(["A", "B"], blocks) == 300

    def test_zero_pop_block(self):
        blocks = self._blocks()
        assert _subcluster_pop(["A", "C"], blocks) == 100

    def test_empty_list(self):
        assert _subcluster_pop([], self._blocks()) == 0

    def test_all_zero(self):
        blocks = self._blocks()
        assert _subcluster_pop(["C"], blocks) == 0


# ---------------------------------------------------------------------------
# _subcluster_centroid
# ---------------------------------------------------------------------------

class TestSubclusterCentroid:
    def _blocks(self):
        return {
            "A": _make_block("A", 100, 42.0, -71.0),
            "B": _make_block("B", 100, 43.0, -72.0),
            "Z": _make_block("Z",   0, 99.0,   0.0),  # zero-pop, should be ignored
        }

    def test_equal_weights_is_midpoint(self):
        blocks = self._blocks()
        latitude, longitude = _subcluster_centroid(["A", "B"], blocks)
        assert latitude  == pytest.approx(42.5)
        assert longitude == pytest.approx(-71.5)

    def test_zero_pop_blocks_ignored(self):
        blocks = self._blocks()
        latitude, longitude = _subcluster_centroid(["A", "Z"], blocks)
        assert latitude  == pytest.approx(42.0)
        assert longitude == pytest.approx(-71.0)

    def test_population_weighted(self):
        blocks = {
            "A": _make_block("A", 100, 40.0, -70.0),
            "B": _make_block("B", 300, 44.0, -74.0),
        }
        latitude, longitude = _subcluster_centroid(["A", "B"], blocks)
        # weighted: (100*40 + 300*44)/400 = 43.0, (100*-70 + 300*-74)/400 = -73.0
        assert latitude  == pytest.approx(43.0)
        assert longitude == pytest.approx(-73.0)

    def test_all_zero_pop_fallback(self):
        blocks = {
            "A": _make_block("A", 0, 42.0, -71.0),
            "B": _make_block("B", 0, 44.0, -73.0),
        }
        latitude, longitude = _subcluster_centroid(["A", "B"], blocks)
        assert latitude  == pytest.approx(43.0)
        assert longitude == pytest.approx(-72.0)


# ---------------------------------------------------------------------------
# _build_subcluster_nodes
# ---------------------------------------------------------------------------

class TestBuildSubclusterNodes:
    def _blocks(self):
        return {
            "A": _make_block("A", 100, 42.0, -71.0),
            "B": _make_block("B", 200, 43.0, -72.0),
            "C": _make_block("C", 300, 44.0, -73.0),
        }

    def test_one_node_per_subcluster(self):
        subclusters = {0: ["A", "B"], 1: ["C"]}
        nodes = _build_subcluster_nodes(subclusters, self._blocks())
        assert len(nodes) == 2

    def test_node_geoids_are_subcluster_ids(self):
        subclusters = {0: ["A"], 1: ["B"]}
        nodes = _build_subcluster_nodes(subclusters, self._blocks())
        geoids = {n["geoid"] for n in nodes}
        assert geoids == {"0", "1"}

    def test_population_summed(self):
        subclusters = {0: ["A", "B"]}
        nodes = _build_subcluster_nodes(subclusters, self._blocks())
        assert nodes[0]["pop"] == 300

    def test_ordered_by_subcluster_id(self):
        subclusters = {5: ["C"], 2: ["A"], 9: ["B"]}
        nodes = _build_subcluster_nodes(subclusters, self._blocks())
        geoids = [n["geoid"] for n in nodes]
        assert geoids == ["2", "5", "9"]


# ---------------------------------------------------------------------------
# _build_subcluster_adj
# ---------------------------------------------------------------------------

class TestBuildSubclusterAdj:
    def _blocks(self):
        # 5 blocks in a rough grid (lat/lon chosen so sphere hull works)
        return {
            "A": _make_block("A", 100, 42.00, -71.00),
            "B": _make_block("B", 100, 42.01, -71.00),
            "C": _make_block("C", 100, 42.02, -71.00),
            "D": _make_block("D", 100, 42.03, -71.00),
            "E": _make_block("E", 100, 42.04, -71.00),
        }

    def test_same_subcluster_pair_excluded(self):
        blocks = self._blocks()
        subclusters = {0: ["A", "B"], 1: ["C", "D", "E"]}
        block_to_subcluster = {"A": 0, "B": 0, "C": 1, "D": 1, "E": 1}
        # A-B are in the same subcluster; only B-C crosses the boundary
        adjacency = {("A", "B"), ("B", "C"), ("C", "D")}
        nodes = _build_subcluster_nodes(subclusters, blocks)
        result = _build_subcluster_adj(subclusters, block_to_subcluster, adjacency, nodes)
        # Only the B-C pair should produce an edge between subclusters 0 and 1
        for geoid_a, geoid_b in result:
            assert geoid_a != geoid_b

    def test_cross_boundary_pair_included(self):
        blocks = self._blocks()
        subclusters = {0: ["A", "B"], 1: ["C", "D", "E"]}
        block_to_subcluster = {"A": 0, "B": 0, "C": 1, "D": 1, "E": 1}
        adjacency = {("B", "C")}
        nodes = _build_subcluster_nodes(subclusters, blocks)
        result = _build_subcluster_adj(subclusters, block_to_subcluster, adjacency, nodes)
        assert ("0", "1") in result or ("1", "0") in result

    def test_no_adjacency_returns_empty(self):
        blocks = self._blocks()
        subclusters = {0: ["A", "B"], 1: ["C"]}
        block_to_subcluster = {"A": 0, "B": 0, "C": 1}
        nodes = _build_subcluster_nodes(subclusters, blocks)
        result = _build_subcluster_adj(subclusters, block_to_subcluster, set(), nodes)
        assert result == set()

    def test_pairs_are_canonically_ordered(self):
        blocks = self._blocks()
        subclusters = {0: ["A"], 1: ["B"]}
        block_to_subcluster = {"A": 0, "B": 1}
        adjacency = {("A", "B")}
        nodes = _build_subcluster_nodes(subclusters, blocks)
        result = _build_subcluster_adj(subclusters, block_to_subcluster, adjacency, nodes)
        for geoid_a, geoid_b in result:
            assert geoid_a <= geoid_b


# ---------------------------------------------------------------------------
# _bisect (integration — requires pymetis)
# ---------------------------------------------------------------------------

class TestBisect:
    def _linear_blocks(self, count=6):
        """count blocks in a line, each 0.01 degrees apart."""
        blocks = {}
        for i in range(count):
            geoid = str(i)
            blocks[geoid] = _make_block(geoid, 100, 42.0 + i * 0.01, -71.0)
        return blocks

    def _chain_nbrs(self, count=6):
        pairs = {(str(i), str(i + 1)) for i in range(count - 1)}
        return _build_rook_nbrs(pairs)

    def test_returns_two_parts(self):
        blocks = self._linear_blocks(6)
        neighbours = self._chain_nbrs(6)
        parts = _bisect(list(blocks.keys()), neighbours, blocks)
        assert len(parts) >= 2

    def test_all_geoids_covered(self):
        blocks = self._linear_blocks(6)
        neighbours = self._chain_nbrs(6)
        parts = _bisect(list(blocks.keys()), neighbours, blocks)
        covered = {geoid for part in parts for geoid in part}
        assert covered == set(blocks.keys())

    def test_single_block_returns_unchanged(self):
        blocks = {"A": _make_block("A", 100, 42.0, -71.0)}
        parts = _bisect(["A"], {}, blocks)
        assert parts == [["A"]]

    def test_parts_are_nonempty(self):
        blocks = self._linear_blocks(6)
        neighbours = self._chain_nbrs(6)
        parts = _bisect(list(blocks.keys()), neighbours, blocks)
        for part in parts:
            assert len(part) > 0
