"""Tests for redistrict.cli (DB-free, using mocks)."""

from unittest.mock import MagicMock, call, patch

import pytest

from redistrict import cli


def _make_nodes(n: int = 6) -> list[dict]:
    """Minimal node list — enough for a Delaunay triangulation."""
    return [
        {"geoid": str(i), "pop": 1000 + i * 100,
         "lat": 41.70 + i * 0.05, "lon": -71.55 + i * 0.03}
        for i in range(n)
    ]


# ---------------------------------------------------------------------------
# run() — stops before METIS, saves edges, status=pending
# ---------------------------------------------------------------------------

class TestRunPendingBehaviour:
    """run() should save a pending stub and NOT call partition.partition."""

    def _patch_db(self, nodes):
        m = MagicMock()
        m.connect.return_value.__enter__ = lambda s: s
        m.connect.return_value.__exit__ = MagicMock(return_value=False)
        conn = MagicMock()
        m.connect.return_value = conn
        m.fetch_nodes.return_value = nodes
        m.fetch_available_states.return_value = ["44"]
        m.get_missing_adjacency_geoids.return_value = []
        m.fetch_adjacency.return_value = set()
        m.write_run.return_value = 99
        return m

    def test_run_does_not_call_partition(self):
        nodes = _make_nodes()
        mock_db = self._patch_db(nodes)

        with patch("redistrict.cli.db", mock_db), \
             patch("redistrict.cli.partition") as mock_partition:
            cli.run("44", "tracts", 2)

        mock_partition.partition.assert_not_called()

    def test_run_saves_pending_status(self):
        nodes = _make_nodes()
        mock_db = self._patch_db(nodes)

        with patch("redistrict.cli.db", mock_db):
            cli.run("44", "tracts", 2)

        call_args = mock_db.write_run.call_args
        params = call_args[0][4]  # positional: conn, geography, statefp, n_districts, params
        assert params["status"] == "pending"

    def test_run_does_not_write_assignments(self):
        nodes = _make_nodes()
        mock_db = self._patch_db(nodes)

        with patch("redistrict.cli.db", mock_db):
            cli.run("44", "tracts", 2)

        mock_db.write_assignments.assert_not_called()

    def test_run_does_not_write_district_geoms(self):
        nodes = _make_nodes()
        mock_db = self._patch_db(nodes)

        with patch("redistrict.cli.db", mock_db):
            cli.run("44", "tracts", 2)

        mock_db.write_district_geoms.assert_not_called()

    def test_run_calls_write_edges(self):
        nodes = _make_nodes()
        mock_db = self._patch_db(nodes)

        with patch("redistrict.cli.db", mock_db):
            cli.run("44", "tracts", 2)

        mock_db.write_edges.assert_called_once()

    def test_run_returns_run_id(self):
        nodes = _make_nodes()
        mock_db = self._patch_db(nodes)

        with patch("redistrict.cli.db", mock_db):
            result = cli.run("44", "tracts", 2)

        assert result == 99

    def test_run_excludes_zero_pop_nodes_from_graph(self):
        """Zero-pop nodes must not appear in the edge set saved to DB."""
        nodes = _make_nodes(6)
        nodes[0]["pop"] = 0  # make first node zero-pop
        mock_db = self._patch_db(nodes)

        with patch("redistrict.cli.db", mock_db):
            cli.run("44", "tracts", 2)

        # write_edges receives active_nodes — zero-pop geoid "0" must not appear
        call_args = mock_db.write_edges.call_args
        active_nodes_arg = call_args[0][2]  # positional: conn, run_id, nodes, edges, adj_pairs
        geoids = {n["geoid"] for n in active_nodes_arg}
        assert "0" not in geoids

    def test_run_records_zero_pop_count_in_params(self):
        nodes = _make_nodes(6)
        nodes[0]["pop"] = 0
        nodes[1]["pop"] = 0
        mock_db = self._patch_db(nodes)

        with patch("redistrict.cli.db", mock_db):
            cli.run("44", "tracts", 2)

        params = mock_db.write_run.call_args[0][4]
        assert params["n_zero_pop_nodes"] == 2
