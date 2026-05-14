"""PyMETIS graph partitioning wrapper."""

import pymetis


def partition(
    adjacency_lists: list[list[int]],
    eweights: list[int],
    nweights: list[int],
    n_parts: int,
) -> list[int]:
    """
    Partition a graph into n_parts using PyMETIS.

    Parameters
    ----------
    adjacency_lists:
        adjacency_lists[i] is the list of neighbor indices for node i.
    eweights:
        Flat edge-weight list aligned with adjacency_lists (same ordering).
    nweights:
        Node weights (population). One value per node.
    n_parts:
        Number of partitions (districts).

    Returns
    -------
    membership : list[int]
        membership[i] is the partition index (0-based) assigned to node i.

    Raises
    ------
    ValueError
        If n_parts < 2 or n_parts >= len(nodes).
    """
    n_nodes = len(adjacency_lists)
    if n_parts < 2:
        raise ValueError(f"n_parts must be >= 2, got {n_parts}")
    if n_parts >= n_nodes:
        raise ValueError(
            f"n_parts ({n_parts}) must be less than the number of nodes ({n_nodes})"
        )

    _cut_count, membership = pymetis.part_graph(
        n_parts,
        adjacency=adjacency_lists,
        eweights=eweights,
        vweights=nweights,
        contiguous=True,
    )
    return list(membership)
