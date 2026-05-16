"""PyMETIS graph partitioning wrapper."""

import pymetis


def partition(
    adjacency_lists: list[list[int]],
    eweights: list[int],
    nweights: list[int],
    n_parts: int,
    ncuts: int = 3,
    niter: int = 20,
) -> tuple[int, list[int]]:
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
    ncuts:
        Independent partitioning attempts; METIS keeps the best result.
        METIS default is 1. We default to 3 for better quality.
    niter:
        Refinement iterations per uncoarsening stage.
        METIS default is 10. We default to 20 for better quality.

    Returns
    -------
    (edge_cut, membership) where edge_cut is the total weighted cut and
    membership[i] is the partition index (0-based) assigned to node i.

    Raises
    ------
    ValueError
        If n_parts < 2 or n_parts >= number of nodes.
    """
    n_nodes = len(adjacency_lists)
    if n_parts < 2:
        raise ValueError(f"n_parts must be >= 2, got {n_parts}")
    if n_parts >= n_nodes:
        raise ValueError(
            f"n_parts ({n_parts}) must be less than the number of nodes ({n_nodes})"
        )

    # contig only works with k-way (recursive=False). Recursive bisection
    # ignores the contig flag entirely, producing disconnected districts.
    options = pymetis.Options(ncuts=ncuts, niter=niter, contig=1, ufactor=8, seed=42)

    result = pymetis.part_graph(
        n_parts,
        adjacency=adjacency_lists,
        eweights=eweights,
        vweights=nweights,
        options=options,
        recursive=False,
    )
    return result.edge_cuts, list(result.vertex_part)
