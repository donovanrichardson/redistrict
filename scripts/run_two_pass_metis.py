#!/usr/bin/env python3
"""
Two-pass METIS redistricting directly on census blocks.

Usage:
    uv run --env-file .env python scripts/run_two_pass_metis.py <statefp> <n_districts>

Example:
    uv run --env-file .env python scripts/run_two_pass_metis.py 44 2   # Rhode Island, 2 districts

Pipeline:
  1. Fetch all census blocks (including zero-pop) from DB
  2. Compute 99th-pct-by-pop-mass threshold
  3. Ensure rook contiguity adjacency is precomputed for all blocks
  4. Pass 1 METIS: nodes=blocks, weights=pop, n_parts=total_pop/threshold
     → produces ~threshold-pop clusters
  5. Derive cluster adjacency from block adjacency
  6. Pass 2 METIS: nodes=clusters, weights=cluster_pop, n_parts=n_districts
     → produces final district assignments
  7. Disaggregate: cluster→district, block→district
  8. Write run record, assignments, district geometries to DB
  9. Export GeoJSON to output/two_pass_runs/
"""
import math
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))

import numpy as np
from shapely.ops import unary_union
from tqdm import tqdm

from redistrict import db, partition
from redistrict.h3_graph import compute_threshold

_FIPS_TO_NAME = {
    "01": "Alabama",       "04": "Arizona",        "05": "Arkansas",
    "06": "California",    "08": "Colorado",        "09": "Connecticut",
    "10": "Delaware",      "11": "DC",              "12": "Florida",
    "13": "Georgia",       "16": "Idaho",           "17": "Illinois",
    "18": "Indiana",       "19": "Iowa",            "20": "Kansas",
    "21": "Kentucky",      "22": "Louisiana",       "23": "Maine",
    "24": "Maryland",      "25": "Massachusetts",   "26": "Michigan",
    "27": "Minnesota",     "28": "Mississippi",     "29": "Missouri",
    "30": "Montana",       "31": "Nebraska",        "32": "Nevada",
    "33": "New Hampshire", "34": "New Jersey",      "35": "New Mexico",
    "36": "New York",      "37": "North Carolina",  "38": "North Dakota",
    "39": "Ohio",          "40": "Oklahoma",        "41": "Oregon",
    "42": "Pennsylvania",  "44": "Rhode Island",    "45": "South Carolina",
    "46": "South Dakota",  "47": "Tennessee",       "48": "Texas",
    "49": "Utah",          "50": "Vermont",         "51": "Virginia",
    "53": "Washington",    "54": "West Virginia",   "55": "Wisconsin",
    "56": "Wyoming",
}

OUTPUT_DIR = os.path.join(os.path.dirname(__file__), "..", "output", "two_pass_runs")
NCUTS_PASS1 = 1   # coarse aggregation — quality matters less, keep memory low
NCUTS_PASS2 = 10  # final district pass — more trials for quality
NITER = 20
_HAVERSINE_SCALE = 10_000  # integer edge-weight ceiling for METIS


def _haversine_km(lat1: float, lon1: float, lat2: float, lon2: float) -> float:
    dlat = math.radians(lat2 - lat1)
    dlon = math.radians(lon2 - lon1)
    a = (math.sin(dlat / 2) ** 2
         + math.cos(math.radians(lat1)) * math.cos(math.radians(lat2))
         * math.sin(dlon / 2) ** 2)
    return 6371.0 * 2.0 * math.asin(math.sqrt(max(0.0, a)))


def _add_bridge_edges(
    blocks: list[dict],
    adjacency: set[tuple[str, str]],
) -> tuple[set[tuple[str, str]], int]:
    """
    Detect disconnected components in the block graph and add bridge edges so
    the graph becomes fully connected.

    Uses a minimum spanning tree over components: each bridge connects the two
    nearest blocks from any two currently-disconnected components, so island
    chains connect to their geographic neighbor rather than all jumping directly
    to the mainland.

    Returns (augmented_adjacency, n_bridges_added).
    """
    geoid_list = [b["geoid"] for b in blocks]
    idx_of = {g: i for i, g in enumerate(geoid_list)}
    n = len(geoid_list)

    # Build index-based adjacency for BFS
    adj: list[set[int]] = [set() for _ in range(n)]
    for a, b_geoid in adjacency:
        ia = idx_of.get(a)
        ib = idx_of.get(b_geoid)
        if ia is not None and ib is not None:
            adj[ia].add(ib)
            adj[ib].add(ia)

    # BFS to find connected components
    comp_of = [-1] * n
    components: list[list[int]] = []
    for start in range(n):
        if comp_of[start] != -1:
            continue
        comp: list[int] = []
        queue = [start]
        comp_of[start] = len(components)
        while queue:
            node = queue.pop()
            comp.append(node)
            for nb in adj[node]:
                if comp_of[nb] == -1:
                    comp_of[nb] = len(components)
                    queue.append(nb)
        components.append(comp)

    k = len(components)
    if k == 1:
        print(f"  Graph is fully connected ({n:,} blocks, 1 component).")
        return adjacency, 0

    block_pops = [int(b["pop"]) for b in blocks]
    comp_pops = [sum(block_pops[i] for i in comp) for comp in components]
    print(f"  Graph has {k} components: "
          + ", ".join(f"{len(components[j]):,} blocks (pop {comp_pops[j]:,})"
                      for j in range(min(k, 5)))
          + ("..." if k > 5 else ""))

    # Pre-build per-component lat/lon arrays, then reduce to convex hull vertices
    from scipy.spatial import ConvexHull

    def _hull_indices(comp: list[int]) -> list[int]:
        """Return the subset of comp indices that lie on the convex hull.
        Falls back to all indices when the component is too small for a hull."""
        if len(comp) < 4:
            return comp
        pts = np.column_stack([
            [float(blocks[i]["lat"]) for i in comp],
            [float(blocks[i]["lon"]) for i in comp],
        ])
        try:
            hull = ConvexHull(pts)
            return [comp[v] for v in hull.vertices]
        except Exception:
            return comp

    hull_comps = [_hull_indices(comp) for comp in components]

    comp_lats = [np.array([float(blocks[i]["lat"]) for i in h]) for h in hull_comps]
    comp_lons = [np.array([float(blocks[i]["lon"]) for i in h]) for h in hull_comps]
    comp_idx  = [np.array(h) for h in hull_comps]

    def _nearest_pair(ci: int, cj: int) -> tuple[int, int, float]:
        """Return (block_idx_in_ci, block_idx_in_cj, sq_dist) using hull vertices."""
        lats_i, lons_i = comp_lats[ci], comp_lons[ci]
        lats_j, lons_j = comp_lats[cj], comp_lons[cj]
        dlat = lats_i[:, None] - lats_j[None, :]
        dlon = lons_i[:, None] - lons_j[None, :]
        sq = dlat ** 2 + dlon ** 2
        flat = int(np.argmin(sq))
        ri, rj = divmod(flat, len(lats_j))
        return int(comp_idx[ci][ri]), int(comp_idx[cj][rj]), float(sq[ri, rj])

    # Prim's MST over the k components (k is tiny, O(k^2) is fine)
    in_tree = {0}
    augmented = set(adjacency)
    n_bridges = 0

    with tqdm(total=k - 1, desc="  Bridging components", unit="bridge") as pbar:
        while len(in_tree) < k:
            best = (float("inf"), -1, -1, -1, -1)  # (dist, ci, cj, bi, bj)
            for ci in in_tree:
                for cj in range(k):
                    if cj in in_tree:
                        continue
                    bi, bj, sq = _nearest_pair(ci, cj)
                    if sq < best[0]:
                        best = (sq, ci, cj, bi, bj)
            sq, ci, cj, bi, bj = best
            ga, gb = geoid_list[bi], geoid_list[bj]
            augmented.add((min(ga, gb), max(ga, gb)))
            in_tree.add(cj)
            n_bridges += 1
            pbar.update(1)

    print(f"  Added {n_bridges} bridge edge(s). Graph now fully connected.")
    return augmented, n_bridges


def _show_cluster_graph(
    cluster_nodes: list[dict],
    adjacency: set[tuple[str, str]],
    title: str,
) -> None:
    """Plot cluster centroids as nodes and adjacency edges as lines."""
    import matplotlib.pyplot as plt

    id_to_node = {n["geoid"]: n for n in cluster_nodes}
    lons = [float(n["lon"]) for n in cluster_nodes]
    lats = [float(n["lat"]) for n in cluster_nodes]

    fig, ax = plt.subplots(figsize=(12, 10))

    # Draw edges
    for ga, gb in adjacency:
        na, nb = id_to_node.get(ga), id_to_node.get(gb)
        if na and nb:
            ax.plot(
                [float(na["lon"]), float(nb["lon"])],
                [float(na["lat"]), float(nb["lat"])],
                color="steelblue", linewidth=0.3, alpha=0.4,
            )

    # Draw nodes sized by population
    pops = [max(1, n["pop"]) for n in cluster_nodes]
    max_pop = max(pops)
    sizes = [2 + 8 * (p / max_pop) for p in pops]
    ax.scatter(lons, lats, s=sizes, c="tomato", zorder=3, linewidths=0)

    ax.set_title(title, fontsize=12)
    ax.set_xlabel("Longitude")
    ax.set_ylabel("Latitude")
    ax.set_aspect("equal")
    plt.tight_layout()
    plt.show(block=False)
    plt.pause(0.1)



def _show_partition(
    geoid_to_part: dict[str, int],
    block_geoms_wkb: dict[str, bytes],
    title: str,
) -> dict[int, object]:
    """
    Union block geometries by partition group, then plot.
    Returns {part_id: unioned_shapely_geometry} for reuse in subsequent plots.
    """
    from shapely import wkb as swkb, get_parts
    from shapely.ops import unary_union
    import geopandas as gpd
    import matplotlib.pyplot as plt
    import matplotlib.patches as mpatches
    import matplotlib.colors as mcolors

    # Group geoms by partition
    groups: dict[int, list] = {}
    for geoid, part_id in geoid_to_part.items():
        raw = block_geoms_wkb.get(geoid)
        if raw:
            groups.setdefault(part_id, []).append(swkb.loads(raw))

    print(f"  Unioning {sum(len(v) for v in groups.values()):,} block polygons "
          f"into {len(groups):,} groups...")
    part_geoms: dict[int, object] = {p: unary_union(gs) for p, gs in groups.items()}

    n_parts = len(part_geoms)
    cmap = plt.cm.get_cmap("tab20")
    part_colour = {p: mcolors.to_hex(cmap((i % 20) / 20))
                   for i, p in enumerate(sorted(part_geoms))}

    rows = []
    for p, geom in part_geoms.items():
        colour = part_colour[p]
        polys = list(get_parts(geom))
        for poly in polys:
            rows.append({"geometry": poly, "colour": colour})
    gdf = gpd.GeoDataFrame(rows, crs="EPSG:4269")

    fig, ax = plt.subplots(figsize=(12, 10))
    gdf.plot(ax=ax, color=gdf["colour"], edgecolor="white", linewidth=0.3, alpha=0.85)
    ax.set_title(f"{title}\n{n_parts:,} groups", fontsize=12)
    ax.set_xlabel("Longitude")
    ax.set_ylabel("Latitude")
    ax.set_aspect("equal")
    plt.tight_layout()
    plt.show(block=False)
    plt.pause(0.1)

    return part_geoms


def _build_metis_graph(
    nodes: list,
    adjacency: set[tuple[str, str]],
    key: str = "geoid",
    pop_key: str = "pop",
    use_haversine: bool = False,
) -> tuple[list[list[int]], list[int], list[int]]:
    """
    Build PyMETIS-ready structures.
    nodes: list of dicts with `key` and `pop_key` (and `lat`/`lon` when use_haversine).
    adjacency: set of (id_a, id_b) string pairs.
    When use_haversine=True, edge weights are 1/distance (nearby = expensive to cut).
    Returns (adjacency_lists, eweights, nweights).
    """
    id_to_idx = {nd[key]: i for i, nd in enumerate(nodes)}
    nn = len(nodes)
    adj_sets: list[set[int]] = [set() for _ in range(nn)]
    raw_w: dict[tuple[int, int], float] = {}

    for a, b in adjacency:
        ia = id_to_idx.get(a)
        ib = id_to_idx.get(b)
        if ia is None or ib is None:
            continue
        adj_sets[ia].add(ib)
        adj_sets[ib].add(ia)
        if use_haversine:
            na, nb = nodes[ia], nodes[ib]
            d = max(0.01, _haversine_km(
                float(na["lat"]), float(na["lon"]),
                float(nb["lat"]), float(nb["lon"]),
            ))
            raw_w[(min(ia, ib), max(ia, ib))] = 1.0 / d

    adjacency_lists = [sorted(s) for s in adj_sets]

    if use_haversine and raw_w:
        max_raw = max(raw_w.values())
        scale = _HAVERSINE_SCALE / max_raw
        eweights = [
            max(1, round(raw_w[(min(i, j), max(i, j))] * scale))
            for i, neighbors in enumerate(adjacency_lists)
            for j in neighbors
        ]
    else:
        eweights = [1 for neighbors in adjacency_lists for _ in neighbors]

    nweights = [max(1, int(nd[pop_key])) for nd in nodes]
    return adjacency_lists, eweights, nweights


def _build_cluster_graph(
    blocks: list[dict],
    membership: list[int],
    block_adjacency: set[tuple[str, str]],
) -> tuple[list[dict], set[tuple[str, str]], list[int]]:
    """
    Aggregate blocks into cluster nodes and derive cluster adjacency.

    Compacts cluster IDs to remove METIS gaps (e.g. IDs 0,1,3 with no cluster 2),
    which would otherwise produce (0,0) centroid phantom nodes.

    Returns (cluster_nodes, cluster_adjacency, membership_compacted).
    cluster_nodes: list of {"geoid": str, "pop": int, "lat": float, "lon": float}
    membership_compacted: membership re-indexed to match cluster_nodes positions.
    """
    # Compact: remove gaps left by METIS in cluster ID space
    used = sorted(set(membership))
    old_to_new = {old: new for new, old in enumerate(used)}
    membership_c = [old_to_new[m] for m in membership]
    n_clusters = len(used)

    cluster_pop: list[int] = [0] * n_clusters
    cluster_lat_sum: list[float] = [0.0] * n_clusters
    cluster_lon_sum: list[float] = [0.0] * n_clusters
    cluster_count: list[int] = [0] * n_clusters

    for i, b in enumerate(blocks):
        c = membership_c[i]
        p = int(b["pop"])
        lat, lon = float(b["lat"]), float(b["lon"])
        cluster_pop[c] += p
        weight = p if p > 0 else 1
        cluster_lat_sum[c] += lat * weight
        cluster_lon_sum[c] += lon * weight
        cluster_count[c] += weight

    cluster_nodes = [
        {
            "geoid": str(c),
            "pop":   cluster_pop[c],
            "lat":   cluster_lat_sum[c] / cluster_count[c],
            "lon":   cluster_lon_sum[c] / cluster_count[c],
        }
        for c in range(n_clusters)
    ]

    geoid_to_idx = {b["geoid"]: i for i, b in enumerate(blocks)}
    cluster_adj: set[tuple[str, str]] = set()
    for a, b in block_adjacency:
        ia = geoid_to_idx.get(a)
        ib = geoid_to_idx.get(b)
        if ia is None or ib is None:
            continue
        ca, cb = membership_c[ia], membership_c[ib]
        if ca != cb:
            sa, sb = str(min(ca, cb)), str(max(ca, cb))
            cluster_adj.add((sa, sb))

    return cluster_nodes, cluster_adj, membership_c


def _repair_pass1_contiguity(
    geoid_to_cluster: dict[str, int],
    block_adjacency: set[tuple[str, str]],
) -> tuple[dict[str, int], int]:
    """
    Split non-contiguous Pass 1 clusters into singlepart clusters.
    Each disconnected fragment gets a new unique cluster ID.
    The largest (by block count) fragment keeps the original cluster ID.
    Returns (updated_geoid_to_cluster, n_fragments_split).
    """
    geoid_set = set(geoid_to_cluster)

    # Build bidirectional adjacency restricted to known geoids
    adj: dict[str, set[str]] = {g: set() for g in geoid_set}
    for a, b in block_adjacency:
        if a in geoid_set and b in geoid_set:
            adj[a].add(b)
            adj[b].add(a)

    # Group geoids by cluster
    cluster_to_geoids: dict[int, list[str]] = {}
    for geoid, c in geoid_to_cluster.items():
        cluster_to_geoids.setdefault(c, []).append(geoid)

    result = dict(geoid_to_cluster)
    next_id = max(geoid_to_cluster.values()) + 1
    n_splits = 0

    for c, geoids in cluster_to_geoids.items():
        if len(geoids) == 1:
            continue
        geoid_set_c = set(geoids)
        visited: set[str] = set()
        components: list[list[str]] = []
        for start in geoids:
            if start in visited:
                continue
            comp: list[str] = []
            queue = [start]
            visited.add(start)
            while queue:
                node = queue.pop()
                comp.append(node)
                for nb in adj[node]:
                    if nb in geoid_set_c and nb not in visited:
                        visited.add(nb)
                        queue.append(nb)
            components.append(comp)

        if len(components) <= 1:
            continue

        # Largest fragment keeps original ID; others get new IDs
        components.sort(key=len, reverse=True)
        for frag in components[1:]:
            for geoid in frag:
                result[geoid] = next_id
            next_id += 1
            n_splits += 1

    return result, n_splits


def main(statefp: str, n_districts: int) -> int:
    os.makedirs(OUTPUT_DIR, exist_ok=True)
    state_name = _FIPS_TO_NAME.get(statefp, statefp)
    conn = db.connect()
    try:
        db.ensure_tables(conn, "blocks")

        print(f"\nFetching census blocks for {state_name} (statefp={statefp})...")
        blocks = db.fetch_nodes(conn, "blocks", statefp)
        if not blocks:
            print("No blocks found. Is blocks_2020 loaded?")
            return
        n_blocks = len(blocks)
        total_pop = sum(int(b["pop"]) for b in blocks)
        active_blocks = [b for b in blocks if int(b["pop"]) > 0]
        zero_pop_blocks = [b for b in blocks if int(b["pop"]) == 0]
        print(f"  {n_blocks:,} blocks  ({len(active_blocks):,} populated,  "
              f"{len(zero_pop_blocks):,} zero-pop)")
        print(f"  Total population: {total_pop:,}")

        # --- Threshold ---
        threshold = compute_threshold(active_blocks)
        n_clusters_target = max(2, min(math.ceil(total_pop / threshold), n_blocks - 1))
        print(f"\nThreshold (99th pct by pop mass): {threshold:,}")
        print(f"  Pass-1 target clusters: {n_clusters_target:,}")

        # --- Block geometries (for visualisation) ---
        print("\nFetching block geometries for visualisation...")
        block_geoms_wkb = db.fetch_block_geoms(conn, statefp)
        print(f"  {len(block_geoms_wkb):,} block polygons loaded.")

        # --- Rook adjacency ---
        geoids = [b["geoid"] for b in blocks]
        print("\nChecking rook adjacency cache...")
        missing = db.get_missing_adjacency_geoids(conn, "blocks", geoids)
        if missing:
            print(f"  Computing adjacency for {len(missing):,} blocks (bulk spatial join)...")
            n_pairs = db.compute_and_store_adjacency_bulk(conn, "blocks", statefp)
            print(f"  {n_pairs:,} new pairs inserted.")
        else:
            print("  Adjacency fully cached.")
        block_adjacency = db.fetch_adjacency(conn, "blocks", geoids)
        print(f"  {len(block_adjacency):,} rook-adjacent block pairs.")

        # --- Bridge disconnected island components (populated only) ---
        print("\nChecking graph connectivity...")
        block_adjacency, n_bridges = _add_bridge_edges(active_blocks, block_adjacency)

        # --- Pass 1: populated blocks → clusters ---
        print(f"\nPass 1 METIS: {len(active_blocks):,} populated blocks "
              f"→ {n_clusters_target:,} clusters...")
        adj_lists_1, ew_1, nw_1 = _build_metis_graph(
            active_blocks, block_adjacency, use_haversine=True,
        )
        cut_1, membership_1 = partition.partition(
            adj_lists_1, ew_1, nw_1, n_clusters_target,
            ncuts=NCUTS_PASS1, niter=NITER, recursive=True,
        )
        print(f"  Edge cut: {cut_1:,}")
        n_clusters_actual = max(membership_1) + 1
        cluster_pops = [0] * n_clusters_actual
        for i, b in enumerate(active_blocks):
            cluster_pops[membership_1[i]] += int(b["pop"])
        ideal_cluster = total_pop / n_clusters_actual
        worst_cluster = max(abs(p - ideal_cluster) / ideal_cluster * 100
                            for p in cluster_pops if p > 0)
        print(f"  {n_clusters_actual:,} clusters, "
              f"ideal pop {ideal_cluster:,.0f}, worst deviation {worst_cluster:.1f}%")

        # --- Assign zero-pop blocks to nearest cluster centroid ---
        geoid_to_cluster: dict[str, int] = {
            active_blocks[i]["geoid"]: membership_1[i]
            for i in range(len(active_blocks))
        }
        if zero_pop_blocks:
            clust_lat: dict[int, list] = {}
            clust_lon: dict[int, list] = {}
            for b in active_blocks:
                c = geoid_to_cluster[b["geoid"]]
                clust_lat.setdefault(c, []).append(float(b["lat"]))
                clust_lon.setdefault(c, []).append(float(b["lon"]))
            clust_ids = sorted(clust_lat)
            cent_lats = np.array([sum(clust_lat[c]) / len(clust_lat[c]) for c in clust_ids])
            cent_lons = np.array([sum(clust_lon[c]) / len(clust_lon[c]) for c in clust_ids])
            for b in zero_pop_blocks:
                lat, lon = float(b["lat"]), float(b["lon"])
                d = (cent_lats - lat) ** 2 + (cent_lons - lon) ** 2
                geoid_to_cluster[b["geoid"]] = clust_ids[int(np.argmin(d))]
            print(f"  {len(zero_pop_blocks):,} zero-pop blocks assigned to nearest cluster.")

        # --- Contiguity repair: split non-contiguous clusters into singleparts ---
        print("Repairing Pass 1 contiguity...")
        geoid_to_cluster, n_splits = _repair_pass1_contiguity(
            geoid_to_cluster, block_adjacency,
        )
        n_clusters_actual = len(set(geoid_to_cluster.values()))
        print(f"  {n_splits:,} fragments split off → {n_clusters_actual:,} clusters total.")

        # Rebuild full membership aligned with `blocks` for _build_cluster_graph
        full_membership_1 = [geoid_to_cluster[b["geoid"]] for b in blocks]

        # --- Pass 2: clusters → districts ---
        print(f"\nPass 2 METIS: {n_clusters_actual:,} clusters → {n_districts} districts...")
        cluster_nodes, cluster_adj, full_membership_1c = _build_cluster_graph(
            blocks, full_membership_1, block_adjacency,
        )

        # --- Pass 1 visualisation (after compaction so cluster IDs are consistent) ---
        geoid_to_cluster_viz = {b["geoid"]: full_membership_1c[i] for i, b in enumerate(blocks)}
        cluster_geoms = _show_partition(
            geoid_to_cluster_viz, block_geoms_wkb,
            title=f"{state_name} — Pass 1: {n_clusters_actual:,} clusters "
                  f"(threshold={threshold:,})",
        )
        print("Checking cluster graph connectivity...")
        cluster_adj, n_cluster_bridges = _add_bridge_edges(cluster_nodes, cluster_adj)

        # --- Pass 2 input graph visualisation ---
        _show_cluster_graph(cluster_nodes, cluster_adj,
                            title=f"{state_name} — Pass 2 input: "
                                  f"{len(cluster_nodes):,} nodes, "
                                  f"{len(cluster_adj):,} edges")

        adj_lists_2, ew_2, nw_2 = _build_metis_graph(
            cluster_nodes, cluster_adj, use_haversine=True,
        )
        cut_2, membership_2 = partition.partition(
            adj_lists_2, ew_2, nw_2, n_districts,
            ncuts=NCUTS_PASS2, niter=NITER,
        )
        print(f"  Edge cut: {cut_2:,}")

        # --- Disaggregate: populated blocks → cluster → district ---
        block_to_district = {
            b["geoid"]: membership_2[full_membership_1c[i]]
            for i, b in enumerate(blocks) if int(b["pop"]) > 0
        }

        # --- Assign zero-pop blocks to nearest district centroid ---
        from collections import defaultdict
        dist_lat: dict[int, list] = defaultdict(list)
        dist_lon: dict[int, list] = defaultdict(list)
        for b in blocks:
            if int(b["pop"]) > 0:
                d = block_to_district[b["geoid"]]
                dist_lat[d].append(float(b["lat"]))
                dist_lon[d].append(float(b["lon"]))
        dist_ids = sorted(dist_lat)
        cent_lats = np.array([sum(dist_lat[d]) / len(dist_lat[d]) for d in dist_ids])
        cent_lons = np.array([sum(dist_lon[d]) / len(dist_lon[d]) for d in dist_ids])
        for b in blocks:
            if int(b["pop"]) == 0:
                lat, lon = float(b["lat"]), float(b["lon"])
                dists = (cent_lats - lat) ** 2 + (cent_lons - lon) ** 2
                block_to_district[b["geoid"]] = dist_ids[int(np.argmin(dists))]

        # --- Pass 2 visualisation (union cluster geoms by district) ---
        import geopandas as gpd
        import matplotlib.pyplot as plt
        import matplotlib.patches as mpatches
        import matplotlib.colors as mcolors
        from shapely.ops import unary_union

        cluster_to_district = {c: membership_2[c] for c in range(len(cluster_nodes))}
        n_d = len(set(cluster_to_district.values()))
        cmap2 = plt.cm.get_cmap("tab10")
        dist_colour = {d: mcolors.to_hex(cmap2(i / max(n_d - 1, 1)))
                       for i, d in enumerate(sorted(set(cluster_to_district.values())))}
        district_geoms: dict[int, list] = {}
        for c, geom in cluster_geoms.items():
            d = cluster_to_district[c]
            district_geoms.setdefault(d, []).append(geom)
        rows2 = [{"geometry": unary_union(gs), "colour": dist_colour[d]}
                 for d, gs in district_geoms.items()]
        gdf2 = gpd.GeoDataFrame(rows2, crs="EPSG:4269")
        fig2, ax2 = plt.subplots(figsize=(12, 10))
        gdf2.plot(ax=ax2, color=gdf2["colour"], edgecolor="white", linewidth=0.5, alpha=0.9)
        legend_elements = [
            mpatches.Patch(facecolor=dist_colour[d], label=f"District {d}")
            for d in sorted(dist_colour)
        ]
        ax2.legend(handles=legend_elements, loc="lower left", fontsize=9, title="District")
        ax2.set_title(
            f"{state_name} — Pass 2: {n_districts} districts  (edge cut {cut_2:,})",
            fontsize=12,
        )
        ax2.set_xlabel("Longitude")
        ax2.set_ylabel("Latitude")
        ax2.set_aspect("equal")
        plt.tight_layout()
        plt.show(block=False)
        plt.pause(0.1)

        # --- Population stats ---
        pop_per_district: dict[int, int] = {}
        for b in active_blocks:
            d = block_to_district[b["geoid"]]
            pop_per_district[d] = pop_per_district.get(d, 0) + int(b["pop"])
        ideal = total_pop / n_districts
        print("\nDistrict populations:")
        for d in sorted(pop_per_district):
            pct = 100 * (pop_per_district[d] - ideal) / ideal
            print(f"  District {d}: {pop_per_district[d]:,}  ({pct:+.1f}%)")
        worst = max(abs(100 * (p - ideal) / ideal) for p in pop_per_district.values())

        # --- Write to DB ---
        params = {
            "status":          "complete",
            "method":          "two_pass_metis",
            "threshold":       threshold,
            "n_blocks":        n_blocks,
            "n_h3_cells":      n_clusters_actual,
            "n_edges":         len(block_adjacency),
            "edge_cut":        cut_2,
            "ncuts":           NCUTS_PASS2,
            "niter":           NITER,
            "recursive":       False,
        }
        print("\nBuilding district geometries from cluster geoms...")
        district_cluster_geoms: dict[int, list] = {}
        for c, geom in cluster_geoms.items():
            d = cluster_to_district.get(c)
            if d is not None:
                district_cluster_geoms.setdefault(d, []).append(geom)
        district_geoms_wkt: dict[int, tuple[str, int]] = {
            dist_id: (unary_union(geoms).wkt, pop_per_district.get(dist_id, 0))
            for dist_id, geoms in district_cluster_geoms.items()
        }

        print("Writing results to database...")
        run_id = db.write_run(conn, "blocks", statefp, n_districts, params)
        db.write_assignments(conn, run_id, block_to_district)
        db.write_district_geoms_wkt(conn, run_id, district_geoms_wkt)
        print(f"  Run ID: {run_id}")

        # --- GeoJSON export ---
        slug = state_name.lower().replace(" ", "_")
        geojson_path = os.path.join(OUTPUT_DIR, f"{slug}_2pass_run{run_id}.geojson")
        db.export_geojson(conn, run_id, geojson_path)
        print(f"  GeoJSON -> {geojson_path}")

        # --- Deviation log ---
        log_path = os.path.join(
            OUTPUT_DIR, f"{slug}_2pass_run{run_id}_deviations.md"
        )
        lines = [
            f"# {state_name} — two-pass METIS — {n_districts} districts — Run {run_id}",
            f"",
            f"method=two_pass_metis  threshold={threshold}  "
            f"clusters={n_clusters_actual}  ncuts_p1={NCUTS_PASS1}  ncuts_p2={NCUTS_PASS2}  niter={NITER}",
            f"blocks={n_blocks:,}  edges={len(block_adjacency):,}  "
            f"districts={n_districts}  ideal={ideal:,.0f}  worst_deviation={worst:.1f}%",
            f"",
            f"| District | Population | Deviation |",
            f"|----------|-----------|-----------|",
        ]
        for d in sorted(pop_per_district):
            pop = pop_per_district[d]
            pct = 100 * (pop - ideal) / ideal
            lines.append(f"| {d} | {pop:,} | {pct:+.1f}% |")
        with open(log_path, "w") as fh:
            fh.write("\n".join(lines) + "\n")
        print(f"  Deviations -> {log_path}")

        return run_id

    finally:
        conn.close()


if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: run_two_pass_metis.py <statefp> <n_districts>")
        sys.exit(1)
    statefp = sys.argv[1]
    try:
        n_districts = int(sys.argv[2])
    except ValueError:
        print("n_districts must be an integer")
        sys.exit(1)
    main(statefp, n_districts)
