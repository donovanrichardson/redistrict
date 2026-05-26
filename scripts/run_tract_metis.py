#!/usr/bin/env python3
"""
Tract-based iterative METIS redistricting.

Usage:
    uv run --env-file .env python scripts/run_tract_metis.py <statefp> <n_districts>

Algorithm:
  1. Fetch all census blocks; group by tract (first 11 chars of GEOID20).
  2. Union block geometries per tract; explode MultiPolygons → singlepart subclusters.
  3. threshold = 2 × median singlepart population.
  4. While any subcluster has pop > threshold:
       For each oversized subcluster, find its rook-connected components.
       METIS-bisect (k=2, haversine-inverse edge weights) each component that
       still exceeds the threshold.  Split any non-contiguous METIS output into
       singlepart subclusters.
  5. Bridge disconnected subcluster components (convex hull + Kruskal's).
  6. METIS k-way (k = n_districts, uniform edge weights) on subclusters.
  7. Write run, assignments, and district geometries to DB; export GeoJSONs.
"""

import json
import math
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))

import numpy as np
from shapely.ops import unary_union

from redistrict import db, partition


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

OUTPUT_DIR = os.path.join(os.path.dirname(__file__), "..", "output", "tract_runs")
NCUTS_BISECT = 1   # coarse bisection — speed over quality
NITER_BISECT = 10
NCUTS_FINAL  = 10  # final district pass
NITER_FINAL  = 20
_HAVERSINE_SCALE = 10_000


def _tract_id(geoid20: str) -> str:
    """First 11 chars of GEOID20 = state(2)+county(3)+tract(6)."""
    return geoid20[:11]


def _haversine_km(lat1: float, lon1: float, lat2: float, lon2: float) -> float:
    dlat = math.radians(lat2 - lat1)
    dlon = math.radians(lon2 - lon1)
    a = (math.sin(dlat / 2) ** 2
         + math.cos(math.radians(lat1)) * math.cos(math.radians(lat2))
         * math.sin(dlon / 2) ** 2)
    return 6371.0 * 2.0 * math.asin(math.sqrt(max(0.0, a)))


def _build_rook_nbrs(
    adjacency: set[tuple[str, str]],
) -> dict[str, list[str]]:
    """Build a neighbour lookup from the full rook adjacency set (built once)."""
    rook_neighbours: dict[str, list[str]] = {}
    for a, b in adjacency:
        rook_neighbours.setdefault(a, []).append(b)
        rook_neighbours.setdefault(b, []).append(a)
    return rook_neighbours


def _rook_components(
    geoid_list: list[str],
    rook_neighbours: dict[str, list[str]],
) -> list[list[str]]:
    """BFS connected components of a block set under rook adjacency."""
    geoid_set = set(geoid_list)
    visited_geoids: set[str] = set()
    components: list[list[str]] = []
    for start_geoid in geoid_list:
        if start_geoid in visited_geoids:
            continue
        component: list[str] = []
        bfs_queue = [start_geoid]
        visited_geoids.add(start_geoid)
        while bfs_queue:
            node = bfs_queue.pop()
            component.append(node)
            for neighbour in rook_neighbours.get(node, []):
                if neighbour in geoid_set and neighbour not in visited_geoids:
                    visited_geoids.add(neighbour)
                    bfs_queue.append(neighbour)
        components.append(component)
    return components


def _bisect(
    geoid_list: list[str],
    rook_neighbours: dict[str, list[str]],
    blocks_by_geoid: dict[str, dict],
) -> list[list[str]]:
    """
    METIS k=2 bisection of a rook-connected block set.
    Returns a list of singlepart parts (>=2 if METIS output is non-contiguous).
    Falls back to the original list as a single part on error.
    """
    if len(geoid_list) < 2:
        return [geoid_list]

    geoid_set = set(geoid_list)
    index_of = {geoid: index for index, geoid in enumerate(geoid_list)}
    adjacency_sets: list[set[int]] = [set() for _ in geoid_list]
    raw_edge_weights: dict[tuple[int, int], float] = {}

    for geoid in geoid_list:
        index_geoid = index_of[geoid]
        for neighbour in rook_neighbours.get(geoid, []):
            if neighbour in geoid_set:
                index_neighbour = index_of[neighbour]
                adjacency_sets[index_geoid].add(index_neighbour)
                key = (min(index_geoid, index_neighbour), max(index_geoid, index_neighbour))
                if key not in raw_edge_weights:
                    distance = max(0.01, _haversine_km(
                        float(blocks_by_geoid[geoid]["lat"]), float(blocks_by_geoid[geoid]["lon"]),
                        float(blocks_by_geoid[neighbour]["lat"]), float(blocks_by_geoid[neighbour]["lon"]),
                    ))
                    raw_edge_weights[key] = 1.0 / distance

    adjacency_lists = [sorted(s) for s in adjacency_sets]
    max_weight = max(raw_edge_weights.values(), default=1.0)
    weight_scale = _HAVERSINE_SCALE / max_weight
    edge_weights: list[int] = []
    for index, neighbours in enumerate(adjacency_lists):
        for neighbour in neighbours:
            edge_weights.append(max(1, int(raw_edge_weights.get((min(index, neighbour), max(index, neighbour)), 1.0) * weight_scale)))
    node_weights = [max(1, int(blocks_by_geoid[geoid]["pop"])) for geoid in geoid_list]

    try:
        _, membership = partition.partition(
            adjacency_lists, edge_weights, node_weights, 2,
            ncuts=NCUTS_BISECT, niter=NITER_BISECT,
        )
    except Exception:
        return [geoid_list]

    partition_part_0 = [geoid_list[index] for index, m in enumerate(membership) if m == 0]
    partition_part_1 = [geoid_list[index] for index, m in enumerate(membership) if m == 1]

    bisect_result: list[list[str]] = []
    for partition_part in (partition_part_0, partition_part_1):
        if partition_part:
            bisect_result.extend(_rook_components(partition_part, rook_neighbours))
    return bisect_result if bisect_result else [geoid_list]


def _subcluster_pop(geoid_list: list[str], blocks_by_geoid: dict[str, dict]) -> int:
    return sum(int(blocks_by_geoid[geoid]["pop"]) for geoid in geoid_list)


def _subcluster_centroid(
    geoid_list: list[str],
    blocks_by_geoid: dict[str, dict],
) -> tuple[float, float]:
    """Population-weighted centroid using populated blocks only."""
    lat_sum = lon_sum = weight = 0.0
    for geoid in geoid_list:
        block = blocks_by_geoid[geoid]
        population = int(block["pop"])
        if population == 0:
            continue
        lat_sum += float(block["lat"]) * population
        lon_sum += float(block["lon"]) * population
        weight  += population
    if weight == 0:
        # fallback for all-zero-pop subclusters (should not arise in normal flow)
        latitudes = [float(blocks_by_geoid[geoid]["lat"]) for geoid in geoid_list]
        longitudes = [float(blocks_by_geoid[geoid]["lon"]) for geoid in geoid_list]
        return sum(latitudes) / len(latitudes), sum(longitudes) / len(longitudes)
    return lat_sum / weight, lon_sum / weight


def _build_subcluster_nodes(
    subclusters: dict[int, list[str]],
    blocks_by_geoid: dict[str, dict],
) -> list[dict]:
    nodes = []
    for subcluster_id in sorted(subclusters):
        geoid_list = subclusters[subcluster_id]
        pop = _subcluster_pop(geoid_list, blocks_by_geoid)
        lat, lon = _subcluster_centroid(geoid_list, blocks_by_geoid)
        nodes.append({"geoid": str(subcluster_id), "pop": pop, "lat": lat, "lon": lon})
    return nodes


def _build_subcluster_adj(
    subclusters: dict[int, list[str]],
    block_to_subcluster: dict[str, int],
    adjacency: set[tuple[str, str]],
    subcluster_nodes: list[dict],
) -> set[tuple[str, str]]:
    """
    Two subclusters are adjacent only if they share a rook block-pair edge AND
    their centroids are connected in the Delaunay triangulation of all subcluster
    centroids.  If a subcluster has rook neighbors but none pass the Delaunay
    test, it is connected to its single closest rook neighbor (by centroid
    distance) so it is not left isolated before bridging.
    """
    import numpy as np

    # --- Step 1: rook subcluster pairs ---
    rook_sub_adj: set[tuple[str, str]] = set()
    rook_neighbours: dict[str, set[str]] = {}
    for geoid_a, geoid_b in adjacency:
        sa = block_to_subcluster.get(geoid_a)
        sb = block_to_subcluster.get(geoid_b)
        if sa is not None and sb is not None and sa != sb:
            pair_a, pair_b = str(min(sa, sb)), str(max(sa, sb))
            rook_sub_adj.add((pair_a, pair_b))
            rook_neighbours.setdefault(str(sa), set()).add(str(sb))
            rook_neighbours.setdefault(str(sb), set()).add(str(sa))

    if len(subcluster_nodes) < 3:
        return rook_sub_adj

    # --- Step 2: sphere convex hull neighbor pairs ---
    # Project centroids onto the unit sphere, run 3D ConvexHull.
    # Each simplex is a triangle; its 3 edges define spatial neighbors.
    # This is equivalent to spherical Delaunay triangulation.
    from scipy.spatial import ConvexHull

    node_geoid_list = [n["geoid"] for n in subcluster_nodes]
    latitudes = np.radians([float(n["lat"]) for n in subcluster_nodes])
    longitudes = np.radians([float(n["lon"]) for n in subcluster_nodes])
    sphere_points = np.column_stack([
        np.cos(latitudes) * np.cos(longitudes),
        np.cos(latitudes) * np.sin(longitudes),
        np.sin(latitudes),
    ])
    try:
        hull = ConvexHull(sphere_points)
    except Exception:
        return rook_sub_adj

    delaunay_pairs: set[tuple[str, str]] = set()
    for simplex in hull.simplices:   # each simplex is a triangle (3 vertices)
        for i in range(3):
            for j in range(i + 1, 3):
                geoid_a = node_geoid_list[simplex[i]]
                geoid_b = node_geoid_list[simplex[j]]
                delaunay_pairs.add((min(geoid_a, geoid_b), max(geoid_a, geoid_b)))

    # --- Step 3: keep only rook edges that are also Delaunay neighbors ---
    subcluster_adjacency = rook_sub_adj & delaunay_pairs

    # --- Step 4: fallback for subclusters whose every rook neighbor was pruned ---
    centroids: dict[str, tuple[float, float]] = {
        n["geoid"]: (float(n["lat"]), float(n["lon"])) for n in subcluster_nodes
    }
    represented: set[str] = set()
    for geoid_a, geoid_b in subcluster_adjacency:
        represented.add(geoid_a)
        represented.add(geoid_b)

    for subcluster_id_str, neighbours in rook_neighbours.items():
        if subcluster_id_str in represented:
            continue
        # All rook neighbors were pruned by Delaunay filter — attach to closest one
        lat0, lon0 = centroids[subcluster_id_str]
        best_subcluster_id = min(
            neighbours,
            key=lambda s: _haversine_km(lat0, lon0, centroids[s][0], centroids[s][1]),
        )
        subcluster_adjacency.add((min(subcluster_id_str, best_subcluster_id), max(subcluster_id_str, best_subcluster_id)))

    return subcluster_adjacency


def _add_bridge_edges(
    blocks: list[dict],
    adjacency: set[tuple[str, str]],
    state_boundary=None,
) -> tuple[set[tuple[str, str]], int, list[int]]:
    """
    Detect disconnected components and add bridge edges.
    External nodes per component are determined by the 2D lat/lon convex hull
    of centroid coordinates within the component.

    Returns (augmented_adjacency, n_bridges_added, comp_of).
    """
    from scipy.spatial import ConvexHull

    geoid_list = [b["geoid"] for b in blocks]
    index_of = {geoid: index for index, geoid in enumerate(geoid_list)}
    node_count = len(geoid_list)

    adjacency_lists: list[set[int]] = [set() for _ in range(node_count)]
    for geoid_a, geoid_b in adjacency:
        index_a = index_of.get(geoid_a)
        index_b = index_of.get(geoid_b)
        if index_a is not None and index_b is not None:
            adjacency_lists[index_a].add(index_b)
            adjacency_lists[index_b].add(index_a)

    component_index_of = [-1] * node_count
    components: list[list[int]] = []
    for start in range(node_count):
        if component_index_of[start] != -1:
            continue
        component: list[int] = []
        bfs_queue = [start]
        component_index_of[start] = len(components)
        while bfs_queue:
            node = bfs_queue.pop()
            component.append(node)
            for neighbour in adjacency_lists[node]:
                if component_index_of[neighbour] == -1:
                    component_index_of[neighbour] = len(components)
                    bfs_queue.append(neighbour)
        components.append(component)

    component_count = len(components)
    if component_count == 1:
        print(f"  Subcluster graph fully connected ({node_count:,} nodes, 1 component).")
        return adjacency, 0, component_index_of

    block_populations = [int(b["pop"]) for b in blocks]
    component_populations = [sum(block_populations[block_index] for block_index in component) for component in components]
    print(f"  Subcluster graph has {component_count} components: "
          + ", ".join(f"{len(components[j]):,} nodes (pop {component_populations[j]:,})"
                      for j in range(min(component_count, 5)))
          + ("..." if component_count > 5 else ""))

    exterior_points: list[tuple[float, float]] = []
    exterior_component_indices: list[int] = []
    exterior_block_indices: list[int] = []

    for component_index, component in enumerate(components):
        if len(component) < 3:
            hull_indices = component
        else:
            component_pts = np.column_stack([
                [float(blocks[block_index]["lat"]) for block_index in component],
                [float(blocks[block_index]["lon"]) for block_index in component],
            ])
            try:
                hull_indices = [component[v] for v in ConvexHull(component_pts).vertices]
            except Exception:
                hull_indices = component

        for block_index in hull_indices:
            exterior_points.append((float(blocks[block_index]["lat"]), float(blocks[block_index]["lon"])))
            exterior_component_indices.append(component_index)
            exterior_block_indices.append(block_index)

    # Spherical convex hull over all exterior centroids projected onto unit sphere
    bridge_candidates: list[tuple[float, int, int, int, int]] = []
    latitudes_rad  = np.radians([point[0] for point in exterior_points])
    longitudes_rad = np.radians([point[1] for point in exterior_points])
    sphere_points = np.column_stack([
        np.cos(latitudes_rad) * np.cos(longitudes_rad),
        np.cos(latitudes_rad) * np.sin(longitudes_rad),
        np.sin(latitudes_rad),
    ])
    try:
        combined = ConvexHull(sphere_points)
        for simplex in combined.simplices:
            for i in range(3):
                for j in range(i + 1, 3):
                    ai, aj = int(simplex[i]), int(simplex[j])
                    ca, cj = exterior_component_indices[ai], exterior_component_indices[aj]
                    if ca != cj:
                        ba, bj = exterior_block_indices[ai], exterior_block_indices[aj]
                        distance = _haversine_km(
                            float(blocks[ba]["lat"]), float(blocks[ba]["lon"]),
                            float(blocks[bj]["lat"]), float(blocks[bj]["lon"]),
                        )
                        bridge_candidates.append((distance, ba, bj, ca, cj))
    except Exception:
        pass

    # Discard bridge candidates whose line segment exits the state boundary
    if state_boundary is not None and bridge_candidates:
        from shapely.geometry import LineString
        n_before = len(bridge_candidates)
        bridge_candidates = [
            (distance, ba, bj, ca, cj)
            for distance, ba, bj, ca, cj in bridge_candidates
            if state_boundary.covers(LineString([
                (float(blocks[ba]["lon"]), float(blocks[ba]["lat"])),
                (float(blocks[bj]["lon"]), float(blocks[bj]["lat"])),
            ]))
        ]
        n_discarded = n_before - len(bridge_candidates)
        if n_discarded:
            print(f"  Discarded {n_discarded} bridge candidate(s) crossing state boundary.")

    # Group candidates by component pair; for each pair add the shortest 1/3 (min 1)
    from tqdm import tqdm

    pair_to_candidates: dict[tuple[int, int], list[tuple[float, int, int]]] = {}
    for distance, block_index_a, block_index_b, ca, cj in bridge_candidates:
        pair_key = (min(ca, cj), max(ca, cj))
        pair_to_candidates.setdefault(pair_key, []).append((distance, block_index_a, block_index_b))

    print(f"  Bridging {len(pair_to_candidates)} component pair(s) "
          f"({sum(len(v) for v in pair_to_candidates.values())} total candidates)...")

    augmented_adjacency = set(adjacency)
    bridge_count = 0

    for (ca, cj), pair_candidates in tqdm(
        pair_to_candidates.items(),
        desc="  bridging pairs",
        unit="pair",
        leave=False,
    ):
        pair_candidates.sort()
        n_to_add = max(1, len(pair_candidates) // 3)
        for distance, block_index_a, block_index_b in pair_candidates[:n_to_add]:
            geoid_a, geoid_b = geoid_list[block_index_a], geoid_list[block_index_b]
            augmented_adjacency.add((min(geoid_a, geoid_b), max(geoid_a, geoid_b)))
            bridge_count += 1
        tqdm.write(
            f"    components ({ca}, {cj}): {len(pair_candidates)} candidate(s) → {n_to_add} bridge(s) added "
            f"(shortest {pair_candidates[0][0]:.1f} km, longest selected {pair_candidates[n_to_add - 1][0]:.1f} km)"
        )

    # Fallback: check remaining connectivity and add minimum bridges for any still-disconnected components
    union_find_parent = list(range(component_count))

    def _find(x: int) -> int:
        while union_find_parent[x] != x:
            union_find_parent[x] = union_find_parent[union_find_parent[x]]
            x = union_find_parent[x]
        return x

    for geoid_a, geoid_b in augmented_adjacency:
        index_a = index_of.get(geoid_a)
        index_b = index_of.get(geoid_b)
        if index_a is not None and index_b is not None:
            ca, cj = component_index_of[index_a], component_index_of[index_b]
            if _find(ca) != _find(cj):
                union_find_parent[_find(ca)] = _find(cj)

    remaining_components = len({_find(c) for c in range(component_count)})
    if remaining_components > 1:
        print(f"  {remaining_components} component(s) still disconnected — running fallback...")
        fallback_candidates: list[tuple[float, int, int, int, int]] = []
        exterior_point_count = len(exterior_block_indices)
        for ai in range(exterior_point_count):
            for aj in range(ai + 1, exterior_point_count):
                ca, cj = exterior_component_indices[ai], exterior_component_indices[aj]
                if _find(ca) == _find(cj):
                    continue
                lat_a, lon_a = exterior_points[ai]
                lat_j, lon_j = exterior_points[aj]
                fallback_candidates.append((
                    _haversine_km(lat_a, lon_a, lat_j, lon_j),
                    exterior_block_indices[ai], exterior_block_indices[aj], ca, cj,
                ))
        fallback_candidates.sort()
        fallback_bridge_count = 0
        for distance, block_index_a, block_index_b, ca, cj in fallback_candidates:
            if _find(ca) != _find(cj):
                geoid_a, geoid_b = geoid_list[block_index_a], geoid_list[block_index_b]
                augmented_adjacency.add((min(geoid_a, geoid_b), max(geoid_a, geoid_b)))
                union_find_parent[_find(ca)] = _find(cj)
                bridge_count += 1
                fallback_bridge_count += 1
                print(f"    fallback bridge ({ca} ↔ {cj}): {distance:.1f} km")
            if len({_find(c) for c in range(component_count)}) == 1:
                break
        print(f"  Fallback added {fallback_bridge_count} bridge(s).")

    print(f"  Added {bridge_count} bridge edge(s) total. Subcluster graph now fully connected.")
    return augmented_adjacency, bridge_count, component_index_of


def _show_subcluster_adjacency(
    subcluster_nodes: list[dict],
    subcluster_rook_adjacency: set[tuple[str, str]],
    subcluster_bridge_adjacency: set[tuple[str, str]],
    subcluster_geometries: dict[int, object],
    title: str,
) -> None:
    """Plot subcluster polygons overlaid with rook edges (blue) and bridge edges (red)."""
    import geopandas as gpd
    import matplotlib.pyplot as plt
    import matplotlib.colors as mcolors
    import matplotlib.lines as mlines

    # Node centroid lookup: geoid str -> (lon, lat) for matplotlib (x, y)
    centroids: dict[str, tuple[float, float]] = {
        n["geoid"]: (float(n["lon"]), float(n["lat"])) for n in subcluster_nodes
    }

    n_parts = len(subcluster_geometries)
    cmap = plt.cm.get_cmap("tab20")
    partition_colour = {p: mcolors.to_hex(cmap((i % 20) / 20))
                        for i, p in enumerate(sorted(subcluster_geometries))}

    subcluster_rows = [{"geometry": geometry, "colour": partition_colour[pid]}
                       for pid, geometry in subcluster_geometries.items()]
    gdf = gpd.GeoDataFrame(subcluster_rows, crs="EPSG:4269")

    fig, ax = plt.subplots(figsize=(14, 11))
    gdf.plot(ax=ax, color=gdf["colour"], edgecolor="white", linewidth=0.3, alpha=0.75)

    for geoid_a, geoid_b in subcluster_rook_adjacency:
        if geoid_a in centroids and geoid_b in centroids:
            x0, y0 = centroids[geoid_a]
            x1, y1 = centroids[geoid_b]
            ax.plot([x0, x1], [y0, y1], color="#3B82F6", linewidth=0.6, alpha=0.5)

    for geoid_a, geoid_b in subcluster_bridge_adjacency:
        if geoid_a in centroids and geoid_b in centroids:
            x0, y0 = centroids[geoid_a]
            x1, y1 = centroids[geoid_b]
            ax.plot([x0, x1], [y0, y1], color="#EF4444", linewidth=2.0, alpha=0.9, zorder=5)

    legend_elems = [
        mlines.Line2D([], [], color="#3B82F6", linewidth=1.5, label=f"Rook ({len(subcluster_rook_adjacency):,})"),
        mlines.Line2D([], [], color="#EF4444", linewidth=2.0, label=f"Bridge ({len(subcluster_bridge_adjacency):,})"),
    ]
    ax.legend(handles=legend_elems, loc="lower left", fontsize=9)
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
    from shapely import wkb as swkb
    from shapely.ops import unary_union
    import geopandas as gpd
    import matplotlib.pyplot as plt
    import matplotlib.patches as mpatches
    import matplotlib.colors as mcolors

    groups: dict[int, list] = {}
    for geoid, part_id in geoid_to_part.items():
        raw = block_geoms_wkb.get(geoid)
        if raw:
            groups.setdefault(part_id, []).append(swkb.loads(raw))

    print(f"  Unioning {sum(len(v) for v in groups.values()):,} block polygons "
          f"into {len(groups):,} groups...")
    partition_geoms: dict[int, object] = {p: unary_union(geoid_list) for p, geoid_list in groups.items()}

    n_parts = len(partition_geoms)
    cmap = plt.cm.get_cmap("tab20")
    partition_colour = {p: mcolors.to_hex(cmap((i % 20) / 20))
                        for i, p in enumerate(sorted(partition_geoms))}

    subcluster_rows = []
    for p, geometry in partition_geoms.items():
        subcluster_rows.append({"geometry": geometry, "colour": partition_colour[p]})
    gdf = gpd.GeoDataFrame(subcluster_rows, crs="EPSG:4269")
    fig, ax = plt.subplots(figsize=(12, 10))
    gdf.plot(ax=ax, color=gdf["colour"], edgecolor="white", linewidth=0.3, alpha=0.85)
    ax.set_title(title, fontsize=12)
    ax.set_xlabel("Longitude")
    ax.set_ylabel("Latitude")
    ax.set_aspect("equal")
    plt.tight_layout()
    plt.show(block=False)
    plt.pause(0.1)
    return partition_geoms


def main(statefp: str, n_districts: int) -> None:
    os.makedirs(OUTPUT_DIR, exist_ok=True)
    state_name = _FIPS_TO_NAME.get(statefp, statefp)
    conn = db.connect()
    try:
        db.ensure_tables(conn, "blocks")

        print(f"\nFetching census blocks for {state_name} (statefp={statefp})...")
        blocks = db.fetch_nodes(conn, "blocks", statefp)
        if not blocks:
            print("No blocks found.")
            return
        total_pop = sum(int(b["pop"]) for b in blocks)
        blocks_by_geoid = {b["geoid"]: b for b in blocks}
        populated_blocks = [b for b in blocks if int(b["pop"]) > 0]
        zero_pop_blocks  = [b for b in blocks if int(b["pop"]) == 0]
        populated_blocks_by_geoid = {b["geoid"]: b for b in populated_blocks}
        print(f"  {len(blocks):,} blocks, total pop {total_pop:,}")

        # --- State boundary ---
        print("\nFetching state boundary...")
        _state_boundary_wkb = db.fetch_state_boundary(conn, statefp)
        if _state_boundary_wkb:
            from shapely import wkb as swkb
            state_boundary = swkb.loads(_state_boundary_wkb)
            print("  State boundary loaded.")
        else:
            state_boundary = None
            print("  No county geometries found — state boundary filter disabled.")

        # --- Block geometries ---
        print("\nFetching block geometries...")
        block_geoms_wkb = db.fetch_block_geoms(conn, statefp)
        print(f"  {len(block_geoms_wkb):,} block polygons loaded.")

        # --- Rook adjacency ---
        all_block_geoids = [b["geoid"] for b in blocks]
        print("\nChecking rook adjacency cache...")
        missing_geoids = db.get_missing_adjacency_geoids(conn, "blocks", all_block_geoids)
        if missing_geoids:
            print(f"  Computing adjacency for {len(missing_geoids):,} blocks...")
            n_pairs = db.compute_and_store_adjacency_bulk(conn, "blocks", statefp)
            print(f"  {n_pairs:,} new pairs inserted.")
        else:
            print("  Adjacency fully cached.")
        rook_adjacency = db.fetch_adjacency(conn, "blocks", all_block_geoids)
        print(f"  {len(rook_adjacency):,} rook-adjacent block pairs.")
        rook_neighbours = _build_rook_nbrs(rook_adjacency)

        # --- Initial subclusters: one singlepart per tract (populated blocks only) ---
        from shapely import wkb as swkb
        from shapely.strtree import STRtree

        print(f"\nGrouping populated blocks by tract "
              f"({len(populated_blocks):,} populated, {len(zero_pop_blocks):,} zero-pop)...")
        tract_to_blocks: dict[str, list[str]] = {}
        for block in populated_blocks:
            tract_to_blocks.setdefault(_tract_id(block["geoid"]), []).append(block["geoid"])

        subclusters: dict[int, list[str]] = {}
        block_to_subcluster: dict[str, int] = {}
        next_subcluster_id = 0
        geometry_split_count = 0

        for tract_geoids in tract_to_blocks.values():
            block_geom_pairs = [
                (geoid, swkb.loads(block_geoms_wkb[geoid]))
                for geoid in tract_geoids if geoid in block_geoms_wkb
            ]
            if not block_geom_pairs:
                continue
            tract_geometry_union = unary_union([geometry for _, geometry in block_geom_pairs])
            if tract_geometry_union.geom_type == "Polygon":
                parts = [tract_geometry_union]
            else:
                parts = [geometry for geometry in tract_geometry_union.geoms
                         if geometry.geom_type in ("Polygon", "MultiPolygon")]

            if len(parts) <= 1:
                subclusters[next_subcluster_id] = tract_geoids
                for geoid in tract_geoids:
                    block_to_subcluster[geoid] = next_subcluster_id
                next_subcluster_id += 1
            else:
                geometry_split_count += len(parts) - 1
                spatial_index = STRtree(parts)
                part_block_geoids: dict[int, list[str]] = {i: [] for i in range(len(parts))}
                for geoid, geometry in block_geom_pairs:
                    spatial_hits = spatial_index.query(geometry, predicate="intersects")
                    part_index = int(spatial_hits[0]) if len(spatial_hits) > 0 else 0
                    part_block_geoids[part_index].append(geoid)
                for part_geoids in part_block_geoids.values():
                    if part_geoids:
                        subclusters[next_subcluster_id] = part_geoids
                        for geoid in part_geoids:
                            block_to_subcluster[geoid] = next_subcluster_id
                        next_subcluster_id += 1

        if geometry_split_count:
            print(f"  {geometry_split_count} tract(s) split into singleparts by geometry explosion.")

        # Threshold: 2x median singlepart population
        all_subcluster_populations = sorted(_subcluster_pop(geoid_list, populated_blocks_by_geoid) for geoid_list in subclusters.values())
        median_pop = all_subcluster_populations[len(all_subcluster_populations) // 2]
        threshold = 2.0 * median_pop
        print(f"  {len(tract_to_blocks):,} tracts → {len(subclusters):,} initial subclusters  "
              f"median pop = {median_pop:,.0f}  threshold = {threshold:,.0f}")

        # --- Iterative bisection ---
        from tqdm import tqdm

        print("\nIterative bisection until no subcluster exceeds threshold...")
        unsplittable_subcluster_ids: set[int] = set()
        iteration = 0
        while True:
            oversized = [
                subcluster_id for subcluster_id, geoid_list in subclusters.items()
                if subcluster_id not in unsplittable_subcluster_ids
                and _subcluster_pop(geoid_list, populated_blocks_by_geoid) > threshold
            ]
            if not oversized:
                break
            iteration += 1
            split_count = 0
            for subcluster_id in tqdm(oversized, desc=f"  Iter {iteration}", unit="sub"):
                geoid_list = subclusters.pop(subcluster_id)
                components = _rook_components(geoid_list, rook_neighbours)
                produced_subcluster_ids: list[int] = []
                for component in components:
                    component_population = _subcluster_pop(component, populated_blocks_by_geoid)
                    if component_population > threshold and len(component) > 1:
                        parts = _bisect(component, rook_neighbours, populated_blocks_by_geoid)
                        split_count += len(parts) - 1
                    else:
                        parts = [component]
                    for partition_part in parts:
                        subclusters[next_subcluster_id] = partition_part
                        for geoid in partition_part:
                            block_to_subcluster[geoid] = next_subcluster_id
                        produced_subcluster_ids.append(next_subcluster_id)
                        next_subcluster_id += 1
                if all(_subcluster_pop(subclusters[new_id], populated_blocks_by_geoid) > threshold
                       for new_id in produced_subcluster_ids):
                    unsplittable_subcluster_ids.update(produced_subcluster_ids)
            give_up_suffix = f" ({len(unsplittable_subcluster_ids)} given up)" if unsplittable_subcluster_ids else ""
            print(f"  Iter {iteration}: {len(oversized)} oversized → "
                  f"{split_count} additional splits → {len(subclusters):,} subclusters total{give_up_suffix}")

        print(f"  Done: {len(subclusters):,} subclusters after {iteration} iteration(s).")

        # --- Visualise subclusters ---
        geoid_to_sub_viz = {geoid: block_to_subcluster[geoid] for geoid in populated_blocks_by_geoid}
        subcluster_geoms = _show_partition(
            geoid_to_sub_viz, block_geoms_wkb,
            title=f"{state_name} — tract subclusters: {len(subclusters):,}",
        )

        # --- Build subcluster nodes and adjacency ---
        subcluster_nodes = _build_subcluster_nodes(subclusters, populated_blocks_by_geoid)
        sub_adj = _build_subcluster_adj(subclusters, block_to_subcluster, rook_adjacency, subcluster_nodes)

        # --- Bridge disconnected subcluster components ---
        print("\nChecking subcluster graph connectivity...")
        subcluster_rook_adjacency = set(sub_adj)
        sub_adj, bridge_count, subcluster_component_of = _add_bridge_edges(
            subcluster_nodes, sub_adj, state_boundary,
        )
        subcluster_bridge_adjacency = sub_adj - subcluster_rook_adjacency

        _show_subcluster_adjacency(
            subcluster_nodes, subcluster_rook_adjacency, subcluster_bridge_adjacency, subcluster_geoms,
            title=f"{state_name} — subcluster adjacency after bridging "
                  f"({len(subcluster_rook_adjacency):,} rook + {len(subcluster_bridge_adjacency)} bridge)",
        )

        # --- Final METIS: equal edge weights, k = n_districts ---
        print(f"\nFinal METIS: {len(subcluster_nodes):,} subclusters → {n_districts} districts...")
        subcluster_id_list = [int(n["geoid"]) for n in subcluster_nodes]
        subcluster_index = {n["geoid"]: i for i, n in enumerate(subcluster_nodes)}

        # Build adj lists with uniform edge weight
        adjacency_sets: list[set[int]] = [set() for _ in subcluster_nodes]
        for geoid_a, geoid_b in sub_adj:
            index_a = subcluster_index.get(geoid_a)
            index_b = subcluster_index.get(geoid_b)
            if index_a is not None and index_b is not None:
                adjacency_sets[index_a].add(index_b)
                adjacency_sets[index_b].add(index_a)
        adjacency_lists = [sorted(s) for s in adjacency_sets]
        edge_weights = [1] * sum(len(s) for s in adjacency_sets)   # uniform edge weights
        node_weights = [max(1, int(n["pop"])) for n in subcluster_nodes]

        _, membership = partition.partition(
            adjacency_lists, edge_weights, node_weights, n_districts,
            ncuts=NCUTS_FINAL, niter=NITER_FINAL,
        )
        print(f"  Done.")

        # --- Disaggregate populated blocks: subcluster → district ---
        subcluster_to_district = {subcluster_id_list[i]: membership[i]
                                   for i in range(len(subcluster_nodes))}

        # --- Build block → district for populated blocks ---
        block_to_district: dict[str, int] = {}
        for geoid in block_to_subcluster:
            subcluster_id = block_to_subcluster[geoid]
            if subcluster_id in subcluster_to_district:
                block_to_district[geoid] = subcluster_to_district[subcluster_id]

        # --- Assign zero-pop blocks by BFS wave-front from populated districts ---
        if zero_pop_blocks:
            print(f"\nAssigning {len(zero_pop_blocks):,} zero-pop blocks by wave-front...")

            # Population-weighted district centroids
            district_lat_sum: dict[int, float] = {}
            district_lon_sum: dict[int, float] = {}
            district_pop_sum: dict[int, float] = {}
            for block in populated_blocks:
                district = block_to_district.get(block["geoid"])
                if district is None:
                    continue
                population = float(block["pop"])
                district_lat_sum[district] = district_lat_sum.get(district, 0.0) + float(block["lat"]) * population
                district_lon_sum[district] = district_lon_sum.get(district, 0.0) + float(block["lon"]) * population
                district_pop_sum[district] = district_pop_sum.get(district, 0.0) + population
            district_centroids: dict[int, tuple[float, float]] = {
                district: (
                    district_lat_sum[district] / district_pop_sum[district],
                    district_lon_sum[district] / district_pop_sum[district],
                )
                for district in district_pop_sum
            }

            unassigned_zero_pop: set[str] = {b["geoid"] for b in zero_pop_blocks}
            wave_changed = True
            while wave_changed:
                wave_changed = False
                for geoid in list(unassigned_zero_pop):
                    adjacent_district_ids = {
                        block_to_district[neighbour]
                        for neighbour in rook_neighbours.get(geoid, [])
                        if neighbour in block_to_district
                    }
                    if not adjacent_district_ids:
                        continue
                    block = blocks_by_geoid[geoid]
                    latitude  = float(block["lat"])
                    longitude = float(block["lon"])
                    best_district = min(
                        adjacent_district_ids,
                        key=lambda district: (
                            (district_centroids[district][0] - latitude) ** 2
                            + (district_centroids[district][1] - longitude) ** 2
                            if district in district_centroids else float("inf")
                        ),
                    )
                    block_to_district[geoid] = best_district
                    unassigned_zero_pop.discard(geoid)
                    wave_changed = True

            # Isolated zero-pop blocks: BFS into components, assign via sphere convex hull
            if unassigned_zero_pop:
                from scipy.spatial import ConvexHull as _ConvexHull

                isolated_components: list[tuple[list[str], float, float]] = []
                visited_zero_pop: set[str] = set()
                for start_geoid in unassigned_zero_pop:
                    if start_geoid in visited_zero_pop:
                        continue
                    component: list[str] = []
                    queue = [start_geoid]
                    visited_zero_pop.add(start_geoid)
                    while queue:
                        current_geoid = queue.pop()
                        component.append(current_geoid)
                        for neighbour in rook_neighbours.get(current_geoid, []):
                            if neighbour in unassigned_zero_pop and neighbour not in visited_zero_pop:
                                visited_zero_pop.add(neighbour)
                                queue.append(neighbour)
                    latitudes  = [float(blocks_by_geoid[g]["lat"]) for g in component]
                    longitudes = [float(blocks_by_geoid[g]["lon"]) for g in component]
                    isolated_components.append((
                        component,
                        sum(latitudes) / len(latitudes),
                        sum(longitudes) / len(longitudes),
                    ))

                subcluster_nodes_by_geoid: dict[str, dict] = {n["geoid"]: n for n in subcluster_nodes}
                hull_nodes = list(subcluster_nodes)
                for index, (component, component_lat, component_lon) in enumerate(isolated_components):
                    hull_nodes.append({"geoid": f"__zp_{index}", "lat": component_lat, "lon": component_lon})

                latitudes_rad  = np.radians([float(n["lat"]) for n in hull_nodes])
                longitudes_rad = np.radians([float(n["lon"]) for n in hull_nodes])
                sphere_points = np.column_stack([
                    np.cos(latitudes_rad) * np.cos(longitudes_rad),
                    np.cos(latitudes_rad) * np.sin(longitudes_rad),
                    np.sin(latitudes_rad),
                ])
                hull_neighbours_by_geoid: dict[str, set[str]] = {}
                try:
                    hull = _ConvexHull(sphere_points)
                    for simplex in hull.simplices:
                        for i in range(3):
                            for j in range(i + 1, 3):
                                geoid_a = hull_nodes[simplex[i]]["geoid"]
                                geoid_b = hull_nodes[simplex[j]]["geoid"]
                                hull_neighbours_by_geoid.setdefault(geoid_a, set()).add(geoid_b)
                                hull_neighbours_by_geoid.setdefault(geoid_b, set()).add(geoid_a)
                except Exception:
                    pass

                new_district_assignments = 0
                for index, (component, component_lat, component_lon) in enumerate(isolated_components):
                    component_geoid_str = f"__zp_{index}"
                    hull_subcluster_neighbours = [
                        geoid_str for geoid_str in hull_neighbours_by_geoid.get(component_geoid_str, set())
                        if geoid_str in subcluster_nodes_by_geoid
                    ]
                    candidates = hull_subcluster_neighbours if hull_subcluster_neighbours \
                        else [n["geoid"] for n in subcluster_nodes]
                    best_subcluster_geoid = min(
                        candidates,
                        key=lambda geoid_str: _haversine_km(
                            component_lat, component_lon,
                            float(subcluster_nodes_by_geoid[geoid_str]["lat"]),
                            float(subcluster_nodes_by_geoid[geoid_str]["lon"]),
                        ),
                    )
                    best_district = subcluster_to_district[int(best_subcluster_geoid)]
                    for geoid in component:
                        block_to_district[geoid] = best_district
                    new_district_assignments += len(component)
                print(f"  {new_district_assignments:,} isolated zero-pop blocks assigned via sphere convex hull.")

            n_zero_pop_assigned = len(zero_pop_blocks) - len(unassigned_zero_pop) \
                if unassigned_zero_pop else len(zero_pop_blocks)
            print(f"  {n_zero_pop_assigned:,} zero-pop blocks assigned by rook adjacency wave-front.")

        # --- Population stats ---
        pop_per_district: dict[int, int] = {}
        for block in blocks:
            if int(block["pop"]) > 0:
                district = block_to_district[block["geoid"]]
                pop_per_district[district] = pop_per_district.get(district, 0) + int(block["pop"])
        ideal_population = total_pop / n_districts
        print("\nDistrict populations:")
        for district in sorted(pop_per_district):
            deviation_pct = 100 * (pop_per_district[district] - ideal_population) / ideal_population
            print(f"  District {district}: {pop_per_district[district]:,}  ({deviation_pct:+.1f}%)")
        worst_deviation = max(abs(100 * (population - ideal_population) / ideal_population) for population in pop_per_district.values())

        # --- Visualise districts ---
        import geopandas as gpd
        import matplotlib.pyplot as plt
        import matplotlib.patches as mpatches
        import matplotlib.colors as mcolors

        district_colour = {
            district: mcolors.to_hex(plt.cm.get_cmap("tab10")(i / max(n_districts - 1, 1)))
            for i, district in enumerate(sorted(pop_per_district))
        }
        district_geoms: dict[int, list] = {}
        for subcluster_id, geometry in subcluster_geoms.items():
            district = subcluster_to_district.get(subcluster_id)
            if district is not None:
                district_geoms.setdefault(district, []).append(geometry)
        if zero_pop_blocks:
            from shapely import wkb as swkb
            for block in zero_pop_blocks:
                geoid = block["geoid"]
                district = block_to_district.get(geoid)
                raw = block_geoms_wkb.get(geoid)
                if district is not None and raw:
                    district_geoms.setdefault(district, []).append(swkb.loads(raw))
        district_rows = [{"geometry": unary_union(geoid_list), "colour": district_colour[district]}
                         for district, geoid_list in district_geoms.items()]
        gdf2 = gpd.GeoDataFrame(district_rows, crs="EPSG:4269")
        district_fig, district_ax = plt.subplots(figsize=(12, 10))
        gdf2.plot(ax=district_ax, color=gdf2["colour"], edgecolor="white", linewidth=0.5, alpha=0.9)
        legend_elems = [
            mpatches.Patch(facecolor=district_colour[district], label=f"District {district}")
            for district in sorted(district_colour)
        ]
        district_ax.legend(handles=legend_elems, loc="lower left", fontsize=9)
        district_ax.set_title(f"{state_name} — {n_districts} districts", fontsize=12)
        district_ax.set_xlabel("Longitude")
        district_ax.set_ylabel("Latitude")
        district_ax.set_aspect("equal")
        plt.tight_layout()
        plt.show(block=False)
        plt.pause(0.1)

        # --- District geometries (union subcluster geoms by district) ---
        district_geoms_wkt: dict[int, tuple[str, int]] = {
            district: (unary_union(geoid_list).wkt, pop_per_district.get(district, 0))
            for district, geoid_list in district_geoms.items()
        }

        # --- Write to DB ---
        run_params = {
            "status":          "complete",
            "method":          "tract_metis",
            "n_tracts":        len(tract_to_blocks),
            "n_geo_splits":    geometry_split_count,
            "median_sub_pop":  int(median_pop),
            "threshold":       int(threshold),
            "n_subclusters":   len(subclusters),
            "bisect_iters":    iteration,
            "n_bridges":       bridge_count,
            "n_blocks":        len(blocks),
            "n_edges":         len(sub_adj),
            "districts":       n_districts,
            "ideal_pop":       int(ideal_population),
            "worst_deviation": round(worst_deviation, 2),
            "ncuts_final":     NCUTS_FINAL,
            "niter_final":     NITER_FINAL,
        }
        run_id = db.write_run(conn, "blocks", statefp, n_districts, run_params)
        print(f"\nRun ID: {run_id}")

        db.write_assignments(conn, run_id, block_to_district)
        print(f"  {len(block_to_district):,} block assignments written.")

        db.write_district_geoms_wkt(conn, run_id, district_geoms_wkt)
        print("  District geometries written.")

        # --- GeoJSON export ---
        state_slug = state_name.lower().replace(" ", "_")
        geojson_path = os.path.join(OUTPUT_DIR, f"{state_slug}_tract_run{run_id}.geojson")
        db.export_geojson(conn, run_id, geojson_path)
        print(f"  GeoJSON → {geojson_path}")

        # --- Subcluster components GeoJSON ---
        import shapely as _shapely
        subcluster_features = []
        for index, node in enumerate(subcluster_nodes):
            subcluster_id = int(node["geoid"])
            geometry = subcluster_geoms.get(subcluster_id)
            if geometry is None:
                continue
            subcluster_features.append({
                "type": "Feature",
                "geometry": json.loads(_shapely.to_geojson(geometry)),
                "properties": {
                    "subcluster_id": subcluster_id,
                    "component_id": subcluster_component_of[index],
                    "pop": node["pop"],
                    "district": subcluster_to_district.get(subcluster_id),
                },
            })
        subcluster_geojson_path = os.path.join(OUTPUT_DIR, f"{state_slug}_tract_run{run_id}_subclusters.geojson")
        with open(subcluster_geojson_path, "w") as file_handle:
            json.dump({"type": "FeatureCollection", "features": subcluster_features}, file_handle)
        print(f"  Subclusters GeoJSON → {subcluster_geojson_path}")

        # --- Deviation log ---
        deviation_log_path = os.path.join(OUTPUT_DIR, f"{state_slug}_tract_run{run_id}_deviations.md")
        lines = [
            f"# {state_name} — tract METIS — {n_districts} districts — Run {run_id}",
            "",
            f"method=tract_metis  threshold={int(threshold)}  "
            f"subclusters={len(subclusters)}  bisect_iters={iteration}  "
            f"ncuts={NCUTS_FINAL}  niter={NITER_FINAL}",
            f"blocks={len(blocks):,}  tracts={len(tract_to_blocks):,}  "
            f"districts={n_districts}  ideal={ideal_population:,.0f}  worst_deviation={worst_deviation:.1f}%",
            "",
            "## District populations",
            "",
            f"{'District':>10}  {'Population':>12}  {'Deviation':>10}",
        ]
        for district in sorted(pop_per_district):
            deviation_pct = 100 * (pop_per_district[district] - ideal_population) / ideal_population
            lines.append(f"{district:>10}  {pop_per_district[district]:>12,}  {deviation_pct:>+10.2f}%")
        with open(deviation_log_path, "w") as file_handle:
            file_handle.write("\n".join(lines) + "\n")
        print(f"  Deviation log → {deviation_log_path}")

        input("\nPress Enter to close plots...")

    finally:
        conn.close()


if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: run_tract_metis.py <statefp> <n_districts>")
        sys.exit(1)
    main(sys.argv[1], int(sys.argv[2]))
