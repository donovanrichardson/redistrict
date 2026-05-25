#!/usr/bin/env python3
"""
Tract-based iterative METIS redistricting.

Usage:
    uv run --env-file .env python scripts/run_tract_metis.py <statefp> <n_districts>

Algorithm:
  1. Fetch all census blocks; group by tract (first 11 chars of GEOID20).
  2. avg_tract_pop = total_pop / n_populated_tracts;  threshold = 2 × avg_tract_pop.
  3. Each tract starts as one subcluster.
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


def _rook_components(
    geoids: list[str],
    adjacency: set[tuple[str, str]],
) -> list[list[str]]:
    """BFS connected components of a block set under rook adjacency."""
    geoid_set = set(geoids)
    adj: dict[str, list[str]] = {g: [] for g in geoid_set}
    for a, b in adjacency:
        if a in geoid_set and b in geoid_set:
            adj[a].append(b)
            adj[b].append(a)
    visited: set[str] = set()
    comps: list[list[str]] = []
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
                if nb not in visited:
                    visited.add(nb)
                    queue.append(nb)
        comps.append(comp)
    return comps


def _bisect(
    geoids: list[str],
    adjacency: set[tuple[str, str]],
    blocks_by_geoid: dict[str, dict],
) -> list[list[str]]:
    """
    METIS k=2 bisection of a rook-connected block set.
    Returns a list of singlepart parts (>=2 if METIS output is non-contiguous).
    Falls back to the original list as a single part on error.
    """
    if len(geoids) < 2:
        return [geoids]

    geoid_set = set(geoids)
    idx = {g: i for i, g in enumerate(geoids)}
    adj_sets: list[set[int]] = [set() for _ in geoids]
    ew_raw: dict[tuple[int, int], float] = {}

    for a, b in adjacency:
        if a in geoid_set and b in geoid_set:
            ia, ib = idx[a], idx[b]
            adj_sets[ia].add(ib)
            adj_sets[ib].add(ia)
            d = max(0.01, _haversine_km(
                float(blocks_by_geoid[a]["lat"]), float(blocks_by_geoid[a]["lon"]),
                float(blocks_by_geoid[b]["lat"]), float(blocks_by_geoid[b]["lon"]),
            ))
            ew_raw[(min(ia, ib), max(ia, ib))] = 1.0 / d

    adj_lists = [sorted(s) for s in adj_sets]
    max_w = max(ew_raw.values(), default=1.0)
    scale = _HAVERSINE_SCALE / max_w
    ew: list[int] = []
    for i, nbrs in enumerate(adj_lists):
        for nb in nbrs:
            ew.append(max(1, int(ew_raw.get((min(i, nb), max(i, nb)), 1.0) * scale)))
    nw = [max(1, int(blocks_by_geoid[g]["pop"])) for g in geoids]

    try:
        _, membership = partition.partition(
            adj_lists, ew, nw, 2,
            ncuts=NCUTS_BISECT, niter=NITER_BISECT, recursive=True,
        )
    except Exception:
        return [geoids]

    part0 = [geoids[i] for i, m in enumerate(membership) if m == 0]
    part1 = [geoids[i] for i, m in enumerate(membership) if m == 1]

    # Split any non-contiguous METIS output into singlepart pieces
    result: list[list[str]] = []
    for part in (part0, part1):
        if part:
            result.extend(_rook_components(part, adjacency))
    return result if result else [geoids]


def _subcluster_pop(geoids: list[str], blocks_by_geoid: dict[str, dict]) -> int:
    return sum(int(blocks_by_geoid[g]["pop"]) for g in geoids)


def _subcluster_centroid(
    geoids: list[str],
    blocks_by_geoid: dict[str, dict],
) -> tuple[float, float]:
    """Population-weighted centroid using populated blocks only."""
    lat_sum = lon_sum = weight = 0.0
    for g in geoids:
        b = blocks_by_geoid[g]
        p = int(b["pop"])
        if p == 0:
            continue
        lat_sum += float(b["lat"]) * p
        lon_sum += float(b["lon"]) * p
        weight  += p
    if weight == 0:
        # fallback for all-zero-pop subclusters (should not arise in normal flow)
        lats = [float(blocks_by_geoid[g]["lat"]) for g in geoids]
        lons = [float(blocks_by_geoid[g]["lon"]) for g in geoids]
        return sum(lats) / len(lats), sum(lons) / len(lons)
    return lat_sum / weight, lon_sum / weight


def _build_subcluster_nodes(
    subclusters: dict[int, list[str]],
    blocks_by_geoid: dict[str, dict],
) -> list[dict]:
    nodes = []
    for sid in sorted(subclusters):
        geoids = subclusters[sid]
        pop = _subcluster_pop(geoids, blocks_by_geoid)
        lat, lon = _subcluster_centroid(geoids, blocks_by_geoid)
        nodes.append({"geoid": str(sid), "pop": pop, "lat": lat, "lon": lon})
    return nodes


def _build_subcluster_adj(
    subclusters: dict[int, list[str]],
    block_to_sub: dict[str, int],
    adjacency: set[tuple[str, str]],
) -> set[tuple[str, str]]:
    sub_adj: set[tuple[str, str]] = set()
    for a, b in adjacency:
        sa = block_to_sub.get(a)
        sb = block_to_sub.get(b)
        if sa is not None and sb is not None and sa != sb:
            sub_adj.add((str(min(sa, sb)), str(max(sa, sb))))
    return sub_adj


def _add_bridge_edges(
    blocks: list[dict],
    adjacency: set[tuple[str, str]],
    node_geoms: dict[str, object] | None = None,
) -> tuple[set[tuple[str, str]], int, list[int]]:
    """
    Detect disconnected components and add minimum bridge edges (Kruskal's on
    combined convex hull of per-component external nodes).  When node_geoms is
    provided, external nodes are those whose geometry touches the component's
    exterior ring (excluding bay-facing nodes).

    Returns (augmented_adjacency, n_bridges_added, comp_of).
    """
    from scipy.spatial import ConvexHull
    from shapely.ops import unary_union as _unary_union

    geoid_list = [b["geoid"] for b in blocks]
    idx_of = {g: i for i, g in enumerate(geoid_list)}
    n = len(geoid_list)

    adj: list[set[int]] = [set() for _ in range(n)]
    for a, b_g in adjacency:
        ia = idx_of.get(a)
        ib = idx_of.get(b_g)
        if ia is not None and ib is not None:
            adj[ia].add(ib)
            adj[ib].add(ia)

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
        print(f"  Subcluster graph fully connected ({n:,} nodes, 1 component).")
        return adjacency, 0, comp_of

    block_pops = [int(b["pop"]) for b in blocks]
    comp_pops = [sum(block_pops[i] for i in comp) for comp in components]
    print(f"  Subcluster graph has {k} components: "
          + ", ".join(f"{len(components[j]):,} nodes (pop {comp_pops[j]:,})"
                      for j in range(min(k, 5)))
          + ("..." if k > 5 else ""))

    def _iter_exteriors(geom):
        if geom.geom_type == "Polygon":
            yield geom.exterior
        elif geom.geom_type in ("MultiPolygon", "GeometryCollection"):
            for part in geom.geoms:
                yield from _iter_exteriors(part)

    ext_pts: list[tuple[float, float]] = []
    ext_comp: list[int] = []
    ext_block: list[int] = []

    for ci, comp in enumerate(components):
        hull_idx: list[int]
        comp_polys = (
            [node_geoms[blocks[bi]["geoid"]]
             for bi in comp
             if blocks[bi]["geoid"] in node_geoms]
            if node_geoms else []
        )
        if comp_polys:
            comp_union = _unary_union(comp_polys)
            exteriors = list(_iter_exteriors(comp_union))
            if exteriors:
                ext_ring_union = _unary_union(exteriors)
                hull_idx = [
                    bi for bi in comp
                    if node_geoms.get(blocks[bi]["geoid"]) is not None
                    and node_geoms[blocks[bi]["geoid"]].intersects(ext_ring_union)
                ]
                if not hull_idx:
                    hull_idx = comp
            else:
                hull_idx = comp
        elif len(comp) < 4:
            hull_idx = comp
        else:
            pts = np.column_stack([
                [float(blocks[i]["lat"]) for i in comp],
                [float(blocks[i]["lon"]) for i in comp],
            ])
            try:
                hull_idx = [comp[v] for v in ConvexHull(pts).vertices]
            except Exception:
                hull_idx = comp

        for bi in hull_idx:
            ext_pts.append((float(blocks[bi]["lat"]), float(blocks[bi]["lon"])))
            ext_comp.append(ci)
            ext_block.append(bi)

    candidates: list[tuple[float, int, int, int, int]] = []
    ext_arr = np.array(ext_pts)
    try:
        combined = ConvexHull(ext_arr)
        for simplex in combined.simplices:
            ai, aj = int(simplex[0]), int(simplex[1])
            ca, cj = ext_comp[ai], ext_comp[aj]
            if ca != cj:
                ba, bj = ext_block[ai], ext_block[aj]
                d = _haversine_km(
                    float(blocks[ba]["lat"]), float(blocks[ba]["lon"]),
                    float(blocks[bj]["lat"]), float(blocks[bj]["lon"]),
                )
                candidates.append((d, ba, bj, ca, cj))
    except Exception:
        pass
    candidates.sort()

    parent = list(range(k))

    def _find(x: int) -> int:
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    augmented = set(adjacency)
    n_bridges = 0

    for d, bi, bj, ca, cj in candidates:
        if _find(ca) != _find(cj):
            ga, gb = geoid_list[bi], geoid_list[bj]
            augmented.add((min(ga, gb), max(ga, gb)))
            parent[_find(ca)] = _find(cj)
            n_bridges += 1
        if n_bridges == k - 1:
            break

    if n_bridges < k - 1:
        fallback: list[tuple[float, int, int, int, int]] = []
        n_ext = len(ext_block)
        for ai in range(n_ext):
            for aj in range(ai + 1, n_ext):
                ca, cj = ext_comp[ai], ext_comp[aj]
                if _find(ca) == _find(cj):
                    continue
                lat_a, lon_a = ext_pts[ai]
                lat_j, lon_j = ext_pts[aj]
                fallback.append((
                    _haversine_km(lat_a, lon_a, lat_j, lon_j),
                    ext_block[ai], ext_block[aj], ca, cj,
                ))
        fallback.sort()
        for d, bi, bj, ca, cj in fallback:
            if _find(ca) != _find(cj):
                ga, gb = geoid_list[bi], geoid_list[bj]
                augmented.add((min(ga, gb), max(ga, gb)))
                parent[_find(ca)] = _find(cj)
                n_bridges += 1
            if n_bridges == k - 1:
                break

    print(f"  Added {n_bridges} bridge edge(s). Subcluster graph now fully connected.")
    return augmented, n_bridges, comp_of


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
    part_geoms: dict[int, object] = {p: unary_union(gs) for p, gs in groups.items()}

    n_parts = len(part_geoms)
    cmap = plt.cm.get_cmap("tab20")
    part_colour = {p: mcolors.to_hex(cmap((i % 20) / 20))
                   for i, p in enumerate(sorted(part_geoms))}

    rows = []
    for p, geom in part_geoms.items():
        rows.append({"geometry": geom, "colour": part_colour[p]})
    gdf = gpd.GeoDataFrame(rows, crs="EPSG:4269")
    fig, ax = plt.subplots(figsize=(12, 10))
    gdf.plot(ax=ax, color=gdf["colour"], edgecolor="white", linewidth=0.3, alpha=0.85)
    ax.set_title(title, fontsize=12)
    ax.set_xlabel("Longitude")
    ax.set_ylabel("Latitude")
    ax.set_aspect("equal")
    plt.tight_layout()
    plt.show(block=False)
    plt.pause(0.1)
    return part_geoms


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
        blocks_by_geoid = {b["geoid"]: b for b in blocks}
        total_pop = sum(int(b["pop"]) for b in blocks)
        print(f"  {len(blocks):,} blocks, total pop {total_pop:,}")

        # --- Block geometries ---
        print("\nFetching block geometries...")
        block_geoms_wkb = db.fetch_block_geoms(conn, statefp)
        print(f"  {len(block_geoms_wkb):,} block polygons loaded.")

        # --- Rook adjacency ---
        geoids = [b["geoid"] for b in blocks]
        print("\nChecking rook adjacency cache...")
        missing = db.get_missing_adjacency_geoids(conn, "blocks", geoids)
        if missing:
            print(f"  Computing adjacency for {len(missing):,} blocks...")
            n_pairs = db.compute_and_store_adjacency_bulk(conn, "blocks", statefp)
            print(f"  {n_pairs:,} new pairs inserted.")
        else:
            print("  Adjacency fully cached.")
        rook_adjacency = db.fetch_adjacency(conn, "blocks", geoids)
        print(f"  {len(rook_adjacency):,} rook-adjacent block pairs.")

        # --- Initial subclusters: one per tract (populated blocks only) ---
        populated_blocks = [b for b in blocks if int(b["pop"]) > 0]
        zero_pop_blocks  = [b for b in blocks if int(b["pop"]) == 0]
        print(f"\nGrouping populated blocks by tract "
              f"({len(populated_blocks):,} populated, {len(zero_pop_blocks):,} zero-pop)...")
        tract_to_blocks: dict[str, list[str]] = {}
        for b in populated_blocks:
            tract_to_blocks.setdefault(_tract_id(b["geoid"]), []).append(b["geoid"])

        n_populated_tracts = sum(
            1 for gs in tract_to_blocks.values()
            if any(int(blocks_by_geoid[g]["pop"]) > 0 for g in gs)
        )
        avg_tract_pop = total_pop / max(n_populated_tracts, 1)
        threshold = 2.0 * avg_tract_pop
        print(f"  {len(tract_to_blocks):,} tracts  "
              f"({n_populated_tracts:,} populated)")
        print(f"  avg tract pop = {avg_tract_pop:,.0f}  "
              f"threshold = {threshold:,.0f}")

        subclusters: dict[int, list[str]] = {}
        block_to_sub: dict[str, int] = {}
        next_id = 0
        for tract_geoids in tract_to_blocks.values():
            subclusters[next_id] = tract_geoids
            for g in tract_geoids:
                block_to_sub[g] = next_id
            next_id += 1

        # --- Iterative bisection ---
        print("\nIterative bisection until no subcluster exceeds threshold...")
        iteration = 0
        while True:
            oversized = [
                sid for sid, gs in subclusters.items()
                if _subcluster_pop(gs, blocks_by_geoid) > threshold
            ]
            if not oversized:
                break
            iteration += 1
            n_split = 0
            for sid in oversized:
                geoids_s = subclusters.pop(sid)
                comps = _rook_components(geoids_s, rook_adjacency)
                for comp in comps:
                    comp_pop = _subcluster_pop(comp, blocks_by_geoid)
                    if comp_pop > threshold and len(comp) > 1:
                        parts = _bisect(comp, rook_adjacency, blocks_by_geoid)
                        n_split += len(parts) - 1
                    else:
                        parts = [comp]
                    for part in parts:
                        subclusters[next_id] = part
                        for g in part:
                            block_to_sub[g] = next_id
                        next_id += 1

            print(f"  Iter {iteration}: {len(oversized)} oversized → "
                  f"{n_split} additional splits → {len(subclusters):,} subclusters total")

        print(f"  Done: {len(subclusters):,} subclusters after {iteration} iteration(s).")

        # --- Assign zero-pop blocks by BFS wave-front from populated subclusters ---
        if zero_pop_blocks:
            zero_pop_geoids = {b["geoid"] for b in zero_pop_blocks}
            zp_by_geoid = {b["geoid"]: b for b in zero_pop_blocks}

            # Neighbour lookup restricted to rook adjacency
            zp_nbrs: dict[str, list[str]] = {g: [] for g in zero_pop_geoids}
            for a, b_g in rook_adjacency:
                if a in zero_pop_geoids:
                    zp_nbrs[a].append(b_g)
                if b_g in zero_pop_geoids:
                    zp_nbrs[b_g].append(a)

            # Population-weighted centroids of final subclusters (populated blocks only)
            sub_cent: dict[int, tuple[float, float]] = {
                sid: _subcluster_centroid(geoids, blocks_by_geoid)
                for sid, geoids in subclusters.items()
            }

            # Wave-front: assign each zero-pop block to the nearest adjacent subcluster
            unassigned_zp: set[str] = set(zero_pop_geoids)
            changed = True
            while changed:
                changed = False
                for geoid in list(unassigned_zp):
                    adj_subs = {
                        block_to_sub[nb]
                        for nb in zp_nbrs[geoid]
                        if nb in block_to_sub
                    }
                    if not adj_subs:
                        continue
                    b = zp_by_geoid[geoid]
                    lat, lon = float(b["lat"]), float(b["lon"])
                    best = min(
                        adj_subs,
                        key=lambda s: (
                            (sub_cent[s][0] - lat) ** 2 + (sub_cent[s][1] - lon) ** 2
                            if s in sub_cent else float("inf")
                        ),
                    )
                    block_to_sub[geoid] = best
                    subclusters[best].append(geoid)
                    unassigned_zp.discard(geoid)
                    changed = True

            # Isolated zero-pop blocks: group contiguous components into new clusters
            if unassigned_zp:
                n_new = 0
                visited_zp: set[str] = set()
                for start in unassigned_zp:
                    if start in visited_zp:
                        continue
                    comp: list[str] = []
                    queue = [start]
                    visited_zp.add(start)
                    while queue:
                        node = queue.pop()
                        comp.append(node)
                        for nb in zp_nbrs[node]:
                            if nb in unassigned_zp and nb not in visited_zp:
                                visited_zp.add(nb)
                                queue.append(nb)
                    subclusters[next_id] = comp
                    for geoid in comp:
                        block_to_sub[geoid] = next_id
                    next_id += 1
                    n_new += 1
                print(f"  {len(unassigned_zp):,} isolated zero-pop blocks "
                      f"→ {n_new} new cluster(s).")

            n_zp_assigned = len(zero_pop_geoids) - len(unassigned_zp) \
                if unassigned_zp else len(zero_pop_geoids)
            print(f"  {n_zp_assigned:,} zero-pop blocks assigned by rook adjacency.")

        # --- Visualise subclusters ---
        geoid_to_sub_viz = {g: block_to_sub[g] for g in blocks_by_geoid}
        sub_geoms = _show_partition(
            geoid_to_sub_viz, block_geoms_wkb,
            title=f"{state_name} — tract subclusters: {len(subclusters):,}",
        )

        # --- Build subcluster nodes and adjacency ---
        sub_nodes = _build_subcluster_nodes(subclusters, blocks_by_geoid)
        sub_adj = _build_subcluster_adj(subclusters, block_to_sub, rook_adjacency)

        # --- Bridge disconnected subcluster components ---
        print("\nChecking subcluster graph connectivity...")
        sub_adj, n_bridges, sub_comp_of = _add_bridge_edges(
            sub_nodes, sub_adj,
            node_geoms={str(sid): sub_geoms[sid]
                        for sid in subclusters if sid in sub_geoms},
        )

        # --- Final METIS: equal edge weights, k = n_districts ---
        print(f"\nFinal METIS: {len(sub_nodes):,} subclusters → {n_districts} districts...")
        sub_id_list = [int(n["geoid"]) for n in sub_nodes]
        sub_idx = {n["geoid"]: i for i, n in enumerate(sub_nodes)}

        # Build adj lists with uniform edge weight
        adj_sets: list[set[int]] = [set() for _ in sub_nodes]
        for ga, gb in sub_adj:
            ia = sub_idx.get(ga)
            ib = sub_idx.get(gb)
            if ia is not None and ib is not None:
                adj_sets[ia].add(ib)
                adj_sets[ib].add(ia)
        adj_lists = [sorted(s) for s in adj_sets]
        ew = [1] * sum(len(s) for s in adj_sets)   # uniform edge weights
        nw = [max(1, int(n["pop"])) for n in sub_nodes]

        _, membership = partition.partition(
            adj_lists, ew, nw, n_districts,
            ncuts=NCUTS_FINAL, niter=NITER_FINAL,
        )
        print(f"  Done.")

        # --- Disaggregate: block → subcluster → district ---
        sub_to_district = {sub_id_list[i]: membership[i]
                           for i in range(len(sub_nodes))}
        block_to_district: dict[str, int] = {
            g: sub_to_district[block_to_sub[g]]
            for g in blocks_by_geoid
        }

        # --- Population stats ---
        pop_per_district: dict[int, int] = {}
        for b in blocks:
            if int(b["pop"]) > 0:
                d = block_to_district[b["geoid"]]
                pop_per_district[d] = pop_per_district.get(d, 0) + int(b["pop"])
        ideal = total_pop / n_districts
        print("\nDistrict populations:")
        for d in sorted(pop_per_district):
            pct = 100 * (pop_per_district[d] - ideal) / ideal
            print(f"  District {d}: {pop_per_district[d]:,}  ({pct:+.1f}%)")
        worst = max(abs(100 * (p - ideal) / ideal) for p in pop_per_district.values())

        # --- Visualise districts ---
        import geopandas as gpd
        import matplotlib.pyplot as plt
        import matplotlib.patches as mpatches
        import matplotlib.colors as mcolors
        from shapely.ops import unary_union as _uu

        dist_colour = {
            d: mcolors.to_hex(plt.cm.get_cmap("tab10")(i / max(n_districts - 1, 1)))
            for i, d in enumerate(sorted(pop_per_district))
        }
        dist_geoms: dict[int, list] = {}
        for sid, geom in sub_geoms.items():
            d = sub_to_district.get(sid)
            if d is not None:
                dist_geoms.setdefault(d, []).append(geom)
        rows2 = [{"geometry": _uu(gs), "colour": dist_colour[d]}
                 for d, gs in dist_geoms.items()]
        gdf2 = gpd.GeoDataFrame(rows2, crs="EPSG:4269")
        fig2, ax2 = plt.subplots(figsize=(12, 10))
        gdf2.plot(ax=ax2, color=gdf2["colour"], edgecolor="white", linewidth=0.5, alpha=0.9)
        legend_elems = [
            mpatches.Patch(facecolor=dist_colour[d], label=f"District {d}")
            for d in sorted(dist_colour)
        ]
        ax2.legend(handles=legend_elems, loc="lower left", fontsize=9)
        ax2.set_title(f"{state_name} — {n_districts} districts", fontsize=12)
        ax2.set_xlabel("Longitude")
        ax2.set_ylabel("Latitude")
        ax2.set_aspect("equal")
        plt.tight_layout()
        plt.show(block=False)
        plt.pause(0.1)

        # --- District geometries (union subcluster geoms by district) ---
        district_geoms_wkt: dict[int, tuple[str, int]] = {
            d: (_uu(gs).wkt, pop_per_district.get(d, 0))
            for d, gs in dist_geoms.items()
        }

        # --- Write to DB ---
        params = {
            "status":          "complete",
            "method":          "tract_metis",
            "n_tracts":        len(tract_to_blocks),
            "n_populated_tracts": n_populated_tracts,
            "avg_tract_pop":   int(avg_tract_pop),
            "threshold":       int(threshold),
            "n_subclusters":   len(subclusters),
            "bisect_iters":    iteration,
            "n_bridges":       n_bridges,
            "n_blocks":        len(blocks),
            "n_edges":         len(sub_adj),
            "districts":       n_districts,
            "ideal_pop":       int(ideal),
            "worst_deviation": round(worst, 2),
            "ncuts_final":     NCUTS_FINAL,
            "niter_final":     NITER_FINAL,
        }
        run_id = db.write_run(conn, "blocks", statefp, n_districts, params)
        print(f"\nRun ID: {run_id}")

        db.write_assignments(conn, run_id, block_to_district)
        print(f"  {len(block_to_district):,} block assignments written.")

        db.write_district_geoms_wkt(conn, run_id, district_geoms_wkt)
        print("  District geometries written.")

        # --- GeoJSON export ---
        slug = state_name.lower().replace(" ", "_")
        geojson_path = os.path.join(OUTPUT_DIR, f"{slug}_tract_run{run_id}.geojson")
        db.export_geojson(conn, run_id, geojson_path)
        print(f"  GeoJSON → {geojson_path}")

        # --- Subcluster components GeoJSON ---
        import shapely as _shapely
        comp_features = []
        for i, node in enumerate(sub_nodes):
            sid = int(node["geoid"])
            geom = sub_geoms.get(sid)
            if geom is None:
                continue
            comp_features.append({
                "type": "Feature",
                "geometry": json.loads(_shapely.to_geojson(geom)),
                "properties": {
                    "subcluster_id": sid,
                    "component_id": sub_comp_of[i],
                    "pop": node["pop"],
                    "district": sub_to_district.get(sid),
                },
            })
        comp_path = os.path.join(OUTPUT_DIR, f"{slug}_tract_run{run_id}_subclusters.geojson")
        with open(comp_path, "w") as fh:
            json.dump({"type": "FeatureCollection", "features": comp_features}, fh)
        print(f"  Subclusters GeoJSON → {comp_path}")

        # --- Deviation log ---
        log_path = os.path.join(OUTPUT_DIR, f"{slug}_tract_run{run_id}_deviations.md")
        lines = [
            f"# {state_name} — tract METIS — {n_districts} districts — Run {run_id}",
            "",
            f"method=tract_metis  threshold={int(threshold)}  "
            f"subclusters={len(subclusters)}  bisect_iters={iteration}  "
            f"ncuts={NCUTS_FINAL}  niter={NITER_FINAL}",
            f"blocks={len(blocks):,}  tracts={len(tract_to_blocks):,}  "
            f"districts={n_districts}  ideal={ideal:,.0f}  worst_deviation={worst:.1f}%",
            "",
            "## District populations",
            "",
            f"{'District':>10}  {'Population':>12}  {'Deviation':>10}",
        ]
        for d in sorted(pop_per_district):
            pct = 100 * (pop_per_district[d] - ideal) / ideal
            lines.append(f"{d:>10}  {pop_per_district[d]:>12,}  {pct:>+10.2f}%")
        with open(log_path, "w") as fh:
            fh.write("\n".join(lines) + "\n")
        print(f"  Deviation log → {log_path}")

        input("\nPress Enter to close plots...")

    finally:
        conn.close()


if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: run_tract_metis.py <statefp> <n_districts>")
        sys.exit(1)
    main(sys.argv[1], int(sys.argv[2]))
