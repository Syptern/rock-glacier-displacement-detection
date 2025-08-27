import numpy as np
import geopandas as gpd
from shapely.geometry import LineString

# ---------------------------
# Circular statistics function
# ---------------------------
def circular_variance(angles_deg):
    """Circular variance of angles in degrees"""
    a = np.deg2rad(angles_deg)
    angles_complex = np.exp(1j * a)
    r = abs(np.sum(angles_complex)) / len(angles_deg)
    return 1 - r

# ---------------------------
# Clustering function
# ---------------------------
def cluster_displacement_points(points_gdf, window_size=16,
                                             min_cluster_size=3, std_threshold=0.2,
                                             disp_threshold=1.0, length_tol_rel=2):
    """
    Cluster points with similar direction and local length consistency.
    
    Parameters
    ----------
    points_gdf : GeoDataFrame
        Must have 'dx', 'dy', 'length', 'angle'
    window_size : int
        Grid step size (map units) for neighbor search
    min_cluster_size : int
        Minimum number of points in a cluster
    std_threshold : float
        Maximum circular variance of cluster angles
    disp_threshold : float
        Minimum displacement length to consider
    length_tol_rel : float
        Maximum allowed relative length difference to direct neighbor
    """
    # Filter points below displacement threshold
    gdf = points_gdf[points_gdf['length'] > disp_threshold].copy()
    gdf = gdf.reset_index(drop=True)

    # Build 2D grid index for neighbor lookup
    step = window_size
    grid_index = {}
    coords = np.array(list(zip(gdf.geometry.x, gdf.geometry.y)))
    grid_coords = np.floor(coords / step).astype(int)
    for idx, (gx, gy) in enumerate(grid_coords):
        grid_index.setdefault((gx, gy), []).append(idx)

    # Prepare cluster tracking
    cluster_ids = np.full(len(gdf), -1, dtype=int)
    cluster_id = 0
    visited = set()

    offsets = [(-1,-1), (-1,0), (-1,1), (0,-1), (0,1), (1,-1), (1,0), (1,1)]

    for i in range(len(gdf)):
        if i in visited:
            continue

        # Seed cluster
        stack = [(i, None)]  # (point_idx, parent_idx)
        cluster_points = []
        cluster_angles = []

        while stack:
            idx, parent_idx = stack.pop()
            if idx in visited:
                continue

            candidate_angle = gdf.loc[idx, 'angle']
            candidate_length = gdf.loc[idx, 'length']

            # Always accept first two points
            if len(cluster_points) >= 2:
                tentative_angles = cluster_angles + [candidate_angle]
                if circular_variance(tentative_angles) > std_threshold:
                    continue

                # Check length relative to parent
                if parent_idx is not None:
                    parent_length = gdf.loc[parent_idx, 'length']
                    if abs(candidate_length - parent_length) > parent_length * length_tol_rel:
                        continue

            # Accept point
            visited.add(idx)
            cluster_points.append(idx)
            cluster_angles.append(candidate_angle)

            # Add neighbors to stack
            gx, gy = grid_coords[idx]
            for ox, oy in offsets:
                neighbor_grid = (gx + ox, gy + oy)
                if neighbor_grid in grid_index:
                    for n_idx in grid_index[neighbor_grid]:
                        if n_idx in visited:
                            continue
                        ang_diff = abs(gdf.loc[n_idx, 'angle'] - candidate_angle)
                        ang_diff = min(ang_diff, 360 - ang_diff)
                        if ang_diff < 50:  # broad directional similarity
                            stack.append((n_idx, idx))

        # Assign cluster_id if cluster large enough
        if len(cluster_points) >= min_cluster_size:
            for idx in cluster_points:
                cluster_ids[idx] = cluster_id
            cluster_id += 1

    gdf['cluster_id'] = cluster_ids

    # Points GeoDataFrame
    points_gdf_clustered = gdf[gdf['cluster_id'] >= 0].copy()

    # Generate lines
    lines = []
    for idx, row in points_gdf_clustered.iterrows():
        x0, y0 = row.geometry.x, row.geometry.y
        x1, y1 = x0 + row.dx, y0 + row.dy
        lines.append(LineString([(x0, y0), (x1, y1)]))

    lines_gdf_clustered = gpd.GeoDataFrame(points_gdf_clustered.drop(columns='geometry'),
                                           geometry=lines,
                                           crs=points_gdf_clustered.crs)

    return points_gdf_clustered, lines_gdf_clustered
