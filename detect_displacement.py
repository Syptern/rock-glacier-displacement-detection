import argparse
import rasterio
import numpy as np
from skimage.feature import match_template
import geopandas as gpd
from shapely.geometry import Point, LineString
import math
import os
    
from directional_clustering import cluster_displacement_points

def local_displacement_points_and_lines(img1, img2, transform, crs, window, search, grid_size):
    h, w = img1.shape
    point_records = []
    line_records = []

    px_x = transform.a
    px_y = -transform.e  

    for y in range(0, h - window, grid_size):
        for x in range(0, w - window, grid_size):
            patch1 = img1[y:y+window, x:x+window]

            # define search area
            y0, x0 = max(0, y-search), max(0, x-search)
            y1, x1 = min(h, y+window+search), min(w, x+window+search)
            search_area = img2[y0:y1, x0:x1]

            if search_area.shape[0] < window or search_area.shape[1] < window:
                continue

            # normalized cross-correlation
            result = match_template(search_area, patch1, pad_input=True)
            ij = np.unravel_index(np.argmax(result), result.shape)
            peak_y, peak_x = ij

            # displacement in pixels (corrected)
            dy = (peak_y + y0) - y - (window // 2)
            dx = (peak_x + x0) - x - (window // 2)

            # convert to map units
            dx_map = dx * px_x
            dy_map = -dy * px_y

            # center coordinates
            cx, cy = x + window // 2, y + window // 2
            X, Y = transform * (cx, cy)

            # line endpoint
            X2 = X + dx_map
            Y2 = Y + dy_map

            # magnitude and angles
            length = math.sqrt(dx_map**2 + dy_map**2)
            angle_rad = math.atan2(dy_map, dx_map)
            angle_deg = math.degrees(angle_rad)
            bearing = (90 - angle_deg) % 360

            # point record
            point_records.append({
                "geometry": Point(X, Y),
                "dx": dx_map,
                "dy": dy_map,
                "length": length,
                "angle": angle_deg,
                "bearing": bearing,
                "corr": result[peak_y, peak_x]
            })

            # line record
            line_records.append({
                "geometry": LineString([(X, Y), (X2, Y2)]),
                "dx": dx_map,
                "dy": dy_map,
                "length": length,
                "angle": angle_deg,
                "bearing": bearing,
                "corr": result[peak_y, peak_x]
            })

    points_gdf = gpd.GeoDataFrame(point_records, crs=crs)
    lines_gdf = gpd.GeoDataFrame(line_records, crs=crs)
    return points_gdf, lines_gdf


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Detect rock glacier displacements from two rasters.")
    parser.add_argument("before", help="Path to the first raster (before).")
    parser.add_argument("after", help="Path to the second raster (after).")
    parser.add_argument("--window_size", type=int, default=16, help="Template (window) size for displacement search.")
    parser.add_argument("--search_size", type=int, default=32, help="Search window size.")
    parser.add_argument("--min_cluster_size", type=int, default=5, help="Minimum number of points per cluster.")
    parser.add_argument("--std_threshold", type=float, default=0.05, help="Maximum circular variance for angles.")
    parser.add_argument("--disp_threshold", type=float, default=1, help="Minimum displacement magnitude to consider.")
    parser.add_argument("--length_tol_rel", type=float, default=0.2, help="Relative tolerance for length differences.")
    parser.add_argument("--outdir", type=str, default="results", help="Output directory.")
    parser.add_argument("--grid_size", type=int, default=16, help="Grid size of displacement vectors")
    args = parser.parse_args()
    
    raster1 = args.before
    raster2 = args.after

    with rasterio.open(raster1) as src1, rasterio.open(raster2) as src2:
        img1 = src1.read(1).astype(np.float32)
        img2 = src2.read(1).astype(np.float32)
        transform = src1.transform
        crs = src1.crs

    points_gdf, lines_gdf = local_displacement_points_and_lines(
        img1, img2, transform, crs, window=args.window_size, search=args.search_size, grid_size=args.grid_size
    )

    points_clustered, lines_clustered = cluster_displacement_points(
        points_gdf,
        window_size=args.window_size,
        min_cluster_size=args.min_cluster_size,
        std_threshold=args.std_threshold,
        disp_threshold=args.disp_threshold,
        length_tol_rel=args.length_tol_rel
    )
    
    outdir = args.outdir
    os.makedirs(outdir, exist_ok=True)
        
    out_points = os.path.join(outdir, "displacement_points.geojson")
    out_lines = os.path.join(outdir, "displacement_lines.geojson")
    
    out_clustered_points = os.path.join(outdir, "clustered_displacement_points.geojson")
    out_clustered_lines = os.path.join(outdir,"clustered_displacement_lines.geojson")
    
    points_gdf.to_file(out_points, driver="GeoJSON")
    lines_gdf.to_file(out_lines, driver="GeoJSON")

    print(f"✅ Displacement points written to {out_points}")
    print(f"✅ Displacement lines written to {out_lines}")
    
    
    points_clustered.to_file(out_clustered_points, driver="GeoJSON")
    lines_clustered.to_file(out_clustered_lines, driver="GeoJSON")
    
    print(f"✅ Displacement clustered points written to {out_points}")
    print(f"✅ Displacement clustered lines written to {out_lines}")
