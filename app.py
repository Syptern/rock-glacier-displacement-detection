import rasterio
import numpy as np
from skimage.feature import match_template
import geopandas as gpd
from shapely.geometry import Point, LineString
import math

from directional_clustering import cluster_displacement_points

def local_displacement_points_and_lines(img1, img2, transform, crs, window=16, search=32):
    h, w = img1.shape
    point_records = []
    line_records = []

    # pixel size in map units
    px_x = transform.a
    px_y = -transform.e  # usually negative in north-up rasters

    for y in range(0, h - window, window):
        for x in range(0, w - window, window):
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
    raster1 = "before.tif"
    raster2 = "after.tif"
    out_points = "displacement_points.geojson"
    out_lines = "displacement_lines.geojson"

    with rasterio.open(raster1) as src1, rasterio.open(raster2) as src2:
        img1 = src1.read(1).astype(np.float32)
        img2 = src2.read(1).astype(np.float32)
        transform = src1.transform
        crs = src1.crs

    points_gdf, lines_gdf = local_displacement_points_and_lines(
        img1, img2, transform, crs, window=16, search=32
    )

    points_clustered, lines_clustered = cluster_displacement_points(
        points_gdf,
        window_size=16,
        min_cluster_size=15,
        std_threshold=0.05,
        disp_threshold=0.8
    )
    
    points_gdf.to_file(out_points, driver="GeoJSON")
    lines_gdf.to_file(out_lines, driver="GeoJSON")

    print(f"✅ Displacement points written to {out_points}")
    print(f"✅ Displacement lines written to {out_lines}")
    
    
    points_clustered.to_file("clustered_points.geojson", driver="GeoJSON")
    lines_clustered.to_file("clustered_lines.geojson", driver="GeoJSON")
    
    print(f"✅ Displacement clustered points written to {out_points}")
    print(f"✅ Displacement clustered lines written to {out_lines}")
