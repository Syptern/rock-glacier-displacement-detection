# Rock Glacier Displacement Detector

A Python tool to **detect and cluster surface displacements** between two georeferenced raster images (e.g., aerial orthophotos, satellite scenes). From my experience DEM derived hillshades cause the least artifacts.
Built for **rock glacier monitoring**, but applicable to any displacement detection task.  
Inspired by **SAGA IMCORR** with added detection of similar movement clusters.

---

## Features

- Compute **local displacement vectors** (per window) between two rasters.  
- Export **points** (window centers) and **lines** (displacement vectors) as GeoJSON.  
- **Cluster** displacements by direction (circular variance) and local neighbor length consistency.  

---

## Installation

Clone the repository and install dependencies:

```bash
git clone https://github.com/yourusername/rock-glacier-displacement.git
cd rock-glacier-displacement
pip install -r requirements.txt
```

---

## Usage

```
python detect_displacement.py before.tif after.tif --window_size 16 --search_size 32 --min_cluster_size 10 --std_threshold 0.2 --disp_threshold 1.0
```

---

## Options

| Flag | Type | Default | Description |
| --- | --- | --- | --- |
| `before` | path | —   | Path to the first (earlier) raster. |
| `after` | path | —   | Path to the second (later) raster. |
| `--window_size` | int | `16` | Template (block) size in pixels for matching. |
| `--search_size` | int | `32` | Half-size (radius) of search window in pixels around each template. |
| `--min_cluster_size` | int | `5` | Minimum number of points to keep a cluster. |
| `--std_threshold` | float | `0.2` | Max **circular variance** of angles per cluster (lower = stricter). |
| `--disp_threshold` | float | `1.0` | Minimum displacement magnitude to consider (map units). |
| `--length_tol_rel` | float | `0.2` | **Local** relative length tolerance vs. neighbor (e.g., `0.2` = ±20%). |
| `--max_length_std` | float | `1.0` | Post-filter: maximum allowed **std** of lengths within a cluster. |
| `--outdir` | path | `results/` | Output directory for GeoJSONs. |

---

## Output examples

### Displacement map

![example displacement vectors output](images\example_displacement_output.png)

### Selected clusters

![example clusters output](images\example_displacement_clusters_output.png)