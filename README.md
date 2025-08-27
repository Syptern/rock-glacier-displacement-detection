# Rock Glacier Displacement Detector

A Python tool to **detect and cluster surface displacements** between two georeferenced raster images (e.g., aerial orthophotos, satellite scenes).  
Built for **rock glacier monitoring**, but applicable to any displacement detection task.  
Inspired by **SAGA IMCORR** with added clustering and filtering.

---

## Features

- Compute **local displacement vectors** (per window) between two rasters.  
- Export **points** (window centers) and **lines** (displacement vectors) as GeoJSON.  
- **Cluster** displacements by direction (circular variance) and local neighbor length consistency.  
- Optional **post-filter** of messy clusters by length standard deviation.  
- Results open directly in **QGIS** (style by `cluster_id`).

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