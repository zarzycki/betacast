#!/usr/bin/env python3
"""
Generate a DATM domain file from any regular rectilinear lat/lon grid.

Reads 1-D lat/lon arrays from a source netCDF file and writes a
DATM-compatible domain file containing xc, yc, area, and mask.

Areas are computed analytically in steradians (radians^2), matching the
NCL function area_global_rectilinear_grid with opt@rearth=1.0.

Usage:
    python generate-datm-domain.py --infile=source.nc --outfile=domain.nc [options]

Examples:
    # ERA5 (0.25-degree, 721x1440; latitude/longitude coords, N->S — auto-flipped)
    python generate-datm-domain.py \
        --infile=/glade/campaign/collections/rda/data/ds633.0/e5.oper.an.sfc/198608/e5.oper.an.sfc.128_165_10u.ll025sc.1986080100_1986083123.nc \
        --outfile=era5-domain.nc \
        --lat_name=latitude --lon_name=longitude

    # CESM LENS (192x288 f09 grid; standard lat/lon coords, S->N)
    python generate-datm-domain.py \
        --infile=/glade/u/home/zarzycki/scratch/CESM_LENS_temp/TREFHT/ens_TREFHT_anom.nc \
        --outfile=lens-domain.nc

    # CR20V3 (256x512 Gaussian ~0.7-degree; latitude/longitude coords, N->S — auto-flipped)
    python generate-datm-domain.py \
        --infile=/glade/u/home/zarzycki/rda/d131003/anl/anl_mean_1901_TMP_2m.nc \
        --outfile=cr20v3-domain.nc \
        --lat_name=latitude --lon_name=longitude
"""

import numpy as np
import xarray as xr
import os
import sys
import argparse
import logging
from datetime import datetime

logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    datefmt='%Y-%m-%d %H:%M:%S',
)
logger = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
# Argument parsing
# ---------------------------------------------------------------------------
parser = argparse.ArgumentParser(description="Generate a DATM domain file from a rectilinear lat/lon grid")
parser.add_argument("--infile",    required=True,
                    help="Source netCDF file containing 1-D lat and lon arrays")
parser.add_argument("--outfile",   required=True,
                    help="Output path for the domain file")
parser.add_argument("--lat_name",  default="lat",
                    help="Name of the latitude variable in --infile (default: lat)")
parser.add_argument("--lon_name",  default="lon",
                    help="Name of the longitude variable in --infile (default: lon)")
parser.add_argument("--title",     default=None,
                    help="Title attribute for the output file (auto-generated if omitted)")
args = parser.parse_args()

# ---------------------------------------------------------------------------
# Read lat/lon
# ---------------------------------------------------------------------------
if not os.path.exists(args.infile):
    logger.error(f"Input file not found: {args.infile}")
    sys.exit(1)

logger.info(f"Reading grid from: {args.infile}")
ds = xr.open_dataset(args.infile, decode_times=False)

if args.lat_name not in ds:
    logger.error(f"Latitude variable '{args.lat_name}' not found. Available: {list(ds.coords) + list(ds.data_vars)}")
    sys.exit(1)
if args.lon_name not in ds:
    logger.error(f"Longitude variable '{args.lon_name}' not found. Available: {list(ds.coords) + list(ds.data_vars)}")
    sys.exit(1)

lat = ds[args.lat_name].values.squeeze().astype(np.float64)
lon = ds[args.lon_name].values.squeeze().astype(np.float64)
ds.close()

if lat.ndim != 1 or lon.ndim != 1:
    logger.error(f"Expected 1-D lat/lon arrays, got shapes lat={lat.shape}, lon={lon.shape}")
    sys.exit(1)

# Ensure S->N ordering (DATM convention)
if lat[0] > lat[-1]:
    logger.info("Latitude is N->S; flipping to S->N")
    lat = lat[::-1]

nlat, nlon = len(lat), len(lon)
logger.info(f"Grid: {nlat} lat x {nlon} lon")
logger.info(f"  lat: {lat[0]:.4f} -> {lat[-1]:.4f}")
logger.info(f"  lon: {lon[0]:.4f} -> {lon[-1]:.4f}")

# ---------------------------------------------------------------------------
# Compute cell areas in steradians (radians^2)
# Matches NCL area_global_rectilinear_grid with opt@rearth=1.0
#
# For a regular rectilinear grid the area of cell (j, i) is:
#   area = dlon_i * (sin(lat_north_j) - sin(lat_south_j))
# where lat bounds are midpoints between adjacent cell centres
# (clamped to ±pi/2 at the poles).
# ---------------------------------------------------------------------------
lat_r = np.deg2rad(lat)
lon_r = np.deg2rad(lon)

# Latitude bounds: midpoints between adjacent centres, poles at ±pi/2
lat_bnds       = np.empty(nlat + 1)
lat_bnds[1:-1] = 0.5 * (lat_r[:-1] + lat_r[1:])
lat_bnds[0]    = -np.pi / 2.0
lat_bnds[-1]   =  np.pi / 2.0

# Longitude widths: midpoints between adjacent centres, wrap at ends
lon_bnds       = np.empty(nlon + 1)
lon_bnds[1:-1] = 0.5 * (lon_r[:-1] + lon_r[1:])
lon_bnds[0]    = lon_r[0]  - 0.5 * (lon_r[1] - lon_r[0])
lon_bnds[-1]   = lon_r[-1] + 0.5 * (lon_r[-1] - lon_r[-2])
dlon = lon_bnds[1:] - lon_bnds[:-1]          # (nlon,)

dlat_sin = np.sin(lat_bnds[1:]) - np.sin(lat_bnds[:-1])   # (nlat,)

# Broadcast to 2-D: area[j, i] = dlat_sin[j] * dlon[i]
area_2d = (dlat_sin[:, None] * dlon[None, :]).astype(np.float32)

logger.info(f"Area sum: {area_2d.sum():.6f} sr  (full sphere = {4*np.pi:.6f} sr)")

# ---------------------------------------------------------------------------
# Build output dataset
# ---------------------------------------------------------------------------
lon_2d = np.tile(lon,          (nlat, 1)).astype(np.float32)
lat_2d = np.tile(lat[:, None], (1, nlon)).astype(np.float32)
mask_2d = np.ones((nlat, nlon), dtype=np.float32)

ds_out = xr.Dataset({
    'xc':   xr.DataArray(lon_2d,  dims=['nj', 'ni'],
                         attrs={'units': 'degrees_east',  'long_name': 'longitude',
                                'mode': 'time-invariant'}),
    'yc':   xr.DataArray(lat_2d,  dims=['nj', 'ni'],
                         attrs={'units': 'degrees_north', 'long_name': 'latitude',
                                'mode': 'time-invariant'}),
    'area': xr.DataArray(area_2d, dims=['nj', 'ni'],
                         attrs={'units': 'radians^2',
                                'long_name': 'area of grid cell in radians squared',
                                'mode': 'time-invariant'}),
    'mask': xr.DataArray(mask_2d, dims=['nj', 'ni'],
                         attrs={'units': 'unitless', 'long_name': 'domain mask',
                                'mode': 'time-invariant'}),
})

title = args.title or f"DATM domain ({nlat}x{nlon}) generated from {os.path.basename(args.infile)}"
ds_out.attrs = {
    'title':         title,
    'Conventions':   'None',
    'source_file':   args.infile,
    'creation_date': datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
}

# ---------------------------------------------------------------------------
# Write output
# ---------------------------------------------------------------------------
os.makedirs(os.path.dirname(os.path.abspath(args.outfile)), exist_ok=True)
if os.path.exists(args.outfile):
    os.remove(args.outfile)
ds_out.to_netcdf(args.outfile)
logger.info(f"Domain file written: {args.outfile}")
