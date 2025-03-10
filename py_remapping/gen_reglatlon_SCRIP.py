import xarray as xr
import numpy as np
import os
import argparse

from ESMF_regridding import *

# Parse command line arguments
parser = argparse.ArgumentParser(description="Generate SCRIP grid file from input data.")
parser.add_argument('--srcfilename', type=str, required=True, help='Source filename for input data')
parser.add_argument('--dstDir', type=str, default="./", help='Directory to save the output SCRIP file')
parser.add_argument('--dstGridName', type=str, default="hwrf_storm_scrip.nc", help='Output SCRIP filename')
parser.add_argument('--metatitle', type=str, default="HWRF", help='Metadata title for the SCRIP file')
parser.add_argument('--appendLL', action='store_true', help='Append top left corner lat/lon to the filename')
args = parser.parse_args()

# Load data
src_file = xr.open_dataset(args.srcfilename, engine='cfgrib', filter_by_keys={'stepType': 'instant', 'typeOfLevel': 'surface'})

# Load lat/lon
lat = src_file.variables['latitude'].values
lon = src_file.variables['longitude'].values
print(lat)
print(lon)

# Figure out top lat/lon corners
top_left_lat = lat[0]
top_left_lon = lon[0]

if args.appendLL:
    args.dstGridName = f"{args.dstGridName}_{top_left_lat}_{top_left_lon}.nc"

print(f"Creating SCRIP file from: {args.srcfilename}")
print(f"Writing SCRIP file to: {os.path.join(args.dstDir, args.dstGridName)}")
print(f"Metaname for SCRIP file is: {args.metatitle}")

# Set options for SCRIP file generation
Opt = {
    'ForceOverwrite': True,
    'PrintTimings': True,
    'Title': args.metatitle
}

# Generate SCRIP file
rectilinear_to_SCRIP(os.path.join(args.dstDir, args.dstGridName), lat, lon, Opt)

print("done")
