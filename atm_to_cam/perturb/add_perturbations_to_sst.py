import os
import sys
import numpy as np
import xarray as xr
import argparse
from datetime import datetime
# Betacast modules
module_paths = [
    ('functions_path', ['../..', 'py_functions']),
]
for path_name, path_parts in module_paths:
    path = os.path.abspath(os.path.join(os.path.dirname(__file__), *path_parts))
    if path not in sys.path:
        sys.path.append(path)
from py_seedfuncs import keyword_values
import horizremap

# Set up argument parser
parser = argparse.ArgumentParser(description='Add SST perturbations for pseudo-global warming')
parser.add_argument('--BEFOREPERTFILE', type=str, required=True,
                    help='Input SST file before perturbation')
parser.add_argument('--AFTERPERTFILE', type=str, required=True,
                    help='Output SST file after perturbation')
parser.add_argument('--pthi', type=str, required=True,
                    help='Path to perturbation namelist file')

args = parser.parse_args()

# Get configuration values from namelist
pthi = args.pthi
warming_case = keyword_values(pthi, "case", "str")
basedir = keyword_values(pthi, "basedir", "str")
adjust_ice = keyword_values(pthi, "adjust_ice", "bool")
output_sst_diag = keyword_values(pthi, "output_sst_diag", "bool")
start_month = keyword_values(pthi, "start_month", "int")
end_month = keyword_values(pthi, "end_month", "int")
current_year = keyword_values(pthi, "current_year", "int")
comp_year = keyword_values(pthi, "comp_year", "int")

# Print configurations
print("************* Running SST perturbation code *************")
print(f"Case: {warming_case}")
print(f"basedir: {basedir}")
print(f"adjust_ice: {adjust_ice}")
print(f"output_sst_diag: {output_sst_diag}")
print(f"start_month: {start_month}")
print(f"end_month: {end_month}")
print(f"current_year: {current_year}")
print(f"comp_year: {comp_year}")
print(f"BEFOREPERTFILE: {args.BEFOREPERTFILE}")
print(f"AFTERPERTFILE: {args.AFTERPERTFILE}")
print("****************************************************")

# Read from input file; .values.copy() loads data into memory
cam_ds = xr.open_dataset(args.BEFOREPERTFILE)
lat = cam_ds['lat'].values
lon = cam_ds['lon'].values
SST = cam_ds['SST_cpl'].values.copy()
ntimes = SST.shape[0]

# Load delta SST based on warming case
if warming_case == "CAMC20C":
    delta_file = xr.open_dataset(
        f"{basedir}/{warming_case}_plev/delta_ts_CAM5-1-1degree_All-Hist_est1_v2.0_1996-2016.nc_Climatology_2016-1996.nc")
    deltaSST_in = delta_file['delta_ts_Climatology_Monthly'].values[start_month-1:end_month, 0, :, :]
    deltaSST = np.mean(deltaSST_in, axis=0)
    delta_lat = delta_file['lat'].values
    delta_lon = delta_file['lon'].values
    delta_file.close()

elif warming_case == "CESMLENS":
    delta_file = xr.open_dataset(f"{basedir}/{warming_case}_plev/ens_SST_anom.nc")
    current_idx_start = current_year * 12 - 1920 * 12 + start_month - 1
    current_idx_end = current_year * 12 - 1920 * 12 + end_month
    deltaSST_in = delta_file['TS'].values[current_idx_start:current_idx_end, :, :]
    deltaSST_current = np.mean(deltaSST_in, axis=0)

    if comp_year < 1920:
        deltaSST_comp = np.zeros_like(deltaSST_current)
    else:
        comp_idx_start = comp_year * 12 - 1920 * 12 + start_month - 1
        comp_idx_end = comp_year * 12 - 1920 * 12 + end_month
        deltaSST_in = delta_file['TS'].values[comp_idx_start:comp_idx_end, :, :]
        deltaSST_comp = np.mean(deltaSST_in, axis=0)

    deltaSST = deltaSST_comp - deltaSST_current
    delta_lat = delta_file['lat'].values
    delta_lon = delta_file['lon'].values
    delta_file.close()

else:
    print(f"Unknown warming case: {warming_case}")
    sys.exit(1)

# Fill missing/NaN values with 0
deltaSST = np.nan_to_num(deltaSST, nan=0.0)

# Bilinear interpolation from delta grid to SST grid (replaces NCL linint2_Wrap)
deltaSST_interp = horizremap.linint2(delta_lon, delta_lat, deltaSST, True, lon, lat)

print(f"max SST delta: {np.nanmax(deltaSST_interp)}    min SST delta: {np.nanmin(deltaSST_interp)}")

# Add delta SST to all timesteps
for ii in range(ntimes):
    SST[ii, :, :] = SST[ii, :, :] + deltaSST_interp.astype(SST.dtype)

print(f"Writing updated SSTs to: {args.AFTERPERTFILE}")
cam_ds['SST_cpl'].values[...] = SST

# Ice adjustment
if adjust_ice:
    ice_thresh = -1.0  # degC, assuming ~35 per mille salinity
    ice = cam_ds['ice_cov'].values.copy()
    ice0 = ice.copy()
    ice = np.where((ice > 0.0) & (SST > ice_thresh), 0.0, ice)
    deltaice_interp = ice - ice0
    print(f"Writing updated ice to: {args.AFTERPERTFILE}")
    cam_ds['ice_cov'].values[...] = ice

# Write modified dataset
cam_ds.to_netcdf(args.AFTERPERTFILE)

# Diagnostics
if output_sst_diag:
    diag_filename = "deltas_sst.nc"
    print(f"outputting diags to {diag_filename}")

    diag_ds = xr.Dataset(
        data_vars={
            'deltaSST_interp': (['lat', 'lon'], deltaSST_interp),
        },
        coords={
            'lat': lat,
            'lon': lon,
        },
        attrs={
            'creation_date': datetime.now().strftime('%c'),
        }
    )

    if adjust_ice:
        diag_ds['deltaice_interp'] = (['lat', 'lon'], deltaice_interp[0, :, :])

    diag_ds.to_netcdf(diag_filename)

print("SST perturbation complete.")
