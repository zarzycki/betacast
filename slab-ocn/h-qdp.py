import numpy as np
import xarray as xr
import logging
import argparse
import sys
import os

module_path = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', 'py_functions'))
if module_path not in sys.path:
    sys.path.append(module_path)
import pyfuncs
import horizremap


def parse_args_hblt():
    """Parse command line arguments for hblt/qdp addition."""
    parser = argparse.ArgumentParser(
        description='Add hblt and qdp variables to SST file'
    )
    parser.add_argument('--sst_file', type=str, required=True,
                        help='Input SST file (output from sst_to_cam.py)')
    parser.add_argument('--output_file', type=str, required=True,
                        help='Output file path')
    parser.add_argument('--which_h', type=str, default='GODAS',
                        choices=['Z16', 'GODAS'],
                        help='Source for boundary layer heights (default: GODAS)')
    parser.add_argument('--which_q', type=str, default='Z16',
                        choices=['Z16', 'NONE'],
                        help='Source for Q flux correction (default: Z16)')
    parser.add_argument('--copy_ice', action='store_true', default=True,
                        help='Copy ice_cov if available (default: True)')
    parser.add_argument('--no_copy_ice', action='store_false', dest='copy_ice',
                        help='Do not copy ice_cov')
    parser.add_argument('--nc_files_dir', type=str, default=None,
                        help='Directory containing hblt/qdp source files')
    parser.add_argument('--verbose', '-v', action='store_true',
                        help='Enable verbose logging')
    return parser.parse_args()


# Parse command line arguments
args = parse_args_hblt()

# Configure logging
pyfuncs.configure_logging(args.verbose)

# Get BETACAST path
BETACAST, PATHTOHERE = pyfuncs.get_betacast_path()

# Unpack arguments
sfile = args.sst_file
output_file = args.output_file
WHICH_H = args.which_h
WHICH_Q = args.which_q
copy_ice_cov_if_avail = args.copy_ice
nc_files_dir = args.nc_files_dir

# Set default nc_files_dir based on BETACAST if not provided
if args.nc_files_dir is None:
    nc_files_dir = os.path.join(BETACAST, 'slab-ocn', 'nc_files')
else:
    nc_files_dir = args.nc_files_dir

logging.info(f"Input SST file: {sfile}")
logging.info(f"Output file: {output_file}")
logging.info(f"nc_files_dir: {nc_files_dir}")
logging.info(f"WHICH_H: {WHICH_H}")
logging.info(f"WHICH_Q: {WHICH_Q}")
logging.info(f"copy_ice_cov_if_avail: {copy_ice_cov_if_avail}")

# Open SST file (decode_times=False to preserve raw time values)
logging.info(f"Opening SST file: {sfile}")
s = xr.open_dataset(sfile, decode_times=False)
SST_cpl = s['SST_cpl'].values
SST_cpl_dims = SST_cpl.shape
ntim = SST_cpl_dims[0]
nlat = SST_cpl_dims[1]
nlon = SST_cpl_dims[2]
lat = s['lat'].values
lon = s['lon'].values

logging.info(f"SST dimensions: ntim={ntim}, nlat={nlat}, nlon={nlon}")

# Check if we should copy ice_cov
copy_ice = copy_ice_cov_if_avail and 'ice_cov' in s.variables
if copy_ice:
    logging.info("Will copy ice_cov from input file")
    ice_cov = s['ice_cov'].values
else:
    logging.info("ice_cov will not be copied")

# Get boundary layer heights
if WHICH_H == "Z16":
    afile = os.path.join(nc_files_dir, "hq_z16.nc")
    logging.info(f"Loading hblt from Z16: {afile}")
    a = xr.open_dataset(afile)
    hblt = a['hblt'].values
    a.close()

elif WHICH_H == "GODAS":
    afile = os.path.join(nc_files_dir, "dbss_obml.mon.ltm.1991-2020.nc")
    logging.info(f"Loading hblt from GODAS: {afile}")
    a = xr.open_dataset(afile, decode_times=False)
    hblt_o = a['dbss_obml'].values
    hblt_o_lon = a['lon'].values
    hblt_o_lat = a['lat'].values
    a.close()

    # Interpolate to SST grid (12 months)
    logging.info("Interpolating GODAS hblt to SST grid")
    nmonths = hblt_o.shape[0]
    hblt = np.empty((nmonths, nlat, nlon), dtype=np.float64)
    for m in range(nmonths):
        hblt[m, :, :] = horizremap.linint2(
            hblt_o_lon, hblt_o_lat, hblt_o[m, :, :], True, lon, lat
        )

else:
    logging.error(f"Unsupported WHICH_H: {WHICH_H}")
    sys.exit(1)

# Get Q flux correction
if WHICH_Q == "Z16":
    qfile = os.path.join(nc_files_dir, "hq_z16.nc")
    logging.info(f"Loading qdp from Z16: {qfile}")
    q = xr.open_dataset(qfile)
    qdp = q['qdp'].values
    q.close()

elif WHICH_Q == "NONE":
    logging.info("Setting qdp to zero")
    qdp = np.zeros_like(hblt)

else:
    logging.error(f"Unsupported WHICH_Q: {WHICH_Q}")
    sys.exit(1)

# Expand hblt and qdp to full time dimension by tiling monthly climatology
nyears = ntim // 12
logging.info(f"Tiling {nyears} years of monthly climatology")

hblt_full = np.empty((ntim, nlat, nlon), dtype=np.float64)
qdp_full = np.empty((ntim, nlat, nlon), dtype=np.float64)

for ii in range(nyears):
    ix = ii * 12
    hblt_full[ix:ix+12, :, :] = hblt
    qdp_full[ix:ix+12, :, :] = qdp

# Fill missing values along longitude (dim 2) then latitude (dim 1)
logging.info("Filling missing values in hblt")
hblt_full = pyfuncs.linmsg_n(hblt_full, -1, 2)
hblt_full = pyfuncs.linmsg_n(hblt_full, -1, 1)
#hblt_full = np.where(np.isnan(hblt_full), 40.0, hblt_full)

# Clamp hblt minimum to 5.0 and maximum to 2500.0
hblt_full = np.where(hblt_full <    5.0,    5.0, hblt_full)
hblt_full = np.where(hblt_full > 2500.0, 2500.0, hblt_full)

# Get time coordinate from input (strip attrs to avoid encoding conflicts)
times = s['time'].values
time_attrs = dict(s['time'].attrs)  # Capture original time attributes

# Close input file
s.close()

# Remove existing output file if present
if os.path.exists(output_file):
    logging.info(f"Removing existing output file: {output_file}")
    os.remove(output_file)

# Write output
logging.info(f"Writing output file: {output_file}")

# NOTE: As of 11/25, it seems if hblt and qdp are not cast as doubles, problems arise
out_vars = {
    'SST_cpl': (['time', 'lat', 'lon'], SST_cpl.astype(np.float64)),
    'hblt': (['time', 'lat', 'lon'], hblt_full.astype(np.float64)),
    'qdp': (['time', 'lat', 'lon'], qdp_full.astype(np.float64)),
}

if copy_ice:
    out_vars['ice_cov'] = (['time', 'lat', 'lon'], ice_cov.astype(np.float64))

out_ds = xr.Dataset(
    out_vars,
    coords={
        'time': times,
        'lat': lat,
        'lon': lon
    }
)

# Set time attributes (use input attrs, fallback to defaults)
out_ds['time'].attrs['units'] = time_attrs.get('units', 'days since 0001-01-01 00:00:00')
out_ds['time'].attrs['calendar'] = time_attrs.get('calendar', '365_day')

# Add global attributes
out_ds.attrs['title'] = "SST file with hblt and qdp added"
out_ds.attrs['source_file'] = sfile
out_ds.attrs['hblt_source'] = WHICH_H
out_ds.attrs['qdp_source'] = WHICH_Q
out_ds.attrs['Conventions'] = "None"

out_ds.to_netcdf(output_file, unlimited_dims=['time'])
logging.info(f"Data written to {output_file}")