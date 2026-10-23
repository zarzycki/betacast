"""
Replace qdp on an existing slab ocean DOCN file with monthly SST bias corrections
(diffMM.nc, e.g., SLAB - FIXEDSST) bilinearly interpolated to the DOCN grid. Months
without a correction file are set to zero. The slab file is modified in place.
Port of add-sst-correction-to-slab.ncl.
"""
import numpy as np
import netCDF4 as nc
import logging
import argparse
import sys
import os

module_path = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', 'py_functions'))
if module_path not in sys.path:
    sys.path.append(module_path)
import pyfuncs
import horizremap


def parse_args_corr():
    """Parse command line arguments for adding SST correction to slab file."""
    parser = argparse.ArgumentParser(
        description='Add monthly SST bias correction to slab DOCN file as qdp'
    )
    parser.add_argument('--slab_file', type=str, required=True,
                        help='Slab DOCN file (with qdp, xc, yc) to modify in place')
    parser.add_argument('--bias_dir', type=str, required=True,
                        help='Directory containing diffMM.nc files with SST bias')
    parser.add_argument('--verbose', '-v', action='store_true',
                        help='Enable verbose logging')
    return parser.parse_args()


# Parse command line arguments
args = parse_args_corr()

# Configure logging
pyfuncs.configure_logging(args.verbose)

SST_write_file = args.slab_file
SST_bias_dir = args.bias_dir

logging.info(f"Slab file: {SST_write_file}")
logging.info(f"Bias dir: {SST_bias_dir}")

inf = nc.Dataset(SST_write_file, "r+")

corr = np.ma.filled(inf.variables['qdp'][:], np.nan)
lat = np.ma.filled(inf.variables['yc'][:], np.nan)
lon = np.ma.filled(inf.variables['xc'][:], np.nan)

logging.info(f"qdp: {corr.dtype} {corr.shape}")

for i in range(12):
    monthstr = f"{i+1:02d}"
    corrFile = os.path.join(SST_bias_dir, f"diff{monthstr}.nc")
    if os.path.isfile(corrFile):
        logging.info(f"Found correction from month: {monthstr}")
        with nc.Dataset(corrFile, "r") as f:
            corrLat = np.ma.filled(f.variables['lat'][:], np.nan)
            corrLon = np.ma.filled(f.variables['lon'][:], np.nan)
            corrVar = np.ma.filled(f.variables['SST'][0, :, :], np.nan)
        corr[i, :, :] = horizremap.linint2(corrLon, corrLat, corrVar, True, lon, lat)
        logging.info(f"corr{monthstr} min: {np.nanmin(corr[i])} max: {np.nanmax(corr[i])}")
    else:
        corr[i, :, :] = 0.0

# No correction where the bias is missing (e.g., land or outside bias grid)
if np.any(np.isnan(corr)):
    logging.warning(f"Setting {np.isnan(corr).sum()} missing qdp values to 0.0")
    corr = np.where(np.isnan(corr), 0.0, corr)

# Write to NetCDF
logging.info("Writing file...")
inf.variables['qdp'][:] = corr
inf.close()
