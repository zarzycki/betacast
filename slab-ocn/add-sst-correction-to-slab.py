"""
Replace qdp on an existing slab ocean DOCN file with monthly SST bias corrections
(diffMM.nc, e.g., SLAB - FIXEDSST) bilinearly interpolated to the DOCN grid. Each diffMM.nc
is applied to every time in calendar month MM (from the decoded time coordinate). Months
without a correction file are set to zero. The slab file is modified in place.
Port of add-sst-correction-to-slab.ncl.
"""
import numpy as np
import netCDF4 as nc
import cftime
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
                        help='Slab DOCN file (with qdp and lat/lon or yc/xc) to modify in place')
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

# Target grid: lat/lon (h-qdp.py output) or yc/xc (older/domain-style slab files)
latname, lonname = ('lat', 'lon') if 'lat' in inf.variables else ('yc', 'xc')
lat = np.ma.filled(inf.variables[latname][:], np.nan)
lon = np.ma.filled(inf.variables[lonname][:], np.nan)
if lat.ndim == 2:
    # 2D (nj, ni) coordinates are only usable if the grid is regular lat-lon
    if not (np.all(lat == lat[:, :1]) and np.all(lon == lon[:1, :])):
        logging.error(f"{latname}/{lonname} in {SST_write_file} are 2D and not a regular lat-lon grid")
        sys.exit(1)
    lat, lon = lat[:, 0], lon[0, :]
logging.info(f"Target grid from {latname}/{lonname}: nlat={lat.size}, nlon={lon.size}")

logging.info(f"qdp: {corr.dtype} {corr.shape}")

# Calendar month of each time in the slab file
ntim = corr.shape[0]
try:
    tvar = inf.variables['time']
    dates = cftime.num2date(tvar[:], tvar.units, calendar=getattr(tvar, 'calendar', 'standard'))
    months = np.array([d.month for d in dates])
    time_desc = f"{dates[0]} to {dates[-1]}"
except Exception as e:
    # Older slab files may lack CF time; whole years are assumed to be repeating Jan-Dec
    if ntim % 12 != 0:
        logging.error(f"Could not decode time in {SST_write_file} (needs CF units/calendar) "
                      f"and it has {ntim} times, not a multiple of 12: {e!r}")
        sys.exit(1)
    logging.warning(f"Could not decode time in {SST_write_file} ({e!r}); file has {ntim} times, "
                    f"ASSUMING {ntim // 12} year(s) of months Jan-Dec in order")
    months = np.arange(ntim) % 12 + 1
    time_desc = f"assumed {ntim // 12} x Jan-Dec"

logging.info("=" * 65)
logging.info(f"Slab file has {ntim} times: {time_desc}")
logging.info("Each diffMM.nc is applied to every time in calendar month MM:")
if ntim != 12:
    logging.warning(f"Slab file is not a single 12-month climatology ({ntim} times), "
                    "so each monthly correction is repeated for every year")

for i in range(12):
    monthstr = f"{i+1:02d}"
    corrFile = os.path.join(SST_bias_dir, f"diff{monthstr}.nc")
    imon = months == i + 1
    if not np.any(imon):
        logging.info(f"  month {monthstr}: not in slab file, skipping")
        continue
    if os.path.isfile(corrFile):
        with nc.Dataset(corrFile, "r") as f:
            corrLat = np.ma.filled(f.variables['lat'][:], np.nan)
            corrLon = np.ma.filled(f.variables['lon'][:], np.nan)
            corrVar = np.ma.filled(f.variables['SST'][0, :, :], np.nan)
        corrOut = horizremap.linint2(corrLon, corrLat, corrVar, True, lon, lat)
        corr[imon, :, :] = corrOut
        logging.info(f"  month {monthstr}: diff{monthstr}.nc -> {imon.sum()} times "
                     f"(min: {np.nanmin(corrOut):.4g}, max: {np.nanmax(corrOut):.4g})")
    else:
        corr[imon, :, :] = 0.0
        logging.warning(f"  month {monthstr}: no diff{monthstr}.nc, qdp set to 0.0 for {imon.sum()} times")
logging.info("=" * 65)

# No correction where the bias is missing (e.g., land or outside bias grid)
if np.any(np.isnan(corr)):
    logging.warning(f"Setting {np.isnan(corr).sum()} missing qdp values to 0.0")
    corr = np.where(np.isnan(corr), 0.0, corr)

# Write to NetCDF
logging.info("Writing file...")
inf.variables['qdp'][:] = corr
inf.close()
