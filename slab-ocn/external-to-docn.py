"""
Take an external file with SST and ice (vars sst and icec, lat N->S) and convert it
to monthly means at the midpoint of each month in the format used for CIME DOCN streams
(lat S->N, with date/datesec). Port of external-to-docn.ncl.
"""
import numpy as np
import xarray as xr
import cftime
import logging
import argparse
import sys
import os

module_path = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', 'py_functions'))
if module_path not in sys.path:
    sys.path.append(module_path)
import pyfuncs


def parse_args_docn():
    """Parse command line arguments for external-to-DOCN conversion."""
    parser = argparse.ArgumentParser(
        description='Convert an external SST/ice file to CIME DOCN format'
    )
    parser.add_argument('--infile', type=str, required=True,
                        help='Input file with sst and icec variables')
    parser.add_argument('--output_file', type=str, default=None,
                        help='Output file path (default: <infile>_CIMEified.nc)')
    parser.add_argument('--verbose', '-v', action='store_true',
                        help='Enable verbose logging')
    return parser.parse_args()


# Parse command line arguments
args = parse_args_docn()

# Configure logging
pyfuncs.configure_logging(args.verbose)

BETACAST, PATHTOHERE = pyfuncs.get_betacast_path()

infile = args.infile
SST_write_file = args.output_file if args.output_file is not None else infile + "_CIMEified.nc"

logging.info(f"Input file: {infile}")
logging.info(f"Output file: {SST_write_file}")

# Open input (decode_times=False to preserve raw time values and attributes)
f = xr.open_dataset(infile, decode_times=False)

# Flip lat to S->N
var_in = f['sst'][:, ::-1, :].load()
ice_in = f['icec'][:, ::-1, :].load()
time_in = f['time']

# get rid of missing data (usually over land areas)
# first go east/west, then north/south for SST and set to 0 for ICEFRAC
var_in.values = pyfuncs.linmsg_n(var_in.values, -1, 2)
var_in.values = pyfuncs.linmsg_n(var_in.values, -1, 1)
ice_in.values = pyfuncs.linmsg_n(ice_in.values, -1, 2)
ice_in.values = np.where(np.isnan(ice_in.values), 0.0, ice_in.values).astype(ice_in.dtype)

units = time_in.attrs.get('units')
calendar = time_in.attrs.get('calendar', 'standard')
try:
    alldates = cftime.num2date(time_in.values, units, calendar=calendar)
except Exception as e:
    logging.error(f"Could not decode time in {infile} (needs CF units/calendar): {e!r}")
    sys.exit(1)

# midpoint of each month in the file's calendar (handles leap years, noleap, 360_day)
middates = []
for d in alldates:
    start = cftime.datetime(d.year, d.month, 1, calendar=calendar)
    end = cftime.datetime(d.year + d.month // 12, d.month % 12 + 1, 1, calendar=calendar)
    middates.append(start + (end - start) / 2)

# get YYYYMMDD and SSSSS, and put time at the midpoint so all three agree
date = np.array([m.year * 10000 + m.month * 100 + m.day for m in middates], dtype=np.int32)
datesec = np.array([m.hour * 3600 + m.minute * 60 + m.second for m in middates], dtype=np.int32)
time_mid = np.asarray(cftime.date2num(middates, units, calendar=calendar), dtype=np.float64)

for d, s in zip(date, datesec):
    logging.info(f"{d} {s}")

#----- write output file with correct format

if os.path.exists(SST_write_file):
    logging.info(f"Removing existing output file: {SST_write_file}")
    os.remove(SST_write_file)

out = xr.Dataset(
    {
        'date': (['time'], date),
        'datesec': (['time'], datesec),
        'SST': var_in,
        'ICEFRAC': ice_in,
    },
)
out = out.assign_coords(time=('time', time_mid, time_in.attrs))
out.attrs = pyfuncs.get_provenance_atts(BETACAST)

out.to_netcdf(SST_write_file, unlimited_dims=['time'],
              encoding={c: {'_FillValue': None} for c in ['time', 'lat', 'lon']})
f.close()

logging.info("all done!")
