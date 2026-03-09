#!/usr/bin/env python3
"""
Generate DATM forcing files from CR20V3 (NOAA 20th Century Reanalysis V3).

Source data: annual per-variable netCDF files in the CR20V3 anl/ directory,
3-hourly at ~0.7-degree Gaussian resolution (256x512).

Output: monthly per-stream netCDF files in TPQW/, Solar/, Precip/ subdirs.

Strategy: all selected timesteps for the year are loaded into memory in a
single pass (8 disk reads total), then sliced by month in numpy.

Usage:
    python gen-forcing-cr20v3.py --year=1901 --rawdatadir=/path/to/anl --outdirbase=/out/ [options]
"""

import numpy as np
import xarray as xr
import os
import sys
import argparse
import logging
import cftime
from datetime import datetime
import subprocess
import shutil
import time as time_mod

import dask
dask.config.set(scheduler='threads', num_workers=4)

logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    datefmt='%Y-%m-%d %H:%M:%S',
)
logger = logging.getLogger(__name__)


def log_stats(name, arr, units=""):
    """Log min/max/mean and NaN count for a numpy array."""
    nnan = int(np.sum(np.isnan(arr)))
    logger.info(f"    {name:8s}: min={arr.min():12.4g}  max={arr.max():12.4g}  "
                f"mean={arr.mean():12.4g}  {units}  NaN={nnan}")


def log_qc(name, arr_before, arr_after):
    """Log how many values were corrected by QC."""
    n = int(np.sum(arr_before != arr_after))
    if n > 0:
        logger.info(f"    QC {name:8s}: {n} values corrected")
    else:
        logger.info(f"    QC {name:8s}: no corrections needed")



# ---------------------------------------------------------------------------
# Argument parsing
# ---------------------------------------------------------------------------
parser = argparse.ArgumentParser(description="Generate DATM forcing from CR20V3")
parser.add_argument("--year",           required=True, type=int,
                    help="Year to process (YYYY); all 12 months are written")
parser.add_argument("--rawdatadir",     required=True,
                    help="Path to CR20V3 anl/ directory")
parser.add_argument("--outdirbase",     required=True,
                    help="Base output directory (TPQW/, Solar/, Precip/ created beneath)")
parser.add_argument("--datafilename",   default="CMZCR20V3.v0.c2025.0.7d",
                    help="Prefix for output filenames")
parser.add_argument("--time_stride",    default=24, type=int,
                    help="Output time stride in hours (3=native, 6=6-hourly, 24=daily 00Z). "
                         "Must be a multiple of 3 (source resolution). Default: 24")
parser.add_argument("--greg_to_noleap", action="store_true",
                    help="Convert proleptic_gregorian calendar to noleap (drop Feb 29)")
parser.add_argument("--convert_nc3",    action="store_true",
                    help="Convert output files to NetCDF3 classic via ncks")

args           = parser.parse_args()
YYYY           = args.year
rawdatadir     = args.rawdatadir
outdirbase     = args.outdirbase
datafilename   = args.datafilename
time_stride    = args.time_stride
greg_to_noleap = args.greg_to_noleap
convert_nc3    = args.convert_nc3

if time_stride % 3 != 0 or not (3 <= time_stride <= 24):
    logger.error(f"--time_stride must be 3, 6, 12, or 24 (got {time_stride})")
    sys.exit(1)

# ---------------------------------------------------------------------------
# Startup banner
# ---------------------------------------------------------------------------
logger.info("=" * 70)
logger.info("  gen-forcing-cr20v3.py")
logger.info(f"  year          : {YYYY}")
logger.info(f"  rawdatadir    : {rawdatadir}")
logger.info(f"  outdirbase    : {outdirbase}")
logger.info(f"  datafilename  : {datafilename}")
logger.info(f"  time_stride   : {time_stride}h")
logger.info(f"  greg_to_noleap: {greg_to_noleap}")
logger.info(f"  convert_nc3   : {convert_nc3}")
logger.info("=" * 70)

# ---------------------------------------------------------------------------
# Variable file mapping:  logical key -> (filename, nc variable name)
# All variables live in anl_mean_{YEAR}_{VARCODE}_{LEVEL}.nc
# ---------------------------------------------------------------------------
VAR_FILES = {
    'tbot': (f'anl_mean_{YYYY}_TMP_2m.nc',    't2m'),    # K
    'u10':  (f'anl_mean_{YYYY}_UGRD_10m.nc',  'u10'),    # m/s
    'v10':  (f'anl_mean_{YYYY}_VGRD_10m.nc',  'v10'),    # m/s
    'qbot': (f'anl_mean_{YYYY}_SPFH_2m.nc',   'q'),      # kg/kg
    'psrf': (f'anl_mean_{YYYY}_PRES_sfc.nc',  'sp'),     # Pa
    'flds': (f'anl_mean_{YYYY}_DLWRF_sfc.nc', 'dlwrf'),  # W/m^2 (instantaneous, no /3600)
    'fsds': (f'anl_mean_{YYYY}_DSWRF_sfc.nc', 'dswrf'),  # W/m^2 (instantaneous, no /3600)
    'prec': (f'anl_mean_{YYYY}_PRATE_sfc.nc', 'prate'),  # kg/m^2/s == mm/s (no conversion)
}

# Leap-year flag (proleptic_gregorian: div-by-4 is always leap)
is_leap = (YYYY % 4 == 0 and YYYY % 100 != 0) or (YYYY % 400 == 0)
logger.info(f"Leap year check: {YYYY} -> is_leap={is_leap}")

# ---------------------------------------------------------------------------
# Open all input datasets (lazy via dask; decode_times=False avoids cftime
# overflow when calling cftime.num2date on already-decoded values)
# ---------------------------------------------------------------------------
logger.info("-" * 70)
logger.info("Opening input files...")
datasets = {}
for key, (fname, vname) in VAR_FILES.items():
    fpath = os.path.join(rawdatadir, fname)
    if not os.path.exists(fpath):
        logger.error(f"Missing required file: {fpath}")
        sys.exit(1)
    fsize_mb = os.path.getsize(fpath) / 1e6
    datasets[key] = xr.open_dataset(fpath, chunks={'time': 'auto'}, decode_times=False)
    shape = datasets[key][vname].shape
    logger.info(f"  {key:6s}: {fname}  ({fsize_mb:.0f} MB)  shape={shape}")

# ---------------------------------------------------------------------------
# Grid info (identical for all files)
# ---------------------------------------------------------------------------
ds0     = datasets['tbot']
lat_raw = ds0['latitude'].values   # (256,)  89.46 -> -89.46  (N->S)
lon_raw = ds0['longitude'].values  # (512,)  0     -> 359.3

# Flip lat to S->N (standard DATM convention)
lat = lat_raw[::-1]
lon = lon_raw
nlat, nlon = len(lat), len(lon)

logger.info(f"Native grid: {nlat} lat x {nlon} lon")
logger.info(f"  lat: {lat[0]:.4f} -> {lat[-1]:.4f} (S->N after flip)")
logger.info(f"  lon: {lon[0]:.4f} -> {lon[-1]:.4f}")

lon_datm_2D = np.tile(lon,          (nlat, 1)).astype(np.float32)
lat_datm_2D = np.tile(lat[:, None], (1, nlon)).astype(np.float32)

# ---------------------------------------------------------------------------
# Parse full-year time coordinate once
# ---------------------------------------------------------------------------
logger.info("-" * 70)
logger.info("Parsing time coordinate...")
time_da       = ds0['time']
time_units    = time_da.attrs.get('units',    'hours since 1800-01-01 00:00:00.0')
time_calendar = time_da.attrs.get('calendar', 'proleptic_gregorian')
time_values   = time_da.values.astype(np.float64)
time_cf       = cftime.num2date(time_values, time_units, time_calendar)

time_month = np.array([t.month for t in time_cf], dtype=np.int32)
time_day   = np.array([t.day   for t in time_cf], dtype=np.int32)
time_hour  = np.array([t.hour  for t in time_cf], dtype=np.int32)

logger.info(f"  Calendar    : {time_calendar}")
logger.info(f"  Time units  : {time_units}")
logger.info(f"  Total steps : {len(time_cf)}  ({len(time_cf)/8:.0f} days * 8 steps/day)")
logger.info(f"  First step  : {time_cf[0]}")
logger.info(f"  Last step   : {time_cf[-1]}")

# ---------------------------------------------------------------------------
# Shared output metadata
# ---------------------------------------------------------------------------
new_time_units    = "days since 1800-01-01 00:00:00"
time_calendar_out = "noleap" if greg_to_noleap else time_calendar

EDGEW = np.array([0.],   dtype=np.float32)
EDGEE = np.array([360.], dtype=np.float32)
EDGES = np.array([-90.], dtype=np.float32)
EDGEN = np.array([90.],  dtype=np.float32)

FILL_FLOAT = np.float32(9.96921e+36)
FILL_OTHER = np.int32(-9999)

# ---------------------------------------------------------------------------
# Helper: write dataset to netCDF, log file size, optionally convert to NC3
# ---------------------------------------------------------------------------
def write_nc(ds, filepath):
    os.makedirs(os.path.dirname(filepath), exist_ok=True)
    if os.path.exists(filepath):
        os.remove(filepath)
    encoding = {}
    for var in ds.variables:
        fv = float(FILL_FLOAT) if ds[var].dtype in (np.float32, np.float64) else int(FILL_OTHER)
        encoding[var] = {'zlib': True, 'complevel': 1, '_FillValue': fv}
    t0 = time_mod.time()
    ds.to_netcdf(filepath, encoding=encoding, format='NETCDF4', unlimited_dims=['time'])
    fsize_mb = os.path.getsize(filepath) / 1e6
    logger.info(f"    Wrote: {os.path.basename(filepath)}  ({fsize_mb:.1f} MB, {time_mod.time()-t0:.1f}s)")
    if convert_nc3:
        if shutil.which('ncks'):
            t0 = time_mod.time()
            result = subprocess.run(['ncks', '-O', '-3', filepath, filepath],
                                    capture_output=True, text=True)
            if result.returncode == 0:
                fsize_mb = os.path.getsize(filepath) / 1e6
                logger.info(f"    NC3 convert OK ({fsize_mb:.1f} MB, {time_mod.time()-t0:.1f}s)")
            else:
                logger.warning(f"    ncks conversion failed: {result.stderr}")
        else:
            logger.warning("    ncks not found; skipping NC3 conversion")


# ---------------------------------------------------------------------------
# Load all selected timesteps for the full year in one pass (8 disk reads),
# then slice by month in numpy (near-instant).
#
# time_stride controls which hours to keep:
#   24h -> 00Z only  (every 8th source step)
#   6h  -> 00/06/12/18Z  (every 2nd source step)
#   3h  -> all steps (native resolution)
# ---------------------------------------------------------------------------
script_t0 = time_mod.time()

annual_sel = (time_hour % time_stride) == 0
if greg_to_noleap and is_leap:
    n_before   = int(annual_sel.sum())
    annual_sel = annual_sel & ~((time_month == 2) & (time_day == 29))
    logger.info(f"Dropped Feb 29 for noleap: {n_before} -> {int(annual_sel.sum())} steps")
annual_indices = np.where(annual_sel)[0].tolist()
annual_cf      = [time_cf[i] for i in annual_indices]
annual_month   = np.array([t.month for t in annual_cf], dtype=np.int32)

steps_per_day = 24 // time_stride
logger.info("-" * 70)
logger.info(f"Loading {len(annual_indices)} steps for {YYYY} "
            f"(stride={time_stride}h, {steps_per_day} steps/day) in one pass...")

def load_var(key):
    vname = VAR_FILES[key][1]
    logger.info(f"  Reading {key} ({vname})...")
    t0  = time_mod.time()
    arr = datasets[key][vname].isel(time=annual_indices).values  # (nsteps,256,512)
    arr = arr[:, ::-1, :].astype(np.float32)                     # flip lat N->S -> S->N
    logger.info(f"    done in {time_mod.time()-t0:.1f}s  shape={arr.shape}  "
                f"min={arr.min():.4g}  max={arr.max():.4g}")
    return arr

tbot_yr = load_var('tbot')
u10_yr  = load_var('u10')
v10_yr  = load_var('v10')
qbot_yr = load_var('qbot')
psrf_yr = load_var('psrf')
flds_yr = load_var('flds')
fsds_yr = load_var('fsds')
prec_yr = load_var('prec')

wind_yr = np.sqrt(u10_yr**2 + v10_yr**2)
del u10_yr, v10_yr
logger.info(f"All variables loaded. Total load time: {time_mod.time()-script_t0:.1f}s")

# ---------------------------------------------------------------------------
# Main loop: slice by month (pure numpy) and write
# ---------------------------------------------------------------------------
for MM in range(1, 13):
    MM_str   = f"{MM:02d}"
    month_t0 = time_mod.time()
    logger.info("=" * 70)
    logger.info(f"  Processing {YYYY}-{MM_str}")
    logger.info("=" * 70)

    mm_mask = annual_month == MM
    if not mm_mask.any():
        logger.warning(f"  No timesteps found for {YYYY}-{MM_str}, skipping")
        continue

    sel_cf = [annual_cf[i] for i in np.where(mm_mask)[0]]
    logger.info(f"  {mm_mask.sum()} timesteps  (numpy slice, no disk I/O)")
    logger.info(f"  First: {sel_cf[0]}   Last: {sel_cf[-1]}")

    # Convert selected times to output calendar
    if greg_to_noleap:
        time_out = np.array([
            cftime.date2num(
                cftime.DatetimeNoLeap(t.year, t.month, t.day, t.hour, t.minute, t.second),
                units=new_time_units, calendar='noleap')
            for t in sel_cf
        ], dtype=np.float64)
    else:
        time_out = np.array([
            cftime.date2num(t, units=new_time_units, calendar=time_calendar)
            for t in sel_cf
        ], dtype=np.float64)

    logger.info(f"  Output time: {time_out[0]:.4f} -> {time_out[-1]:.4f}  ({new_time_units}, {time_calendar_out})")

    time_attrs = {
        'long_name': 'observation time',
        'units':     new_time_units,
        'calendar':  time_calendar_out,
    }

    # Slice this month from pre-loaded annual arrays
    tbot = tbot_yr[mm_mask]
    wind = wind_yr[mm_mask]
    qbot = qbot_yr[mm_mask]
    psrf = psrf_yr[mm_mask]
    flds = flds_yr[mm_mask]
    fsds = fsds_yr[mm_mask]
    prec = prec_yr[mm_mask]
    zbot = np.full_like(tbot, 10.0, dtype=np.float32)

    logger.info("  --- Raw stats (pre-QC) ---")
    log_stats('tbot', tbot, 'K')
    log_stats('wind', wind, 'm/s')
    log_stats('psrf', psrf, 'Pa')
    log_stats('qbot', qbot, 'kg/kg')
    log_stats('flds', flds, 'W/m2')
    log_stats('fsds', fsds, 'W/m2')
    log_stats('prec', prec, 'mm/s')

    # ------------------------------------------------------------------
    # QC / bounds enforcement
    # ------------------------------------------------------------------
    logger.info("  --- Applying QC / bounds ---")

    prec_orig = prec.copy(); prec = np.where(prec < 0.0, 0.0, prec);   log_qc('prec<0',  prec_orig, prec)
    fsds_orig = fsds.copy(); fsds = np.where(fsds <= 0.0, 1e-8, fsds); log_qc('fsds<=0', fsds_orig, fsds)
    flds_orig = flds.copy(); flds = np.where(flds <= 0.0, 1e-8, flds); log_qc('flds<=0', flds_orig, flds)
    qbot_orig = qbot.copy(); qbot = np.clip(qbot, 1e-8, 0.2);          log_qc('qbot',    qbot_orig, qbot)
    wind_orig = wind.copy(); wind = np.clip(wind, 0.0,  45.0);         log_qc('wind',    wind_orig, wind)
    psrf_orig = psrf.copy(); psrf = np.clip(psrf, 30000.0, 120.0e3);   log_qc('psrf',    psrf_orig, psrf)
    tbot_orig = tbot.copy(); tbot = np.clip(tbot, 175.0,   340.0);     log_qc('tbot',    tbot_orig, tbot)
    del prec_orig, fsds_orig, flds_orig, qbot_orig, wind_orig, psrf_orig, tbot_orig

    logger.info("  --- Post-QC stats ---")
    log_stats('tbot', tbot, 'K')
    log_stats('wind', wind, 'm/s')
    log_stats('psrf', psrf, 'Pa')
    log_stats('qbot', qbot, 'kg/kg')
    log_stats('flds', flds, 'W/m2')
    log_stats('fsds', fsds, 'W/m2')
    log_stats('prec', prec, 'mm/s')

    # ------------------------------------------------------------------
    # Build base dataset (coordinates shared by all three streams)
    # ------------------------------------------------------------------
    def base_ds():
        ds = xr.Dataset()
        ds['time']   = xr.DataArray(time_out, dims=['time'], attrs=time_attrs)
        ds['LONGXY'] = xr.DataArray(lon_datm_2D, dims=['lat', 'lon'],
                                    attrs={'long_name': 'longitude', 'units': 'degrees_east',  'mode': 'time-invariant'})
        ds['LATIXY'] = xr.DataArray(lat_datm_2D, dims=['lat', 'lon'],
                                    attrs={'long_name': 'latitude',  'units': 'degrees_north', 'mode': 'time-invariant'})
        ds['EDGEW']  = xr.DataArray(EDGEW, dims=['scalar'], attrs={'long_name': 'western edge in atmospheric data', 'mode': 'time-invariant'})
        ds['EDGEE']  = xr.DataArray(EDGEE, dims=['scalar'], attrs={'long_name': 'eastern edge in atmospheric data', 'mode': 'time-invariant'})
        ds['EDGES']  = xr.DataArray(EDGES, dims=['scalar'], attrs={'long_name': 'southern edge in atmospheric data', 'mode': 'time-invariant'})
        ds['EDGEN']  = xr.DataArray(EDGEN, dims=['scalar'], attrs={'long_name': 'northern edge in atmospheric data', 'mode': 'time-invariant'})
        return ds

    global_attrs = {
        'title':         f"DATM forcing from CR20V3: {YYYY} {MM_str}",
        'source':        rawdatadir,
        'Conventions':   "None",
        'creation_date': datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
        'notes':         f"Generated with Betacast toolkit; {time_stride}h stride from 3-hourly CR20V3",
    }

    # Stream 1: TPQWL
    logger.info("  --- Writing TPQWL ---")
    ds_tpqw = base_ds()
    ds_tpqw['TBOT'] = xr.DataArray(tbot, dims=['time', 'lat', 'lon'],
        attrs={'mode': 'time-dependent', 'units': 'K',           'long_name': 'temperature at the lowest model level'})
    ds_tpqw['WIND'] = xr.DataArray(wind, dims=['time', 'lat', 'lon'],
        attrs={'mode': 'time-dependent', 'units': 'm/s',         'long_name': 'wind magnitude at the lowest model level'})
    ds_tpqw['PSRF'] = xr.DataArray(psrf, dims=['time', 'lat', 'lon'],
        attrs={'mode': 'time-dependent', 'units': 'Pa',          'long_name': 'surface pressure'})
    ds_tpqw['QBOT'] = xr.DataArray(qbot, dims=['time', 'lat', 'lon'],
        attrs={'mode': 'time-dependent', 'units': 'kg/kg',       'long_name': 'specific humidity at lowest model level'})
    ds_tpqw['ZBOT'] = xr.DataArray(zbot, dims=['time', 'lat', 'lon'],
        attrs={'mode': 'time-dependent', 'units': 'm',           'long_name': 'reference height'})
    ds_tpqw['FLDS'] = xr.DataArray(flds, dims=['time', 'lat', 'lon'],
        attrs={'mode': 'time-dependent', 'units': 'W/m**2',      'long_name': 'downward longwave flux at surface'})
    ds_tpqw.attrs = global_attrs
    write_nc(ds_tpqw, os.path.join(outdirbase, 'TPQW',   f"{datafilename}.TPQWL.{YYYY}-{MM_str}.nc"))
    del ds_tpqw, tbot, wind, psrf, qbot, zbot, flds

    # Stream 2: Solar
    logger.info("  --- Writing Solar ---")
    ds_solar = base_ds()
    ds_solar['FSDS'] = xr.DataArray(fsds, dims=['time', 'lat', 'lon'],
        attrs={'mode': 'time-dependent', 'units': 'W/m**2',      'long_name': 'downward solar flux at surface'})
    ds_solar.attrs = global_attrs
    write_nc(ds_solar, os.path.join(outdirbase, 'Solar',  f"{datafilename}.Solar.{YYYY}-{MM_str}.nc"))
    del ds_solar, fsds

    # Stream 3: Precip
    logger.info("  --- Writing Precip ---")
    ds_prec = base_ds()
    ds_prec['PRECTmms'] = xr.DataArray(prec, dims=['time', 'lat', 'lon'],
        attrs={'mode': 'time-dependent', 'units': 'mm H2O / sec', 'long_name': 'PRECTmms total precipitation'})
    ds_prec.attrs = global_attrs
    write_nc(ds_prec, os.path.join(outdirbase, 'Precip', f"{datafilename}.Prec.{YYYY}-{MM_str}.nc"))
    del ds_prec, prec

    logger.info(f"  Month {MM_str} complete in {time_mod.time()-month_t0:.1f}s")

# ---------------------------------------------------------------------------
# Close input files
# ---------------------------------------------------------------------------
logger.info("=" * 70)
for ds in datasets.values():
    ds.close()

elapsed = time_mod.time() - script_t0
logger.info(f"Year {YYYY} complete. Total elapsed: {elapsed:.1f}s ({elapsed/60:.1f} min)")
