#!/usr/bin/env python3
import numpy as np
import xarray as xr
import os
import sys
import argparse
import logging
import cftime
from datetime import datetime

# Import provided modules
import meteo
import horizremap
from constants import grav, Rd, gamma_s, p0

# Set up logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

# Parse command line arguments
parser = argparse.ArgumentParser(description="Process ERA5 data for DATM forcing")
parser.add_argument("--era5_file", required=True, help="Path to ERA5 input file")
parser.add_argument("--year", required=True, help="Year (YYYY)")
parser.add_argument("--month", required=True, help="Month (MM)")
parser.add_argument("--domain_file", required=True, help="Path to domain file")
parser.add_argument("--newgrid", action="store_true", help="Regrid to new grid")
parser.add_argument("--wgt_filename", default="", help="Path to ESMF weight file")
parser.add_argument("--do_q", action="store_true", help="Calculate specific humidity")
parser.add_argument("--do_flds", action="store_true", help="Include longwave downward radiation")
parser.add_argument("--outdirbase", required=True, help="Base output directory")
parser.add_argument("--datafilename", default="CMZERA5.v0.c2021.0.5d", help="Base name for output files")
parser.add_argument("--greg_to_noleap", action="store_true", help="Convert from Gregorian to no-leap calendar")

args = parser.parse_args()

RAWERA5FILE = args.era5_file
YYYY = args.year
MM = args.month
newgrid = args.newgrid
wgt_filename = args.wgt_filename
do_q = args.do_q
do_flds = args.do_flds
outdirbase = args.outdirbase
datafilename = args.datafilename
greg_to_noleap = args.greg_to_noleap

# Open the domain file
try:
    logger.info(f"Opening domain file: {args.domain_file}")
    d = xr.open_dataset(args.domain_file)
except Exception as e:
    logger.error(f"Error opening domain file: {e}")
    sys.exit(1)

# Determine if it's a leap month
is_leap_month = False
if ((int(YYYY) % 4 == 0 and int(YYYY) % 100 != 0) or int(YYYY) % 400 == 0) and int(MM) == 2:
    logger.info(f"We have a February leap year for {YYYY} {MM}")
    is_leap_month = True

try:
    logger.info(f"Adding {RAWERA5FILE} ...")
    f = xr.open_dataset(RAWERA5FILE)
except Exception as e:
    logger.error(f"Error opening ERA5 file: {e}")
    sys.exit(1)

# Get coordinates from domain file
logger.info("Getting coordinates from domain file...")
lon_datm_2D = d.xc
lat_datm_2D = d.yc

# Assuming domain is regular lat-lon grid, pull first lat/lon to get 1D lat/lon arrays
lon_datm = lon_datm_2D.isel(lat=0).values
lat_datm = lat_datm_2D.isel(lon=0).values

# Get coordinates from ERA5 file
logger.info("Reading relevant coordinates from reanalysis...")
lon_era5 = f.longitude.values
lat_era5 = f.latitude.values[::-1]  # reversed to match NCL [::-1]

# Determine ENIX (end index)
logger.info("Figuring out ENIX")
STIX = 0  # Start index
ENIX = len(f.time) - 1
timestride = 1

if is_leap_month and greg_to_noleap:
    # If this is Feb/leap and we want no leap, find days with DD "29" and cut off
    logger.info(f"Correcting ENIX from {ENIX} to...")

    # Try to decode time values
    try:
        tmp_times = f.time.to_index()
        tmp_days = tmp_times.day.values
    except:
        # Fallback: Extract day from time value string representation
        tmp_times_str = [str(t) for t in f.time.values]
        tmp_days = [int(t.split('-')[2].split('T')[0]) if 'T' in t else int(t.split('-')[2]) for t in tmp_times_str]

    ENIX = sum(np.array(tmp_days) <= 28) - 1
    logger.info(f"... to {ENIX}")

# Get full time
time_era5 = f.time[STIX:ENIX+1:timestride]
logger.info(f"Time range: {time_era5.values}")

# Extract variables from ERA5 file
logger.info("Pull required variables off files...")
u10_era5 = f.u10[STIX:ENIX+1:timestride, ::-1, :].values.astype(np.float32)
v10_era5 = f.v10[STIX:ENIX+1:timestride, ::-1, :].values.astype(np.float32)
tbot_era5 = f.t2m[STIX:ENIX+1:timestride, ::-1, :].values.astype(np.float32)
tdew_era5 = f.d2m[STIX:ENIX+1:timestride, ::-1, :].values.astype(np.float32)
psrf_era5 = f.sp[STIX:ENIX+1:timestride, ::-1, :].values.astype(np.float32)
prec_era5 = f.mtpr[STIX:ENIX+1:timestride, ::-1, :].values.astype(np.float32)
fsds_era5 = f.ssrd[STIX:ENIX+1:timestride, ::-1, :].values.astype(np.float32)
if do_flds:
    flds_era5 = f.strd[STIX:ENIX+1:timestride, ::-1, :].values.astype(np.float32)

# Print stats
logger.info("READ FROM ERA5 STATS:")
logger.info(f"u10_era5 max: {np.max(u10_era5)}   min: {np.min(u10_era5)}")
logger.info(f"v10_era5 max: {np.max(v10_era5)}   min: {np.min(v10_era5)}")
logger.info(f"tbot_era5 max: {np.max(tbot_era5)}   min: {np.min(tbot_era5)}")
logger.info(f"tdew_era5 max: {np.max(tdew_era5)}   min: {np.min(tdew_era5)}")
logger.info(f"psrf_era5 max: {np.max(psrf_era5)}   min: {np.min(psrf_era5)}")
logger.info(f"prec_era5 max: {np.max(prec_era5)}   min: {np.min(prec_era5)}")
logger.info(f"fsds_era5 max: {np.max(fsds_era5)}   min: {np.min(fsds_era5)}")
if do_flds:
    logger.info(f"flds_era5 max: {np.max(flds_era5)}   min: {np.min(flds_era5)}")

# Process wind
logger.info("Processing wind...")
wind_era5 = np.sqrt(u10_era5**2 + v10_era5**2)

# Process specific humidity if requested
if do_q:
    logger.info("Processing Q...")
    qbot_era5 = meteo.mixhum_ptd(psrf_era5, tdew_era5, 2)

# Interpolate to new grid if requested
if newgrid:
    if wgt_filename:
        logger.info("Using ESMF mapping")

        # Use the provided remap_with_weights_wrapper function for regridding
        tbot_datm, lat_datm, lon_datm = horizremap.remap_with_weights_wrapper(tbot_era5, wgt_filename)
        if do_q:
            qbot_datm, _, _ = horizremap.remap_with_weights_wrapper(qbot_era5, wgt_filename)
        else:
            qbot_datm, _, _ = horizremap.remap_with_weights_wrapper(tdew_era5, wgt_filename)
        wind_datm, _, _ = horizremap.remap_with_weights_wrapper(wind_era5, wgt_filename)
        psrf_datm, _, _ = horizremap.remap_with_weights_wrapper(psrf_era5, wgt_filename)
        prec_datm, _, _ = horizremap.remap_with_weights_wrapper(prec_era5, wgt_filename)
        fsds_datm, _, _ = horizremap.remap_with_weights_wrapper(fsds_era5, wgt_filename)
        if do_flds:
            flds_datm, _, _ = horizremap.remap_with_weights_wrapper(flds_era5, wgt_filename)
    else:
        logger.info("Bilinear interpolation")

        # Use bilinear interpolation
        tbot_datm = horizremap.linint2(lon_era5, lat_era5, tbot_era5, True, lon_datm, lat_datm)
        if do_q:
            qbot_datm = horizremap.linint2(lon_era5, lat_era5, qbot_era5, True, lon_datm, lat_datm)
        else:
            qbot_datm = horizremap.linint2(lon_era5, lat_era5, tdew_era5, True, lon_datm, lat_datm)
        wind_datm = horizremap.linint2(lon_era5, lat_era5, wind_era5, True, lon_datm, lat_datm)
        psrf_datm = horizremap.linint2(lon_era5, lat_era5, psrf_era5, True, lon_datm, lat_datm)
        prec_datm = horizremap.linint2(lon_era5, lat_era5, prec_era5, True, lon_datm, lat_datm)
        fsds_datm = horizremap.linint2(lon_era5, lat_era5, fsds_era5, True, lon_datm, lat_datm)
        if do_flds:
            flds_datm = horizremap.linint2(lon_era5, lat_era5, flds_era5, True, lon_datm, lat_datm)
else:
    logger.info("Copying vars")
    tbot_datm = tbot_era5.copy()
    del tbot_era5
    if do_q:
        qbot_datm = qbot_era5.copy()
        del qbot_era5
    else:
        qbot_datm = tdew_era5.copy()
        del tdew_era5
    wind_datm = wind_era5.copy()
    del wind_era5
    psrf_datm = psrf_era5.copy()
    del psrf_era5
    prec_datm = prec_era5.copy()
    del prec_era5
    fsds_datm = fsds_era5.copy()
    del fsds_era5
    if do_flds:
        flds_datm = flds_era5.copy()
        del flds_era5

# Enforce zero checks
logger.info("Enforce zero checks")
eps = 1.0e-8
prec_datm = np.where(prec_datm < 0.0, 0.0, prec_datm)
fsds_datm = np.where(fsds_datm <= 0.0, eps, fsds_datm)
if do_flds:
    flds_datm = np.where(flds_datm <= 0.0, eps, flds_datm)
if do_q:
    qbot_datm = np.where(qbot_datm <= 0.0, eps, qbot_datm)

# Enforce specified capping checks
logger.info("Enforce specified capping checks")

# Cap wind
maxwind = 45.0
minwind = 0.0
logger.info(f"wind_datm max: {np.max(wind_datm)}   min: {np.min(wind_datm)}")
logger.info(f"min/max set by user: {minwind} {maxwind}")
logger.info(f"Number of wind_datm over max to be corrected: {np.sum(wind_datm > maxwind)}")
logger.info(f"Number of wind_datm under min to be corrected: {np.sum(wind_datm < minwind)}")
wind_datm = np.where(wind_datm > maxwind, maxwind, wind_datm)
wind_datm = np.where(wind_datm < minwind, minwind, wind_datm)

# Cap pressure
maxpsrf = 120000.0
minpsrf = 30000.0
logger.info(f"psrf_datm max: {np.max(psrf_datm)}   min: {np.min(psrf_datm)}")
logger.info(f"min/max set by user: {minpsrf} {maxpsrf}")
logger.info(f"Number of psrf_datm over max to be corrected: {np.sum(psrf_datm > maxpsrf)}")
logger.info(f"Number of psrf_datm under min to be corrected: {np.sum(psrf_datm < minpsrf)}")
psrf_datm = np.where(psrf_datm > maxpsrf, maxpsrf, psrf_datm)
psrf_datm = np.where(psrf_datm < minpsrf, minpsrf, psrf_datm)

# Cap temperature
maxtbot = 340.0
mintbot = 175.0
logger.info(f"tbot_datm max: {np.max(tbot_datm)}   min: {np.min(tbot_datm)}")
logger.info(f"min/max set by user: {mintbot} {maxtbot}")
logger.info(f"Number of tbot_datm over max to be corrected: {np.sum(tbot_datm > maxtbot)}")
logger.info(f"Number of tbot_datm under min to be corrected: {np.sum(tbot_datm < mintbot)}")
tbot_datm = np.where(tbot_datm > maxtbot, maxtbot, tbot_datm)
tbot_datm = np.where(tbot_datm < mintbot, mintbot, tbot_datm)

# Cap specific humidity if requested
if do_q:
    maxqbot = 0.2
    logger.info(f"qbot_datm max: {np.max(qbot_datm)}   min: {np.min(qbot_datm)}")
    logger.info(f"max set by user: {maxqbot}")
    logger.info(f"Number of qbot_datm over max to be corrected: {np.sum(qbot_datm > maxqbot)}")
    qbot_datm = np.where(qbot_datm > maxqbot, maxqbot, qbot_datm)

# Create dummy variables
logger.info("Dummy variables")
zbot_datm = np.full_like(tbot_datm, 10.0)
zbot_datm_attrs = {"long_name": "reference height", "mode": "time-dependent", "units": "m"}

# Convert units
logger.info("Convert units")
fsds_datm = fsds_datm / 3600.  # convert from J/s to W/m2 (over 1 hour)
if do_flds:
    flds_datm = flds_datm / 3600.

# Set time coordinates
logger.info("Set time units")
if greg_to_noleap:
    new_time_units = "days since 1900-01-01 00:00:00"

    # Get dates from time index
    try:
        tmp_dates = time_era5.to_index()
        utc_date = [(d.year, d.month, d.day, d.hour, d.minute, d.second) for d in tmp_dates]
    except:
        # Fallback: try to parse the time strings
        logger.warning("Failed to decode time with to_index(), trying cftime")
        utc_date = []
        time_units = getattr(time_era5, 'units', 'days since 1900-01-01')
        time_calendar = getattr(time_era5, 'calendar', 'standard')
        try:
            dates = cftime.num2date(time_era5.values, units=time_units, calendar=time_calendar)
            utc_date = [(d.year, d.month, d.day, d.hour, d.minute, d.second) for d in dates]
        except:
            logger.error("Failed to decode time, using fallback values")
            # Last resort: use fixed values for demonstration
            utc_date = [(int(YYYY), int(MM), 1, 0, 0, 0)]

    # Convert to noleap calendar
    time_datm = []
    for year, month, day, hour, minute, second in utc_date:
        date_noleap = cftime.DatetimeNoLeap(year, month, day, hour, minute, second)
        days_since = cftime.date2num(date_noleap, units=new_time_units, calendar="noleap")
        time_datm.append(days_since)
    time_datm = np.array(time_datm, dtype=np.float64)
else:
    new_time_units = "days since 1986-08-01 00:00:00"

    # Convert using standard calendar
    try:
        tmp_dates = time_era5.to_index()
        # Convert to new reference date
        time_datm = [(d - np.datetime64('1986-08-01')).total_seconds() / 86400 for d in tmp_dates]
        time_datm = np.array(time_datm, dtype=np.float64)
    except:
        # Fallback
        logger.warning("Failed to decode time with standard methods, using cftime")
        time_units = getattr(time_era5, 'units', 'days since 1900-01-01')
        time_calendar = getattr(time_era5, 'calendar', 'standard')
        try:
            dates = cftime.num2date(time_era5.values, units=time_units, calendar=time_calendar)
            time_datm = cftime.date2num(dates, units=new_time_units, calendar="standard")
        except:
            logger.error("Failed to decode time, using original values")
            time_datm = time_era5.values.astype(np.float64)

# Set metadata for coordinates
lat_datm_float = lat_datm.astype(np.float32)
lon_datm_float = lon_datm.astype(np.float32)

lat_datm_float_attrs = {
    "long_name": "latitude",
    "units": "degrees_north",
    "mode": "time-invariant"
}
lon_datm_float_attrs = {
    "long_name": "longitude",
    "units": "degrees_east",
    "mode": "time-invariant"
}

# 2D lats and lons
lat_datm_2D_float = lat_datm_2D.values.astype(np.float32)
lon_datm_2D_float = lon_datm_2D.values.astype(np.float32)

lat_datm_2D_float_attrs = {
    "long_name": "latitude",
    "units": "degrees_north",
    "mode": "time-invariant"
}
lon_datm_2D_float_attrs = {
    "long_name": "longitude",
    "units": "degrees_east",
    "mode": "time-invariant"
}

# Set up edge variables and their attributes
EDGEW = np.array([0.])
EDGEE = np.array([360.])
EDGES = np.array([-90.])
EDGEN = np.array([90.])

EDGEW_attrs = {
    "long_name": "western edge in atmospheric data",
    "mode": "time-invariant"
}
EDGEE_attrs = {
    "long_name": "eastern edge in atmospheric data",
    "mode": "time-invariant"
}
EDGES_attrs = {
    "long_name": "southern edge in atmospheric data",
    "mode": "time-invariant"
}
EDGEN_attrs = {
    "long_name": "northern edge in atmospheric data",
    "mode": "time-invariant"
}

# Set metadata for other variables
tbot_datm_attrs = {
    "mode": "time-dependent",
    "units": "K",
    "long_name": "temperature at the lowest model level"
}
wind_datm_attrs = {
    "mode": "time-dependent",
    "units": "m/s",
    "long_name": "wind magnitude at the lowest model level"
}
psrf_datm_attrs = {
    "mode": "time-dependent",
    "units": "Pa",
    "long_name": "surface pressure"
}
qbot_datm_attrs = {
    "mode": "time-dependent",
    "units": "kg/kg" if do_q else "K",
    "long_name": "specific humidity at lowest model level" if do_q else "dew point temperature at 2m"
}
prec_datm_attrs = {
    "mode": "time-dependent",
    "units": "mm H2O / sec",
    "long_name": "PRECTmms total precipitation"
}
fsds_datm_attrs = {
    "mode": "time-dependent",
    "units": "W/m**2",
    "long_name": "downward solar flux at surface"
}
time_datm_attrs = {
    "long_name": "observation time",
    "units": new_time_units,
    "calendar": "noleap" if greg_to_noleap else "standard"
}

if do_flds:
    flds_datm_attrs = {
        "mode": "time-dependent",
        "units": "W/m**2",
        "long_name": "downward longwave flux at surface"
    }

# Define output folders and file types
filenames = ["TPQWL", "Solar", "Prec"]
outfolders = ["TPQW", "Solar", "Precip"]

# Create output files
for ii in range(len(filenames)):
    fullfilename = os.path.join(outdirbase, outfolders[ii],
                              f"{datafilename}.{filenames[ii]}.{YYYY}-{MM}.nc")

    logger.info(f"Writing {fullfilename} output!")

    # Create output directory if it doesn't exist
    os.makedirs(os.path.join(outdirbase, outfolders[ii]), exist_ok=True)

    # Remove any pre-existing file
    if os.path.exists(fullfilename):
        os.remove(fullfilename)

    # Create datasets with dimensions
    # Extract shapes for dimension setup
    if newgrid:
        time_dim = len(time_datm)
        lat_dim = len(lat_datm)
        lon_dim = len(lon_datm)
        scalar_dim = 1
    else:
        time_shape = tbot_datm.shape[0]
        lat_shape = tbot_datm.shape[1]
        lon_shape = tbot_datm.shape[2]
        time_dim = time_shape
        lat_dim = lat_shape
        lon_dim = lon_shape
        scalar_dim = 1

    # Create a new dataset
    ds = xr.Dataset()

    # Add time dimension and coordinate
    ds['time'] = xr.DataArray(time_datm, dims=['time'])
    ds['time'].attrs = time_datm_attrs

    # Add lat/lon coordinates
    ds['LONGXY'] = xr.DataArray(lon_datm_2D_float, dims=['lat', 'lon'])
    ds['LONGXY'].attrs = lon_datm_2D_float_attrs

    ds['LATIXY'] = xr.DataArray(lat_datm_2D_float, dims=['lat', 'lon'])
    ds['LATIXY'].attrs = lat_datm_2D_float_attrs

    # Add edge variables
    ds['EDGEW'] = xr.DataArray(EDGEW, dims=['scalar'])
    ds['EDGEW'].attrs = EDGEW_attrs

    ds['EDGEE'] = xr.DataArray(EDGEE, dims=['scalar'])
    ds['EDGEE'].attrs = EDGEE_attrs

    ds['EDGES'] = xr.DataArray(EDGES, dims=['scalar'])
    ds['EDGES'].attrs = EDGES_attrs

    ds['EDGEN'] = xr.DataArray(EDGEN, dims=['scalar'])
    ds['EDGEN'].attrs = EDGEN_attrs

    # Add data variables based on file type
    if ii == 0:  # TPQWL file
        ds['TBOT'] = xr.DataArray(tbot_datm, dims=['time', 'lat', 'lon'])
        ds['TBOT'].attrs = tbot_datm_attrs

        ds['WIND'] = xr.DataArray(wind_datm, dims=['time', 'lat', 'lon'])
        ds['WIND'].attrs = wind_datm_attrs

        ds['PSRF'] = xr.DataArray(psrf_datm, dims=['time', 'lat', 'lon'])
        ds['PSRF'].attrs = psrf_datm_attrs

        var_name = 'QBOT' if do_q else 'TDEW'
        ds[var_name] = xr.DataArray(qbot_datm, dims=['time', 'lat', 'lon'])
        ds[var_name].attrs = qbot_datm_attrs

        ds['ZBOT'] = xr.DataArray(zbot_datm, dims=['time', 'lat', 'lon'])
        ds['ZBOT'].attrs = zbot_datm_attrs

        if do_flds:
            ds['FLDS'] = xr.DataArray(flds_datm, dims=['time', 'lat', 'lon'])
            ds['FLDS'].attrs = flds_datm_attrs

    elif ii == 1:  # Solar file
        ds['FSDS'] = xr.DataArray(fsds_datm, dims=['time', 'lat', 'lon'])
        ds['FSDS'].attrs = fsds_datm_attrs

    elif ii == 2:  # Prec file
        ds['PRECTmms'] = xr.DataArray(prec_datm, dims=['time', 'lat', 'lon'])
        ds['PRECTmms'].attrs = prec_datm_attrs

    # Add global attributes
    ds.attrs = {
        'title': f"DATM forcing using ERA5: {YYYY} {MM}",
        'source_file': RAWERA5FILE,
        'Conventions': "None",
        'creation_date': datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
        'notes': "Generated with Betacast toolkit"
    }

    # Set up encoding for compression if desired
    encoding = {}
    for var in ds.variables:
        encoding[var] = {'zlib': True, 'complevel': 1}

    # Write to netCDF file
    ds.to_netcdf(
        fullfilename,
        encoding=encoding,
        format='NETCDF4'
    )

    logger.info(f"Successfully wrote {fullfilename}")

logger.info("Processing completed successfully!")
return 0