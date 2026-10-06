import numpy as np
import xarray as xr
from pathlib import Path
import glob
import os
import sys
from datetime import datetime, timedelta
import pandas as pd
import argparse

module_path = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', 'py_functions'))
if module_path not in sys.path:
    sys.path.append(module_path)
import pyfuncs

# This is how you'd invoke a single year
# python gen-rda-era5sst-year.py --year 2024
#
# This is how you'd invoke a loop!
# for year in $(seq 1982 2024); do python gen-rda-era5sst-year.py --year "$year" ; done

# Parse command-line options
parser = argparse.ArgumentParser(
    description="Reformat ERA5 hourly SST/ice (NCAR GDEX d633000) into CIME-friendly daily time series."
)
parser.add_argument(
    "--year", type=int, required=True,
    help="4-digit year to process (e.g., 2010)"
)
args = parser.parse_args()

def create_date_variables(time_days, reference_year=1970):
    """Create date and datesec variables from time"""
    # Convert to actual dates
    ref_date = datetime(reference_year, 1, 1)
    dates = []
    datesecs = []

    for t in time_days:
        actual_date = ref_date + timedelta(days=float(t))
        # date: YYYYMMDD format
        date_int = int(actual_date.strftime('%Y%m%d'))
        # datesec: seconds since start of day (should be 0 for daily data)
        datesec_int = 0  # Daily data centered at start of day

        dates.append(date_int)
        datesecs.append(datesec_int)

    return np.array(dates), np.array(datesecs)

# ==================== USER CONFIGURATION ====================

ref_year=1970    # This is the "years since" identifier
TTHRESH = 271.9  # Initial cut for ice vs. open ocean
KtoC = 273.15    # Kelvin to Celsius conversion
SMOOTH_ICE = True
SMOOTH_ITER = 3
year = args.year
input_folder = "/glade/campaign/collections/gdex/data/d633000/e5.oper.an.sfc/"
output_folder = "/glade/derecho/scratch/zarzycki/"
nc_format = 'NETCDF4'  # use NETCDF3_CLASSIC for E3SM, NETCDF4 otherwise

# ==================== DERIVED CONFIGURATION ====================

output_file = f"{output_folder.rstrip('/')}/sst_era5_daily_0.25deg_yr{year}_timeseries.nc"

# Create output directory if it doesn't exist
Path(output_folder).mkdir(parents=True, exist_ok=True)

print(f"Input folder: {input_folder}")
print(f"Year: {year}")
print(f"Output folder: {output_folder}")

# Find the monthly SST (128_034_sstk) and ice (128_031_ci) files for specified year
sst_files = []
ice_files = []
for month in range(1, 13):
    month_folder = f"{input_folder.rstrip('/')}/{year}{month:02d}"
    sst_match = sorted(glob.glob(f"{month_folder}/e5.oper.an.sfc.128_034_sstk.ll025sc.{year}{month:02d}*.nc"))
    ice_match = sorted(glob.glob(f"{month_folder}/e5.oper.an.sfc.128_031_ci.ll025sc.{year}{month:02d}*.nc"))
    if len(sst_match) != 1 or len(ice_match) != 1:
        raise ValueError(f"Expected 1 sstk and 1 ci file in {month_folder}, found {len(sst_match)} and {len(ice_match)}")
    sst_files.append(sst_match[0])
    ice_files.append(ice_match[0])

print(f"Found {len(sst_files)} monthly sstk files and {len(ice_files)} monthly ci files")

# Get coordinates from first file
# ERA5 is on 0.25deg nodes (lat 90 -> -90, 721 pts, includes poles; lon 0 -> 359.75)
print("Reading coordinates from first file...")
ds_first = xr.open_dataset(sst_files[0])
era5_lat = ds_first.latitude.values
era5_lon = ds_first.longitude.values
ds_first.close()

print(f"Native ERA5 grid: {len(era5_lat)} x {len(era5_lon)} (lat x lon)")
print(f"Latitude range: {era5_lat[0]:.3f} to {era5_lat[-1]:.3f}")
print(f"Longitude range: {era5_lon[0]:.3f} to {era5_lon[-1]:.3f}")

# Initialize lists to store data
sst_list = []
ice_list = []
time_list = []

# Process each month, then each day within the month
# Each daily value is the mean of the 24 hourly analyses (00Z-23Z) on that day
print("\nProcessing files...")
for i in range(len(sst_files)):
    print(f"Processing {i+1}/{len(sst_files)}: {Path(sst_files[i]).name}")

    # Load files (xarray converts the 9.999e20 land fill to NaN)
    ds_sst = xr.open_dataset(sst_files[i])
    ds_ice = xr.open_dataset(ice_files[i])

    # Get the unique days in this month
    days_in_month = pd.to_datetime(ds_sst.time.values).normalize().unique()

    for day in days_in_month:
        day_start = day
        day_end = day + pd.Timedelta(hours=23)

        sst_hourly = ds_sst.SSTK.sel(time=slice(day_start, day_end))
        ice_hourly = ds_ice.CI.sel(time=slice(day_start, day_end))

        if sst_hourly.sizes['time'] != 24 or ice_hourly.sizes['time'] != 24:
            print(f"Warning: {day.strftime('%Y-%m-%d')} has {sst_hourly.sizes['time']} sstk and {ice_hourly.sizes['time']} ci hours (expected 24)")

        # Daily mean, convert SST K -> C
        sst_data = sst_hourly.mean(dim='time').values - KtoC
        ice_data = ice_hourly.mean(dim='time').values
        sst_list.append(sst_data)
        ice_list.append(ice_data)

        # Process time - convert datetime to days since reference
        reference_dt = pd.to_datetime(f'{ref_year}-01-01 00:00:00')
        time_days = (day - reference_dt).days + (day - reference_dt).seconds / 86400.0
        time_list.append(time_days)

    ds_sst.close()
    ds_ice.close()

print(f"Processed {len(time_list)} days")

# Validate that we have a reasonable number of days for a full year
expected_days = 366 if year % 4 == 0 and (year % 100 != 0 or year % 400 == 0) else 365
if len(time_list) != expected_days:
    print(f"Warning: Only found {len(time_list)} days, expected {expected_days} for year {year}")
else:
    print(f"✓ Found complete year data ({len(time_list)} days)")

# Convert lists to arrays
print("\nConverting data to arrays...")
sst_array = np.array(sst_list)  # Shape: (time, lat, lon)
ice_array = np.array(ice_list)  # Shape: (time, lat, lon)
time_array = np.array(time_list)

print(f"ERA5 data arrays shape: {sst_array.shape}")

# ==================== MOVE TO OISST CELL CENTERS ====================

# ERA5 points sit on the 0.25deg nodes, but OISST (and domains/domain.ocn.720x1440.nc)
# are cell centered (lat -89.875 -> 89.875, lon 0.125 -> 359.875). Here we average
# each pair of neighboring ERA5 points so we can reuse the existing 720x1440 domain.
# NOTE: if either neighbor is NaN (i.e., land) the average is NaN, which then gets
# filled below just like OISST land points.
print("Averaging ERA5 nodes to 720x1440 cell centers...")

# Flip lat so it goes S -> N like OISST
era5_lat = era5_lat[::-1]
sst_array = sst_array[:, ::-1, :]
ice_array = ice_array[:, ::-1, :]

# Lat: average row j and j+1 (721 -> 720)
lat_coords = 0.5 * (era5_lat[:-1] + era5_lat[1:])
sst_array = 0.5 * (sst_array[:, :-1, :] + sst_array[:, 1:, :])
ice_array = 0.5 * (ice_array[:, :-1, :] + ice_array[:, 1:, :])

# Lon: average column i and i+1, wrapping 359.75 with 0 (1440 -> 1440)
lon_coords = era5_lon + 0.125
sst_array = 0.5 * (sst_array + np.roll(sst_array, -1, axis=2))
ice_array = 0.5 * (ice_array + np.roll(ice_array, -1, axis=2))

print(f"Data arrays shape: {sst_array.shape}")
print(f"Latitude range: {lat_coords[0]:.3f} to {lat_coords[-1]:.3f}")
print(f"Longitude range: {lon_coords[0]:.3f} to {lon_coords[-1]:.3f}")

# This will first essentially using linear interpolation to fill nans along lats
print(f"Filling missing data...")
sst_array = pyfuncs.linmsg_n(sst_array, -1, 2, fill_all_nans=(TTHRESH-KtoC))
sst_array = pyfuncs.linmsg_n(sst_array, -1, 1, fill_all_nans=(TTHRESH-KtoC))

# Smooth ice if requested
if SMOOTH_ICE:
    print(f"smoothing derived ice field {SMOOTH_ITER} times!")
    ice_array = pyfuncs.smooth_with_smth9(ice_array, SMOOTH_ITER, p=0.50, q=0.25)

# Handling missing values (there really should be much/any, but let's just do it)
print(f"Final final missing values!")
sst_array = np.where(np.isnan(sst_array), (TTHRESH-KtoC), sst_array)
sst_array = np.where(sst_array > 500, (TTHRESH-KtoC), sst_array)
ice_array = np.where(np.isnan(ice_array), 0.0, ice_array)
ice_array = np.where(ice_array > 500, 1.0, ice_array)

# Create date variables
print("Creating date variables...")
date_array, datesec_array = create_date_variables(time_array, reference_year=ref_year)

# Create output dataset
# Admittedly, this is done to mimic some existing CESM SST files
# NOTE: The calendar attribute is very important!
# "standard" includes leap years, "noleap" is the CESM default without the calendar
# attribute and will cause Feb 29 to drive out-of-sync conditions
print(f"Creating output dataset: {output_file}")
ds_out = xr.Dataset(
    {
        'SST_cpl': (['time', 'lat', 'lon'], sst_array, {
            '_FillValue': -999.0,
            '_FillValue_original': -999,
            'units': 'degrees C',
            'long_name': 'Daily sea surface temperature'
        }),
        'ice_cov': (['time', 'lat', 'lon'], ice_array, {
            '_FillValue': -999.0,
            '_FillValue_original': -999,
            'units': 'percentage',
            'long_name': 'Sea ice concentration'
        }),
        'time': (['time'], time_array, {
            'long_name': 'Center time of the day',
            'calendar': 'standard',
            'units': f'days since {ref_year}-01-01 00:00:00'
        }),
        'date': (['time'], date_array, {
            'long_name': 'current date (YYYYMMDD)'
        }),
        'datesec': (['time'], datesec_array, {
            'long_name': 'current seconds of current date'
        })
    },
    coords={
        'lat': (['lat'], lat_coords, {
            'units': 'degrees_north',
            'long_name': 'latitude',
            '_FillValue': -900.0
        }),
        'lon': (['lon'], lon_coords, {
            'units': 'degrees_east',
            'long_name': 'longitude',
            '_FillValue': -900.0
        })
    }
)

# Add global attributes
ds_out.attrs = {
    'title': f'ERA5 daily mean SST/ice reformatted for CESM/E3SM - Year {year}',
    'source_data': 'ECMWF ERA5 hourly surface analysis (NCAR GDEX d633000, 128_034_sstk and 128_031_ci)',
    'source_folder': str(Path(input_folder).resolve()),
    'native_grid': '0.25 degree global grid (721 x 1440 nodes)',
    'output_grid': '0.25 degree global grid (720 x 1440 cell centers, matches OISST)',
    'time_range': f'{Path(sst_files[0]).name[-24:-3]} to {Path(sst_files[-1]).name[-24:-3]}',
    'creation_date': datetime.now().strftime('%Y-%m-%d %H:%M:%S'),
    'processing': 'Daily mean of hourly data, 2-point average to cell centers, variable renaming and format conversion',
    'year_processed': year,
    'days_processed': len(time_list)
}

# Write to file with compression
print("Writing to NetCDF file...")
encoding = {
    'time': {'dtype': 'float32'},
    'lat': {'dtype': 'float64'},
    'lon': {'dtype': 'float64'},
    'date': {'dtype': 'int32'},
    'datesec': {'dtype': 'int32'}
}
# update encoding to support compression/chunking in netcdf4
if nc_format[:7].upper() == "NETCDF4":
    compression = {'zlib': True, 'complevel': 1}
    encoding.update({
        'SST_cpl': compression,
        'ice_cov': compression
    })

ds_out.to_netcdf(output_file, format=nc_format, encoding=encoding, unlimited_dims=['time'])
ds_out.close()

print(f"\n✓ Successfully created: {output_file}")
print(f"✓ Output dimensions: time={len(time_array)}, lat={len(lat_coords)}, lon={len(lon_coords)}")
print(f"✓ Year {year}: {len(time_list)} days from {date_array[0]} to {date_array[-1]}")
print("✓ Conversion completed!")
print(f"✓ Output file size: {Path(output_file).stat().st_size / (1024**3):.2f} GB")
