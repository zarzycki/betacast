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
# python rda-oi-years.py --year 2024
#
# This is how you'd invoke a loop!
# for year in $(seq 1982 2024); do python rda-oi-years.py --year "$year" ; done

# Parse command-line options
parser = argparse.ArgumentParser(
    description="Reformat OISST AVHRR daily files into CIME-friendly time series."
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
input_folder = f"/glade/u/home/zarzycki/rda/d277007/avhrr_v2.1/{year}/"
output_folder = "/glade/derecho/scratch/zarzycki/"
DEBUG=False  # This should be false pretty much all the time (see below)
nc_format = 'NETCDF4'  # use NETCDF3_CLASSIC for E3SM, NETCDF4 otherwise

# ==================== DERIVED CONFIGURATION ====================

input_pattern = f"{input_folder.rstrip('/')}/oisst-avhrr-v02r01.{year}*.nc"
output_file = f"{output_folder.rstrip('/')}/sst_reynolds_daily_0.25deg_yr{year}_timeseries.nc"

# Create output directory if it doesn't exist
Path(output_folder).mkdir(parents=True, exist_ok=True)

print(f"Input folder: {input_folder}")
print(f"Year: {year}")
print(f"Output folder: {output_folder}")
print(f"Looking for files: {input_pattern}")

# Find all files for specified year
files = sorted(glob.glob(input_pattern))
print(f"Found {len(files)} files to process")

if not files:
    raise ValueError(f"No files found matching pattern: {input_pattern}")

# Validate that we have a reasonable number of files for a full year
expected_days = 366 if year % 4 == 0 and (year % 100 != 0 or year % 400 == 0) else 365
if len(files) < expected_days - 1:
    print(f"Warning: Only found {len(files)} files, expected ~{expected_days} for year {year}")
else:
    print(f"✓ Found complete year data ({len(files)} files)")

# Get coordinates from first file
print("Reading coordinates from first file...")
ds_first = xr.open_dataset(files[0])
lat_coords = ds_first.lat.values
lon_coords = ds_first.lon.values
ds_first.close()

print(f"Native grid: {len(lat_coords)} x {len(lon_coords)} (lat x lon)")
print(f"Latitude range: {lat_coords[0]:.3f} to {lat_coords[-1]:.3f}")
print(f"Longitude range: {lon_coords[0]:.3f} to {lon_coords[-1]:.3f}")

# Initialize lists to store data
sst_list = []
ice_list = []
time_list = []

# Process each file
print("\nProcessing files...")
for i, file in enumerate(files):
    if (i + 1) % 10 == 0 or i == 0:  # Print progress every 10 files
        print(f"Processing {i+1}/{len(files)}: {Path(file).name}")

    # Load file
    ds = xr.open_dataset(file)

    # Extract SST data (xarray automatically applies scaling)
    sst_data = ds.sst.values.squeeze()  # Remove time and zlev dimensions
    sst_list.append(sst_data)

    # Extract ice data (xarray automatically applies scaling)
    ice_data = ds.ice.values.squeeze()  # Remove time and zlev dimensions
    ice_list.append(ice_data)

    # Process time - convert datetime to days since reference
    time_dt = pd.to_datetime(ds.time.values[0])
    reference_dt = pd.to_datetime(f'{ref_year}-01-01 00:00:00')
    time_days = (time_dt - reference_dt).days + (time_dt - reference_dt).seconds / 86400.0
    time_list.append(time_days)

    ds.close()

print(f"Processed {len(files)} files")

# Convert lists to arrays
print("\nConverting data to arrays...")
sst_array = np.array(sst_list)  # Shape: (time, lat, lon)
ice_array = np.array(ice_list)  # Shape: (time, lat, lon)
time_array = np.array(time_list)

print(f"Data arrays shape: {sst_array.shape}")

# ==================== DEBUG MODIFICATIONS ====================

# This block of code will impose a square box of warming on 8/2 over the NATL
# In 9/2025, this was used to verify calendar attributes and stream reading
if DEBUG:
    print("DEBUG mode enabled - applying test modifications...")

    # Find August 2nd
    ref_date = datetime(ref_year, 1, 1)
    target_date = datetime(year, 8, 2)  # August 2nd of the specified year

    perturb_date_found = False
    for i, t in enumerate(time_array):
        actual_date = ref_date + timedelta(days=float(t))
        if actual_date.month == 8 and actual_date.day == 2:
            print(f"Found August 2nd at time index {i}: {actual_date.strftime('%Y-%m-%d')}")

            # Define North Atlantic patch coordinates (lat: 20-32°N, lon: 297-320°E)
            lat_min_idx = np.argmin(np.abs(lat_coords - 20))    # ~20°N
            lat_max_idx = np.argmin(np.abs(lat_coords - 32))    # ~32°N
            lon_min_idx = np.argmin(np.abs(lon_coords - 297))   # ~297°E
            lon_max_idx = np.argmin(np.abs(lon_coords - 320))   # ~320°E

            print(f"  Applying -3K SST patch to:")
            print(f"  Lat indices {lat_min_idx}-{lat_max_idx} ({lat_coords[lat_min_idx]:.2f}°-{lat_coords[lat_max_idx]:.2f}°)")
            print(f"  Lon indices {lon_min_idx}-{lon_max_idx} ({lon_coords[lon_min_idx]:.2f}°-{lon_coords[lon_max_idx]:.2f}°)")

            # Apply -3K cooling to the patch
            sst_array[i, lat_min_idx:lat_max_idx+1, lon_min_idx:lon_max_idx+1] -= 5.0

            perturb_date_found = True
            break

    if not perturb_date_found:
        print("Warning: August 2nd not found in the time series")

# ==================== END DEBUG MODIFICATIONS ====================

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
    'title': f'OISST AVHRR daily data reformatted for CESM/E3SM - Year {year}',
    'source_data': 'NOAA/NCEI 1/4 Degree Daily OISST Analysis v2.1',
    'source_folder': str(Path(input_folder).resolve()),
    'native_grid': '0.25 degree global grid',
    'time_range': f'{Path(files[0]).name[-14:-3]} to {Path(files[-1]).name[-14:-3]}',
    'creation_date': datetime.now().strftime('%Y-%m-%d %H:%M:%S'),
    'processing': 'Variable renaming and format conversion only',
    'year_processed': year,
    'files_processed': len(files)
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
print(f"✓ Year {year}: {len(files)} days from {Path(files[0]).name[-14:-3]} to {Path(files[-1]).name[-14:-3]}")
print("✓ Conversion completed!")
print(f"✓ Output file size: {Path(output_file).stat().st_size / (1024**3):.2f} GB")
