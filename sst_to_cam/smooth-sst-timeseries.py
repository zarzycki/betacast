import numpy as np
from netCDF4 import Dataset
from pathlib import Path
from datetime import datetime
import argparse

# This is how you'd invoke this script
# python smooth-sst-timeseries.py file.nc --smooth_window 5
#
# Applies a centered running mean in time to SST_cpl and ice_cov at each (lat, lon).
# Near the start/end of the time series the window shrinks so it stays centered
# (e.g., with --smooth_window 5, day 2 is a 3-day mean and days 1 and N are untouched).
# Output is written next to the input file as file_smooth5.nc
#
# NOTE: The data is processed in blocks of BLOCK_SIZE days so that long (e.g., 1979-2025)
# 0.25deg files don't need to fit in memory. Each block reads +/- half_window extra
# days on either side so the result is identical to smoothing everything at once.

# Parse command-line options
parser = argparse.ArgumentParser(
    description="Apply a centered running mean in time to SST_cpl and ice_cov."
)
parser.add_argument(
    "input_file", type=str,
    help="SST/ice timeseries file (e.g., sst_era5_daily_0.25deg_yr2005_timeseries.nc)"
)
parser.add_argument(
    "--smooth_window", type=int, required=True,
    help="Number of days in running mean (must be odd, e.g., 5)"
)
args = parser.parse_args()

# ==================== CONFIGURATION ====================

smooth_window = args.smooth_window
input_file = args.input_file
vars_to_smooth = ['SST_cpl', 'ice_cov']
BLOCK_SIZE = 365  # days per block, ~1.5GB per var per block at 0.25deg

if smooth_window < 1 or smooth_window % 2 == 0:
    raise ValueError(f"--smooth_window must be a positive odd integer, got {smooth_window}")

half_window = smooth_window // 2

input_path = Path(input_file)
output_file = str(input_path.with_name(f"{input_path.stem}_smooth{smooth_window}{input_path.suffix}"))

print(f"Input file: {input_file}")
print(f"Output file: {output_file}")
print(f"Smoothing window: {smooth_window} days (+/- {half_window})")

# ==================== CREATE OUTPUT FILE ====================

# Output is an exact copy of the input file structure (same format, dims, attributes,
# compression, chunking), then we fill in the variables below
src = Dataset(input_file, 'r')
src.set_auto_mask(False)

print(f"Input format: {src.data_model}")
dst = Dataset(output_file, 'w', format=src.data_model)

# Copy dimensions
for name, dim in src.dimensions.items():
    dst.createDimension(name, None if dim.isunlimited() else len(dim))

# Copy variable definitions and attributes (_FillValue has to be set at creation)
for name, var in src.variables.items():
    fill_value = var.getncattr('_FillValue') if '_FillValue' in var.ncattrs() else None
    if src.data_model.startswith('NETCDF4'):
        filters = var.filters()
        chunking = var.chunking()
        dst_var = dst.createVariable(name, var.dtype, var.dimensions, fill_value=fill_value,
                                     zlib=filters['zlib'], complevel=filters['complevel'],
                                     shuffle=filters['shuffle'],
                                     chunksizes=None if chunking == 'contiguous' else chunking)
    else:
        dst_var = dst.createVariable(name, var.dtype, var.dimensions, fill_value=fill_value)
    dst_var.setncatts({att: var.getncattr(att) for att in var.ncattrs() if att != '_FillValue'})

# Copy global attributes and note what we did
dst.setncatts({att: src.getncattr(att) for att in src.ncattrs()})
dst.smoothing = f'Centered {smooth_window}-day running mean applied to {", ".join(vars_to_smooth)} (window shrinks at ends)'
dst.smoothing_source_file = str(input_path.resolve())
dst.smoothing_date = datetime.now().strftime('%Y-%m-%d %H:%M:%S')

# Copy everything we are NOT smoothing as-is (time, date, datesec, lat, lon, ...)
for name, var in src.variables.items():
    if name not in vars_to_smooth:
        print(f"Copying {name}...")
        dst[name][:] = var[:]

# ==================== SMOOTH ====================

nt = len(src.dimensions['time'])
print(f"Found {nt} times")

for var in vars_to_smooth:
    print(f"Smoothing {var}...")
    nlat = src[var].shape[1]
    nlon = src[var].shape[2]

    for start in range(0, nt, BLOCK_SIZE):
        end = min(start + BLOCK_SIZE, nt)   # this block writes times start to end-1

        # Read the block plus +/- half_window days on either side (bounded by the file)
        read_start = max(start - half_window, 0)
        read_end = min(end + half_window, nt)
        print(f"  times {start+1}-{end} of {nt} (reading {read_start+1}-{read_end})")
        data = src[var][read_start:read_end, :, :]  # Shape: (time, lat, lon)

        smoothed = np.zeros((end - start, nlat, nlon), dtype=data.dtype)
        for t in range(start, end):
            # Shrink the window near the ends so it stays centered on day t
            this_half = min(half_window, t, nt - 1 - t)
            # Indices in data are offset by read_start
            lo = t - this_half - read_start
            hi = t + this_half + 1 - read_start
            smoothed[t - start, :, :] = np.mean(data[lo:hi, :, :], axis=0)

        dst[var][start:end, :, :] = smoothed

src.close()
dst.close()

print(f"\n✓ Successfully created: {output_file}")
print(f"✓ Output file size: {Path(output_file).stat().st_size / (1024**3):.2f} GB")
