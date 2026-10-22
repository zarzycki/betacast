import xarray as xr
import numpy as np
import os
import subprocess
import datetime
import argparse

def parse_args():
    parser = argparse.ArgumentParser(description='Normalize DATM deltas for climate simulations')
    # Required arguments
    parser.add_argument('--current_year', type=int, required=True,
                        help='Reference year for normalization')
    parser.add_argument('--basedir', type=str, required=True,
                        help='Base directory containing anomaly files')

    # Optional configuration parameters with defaults
    parser.add_argument('--start_year', type=int, default=1920,
                        help='Start year for calculations (default: 1920)')
    parser.add_argument('--start_month', type=int, default=1,
                        help='Start month (default: 1)')
    parser.add_argument('--end_month', type=int, default=12,
                        help='End month (default: 12)')
    parser.add_argument('--shouldermon', type=int, default=1,
                        help='Shoulder month (default: 1)')
    parser.add_argument('--debug', action='store_true',
                        help='Enable debug mode (default: False)')
    parser.add_argument('--no-compress', dest='do_compress', action='store_false',
                        help='Disable compression of output files (default: compression enabled)')
    parser.set_defaults(do_compress=True)

    return parser.parse_args()

# Parse command-line arguments
args = parse_args()
basedir = args.basedir
current_year = args.current_year

# Configuration parameters from command-line arguments
# These rarely get updated and were "hardcoded" in the old NCL version, but here
# it makes sense to give the user more control (see defaults in the function above)
start_year = args.start_year
start_month = args.start_month
end_month = args.end_month
shouldermon = args.shouldermon
debug = args.debug
do_compress = args.do_compress

print(f"basedir: {basedir}")
print(f"current_year: {current_year}")
print(f"Configuration: start_year={start_year}, start_month={start_month}, end_month={end_month}, shouldermon={shouldermon}")
print(f"Debug mode: {debug}, Compression: {do_compress}")

# These are the vars we need to loop over and how to process them
# ens_vars are what are from LENS
# internal_vars are the "target" for the output ANOM files
# operators means are the anomalies a multiplier (e.g., 1.05) or a straight diff (-1).
ens_vars = ["FLDS", "PRECT", "QBOT", "TBOT"]
internal_vars = ["FLDS", "PRECT", "Q", "T"]
operators = ["ratio", "ratio", "diff", "diff"]

for ii in range(len(ens_vars)):
    newfile = f"{basedir}/ens_{ens_vars[ii]}_{current_year}ref_anom.nc"
    # Check if newfile exists
    if not os.path.exists(newfile):
        # Print diagnostics
        print(f"ANOM_DIAG: {current_year} {ens_vars[ii]} {internal_vars[ii]} {operators[ii]}")
        # Load data
        # decode_times=False so time/time_bnds are written back exactly as read (matches NCL)
        f = xr.open_dataset(f"{basedir}/ens_{ens_vars[ii]}_anom.nc", decode_times=False)
        # Load the reference year
        refyr_stix = current_year*12 - start_year*12 + start_month-1 - shouldermon
        refyr_enix = current_year*12 - start_year*12 + end_month-1
        print(f"ANOM_DIAG: refyr_stix: {refyr_stix}")
        print(f"ANOM_DIAG: refyr_enix: {refyr_enix}")
        # Select reference year data
        ref_year = f[internal_vars[ii]].isel(time=slice(refyr_stix, refyr_enix+1))
        # Average across time for the reference year
        # Accumulate in float64 like NCL dim_avg_n, then cast back to the input type
        ref_year_avg = ref_year.mean(dim="time", skipna=True, dtype="float64").astype(ref_year.dtype)
        # Load the existing anom for the whole datastream
        new_anom = f[internal_vars[ii]].copy(deep=True)
        # Update the anom based on the ref year avg
        if operators[ii] == "diff":
            # Equivalent to NCL subtract with conform
            new_anom = (new_anom - ref_year_avg).assign_attrs(f[internal_vars[ii]].attrs)
        else:
            # Equivalent to NCL division with conform
            new_anom = (new_anom / ref_year_avg).assign_attrs(f[internal_vars[ii]].attrs)
        if debug:
            debug_file = f"simple_{ens_vars[ii]}.nc"
            if os.path.exists(debug_file):
                os.remove(debug_file)
            # Create debug NetCDF file
            debug_ds = xr.Dataset({
                "new_anom": new_anom,
                "ref_year_avg": ref_year_avg,
                internal_vars[ii]: f[internal_vars[ii]]
            })
            debug_ds.attrs["creation_date"] = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
            debug_ds.to_netcdf(debug_file)
        print(f"Writing {newfile}")
        # Remove existing file if it exists
        if os.path.exists(newfile):
            os.remove(newfile)
        # Create new output dataset
        out_ds = xr.Dataset({
            internal_vars[ii]: new_anom
        })
        # Copy time_bnds from original file
        out_ds["time_bnds"] = f["time_bnds"]
        # Set creation date attribute
        out_ds.attrs["creation_date"] = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        # No NaN _FillValue on coords/time_bnds; keep the input fill/missing value on the data var
        enc = {k: {"_FillValue": None} for k in out_ds.variables}
        enc[internal_vars[ii]] = {k: f[internal_vars[ii]].encoding[k] for k in ("_FillValue", "missing_value") if k in f[internal_vars[ii]].encoding}
        # Write to NetCDF file with unlimited time dimension
        out_ds.to_netcdf(newfile, unlimited_dims=["time"], encoding=enc)
        if do_compress:
            print("... compressing")
            systemcall = f"ncks -4 -L 1 -O {newfile} {newfile}"
            print(f" {systemcall}")
            subprocess.run(systemcall, shell=True, check=True)
            print("... success!")
    else:
        print(f"{newfile} already exists!")