#!/usr/bin/env python3
# Python port of add_perturbations_to_datm.ncl
#
# Adds CESM-LENS ensemble-mean anomalies (T, Q, and optionally FLDS/FSDS) to existing DATM
# forcing files in place. T/Q are additive, FLDS/FSDS/PRECT are multiplicative (PRECT is
# scaled by 1 + 0.03*dT, i.e., ~3%/K).
#
# python add_perturbations_to_datm.py --datm_file_name TPQWL.nc --datm2_file_name Prec.nc \
#     --datm3_file_name Solar.nc --current_year 1996 --comp_year 2099

import numpy as np
import xarray as xr
import netCDF4
import cftime
import argparse
from datetime import datetime
from scipy.interpolate import interp1d, RegularGridInterpolator

# Parse command line arguments
parser = argparse.ArgumentParser(description="Add LENS anomalies to DATM forcing files")
parser.add_argument("--datm_file_name", required=True, help="DATM TPQWL file (TBOT, QBOT, FLDS) to modify in place")
parser.add_argument("--datm2_file_name", required=True, help="DATM precip file (PRECTmms) to modify in place")
parser.add_argument("--datm3_file_name", default=None, help="DATM solar file (FSDS) to modify in place (required with --do_sw_fluxes)")
parser.add_argument("--current_year", type=int, default=1996, help="Current year (valid 1921-onwards)")
parser.add_argument("--comp_year", type=int, default=2099, help="Comparison year (< start_year for no perturbation)")
parser.add_argument("--delta_basedir", default="/glade/u/home/zarzycki/scratch/CESM_LENS_temp",
                    help="Base dir containing {T,Q,FSDS,FLDS}/ens_{VAR}_anom.nc")
parser.add_argument("--do_sw_fluxes", action="store_true", help="Perturb shortwave (FSDS)")
parser.add_argument("--no_lw_fluxes", dest="do_lw_fluxes", action="store_false", help="Do not perturb longwave (FLDS)")
parser.add_argument("--output_diag", action="store_true", help="Write first 5 times of orig/perturbed fields to datm.nc")
parser.add_argument("--lev_index", type=int, default=29, help="Level index of T/Q anomalies to use (default: 29, lowest level)")
parser.add_argument("--start_year", type=int, default=1920, help="First year of anomaly files (default: 1920)")
parser.add_argument("--start_month", type=int, default=1, help="Start month (default: 1)")
parser.add_argument("--end_month", type=int, default=12, help="End month (default: 12)")
parser.add_argument("--shouldermon", type=int, default=1, help="Shoulder months before start_month (default: 1)")
parser.set_defaults(do_lw_fluxes=True)

args = parser.parse_args()

datm_file_name = args.datm_file_name
datm2_file_name = args.datm2_file_name
datm3_file_name = args.datm3_file_name
current_year = args.current_year
comp_year = args.comp_year
delta_basedir = args.delta_basedir
do_sw_fluxes = args.do_sw_fluxes
do_lw_fluxes = args.do_lw_fluxes
output_diag = args.output_diag
lev_index = args.lev_index
start_year = args.start_year
start_month = args.start_month
end_month = args.end_month
shouldermon = args.shouldermon

if do_sw_fluxes and datm3_file_name is None:
    parser.error("--datm3_file_name is required with --do_sw_fluxes")

# The anomalies actually applied come from comp_year (or current_year if comp_year predates
# the anomaly files, which is how a "no perturbation" run is configured)
if comp_year < start_year:
    anom_year = current_year
else:
    anom_year = comp_year
stix = anom_year*12 - start_year*12 + start_month-1 - shouldermon
enix = anom_year*12 - start_year*12 + end_month-1

# Variables to perturb: T/Q are 3D (time,lev,lat,lon), fluxes are 2D (time,lat,lon)
delta_vars = ["T", "Q"]
if do_sw_fluxes:
    delta_vars.append("FSDS")
if do_lw_fluxes:
    delta_vars.append("FLDS")

print("Reading deltas...")
deltas = {}
for var in delta_vars:
    f = xr.open_dataset(f"{delta_basedir}/{var}/ens_{var}_anom.nc")
    delta = f[var].isel(time=slice(stix, enix+1)).astype("float64").load()
    if delta.ndim == 4:
        delta = delta.isel({delta.dims[1]: lev_index})
    deltas[var] = delta
    f.close()

delta_lat = deltas["T"]["lat"].values
delta_lon = deltas["T"]["lon"].values

# Day of year for the anomalies. Last anomaly is the end of the year, which keeps the
# time axis monotonic for the interpolation below.
time_deltas_doy = np.array([d.dayofyr for d in deltas["T"]["time"].values], dtype="float64")
time_deltas_doy[-1] = 366
print(time_deltas_doy)

print(datm_file_name)
print(datm2_file_name)
print(datm3_file_name)
datm_file = netCDF4.Dataset(datm_file_name, "r+")
datm2_file = netCDF4.Dataset(datm2_file_name, "r+")
if do_sw_fluxes:
    datm3_file = netCDF4.Dataset(datm3_file_name, "r+")

datm_time = datm_file.variables["time"]
datm_dates = cftime.num2date(datm_time[:], datm_time.units, calendar=getattr(datm_time, "calendar", "noleap"))
datm_deltas_doy = np.array([d.dayofyr + d.hour/24. for d in datm_dates])

print("Time interpolation")
# interp LENS to time deltas in DATM
interp = {}
for var in delta_vars:
    interp[var] = interp1d(time_deltas_doy, deltas[var].values, axis=0, bounds_error=False)(datm_deltas_doy)
interp["PRECT"] = 1.0 + (0.03 * interp["T"])

print("Spatial interpolation output")
# now interp the time interpolated data to the DATM grid
datm_lat = datm_file.variables["LATIXY"][:, 0]
datm_lon = datm_file.variables["LONGXY"][0, :]

# Pad one column on each side of the source grid so the interpolation is cyclic in lon
lon_ext = np.concatenate([[delta_lon[-1]-360.], delta_lon, [delta_lon[0]+360.]])
query_lon = delta_lon[0] + np.mod(datm_lon - delta_lon[0], 360.)
query_lat2d, query_lon2d = np.meshgrid(datm_lat, query_lon, indexing="ij")
query_pts = np.column_stack([query_lat2d.ravel(), query_lon2d.ravel()])

interp_DATM = {}
for var in interp:
    data_ext = np.concatenate([interp[var][:, :, -1:], interp[var], interp[var][:, :, :1]], axis=2)
    interpolator = RegularGridInterpolator((delta_lat, lon_ext), np.moveaxis(data_ext, 0, -1),
                                           method="linear", bounds_error=False, fill_value=np.nan)
    out = interpolator(query_pts).reshape(len(datm_lat), len(datm_lon), -1)
    interp_DATM[var] = np.moveaxis(out, -1, 0)

print("Read and update DATM data")
# (file, DATM variable, interp variable, operator)
updates = [(datm_file, "TBOT", "T", "add"),
           (datm_file, "QBOT", "Q", "add"),
           (datm2_file, "PRECTmms", "PRECT", "mult")]
if do_sw_fluxes:
    updates.append((datm3_file, "FSDS", "FSDS", "mult"))
if do_lw_fluxes:
    updates.append((datm_file, "FLDS", "FLDS", "mult"))

orig = {}
new = {}
for ncfile, datm_var, interp_var, operator in updates:
    orig[datm_var] = ncfile.variables[datm_var][:]
    if operator == "add":
        new[datm_var] = orig[datm_var] + interp_DATM[interp_var].astype(orig[datm_var].dtype)
    else:
        new[datm_var] = orig[datm_var] * interp_DATM[interp_var].astype(orig[datm_var].dtype)

# write new fields to perturbation file
print("Writing output")
for ncfile, datm_var, interp_var, operator in updates:
    ncfile.variables[datm_var][:] = new[datm_var]

datm_file.close()
datm2_file.close()
if do_sw_fluxes:
    datm3_file.close()

if output_diag:
    diag_filename = "datm.nc"
    print(f"outputting diags to {diag_filename}")
    diag_vars = {}
    for datm_var in new:
        diag_vars[f"{datm_var}_orig"] = (("time", "lat", "lon"), np.ma.filled(orig[datm_var][0:5], np.nan))
        diag_vars[datm_var] = (("time", "lat", "lon"), np.ma.filled(new[datm_var][0:5], np.nan))
    for interp_var in interp_DATM:
        diag_vars[f"interp_{interp_var}_DATM"] = (("time", "lat", "lon"), interp_DATM[interp_var][0:5])
    diag = xr.Dataset(diag_vars, attrs={"creation_date": datetime.now().strftime("%c")})
    diag.to_netcdf(diag_filename, unlimited_dims=["time"])
