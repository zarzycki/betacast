import os
import sys
import xarray as xr
import numpy as np
import shutil
import argparse

# Betacast modules
module_paths = [
    ('functions_path', ['../..', 'py_functions']),
    ('functions_path', ['../..', 'py_remapping'])
]
for path_name, path_parts in module_paths:
    path = os.path.abspath(os.path.join(os.path.dirname(__file__), *path_parts))
    print(f"{path_name}: {path}")
    if path not in sys.path:
        sys.path.append(path)
from py_seedfuncs import keyword_values, gc_latlon
import pyfuncs
import meteo
import cam2cam
import horizremap
import vertremap
from ESMF_regridding import *

# Set up argument parser
parser = argparse.ArgumentParser(description='Add perturbations to CAM file')
parser.add_argument('--BEFOREPERTFILE', type=str, required=True,
                    help='Input CAM file before perturbation')
parser.add_argument('--AFTERPERTFILE', type=str, required=True,
                    help='Output CAM file after perturbation')
parser.add_argument('--gridfile', type=str, required=True,
                    help='Grid file for remapping')
parser.add_argument('--MAPFILEPATH', type=str, required=True,
                    help='Path for map files')
parser.add_argument('--pthi', type=str, required=True,
                    help='Path to perturbation namelist file')

args = parser.parse_args()

# Set variables from command line arguments
beforefile = args.BEFOREPERTFILE
afterfile = args.AFTERPERTFILE
gridfile = args.gridfile
mapfilepath = args.MAPFILEPATH
if mapfilepath != './':
    mapfilepath += '/'
pthi = args.pthi
print(f"mapfilepath: {mapfilepath}")
tmpfilepath = mapfilepath

print("Command line arguments:")
print(f"BEFOREPERTFILE: {args.BEFOREPERTFILE}")
print(f"AFTERPERTFILE: {args.AFTERPERTFILE}")
print(f"gridfile: {args.gridfile}")
print(f"MAPFILEPATH: {args.MAPFILEPATH}")
print(f"pthi: {args.pthi}")

# Get configuration values
warming_case = keyword_values(pthi, "case", "str")
basedir = keyword_values(pthi, "basedir", "str")
start_month = keyword_values(pthi, "start_month", "int")
end_month = keyword_values(pthi, "end_month", "int")
current_year = keyword_values(pthi, "current_year", "int")
comp_year = keyword_values(pthi, "comp_year", "int")

correct_sfc = keyword_values(pthi, "correct_sfc", "bool")
plevs = keyword_values(pthi, "plevs", "bool")
update_pressure = keyword_values(pthi, "update_pressure", "bool")
update_winds = keyword_values(pthi, "update_winds", "bool")
do_ps_corr = keyword_values(pthi, "do_ps_corr", "bool")
smooth_deltas = keyword_values(pthi, "smooth_deltas", "bool")

if smooth_deltas:
    smooth_delta_iter = keyword_values(pthi, "smooth_delta_iter", "int")

output_atm_diag = keyword_values(pthi, "output_atm_diag", "bool")
keep_esmf = keyword_values(pthi, "keep_esmf", "bool")
extra_diags_atm = keyword_values(pthi, "extra_diags_atm", "bool")

# Print configurations
print("************* Running ATM perturbation code *************")
print(f"Case: {warming_case}")
print(f"basedir: {basedir}")
print(f"start_month: {start_month}")
print(f"end_month: {end_month}")
print(f"current_year: {current_year}")
print(f"comp_year: {comp_year}")
print(f"BEFOREPERTFILE: {args.BEFOREPERTFILE}")
print(f"AFTERPERTFILE: {args.AFTERPERTFILE}")

print(f"plevs: {plevs}")
print(f"correct_sfc: {correct_sfc}")
print(f"update_pressure: {update_pressure}")
print(f"update_winds: {update_winds}")
print(f"do_ps_corr: {do_ps_corr}")

print(f"keep_esmf: {keep_esmf}")

print(f"smooth_deltas: {smooth_deltas}")
if smooth_deltas:
    print(f"smooth_delta_iter: {smooth_delta_iter}")

print(f"output_atm_diag: {output_atm_diag}")
print(f"extra_diags_atm: {extra_diags_atm}")

print("****************************************************")

# Set up mapfile path
mapfilepath = os.environ.get('MAPFILEPATH', './')
if mapfilepath != './':
    mapfilepath += '/'
print(f"mapfilepath: {mapfilepath}")
tmpfilepath = mapfilepath

# Copy before to after file
beforefile = args.BEFOREPERTFILE
afterfile = args.AFTERPERTFILE

# Open the after CAM file for modification
cam_ds = xr.open_dataset(beforefile)

# Extract variables from target mesh
T = cam_ds.T.values
Q = cam_ds.Q.values
U = cam_ds.U.values
V = cam_ds.V.values
lat = cam_ds.lat.values
lon = cam_ds.lon.values
camlev = cam_ds.lev.values
hyam = cam_ds.hyam.values
hybm = cam_ds.hybm.values
hyai = cam_ds.hyai.values
hybi = cam_ds.hybi.values
P0 = 100000.0  # Pa
PS = cam_ds.PS.values
ncol = len(lat)
nlev = len(hyam)

# Handle different warming_cases for delta files
if warming_case == "CAMC20C":  # Can only be used for attribution
    # Load delta files
    deltaFilePS = xr.open_dataset(f"{basedir}/{warming_case}_plev/delta_ps_CAM5-1-1degree_All-Hist_est1_v2.0_1996-2016.nc_Climatology_2016-1996.nc")
    deltaFileT = xr.open_dataset(f"{basedir}/{warming_case}_plev/delta_ta_CAM5-1-1degree_All-Hist_est1_v2.0_1996-2016.nc_Climatology_2016-1996.nc")
    deltaFileQ = xr.open_dataset(f"{basedir}/{warming_case}_plev/delta_hus_CAM5-1-1degree_All-Hist_est1_v2.0_1996-2016.nc_Climatology_2016-1996.nc")

    lat_in = deltaFileT.lat.values
    lon_in = deltaFileT.lon.values

    # Extract and process deltas
    deltaPS_in = deltaFilePS.delta_ps_Climatology_Monthly[start_month-1:end_month+1, 0, :, :]
    deltaT_in = deltaFileT.delta_ta_Climatology_Monthly[start_month-1:end_month+1, ::-1, :, :]
    deltaQ_in = deltaFileQ.delta_hus_Climatology_Monthly[start_month-1:end_month+1, ::-1, :, :]

    # Average over months (dim_avg_n_Wrap equivalent)
    deltaPS = deltaPS_in.mean(dim='time').values
    deltaT = deltaT_in.mean(dim='time').values
    deltaQ = deltaQ_in.mean(dim='time').values

elif warming_case == "CESMLENS":  # Can be used for attribution and future projections
    if plevs:
        deltaFilePS = xr.open_dataset(f"{basedir}/{warming_case}_plev/ens_PS_anom.nc")
        deltaFileT = xr.open_dataset(f"{basedir}/{warming_case}_plev/ens_T_anom.nc")
        deltaFileQ = xr.open_dataset(f"{basedir}/{warming_case}_plev/ens_Q_anom.nc")
        plev_in = deltaFileT.plev.values
    else:
        deltaFilePS = xr.open_dataset(f"{basedir}/{warming_case}_mlev/PS/ens_PS_anom.nc")
        deltaFileT = xr.open_dataset(f"{basedir}/{warming_case}_mlev/T/ens_T_anom.nc")
        deltaFileQ = xr.open_dataset(f"{basedir}/{warming_case}_mlev/Q/ens_Q_anom.nc")
        hyam_in = deltaFileT.hyam.values
        hybm_in = deltaFileT.hybm.values
        lev_in = deltaFileT.lev.values

    lat_in = deltaFileT.lat.values
    lon_in = deltaFileT.lon.values

    # Calculate indices and extract current year deltas
    current_start_idx = (current_year * 12 - 1920 * 12) + start_month - 1
    current_end_idx = (current_year * 12 - 1920 * 12) + end_month - 1

    deltaPS_in = deltaFilePS.PS[current_start_idx:current_end_idx+1]
    deltaT_in = deltaFileT.T[current_start_idx:current_end_idx+1]
    deltaQ_in = deltaFileQ.Q[current_start_idx:current_end_idx+1]

    # Average over months
    deltaPS_current = deltaPS_in.mean(dim='time')
    deltaT_current = deltaT_in.mean(dim='time')
    deltaQ_current = deltaQ_in.mean(dim='time')

    if comp_year < 1920:
        # Initialize with zeros while keeping metadata
        deltaPS_comp = xr.zeros_like(deltaPS_current)
        deltaT_comp = xr.zeros_like(deltaT_current)
        deltaQ_comp = xr.zeros_like(deltaQ_current)
    else:
        comp_start_idx = (comp_year * 12 - 1920 * 12) + start_month - 1
        comp_end_idx = (comp_year * 12 - 1920 * 12) + end_month - 1

        deltaPS_in = deltaFilePS.PS[comp_start_idx:comp_end_idx+1]
        deltaT_in = deltaFileT.T[comp_start_idx:comp_end_idx+1]
        deltaQ_in = deltaFileQ.Q[comp_start_idx:comp_end_idx+1]

        deltaPS_comp = deltaPS_in.mean(dim='time')
        deltaT_comp = deltaT_in.mean(dim='time')
        deltaQ_comp = deltaQ_in.mean(dim='time')

    # Calculate differences
    deltaPS = deltaPS_comp.values - deltaPS_current.values
    deltaT = deltaT_comp.values - deltaT_current.values
    deltaQ = deltaQ_comp.values - deltaQ_current.values
else:
    print("Unknown warming case...")
    print(f"----{warming_case}----")

# Fill missing values
deltaT = np.nan_to_num(deltaT, 0)
deltaQ = np.nan_to_num(deltaQ, 0)
deltaPS = np.nan_to_num(deltaPS, 0)

print("Generating weights!")
if os.path.isfile(gridfile):
    gridfilebasename = os.path.basename(gridfile)
    gridfilebasename = os.path.splitext(gridfilebasename)[0]
    print(f"Using existing SCRIP grid: {gridfilebasename}")
else:
    print("No gridfile available, generating one using unstructured_to_scrip")
    # gen SE grid
    Opt_se = {
        "ForceOverwrite": True,
        "PrintTimings": True,
        "Title": "SE Grid"
    }
    seGridName = f"grid_se_{ncol}.nc"
    if not keep_esmf or not os.path.isfile(os.path.join(mapfilepath, seGridName)):
        unstructured_to_scrip(os.path.join(mapfilepath, seGridName), lat, lon, Opt_se)

# gen RLL grid
Opt_ll = {
   "ForceOverwrite": True,
   "PrintTimings": True,
   "Title": "Deltas grid"
}
llGridName = "grid_deltas.nc"
if not keep_esmf or not os.path.isfile(os.path.join(mapfilepath, llGridName)):
    rectilinear_to_SCRIP(os.path.join(mapfilepath, llGridName), lat_in, lon_in, Opt_ll)

# Create weights
Opt = {
   "RemovePETLog": True,
   "InterpMethod": "patch",
   "ForceOverwrite": True,
   "PrintTimings": True,
   "DstESMF": False
}

# Forward mapping (source->ll)
if os.path.isfile(gridfile):
    Opt["SrcESMF"] = False
    src_path = gridfile
    wgtFileName1 = f"map_{gridfilebasename}_to_ll.nc"
else:
    Opt["SrcESMF"] = True
    src_path = os.path.join(mapfilepath, seGridName)
    wgtFileName1 = f"map_se_{ncol}_to_ll.nc"

if not keep_esmf or not os.path.isfile(os.path.join(mapfilepath, wgtFileName1)):
    esmf_regrid_gen_weights(
                          src_path,
                          os.path.join(mapfilepath, llGridName),
                          os.path.join(mapfilepath, wgtFileName1),
                          Opt)

# Backward mapping (ll->source)
Opt["SrcESMF"] = False
if os.path.isfile(gridfile):
    Opt["DstESMF"] = False
    dst_path = gridfile
    wgtFileName2 = f"map_ll_to_{gridfilebasename}.nc"
else:
    Opt["DstESMF"] = True
    dst_path = os.path.join(mapfilepath, seGridName)
    wgtFileName2 = f"map_ll_to_se_{ncol}.nc"

if not keep_esmf or not os.path.isfile(os.path.join(mapfilepath, wgtFileName2)):
    esmf_regrid_gen_weights(
                          os.path.join(mapfilepath, llGridName),
                          dst_path,
                          os.path.join(mapfilepath, wgtFileName2),
                          Opt)

print("Beginning PS to lat-lon interp...")
print(f"PS.shape: {PS.shape}")
# NOTE, PS[0] squeezes off the singleton dimension "time" -- it works for ND arrays (so both struct and unstruct)
# versions PS[0,:], PS[0,:,:], etc.
PS_deltaGrid, _, _ = horizremap.remap_with_weights_wrapper(
                                        PS[0],
                                        os.path.join(mapfilepath, wgtFileName1))

# Perform vertical interpolation
print("Beginning vertical interpolation...")
if warming_case == "CAMC20C":
    #CMZ interface def pressure_to_hybrid(p_levels, data_on_p_levels, ps, a_coeff, b_coeff, level_dim=0, p0=100000, kflag=1):
    deltaTCAM = vertremap.pressure_to_hybrid(plev_in, deltaT, PS_deltaGrid, hyam, hybm)
    deltaQCAM = vertremap.pressure_to_hybrid(plev_in, deltaQ, PS_deltaGrid, hyam, hybm)
elif warming_case == "CESMLENS":
    if plevs:
        deltaTCAM = vertremap.pressure_to_hybrid(plev_in, deltaT, PS_deltaGrid, hyam, hybm)
        deltaQCAM = vertremap.pressure_to_hybrid(plev_in, deltaQ, PS_deltaGrid, hyam, hybm)
    else:
        deltaTCAM = vertremap.hyi2hyo(P0, hyam_in, hybm_in, PS_deltaGrid, deltaT, hyam, hybm, intflg=1, unstructured=False)
        deltaQCAM = vertremap.hyi2hyo(P0, hyam_in, hybm_in, PS_deltaGrid, deltaQ, hyam, hybm, intflg=1, unstructured=False)

deltaPSCAM = deltaPS  # no vertical interpolation needed for surface pressure

# Apply smoothing if requested
if smooth_deltas:
    print(f"Applying smoothing with {smooth_delta_iter} iterations...")
    deltaQCAM = pyfuncs.smooth_with_smth9(deltaQCAM, smooth_delta_iter, p=0.5, q=0.25)
    deltaTCAM = pyfuncs.smooth_with_smth9(deltaTCAM, smooth_delta_iter, p=0.5, q=0.25)
    deltaPSCAM = pyfuncs.smooth_with_smth9(deltaPSCAM, smooth_delta_iter, p=0.5, q=0.25)

# Calculate extra diagnostics if requested
if extra_diags_atm:
    dp = cam2cam.dpres_hybrid_ccm(deltaPSCAM, p0, hyai, hybi)
    pw = meteo.prcwater_dp(deltaQCAM, dp)

# Interpolate deltas back to original grid
print("Interpolating deltas back to original grid...")
deltaTCAM_interp, _, _ = horizremap.remap_with_weights_wrapper(
                                deltaTCAM,
                                os.path.join(mapfilepath, wgtFileName2))
deltaQCAM_interp, _, _ = horizremap.remap_with_weights_wrapper(
                                deltaQCAM,
                                os.path.join(mapfilepath, wgtFileName2))
deltaPSCAM_interp, _, _ = horizremap.remap_with_weights_wrapper(
                                deltaPSCAM,
                                os.path.join(mapfilepath, wgtFileName2))

# Fill any missing values
deltaTCAM_interp = np.nan_to_num(deltaTCAM_interp, 0)
deltaQCAM_interp = np.nan_to_num(deltaQCAM_interp, 0)
deltaPSCAM_interp = np.nan_to_num(deltaPSCAM_interp, 0)

if do_ps_corr:

    # Perform empirical correction to PS with emphasis over low PS areas
    anom_scaling = 3.0  # vertical average reference Tanom for scaling
    print(f"Doing empirical ps_corr with anom_scaling: {anom_scaling}")

    print("Calculate pint and dpint")
    # Do weighted integral of deltaTCAM_interp
    dpint = np.zeros_like(deltaTCAM_interp)
    pint = np.zeros((nlev+1, ncol), dtype=hyai.dtype)

    # Calculate interface pressures
    pint = np.expand_dims(hyai, 1) * P0 + np.expand_dims(hybi, 1) * np.expand_dims(PS[0], 0)
    dpint = pint[1:nlev+1, :] - pint[0:nlev, :]

    print("Calculate weighted T anomaly and sign")
    # Calculate vertical weighted average
    Tanom = np.sum(deltaTCAM_interp * dpint, axis=0) / np.sum(dpint, axis=0)
    # Find where column is warmer and where is colder
    Tsign = np.where(Tanom >= 0.0, 1.0, -1.0)

    print("Remap to RLL")
    Tsign_deltaGrid, _, _ = horizremap.remap_with_weights_wrapper(
                                    Tsign,
                                    os.path.join(mapfilepath, wgtFileName1))

    print("Smooth")
    smoothiter = 50
    Tsign_deltaGrid = pyfuncs.smooth_with_smth9(Tsign_deltaGrid, smoothiter, p=0.5, q=0.25)

    print(f"Missing values: {np.sum(np.isnan(Tsign_deltaGrid))}")

    print("Remap to SE")
    Tsign, _, _ = horizremap.remap_with_weights_wrapper(
                                    Tsign_deltaGrid,
                                    os.path.join(mapfilepath, wgtFileName2))

    print(f"before {np.sum(np.isnan(Tsign))}")

    # Handle missing values
    if np.any(np.isnan(Tsign)):
        print("Found some missing Tsign data")
        print(f"Number of missing values: {np.sum(np.isnan(Tsign))}")
        Tsignorig = Tsign.copy()

        for ii in range(ncol):
            if np.isnan(Tsign[ii]):
                thisLat = lat[ii]
                thisLon = lon[ii]
                # Need to implement gc_latlon equivalent
                gcdist = gc_latlon(thisLat, thisLon, lat, lon, 2, 2)
                gcdist = np.where(np.isnan(Tsignorig), 999999., gcdist)
                Tsign[ii] = Tsign[np.argmin(gcdist)]
                print(f"Replacing Tsign at {ii} with {np.argmin(gcdist)}")
                print(f"this lat/lon {lat[ii]} {lon[ii]} nearest: {lat[np.argmin(gcdist)]} {lon[np.argmin(gcdist)]}")

    print(f"after {np.sum(np.isnan(Tsign))}")

    # Scale the multiplier based on the magnitude of the column anomaly
    Tsign = Tsign * (np.abs(Tanom)/anom_scaling)

    print(f"Missing values: {np.sum(np.isnan(Tsign))}")

    # Correction coefficients derived from ne30 run
    rc = -0.007590353
    rc_intercept = 771.9941
    print(f"CORR: using dPSL = {rc}*PS + {rc_intercept}")

    # Generate PS corr based on PS and constants
    PScorr = PS[0,:] * rc + rc_intercept

    print(f"Missing values: {np.sum(np.isnan(PScorr))}")

    # Scale correction by Tanom + heating/cooling
    PScorr = PScorr * Tsign

    print(f"Missing values: {np.sum(np.isnan(PScorr))}")

    PS[0,:] = PS[0,:] + PScorr

# Update arrays
print(f"T before adding deltas - min: {T.min():.2f}, max: {T.max():.2f}, mean: {T.mean():.2f}")
T[0,:,:] = T[0,:,:] + deltaTCAM_interp
print(f"T after adding deltas - min: {T.min():.2f}, max: {T.max():.2f}, mean: {T.mean():.2f}")
print(f"Q before adding deltas - min: {Q.min():.6f}, max: {Q.max():.6f}, mean: {Q.mean():.6f}")
Q[0,:,:] = Q[0,:,:] + deltaQCAM_interp
print(f"Q after adding deltas - min: {Q.min():.6f}, max: {Q.max():.6f}, mean: {Q.mean():.6f}")
if update_pressure:
    print(f"PS before adding deltas - min: {PS.min():.2f}, max: {PS.max():.2f}, mean: {PS.mean():.2f}")
    PS[0,:] = PS[0,:] + deltaPSCAM_interp
    print(f"PS after adding deltas - min: {PS.min():.2f}, max: {PS.max():.2f}, mean: {PS.mean():.2f}")
if update_winds:
    print(f"U before adding deltas - min: {U.min():.2f}, max: {U.max():.2f}, mean: {U.mean():.2f}")
    U[0,:,:] = U[0,:,:] + deltaUCAM_interp
    print(f"U after adding deltas - min: {U.min():.2f}, max: {U.max():.2f}, mean: {U.mean():.2f}")
    print(f"V before adding deltas - min: {V.min():.2f}, max: {V.max():.2f}, mean: {V.mean():.2f}")
    V[0,:,:] = V[0,:,:] + deltaVCAM_interp
    print(f"V after adding deltas - min: {V.min():.2f}, max: {V.max():.2f}, mean: {V.mean():.2f}")

if correct_sfc:
    for ii in range(nlev):
        T[0,ii,:] = T[0,ii,:] - deltaTCAM_interp[nlev-1,:]
        Q[0,ii,:] = Q[0,ii,:] - deltaQCAM_interp[nlev-1,:]

# Where moisture is negative due to deltas, reset
print(f"Reset {np.sum(Q <= 0)} Q elements for being negative")
Q = np.where(Q <= 0, 1.0e-9, Q)

# Check for missing values
for var_name, var in [('Q', Q), ('T', T), ('PS', PS)]:
    if np.any(np.isnan(var)):
        print(f"{var_name} data contains some missing values. Beware.")
        print(f"Number of missing values: {np.sum(np.isnan(var))}")

# Write variables back to file
print("Writing variables back to file...")
print("Writing T and Q")
cam_ds['T'].values[...] = T
cam_ds['Q'].values[...] = Q
if update_pressure or do_ps_corr:
    print("Writing PS")
    cam_ds['PS'].values[...] = PS
if update_winds:
    print("Writing U and V")
    cam_ds['U'].values[...] = U
    cam_ds['V'].values[...] = V

cam_ds.to_netcdf(afterfile, mode='w')
cam_ds.close()

if output_atm_diag:
    diag_filename = os.path.join(tmpfilepath, "deltas_atm.nc")
    print(f"outputting diags to {diag_filename}")

    if os.path.exists(diag_filename):
        os.remove(diag_filename)

    # Create diagnostic file with explicit dimensions
    diag_ds = xr.Dataset(
        data_vars={
            'PS_deltaGrid': (('lat', 'lon'), PS_deltaGrid),
            'deltaPSCAM': (('lat', 'lon'), deltaPSCAM),
            'deltaTCAM': (('lev', 'lat', 'lon'), deltaTCAM),
            'deltaQCAM': (('lev', 'lat', 'lon'), deltaQCAM)
        },
        coords={
            'lat': lat_in,
            'lon': lon_in,
            'lev': camlev
        },
        attrs={
            'creation_date': datetime.now().strftime("%a %b %d %H:%M:%S %Y")
        }
    )

    if extra_diags_atm:
        diag_ds['pw'] = pw

    if do_ps_corr:
        Tanom_remap, _, _ = horizremap.remap_with_weights_wrapper(
                                    Tanom,
                                    os.path.join(mapfilepath, wgtFileName1))
        Tsign_remap, _, _ = horizremap.remap_with_weights_wrapper(
                                    Tsign,
                                    os.path.join(mapfilepath, wgtFileName1))
        PScorr_remap, _, _ = horizremap.remap_with_weights_wrapper(
                                    PScorr,
                                    os.path.join(mapfilepath, wgtFileName1))
        diag_ds['Tanom'] = (('lat', 'lon'), Tanom_remap)
        diag_ds['Tsign'] = (('lat', 'lon'), Tsign_remap)
        diag_ds['PScorr'] = (('lat', 'lon'), PScorr_remap)

    diag_ds.to_netcdf(diag_filename)
