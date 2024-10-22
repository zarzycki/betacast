import os
import sys
import xarray as xr
import numpy as np
from pathlib import Path
import shutil
from scipy.interpolate import griddata

# Betacast modules
module_path = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', 'tcseed'))
print(module_path)
if module_path not in sys.path:
    sys.path.append(module_path)
from py_seedfuncs import keyword_values
import meteo
import cam2cam

# Get config file path from environment or use default
pthi = "../../namelists/perturb.sample.nl"

# Get configuration values
case = keyword_values(pthi, "case", "str")
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
esmf_remap = keyword_values(pthi, "esmf_remap", "bool")
smooth_deltas = keyword_values(pthi, "smooth_deltas", "bool")

if smooth_deltas:
    smooth_delta_iter = keyword_values(pthi, "smooth_delta_iter", "int")

output_atm_diag = keyword_values(pthi, "output_atm_diag", "bool")
keep_esmf = keyword_values(pthi, "keep_esmf", "bool")
extra_diags_atm = keyword_values(pthi, "extra_diags_atm", "bool")

# Print configurations
print("************* Running ATM perturbation code *************")
print(f"Case: {case}")
print(f"basedir: {basedir}")
print(f"start_month: {start_month}")
print(f"end_month: {end_month}")
print(f"current_year: {current_year}")
print(f"comp_year: {comp_year}")
print(f"BEFOREPERTFILE: {os.environ.get('BEFOREPERTFILE')}")
print(f"AFTERPERTFILE: {os.environ.get('AFTERPERTFILE')}")

print(f"plevs: {plevs}")
print(f"correct_sfc: {correct_sfc}")
print(f"update_pressure: {update_pressure}")
print(f"update_winds: {update_winds}")
print(f"do_ps_corr: {do_ps_corr}")

print(f"esmf_remap: {esmf_remap}")
print(f"keep_esmf: {keep_esmf}")

print(f"smooth_deltas: {smooth_deltas}")
if smooth_deltas:
    print(f"smooth_delta_iter: {smooth_delta_iter}")

print(f"output_atm_diag: {output_atm_diag}")
print(f"extra_diags_atm: {extra_diags_atm}")

print("****************************************************")

# Set up mapfile path
mapfilepath = os.environ.get('MAPFILEPATH', './')
print(f"mapfilepath: {mapfilepath}")

# Copy before to after file
beforefile = os.environ.get('BEFOREPERTFILE')
afterfile = os.environ.get('AFTERPERTFILE')
shutil.copy2(beforefile, afterfile)

# Open the CAM file for modification
cam_ds = xr.open_dataset(afterfile, mode='a')

# Extract variables
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

# Handle different cases for delta files
if case == "CAMC20C":
    # Load delta files
    delta_ps = xr.open_dataset(f"{basedir}/{case}_plev/delta_ps_CAM5-1-1degree_All-Hist_est1_v2.0_1996-2016.nc_Climatology_2016-1996.nc")
    delta_t = xr.open_dataset(f"{basedir}/{case}_plev/delta_ta_CAM5-1-1degree_All-Hist_est1_v2.0_1996-2016.nc_Climatology_2016-1996.nc")
    delta_q = xr.open_dataset(f"{basedir}/{case}_plev/delta_hus_CAM5-1-1degree_All-Hist_est1_v2.0_1996-2016.nc_Climatology_2016-1996.nc")

    # Extract and process deltas
    deltaPS_in = delta_ps.delta_ps_Climatology_Monthly[start_month-1:end_month, 0, :, :].values
    deltaT_in = delta_t.delta_ta_Climatology_Monthly[start_month-1:end_month, ::-1, :, :].values
    deltaQ_in = delta_q.delta_hus_Climatology_Monthly[start_month-1:end_month, ::-1, :, :].values

    # Average over months
    deltaPS = np.mean(deltaPS_in, axis=0)
    deltaT = np.mean(deltaT_in, axis=0)
    deltaQ = np.mean(deltaQ_in, axis=0)

elif case == "CESMLENS":
    plev_suffix = "_plev" if plevs else "_mlev"

    # Load delta files
    if plevs:
        delta_ps = xr.open_dataset(f"{basedir}/{case}{plev_suffix}/ens_PS_anom.nc")
        delta_t = xr.open_dataset(f"{basedir}/{case}{plev_suffix}/ens_T_anom.nc")
        delta_q = xr.open_dataset(f"{basedir}/{case}{plev_suffix}/ens_Q_anom.nc")
    else:
        delta_ps = xr.open_dataset(f"{basedir}/{case}{plev_suffix}/PS/ens_PS_anom.nc")
        delta_t = xr.open_dataset(f"{basedir}/{case}{plev_suffix}/T/ens_T_anom.nc")
        delta_q = xr.open_dataset(f"{basedir}/{case}{plev_suffix}/Q/ens_Q_anom.nc")
        hyam_in = delta_t.hyam.values
        hybm_in = delta_t.hybm.values

    # Calculate indices for time selection
    current_start_idx = (current_year * 12 - 1920 * 12) + start_month - 1
    current_end_idx = (current_year * 12 - 1920 * 12) + end_month - 1

    # Extract current year deltas
    deltaPS_in = delta_ps.PS[current_start_idx:current_end_idx].values.astype(np.float64)
    deltaT_in = delta_t.T[current_start_idx:current_end_idx].values.astype(np.float64)
    deltaQ_in = delta_q.Q[current_start_idx:current_end_idx].values.astype(np.float64)

    deltaPS_current = np.mean(deltaPS_in, axis=0)
    deltaT_current = np.mean(deltaT_in, axis=0)
    deltaQ_current = np.mean(deltaQ_in, axis=0)

    if comp_year < 1920:
        deltaPS = np.zeros_like(deltaPS_current)
        deltaT = np.zeros_like(deltaT_current)
        deltaQ = np.zeros_like(deltaQ_current)
    else:
        comp_start_idx = (comp_year * 12 - 1920 * 12) + start_month - 1
        comp_end_idx = (comp_year * 12 - 1920 * 12) + end_month - 1

        deltaPS_in = delta_ps.PS[comp_start_idx:comp_end_idx].values.astype(np.float64)
        deltaT_in = delta_t.T[comp_start_idx:comp_end_idx].values.astype(np.float64)
        deltaQ_in = delta_q.Q[comp_start_idx:comp_end_idx].values.astype(np.float64)

        deltaPS_comp = np.mean(deltaPS_in, axis=0)
        deltaT_comp = np.mean(deltaT_in, axis=0)
        deltaQ_comp = np.mean(deltaQ_in, axis=0)

        deltaPS = deltaPS_comp - deltaPS_current
        deltaT = deltaT_comp - deltaT_current
        deltaQ = deltaQ_comp - deltaQ_current

# Fill missing values
deltaT = np.nan_to_num(deltaT, 0)
deltaQ = np.nan_to_num(deltaQ, 0)
deltaPS = np.nan_to_num(deltaPS, 0)

# ESMF remapping section
if esmf_remap:
    import ESMF

    def create_esmf_grid(lats, lons, is_sphere=True):
        grid = ESMF.Grid(np.array([len(lats), len(lons)]),
                        staggering=[False, False],
                        coord_sys=ESMF.CoordSys.SPH_DEG if is_sphere else ESMF.CoordSys.CART,
                        num_peri_dims=0, periodic_dim=None)

        [lon_grid, lat_grid] = [0, 1] if is_sphere else [1, 0]
        lon_ptr = grid.get_coords(lon_grid)
        lat_ptr = grid.get_coords(lat_grid)

        lon_ptr[...] = lons.reshape((len(lons), 1))
        lat_ptr[...] = lats.reshape((1, len(lats)))

        return grid

    # Create source and destination grids
    if 'gridfile' in locals():
        grid_src = xr.open_dataset(gridfile)
        src_grid = create_esmf_grid(grid_src.lat.values, grid_src.lon.values)
    else:
        src_grid = create_esmf_grid(lat, lon)

    dst_grid = create_esmf_grid(deltaT.lat.values, deltaT.lon.values)

    # Create ESMF regridder
    regridder = ESMF.Regridder(src_grid, dst_grid, "bilinear")

    # Perform regridding
    PS_deltaGrid = regridder(PS[0])
else:
    # Use scipy's griddata for interpolation
    points = np.column_stack((lat.flatten(), lon.flatten()))
    xi = np.column_stack((deltaT.lat.values.flatten(), deltaT.lon.values.flatten()))
    PS_deltaGrid = griddata(points, PS[0].flatten(), xi, method='linear')

PS_deltaGrid = xr.DataArray(PS_deltaGrid, dims=['lat', 'lon'],
                           coords={'lat': deltaT.lat, 'lon': deltaT.lon})
PS_deltaGrid.attrs['units'] = 'Pa'

# Perform vertical interpolation
print("Beginning vertical interpolation...")
if case == "CAMC20C":
    #CMZ interface def pressure_to_hybrid(p_levels, data_on_p_levels, ps, a_coeff, b_coeff, level_dim=0, p0=100000, kflag=1):
    deltaTCAM = vertremap.pressure_to_hybrid(deltaT.plev.values, deltaT, PS_deltaGrid, hyam, hybm)
    deltaQCAM = vertremap.pressure_to_hybrid(deltaQ.plev.values, deltaQ, PS_deltaGrid, hyam, hybm)
elif case == "CESMLENS":
    if plevs:
        deltaTCAM = vertremap.pressure_to_hybrid(deltaT.plev.values, deltaT, PS_deltaGrid, hyam, hybm)
        deltaQCAM = vertremap.pressure_to_hybrid(deltaQ.plev.values, deltaQ, PS_deltaGrid, hyam, hybm)
    else:
        # CMZ TODO
        #deltaTCAM = hyi2hyo(P0, hyam_in, hybm_in, PS_deltaGrid, deltaT, hyam, hybm)
        #deltaQCAM = hyi2hyo(P0, hyam_in, hybm_in, PS_deltaGrid, deltaQ, hyam, hybm)

deltaPSCAM = deltaPS  # no vertical interpolation needed for surface pressure

# Apply smoothing if requested
if smooth_deltas:
    print(f"Applying smoothing with {smooth_delta_iter} iterations...")
    deltaQCAM = pyfuncs.smooth_with_smth9(deltaQCAM, smooth_delta_iter, p=0.5, q=0.25)
    deltaTCAM = pyfuncs.smooth_with_smth9(deltaTCAM, smooth_delta_iter, p=0.5, q=0.25)
    deltaPSCAM = pyfuncs.smooth_with_smth9(deltaPSCAM, smooth_delta_iter, p=0.5, q=0.25)

# Calculate extra diagnostics if requested
if extra_diags_atm:
    dp = cam2cam.dpres_hybrid_ccm(data_vars['ps'], p0, hyai, hybi)
    pw = meteo.prcwater_dp(deltaQCAM, dp)

# Interpolate deltas back to original grid
print("Interpolating deltas back to original grid...")
if esmf_remap:
    # Using existing ESMF regridder
    deltaTCAM_interp = regridder(deltaTCAM)
    deltaQCAM_interp = regridder(deltaQCAM)
    deltaPSCAM_interp = regridder(deltaPSCAM)
else:
    # Using scipy's griddata
    points = np.column_stack((deltaT.lat.values.flatten(), deltaT.lon.values.flatten()))
    xi = np.column_stack((lat.flatten(), lon.flatten()))

    deltaTCAM_interp = np.array([
        griddata(points, deltaTCAM[k].flatten(), xi, method='linear')
        for k in range(nlev)
    ])
    deltaQCAM_interp = np.array([
        griddata(points, deltaQCAM[k].flatten(), xi, method='linear')
        for k in range(nlev)
    ])
    deltaPSCAM_interp = griddata(points, deltaPSCAM.flatten(), xi, method='linear')

# Fill any missing values
deltaTCAM_interp = np.nan_to_num(deltaTCAM_interp, 0)
deltaQCAM_interp = np.nan_to_num(deltaQCAM_interp, 0)
deltaPSCAM_interp = np.nan_to_num(deltaPSCAM_interp, 0)

# Continue with the next part...


