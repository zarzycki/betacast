import xarray as xr
import numpy as np
from datetime import datetime
import argparse
import os
import sys
import netCDF4 as nc4

# Betacast modules
module_path = os.path.abspath(os.path.join(os.path.dirname(__file__), '../..', 'py_functions'))
if module_path not in sys.path:
    sys.path.append(module_path)
from py_seedfuncs import *

# Parse args
parser = argparse.ArgumentParser(description='Process vortex seed file and namelist')
parser.add_argument('--se_inic', type=str, required=True, help='Path to the seed file (netCDF format)')
parser.add_argument('--vortex_namelist', type=str, required=True, help='Path to the vortex namelist file')
args = parser.parse_args()

# Read relevant settings
seedfile = args.se_inic
pthi = args.vortex_namelist
doplot = False
deg_bnd = 15.0

# Constants
P0 = 100000.0

# Get values from namelist
invert_vortex = keyword_values(pthi, "invert_vortex", "bool", default=False, verbose=True)
modify_q = keyword_values(pthi, "modify_q", "bool", verbose=True)
modify_q_mult = 1.0 if not modify_q else keyword_values(pthi, "modify_q_mult", "float", verbose=True)
gamma_ = keyword_values(pthi, "gamma", "float", verbose=True)
cen_lat = keyword_values(pthi, "psminlat", "float", verbose=True)
cen_lon = keyword_values(pthi, "psminlon", "float", verbose=True)
zp = keyword_values(pthi, "zp", "float", verbose=True)
exppr = keyword_values(pthi, "exppr", "float", verbose=True)
restart_file = keyword_values(pthi, "restart_file", "bool", verbose=True)
adjust_vortex = keyword_values(pthi, "adjust_vortex", "bool", default=False, verbose=True)

if adjust_vortex:
    print("THEARR: ADJUSTING VORTEX!")
    # These are the TARGET values
    minp_obsv = keyword_values(pthi, "minp", "float", verbose=True)
    target_rmw_obsv = keyword_values(pthi, "target_rmw", "float", verbose=True)
    # These are the inversion values
    rp_data = keyword_values(pthi, "rp", "float", verbose=True)
    dp_data = keyword_values(pthi, "dp", "float", verbose=True)
# If we are seeding, read target pressure and RMW
# If we are NOT seeding, assume these have been precalculated or set
else:
    if not invert_vortex:
        minp = keyword_values(pthi, "minp", "float", verbose=True)
        target_rmw = keyword_values(pthi, "target_rmw", "float", verbose=True)

        # Handle default values for minp and target_rmw
        if minp < -19999.0:
            minp = 995.0
        elif 0.0 < minp <= -19999.0:
            minp = abs(minp)

        if target_rmw < 0.0:
            target_rmw = 200000.0
    else:
        rp = keyword_values(pthi, "rp", "float", verbose=True)
        dp = keyword_values(pthi, "dp", "float", verbose=True)

# Start the script timer
wcSeedStrt = datetime.now()

# Load inputFile as xarray now -- using with it will auto close when indent returns
with xr.open_dataset(seedfile) as inputFile:
    # All the code that uses inputFile goes here, indented one level
    if restart_file:
        if 'lat_d' in inputFile.variables:
            lat = inputFile['lat_d'].values
            lon = inputFile['lon_d'].values
            area = inputFile['area_d'].values
        else:
            if 'lat' in inputFile.variables:
                lat = inputFile['lat'].values
                lon = inputFile['lon'].values
                area = inputFile['area'].values
            else:
                print("Cannot find lat or lat_d on input file")
                sys.exit(1)
        u = inputFile['U'].values
        v = inputFile['V'].values
        ps = inputFile['PSDRY'].values
        t = inputFile['T'].values
        dpq = inputFile['dpQ'].values
        dpcldice = inputFile['dpCLDICE'].values
        dpcldliq = inputFile['dpCLDLIQ'].values
    else:
        lat = inputFile['lat'].values
        lon = inputFile['lon'].values
        u = inputFile['U'].values
        v = inputFile['V'].values
        ps = inputFile['PS'].values
        t = inputFile['T'].values
        q = inputFile['Q'].values

    # Common to both restart and non restart file
    hyai = inputFile['hyai'].values
    hybi = inputFile['hybi'].values
    hyam = inputFile['hyam'].values
    hybm = inputFile['hybm'].values

    # Done reading things via xarray...

# Get array shapes
ncol = lat.shape[0]
nlev = hyam.shape[0]
nlevi = hybi.shape[0]

if restart_file:
    # Calculate 3-D delta pressure values
    dp3d = np.zeros((1, nlev, ncol))
    for ii in range(ncol):
        dp3d[0, :, ii] = ((hyai[1:nlevi] - hyai[:nlevi-1]) * P0) + ((hybi[1:nlevi] - hybi[:nlevi-1]) * ps[0, ii])

    # Get q from its dp fraction --- q (kg/kg) = dpq (Pa * kg/kg) / dp3d (Pa)
    q = dpq / dp3d

    # Pressure correction
    avg_ps_in = np.sum(ps * area) / np.sum(area)
    print(f"Average PSDRY in: {avg_ps_in}")

    corrMass = np.sum(dpq, axis=1) + np.sum(dpcldice, axis=1) + np.sum(dpcldliq, axis=1)
    print(f"Initial wet pressure for mass correction: {np.mean(corrMass)}")
    ps = ps + corrMass
    print(f"Total PS after adding condensates: {np.sum(ps * area) / np.sum(area)}")

    avg_ps_out = np.sum(ps * area) / np.sum(area)

    # Before seeding - store original integral values
    T_global_orig = calc_mass_weighted_integral(t, dp3d, area)
    Q_global_orig = calc_mass_weighted_integral(q, dp3d, area)
    U_global_orig = calc_mass_weighted_integral(u, dp3d, area)
    V_global_orig = calc_mass_weighted_integral(v, dp3d, area)
    print(f"INTEGRAL: BEFORE  -- T: {T_global_orig:25.17e}, Q: {Q_global_orig:25.17e}, U: {U_global_orig:25.17e}, V: {V_global_orig:25.17e}")

# Create copies of unmodified fields to query later if needed
u_orig = np.copy(u)
v_orig = np.copy(v)
ps_orig = np.copy(ps)
t_orig = np.copy(t)
q_orig = np.copy(q)

# Calculate great circle distances
gcdist, _ = gc_latlon(cen_lat, cen_lon, lat, lon, 2, 2)

# Find relevant parameters for TC seed
if adjust_vortex:
    print("GETRPDP: ADJUSTING VORTEX!")
    minix = np.argmin(gcdist)
    ambps = ps[0, minix] + dp_data  # ambient PS is "minimum PS" + dp_data calculated before
    print(f"ambient ps: {ambps}")
    dp_obsv = ambps - minp_obsv * 100.0 if minp_obsv > 0.0 else -minp_obsv
    rp_obsv = get_rp_from_dp_rmw(cen_lat, dp_obsv, target_rmw_obsv)
    print(f"rp: {rp_obsv}  -- dp: {dp_obsv}")
elif not invert_vortex:
    minix = np.argmin(gcdist)
    ambps = ps[0, minix]
    print(f"ambient ps: {ambps}")
    dp = ambps - minp * 100.0 if minp > 0.0 else -minp
    rp = get_rp_from_dp_rmw(cen_lat, dp, target_rmw)
    print(f"rp: {rp}  -- dp: {dp}")

# Add the vortex seed to the state fields
for ii in range(ncol):
    if gcdist[ii] <= deg_bnd:
        for kk in range(nlev):
            p = hyam[kk] * P0 + hybm[kk] * ps[0, ii]
            if not adjust_vortex:
                # Note, -999 is a failsafe for "z" since zcoords = 0, should be pressure.
                theArr = tctestcase(cen_lon, cen_lat, dp, rp, zp, exppr, gamma_, lon[ii], lat[ii], p, -999, 0, ps[0, ii], u[0, kk, ii], v[0, kk, ii], t[0, kk, ii], q[0, kk, ii], invert_vortex, modify_q, modify_q_mult)
                v[0, kk, ii] = theArr[0]
                u[0, kk, ii] = theArr[1]
                q[0, kk, ii] = theArr[2]
                t[0, kk, ii] = theArr[3]
                if kk == 0:
                    ps[0, ii] = theArr[4]
            else:
                anom_data = tctestcase(cen_lon, cen_lat, dp_data, rp_data, zp, exppr, gamma_, lon[ii], lat[ii], p, -999, 0, ps[0, ii], u[0, kk, ii], v[0, kk, ii], t[0, kk, ii], q[0, kk, ii], True, modify_q, modify_q_mult, return_anoms=True)
                anom_obsv = tctestcase(cen_lon, cen_lat, dp_obsv, rp_obsv, zp, exppr, gamma_, lon[ii], lat[ii], p, -999, 0, ps[0, ii], u[0, kk, ii], v[0, kk, ii], t[0, kk, ii], q[0, kk, ii], False, modify_q, modify_q_mult, return_anoms=True)
                v[0, kk, ii] = v[0, kk, ii] + anom_obsv[0] + anom_data[0]
                u[0, kk, ii] = u[0, kk, ii] + anom_obsv[1] + anom_data[1]
                q[0, kk, ii] = q[0, kk, ii] + anom_obsv[2] + anom_data[2]
                t[0, kk, ii] = t[0, kk, ii] + anom_obsv[3] + anom_data[3]
                if kk == 0:
                    ps[0, ii] = ps[0, ii] + anom_obsv[4] + anom_data[4]

# Correct mass
if restart_file:
    q_perc = q / q_orig
    print(f"qperc min/max: {np.min(q_perc)} {np.max(q_perc)}")
    dpq = q_perc * dpq
    corrMass = np.sum(dpq, axis=1) + np.sum(dpcldice, axis=1) + np.sum(dpcldliq, axis=1)
    print(f"New wet pressure for mass correction: {np.mean(corrMass)}")
    ps = ps - corrMass
    avg_ps_out = np.sum(ps * area) / np.sum(area)
    print(f"Average PSDRY out: {avg_ps_out}")
    ps = ps - (avg_ps_out - avg_ps_in)
    avg_ps_final = np.sum(ps * area) / np.sum(area)
    print(f"Final PSDRY out: {avg_ps_final}")

    # Update dp3d to reflect new psdry
    for ii in range(ncol):
        dp3d[0, :, ii] = ((hyai[1:nlevi] - hyai[:nlevi-1]) * P0) + ((hybi[1:nlevi] - hybi[:nlevi-1]) * ps[0, ii])

    # Calculate new global values
    T_global_new = calc_mass_weighted_integral(t, dp3d, area)
    Q_global_new = calc_mass_weighted_integral(q, dp3d, area)
    U_global_new = calc_mass_weighted_integral(u, dp3d, area)
    V_global_new = calc_mass_weighted_integral(v, dp3d, area)
    print(f"INTEGRAL: AFTER   -- T: {T_global_new:25.17e}, Q: {Q_global_new:25.17e}, U: {U_global_new:25.17e}, V: {V_global_new:25.17e}")

    # Calculate global integral diffs
    T_global_diff = T_global_orig - T_global_new
    Q_global_diff = Q_global_orig - Q_global_new
    U_global_diff = U_global_orig - U_global_new
    V_global_diff = V_global_orig - V_global_new
    print(f"INTEGRAL: DIFF   -- T: {T_global_diff:25.17e}, Q: {Q_global_diff:25.17e}, U: {U_global_diff:25.17e}, V: {V_global_diff:25.17e}")

    # Update state variables by adding a correction factor to enforce integrals are the same
    t[0,:,:] = t[0,:,:] + calc_inverse_mass_weighted_field(T_global_diff, dp3d, area)
    q[0,:,:] = q[0,:,:] + calc_inverse_mass_weighted_field(Q_global_diff, dp3d, area)
    u[0,:,:] = u[0,:,:] + calc_inverse_mass_weighted_field(U_global_diff, dp3d, area)
    v[0,:,:] = v[0,:,:] + calc_inverse_mass_weighted_field(V_global_diff, dp3d, area)

    # After modifications - calculate new values
    T_global_new = calc_mass_weighted_integral(t, dp3d, area)
    Q_global_new = calc_mass_weighted_integral(q, dp3d, area)
    U_global_new = calc_mass_weighted_integral(u, dp3d, area)
    V_global_new = calc_mass_weighted_integral(v, dp3d, area)
    print(f"INTEGRAL: FINAL   -- T: {T_global_new:25.17e}, Q: {Q_global_new:25.17e}, U: {U_global_new:25.17e}, V: {V_global_new:25.17e}")

    # Verify differences
    T_global_diff = T_global_orig - T_global_new
    Q_global_diff = Q_global_orig - Q_global_new
    U_global_diff = U_global_orig - U_global_new
    V_global_diff = V_global_orig - V_global_new
    print(f"INTEGRAL: FINAL_DIFF   -- T: {T_global_diff:25.17e}, Q: {Q_global_diff:25.17e}, U: {U_global_diff:25.17e}, V: {V_global_diff:25.17e}")

# Writing back to file
# Note, we use nc4 here to write instead of xarray because it's easier to preserve
# netcdf type and other encodings. Also faster.
print(f"Printing file {seedfile}")

if restart_file:
    # Write using netCDF4
    with nc4.Dataset(seedfile, 'r+') as nc_file:
        nc_file.variables['PSDRY'][:] = ps
        nc_file.variables['U'][:] = u
        nc_file.variables['V'][:] = v
        nc_file.variables['T'][:] = t
        nc_file.variables['dpQ'][:] = dpq
else:
    with nc4.Dataset(seedfile, 'r+') as nc_file:
        # Update main variables
        nc_file.variables['PS'][:] = ps
        nc_file.variables['U'][:] = u
        nc_file.variables['V'][:] = v
        nc_file.variables['T'][:] = t
        nc_file.variables['Q'][:] = q

        # Helper function for copying variables with netCDF4
        def copy_if_not_exists(nc_file, var_name, new_var_name, data):
            if new_var_name not in nc_file.variables:
                print(f"copying {var_name} -> {new_var_name}")
                orig_var = nc_file.variables[var_name]
                new_var = nc_file.createVariable(new_var_name, orig_var.datatype, orig_var.dimensions)
                new_var[:] = data
            else:
                print(f"{new_var_name} already exists. Not copying.")

        # Add orig variables if needed
        if 'U_orig' not in nc_file.variables:
            copy_if_not_exists(nc_file, 'U', 'U_orig', u_orig)
            copy_if_not_exists(nc_file, 'V', 'V_orig', v_orig)
            copy_if_not_exists(nc_file, 'PS', 'PS_orig', ps_orig)
            copy_if_not_exists(nc_file, 'T', 'T_orig', t_orig)
            copy_if_not_exists(nc_file, 'Q', 'Q_orig', q_orig)
        elif 'U_orig2' not in nc_file.variables:
            copy_if_not_exists(nc_file, 'U', 'U_orig2', u_orig)
            copy_if_not_exists(nc_file, 'V', 'V_orig2', v_orig)
            copy_if_not_exists(nc_file, 'PS', 'PS_orig2', ps_orig)
            copy_if_not_exists(nc_file, 'T', 'T_orig2', t_orig)
            copy_if_not_exists(nc_file, 'Q', 'Q_orig2', q_orig)
        else:
            print("We already have U_orig + U_orig2 on file, not adding any more!")

print(f"Successfully written changes to {seedfile}")

elapsed_time = datetime.now() - wcSeedStrt
print(f"Time to seed: {elapsed_time}")
