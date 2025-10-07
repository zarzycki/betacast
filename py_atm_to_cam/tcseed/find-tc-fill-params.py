import argparse
import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
import shutil
import os
import sys

# Betacast modules
module_path = os.path.abspath(os.path.join(os.path.dirname(__file__), '../..', 'py_functions'))
if module_path not in sys.path:
    sys.path.append(module_path)
from py_seedfuncs import keyword_values, gc_latlon, tctestcase, radialAvg2D_unstruc, radialAvg3D_unstruc, replace_or_add_variable

# Set up argument parsing
parser = argparse.ArgumentParser(description="Process vortex with seeding.")
parser.add_argument('--inic_file', required=True, help='Path to the initialization netCDF file')
parser.add_argument('--vortex_namelist', required=True, help='Path to the namelist file')
args = parser.parse_args()

# Initialize variables from namelist
pthi = args.vortex_namelist
inic_file = args.inic_file

invert_vortex = keyword_values(pthi, "invert_vortex", "bool")

if not invert_vortex:
    print("Not inverting vortex, we can exit gracefully!")
    exit()

deltaMax = keyword_values(pthi, "deltaMax", "float")
psminlat = keyword_values(pthi, "psminlat", "float")
psminlon = keyword_values(pthi, "psminlon", "float")
modify_q = keyword_values(pthi, "modify_q", "bool")
modify_q_mult = keyword_values(pthi, "modify_q_mult", "float")
gamma_ = keyword_values(pthi, "gamma", "float")
restart_file = keyword_values(pthi, "restart_file", "bool")

if gamma_ < 0.0:
    gamma_ = 0.0065
    print(f"setting gamma_ to default of {gamma_}")

# ================== INTERNAL SETTINGS ===================
truePS_scan_radius = 300.0  # km
rad_for_corr = 800.0  # km

n_p_steps = 18
rp_arr = np.linspace(80000., 300000., n_p_steps)
dp_arr = np.linspace(200., 6000., n_p_steps)
exppr_arr = np.linspace(1.1, 1.9, n_p_steps)

n_t_steps = 60
zp_arr = np.linspace(6000., 16000., n_t_steps)

pi = np.pi
convert = 180.0 / pi

# Load a regridded state field with the TC in it.
f = xr.open_dataset(inic_file)

if restart_file:
    lat, lon = (f['lat_d'].values, f['lon_d'].values) if 'lat_d' in f.variables else (f['lat'].values, f['lon'].values)
    hybm = f['hybm'].values
    hyam = f['hyam'].values
    P0 = 100000.0
    lev = hyam * P0 + hybm * np.mean(f['PSDRY'].values)
    lev = lev / 100.0  # Convert to mb
else:
    lat, lon = f['lat'].values, f['lon'].values
    lev = f['lev'].values

ncol = len(lat)
nlev = len(lev)

if restart_file:
    ps = f['PSDRY'].values
    corrMass = ps + np.sum(f['dpQ'].values, axis=1) + np.sum(f['dpCLDICE'].values, axis=1) + np.sum(f['dpCLDLIQ'].values, axis=1)
    ps += corrMass
    mps = ps[0, :]
    T = f['T'][0, :, :].values
else:
    mps = f['PS'][0, :].values
    T = f['T'][0, :, :].values

# Find the minimum location of psminlat
gcdist, _ = gc_latlon(psminlat, psminlon, lat, lon, 2, 4)
tmpps = np.where(gcdist > truePS_scan_radius, np.nan, mps)
minix = np.nanargmin(tmpps)
psminlat, psminlon = lat[minix], lon[minix]

print(lat)
print(lon)
print(mps)
print(T[0,:])

# Radial profiles of pressure and temperature
rad_psl = radialAvg2D_unstruc(mps, lat, lon, deltaMax, psminlat, psminlon, rad_for_corr, False)
Tavg = radialAvg3D_unstruc(T, lat, lon, lev, deltaMax, psminlat, psminlon, rad_for_corr, False)
print(Tavg)
print("---")

print(f"psminlat {psminlat}, psminlon {psminlon}")
print(rad_psl['radial_average'])
print(rad_psl['radius'])
print(rad_psl['hit_count'])

nrad_T = Tavg['radial_average'].shape[1]
Tenv = Tavg['radial_average'][:, nrad_T - 1]
Tanom = Tavg['radial_average'] - Tenv[:, np.newaxis]
Tavg = Tanom
print(Tavg)
print(lev)
print(Tenv)
print("***")

# ===================OPTIMIZATION ROUND 1==============
wcOptimStrt = "now"  # Placeholder for the start time of optimization

# Convert radius associated with radial averaging to degrees along a line of constant latitude
mylon = rad_psl['radius'] / (111.0 * np.cos(np.radians(psminlat)))
psl_amb = rad_psl['radial_average'][nrad_T - 1]
print(psl_amb)

anlPS = np.full_like(mylon, 100000.0)
anlT = np.zeros((nlev, len(mylon)))

# Define the center of the TC
cen_lat = psminlat
cen_lon = 0.0
zp = 10000.0  # Initial guess for zp

invert_vortex = False

# Determine the range of rp, dp, and exppr for optimization
rp_ext = rp_arr[-1] - rp_arr[0]
dp_ext = dp_arr[-1] - dp_arr[0]
exppr_ext = exppr_arr[-1] - exppr_arr[0]
print(f"orig rp range: {rp_ext}    orig dp range: {dp_ext}     orig exppr range: {exppr_ext}")

# Initialize "skill" array for optimization
rmsd = np.zeros((1 + 3, int(n_p_steps ** 3)))

# Initialize linear counter
counter = 0

print(rp_arr)
print(dp_arr)
print(exppr_arr)
print(mylon)
print(psl_amb)

# Optimization loop over PS variables
for aa in range(len(rp_arr)):
    for bb in range(len(dp_arr)):
        for cc in range(len(exppr_arr)):
            rp, dp, exppr = rp_arr[aa], dp_arr[bb], exppr_arr[cc]

            # Use the "mylon" strip to calculate the analytic TC profile at each longitude along the cen_lat line
            for ii in range(len(mylon)):
                tmp = tctestcase(cen_lon, cen_lat, dp, rp, zp, exppr, gamma_, mylon[ii], cen_lat, lev[-1] * 100.0, -999, 0, psl_amb, 0.0, 0.0, Tenv[-1], 0.0, invert_vortex, modify_q, modify_q_mult)
                anlPS[ii] = tmp[4]

            # Root Mean Square Deviation (RMSD)
            rmsd[0, counter] = np.sqrt(np.mean((anlPS - rad_psl['radial_average']) ** 2))
            rmsd[1, counter] = rp
            rmsd[2, counter] = dp
            rmsd[3, counter] = exppr
            counter += 1

# Find "best" configuration from round 1
best_config = np.argmin(rmsd[0, :])
print(f"BEST: {rmsd[0, best_config]} {rmsd[1, best_config]} {rmsd[2, best_config]} {rmsd[3, best_config]}")

# ===================OPTIMIZATION ROUND 2==============
# Shrink the space by reducing the search range and recalculating rp, dp, and exppr

mult_factor = 0.52
rp_ext = (rp_arr[1] - rp_arr[0]) * mult_factor
dp_ext = (dp_arr[1] - dp_arr[0]) * mult_factor
exppr_ext = (exppr_arr[1] - exppr_arr[0]) * mult_factor
print(f"new rp range: {2. * rp_ext}    new dp range: {2. * dp_ext}     new exppr range: {2. * exppr_ext}")

# Create updated arrays based on best configuration from round 1
rp_arr = np.linspace(rmsd[1, best_config] - rp_ext, rmsd[1, best_config] + rp_ext, n_p_steps)
dp_arr = np.linspace(rmsd[2, best_config] - dp_ext, rmsd[2, best_config] + dp_ext, n_p_steps)
exppr_arr = np.linspace(rmsd[3, best_config] - exppr_ext, rmsd[3, best_config] + exppr_ext, n_p_steps)
rmsd = np.zeros((1 + 3, int(n_p_steps ** 3)))  # Reinitialize skill array
counter = 0

# Re-run optimization loop with updated ranges
for aa in range(len(rp_arr)):
    for bb in range(len(dp_arr)):
        for cc in range(len(exppr_arr)):
            rp, dp, exppr = rp_arr[aa], dp_arr[bb], exppr_arr[cc]

            for ii in range(len(mylon)):
                tmp = tctestcase(cen_lon, cen_lat, dp, rp, zp, exppr, gamma_, mylon[ii], cen_lat, lev[-1] * 100.0, -999, 0, psl_amb, 0.0, 0.0, Tenv[-1], 0.0, invert_vortex, modify_q, modify_q_mult)
                anlPS[ii] = tmp[4]

            rmsd[0, counter] = np.sqrt(np.mean((anlPS - rad_psl['radial_average']) ** 2))
            rmsd[1, counter] = rp
            rmsd[2, counter] = dp
            rmsd[3, counter] = exppr
            counter += 1

# Find "best" configuration from round 2
best_config = np.argmin(rmsd[0, :])
print(f"BEST: {rmsd[0, best_config]} {rmsd[1, best_config]} {rmsd[2, best_config]} {rmsd[3, best_config]}")

# Save settings from the "best" config for the T optimization
rp = rmsd[1, best_config]
dp = rmsd[2, best_config]
exppr = rmsd[3, best_config]

# ===================OPTIMIZATION T ANOM==============
# Initialize skill array for T anomaly optimization
corr = np.zeros((1 + 1, int(n_t_steps ** 1)))  # For zp only
counter = 0

# Loop over T variables
for dd in range(len(zp_arr)):
    zp = zp_arr[dd]

    for ii in range(len(mylon)):
        for jj in range(nlev):
            tmp = tctestcase(cen_lon, cen_lat, dp, rp, zp, exppr, gamma_, mylon[ii], cen_lat, lev[jj] * 100.0, -999, 0, psl_amb, 0.0, 0.0, Tenv[jj], 0.0, invert_vortex, modify_q, modify_q_mult)
            anlT[jj, ii] = tmp[3]

    # Calculate the anomaly from the TC test case by taking the T at the "center" of the vortex and T at the end of the lat strip
    anlT_anom = anlT[:, 0] - anlT[:, -1]

    # Correlation skill calculation
    print(anlT_anom.shape)
    print(Tavg[:, 0].shape)
    corr[0, counter] = np.corrcoef(anlT_anom, Tavg[:, 0])[0, 1]
    corr[1, counter] = zp
    counter += 1

# Find best configuration for T anomaly
nlev = corr.shape[1]  # Get the number of levels
for i in range(nlev):
    print(f"{corr[0, i]:.6f} {corr[1, i]:.2f}")

best_config = np.argmax(corr[0, :])
print(f"BEST: {corr[0, best_config]} {corr[1, best_config]}")
zp = corr[1, best_config]

# Final best settings:
print(f"Final BEST SETTINGS: rp = {rp}, dp = {dp}, zp = {zp}, exppr = {exppr}")

nlev = len(lev)  # Number of levels
ncol = len(mylon)  # Number of radial points (mylon)
anlT = np.zeros((nlev, ncol))  # Size: (nlev, len(mylon))
anlPS = np.zeros(ncol)  # Size: (len(mylon))

# Iterate over radial points and levels
for ii in range(ncol):
    for jj in range(nlev):
        # Call the tctestcase function and unpack the result into anlT and anlPS
        tmp = tctestcase(cen_lon, cen_lat, dp, rp, zp, exppr, gamma_, mylon[ii], cen_lat, lev[jj] * 100.0,
                         -999, 0, psl_amb, 0.0, 0.0, Tenv[jj], 0.0, invert_vortex, modify_q, modify_q_mult)

        # Assign the temperature at level jj and radius ii
        anlT[jj, ii] = tmp[3]  # Equivalent to tmp(3) in NCL (temperature)

    # Assign the pressure at radius ii
    anlPS[ii] = tmp[4]  # Equivalent to tmp(4) in NCL (pressure)

# Compute the temperature anomaly
anlT_anom = anlT[:, 0] - anlT[:, -1]  # Subtract the first column from the last column

# Save settings from the "best" config for the T optimization
print("BEST SETTINGS --------------------------")
print(f"rp: {rp}")
print(f"dp: {dp}")
print(f"zp: {zp}")
print(f"gamma_: {gamma_}")
print(f"exppr: {exppr}")
print(f"cen_lat: {psminlat}")
print(f"cen_lon: {psminlon}")
print(f"modify_q: {modify_q}")
print(f"modify_q_mult: {modify_q_mult}")
print("--------------------------")

# Create a backup of the namelist file
shutil.copy(pthi, f"{pthi}.BAK")

# Replace or add new variables in the namelist file
replace_or_add_variable(pthi, "rp", rp)
replace_or_add_variable(pthi, "dp", dp)
replace_or_add_variable(pthi, "zp", zp)
replace_or_add_variable(pthi, "exppr", exppr)
replace_or_add_variable(pthi, "psminlat", psminlat)
replace_or_add_variable(pthi, "psminlon", psminlon)

print("Namelist variables updated successfully!")

# ===================PLOTTING==============
fig, (ax1, ax2) = plt.subplots(nrows=2, figsize=(8, 12))  # Two subplots, one above the other

# First panel: Radial Pressure Profile
ax1.plot(rad_psl['radius'], rad_psl['radial_average'], label='Radial Avg Pressure', color='black')
ax1.plot(rad_psl['radius'], anlPS, label='Analytic Pressure', linestyle='--', color='red')
ax1.set_title('Radial Pressure Profile')
ax1.set_xlabel('Radius (km)')
ax1.set_ylabel('Surface pressure')
ax1.legend()
ax1.set_xlim(0, 900)
ax1.set_ylim(95000, 102000)

# Second panel: Temperature Anomaly Profile
ax2.plot(Tavg[:,0], lev, label='Temperature Anomaly', color='black')
ax2.plot(anlT_anom, lev, color='red')  # Assuming this is the comparison line
ax2.invert_yaxis()
ax2.set_title('Temperature Anomaly Profile')
ax2.set_xlabel('Temperature Anomaly (K)')
ax2.set_ylabel('hybrid level at midpoints (1000*(A+B))')
ax2.set_xlim(-10, 10)
ax2.set_ylim(1010, 100)

# Adjust the layout for better spacing between the subplots
plt.tight_layout()

# Save the figure as PNG
output_path = 'best_fit_profiles.png'
plt.savefig(output_path)