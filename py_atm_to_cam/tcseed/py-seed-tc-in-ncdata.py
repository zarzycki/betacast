import xarray as xr
import numpy as np
from datetime import datetime
import matplotlib.pyplot as plt
from matplotlib import cm
import matplotlib.tri as tri
import argparse

# Betacast modules
module_path = os.path.abspath(os.path.join(os.path.dirname(__file__), '../..', 'py_functions'))
if module_path not in sys.path:
    sys.path.append(module_path)
from py_seedfuncs import *

def wallClockElapseTime(start_time, message):
    elapsed_time = datetime.now() - start_time
    print(f"{message}: {elapsed_time}")

def copy_variable(inputFile, var_name, new_var_name, data, overwrite=False):
    if new_var_name not in inputFile.variables or overwrite:
        print(f"copying {var_name} -> {new_var_name}")
        dims = inputFile[var_name].dims  # Get the dimensions of the original variable
        inputFile[new_var_name] = (dims, data)  # Use passed data
    else:
        print(f"{new_var_name} already exists. Not copying.")

def process_vortex(seedfile, pthi, doplot=False, deg_bnd=15.0):

    invert_vortex  = keyword_values( pthi, "invert_vortex", "bool")
    modify_q = keyword_values( pthi, "modify_q", "bool")
    modify_q_mult = 1.0 if not modify_q else keyword_values(pthi, "modify_q_mult", "float")
    gamma_ = keyword_values(pthi, "gamma", "float")
    cen_lat = keyword_values(pthi, "psminlat", "float")
    cen_lon = keyword_values(pthi, "psminlon", "float")
    zp = keyword_values(pthi, "zp", "float")
    exppr = keyword_values(pthi, "exppr", "float")
    restart_file = keyword_values(pthi, "restart_file", "bool")

    if not invert_vortex:
        minp = keyword_values(pthi, "minp", "float")
        target_rmw = keyword_values(pthi, "target_rmw", "float")

        # Handle default values for minp and target_rmw
        if minp < -19999.0:
            minp = 995.0
        elif 0.0 < minp <= -19999.0:
            minp = abs(minp)

        if target_rmw < 0.0:
            target_rmw = 200000.0

    else:
        rp = keyword_values(pthi, "rp", "float")
        dp = keyword_values(pthi, "dp", "float")

    wcSeedStrt = datetime.now()

    # Open netCDF file
    inputFile = xr.open_dataset(seedfile, mode='a')

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
                return
    else:
        lat = inputFile['lat'].values
        lon = inputFile['lon'].values

    hyai = inputFile['hyai'].values
    hybi = inputFile['hybi'].values
    hyam = inputFile['hyam'].values
    hybm = inputFile['hybm'].values
    P0 = 100000.0

    ncol = lat.shape[0]
    nlev = hyam.shape[0]
    nlevi = hybi.shape[0]

    if restart_file:
        u = inputFile['U'].values
        v = inputFile['V'].values
        ps = inputFile['PSDRY'].values
        t = inputFile['T'].values
        dpq = inputFile['dpQ'].values

        # Calculate dp3d
        dp3d = np.zeros((nlev, ncol))
        for ii in range(ncol):
            dp3d[:, ii] = ((hyai[1:nlevi] - hyai[:nlevi-1]) * P0) + ((hybi[1:nlevi] - hybi[:nlevi-1]) * ps[ii])

        q = dpq / dp3d

        # Pressure correction
        avg_ps_in = np.sum(ps * area) / np.sum(area)
        corrMass = ps + np.sum(inputFile['dpQ'].values, axis=0) + np.sum(inputFile['dpCLDICE'].values, axis=0) + np.sum(inputFile['dpCLDLIQ'].values, axis=0)
        ps = ps + corrMass
        avg_ps_out = np.sum(ps * area) / np.sum(area)

    else:
        u = inputFile['U'].values
        v = inputFile['V'].values
        ps = inputFile['PS'].values
        t = inputFile['T'].values
        q = inputFile['Q'].values

    print(ps.shape)

    u_orig = np.copy(u)
    v_orig = np.copy(v)
    ps_orig = np.copy(ps)
    t_orig = np.copy(t)
    q_orig = np.copy(q)

    # Calculate gcdist
    gcdist, _ = gc_latlon(cen_lat, cen_lon, lat, lon, 2, 2)

    if not invert_vortex:
        minix = np.argmin(gcdist)
        ambps = ps[0,minix]
        print(f"ambient ps: {ambps}")
        dp = ambps - minp * 100.0 if minp > 0.0 else -minp
        rp = get_rp_from_dp_rmw(cen_lat, dp, target_rmw)
        print(rp)

    for ii in range(ncol):
        #if ii % 1000 == 0:
        print(f"At ncol: {ii} of {ncol}")
        if gcdist[ii] <= deg_bnd:
            for kk in range(nlev):
                p = hyam[kk] * P0 + hybm[kk] * ps[0,ii]
                theArr = tctestcase(cen_lon, cen_lat, dp, rp, zp, exppr, gamma_, lon[ii], lat[ii], p, -999, 0, ps[0,ii], u[0,kk, ii], v[0,kk, ii], t[0,kk, ii], q[0,kk, ii], invert_vortex, modify_q, modify_q_mult)
                v[0, kk, ii] = theArr[0]
                u[0, kk, ii] = theArr[1]
                q[0, kk, ii] = theArr[2]
                t[0, kk, ii] = theArr[3]

            ps[0,ii] = theArr[4]

    # Writing back to file
    print(f"Printing file {seedfile}")
    if restart_file:
        q_perc = q / q_orig
        dpq = q_perc * dpq
        corrMass = np.sum(dpq, axis=0) + np.sum(inputFile['dpCLDICE'].values, axis=0) + np.sum(inputFile['dpCLDLIQ'].values, axis=0)
        ps = ps - corrMass
        avg_ps_out = np.sum(ps * area) / np.sum(area)
        ps = ps - (avg_ps_out - avg_ps_in)

        inputFile['PSDRY'].values = ps
        inputFile['U'].values = u
        inputFile['V'].values = v
        inputFile['T'].values = t
        inputFile['dpQ'].values = dpq
    else:
        if 'U_orig' not in inputFile.variables:
            print("copying u -> u_orig")
            copy_variable(inputFile, 'U', 'U_orig', u_orig)
            copy_variable(inputFile, 'V', 'V_orig', v_orig)
            copy_variable(inputFile, 'PS', 'PS_orig', ps_orig)
            copy_variable(inputFile, 'T', 'T_orig', t_orig)
            copy_variable(inputFile, 'Q', 'Q_orig', q_orig)
        elif 'U_orig2' not in inputFile.variables:
            print("copying u -> u_orig2")
            copy_variable(inputFile, 'U', 'U_orig2', u_orig)
            copy_variable(inputFile, 'V', 'V_orig2', v_orig)
            copy_variable(inputFile, 'PS', 'PS_orig2', ps_orig)
            copy_variable(inputFile, 'T', 'T_orig2', t_orig)
            copy_variable(inputFile, 'Q', 'Q_orig2', q_orig)
        else:
            print("We already have U_orig + U_orig2 on file, not adding any more!")

        inputFile['PS'].values = ps
        inputFile['U'].values = u
        inputFile['V'].values = v
        inputFile['T'].values = t
        inputFile['Q'].values = q

    inputFile.to_netcdf(seedfile, mode='a')
    print(f"Successfully written changes to {seedfile}")

    wallClockElapseTime(wcSeedStrt, "Time to seed")

    if doplot:
        plot_diagnostic_map(lon, lat, ps, u)

def plot_diagnostic_map(lon, lat, ps, u):

    triang = tri.Triangulation(lon, lat)

    fig, ax = plt.subplots(figsize=(8, 8))
    contour = ax.tricontourf(triang, u[0,9,:], cmap=cm.viridis)
    plt.colorbar(contour)
    plt.title("U Wind Contour")
    plt.show()

# Add argparse to handle command-line inputs
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Process vortex seed file and namelist')
    parser.add_argument('--se_inic', type=str, required=True, help='Path to the seed file (netCDF format)')
    parser.add_argument('--vortex_namelist', type=str, required=True, help='Path to the vortex namelist file')

    args = parser.parse_args()

    process_vortex(seedfile=args.se_inic, pthi=args.vortex_namelist)

