import xarray as xr
import numpy as np
from datetime import datetime
import matplotlib.pyplot as plt
from matplotlib import cm

# Assuming you have equivalent Python functions for `keyword_values`, `gc_latlon`, `get_rp_from_dp_rmw`, `tctestcase`, etc.

def wallClockElapseTime(start_time, message):
    elapsed_time = datetime.now() - start_time
    print(f"{message}: {elapsed_time}")

def process_vortex(invert_vortex, modify_q, seedfile, pthi, doplot=False, deg_bnd=15.0):
    # Load namelist values
    modify_q_mult = 1.0 if not modify_q else keyword_values(pthi, "modify_q_mult", "float")
    gamma_ = keyword_values(pthi, "gamma", "float")
    cen_lat = keyword_values(pthi, "psminlat", "float")
    cen_lon = keyword_values(pthi, "psminlon", "float")
    zp = keyword_values(pthi, "zp", "float")
    exppr = keyword_values(pthi, "exppr", "float")
    restart_file = keyword_values(pthi, "restart_file", "logical")
    
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
    inputFile = xr.open_dataset(seedfile)

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

    u_orig = np.copy(u)
    v_orig = np.copy(v)
    ps_orig = np.copy(ps)
    t_orig = np.copy(t)
    q_orig = np.copy(q)

    # Calculate gcdist
    gcdist = gc_latlon(cen_lat, cen_lon, lat, lon, 2, 2)

    if not invert_vortex:
        minix = np.argmin(gcdist)
        ambps = ps[minix]
        dp = ambps - minp * 100.0 if minp > 0.0 else -minp
        rp = get_rp_from_dp_rmw(cen_lat, dp, target_rmw)

    for ii in range(ncol):
        if ii % 1000 == 0:
            print(f"At ncol: {ii} of {ncol}")
        if gcdist[ii] <= deg_bnd:
            for kk in range(nlev):
                p = hyam[kk] * P0 + hybm[kk] * ps[ii]
                theArr = tctestcase(cen_lon, cen_lat, dp, rp, zp, exppr, gamma_, lon[ii], lat[ii], p, -999, 0, ps[ii], u[kk, ii], v[kk, ii], t[kk, ii], q[kk, ii], invert_vortex, modify_q, modify_q_mult)
                v[kk, ii] = theArr[0]
                u[kk, ii] = theArr[1]
                q[kk, ii] = theArr[2]
                t[kk, ii] = theArr[3]

            ps[ii] = theArr[4]

    # Writing back to file
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
        inputFile['PS'].values = ps
        inputFile['U'].values = u
        inputFile['V'].values = v
        inputFile['T'].values = t
        inputFile['Q'].values = q

    wallClockElapseTime(wcSeedStrt, "Time to seed")

    # Optional plotting
    if doplot:
        plot_diagnostic_map(lon, lat, ps, u)

def plot_diagnostic_map(lon, lat, ps, u):
    fig, ax = plt.subplots(figsize=(8, 8))
    contour = ax.contourf(lon, lat, ps, cmap=cm.viridis)
    plt.colorbar(contour)
    plt.title("PS Contour")
    plt.show()

    fig, ax = plt.subplots(figsize=(8, 8))
    contour = ax.contourf(lon, lat, u[26], cmap=cm.viridis)
    plt.colorbar(contour)
    plt.title("U Wind Contour at Level 26")
    plt.show()

# Run the function
pthi = "../../namelists/unseed.default.nl"
seedfile = "/path/to/seedfile.nc"
process_vortex(invert_vortex=True, modify_q=True, seedfile=seedfile, pthi=pthi)

