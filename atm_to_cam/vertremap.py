import numpy as np
#from scipy.interpolate import interp1d
from numba import jit
#from tqdm import tqdm
import time

# def interpolate_to_hybrid_levels(p_levels, data_on_p_levels, ps, a_coeff, b_coeff, p0=100000):
#     """
#     Interpolate data on constant pressure levels to hybrid sigma levels.
#
#     Parameters:
#     p_levels: 1D array of pressure levels (Pa) on which data is currently defined
#     data_on_p_levels: 3D array of data defined on pressure levels (e.g., T, U, V)
#     ps: 2D array of surface pressure (Pa)
#     a_coeff: 1D array of a(k) hybrid coefficients
#     b_coeff: 1D array of b(k) hybrid coefficients
#     p0: Reference pressure, typically 100000 Pa (1000 hPa)
#
#     Returns:
#     data_on_hybrid_levels: 4D array of data interpolated to hybrid sigma levels
#     """
#     start_time = time.time()
#
#     # Compute hybrid pressure levels
#     hybrid_p = np.array([a * p0 + b * ps for a, b in zip(a_coeff, b_coeff)])
#
#     # Prepare an output array for the interpolated data
#     nlevs = len(a_coeff)
#     nlat, nlon = ps.shape
#     data_on_hybrid_levels = np.zeros((nlevs, nlat, nlon))
#
#     # Interpolate for each lat/lon point with progress bar
#     for i in tqdm(range(nlat), desc="Interpolating", unit="lat"):
#         for j in range(nlon):
#             f_interp = interp1d(p_levels, data_on_p_levels[:, i, j], fill_value="extrapolate")
#             data_on_hybrid_levels[:, i, j] = f_interp(hybrid_p[:, i, j])
#
#     # Print elapsed time
#     elapsed_time = time.time() - start_time
#     print(f"Interpolation completed in {elapsed_time:.2f} seconds.")
#
#     return data_on_hybrid_levels

#### NUMBA

@jit(nopython=True)
def interpolate_to_hybrid_levels_numba(p_levels, data_on_p_levels, ps, hybrid_p):
    """
    JIT-compiled function to interpolate data on constant pressure levels to hybrid sigma levels.
    This function is optimized with numba.
    """
    nlevs, nlat, nlon = hybrid_p.shape
    data_on_hybrid_levels = np.zeros((nlevs, nlat, nlon))

    # Interpolate for each lat/lon point
    for i in range(nlat):
        for j in range(nlon):
            data_on_hybrid_levels[:, i, j] = np.interp(hybrid_p[:, i, j], p_levels, data_on_p_levels[:, i, j])

    return data_on_hybrid_levels

@jit(nopython=True)
def p2hyo(p_levels, data_on_p_levels, ps, a_coeff, b_coeff, p0, kflag):
    """
    Python equivalent of the P2HYO Fortran routine from NCL
    """
    nlevi = p_levels.size
    nlat, nlon = ps.shape
    nlevo = a_coeff.size
    data_on_hybrid_levels = np.full((nlevo, nlat, nlon), np.nan)

    # Compute hybrid pressure levels
    hybrid_p = np.zeros((nlevo, nlat, nlon))
    for i in range(nlevo):
        hybrid_p[i, :, :] = a_coeff[i] * p0 + b_coeff[i] * ps

    pimin = np.log(p_levels[0])
    pimax = np.log(p_levels[-1])

    for i in range(nlat):
        for j in range(nlon):
            for ko in range(nlevo):
                po_log = np.log(hybrid_p[ko, i, j])
                if po_log < pimin:
                    if kflag in [0, 1, 3]:
                        if kflag == 0:
                            continue  # Set to missing
                        elif kflag == 1:
                            data_on_hybrid_levels[ko, i, j] = data_on_p_levels[0, i, j]
                        elif kflag == 3:
                            dxdp = (data_on_p_levels[1, i, j] - data_on_p_levels[0, i, j]) / (np.log(p_levels[1]) - pimin)
                            data_on_hybrid_levels[ko, i, j] = data_on_p_levels[0, i, j] + (po_log - pimin) * dxdp
                    else:
                        dxdp = (data_on_p_levels[1, i, j] - data_on_p_levels[0, i, j]) / (np.log(p_levels[1]) - pimin)
                        data_on_hybrid_levels[ko, i, j] = data_on_p_levels[0, i, j] + (po_log - pimin) * dxdp
                elif po_log > pimax:
                    if kflag in [0, 1, 2]:
                        if kflag == 0:
                            continue  # Set to missing
                        elif kflag == 1:
                            data_on_hybrid_levels[ko, i, j] = data_on_p_levels[-1, i, j]
                        elif kflag == 2:
                            dxdp = (data_on_p_levels[-1, i, j] - data_on_p_levels[-2, i, j]) / (pimax - np.log(p_levels[-2]))
                            data_on_hybrid_levels[ko, i, j] = data_on_p_levels[-1, i, j] + (po_log - pimax) * dxdp
                    else:
                        dxdp = (data_on_p_levels[-1, i, j] - data_on_p_levels[-2, i, j]) / (pimax - np.log(p_levels[-2]))
                        data_on_hybrid_levels[ko, i, j] = data_on_p_levels[-1, i, j] + (po_log - pimax) * dxdp
                else:
                    for ki in range(nlevi - 1):
                        if np.log(p_levels[ki]) <= po_log < np.log(p_levels[ki + 1]):
                            dxdp = (data_on_p_levels[ki + 1, i, j] - data_on_p_levels[ki, i, j]) / (np.log(p_levels[ki + 1]) - np.log(p_levels[ki]))
                            data_on_hybrid_levels[ko, i, j] = data_on_p_levels[ki, i, j] + (po_log - np.log(p_levels[ki])) * dxdp
                            break

    return data_on_hybrid_levels

def interpolate_to_hybrid_levels(p_levels, data_on_p_levels, ps, a_coeff, b_coeff, p0=100000, kflag=1):
    """
    Interpolate data on constant pressure levels to hybrid sigma levels.

    This function is a wrapper that prepares the input for the JIT-compiled version.
    """

    # Start timing
    start_time = time.time()

    # Check if p_levels is monotonically increasing
    if not np.all(np.diff(p_levels) > 0):
       raise ValueError("p_levels must be monotonically increasing")

    # Call the Python-equivalent Fortran function
    data_on_hybrid_levels = p2hyo(p_levels, data_on_p_levels, ps, a_coeff, b_coeff, p0, kflag)

    # Print elapsed time
    elapsed_time = time.time() - start_time
    print(f"Interpolation completed in {elapsed_time:.2f} seconds.")

    return data_on_hybrid_levels

# def interpolate_to_hybrid_levels(p_levels, data_on_p_levels, ps, a_coeff, b_coeff, p0=100000):
#     """
#     Interpolate data on constant pressure levels to hybrid sigma levels.
#
#     This function is a wrapper that prepares the input for the JIT-compiled version.
#     """
#     # Start timing
#     start_time = time.time()
#
#     # Compute hybrid pressure levels
#     hybrid_p = np.array([a * p0 + b * ps for a, b in zip(a_coeff, b_coeff)])
#
#     # Check if p_levels is monotonically increasing
#     if not np.all(np.diff(p_levels) > 0):
#         raise ValueError("p_levels must be monotonically increasing")
#
#     # Call the JIT-compiled function
#     data_on_hybrid_levels = interpolate_to_hybrid_levels_numba(p_levels, data_on_p_levels, ps, hybrid_p)
#
#     # Print elapsed time
#     elapsed_time = time.time() - start_time
#     print(f"Interpolation completed in {elapsed_time:.2f} seconds.")
#
#     return data_on_hybrid_levels
#
#