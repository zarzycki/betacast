import numpy as np
#from scipy.interpolate import interp1d
from numba import jit
#from tqdm import tqdm
import time

# def pressure_to_hybrid(p_levels, data_on_p_levels, ps, a_coeff, b_coeff, p0=100000):
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

# @jit(nopython=True)
# def pressure_to_hybrid_numba(p_levels, data_on_p_levels, ps, hybrid_p):
#     """
#     JIT-compiled function to interpolate data on constant pressure levels to hybrid sigma levels.
#     This function is optimized with numba.
#     """
#     nlevs, nlat, nlon = hybrid_p.shape
#     data_on_hybrid_levels = np.zeros((nlevs, nlat, nlon))
#
#     # Interpolate for each lat/lon point
#     for i in range(nlat):
#         for j in range(nlon):
#             data_on_hybrid_levels[:, i, j] = np.interp(hybrid_p[:, i, j], p_levels, data_on_p_levels[:, i, j])
#
#     return data_on_hybrid_levels

@jit(nopython=True)
def int2p(pin, xin, pout, linlog):
    """
    Interpolates data from one set of pressure levels to another.

    Parameters:
    -----------
    pin : array-like
        Array of input pressure levels. If multi-dimensional, the level dimension must be in the rightmost position.

    xin : array-like
        Array of data to be interpolated. Should have the same shape as `pin`.

    pout : array-like
        Array of output pressure levels. If multi-dimensional, the level dimension must be in the rightmost position
        and all other dimensions must match `xin`. If one-dimensional, all of `xin` will be interpolated to the same levels.

    linlog : int
        Integer indicating the type of interpolation:
        - abs(linlog) == 1 --> linear interpolation
        - abs(linlog) != 1 --> logarithmic interpolation
        - If linlog is negative, extrapolation outside the range of `pin` will occur.

    Returns:
    --------
    array-like
        Interpolated data array of the same shape as `pout`.
    """

    # Initialize some things
    xflipped = False

    # Initialize the output array
    xout = np.full_like(pout, np.nan)

    # Determine whether to perform linear or logarithmic interpolation
    is_linear = abs(linlog) == 1

    # Reverse arrays if necessary to ensure they are in descending order
    if pin[0] < pin[-1]:
        pin = pin[::-1]
        xin = xin[::-1]
        xflipped = True

    if pout[0] < pout[-1]:
        pout = pout[::-1]

    # Remove missing data
    valid_mask = ~np.isnan(xin) & ~np.isnan(pin)
    pin = pin[valid_mask]
    xin = xin[valid_mask]

    if len(pin) < 2:
        raise ValueError("Not enough valid data points for interpolation")

    # Interpolation
    for np_idx in range(len(pout)):
        for nl in range(len(pin) - 1):
            if pin[nl] >= pout[np_idx] >= pin[nl + 1]:
                if is_linear:
                    slope = (xin[nl] - xin[nl + 1]) / (pin[nl] - pin[nl + 1])
                    xout[np_idx] = xin[nl + 1] + slope * (pout[np_idx] - pin[nl + 1])
                else:
                    slope = (xin[nl] - xin[nl + 1]) / (np.log(pin[nl]) - np.log(pin[nl + 1]))
                    xout[np_idx] = xin[nl + 1] + slope * (np.log(pout[np_idx]) - np.log(pin[nl + 1]))
                break

    # Extrapolation
    if linlog < 0:
        if is_linear:
            slope_high = (xin[1] - xin[0]) / (pin[1] - pin[0])
            slope_low = (xin[-1] - xin[-2]) / (pin[-1] - pin[-2])
        else:
            slope_high = (xin[1] - xin[0]) / (np.log(pin[1]) - np.log(pin[0]))
            slope_low = (xin[-1] - xin[-2]) / (np.log(pin[-1]) - np.log(pin[-2]))

        for np_idx in range(len(pout)):
            if pout[np_idx] > pin[0]:  # Above range, extrapolate
                if is_linear:
                    xout[np_idx] = xin[0] + slope_high * (pout[np_idx] - pin[0])
                else:
                    xout[np_idx] = xin[0] + slope_high * (np.log(pout[np_idx]) - np.log(pin[0]))
            elif pout[np_idx] < pin[-1]:  # Below range, extrapolate
                if is_linear:
                    xout[np_idx] = xin[-1] + slope_low * (pout[np_idx] - pin[-1])
                else:
                    xout[np_idx] = xin[-1] + slope_low * (np.log(pout[np_idx]) - np.log(pin[-1]))

    # If xin was flipped, we need to flip xout for consistency
    if xflipped:
        xout = xout[::-1]

    return xout




def int2p_n(pin, xin, pout, linlog, dim=-1):
    """
    Wrapper for int2p to handle multi-dimensional arrays along a specific dimension.

    Parameters:
    -----------
    pin : array-like
        Array of input pressure levels.

    xin : array-like
        Multi-dimensional array of data to be interpolated.

    pout : array-like
        Array of output pressure levels.

    linlog : int
        Integer indicating the type of interpolation:
        - abs(linlog) == 1 --> linear interpolation
        - abs(linlog) != 1 --> logarithmic interpolation
        - If linlog is negative, extrapolation outside the range of `pin` will occur.

    dim : int, optional (default=-1)
        The dimension index corresponding to pressure levels in `xin`.

    Returns:
    --------
    xout : array-like
        Interpolated data array with the same shape as `xin`, except along the `dim` dimension, which matches the size of `pout`.
    """

    # Move the specified dimension to the last position for easier iteration
    xin_moved = np.moveaxis(xin, dim, -1)

    # Prepare the output array with the new shape
    new_shape = list(xin_moved.shape)
    new_shape[-1] = len(pout)  # Update the last dimension to match pout
    xout = np.empty(new_shape, dtype=xin.dtype)

    # Iterate over all but the last dimension (which is now the pressure dimension)
    for index in np.ndindex(xin_moved.shape[:-1]):
        xout[index] = int2p(pin, xin_moved[index], pout, linlog)

    # Move the dimension back to its original position
    xout = np.moveaxis(xout, -1, dim)

    return xout





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

def pressure_to_hybrid(p_levels, data_on_p_levels, ps, a_coeff, b_coeff, level_dim=0, p0=100000, kflag=1):
    """
    Interpolate data on constant pressure levels to hybrid sigma levels.

    This function is a wrapper that prepares the input for the JIT-compiled version.
    """

    # Start timing
    start_time = time.time()

    # Check if p_levels is monotonically increasing
    if not np.all(np.diff(p_levels) > 0):
        print("p_levels is not monotonically increasing. Attempting to flip arrays...")

        # Flip the p_levels and data_on_p_levels arrays
        p_levels = p_levels[::-1]
        data_on_p_levels = np.flip(data_on_p_levels, axis=level_dim)

        # Recheck if p_levels is now monotonically increasing
        if not np.all(np.diff(p_levels) > 0):
            raise ValueError("p_levels must be monotonically increasing. Flipping arrays did not resolve the issue.")


    # Call the Python-equivalent Fortran function
    data_on_hybrid_levels = p2hyo(p_levels, data_on_p_levels, ps, a_coeff, b_coeff, p0, kflag)

    # Print elapsed time
    elapsed_time = time.time() - start_time
    print(f"Interpolation completed in {elapsed_time:.2f} seconds.")

    return data_on_hybrid_levels

# def pressure_to_hybrid(p_levels, data_on_p_levels, ps, a_coeff, b_coeff, p0=100000):
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
#     data_on_hybrid_levels = pressure_to_hybrid_numba(p_levels, data_on_p_levels, ps, hybrid_p)
#
#     # Print elapsed time
#     elapsed_time = time.time() - start_time
#     print(f"Interpolation completed in {elapsed_time:.2f} seconds.")
#
#     return data_on_hybrid_levels
#
#



###### NEEDS TO BE TESTED

@jit(nopython=True)
def ds2hbd(dati, sigi, sigo, intyp):
    """
    Python equivalent of the DS2HBD Fortran subroutine.

    Parameters:
    - dati: Input data on sigma coordinates (1D array)
    - sigi: Sigma levels (1D array)
    - sigo: Output sigma levels for hybrid coordinates (1D array)
    - intyp: Interpolation type (1: linear, 2: log, 3: log-log)

    Returns:
    - dato: Interpolated data on hybrid coordinates (1D array)
    """
    nli = sigi.size
    nlo = sigo.size
    dato = np.zeros(nlo)

    def a2ln(a1):
        return np.log(np.log(a1 + 1.001))

    for k in range(nlo):
        if sigo[k] <= sigi[1]:
            kp = 0
        elif sigo[k] >= sigi[nli-2]:
            kp = nli - 2
        else:
            kp = 0
            while sigo[k] > sigi[kp+1]:
                kp += 1

        if intyp == 1:  # Linear interpolation
            dato[k] = dati[kp] + (dati[kp+1] - dati[kp]) * (sigo[k] - sigi[kp]) / (sigi[kp+1] - sigi[kp])
        elif intyp == 2:  # Log interpolation
            dato[k] = dati[kp] + (dati[kp+1] - dati[kp]) * np.log(sigo[k] / sigi[kp]) / np.log(sigi[kp+1] / sigi[kp])
        elif intyp == 3:  # Log-log interpolation
            dato[k] = dati[kp] + (dati[kp+1] - dati[kp]) * (a2ln(sigo[k]) - a2ln(sigi[kp])) / (a2ln(sigi[kp+1]) - a2ln(sigi[kp]))

    return dato

@jit(nopython=True)
def dh2sdrv(dati, hya, hyb, p0, psfc, sigi, intyp):
    """
    Python equivalent of the DH2SDRV Fortran subroutine.

    Parameters:
    - dati: Input data on sigma levels (1D array)
    - hya: Hybrid A coefficients (1D array)
    - hyb: Hybrid B coefficients (1D array)
    - p0: Reference pressure
    - psfc: Surface pressure
    - sigi: Sigma levels (1D array)
    - intyp: Interpolation type (1: linear, 2: log, 3: log-log)

    Returns:
    - dato: Interpolated data on hybrid coordinates (1D array)
    """
    nlvi = sigi.size
    nlvo = hya.size
    sigo = np.zeros(nlvo)

    # Convert hybrid to sigma coordinates
    for n in range(nlvo):
        sigo[n] = (hya[n] * p0) / psfc + hyb[n]

    # Perform interpolation using ds2hbd
    dato = ds2hbd(dati, sigi, sigo, intyp)

    return dato

def sigma_to_hybrid(dati, sigma, hya, hyb, psfc, p0=100000, intyp=1):
    """
    Interpolate data from sigma levels to hybrid sigma levels.

    This function is a wrapper that prepares the input for the JIT-compiled version.
    """

    # Start timing
    start_time = time.time()

    # Check if sigma is monotonically increasing
    if not np.all(np.diff(sigma) > 0):
       raise ValueError("Sigma levels must be monotonically increasing")

    # Call the Python-equivalent Fortran function
    dato = dh2sdrv(dati, hya, hyb, p0, psfc, sigma, intyp)

    # Print elapsed time
    elapsed_time = time.time() - start_time
    print(f"Sigma to hybrid interpolation completed in {elapsed_time:.2f} seconds.")

    return dato
