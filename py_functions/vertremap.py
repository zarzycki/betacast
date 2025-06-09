import time
import logging
import warnings

import numpy as np
from numba import jit
from scipy import interpolate

logger = logging.getLogger(__name__)

_VERT_REMAP_VARS = [
    't', 'u', 'v', 'q', 'cldice', 'cldliq',
    'z', 'theta', 'rho', 'w', 'o3'
]

__pres_lev_mandatory__ = np.array([
    1000, 925, 850, 700, 500, 400, 300, 250, 200, 150, 100, 70, 50, 30, 20, 10,
    7, 5, 3, 2, 1
]).astype(np.float64) * 100.0  # Convert mb to Pa

# 7/2025, CMZ: These are nominal levels taken from a HRRRv3 wrfnat HRRR file
__pres_lev_hrrr__ = np.array([
    101500.0, 100800.0, 99900.0, 98500.0, 96710.0, 94540.0, 92170.0, 89610.0, 86750.0,
    83590.0, 80140.0, 76380.0, 72330.0, 67980.0, 63340.0, 58450.0, 53710.0,
    49450.0, 45640.0, 42210.0, 39110.0, 36320.0, 33820.0, 31560.0, 29530.0,
    27700.0, 26000.0, 24472.0, 23069.0, 21666.0, 20263.0, 18860.0, 17457.0,
    16054.0, 14651.0, 13248.0, 11845.0, 10443.0, 9139.0, 8094.0, 7244.0,
    6453.0, 5711.0, 5025.0, 4388.0, 3791.0, 3233.8, 2716.4, 2500.0
]).astype(np.float64)


# These HRRRv1 and v2 eta levels, HRRRv3 switched to hybrid coord with lid ~10mb
def hrrr_eta_to_plev(nominal_surface_pressure_Pa=1000.0):
    eta_levels = [1.0000, 0.9980, 0.9940, 0.9870, 0.9750, 0.9590, 0.9390, 0.9160, 0.8920, 0.8650,
                  0.8350, 0.8020, 0.7660, 0.7270, 0.6850, 0.6400, 0.5920, 0.5420, 0.4970, 0.4565,
                  0.4205, 0.3877, 0.3582, 0.3317, 0.3078, 0.2863, 0.2670, 0.2496, 0.2329, 0.2188,
                  0.2047, 0.1906, 0.1765, 0.1624, 0.1483, 0.1342, 0.1201, 0.1060, 0.0919, 0.0778,
                  0.0657, 0.0568, 0.0486, 0.0409, 0.0337, 0.0271, 0.0209, 0.0151, 0.0097, 0.0047, 0.0000]
    p_interfaces = np.array(eta_levels) * nominal_surface_pressure_Pa * 100.
    p_midpoints = 0.5 * (p_interfaces[:-1] + p_interfaces[1:])
    return p_midpoints


@jit(nopython=True)
def int2p(pin, xin, pout, linlog):
    """
    Interpolates data from one set of pressure levels to another.

    Parameters:
    -----------
    pin : array-like
        1-D array of input pressure levels.

    xin : array-like
        1-D array of data to be interpolated. Should have the same length as `pin`.

    pout : array-like
        1-D array of output pressure levels.

    linlog : int
        Integer indicating the type of interpolation:
        - abs(linlog) == 1 --> linear interpolation
        - abs(linlog) != 1 --> logarithmic interpolation
        - If linlog is negative, extrapolation outside the range of `pin` will occur.

    Returns:
    --------
    array-like
        Interpolated data array of the same length as `pout`.
    """

    # Initialize needed bools
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
        Array of input pressure levels. Can be:
        - 1D array: same pressure levels for all spatial points
        - Multi-dimensional: different pressure levels for each spatial point
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
        Interpolated data array with the same shape as `xin`, except along the `dim` dimension.
    """

    # Determine if pin is 1D or multi-dimensional
    pin_is_1d = (pin.ndim == 1)

    # Move the level dimension to the last position for easier iteration
    xin_moved = np.moveaxis(xin, dim, -1)

    # If a 3D pressure field (for example) move lev to the end
    if not pin_is_1d:
        pin_moved = np.moveaxis(pin, dim, -1)

    # Prepare the output array
    new_shape = list(xin_moved.shape)
    new_shape[-1] = len(pout)
    xout = np.empty(new_shape, dtype=xin.dtype)

    # Perform interpolation for each spatial point
    for index in np.ndindex(xin_moved.shape[:-1]):
        if pin_is_1d:
            # Use the same 1D pressure profile for all spatial points
            xout[index] = int2p(pin, xin_moved[index], pout, linlog)
        else:
            # Use spatially-varying pressure profiles
            xout[index] = int2p(pin_moved[index], xin_moved[index], pout, linlog)

    # Move the level dimension back to its original position
    xout = np.moveaxis(xout, -1, dim)

    return xout


def int2p_n_extrap(pin, xin, pout, linlog, dim=-1, flip_p=False,
                     extrapolate=False, variable=None, ps=None, t_bot=None, phi_sfc=None):
    """
    Enhanced wrapper for int2p with extrapolation options.

    Parameters:
    -----------
    pin : array-like
        Array of input pressure levels. Can be:
        - 1D array: same pressure levels for all spatial points
        - Multi-dimensional: different pressure levels for each spatial point
        - Should be organized top-to-bottom (see flip_p)
    xin : array-like
        Multi-dimensional array of data to be interpolated.
    pout : array-like
        Array of output pressure levels.
    linlog : int
        Integer indicating the type of interpolation.
    dim : int, optional (default=-1)
        The dimension index corresponding to pressure levels.
    flip_p : bool, optional (default=False)
        If True, flip arrays before interpolation and flip output back.
        This is needed if p dimension is natively bottom->top
    extrapolate : bool, optional (default=False)
        If True, use sophisticated extrapolation for below-ground points.
    variable : str, optional
        Variable type for sophisticated extrapolation:
        - 'temperature': Use ECMWF temperature extrapolation
        - 'geopotential': Use ECMWF geopotential extrapolation
        - 'other' or None: Use simple persistence
    ps : array-like, optional
        Surface pressure (2D array). Required if extrapolate=True.
    t_bot : array-like, optional
        Bottom-level temperature. Required for geopotential extrapolation.
    phi_sfc : array-like, optional
        Surface geopotential. Required for temperature extrapolation.

    Returns:
    --------
    xout : array-like
        Interpolated data array with sophisticated extrapolation applied.
    """

    # First do the basic interpolation
    xout = int2p_n(pin, xin, pout, linlog, dim)

    # If no extrapolation requested, return basic result
    if not extrapolate:
        return xout

    # Validate inputs for extrapolation
    if ps is None:
        raise ValueError("Surface pressure (ps) required for extrapolation")
    if variable == 'temperature' and phi_sfc is None:
        raise ValueError("Surface geopotential (phi_sfc) required for temperature extrapolation")
    if variable == 'geopotential' and (t_bot is None or phi_sfc is None):
        raise ValueError("Both t_bot and phi_sfc required for geopotential extrapolation")

    # Sophisticated extrapolation only works with multi-dimensional pressure fields
    if pin.ndim == 1:
        print("Warning: Sophisticated extrapolation requires 3D pressure field. Using basic result.")
        return xout

    # Prepare data for extrapolation (handle flipping if needed)
    pin_for_extrap = pin
    xin_for_extrap = xin
    pout_for_extrap = pout
    xout_for_extrap = xout

    # Flip everything if needed
    if flip_p:
        pin_for_extrap = np.flip(pin, axis=dim)
        xin_for_extrap = np.flip(xin, axis=dim)
        pout_for_extrap = pout[::-1]
        xout_for_extrap = np.flip(xout, axis=dim)

    # Apply extrapolation
    xout_extrap = _vertical_remap_extrap(
        new_levels=pout_for_extrap,
        lev_dim=dim,
        data=xin_for_extrap,
        output=xout_for_extrap,
        pressure=pin_for_extrap,
        ps=ps,
        variable=variable or 'other',
        t_bot=t_bot,
        phi_sfc=phi_sfc
    )

    # Apply final flipping if needed
    if flip_p:
        xout_extrap = np.flip(xout_extrap, axis=dim)

    return xout_extrap


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


def pressure_to_hybrid(p_levels, data_on_p_levels, ps, a_coeff, b_coeff, level_dim=0, p0=100000, kflag=1, name=None):
    """
    Interpolate data on constant pressure levels to hybrid sigma levels.

    This function is a wrapper that prepares the input for the JIT-compiled version.
    """

    # Start timing
    start_time = time.time()

    # Check if p_levels is monotonically increasing
    if not np.all(np.diff(p_levels) > 0):
        logging.info("p_levels is not monotonically increasing. Attempting to flip arrays...")

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
    if name:
        logging.info(f"Interpolation for {name} completed in {elapsed_time:.2f} seconds.")
    else:
        logging.info(f"Interpolation completed in {elapsed_time:.2f} seconds.")

    return data_on_hybrid_levels


def pres2hyb_all(data_vars, ps, hya, hyb):
    """
    Interpolate data on constant pressure levels to hybrid sigma levels.

    This function is a wrapper that prepares the input for the JIT-compiled version.
    """

    data_vint = {}

    allowable_interp_vars = _VERT_REMAP_VARS

    # Create non-interpolated vars
    data_vint['hya'] = hya
    data_vint['hyb'] = hyb

    # Loop ovr the keys in data_in and interpolate if the key is in allowable_interp_vars
    for key in data_vars:
        if key in allowable_interp_vars:
            data_vint[key] = pressure_to_hybrid(data_vars['lev'], data_vars[key], ps, data_vint['hya'], data_vint['hyb'], name=key)
        else:
            data_vint[key] = data_vars[key]

    return data_vint


def interp_hybrid_to_pressure_wrapper(data_vars, ps, hyam, hybm, new_levels, lev_dim=0, method='log', extrapolate=True):

    allowable_interp_vars = _VERT_REMAP_VARS

    for var_name, data in data_vars.items():
        if isinstance(data, np.ndarray) and data.ndim == 3:
            variable_type = 'temperature' if var_name == 't' else 'geopotential' if var_name == 'z' else 'other'
            logging.info(f"interp_hybrid_to_pressure for variable {var_name}, using {variable_type}")
            data_vars[var_name] = interp_hybrid_to_pressure(
                data=data,
                ps=ps,
                hyam=hyam,
                hybm=hybm,
                new_levels=new_levels,
                lev_dim=lev_dim,
                method=method,
                extrapolate=extrapolate,
                variable=variable_type,
                t_bot=data_vars['ts'],
                phi_sfc=data_vars['phis']
            )
    return data_vars


def _temp_extrapolate(data, lev_dim, lev, p_sfc, ps, phi_sfc):
    """Extrapolates temperature below ground using the ECMWF formulation."""
    R_d = 287.04  # dry air gas constant
    g_inv = 1 / 9.80616  # inverse of gravity
    alpha = 0.0065 * R_d * g_inv

    tstar = data.take(indices=-1, axis=lev_dim) * (1 + alpha * (ps / p_sfc - 1))
    hgt = phi_sfc * g_inv
    t0 = tstar + 0.0065 * hgt
    tplat = np.minimum(298, t0)

    tprime0 = np.where((2000 <= hgt) & (hgt <= 2500),
                       0.002 * ((2500 - hgt) * t0 + (hgt - 2000) * tplat),
                       np.nan)
    tprime0 = np.where(2500 < hgt, tplat, tprime0)

    alnp = np.where(hgt < 2000, alpha * np.log(lev / ps),
                    R_d * (tprime0 - tstar) / phi_sfc * np.log(lev / ps))
    alnp = np.where(tprime0 < tstar, 0, alnp)

    result = tstar * (1 + alnp + 0.5 * (alnp**2) + (1 / 6) * (alnp**3))

    return result


def _geo_height_extrapolate(t_bot, lev, p_sfc, ps, phi_sfc):
    """Extrapolates geopotential height below ground using the ECMWF formulation."""
    R_d = 287.04  # dry air gas constant
    g_inv = 1 / 9.80616  # inverse of gravity
    alpha = 0.0065 * R_d * g_inv

    tstar = t_bot * (1 + alpha * (ps / p_sfc - 1))
    hgt = phi_sfc * g_inv
    t0 = tstar + 0.0065 * hgt

    epsilon = 1e-12
    # Added epsilon so we don't divide by zero
    alph = np.where((tstar <= 290.5) & (t0 > 290.5),
                    R_d / (phi_sfc + epsilon) * (290.5 - tstar), alpha)

    alph = np.where((tstar > 290.5) & (t0 > 290.5), 0, alph)
    tstar = np.where((tstar > 290.5) & (t0 > 290.5), 0.5 * (290.5 + tstar), tstar)
    tstar = np.where(tstar < 255, 0.5 * (tstar + 255), tstar)

    alnp = alph * np.log(lev / ps)
    return hgt - R_d * tstar * g_inv * np.log(lev / ps) * (1 + 0.5 * alnp + (1 / 6) * alnp**2)


def _linear_interpolate_1d(new_levels, p_levels, data_on_p_levels, axis=0, fill_value=np.nan):
    """Custom linear interpolation along a specified axis using a linear x-scale."""

    # Ensure the axis of interpolation is the first axis for easier manipulation
    p_levels = np.moveaxis(p_levels, axis, 0)
    data_on_p_levels = np.moveaxis(data_on_p_levels, axis, 0)
    new_levels = np.moveaxis(new_levels, axis, 0)

    # Initialize output array with the shape of data_on_p_levels, but with the new_levels size
    interp_shape = list(data_on_p_levels.shape)
    interp_shape[0] = len(new_levels)  # Update size of the interpolation axis to match new_levels
    interp_data = np.full(interp_shape, fill_value, dtype=data_on_p_levels.dtype)

    # Perform interpolation for each slice along the axis
    for idx in np.ndindex(data_on_p_levels.shape[1:]):
        # Extract 1D slices for interpolation
        p_slice = p_levels[(slice(None),) + idx]
        dp_slice = data_on_p_levels[(slice(None),) + idx]

        # Interpolate over the 1D slice
        interp_slice = np.interp(new_levels, p_slice, dp_slice, left=fill_value, right=fill_value)

        # Assign the interpolated slice back to the corresponding location in interp_data
        interp_data[(slice(None),) + idx] = interp_slice

    # Move the interpolated data back to its original axis
    return np.moveaxis(interp_data, 0, axis)


def _log_interpolate_1d(new_levels, p_levels, data_on_p_levels, axis=0, fill_value=np.nan):
    """Custom log-linear interpolation along a specified axis using a logarithmic x-scale."""

    # Convert pressure levels and new levels to their logarithmic scales
    log_new_levels = np.log(new_levels)
    log_p_levels = np.log(p_levels)

    # Initialize output array with the shape of data_on_p_levels, but with the new_levels size
    interp_shape = list(data_on_p_levels.shape)
    interp_shape[0] = len(new_levels)  # Update size of the interpolation axis to match new_levels
    interp_data = np.full(interp_shape, fill_value, dtype=data_on_p_levels.dtype)

    # Perform interpolation for each slice along the axis
    for idx in np.ndindex(data_on_p_levels.shape[1:]):
        # Extract 1D slices for interpolation
        lp_slice = log_p_levels[(slice(None),) + idx]
        dp_slice = data_on_p_levels[(slice(None),) + idx]

        #logging.info("------------")
        #logging.info(np.exp(log_new_levels))
        #logging.info(np.exp(lp_slice))
        #logging.info(dp_slice)
        # Interpolate over the 1D slice
        interp_slice = np.interp(log_new_levels, lp_slice, dp_slice, left=-999.9, right=fill_value)
        #logging.info(interp_slice)
        if np.any(interp_slice == -999.9):
            # Find the indices of the first two valid values (not equal to -999.9)
            valid_indices = np.where(interp_slice != -999.9)[0][:2]

            # Get the first two valid values
            first_two_valid_values = interp_slice[valid_indices]

            # Calculate the corresponding log pressure levels for the first two valid values
            lp_valid = log_new_levels[valid_indices]

            # Calculate the derivative (slope) between the first two valid points
            slope = (first_two_valid_values[1] - first_two_valid_values[0]) / (lp_valid[1] - lp_valid[0])

            # Extrapolate the -999.9 values using the slope
            for i in np.where(interp_slice == -999.9)[0]:
                interp_slice[i] = first_two_valid_values[0] + slope * (log_new_levels[i] - lp_valid[0])
        #logging.info(interp_slice)
        # Assign the interpolated slice back to the corresponding location in interp_data
        interp_data[(slice(None),) + idx] = interp_slice

    # Move the interpolated data back to its original axis
    return interp_data


def _func_interpolate(method='linear'):
    """Define custom interpolation function."""
    if method == 'linear':
        return _linear_interpolate_1d
    elif method == 'log':
        return _log_interpolate_1d
    else:
        raise ValueError(f'Unknown interpolation method: {method}. Supported methods are: "log" and "linear".')


def _pressure_from_hybrid(psfc, hya, hyb, p0=100000.):
    """Calculate pressure at the hybrid levels."""
    # Dynamically expand hya and hyb to match the shape of psfc
    expanded_shape = (len(hya),) + (1,) * (psfc.ndim)
    hya_expanded = hya.reshape(expanded_shape)
    hyb_expanded = hyb.reshape(expanded_shape)
    return hya_expanded * p0 + hyb_expanded * psfc


def _vertical_remap(func_interpolate, new_levels, xcoords, data, interp_axis=0):
    """Execute the defined interpolation function on data."""
    with warnings.catch_warnings():
        warnings.filterwarnings(
            "ignore", r"Interpolation point out of data bounds encountered")
        return func_interpolate(new_levels, xcoords, data, axis=interp_axis)


def _vertical_remap_extrap(new_levels, lev_dim, data, output, pressure, ps, variable, t_bot, phi_sfc):
    """A helper function to call the appropriate extrapolation function based on the user's inputs."""
    ## NEEDS TO BE TOP-TO-BOTTOM (i.e., small pressures in the 0th index, large numbers in -1!

    # Get indices for the top and bottom of the vertical dimension
    bottom_index = -1

    bottom_level = np.take_along_axis(pressure, np.expand_dims(np.argmax(pressure, axis=lev_dim), axis=lev_dim), axis=lev_dim).squeeze()
    bottom_value = data.take(indices=bottom_index, axis=lev_dim)

    # Prepare new_levels for broadcasting
    new_levels_broadcasted = np.expand_dims(new_levels, axis=tuple(range(1, pressure.ndim)))

    # Downward Extrapolation
    if variable == 'temperature':
        extrapolated = _temp_extrapolate(data, lev_dim, new_levels_broadcasted, bottom_level, ps, phi_sfc)
    elif variable == 'geopotential':
        extrapolated = _geo_height_extrapolate(t_bot, new_levels_broadcasted, bottom_level, ps, phi_sfc)
    else:
        extrapolated = bottom_value

    # Apply downward extrapolation
    output = np.where(new_levels_broadcasted > (bottom_level - 0.1), extrapolated, output)

    # Debugging: Count and print number of NaNs after downward extrapolation
    nan_count_downward = np.isnan(output).sum()
    if nan_count_downward > 0:
        logging.info(f"NaNs detected after downward extrapolation: {nan_count_downward}")
        nan_indices_downward = np.argwhere(np.isnan(output))
        logging.info("NaNs found at indices (downward extrapolation):")
        logging.info(nan_indices_downward)
        logging.info("Inspecting inputs at NaN indices for downward extrapolation:")
        for idx in nan_indices_downward:
            logging.info(f"Index {idx}:")

            # Ensure the index is within bounds for the data array
            if idx[0] < data.shape[0] and all(i < s for i, s in zip(idx[1:], data.shape[1:])):
                logging.info(f"data at index: {data[tuple(idx)]}")
            else:
                logging.info("Data at index: Index out of bounds for data array.")

            # Ensure the index is within bounds for the bottom_level, ps, and phi_sfc arrays
            if all(i < s for i, s in zip(idx[1:], bottom_level.shape)):
                logging.info(f"bottom_level at index: {bottom_level[tuple(idx[1:])]}")
                logging.info(f"ps at index: {ps[tuple(idx[1:])]}")
                logging.info(f"phi_sfc at index: {phi_sfc[tuple(idx[1:])]}")
            else:
                logging.info("Bottom_level, ps, or phi_sfc index out of bounds.")

            # Print the corresponding new_levels value
            if idx[0] < len(new_levels):
                logging.info(f"new_levels at level index: {new_levels[idx[0]]}")
            else:
                logging.info("new_levels index out of bounds.")

    return output


def interp_hybrid_to_pressure(data: np.ndarray,
                              ps: np.ndarray,
                              hyam: np.ndarray,
                              hybm: np.ndarray,
                              p0: float = 100000.,
                              new_levels: np.ndarray = __pres_lev_mandatory__,
                              lev_dim: int = 0,
                              method: str = 'linear',
                              extrapolate: bool = False,
                              variable: str = None,
                              t_bot: np.ndarray = None,
                              phi_sfc: np.ndarray = None) -> np.ndarray:
    """Interpolate and extrapolate data from hybrid-sigma levels to isobaric levels."""

    # Check inputs
    if extrapolate and variable is None:
        raise ValueError("If extrapolate is True, variable must be provided.")
    if variable in ['geopotential', 'temperature'] and (t_bot is None or phi_sfc is None):
        raise ValueError("If variable is 'geopotential' or 'temperature', both t_bot and phi_sfc must be provided")
    if variable not in ['geopotential', 'temperature', 'other', None]:
        raise ValueError(f"The value of variable is {variable}, but the accepted values are 'temperature', 'geopotential', 'other', or None.")

    func_interpolate = _func_interpolate(method)

    # Calculate pressure levels at the hybrid levels
    pressure = _pressure_from_hybrid(ps, hyam, hybm, p0)  # Pa

    # Prepare output data structure
    out_shape = list(data.shape)
    out_shape[lev_dim] = new_levels.size
    output = np.empty(out_shape, dtype=data.dtype)

    # Perform the vertical interpolation
    output = _vertical_remap(
        func_interpolate, new_levels, pressure, data, interp_axis=lev_dim)

    if extrapolate:
        output = _vertical_remap_extrap(new_levels, lev_dim, data, output, pressure, ps, variable, t_bot, phi_sfc)

    return output


def interpolate_mpas_columns_wrapper(mpas_data, data_horiz):

    t_wrf = np.zeros((mpas_data['nlev'], mpas_data['ncell']))
    theta_wrf = np.zeros((mpas_data['nlev'], mpas_data['ncell']))
    rho_wrf = np.zeros((mpas_data['nlev'], mpas_data['ncell']))
    w_wrf = np.zeros((mpas_data['nlevi'], mpas_data['ncell']))
    q_wrf = np.zeros((mpas_data['nlev'], mpas_data['ncell']))
    u_wrf = np.zeros((mpas_data['nlev'], mpas_data['ncell']))
    v_wrf = np.zeros((mpas_data['nlev'], mpas_data['ncell']))

    for ix in range(mpas_data['ncell']):
        if ix % 10000 == 0:
            logging.info(f"... MPAS_VERT_INTERP: {100. * ix / (mpas_data['ncell'] - 1):.2f}%")

        if ix == 0:
            z_src_orientation = 'bottom_to_top' if data_horiz['z'][-1, ix] > data_horiz['z'][0, ix] else 'top_to_bottom'
            z_mpas_orientation = 'bottom_to_top' if mpas_data['z'][ix, -1] > mpas_data['z'][ix, 0] else 'top_to_bottom'

            logger.debug(f"z_src_orientation is {z_src_orientation}")
            logger.debug(f"z_mpas_orientation is {z_mpas_orientation}")

            # Check if the orientations are the same or different
            if z_src_orientation == z_mpas_orientation:
                logger.debug("z_fv and mpas_z have the same orientation")
                flip_array = False
            else:
                logger.debug("z_fv and mpas_z have different orientations")
                logger.debug("we will flip arrays in interpolate_mpas_columns_wrapper")
                flip_array = True

        # We want to priorize the mpas grid orientation, so flip the src arrays if z
        # orientation differs between the two
        if flip_array:
            t_wrf[:, ix], theta_wrf[:, ix], rho_wrf[:, ix], w_wrf[:, ix], \
            q_wrf[:, ix], u_wrf[:, ix], v_wrf[:, ix] = interpolate_single_mpas_column(
                ix,
                mpas_data['nlev'], mpas_data['nlevi'], mpas_data['z'],
                data_horiz['t'][::-1, ix], data_horiz['z'][::-1, ix],
                data_horiz['theta'][::-1, ix], data_horiz['rho'][::-1, ix],
                data_horiz['w'][::-1, ix], data_horiz['q'][::-1, ix],
                data_horiz['u'][::-1, ix], data_horiz['v'][::-1, ix]
            )
        else:
            t_wrf[:, ix], theta_wrf[:, ix], rho_wrf[:, ix], w_wrf[:, ix], \
            q_wrf[:, ix], u_wrf[:, ix], v_wrf[:, ix] = interpolate_single_mpas_column(
                ix,
                mpas_data['nlev'], mpas_data['nlevi'], mpas_data['z'],
                data_horiz['t'][:, ix], data_horiz['z'][:, ix],
                data_horiz['theta'][:, ix], data_horiz['rho'][:, ix],
                data_horiz['w'][:, ix], data_horiz['q'][:, ix],
                data_horiz['u'][:, ix], data_horiz['v'][:, ix]
            )

    data_horiz['t'] = t_wrf
    data_horiz['theta'] = theta_wrf
    data_horiz['rho'] = rho_wrf
    data_horiz['w'] = w_wrf
    data_horiz['q'] = q_wrf
    data_horiz['u'] = u_wrf
    data_horiz['v'] = v_wrf

    return data_horiz


def standard_atmosphere(z):
    """
    Generate a US Standard Atmosphere 1976 profile for given heights.

    Parameters:
    z (array): Heights in meters

    Returns:
    tuple: (temperature, pressure, density)
    """
    # Constants
    g0 = 9.80665  # m/s^2
    R = 287.053  # J/(kg*K)
    M = 0.0289644  # kg/mole (mean molecular mass of air)

    # Layer definitions
    h = np.array([0, 11000, 20000, 32000, 47000, 51000, 71000, 84852])  # m
    P = np.array([101325, 22632.1, 5474.89, 868.019, 110.906, 66.9389, 3.95642, 0.3734])  # Pa
    T = np.array([288.15, 216.65, 216.65, 228.65, 270.65, 270.65, 214.65, 186.95])  # K
    L = np.array([-6.5, 0, 1.0, 2.8, 0, -2.8, -2.0, 0]) / 1000  # K/m

    # Find the appropriate layer for each altitude
    layer = np.searchsorted(h, z) - 1
    layer = np.clip(layer, 0, len(h) - 2)

    # Calculate temperature
    dh = z - h[layer]
    T_z = T[layer] + L[layer] * dh

    # Calculate pressure
    P_z = np.zeros_like(z)
    for i in range(len(h) - 1):
        mask = layer == i
        if L[i] != 0:  # Non-isothermal layer
            P_z[mask] = P[i] * (T_z[mask] / T[i]) ** (-g0 / (R * L[i]))
        else:  # Isothermal layer
            P_z[mask] = P[i] * np.exp(-g0 * dh[mask] / (R * T[i]))

    # Calculate density
    rho_z = P_z / (R * T_z)

    return T_z, P_z, rho_z



def taper_to_standard(valid_data, valid_z, full_z, std_data, taper_rate=1.0):
    """
    Taper between valid data and standard atmosphere data with controllable taper rate.

    Parameters:
    valid_data (array): Valid data points
    valid_z (array): Heights corresponding to valid data
    full_z (array): Full height array
    std_data (array): Standard atmosphere data for full height array
    taper_rate (float): Controls the rate of tapering. Default is 1.0.
                        Higher values make the taper more abrupt,
                        lower values make it more gradual.

    Returns:
    array: Tapered data
    """
    nan_mask = np.isnan(valid_data)

    if np.all(nan_mask):
        return std_data
    if not np.any(nan_mask):
        return valid_data

    first_valid = np.argmin(nan_mask)
    last_valid = len(nan_mask) - np.argmin(nan_mask[::-1]) - 1

    tapered_data = np.copy(valid_data)

    def taper_weight(x):
        return np.clip(x ** taper_rate, 0, 1)

    if first_valid > 0:
        w = taper_weight((full_z[:first_valid] - full_z[0]) / (full_z[first_valid] - full_z[0]))
        tapered_data[:first_valid] = w * valid_data[first_valid] + (1 - w) * std_data[:first_valid]

    if last_valid < len(full_z) - 1:
        w = taper_weight((full_z[last_valid+1:] - full_z[last_valid]) / (full_z[-1] - full_z[last_valid]))
        tapered_data[last_valid+1:] = (1 - w) * valid_data[last_valid] + w * std_data[last_valid+1:]

    tapered_data[np.isnan(tapered_data)] = std_data[np.isnan(tapered_data)]

    return tapered_data






def interpolate_single_mpas_column(ix, mpas_nlev, mpas_nlevi, mpas_z, t_fv, z_fv, theta_fv, rho_fv, w_fv, q_fv, u_fv, v_fv):

    zmid = (mpas_z[ix, 1:mpas_nlevi] + mpas_z[ix, 0:mpas_nlevi - 1]) / 2.0
    zint = mpas_z[ix, :]

    t_wrf_col = np.zeros(mpas_nlev)
    theta_wrf_col = np.zeros(mpas_nlev)
    rho_wrf_col = np.zeros(mpas_nlev)
    w_wrf_col = np.zeros(mpas_nlevi)
    q_wrf_col = np.zeros(mpas_nlev)
    u_wrf_col = np.zeros(mpas_nlev)
    v_wrf_col = np.zeros(mpas_nlev)

    if ix == 0:
        logger.debug("interpolate_single_mpas_column T metrics at ix = 0:")
        logger.debug(f"z_fv[0]: {z_fv[0]}, z_fv[nlev-1]: {z_fv[-1]}")
        logger.debug(f"zmid[0]: {zmid[0]}, zmid[nlev-1]: {zmid[-1]}")
        logger.debug(f"t_fv[0]: {t_fv[0]}, t_fv[nlev-1]: {t_fv[-1]}")
        logger.debug(f"t_fv range: {t_fv.min()} to {t_fv.max()}")
        logger.debug(f"t_fv size: {t_fv.size}, z_fv size: {z_fv.size}, zmid size: {zmid.size}")

    t_wrf_col = z_to_z_interp_wrapper(t_fv, z_fv, zmid, extrap_low="linear", extrap_high="linear")
    theta_wrf_col = z_to_z_interp_wrapper(theta_fv, z_fv, zmid, extrap_low="linear", extrap_high="linear")
    rho_wrf_col = z_to_z_interp_wrapper(rho_fv, z_fv, zmid, extrap_low="linear", extrap_high="log")
    w_wrf_col = z_to_z_interp_wrapper(w_fv, z_fv, zint, extrap_low="linear", extrap_high="fade")
    q_wrf_col = z_to_z_interp_wrapper(q_fv, z_fv, zmid, extrap_low="linear", extrap_high="linear")

    # u and v we don't extrapolate aloft to prevent wind speeds from getting too high
    u_wrf_col = z_to_z_interp_wrapper(u_fv, z_fv, zmid, extrap_low="linear", extrap_high="persist")
    v_wrf_col = z_to_z_interp_wrapper(v_fv, z_fv, zmid, extrap_low="linear", extrap_high="persist")

    # Get standard atmosphere for this column
    #std_T, std_P, std_rho = standard_atmosphere(zmid)
    #std_theta = std_T * (100000 / std_P) ** (2/7)  # Potential temperature
    #t_wrf_col = taper_to_standard(t_wrf_col, zmid, zmid, std_T, taper_rate=2.5)
    #theta_wrf_col = taper_to_standard(theta_wrf_col, zmid, zmid, std_theta, taper_rate=2.5)
    #rho_wrf_col = taper_to_standard(rho_wrf_col, zmid, zmid, std_rho, taper_rate=2.5)

    # Initialize all adjustment counters to 0
    adjustments = {var + '_adjustments': 0 for var in ['t_pos', 'theta_pos', 'rho_pos', 'rho_mono', 'q_pos']}
    t_wrf_col, adjustments['t_pos_adjustments'] = enforce_positive_values(t_wrf_col, replacement_value=180.0)
    theta_wrf_col, adjustments['theta_pos_adjustments'] = enforce_positive_values(theta_wrf_col, replacement_value=500.0)
    rho_wrf_col, adjustments['rho_pos_adjustments'] = enforce_positive_values(rho_wrf_col, replacement_value=1.0e-7)
    # rho_wrf_col, adjustments['rho_mono_adjustments'] = enforce_monotonic_values(rho_wrf_col)
    q_wrf_col, adjustments['q_pos_adjustments'] = enforce_positive_values(q_wrf_col, replacement_value=0.0)
    if ix == 0 or any(adjustments.values()):
        logger.debug(f"Column {ix} adjustments:")
        for adjustment_type, count in adjustments.items():
            logger.debug(f"  {adjustment_type}: {count}")

    return t_wrf_col, theta_wrf_col, rho_wrf_col, w_wrf_col, q_wrf_col, u_wrf_col, v_wrf_col


@jit(nopython=True)
def z_to_z_interp(data_src, z_src, z_target, extrap_low, extrap_high):
    nlev = len(z_target)
    data_target = np.full(nlev, np.nan, dtype=data_src.dtype)

    # Interpolation
    for i in range(nlev):
        if z_src[0] <= z_target[i] <= z_src[-1]:
            for j in range(len(z_src) - 1):
                if z_src[j] <= z_target[i] <= z_src[j + 1]:
                    data_target[i] = data_src[j] + (data_src[j + 1] - data_src[j]) * (z_target[i] - z_src[j]) / (z_src[j + 1] - z_src[j])
                    break

    # Extrapolation
    ixvalid = np.where(~np.isnan(data_target))[0]
    if len(ixvalid) > 0:
        lowest_valid, highest_valid = ixvalid[0], ixvalid[-1]

        # Low extrapolation
        if extrap_low == -90:  # nan/None
            pass
        elif extrap_low == -91:  # Linear
            deriv = (data_target[lowest_valid + 1] - data_target[lowest_valid]) / (z_target[lowest_valid + 1] - z_target[lowest_valid])
            for i in range(lowest_valid):
                data_target[i] = data_target[lowest_valid] + deriv * (z_target[i] - z_target[lowest_valid])
        elif extrap_low == -92:  # Log
            log_ratio = np.log(data_target[lowest_valid + 1] / data_target[lowest_valid]) / (z_target[lowest_valid + 1] - z_target[lowest_valid])
            for i in range(lowest_valid):
                data_target[i] = data_target[lowest_valid] * np.exp(log_ratio * (z_target[i] - z_target[lowest_valid]))
        elif extrap_low == -93:  # Persist
            data_target[:lowest_valid] = data_target[lowest_valid]
        elif extrap_low == -94:  # Fade
            for i in range(lowest_valid):
                w = i / lowest_valid
                data_target[i] = w * data_target[lowest_valid]
        else:  # Constant value
            data_target[:lowest_valid] = extrap_low

        # High extrapolation
        if extrap_high == -90:  # nan/None
            pass
        elif extrap_high == -91:  # Linear
            deriv = (data_target[highest_valid] - data_target[highest_valid - 1]) / (z_target[highest_valid] - z_target[highest_valid - 1])
            for i in range(highest_valid + 1, nlev):
                data_target[i] = data_target[highest_valid] + deriv * (z_target[i] - z_target[highest_valid])
        elif extrap_high == -92:  # Log
            log_ratio = np.log(data_target[highest_valid] / data_target[highest_valid - 1]) / (z_target[highest_valid] - z_target[highest_valid - 1])
            for i in range(highest_valid + 1, nlev):
                data_target[i] = data_target[highest_valid] * np.exp(log_ratio * (z_target[i] - z_target[highest_valid]))
        elif extrap_high == -93:  # Persist
            data_target[highest_valid + 1:] = data_target[highest_valid]
        elif extrap_high == -94:  # Fade
            for i in range(highest_valid + 1, nlev):
                w = (nlev - i - 1) / (nlev - highest_valid - 1)
                data_target[i] = w * data_target[highest_valid]
        else:  # Constant value
            data_target[highest_valid + 1:] = extrap_high

    return data_target


# WARNING: Using -90, -91, -92, -93, -94 as "protected" integers. If you set the extrap option to one of these values it will not behave correctly!
@jit(nopython=True)
def get_extrap_option(option):
    if isinstance(option, (int, float)):
        return int(option) if option in [0, 1, 2, 3] else float(option)
    elif isinstance(option, str):
        option_lower = option.lower()
        if option_lower == "nan":
            return -90
        elif option_lower == "linear":
            return -91
        elif option_lower == "log":
            return -92
        elif option_lower == "persist":
            return -93
        elif option_lower == "fade":
            return -94
    return -90  # Default to nan if invalid option


@jit(nopython=True)
def enforce_positive_values(data, eps=1e-8, replacement_value=None):
    """
    Enforce positive values in the data array.

    Parameters:
    -----------
    data : array-like
        Input data array.
    eps : float, optional (default 1e-8)
        Minimum positive value to enforce if replacement_value is None.
    replacement_value : float, optional (default None)
        Value to use for replacing negative values. If None, use max(data[i], eps).

    Returns:
    --------
    tuple: (modified_data, num_replacements)
        modified_data: Array with enforced positive values.
        num_replacements: Number of values that were replaced.
    """
    modified_data = data.copy()
    num_replacements = 0

    for i in range(len(data)):
        if data[i] <= 0:
            if replacement_value is not None:
                modified_data[i] = replacement_value
            else:
                modified_data[i] = max(data[i], eps)
            num_replacements += 1

    return modified_data, num_replacements


@jit(nopython=True)
def enforce_monotonic_values(data):
    """
    Enforce monotonically decreasing values in the data array.

    Parameters:
    -----------
    data : array-like
        Input data array.

    Returns:
    --------
    tuple: (modified_data, num_adjustments)
        modified_data: Array with enforced monotonically decreasing values.
        num_adjustments: Number of values that were adjusted.
    """
    modified_data = data.copy()
    num_adjustments = 0

    for i in range(1, len(data)):
        if modified_data[i] >= modified_data[i - 1]:
            modified_data[i] = modified_data[i - 1] * 0.98
            num_adjustments += 1

    return modified_data, num_adjustments


@jit(nopython=True)
def z_to_z_interp_wrapper(data_src, z_src, z_target, extrap_low="nan", extrap_high="nan"):
    """
    Interpolate data from one set of z-coordinates to another, with options for extrapolation.

    Parameters:
    -----------
    data_src : array-like
        Source data values corresponding to z_src coordinates.
    z_src : array-like
        Source z-coordinates.
    z_target : array-like
        Target z-coordinates to interpolate onto.
    extrap_low : str or float, optional (default "nan")
        Extrapolation method for values below the lowest z_src:
        - "nan": No extrapolation (returns NaN)
        - "linear": Linear extrapolation
        - "log": Logarithmic extrapolation
        - "persist": Use the lowest valid value
        - "fade": Linear fade to zero
        - float: Use this constant value
    extrap_high : str or float, optional (default "nan")
        Extrapolation method for values above the highest z_src:
        - "nan": No extrapolation (returns NaN)
        - "linear": Linear extrapolation
        - "log": Logarithmic extrapolation
        - "persist": Use the highest valid value
        - "fade": Linear fade to zero
        - float: Use this constant value

    Returns:
    --------
    array-like
        Interpolated data values corresponding to z_target coordinates.

    Notes:
    ------
    - The function uses linear interpolation within the range of z_src.
    - Extrapolation methods are applied outside the range of z_src.
    - The result is returned with the same data type as data_src.
    - This function is optimized with Numba's @jit decorator for performance.

    Warning:
    --------
    Do not use -90, -91, -92, -93, or -94 as extrap_low or extrap_high values, as these are
    reserved for internal use and will not behave as expected.
    """

    # Ensure all inputs are float64
    data_src_64 = data_src.astype(np.float64)
    z_src_64 = z_src.astype(np.float64)
    z_target_64 = z_target.astype(np.float64)

    # Get "integer" option for handling extrapolation
    extrap_low_option = get_extrap_option(extrap_low)
    extrap_high_option = get_extrap_option(extrap_high)

    result = z_to_z_interp(data_src_64, z_src_64, z_target_64, extrap_low_option, extrap_high_option)

    # Convert result back to original dtype of data_src
    return result.astype(data_src.dtype)




# @jit(nopython=True)
# def _hyi2hyo_core_1d(p0, hyai, hybi, psfc, xi, hyao, hybo, intflg):
#     """
#     Core computation for single column or 1D unstructured grid case.
#     """
#     ncol = psfc.shape[0]
#     klevi = len(hyai)
#     klevo = len(hyao)
#     xo = np.full((klevo, ncol), np.nan)
#
#     peps = 1.0e-3
#
#     # Loop over each column
#     for ic in range(ncol):
#         # Calculate input pressure levels
#         pi = hyai * p0 + hybi * psfc[ic]
#         pilow = pi[-1] - peps
#
#         # Calculate output pressure levels
#         po = hyao * p0 + hybo * psfc[ic]
#
#         # Interpolate to each output level
#         for ko in range(klevo):
#             if po[ko] < pi[0] or po[ko] > pilow:
#                 if intflg == 0:
#                     continue  # Leave as NaN
#                 elif intflg == 1:
#                     # Set to nearest input value
#                     if po[ko] < pi[0]:
#                         xo[ko, ic] = xi[0, ic]
#                     else:
#                         xo[ko, ic] = xi[-1, ic]
#             else:
#                 for ki in range(klevi - 1):
#                     if pi[ki] <= po[ko] < pi[ki + 1]:
#                         xo[ko, ic] = xi[ki, ic] + \
#                             (xi[ki + 1, ic] - xi[ki, ic]) * \
#                             (np.log(po[ko]) - np.log(pi[ki])) / \
#                             (np.log(pi[ki + 1]) - np.log(pi[ki]))
#                         break
#
#     return xo
#
# @jit(nopython=True)
# def _hyi2hyo_core_2d(p0, hyai, hybi, psfc, xi, hyao, hybo, intflg):
#     """
#     Core computation for 2D structured grid case (no time dimension).
#     """
#     nlat, nlon = psfc.shape
#     klevi = len(hyai)
#     klevo = len(hyao)
#     xo = np.full((klevo, nlat, nlon), np.nan)
#
#     peps = 1.0e-3
#
#     for j in range(nlat):
#         for i in range(nlon):
#             # Calculate input pressure levels
#             pi = hyai * p0 + hybi * psfc[j, i]
#             pilow = pi[-1] - peps
#
#             # Calculate output pressure levels
#             po = hyao * p0 + hybo * psfc[j, i]
#
#             # Interpolate to each output level
#             for ko in range(klevo):
#                 if po[ko] < pi[0] or po[ko] > pilow:
#                     if intflg == 0:
#                         continue  # Leave as NaN
#                     elif intflg == 1:
#                         # Set to nearest input value
#                         if po[ko] < pi[0]:
#                             xo[ko, j, i] = xi[0, j, i]
#                         else:
#                             xo[ko, j, i] = xi[-1, j, i]
#                 else:
#                     for ki in range(klevi - 1):
#                         if pi[ki] <= po[ko] < pi[ki + 1]:
#                             xo[ko, j, i] = xi[ki, j, i] + \
#                                 (xi[ki + 1, j, i] - xi[ki, j, i]) * \
#                                 (np.log(po[ko]) - np.log(pi[ki])) / \
#                                 (np.log(pi[ki + 1]) - np.log(pi[ki]))
#                             break
#
#     return xo
#
# @jit(nopython=True)
# def _hyi2hyo_core_4d(p0, hyai, hybi, psfc, xi, hyao, hybo, intflg):
#     """Time, lat, lon structured grid case."""
#     ntime, nlat, nlon = psfc.shape
#     klevi = len(hyai)
#     klevo = len(hyao)
#     xo = np.full((ntime, klevo, nlat, nlon), np.nan)
#
#     peps = 1.0e-3
#
#     for t in range(ntime):
#         for j in range(nlat):
#             for i in range(nlon):
#                 pi = hyai * p0 + hybi * psfc[t, j, i]
#                 pilow = pi[-1] - peps
#                 po = hyao * p0 + hybo * psfc[t, j, i]
#
#                 for ko in range(klevo):
#                     if po[ko] < pi[0] or po[ko] > pilow:
#                         if intflg == 0:
#                             continue
#                         elif intflg == 1:
#                             if po[ko] < pi[0]:
#                                 xo[t, ko, j, i] = xi[t, 0, j, i]
#                             else:
#                                 xo[t, ko, j, i] = xi[t, -1, j, i]
#                     else:
#                         for ki in range(klevi - 1):
#                             if pi[ki] <= po[ko] < pi[ki + 1]:
#                                 xo[t, ko, j, i] = xi[t, ki, j, i] + \
#                                     (xi[t, ki + 1, j, i] - xi[t, ki, j, i]) * \
#                                     (np.log(po[ko]) - np.log(pi[ki])) / \
#                                     (np.log(pi[ki + 1]) - np.log(pi[ki]))
#                                 break
#     return xo
#
# @jit(nopython=True)
# def _hyi2hyo_core_time_col(p0, hyai, hybi, psfc, xi, hyao, hybo, intflg):
#     """Time, column unstructured case."""
#     ntime, ncol = psfc.shape
#     klevi = len(hyai)
#     klevo = len(hyao)
#     xo = np.full((ntime, klevo, ncol), np.nan)
#
#     peps = 1.0e-3
#
#     for t in range(ntime):
#         for ic in range(ncol):
#             pi = hyai * p0 + hybi * psfc[t, ic]
#             pilow = pi[-1] - peps
#             po = hyao * p0 + hybo * psfc[t, ic]
#
#             for ko in range(klevo):
#                 if po[ko] < pi[0] or po[ko] > pilow:
#                     if intflg == 0:
#                         continue
#                     elif intflg == 1:
#                         if po[ko] < pi[0]:
#                             xo[t, ko, ic] = xi[t, 0, ic]
#                         else:
#                             xo[t, ko, ic] = xi[t, -1, ic]
#                 else:
#                     for ki in range(klevi - 1):
#                         if pi[ki] <= po[ko] < pi[ki + 1]:
#                             xo[t, ko, ic] = xi[t, ki, ic] + \
#                                 (xi[t, ki + 1, ic] - xi[t, ki, ic]) * \
#                                 (np.log(po[ko]) - np.log(pi[ki])) / \
#                                 (np.log(pi[ki + 1]) - np.log(pi[ki]))
#                             break
#     return xo
#
# def hyi2hyo(p0, hyai, hybi, psfc, xi, hyao, hybo, intflg=0, unstructured=False):
#     """
#     Interpolate data from one set of hybrid levels to another.
#
#     Parameters
#     ----------
#     p0 : float
#         Reference pressure in Pa.
#     hyai : ndarray
#         Input hybrid coefficients A, must be ordered top-to-bottom.
#     hybi : ndarray
#         Input hybrid coefficients B, must be ordered top-to-bottom.
#     psfc : float or ndarray
#         Surface pressure in Pa. Can be:
#         - scalar: single column
#         - 1D array: (ncol,) for unstructured or (nlat,) for structured
#         - 2D array: (time, ncol) for unstructured or (nlat, nlon) for structured
#         - 3D array: (time, nlat, nlon) for structured
#     xi : ndarray
#         Input data to interpolate. Shape must be:
#         - (nlev,) if psfc is scalar
#         - (nlev, ncol) for 1D unstructured or (nlev, nlat) for 1D structured
#         - (time, nlev, ncol) for 2D unstructured or (nlev, nlat, nlon) for 2D structured
#         - (time, nlev, nlat, nlon) for 3D structured
#     hyao : ndarray
#         Output hybrid coefficients A, must be ordered top-to-bottom.
#     hybo : ndarray
#         Output hybrid coefficients B, must be ordered top-to-bottom.
#     intflg : int, optional
#         Interpolation flag:
#         - 0: Set values outside input pressure range to NaN (default)
#         - 1: Set values to nearest input pressure level
#     unstructured : bool, optional
#         If True, treats non-level dimensions as unstructured columns.
#         If False, assumes structured lat/lon grid. Default is False.
#
#     Returns
#     -------
#     ndarray
#         Interpolated data on new hybrid levels.
#         Shape matches input but with lev dimension replaced by new length.
#     """
#     # Basic input validation
#     if not isinstance(p0, (int, float)):
#         raise ValueError("p0 must be a scalar value")
#
#     for arr, name in [(hyai, 'hyai'), (hybi, 'hybi'), (hyao, 'hyao'), (hybo, 'hybo')]:
#         if not isinstance(arr, np.ndarray):
#             raise ValueError(f"{name} must be a numpy array")
#
#     if hyai.shape != hybi.shape or hyao.shape != hybo.shape:
#         raise ValueError("hybrid coefficients must have matching shapes")
#
#     if intflg not in [0, 1]:
#         raise ValueError("intflg must be 0 or 1")
#
#     # Convert scalar psfc to array for consistency
#     if isinstance(psfc, (int, float)):
#         psfc = np.asarray([psfc])
#         xi = np.atleast_2d(xi)
#     else:
#         psfc = np.asarray(psfc)
#
#     # Handle different dimensionality cases
#     if unstructured:
#         if psfc.ndim == 1:  # Single time, columns
#             if xi.ndim != 2 or xi.shape[0] != len(hyai):
#                 raise ValueError("For 1D unstructured psfc, xi must be 2D (nlev, ncol)")
#             return _hyi2hyo_core_1d(p0, hyai, hybi, psfc, xi, hyao, hybo, intflg)
#         elif psfc.ndim == 2:  # Time, columns
#             if xi.ndim != 3 or xi.shape[1] != len(hyai):
#                 raise ValueError("For 2D unstructured psfc, xi must be 3D (time, nlev, ncol)")
#             return _hyi2hyo_core_time_col(p0, hyai, hybi, psfc, xi, hyao, hybo, intflg)
#         else:
#             raise ValueError("For unstructured grids, psfc must be 1D or 2D")
#     else:  # Structured grid
#         if psfc.ndim == 2:  # Single time, lat/lon
#             if xi.ndim != 3 or xi.shape[0] != len(hyai):
#                 raise ValueError("For 2D structured psfc, xi must be 3D (nlev, nlat, nlon)")
#             return _hyi2hyo_core_2d(p0, hyai, hybi, psfc, xi, hyao, hybo, intflg)
#         elif psfc.ndim == 3:  # Time, lat, lon
#             if xi.ndim != 4 or xi.shape[1] != len(hyai):
#                 raise ValueError("For 3D structured psfc, xi must be 4D (time, nlev, nlat, nlon)")
#             return _hyi2hyo_core_4d(p0, hyai, hybi, psfc, xi, hyao, hybo, intflg)
#         else:
#             raise ValueError("For structured grids, psfc must be 2D or 3D")



@jit(nopython=True)
def _hyi2hyo_column(p0, hyai, hybi, ps, xi, hyao, hybo, intflg):
    """
    Core computation for a single vertical column interpolation.
    """
    klevo = len(hyao)
    xo = np.full(klevo, np.nan)

    # Calculate pressures for this column
    pi = hyai * p0 + hybi * ps
    po = hyao * p0 + hybo * ps
    pilow = pi[-1] - 1.0e-3

    # Interpolate each output level
    for ko in range(klevo):
        if po[ko] < pi[0] or po[ko] > pilow:
            if intflg == 0:
                continue
            elif intflg == 1:
                xo[ko] = xi[0] if po[ko] < pi[0] else xi[-1]
        else:
            for ki in range(len(hyai) - 1):
                if pi[ki] <= po[ko] < pi[ki + 1]:
                    xo[ko] = xi[ki] + \
                        (xi[ki + 1] - xi[ki]) * \
                        (np.log(po[ko]) - np.log(pi[ki])) / \
                        (np.log(pi[ki + 1]) - np.log(pi[ki]))
                    break
    return xo

def hyi2hyo(p0, hyai, hybi, psfc, xi, hyao, hybo, intflg=0, unstructured=False):
    """
    Interpolate data from one set of hybrid levels to another.

    Parameters
    ----------
    [same as before]
    """
    # Basic input validation
    if not isinstance(p0, (int, float)):
        raise ValueError("p0 must be a scalar value")

    for arr, name in [(hyai, 'hyai'), (hybi, 'hybi'), (hyao, 'hyao'), (hybo, 'hybo')]:
        if not isinstance(arr, np.ndarray):
            raise ValueError(f"{name} must be a numpy array")

    if hyai.shape != hybi.shape or hyao.shape != hybo.shape:
        raise ValueError("hybrid coefficients must have matching shapes")

    if intflg not in [0, 1]:
        raise ValueError("intflg must be 0 or 1")

    # Handle scalar case
    if isinstance(psfc, (int, float)):
        psfc = np.array([psfc])
        xi = np.atleast_2d(xi)
    else:
        psfc = np.asarray(psfc)

    # Determine if we have a time dimension and level axis
    if unstructured:
        has_time = psfc.ndim == 2  # (time, ncol) vs (ncol,)
        lev_axis = 1 if has_time else 0
    else:
        has_time = psfc.ndim == 3  # (time, nlat, nlon) vs (nlat, nlon)
        lev_axis = 1 if has_time else 0

    print(f"unstructured: {unstructured}")
    print(f"has_time: {has_time}")

    # Verify dimensions match expectations
    if has_time:
        if unstructured and xi.ndim != 3:  # (time, nlev, ncol)
            raise ValueError("For time + unstructured, xi must be 3D (time, nlev, ncol)")
        elif not unstructured and xi.ndim != 4:  # (time, nlev, nlat, nlon)
            raise ValueError("For time + structured, xi must be 4D (time, nlev, nlat, nlon)")
    else:
        if unstructured and xi.ndim != 2:  # (nlev, ncol)
            raise ValueError("For unstructured, xi must be 2D (nlev, ncol)")
        elif not unstructured and xi.ndim != 3:  # (nlev, nlat, nlon)
            raise ValueError("For structured, xi must be 3D (nlev, nlat, nlon)")

    if xi.shape[lev_axis] != len(hyai):
        raise ValueError("Level dimension of xi must match length of hybrid coefficients")

    # Prepare output array
    out_shape = list(xi.shape)
    out_shape[lev_axis] = len(hyao)
    xo = np.full(out_shape, np.nan)

    # Handle different dimensionality cases through iteration
    if has_time:
        for t in range(xi.shape[0]):
            if unstructured:
                # Time, level, column case
                for c in range(psfc.shape[1]):
                    xo[t, :, c] = _hyi2hyo_column(p0, hyai, hybi, psfc[t, c],
                                               xi[t, :, c], hyao, hybo, intflg)
            else:
                # Time, level, lat, lon case
                for j in range(psfc.shape[1]):
                    for i in range(psfc.shape[2]):
                        xo[t, :, j, i] = _hyi2hyo_column(p0, hyai, hybi, psfc[t, j, i],
                                                      xi[t, :, j, i], hyao, hybo, intflg)
    else:
        if unstructured:
            # Level, column case
            for c in range(psfc.shape[0]):
                xo[:, c] = _hyi2hyo_column(p0, hyai, hybi, psfc[c],
                                        xi[:, c], hyao, hybo, intflg)
        else:
            # Level, lat, lon case
            for j in range(psfc.shape[0]):
                for i in range(psfc.shape[1]):
                    xo[:, j, i] = _hyi2hyo_column(p0, hyai, hybi, psfc[j, i],
                                               xi[:, j, i], hyao, hybo, intflg)

    return xo


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
        elif sigo[k] >= sigi[nli - 2]:
            kp = nli - 2
        else:
            kp = 0
            while sigo[k] > sigi[kp + 1]:
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
    logging.info(f"Sigma to hybrid interpolation completed in {elapsed_time:.2f} seconds.")

    return dato