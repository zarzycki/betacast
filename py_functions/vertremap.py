import time
import logging
import warnings

import numpy as np
from numba import jit

logger = logging.getLogger(__name__)

__pres_lev_mandatory__ = np.array([
    1000, 925, 850, 700, 500, 400, 300, 250, 200, 150, 100, 70, 50, 30, 20, 10,
    7, 5, 3, 2, 1
]).astype(np.float64) * 100.0  # Convert mb to Pa


def interpolate_to_levels(ppin, xxin, ppout, linlog=1, xmsg=np.nan):
    """
    Interpolate or extrapolate data from one set of pressure levels to another.

    Parameters:
    -----------
    ppin : numpy.ndarray
        Input pressure levels (1D array).
    xxin : numpy.ndarray
        Input variable values at ppin levels (1D array).
    ppout : numpy.ndarray
        Output pressure levels (1D array).
    linlog : int, optional
        Flag indicating interpolation method:
        - 1: linear interpolation (default)
        - 0: logarithmic interpolation
        - Negative: extrapolation allowed using the nearest valid slope.
    xmsg : float, optional
        Missing data marker. Default is np.nan.

    Returns:
    --------
    xxout : numpy.ndarray
        Interpolated variable values at ppout levels.
    """

    # Initialize output array with missing values
    xxout = np.full_like(ppout, xmsg)

    # Error check: make sure we have enough points
    if len(ppin) < 2 or len(ppout) < 1:
        return xxout  # Return all missing values

    # Determine if input and output pressures need to be reordered
    ppin_reordered = ppin[::-1] if ppin[0] < ppin[1] else ppin
    xxin_reordered = xxin[::-1] if ppin[0] < ppin[1] else xxin
    ppout_reordered = ppout[::-1] if ppout[0] < ppout[1] else ppout

    # Filter out missing data in xxin
    valid_mask = (xxin_reordered != xmsg) & (ppin_reordered != xmsg)
    pin_valid = ppin_reordered[valid_mask]
    xin_valid = xxin_reordered[valid_mask]

    if len(pin_valid) < 2:
        return xxout  # Not enough valid data points

    # Perform interpolation
    if linlog == 1:
        # Linear interpolation
        xxout = np.interp(ppout_reordered, pin_valid, xin_valid)
    else:
        # Logarithmic interpolation
        xxout = np.interp(np.log(ppout_reordered), np.log(pin_valid), xin_valid)

    # Extrapolation if linlog < 0
    if linlog < 0:
        if ppout_reordered[0] > pin_valid[0]:  # Extrapolate above range
            slope = (xin_valid[1] - xin_valid[0]) / (pin_valid[1] - pin_valid[0])
            xxout[ppout_reordered > pin_valid[0]] = xin_valid[0] + slope * (ppout_reordered[ppout_reordered > pin_valid[0]] - pin_valid[0])
        if ppout_reordered[-1] < pin_valid[-1]:  # Extrapolate below range
            slope = (xin_valid[-1] - xin_valid[-2]) / (pin_valid[-1] - pin_valid[-2])
            xxout[ppout_reordered < pin_valid[-1]] = xin_valid[-1] + slope * (ppout_reordered[ppout_reordered < pin_valid[-1]] - pin_valid[-1])

    # Reorder output back to original if necessary
    if ppout[0] < ppout[1]:
        xxout = xxout[::-1]

    return xxout


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
    logging.info(f"Interpolation completed in {elapsed_time:.2f} seconds.")

    return data_on_hybrid_levels


def pres2hyb_all(data_vars, ps, hya, hyb):
    """
    Interpolate data on constant pressure levels to hybrid sigma levels.

    This function is a wrapper that prepares the input for the JIT-compiled version.
    """

    data_vint = {}

    allowable_interp_vars = ['t', 'u', 'v', 'q', 'cldice', 'cldliq', 'z', 'theta', 'rho', 'w']

    # Create non-interpolated vars
    data_vint['hya'] = hya
    data_vint['hyb'] = hyb

    # Loop ovr the keys in data_in and interpolate if the key is in allowable_interp_vars
    for key in data_vars:
        if key in allowable_interp_vars:
            data_vint[key] = pressure_to_hybrid(data_vars['lev'], data_vars[key], ps, data_vint['hya'], data_vint['hyb'])
        else:
            data_vint[key] = data_vars[key]

    return data_vint


def interp_hybrid_to_pressure_wrapper(data_vars, ps, hyam, hybm, new_levels, lev_dim=0, method='log', extrapolate=True):

    allowable_interp_vars = ['t', 'u', 'v', 'q', 'cldice', 'cldliq', 'z', 'theta', 'rho', 'w']

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

    t_wrf_col = z_to_z_interp(t_fv, z_fv, zmid, extrapLow=True, extrapHigh=True)

    theta_wrf_col = z_to_z_interp(theta_fv, z_fv, zmid, extrapLow=True, extrapHigh=True)
    rho_wrf_col = z_to_z_interp(rho_fv, z_fv, zmid, extrapLow=True, extrapHigh=True)
    w_wrf_col = z_to_z_interp(w_fv, z_fv, zint, extrapLow=True, extrapHigh=True)
    q_wrf_col = z_to_z_interp(q_fv, z_fv, zmid, extrapLow=True, extrapHigh=True)

    # u and v we don't extrapolate aloft to prevent wind speeds from getting too high
    u_wrf_col = z_to_z_interp(u_fv, z_fv, zmid, extrapLow=True, extrapHigh=False)
    v_wrf_col = z_to_z_interp(v_fv, z_fv, zmid, extrapLow=True, extrapHigh=False)

#     if ix == 30:
#       logging.info(u_wrf_col)
#       logging.info(u_fv[::-1])
#       logging.info(zmid)
#       logging.info(z_fv[::-1])

    return t_wrf_col, theta_wrf_col, rho_wrf_col, w_wrf_col, q_wrf_col, u_wrf_col, v_wrf_col


@jit(nopython=True)
def z_to_z_interp(theta_fv, z_fv, thisCol, extrapLow=False, extrapHigh=False):
    nlev = len(thisCol)
    t_wrf = np.zeros(nlev)

    for i in range(nlev):
        # logger.debug(f"thisCol[i]: {thisCol[i]}, z_fv[0] {z_fv[0]}, z_fv[-1], {z_fv[-1]}")
        if thisCol[i] <= z_fv[0] or thisCol[i] >= z_fv[-1]:
            t_wrf[i] = np.nan
        else:
            for j in range(len(z_fv) - 1):
                if z_fv[j] <= thisCol[i] <= z_fv[j + 1]:
                    t_wrf[i] = theta_fv[j] + (theta_fv[j + 1] - theta_fv[j]) * (thisCol[i] - z_fv[j]) / (z_fv[j + 1] - z_fv[j])
                    break

#     if np.all(np.isnan(t_wrf)):
#         logging.error("Interpolation failed: t_wrf array is entirely NaN. Exiting program.")
#         logging.error("This *may* be because the arrays you have sent in are flipped...")
#         sys.exit(1)

    if np.any(np.isnan(t_wrf)):
        ixvalid = np.where(~np.isnan(t_wrf))[0]

        if np.isnan(t_wrf[0]):
            lowest_valid = np.min(ixvalid)
            if extrapLow:
                deriv = (t_wrf[lowest_valid] - t_wrf[lowest_valid + 1]) / (thisCol[lowest_valid] - thisCol[lowest_valid + 1])
            else:
                deriv = 0.0
            for ff in range(0, lowest_valid):
                t_wrf[ff] = t_wrf[lowest_valid] + deriv * (thisCol[ff] - thisCol[lowest_valid])

        if np.isnan(t_wrf[-1]):
            highest_valid = np.max(ixvalid)
            if extrapHigh:
                deriv = (t_wrf[highest_valid - 1] - t_wrf[highest_valid]) / (thisCol[highest_valid - 1] - thisCol[highest_valid])
            else:
                deriv = 0.0
            for ff in range(highest_valid + 1, nlev):
                t_wrf[ff] = t_wrf[highest_valid] + deriv * (thisCol[ff] - thisCol[highest_valid])

    return t_wrf


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