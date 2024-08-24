import typing
import warnings
import metpy.interpolate
import numpy as np
from numba import jit

supported_types = typing.Union[np.ndarray]

__pres_lev_mandatory__ = np.array([
    1000, 925, 850, 700, 500, 400, 300, 250, 200, 150, 100, 70, 50, 30, 20, 10,
    7, 5, 3, 2, 1
]).astype(np.float64) * 100.0  # Convert mb to Pa

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

        #print("------------")
        #print(np.exp(log_new_levels))
        #print(np.exp(lp_slice))
        #print(dp_slice)
        # Interpolate over the 1D slice
        interp_slice = np.interp(log_new_levels, lp_slice, dp_slice, left=-999.9, right=fill_value)
        #print(interp_slice)
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
        #print(interp_slice)
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
    output = np.where(new_levels_broadcasted > bottom_level-0.1, extrapolated, output)

    # Debugging: Count and print number of NaNs after downward extrapolation
    nan_count_downward = np.isnan(output).sum()
    if nan_count_downward > 0:
        print(f"NaNs detected after downward extrapolation: {nan_count_downward}")
        nan_indices_downward = np.argwhere(np.isnan(output))
        print("NaNs found at indices (downward extrapolation):")
        print(nan_indices_downward)
        print("Inspecting inputs at NaN indices for downward extrapolation:")
        for idx in nan_indices_downward:
            print(f"Index {idx}:")

            # Ensure the index is within bounds for the data array
            if idx[0] < data.shape[0] and all(i < s for i, s in zip(idx[1:], data.shape[1:])):
                print(f"data at index: {data[tuple(idx)]}")
            else:
                print(f"data at index: Index out of bounds for data array.")

            # Ensure the index is within bounds for the bottom_level, ps, and phi_sfc arrays
            if all(i < s for i, s in zip(idx[1:], bottom_level.shape)):
                print(f"bottom_level at index: {bottom_level[tuple(idx[1:])]}")
                print(f"ps at index: {ps[tuple(idx[1:])]}")
                print(f"phi_sfc at index: {phi_sfc[tuple(idx[1:])]}")
            else:
                print(f"bottom_level, ps, or phi_sfc index out of bounds.")

            # Print the corresponding new_levels value
            if idx[0] < len(new_levels):
                print(f"new_levels at level index: {new_levels[idx[0]]}")
            else:
                print(f"new_levels index out of bounds.")

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