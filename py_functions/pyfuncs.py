import numpy as np
import xarray as xr
import argparse
import os
import sys
import cftime
from scipy.ndimage import gaussian_filter
from numba import jit
import logging
logger = logging.getLogger(__name__)


def configure_logging(verbose=False):
    logging_level = logging.DEBUG if verbose else logging.INFO
    logging.basicConfig(
        level=logging_level,
        datefmt='%Y-%m-%d %H:%M:%S',
        format='%(asctime)s - %(levelname)s - %(message)s',
        handlers=[
            logging.FileHandler("atm_to_cam.log"),
            logging.StreamHandler(sys.stdout)
        ]
    )
    numba_logger = logging.getLogger('numba')
    numba_logger.setLevel(logging.WARNING)
    mpl_logger = logging.getLogger('matplotlib')
    mpl_logger.setLevel(logging.WARNING)
    if verbose:
        logging.info("Verbose mode is on. More detailed logging information will be shown.")


def get_betacast_path():
    BETACAST = os.getenv("BETACAST")
    if BETACAST is None:
        logging.info("We are running local-only atm_to_cam. Set export BETACAST env to run elsewhere")
        local_only = True
        PATHTOHERE = "./"
    else:
        logging.info("Not local only!")
        local_only = False
        PATHTOHERE = os.path.join(BETACAST, "atm_to_cam")

    return local_only, PATHTOHERE


def parse_args():
    parser = argparse.ArgumentParser(description="Process command-line arguments for climate data processing.")

    # Required arguments
    parser.add_argument('--datasource', type=str, required=True,
                        help='Data source: CFSR, ERA5, ERAI, GFS, CAM, RAP')
    parser.add_argument('--numlevels', type=int, required=True,
                        help='Number of vertical levels')
    parser.add_argument('--YYYYMMDDHH', type=int, required=True,
                        help='Initialization date in YYYYMMDDHH format')
    parser.add_argument('--data_filename', type=str, required=True,
                        help='Full path to file containing initial information')
    parser.add_argument('--wgt_filename', type=str, required=True,
                        help='Full path to ESMF weight file from ANL -> MOD')
    parser.add_argument('--se_inic', type=str, required=True,
                        help='Full path of file to write as output')

    # Optional arguments with defaults
    parser.add_argument('--dycore', type=str, default='se',
                        help='Dycore option: fv, se, mpas (default: se)')
    parser.add_argument('--RDADIR', type=str, default='',
                        help='Base path to RDA folder (default: "")')
    parser.add_argument('--mpas_as_cam', action='store_true',
                        help='If set, write MPAS in CAM physics output (default: False)')
    parser.add_argument('--compress_file', action='store_true',
                        help='If set, will attempt NetCDF "chunking" compression (default: False)')
    parser.add_argument('--write_floats', action='store_true',
                        help='If set, write outputs as single precision (default: False)')
    parser.add_argument('--add_cloud_vars', action='store_true', default=True,
                        help='If set, add CLDICE and CLDLIQ to output file (default: True)')
    parser.add_argument('--adjust_config', type=str, default=' ',
                        help='String defining how to perform hydro adjustment (default: "")')
    parser.add_argument('--model_topo_file', type=str, default='',
                        help='File containing PHIS for FV or SE, or an MPAS inic file (default: "")')
    parser.add_argument('--mod_in_topo', type=str, default='',
                        help='Full path to PHIS field from downscaling MOD for pressure surface calculation (default: "")')
    parser.add_argument('--mod_remap_file', type=str, default='',
                        help='Full path to ESMF weight file that goes downscaling MOD -> ANL (default: "")')
    parser.add_argument('--mpasfile', type=str, default='',
                        help='XXXXXXX')
    parser.add_argument('-v', '--verbose', action='store_true', help='Increase output verbosity for debugging')

    return parser.parse_args()


def split_by_lengths(s, lengths):
    """
    Splits a string `s` into parts defined by `lengths`.
    """
    parts = []
    index = 0
    for length in lengths:
        parts.append(s[index:index + length])
        index += length
    return parts


def find_closest_time(times, yearstr, monthstr, daystr, cyclestr, return_isel=False):
    if isinstance(times, np.ndarray):
        # Convert the strings to an ISO 8601 datetime format
        datetime_str = f"{yearstr}-{monthstr.zfill(2)}-{daystr.zfill(2)}T00:00:00"
        datetime_base = np.datetime64(datetime_str)

        # Add the number of seconds from cyclestr to get the exact datetime
        datetime_target = datetime_base + np.timedelta64(int(cyclestr), 's')
        logging.info(f"Target datetime: {datetime_target}")

        closest_index = np.argmin(np.abs(times - datetime_target))
        if return_isel:
            closest_time = closest_index
        else:
            closest_time = times[closest_index]
    else:
        target_time = np.datetime64(f"{yearstr}-{monthstr}-{daystr}T{cyclestr[:2]}:00")
        closest_time = min(times.values, key=lambda x: abs(x - target_time))

    return closest_time


def numpy_to_dataarray(numpy_array, dims, coords=None, attrs=None, name=None):
    """
    Converts a numpy array to an xarray DataArray.

    Parameters:
    - numpy_array: The numpy array to convert.
    - dims: A list of dimension names for the DataArray.
    - coords: Optional dictionary of coordinate arrays.
    - attrs: Optional dictionary of attributes to add to the DataArray.
    - name: Optional name for the DataArray.

    Returns:
    - data_array: The resulting xarray DataArray.
    """
    data_array = xr.DataArray(data=numpy_array, dims=dims, coords=coords, attrs=attrs, name=name)
    return data_array


def add_time_define_precision(var_in, precision, isncol, lat_dim="lat", lon_dim="lon"):
    """
    Adds a time dimension and converts the precision of the input variable.

    Parameters:
    var_in (xarray.DataArray or numpy.ndarray): Input variable.
    precision (str): Desired precision ("float", "single", "double").
    isncol (bool): True if the variable is ncol; False otherwise.
    lat_dim (str): Name of the latitude dimension (default is "lat").
    lon_dim (str): Name of the longitude dimension (default is "lon").

    Returns:
    xarray.DataArray: Output variable with time dimension and specified precision.
    """
    var_dims = var_in.shape

    if isncol:
        if len(var_dims) == 1:  # ncol -> time, ncol
            logging.info(f"Converting a ncol -> ncol {precision} with time attached!")
            ncol = var_dims[0]
            if precision in ["float", "single"]:
                var_out = xr.DataArray(np.zeros((1, ncol), dtype=np.float32), dims=["time", "ncol"])
            elif precision == "double":
                var_out = xr.DataArray(np.zeros((1, ncol), dtype=np.float64), dims=["time", "ncol"])
            else:
                raise ValueError("Invalid precision specified")
            var_out[0, :] = var_in.astype(var_out.dtype)

        elif len(var_dims) == 2:  # lev, ncol -> time, lev, ncol
            logging.info(f"Converting a lev, ncol -> lev, ncol {precision} with time attached!")
            nlev, ncol = var_dims
            if precision in ["float", "single"]:
                var_out = xr.DataArray(np.zeros((1, nlev, ncol), dtype=np.float32), dims=["time", "lev", "ncol"])
            elif precision == "double":
                var_out = xr.DataArray(np.zeros((1, nlev, ncol), dtype=np.float64), dims=["time", "lev", "ncol"])
            else:
                raise ValueError("Invalid precision specified")
            var_out[0, :, :] = var_in.astype(var_out.dtype)

    else:  # Not ncol
        if len(var_dims) == 2:  # lat_dim, lon_dim -> time, lat_dim, lon_dim
            logging.info(f"Converting a {lat_dim}, {lon_dim} -> {lat_dim}, {lon_dim} {precision} with time attached!")
            nlat, nlon = var_dims
            if precision in ["float", "single"]:
                var_out = xr.DataArray(np.zeros((1, nlat, nlon), dtype=np.float32), dims=["time", lat_dim, lon_dim])
            elif precision == "double":
                var_out = xr.DataArray(np.zeros((1, nlat, nlon), dtype=np.float64), dims=["time", lat_dim, lon_dim])
            else:
                raise ValueError("Invalid precision specified")
            var_out[0, :, :] = var_in.astype(var_out.dtype)

        elif len(var_dims) == 3:  # lev, lat_dim, lon_dim -> time, lev, lat_dim, lon_dim
            logging.info(f"Converting a lev, {lat_dim}, {lon_dim} -> lev, {lat_dim}, {lon_dim} {precision} with time attached!")
            nlev, nlat, nlon = var_dims
            if precision in ["float", "single"]:
                var_out = xr.DataArray(np.zeros((1, nlev, nlat, nlon), dtype=np.float32), dims=["time", "lev", lat_dim, lon_dim])
            elif precision == "double":
                var_out = xr.DataArray(np.zeros((1, nlev, nlat, nlon), dtype=np.float64), dims=["time", "lev", lat_dim, lon_dim])
            else:
                raise ValueError("Invalid precision specified")
            var_out[0, :, :, :] = var_in.astype(var_out.dtype)

    # Copy any metadata over to the output variable
    var_out.attrs = var_in.attrs

    return var_out

def replace_nans_with_fill_value(data_horiz, variables, NC_FLOAT_FILL):
    """
    Replace NaN values with _FillValue in specified numpy arrays within a dictionary.

    Parameters:
    -----------
    NC_FLOAT_FILL : float
        The fill value to replace NaN values with.

    data_horiz : dict
        A dictionary containing numpy arrays where NaNs should be replaced.

    variables : list of str
        A list of keys in data_horiz specifying which arrays to modify.

    Returns:
    --------
    dict
        The updated data_horiz dictionary with NaNs replaced in the specified variables.
    """
    for var in variables:
        if var in data_horiz:
            data_horiz[var] = np.where(np.isnan(data_horiz[var]), NC_FLOAT_FILL, data_horiz[var])
        else:
            raise KeyError(f"Variable {var} not found in data_horiz.")
    return data_horiz


def clip_and_count(arr, min_thresh=None, max_thresh=None, var_name="Variable"):
    """
    Clips the values in the array to the specified minimum and maximum thresholds,
    and returns the clipped array along with the number of adjustments made.

    Parameters:
    - arr: numpy.ndarray
        The array to be clipped.
    - min_thresh: float or None
        The minimum threshold value. Values below this will be set to min_thresh.
    - max_thresh: float or None
        The maximum threshold value. Values above this will be set to max_thresh.
    - var_name: str
        The name of the variable for diagnostic output.

    Returns:
    - clipped_arr: numpy.ndarray
        The clipped array.
    """

    # Apply the maximum threshold if specified
    if max_thresh is not None:
        max_adjustments = np.sum(arr > max_thresh)
        arr = np.where(arr > max_thresh, max_thresh, arr)
        logging.info(f"Number of {var_name} adjustments for values above {max_thresh}: {max_adjustments}")

    # Apply the minimum threshold if specified
    if min_thresh is not None:
        min_adjustments = np.sum(arr < min_thresh)
        arr = np.where(arr < min_thresh, min_thresh, arr)
        logging.info(f"Number of {var_name} adjustments for values below {min_thresh}: {min_adjustments}")

    if np.any(np.isnan(arr)):
        logging.info(f"WARNING: {var_name} is missing data...")
        logging.info("This is expected for a regional dataset, but not desirable for a global one")

    return arr


def create_cf_time(year, month, day, hour, base_time="1850-01-01 00:00:00", time_units="days since", calendar="standard"):
    """
    Create a CF-compliant time array from year, month, day, and hour.

    Parameters:
    - year, month, day, hour: Arrays or scalars of year, month, day, and hour
    - base_time: The reference time in the format "YYYY-MM-DD HH:MM:SS"
    - time_units: The time units for the output (e.g., "days since" or "hours since")

    Returns:
    - time: A numpy array of CF-compliant time values
    - time_attrs: A dictionary with attributes to assign to the time variable in xarray
    """

    base_year, base_month, base_day, base_hour = [int(val) for val in base_time.split()[0].split('-')] + [int(base_time.split()[1].split(':')[0])]

    nptime = cftime.date2num(cftime.DatetimeNoLeap(year, month, day, hour), units=f"{time_units} {base_time}")

    nptime = np.array([nptime], dtype=np.float64)

    time_attrs = {
        "units": f"{time_units} {base_time}",
        "calendar": f"{calendar}"
    }

    return nptime, time_attrs


def print_all_variables_info(name_width=22, type_width=40, shape_width=18, value_width=50, size_width=15):
    # Get the current frame
    frame = sys._getframe(1)

    # Get all local variables in the current frame
    local_vars = frame.f_locals

    logging.info(f"{'Variable Name':<{name_width}} {'Type':<{type_width}} {'Shape':<{shape_width}} {'Value':<{value_width}} {'Size (bytes)':<{size_width}}")
    logging.info("=" * (name_width + type_width + shape_width + value_width + size_width))

    for var_name, var_value in local_vars.items():
        var_type = str(type(var_value))[:type_width]
        var_shape = ''
        var_size = sys.getsizeof(var_value)

        # Determine the shape and size for numpy arrays and xarray DataArrays
        if isinstance(var_value, (np.ndarray, xr.DataArray)):
            var_shape = str(var_value.shape)[:shape_width]
            var_size = var_value.nbytes if isinstance(var_value, np.ndarray) else var_value.nbytes
        elif isinstance(var_value, xr.Dataset):
            var_shape = str({dim: size for dim, size in var_value.sizes.items()})[:shape_width]
            var_size = sum([var.nbytes for var in var_value.data_vars.values()])

        # Truncate the value for xarray.Dataset to avoid multi-line printing
        if isinstance(var_value, xr.Dataset):
            display_value = "<xarray.Dataset>"[:value_width]
        elif isinstance(var_value, xr.DataArray):
            flattened_value = var_value.values.flatten()
            display_value = f"{flattened_value[0]}..." if flattened_value.size > 0 else "Empty"
        elif isinstance(var_value, np.ndarray):
            display_value = f"{var_value.flatten()[0]}..." if var_value.size > 0 else "Empty"
        else:
            display_value = str(var_value)[:value_width]

        logging.info(f"{var_name:<{name_width}} {var_type:<{type_width}} {var_shape:<{shape_width}} {display_value:<{value_width}} {var_size:<{size_width}}")


def print_debug_file(output_filename, **kwargs):
    """
    Create an xarray Dataset from the given variables and save it to a NetCDF file.

    Parameters:
    - output_filename: Name of the output NetCDF file
    - **kwargs: Arbitrary number of variables to include in the Dataset.
                Each key should be the variable name.
                Each value should be a tuple where the first element is a list of dimension names,
                and the second element is the variable data (numpy array).
    """
    # Prepare the data dictionary for the Dataset
    data_dict = {}
    for var_name, (dims, var_data) in kwargs.items():
        if isinstance(var_data, np.ndarray):
            # Add the variable to the Dataset dictionary
            data_dict[var_name] = (dims, var_data.astype(np.float32))
        else:
            raise ValueError(f"Unsupported data type for variable '{var_name}': {type(var_data)}")

    # Create an xarray.Dataset
    ds = xr.Dataset(data_dict)

    # Save the dataset to a NetCDF file
    ds.to_netcdf(output_filename)


@jit(nopython=True)
def dsmth9(X, P, Q, LWRAP, XMSG=np.nan):
    """
    Perform 9-point smoothing on a 2D array X.

    Parameters:
    - X: 2D input/output array (numpy array)
    - P: First weight (suggested value 0.50)
    - Q: Second weight (suggested value 0.25)
    - XMSG: Value of missing points
    - LWRAP: Boolean flag to include wraparound points in smoothing

    Returns:
    - X: Smoothed 2D array
    - IER: Error code (0 for success, 1 if NI or NJ is less than 3)
    """
    NI, NJ = X.shape

    if NI < 3 or NJ < 3:
        print(f"Too few points in smth9: error {NI} {NJ}")
        return X

    # Precompute constants
    PO4 = P / 4.0
    QO4 = Q / 4.0

    # Initialize the work array WRK
    WRK = np.full_like(X, XMSG)

    # Set the loop bounds depending on the wraparound flag
    if LWRAP:
        NIB = 0
        NIE = NI
    else:
        NIB = 1
        NIE = NI - 1
    NJB = 1
    NJE = NJ - 1

    # Perform smoothing
    for J in range(NJB, NJE):
        for I in range(NIB, NIE):
            JM1 = J - 1
            JP1 = J + 1
            IM1 = I - 1 if I > 0 else NI - 1
            IP1 = I + 1 if I < NI - 1 else 0

            if (
                X[I, J] == XMSG or X[IM1, JP1] == XMSG or X[IM1, J] == XMSG or
                X[IM1, JM1] == XMSG or X[I, JM1] == XMSG or X[IP1, JM1] == XMSG or
                X[IP1, J] == XMSG or X[IP1, JP1] == XMSG or X[I, JP1] == XMSG
            ):
                WRK[I, J] = X[I, J]
            else:
                TERM1 = PO4 * (X[IM1, J] + X[I, JM1] + X[IP1, J] + X[I, JP1] - 4.0 * X[I, J])
                TERM2 = QO4 * (X[IM1, JP1] + X[IM1, JM1] + X[IP1, JM1] + X[IP1, JP1] - 4.0 * X[I, J])
                WRK[I, J] = X[I, J] + TERM1 + TERM2

    # Transfer back to original array
    X[NIB:NIE, NJB:NJE] = WRK[NIB:NIE, NJB:NJE]

    return X


def latRegWgt(lat, nType="float", opt=0):
    """
    Generates [sin(lat+dlat/2) - sin(lat-dlat/2)] weights for equally spaced (regular) global grids that will sum to 2.0.

    Parameters:
    - lat: 1D array-like of latitudes of the global grid (in degrees).
    - nType: The type of variable to be returned ("float" or "double").
    - opt: Not used. Set to 0.

    Returns:
    - weights: A 1D array of weights of the same size as lat.
    """

    # Convert lat to radians
    lat_rad = np.radians(lat)
    dlat = np.abs(np.diff(lat_rad).mean())  # Assume uniform spacing and get the mean difference

    # Calculate weights
    weights = np.sin(lat_rad + dlat / 2) - np.sin(lat_rad - dlat / 2)

    # Ensure the weights sum to 2.0
    weights *= 2.0 / np.sum(weights)

    # Convert to the specified type
    if nType == "float":
        weights = weights.astype(np.float32)
    elif nType == "double":
        weights = weights.astype(np.float64)
    else:
        raise ValueError("nType must be 'float' or 'double'.")

    return weights


def smooth_with_smth9(var, numiter, p=0.5, q=0.25):
    if var.ndim == 3:
        for ii in range(numiter):
            logging.info(f"SMOOTH (SMTH9) ITER: {ii+1}")
            for level in range(var.shape[0]):
                smoothed = dsmth9(var[level, :, :], p, q, True)
                if smoothed.shape != var[level, :, :].shape:
                    raise ValueError(f"Shape mismatch after smoothing. Expected {var[level, :, :].shape}, but got {smoothed.shape}.")
                var[level, :, :] = smoothed
    elif var.ndim == 2:
        for ii in range(numiter):
            logging.info(f"SMOOTH (SMTH9) ITER: {ii+1}")
            smoothed = dsmth9(var, p, q, True)
            if smoothed.shape != var.shape:
                raise ValueError(f"Shape mismatch after smoothing. Expected {var.shape}, but got {smoothed.shape}.")
            var = smoothed
    else:
        raise ValueError("Input array must be either 2D or 3D.")

    return var


def smooth_with_gaussian(var, numiter, sigma=1, truncate=4.0):

    for ii in range(numiter):
        logging.info(f"SMOOTH (SCIPY GAUS) ITER: {ii+1}")
        var = gaussian_filter(var, sigma=sigma, truncate=truncate)

    return var


def print_and_return_varsize(ps_fv, u_fv, v_fv, t_fv, q_fv):
    """
    Function to print and return the maximum size of xarray DataArray variables in bytes.

    Parameters:
    -----------
    ps_fv, u_fv, v_fv, t_fv, q_fv : xarray.DataArray
        Variables whose sizes are to be determined.

    Returns:
    --------
    max_size : int
        The maximum size among the provided variables in bytes.
    """

    ps_sz = ps_fv.nbytes
    u_sz = u_fv.nbytes
    v_sz = v_fv.nbytes
    t_sz = t_fv.nbytes
    q_sz = q_fv.nbytes

    logging.info("~~~ VAR SIZES (bytes)")
    logging.info(f"ps: {ps_sz}")
    logging.info(f"u: {u_sz}")
    logging.info(f"v: {v_sz}")
    logging.info(f"t: {t_sz}")
    logging.info(f"q: {q_sz}")

    max_size = max(ps_sz, u_sz, v_sz, t_sz, q_sz)
    logging.info(f"Var max size: {max_size}")

    return max_size


def print_min_max_dict(data_vars):
    """
    Prints the minimum and maximum values for each key in the data_vars dictionary.

    Parameters:
    -----------
    data_vars : dict
        A dictionary where each key corresponds to a variable name and its value is a NumPy array or xarray DataArray.
    """
    logging.info("=" * 65)
    for key, data in data_vars.items():
        if isinstance(data, (np.ndarray, xr.DataArray)):
            logging.info(f"{key.upper()} shape: {data.shape} | Max: {data.max()} | Min: {data.min()}")
    logging.info("=" * 65)