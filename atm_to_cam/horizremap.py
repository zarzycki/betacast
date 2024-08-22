import numpy as np
import xarray as xr
#from numba import jit
import time
import scipy.sparse as sp

def remap_with_weights(src_data, sparse_map, dst_grid_dims, src_grid_type, dst_grid_type, **kwargs):
    """
    Regrid data using ESMF weights.

    Parameters:
    - src_data: 2D numpy array of the source data (e.g., t_cam).
    - wgt_filename: Path to the ESMF weight file.
    - datafile: Path to the netCDF file containing the source data.
    - var_name: Variable name in the source datafile to apply the map onto.
    - kwargs: Additional keyword arguments.

    Returns:
    - data_out: Regridded data as a 2D numpy array.
    """

    if src_grid_type == "structured":
        FLATTENTYPE="C"
    elif src_grid_type == "unstructured":
        FLATTENTYPE="F"
    else:
        raise ValueError(f"Unknown grid type: {src_grid_type}")

    # Flatten the source data to create a vector of length n_a
    # It appears ESMF uses "C" ordering, although we may need to use "F" if src is unstructured?
    src_vector = src_data.flatten(order=FLATTENTYPE)

    # Apply map onto the src column vector length n_a to compute vector length n_b
    field_target = sparse_map @ src_vector

    # Reshape 1-D vector returned to dst_grid_dims

    data_out = np.reshape(field_target, dst_grid_dims, order="F")

    return data_out

def remap_with_weights_wrapper(src_data, wgt_filename, **kwargs):
    """
    Wrapper to handle regridding of multi-dimensional data using ESMF weights.

    Parameters:
    - src_data: Multi-dimensional numpy array of the source data.
    - wgt_filename: Path to the ESMF weight file.

    Returns:
    - regridded_data: Multi-dimensional numpy array after regridding.
    """

    start_time = time.time()

    # Open the weight file
    xwgt = xr.open_dataset(wgt_filename)
    srclat = xwgt['yc_a']
    srclon = xwgt['xc_a']
    dstlat = xwgt['yc_b']
    dstlon = xwgt['xc_b']
    src_grid_dims = xwgt['src_grid_dims'].values
    dst_grid_dims = xwgt['dst_grid_dims'].values

    print(f"Src grid dims: {src_grid_dims}, dst grid dims: {dst_grid_dims}")
    src_grid_type = "structured" if src_grid_dims.size == 2 else "unstructured"
    dst_grid_type = "structured" if dst_grid_dims.size == 2 else "unstructured"
    print(f"Src grid type: {src_grid_type}, Dst grid type: {dst_grid_type}")

    n_a = xwgt['n_a'].size  # col dimension
    n_b = xwgt['n_b'].size  # row dimension
    n_s = xwgt['n_s'].size  # nnz dimension

    print("Map contains {0} rows, {1} cols and {2} nnz values".format(n_b, n_a, n_s))

    rows = xwgt['row'][:] - 1  # row indices (1-based to 0-based)
    cols = xwgt['col'][:] - 1  # col indices (1-based to 0-based)
    nnzvals = xwgt['S'][:]  # nnz map values

    # Create sparse matrix map
    sparse_map = sp.coo_matrix((nnzvals, (rows, cols)), shape=(n_b, n_a))

    src_dims = src_data.shape
    if src_grid_type == "structured":
        # Structured grid: last two dimensions are nlat, nlon
        nlat, nlon = src_dims[-2], src_dims[-1]
        extra_dims = src_dims[:-2]
    elif src_grid_type == "unstructured":
        # Unstructured grid: only last dimension
        npoints = src_dims[-1]
        extra_dims = src_dims[:-1]
    else:
        raise ValueError(f"Unknown grid type: {src_grid_type}")

    # Prepare an array to hold the regridded data
    dst_grid_dims = xr.open_dataset(wgt_filename)['dst_grid_dims'].values
    regridded_shape = extra_dims + tuple(dst_grid_dims)
    regridded_data = np.zeros(regridded_shape)

    # Iterate over all combinations of the extra dimensions
    for idx in np.ndindex(*extra_dims):

        # Extract horizontal slice
        slice_2d = src_data[idx]

        # Apply the regridding to this slice
        if src_grid_type == "structured":
            regridded_slice = remap_with_weights(slice_2d[np.newaxis, :, :], sparse_map, dst_grid_dims, src_grid_type, dst_grid_type, **kwargs)
        elif src_grid_type == "unstructured":
            regridded_slice = remap_with_weights(slice_2d[np.newaxis, :]   , sparse_map, dst_grid_dims, src_grid_type, dst_grid_type, **kwargs)
        else:
            raise ValueError(f"Unknown grid type: {src_grid_type}")

        # Store the regridded slice in the appropriate location in the output array
        regridded_data[idx] = regridded_slice

    elapsed_time = time.time() - start_time
    print(f"Regridding completed in {elapsed_time:.2f} seconds.")

    if dst_grid_type == "structured":
        # Swap the last two axes to go from lon, lat axes to lat, lon.
        regridded_data = regridded_data.swapaxes(-1, -2)
        # Reshape dstlat and dstlon arrays to oneD lat/lon
        # This probably breaks for curvilinear grids
        dstlat = np.reshape(dstlat.values, dst_grid_dims, order="F")
        dstlon = np.reshape(dstlon.values, dst_grid_dims, order="F")
        dstlat = dstlat[0,:]
        dstlon = dstlon[:,0]
    else:
        # Otherwise, just return 1D ncol
        dstlat = dstlat.values
        dstlon = dstlon.values

    return regridded_data, dstlat, dstlon


def remap_all(data_in, wgt_filename, dycore='se'):

    allowable_interp_vars = ['ps', 't', 'u', 'v', 'q', 'cldice', 'cldliq', 'z', 'theta', 'rho', 'w', 'phis', 'ts']

    data_out = {}

    # Loop ovr the keys in data_in and interpolate if the key is in allowable_interp_vars
    for key in data_in:
        if key in allowable_interp_vars:
            data_out[key], _, _ = remap_with_weights_wrapper(data_in[key], wgt_filename)
        else:
            data_out[key] = data_in[key]

    data_out['ps'], data_out['lat'], data_out['lon'] = remap_with_weights_wrapper(data_in['ps'], wgt_filename)

    return data_out
