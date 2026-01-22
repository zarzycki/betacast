import numpy as np
import xarray as xr
import time
import scipy.sparse as sp
import logging
from scipy.interpolate import RegularGridInterpolator
logger = logging.getLogger(__name__)

_HORIZ_REMAP_VARS = [
    'ps', 't', 'u', 'v', 'q', 'cldice', 'cldliq', 'numliq', 'numice', 'numcld',
    'z', 'theta', 'rho', 'w', 'phis', 'ts', 'o3'
]

def remap_with_weights(src_data, sparse_map, dst_grid_dims, src_grid_type, dst_grid_type, dest_frac=None, dest_frac_thresh=0.999, **kwargs):
    """
    Regrid data using ESMF weights.

    Parameters:
    - src_data: 2D numpy array of the source data.
    - sparse_map: Sparse matrix containing the regridding weights.
    - dst_grid_dims: Array specifying the destination grid dimensions.
    - src_grid_type: String indicating source grid type ("structured" or "unstructured").
    - dst_grid_type: String indicating destination grid type ("structured" or "unstructured").
    - dest_frac: Optional array of destination cell fractions participating in regridding.
    - dest_frac_thresh: Threshold for dest_frac; cells below this are set to NaN (default 0.999).
    - kwargs: Additional keyword arguments (currently unused).

    Returns:
    - data_out: Regridded data as a 2D numpy array with shape dst_grid_dims.
    """

    if src_grid_type == "structured":
        FLATTENTYPE = "C"
    elif src_grid_type == "unstructured":
        FLATTENTYPE = "F"
    else:
        raise ValueError(f"Unknown grid type: {src_grid_type}")

    # Flatten the source data to create a vector of length n_a
    # It appears ESMF uses "C" ordering, although we may need to use "F" if src is unstructured?
    src_vector = src_data.flatten(order=FLATTENTYPE)

    # Apply map onto the src column vector length n_a to compute vector length n_b
    field_target = sparse_map @ src_vector

    # If we've passed in fraction of dest cells participating in regridding, set
    # all points < dest_frac_thresh to nan so we don't get a lot of 0's from
    # the matrix solve
    if dest_frac is not None:
        field_target = np.where(dest_frac > dest_frac_thresh, field_target, np.nan)

    # Reshape 1-D vector returned to dst_grid_dims
    data_out = np.reshape(field_target, dst_grid_dims, order="F")

    return data_out


def remap_with_weights_wrapper(src_data, wgt_filename, return_xarray=False, **kwargs):
    """
    Wrapper to handle regridding of multi-dimensional data using ESMF weights.

    Parameters:
    - src_data: Multi-dimensional numpy array or xarray DataArray of the source data.
    - wgt_filename: Path to the ESMF weight file.
    - return_xarray: Boolean, if True, returns an xarray DataArray with coordinates
      and dimensions (default False).
    - kwargs: Additional keyword arguments passed to remap_with_weights.

    Returns:
    - regridded_data: Multi-dimensional numpy array or xarray DataArray after regridding.
    - dstlat: 1D array of destination latitude values.
    - dstlon: 1D array of destination longitude values.
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

    logging.info(f"Src grid dims: {src_grid_dims}, dst grid dims: {dst_grid_dims}")
    src_grid_type = "structured" if src_grid_dims.size == 2 else "unstructured"
    dst_grid_type = "structured" if dst_grid_dims.size == 2 else "unstructured"
    logging.info(f"Src grid type: {src_grid_type}, Dst grid type: {dst_grid_type}")

    n_a = xwgt['n_a'].size  # col dimension
    n_b = xwgt['n_b'].size  # row dimension
    n_s = xwgt['n_s'].size  # nnz dimension

    logging.info("Map contains {0} rows, {1} cols and {2} nnz values".format(n_b, n_a, n_s))

    rows = xwgt['row'][:] - 1  # row indices (1-based to 0-based)
    cols = xwgt['col'][:] - 1  # col indices (1-based to 0-based)
    nnzvals = xwgt['S'][:]  # nnz map values

    frac_b = xwgt['frac_b'][:]

    # Create sparse matrix map
    sparse_map = sp.coo_matrix((nnzvals, (rows, cols)), shape=(n_b, n_a))

    # Check if src_data is xarray DataArray and extract dimensions and coordinates
    is_xarray_input = isinstance(src_data, xr.DataArray)
    src_dims = src_data.shape
    src_coords = src_data.coords if is_xarray_input else {}

    # Get src spatial information and then back out non-spatial dimensions
    if src_grid_type == "structured":
        # Structured grid: last two dimensions are nlat, nlon
        nlat, nlon = src_dims[-2], src_dims[-1]
        non_spatial_dims = src_data.dims[:-2] if is_xarray_input else range(len(src_dims) - 2)
    elif src_grid_type == "unstructured":
        # Unstructured grid: last dimension is ncol
        npoints = src_dims[-1]
        non_spatial_dims = src_data.dims[:-1] if is_xarray_input else range(len(src_dims) - 1)
    else:
        raise ValueError(f"Unknown grid type: {src_grid_type}")

    # If structured, cut last two dims off -- if unstructured, cut last one off.
    # This gives aux dims like time, lev, etc.
    # Then add the destination dimensions on the back. This way you get merger of
    # src aux dims and dest spatial dims
    regridded_shape = src_dims[:-2 if src_grid_type == "structured" else -1] + tuple(dst_grid_dims)
    regridded_data = np.zeros(regridded_shape)

    logging.debug(f"is_xarray_input: {is_xarray_input}")
    logging.debug(f"non_spatial_dims: {non_spatial_dims}")
    logging.debug(f"regridded_shape: {tuple(map(int, regridded_shape))}")

    # Use .sizes for xarray DataArray, and .shape for numpy arrays
    if is_xarray_input:
        dim_sizes = [src_data.sizes[dim] for dim in non_spatial_dims]
    else:
        dim_sizes = [src_data.shape[i] for i in range(len(non_spatial_dims))]

    # Calculate total number of iterations for the loop
    total_iterations = np.prod(dim_sizes)

    # Iterate over all combinations of the non-spatial dimensions
    for iteration, idx in enumerate(np.ndindex(*dim_sizes), start=1):
        # Logging message with detailed information
        logging.debug(f"Horizontally remapping index {idx} ({iteration}/{total_iterations}) "
                     f"of shape {tuple(dim_sizes)} with destination grid {tuple(dst_grid_dims)}")

        # Extract horizontal slice
        slice_2d = src_data[idx] if is_xarray_input else src_data[idx]

        if isinstance(slice_2d, xr.DataArray):
            logging.debug("Converting slice_2d to numpy")
            slice_2d = slice_2d.values  # Convert to numpy.ndarray

        logging.debug(f"slice_2d shape before newaxis: {slice_2d.shape} and type: {type(slice_2d)}")

        # Apply the regridding to this slice
        if src_grid_type == "structured":
            regridded_slice = remap_with_weights(slice_2d[np.newaxis, :, :], sparse_map, dst_grid_dims, src_grid_type, dst_grid_type, dest_frac=frac_b, **kwargs)
        elif src_grid_type == "unstructured":
            regridded_slice = remap_with_weights(slice_2d[np.newaxis, :], sparse_map, dst_grid_dims, src_grid_type, dst_grid_type, dest_frac=frac_b, **kwargs)
        else:
            raise ValueError(f"Unknown grid type: {src_grid_type}")

        # Store the regridded slice in the appropriate location in the output array
        regridded_data[idx] = regridded_slice

    elapsed_time = time.time() - start_time
    logging.info(f"Horizointal regridding completed in {elapsed_time:.2f} seconds.")

    if dst_grid_type == "structured":
        # Swap the last two axes to go from lon, lat axes to lat, lon.
        regridded_data = regridded_data.swapaxes(-1, -2)
        # Reshape dstlat and dstlon arrays to oneD lat/lon
        # This probably breaks for curvilinear grids
        dstlat = np.reshape(dstlat.values, dst_grid_dims, order="F")
        dstlon = np.reshape(dstlon.values, dst_grid_dims, order="F")
        dstlat = dstlat[0, :]
        dstlon = dstlon[:, 0]
    else:
        # Otherwise, just return 1D ncol
        dstlat = dstlat.values
        dstlon = dstlon.values

    if return_xarray:
        # Construct xarray DataArray with appropriate coordinates and dimensions
        coords = {}
        dims = list(non_spatial_dims)
        if dst_grid_type == "structured":
            coords.update({dim: src_coords[dim].values for dim in non_spatial_dims})
            coords['lat'] = dstlat
            coords['lon'] = dstlon
            dims.extend(['lat', 'lon'])
        else:
            coords.update({dim: src_coords[dim].values for dim in non_spatial_dims})
            coords['ncol'] = np.arange(dst_grid_dims[0])
            dims.append('ncol')

        logging.debug(f"return_xarray: regridded_data shape: {regridded_data.shape}")
        logging.debug(f"return_xarray: dims: {dims}")

        regridded_data = xr.DataArray(regridded_data, dims=dims, coords=coords)

    return regridded_data, dstlat, dstlon


def remap_all(data_in, wgt_filename, dycore='se'):
    """
    Regrid all applicable variables in a data dictionary.

    Parameters:
    - data_in: Dictionary of variable names to numpy arrays or xarray DataArrays.
    - wgt_filename: Path to the ESMF weight file.
    - dycore: Dynamical core identifier (default 'se'); currently unused.

    Returns:
    - data_out: Dictionary containing regridded data. Variables in _HORIZ_REMAP_VARS
      are interpolated; others are passed through unchanged. The 'ps' variable is
      always regridded, and 'lat'/'lon' coordinates are added from its output.
    """

    allowable_interp_vars = _HORIZ_REMAP_VARS

    data_out = {}

    # Loop ovr the keys in data_in and interpolate if the key is in allowable_interp_vars
    for key in data_in:
        logging.info(key)
        if key in allowable_interp_vars:
            data_out[key], _, _ = remap_with_weights_wrapper(data_in[key], wgt_filename)
        else:
            data_out[key] = data_in[key]

    data_out['ps'], data_out['lat'], data_out['lon'] = remap_with_weights_wrapper(data_in['ps'], wgt_filename)

    return data_out


def uv_cell_to_edge(uZonal, uMerid, nlev, lonEdge, latEdge, lonCell, latCell, edgeNormalVectors, cellsOnEdge):
    """
    Converts zonal and meridional wind components from cell centers to edge centers.

    Parameters:
    -----------
    uZonal : numpy.ndarray
        Zonal wind component at cell centers, shape (nlev, nCells).
    uMerid : numpy.ndarray
        Meridional wind component at cell centers, shape (nlev, nCells).
    nlev : int
        Number of vertical levels.
    lonEdge : numpy.ndarray
        Longitudes of edge centers, shape (nEdges,).
    latEdge : numpy.ndarray
        Latitudes of edge centers, shape (nEdges,).
    lonCell : numpy.ndarray
        Longitudes of cell centers, shape (nCells,).
    latCell : numpy.ndarray
        Latitudes of cell centers, shape (nCells,).
    edgeNormalVectors : numpy.ndarray
        Edge normal vectors, shape (nEdges, 3).
    cellsOnEdge : numpy.ndarray
        Cell indices on edges, shape (nEdges, 2).

    Returns:
    --------
    uNormal : numpy.ndarray
        Normal wind component at edge centers, shape (nlev, nEdges).
    """

    nCells = len(lonCell)
    nEdges = len(lonEdge)

    east = np.zeros((3, nCells))
    north = np.zeros((3, nCells))

    for iCell in range(nCells):
        east[0, iCell] = -np.sin(lonCell[iCell])
        east[1, iCell] = np.cos(lonCell[iCell])
        east[2, iCell] = 0.0

        # Normalize
        east[:, iCell] /= np.sqrt(np.sum(east[:, iCell] ** 2))

        north[0, iCell] = -np.cos(lonCell[iCell]) * np.sin(latCell[iCell])
        north[1, iCell] = -np.sin(lonCell[iCell]) * np.sin(latCell[iCell])
        north[2, iCell] = np.cos(latCell[iCell])

        # Normalize
        north[:, iCell] /= np.sqrt(np.sum(north[:, iCell] ** 2))

    uNormal = np.zeros((nlev, nEdges), dtype=uZonal.dtype)

    for iEdge in range(nEdges):
        if iEdge % 10000 == 0:
            logging.info(f"... UV_CELL_TO_EDGE: {100. * iEdge / (nEdges - 1):.2f}%")

        cell1 = cellsOnEdge[iEdge, 0] - 1
        cell2 = cellsOnEdge[iEdge, 1] - 1

        uNormal[:, iEdge] = (
            uZonal[:, cell1] * 0.5 * (edgeNormalVectors[iEdge, 0] * east[0, cell1] +
                                      edgeNormalVectors[iEdge, 1] * east[1, cell1] +
                                      edgeNormalVectors[iEdge, 2] * east[2, cell1])
            +
            uMerid[:, cell1] * 0.5 * (edgeNormalVectors[iEdge, 0] * north[0, cell1] +
                                      edgeNormalVectors[iEdge, 1] * north[1, cell1] +
                                      edgeNormalVectors[iEdge, 2] * north[2, cell1])
            +
            uZonal[:, cell2] * 0.5 * (edgeNormalVectors[iEdge, 0] * east[0, cell2] +
                                      edgeNormalVectors[iEdge, 1] * east[1, cell2] +
                                      edgeNormalVectors[iEdge, 2] * east[2, cell2])
            +
            uMerid[:, cell2] * 0.5 * (edgeNormalVectors[iEdge, 0] * north[0, cell2] +
                                      edgeNormalVectors[iEdge, 1] * north[1, cell2] +
                                      edgeNormalVectors[iEdge, 2] * north[2, cell2])
        )

    return uNormal


def linint2(xi, yi, fi, fiCyclicX, xo, yo):
    """
    Bilinear interpolation from one rectilinear grid to another.

    Parameters:
    - xi: 1D or 2D array of X coordinates of the input data.
    - yi: 1D or 2D array of Y coordinates of the input data.
    - fi: 2D or 3D array of input data to be interpolated.
    - fiCyclicX: Boolean indicating if the X dimension is cyclic.
    - xo: 1D array of X coordinates of the output data.
    - yo: 1D array of Y coordinates of the output data.

    Returns:
    - Interpolated array of the same dimensions as fi, except the last two dimensions.
    """

    # Handle cyclic conditions
    if fiCyclicX:
        # Extend the fi array and xi array for cyclic boundary
        xi_extended = np.concatenate((xi, xi[:1] + 360.0), axis=-1)
        fi_extended = np.concatenate((fi, fi[..., :1]), axis=-1)
        xi = xi_extended
        fi = fi_extended

    # Create the interpolator object
    interpolator = RegularGridInterpolator((yi, xi), fi, method='linear', bounds_error=False, fill_value=np.nan)

    # Create a meshgrid for target coordinates
    xo_mesh, yo_mesh = np.meshgrid(xo, yo, indexing='xy')

    # Prepare the points for interpolation
    target_points = np.vstack((yo_mesh.ravel(), xo_mesh.ravel())).T

    # Perform the interpolation
    interpolated_data = interpolator(target_points)

    # Reshape the interpolated data to match the target grid shape
    interpolated_data = interpolated_data.reshape(yo_mesh.shape)

    return interpolated_data