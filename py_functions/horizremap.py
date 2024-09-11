import numpy as np
import xarray as xr
import time
import scipy.sparse as sp
import logging
from scipy.interpolate import RegularGridInterpolator
logger = logging.getLogger(__name__)


def remap_with_weights(src_data, sparse_map, dst_grid_dims, src_grid_type, dst_grid_type, dest_frac=None, dest_frac_thresh=0.999, **kwargs):
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

        if isinstance(slice_2d, xr.DataArray):
            logging.info("converting slice_2d to numpy")
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
    logging.info(f"Regridding completed in {elapsed_time:.2f} seconds.")

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

    return regridded_data, dstlat, dstlon


def remap_all(data_in, wgt_filename, dycore='se'):

    allowable_interp_vars = ['ps', 't', 'u', 'v', 'q', 'cldice', 'cldliq', 'z', 'theta', 'rho', 'w', 'phis', 'ts']

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