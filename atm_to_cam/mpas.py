import numpy as np
from scipy.interpolate import interp1d
from numba import jit

# def z_to_z_interp(theta_fv, z_fv, thisCol, extrapLow=False, extrapHigh=False):
#     """
#     Interpolates data from one vertical grid to another.
#
#     Parameters:
#     -----------
#     theta_fv : numpy.ndarray
#         The input field to be interpolated (e.g., potential temperature).
#     z_fv : numpy.ndarray
#         The vertical coordinate of the input field (e.g., heights or pressure levels).
#     thisCol : numpy.ndarray
#         The target vertical coordinate to interpolate onto.
#     extrapLow : bool, optional
#         Whether to extrapolate below the lowest valid level. Default is False.
#     extrapHigh : bool, optional
#         Whether to extrapolate above the highest valid level. Default is False.
#
#     Returns:
#     --------
#     t_wrf : numpy.ndarray
#         The interpolated field.
#     """
#     nlev = len(thisCol)
#
#     # Interpolate using 1D interpolation
#     f_interp = interp1d(z_fv, theta_fv, bounds_error=False, fill_value=np.nan)
#     t_wrf = f_interp(thisCol)
#
#     # Check if there are any missing values and handle extrapolation
#     if np.any(np.isnan(t_wrf)):
#         ixvalid = np.where(~np.isnan(t_wrf))[0]
#
#         # Extrapolate below the lowest valid level if necessary
#         if np.isnan(t_wrf[0]):
#             lowest_valid = np.min(ixvalid)
#             if extrapLow:
#                 deriv = (t_wrf[lowest_valid] - t_wrf[lowest_valid+1]) / (thisCol[lowest_valid] - thisCol[lowest_valid+1])
#             else:
#                 deriv = 0.0
#             for ff in range(0, lowest_valid):
#                 t_wrf[ff] = t_wrf[lowest_valid] + deriv * (thisCol[ff] - thisCol[lowest_valid])
#
#         # Extrapolate above the highest valid level if necessary
#         if np.isnan(t_wrf[-1]):
#             highest_valid = np.max(ixvalid)
#             if extrapHigh:
#                 deriv = (t_wrf[highest_valid-1] - t_wrf[highest_valid]) / (thisCol[highest_valid-1] - thisCol[highest_valid])
#             else:
#                 deriv = 0.0
#             for ff in range(highest_valid+1, nlev):
#                 t_wrf[ff] = t_wrf[highest_valid] + deriv * (thisCol[ff] - thisCol[highest_valid])
#
#     return t_wrf
#

@jit(nopython=True)
def z_to_z_interp(theta_fv, z_fv, thisCol, extrapLow=False, extrapHigh=False):
    nlev = len(thisCol)
    t_wrf = np.zeros(nlev)

    for i in range(nlev):
        if thisCol[i] <= z_fv[0]:
            t_wrf[i] = theta_fv[0]
        elif thisCol[i] >= z_fv[-1]:
            t_wrf[i] = theta_fv[-1]
        else:
            for j in range(len(z_fv) - 1):
                if z_fv[j] <= thisCol[i] <= z_fv[j + 1]:
                    t_wrf[i] = theta_fv[j] + (theta_fv[j + 1] - theta_fv[j]) * (thisCol[i] - z_fv[j]) / (z_fv[j + 1] - z_fv[j])
                    break

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

@jit(nopython=True)
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
        #if iEdge % 10000 == 0:
        #    print(f"... UV_CELL_TO_EDGE: {100. * iEdge / (nEdges - 1):.2f}%")

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




@jit(nopython=True)
def interpolate_mpas_columns_wrapper(mpas_ncell, mpas_nlev, mpas_nlevi, mpas_z, t_fv, z_fv, theta_fv, rho_fv, w_fv, q_fv, u_fv, v_fv, mpas_as_cam):
    t_wrf = np.zeros((mpas_nlev, mpas_ncell))
    theta_wrf = np.zeros((mpas_nlev, mpas_ncell))
    rho_wrf = np.zeros((mpas_nlev, mpas_ncell))
    w_wrf = np.zeros((mpas_nlevi, mpas_ncell))
    q_wrf = np.zeros((mpas_nlev, mpas_ncell))
    u_wrf = np.zeros((mpas_nlev, mpas_ncell))
    v_wrf = np.zeros((mpas_nlev, mpas_ncell))

    for ix in range(mpas_ncell):
        #if ix % 10000 == 0:
        #    print(f"... MPAS_VERT_INTERP: {100. * ix / (mpas_ncell - 1):.2f}%")

        t_wrf[:, ix], theta_wrf[:, ix], rho_wrf[:, ix], w_wrf[:, ix], q_wrf[:, ix], u_wrf[:, ix], v_wrf[:, ix] = interpolate_single_mpas_column(
            ix, mpas_nlev, mpas_nlevi, mpas_z, t_fv, z_fv, theta_fv, rho_fv, w_fv, q_fv, u_fv, v_fv, mpas_as_cam
        )

    return t_wrf, theta_wrf, rho_wrf, w_wrf, q_wrf, u_wrf, v_wrf




@jit(nopython=True)
def interpolate_single_mpas_column(ix, mpas_nlev, mpas_nlevi, mpas_z, t_fv, z_fv, theta_fv, rho_fv, w_fv, q_fv, u_fv, v_fv, mpas_as_cam):
    zmid = (mpas_z[ix, 1:mpas_nlevi] + mpas_z[ix, 0:mpas_nlevi-1]) / 2.0
    zint = mpas_z[ix, :]

    t_wrf_col = np.zeros(mpas_nlev)
    theta_wrf_col = np.zeros(mpas_nlev)
    rho_wrf_col = np.zeros(mpas_nlev)
    w_wrf_col = np.zeros(mpas_nlevi)
    q_wrf_col = np.zeros(mpas_nlev)
    u_wrf_col = np.zeros(mpas_nlev)
    v_wrf_col = np.zeros(mpas_nlev)

    if mpas_as_cam:
        t_wrf_col = z_to_z_interp(t_fv[::-1, ix], z_fv[::-1, ix], zmid, extrapLow=True, extrapHigh=True)

    theta_wrf_col = z_to_z_interp(theta_fv[::-1, ix], z_fv[::-1, ix], zmid, extrapLow=True, extrapHigh=True)
    rho_wrf_col = z_to_z_interp(rho_fv[::-1, ix], z_fv[::-1, ix], zmid, extrapLow=True, extrapHigh=True)
    w_wrf_col = z_to_z_interp(w_fv[::-1, ix], z_fv[::-1, ix], zint, extrapLow=True, extrapHigh=True)
    q_wrf_col = z_to_z_interp(q_fv[::-1, ix], z_fv[::-1, ix], zmid, extrapLow=True, extrapHigh=True)

    # u and v we don't extrapolate aloft to prevent wind speeds from getting too high
    u_wrf_col = z_to_z_interp(u_fv[::-1, ix], z_fv[::-1, ix], zmid, extrapLow=True, extrapHigh=False)
    v_wrf_col = z_to_z_interp(v_fv[::-1, ix], z_fv[::-1, ix], zmid, extrapLow=True, extrapHigh=False)

    return t_wrf_col, theta_wrf_col, rho_wrf_col, w_wrf_col, q_wrf_col, u_wrf_col, v_wrf_col