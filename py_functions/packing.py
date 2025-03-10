import numpy as np
from scipy.interpolate import RegularGridInterpolator

def latlon_to_ncol(var_in):
    vardims = var_in.shape
    dims = len(vardims)

    if dims == 2:
        var_out = var_in.flatten(order='F')
    elif dims == 3:
        nlev = vardims[0]
        nlat = vardims[1]
        nlon = vardims[2]
        ncol = nlat * nlon
        print(f"unpacking -> nlev: {nlev}    nlat: {nlat}    nlon: {nlon}    ncol: {ncol}")
        var_out = np.empty((nlev, ncol), dtype=var_in.dtype)
        for ii in range(nlev):
            var_out[ii, :] = var_in[ii, :, :].flatten(order='F')
    else:
        raise ValueError(f"{vardims} dims not supported")

    return var_out



def ncol_to_latlon(var_out, nlat, nlon):
    vardims = var_out.shape

    if len(vardims) == 1:
        var_in = var_out.reshape((nlat, nlon), order='F')
    elif len(vardims) == 2:
        print(vardims)
        nlev = vardims[0]
        ncol = vardims[1]
        print(f"repacking -> nlev: {nlev}    nlat: {nlat}    nlon: {nlon}    ncol: {ncol}")
        var_in = np.empty((nlev, nlat, nlon), dtype=var_out.dtype)
        for ii in range(nlev):
            var_in[ii, :, :] = var_out[ii, :].reshape((nlat, nlon), order='F')
    else:
        raise ValueError(f"{vardims} dims not supported")

    return var_in


def repack_fv(data_horiz, fv_dims):

    nlat, nlon = fv_dims

    print("Repacking FV variables...")

    keys_to_repack = ['ps', 't', 'q', 'u', 'v', 'cldice', 'cldliq', 'ts', 'phis', 'correct_or_not']

    for key in keys_to_repack:
        if key in data_horiz:
            data_horiz[key] = ncol_to_latlon(data_horiz[key], nlat, nlon)

    return data_horiz




def unpack_fv(data_horiz):

    print("Unpacking FV variables...")

    keys_to_unpack = ['ps', 't', 'q', 'u', 'v', 'cldice', 'cldliq', 'ts', 'phis']

    for key in keys_to_unpack:
        if key in data_horiz:
            data_horiz[key] = latlon_to_ncol(data_horiz[key])

    return data_horiz




def initialize_fv_grid(nfvlat, nfvlon):
    """Initialize the FV grid and derived quantities."""
    nfvslat = nfvlat - 1
    nfvslon = nfvlon

    fvlat = np.linspace(-90, 90, nfvlat)
    delfvlon = 360.0 / nfvlon
    fvlon = np.linspace(0, 360.0 - delfvlon, nfvlon)

    fvslat = (fvlat[:-1] + fvlat[1:]) / 2.0
    fvslon = fvlon - (delfvlon / 2.0)

    return fvlat, fvlon, fvslat, fvslon


def interpolate_uv_to_slat_slon(data_horiz, numlevels, fvlat, fvlon, fvslat, fvslon):

    nfvlon=len(fvlon)
    nfvlat=len(fvlat)
    nfvslat=len(fvslat)
    nfvslon=len(fvslon)

    slat_lon_mesh = np.array(np.meshgrid(fvslat, fvlon)).T.reshape(-1, 2)
    lat_slon_mesh = np.array(np.meshgrid(fvlat, fvslon)).T.reshape(-1, 2)

    # Create a cyclic ghost point for slon interpolation
    ghost_lon = fvlon[-1] - 360.0
    tmp_fvlon = np.hstack(([ghost_lon], fvlon))

    us = np.zeros((numlevels, nfvslat, nfvlon))
    vs = np.zeros((numlevels, nfvlat, nfvlon))

    for level in range(numlevels):
        # Extend v to match ghost point since this var is interpolated in slon
        v_fv_extended = np.hstack((data_horiz['v'][level, :, -1][:, np.newaxis], data_horiz['v'][level, :, :]))
        u_interpolator = RegularGridInterpolator((fvlat, fvlon), data_horiz['u'][level, :, :], method='linear')
        v_interpolator = RegularGridInterpolator((fvlat, tmp_fvlon), v_fv_extended, method='linear')

        us[level, :, :] = u_interpolator(slat_lon_mesh).reshape(nfvslat, nfvlon)
        vs[level, :, :] = v_interpolator(lat_slon_mesh).reshape(nfvlat, nfvlon)

    del ghost_lon, tmp_fvlon, v_interpolator, u_interpolator, v_fv_extended

    return us, vs