import numpy as np

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
