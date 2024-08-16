import os
import numpy as np
import xarray as xr
import datetime
import argparse
import sys
import glob
from scipy.interpolate import interp1d
from numba import jit
from tqdm import tqdm
import time
import scipy.sparse as sp

from pyfuncs import *

data = xr.open_dataset('ncl_era5_regrid.nc')
wgt_filename = './test_files/ne30_to_1x1_patch.nc'

t_fv = data["t_fv"].values

t_cam, selat, selon = esmf_regrid_with_weights_wrapper(t_fv, wgt_filename)

print(t_cam.shape)
print(selat.shape)
print(selon.shape)

# Create an xarray.Dataset
ds = xr.Dataset(
    {
        "t_cam": (["level", "lat", "lon"], t_cam.astype(np.float32))
    },
    coords={
        "lat": (["lat"], selat.astype(np.float32)),
        "lon": (["lon"], selon.astype(np.float32)),
    }
)
output_filename = "tmp.nc"
ds.to_netcdf(output_filename)