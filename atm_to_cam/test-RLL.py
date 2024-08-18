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

from horizremap import *

# data = xr.open_dataset('./test_files/ne30np4_nc3000_Co060_Fi001_PF_nullRR_Nsw042_20171020.nc')
# wgt_filename = './test_files/ne30_to_1x1_patch.nc'
# t_fv = data["PHIS"].values
# t_cam, selat, selon = remap_with_weights_wrapper(t_fv, wgt_filename)
# print(t_cam.shape)
# print(selat.shape)
# print(selon.shape)
# # Create an xarray.Dataset
# ds = xr.Dataset(
#     {
#         "t_cam": (["lat", "lon"], t_cam.astype(np.float32))
#     },
#     coords={
#         "lat": (["lat"], selat.astype(np.float32)),
#         "lon": (["lon"], selon.astype(np.float32)),
#     }
# )
# output_filename = "ncol_to_RLL.nc"
# ds.to_netcdf(output_filename)
#


data = xr.open_dataset('./test_files/ds633.0/e5.oper.invariant/197901/e5.oper.invariant.128_129_z.ll025sc.1979010100_1979010100.nc')
wgt_filename = './test_files/map_era5_0.25x0.25_TO_ne30_patc.nc'
t_fv = data["Z"].isel(time=0).values
t_cam, selat, selon = remap_with_weights_wrapper(t_fv, wgt_filename)
print(t_cam.shape)
print(selat.shape)
print(selon.shape)
# Create an xarray.Dataset
ds = xr.Dataset(
    {
        "t_cam": (["ncol"], t_cam.astype(np.float32))
    },
    coords={
        "lat": (["ncol"], selat.astype(np.float32)),
        "lon": (["ncol"], selon.astype(np.float32)),
    }
)
output_filename = "RLL_to_ncol.nc"
ds.to_netcdf(output_filename)