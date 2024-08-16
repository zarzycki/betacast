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
import vertremap
import horizremap

def main():

    # Check environment variable
    BETACAST = os.getenv("BETACAST")
    if BETACAST is None:
        print("We are running local-only atm_to_cam. Set export BETACAST env to run elsewhere")
        local_only = True
        PATHTOHERE = "./"
    else:
        print("Not local only!")
        local_only = False
        PATHTOHERE = os.path.join(BETACAST, "atm_to_cam")

    define_constants()

    args = parse_args()

    # Example of using the arguments
    print(f"Using datasource: {args.datasource}")
    print(f"Number of levels: {args.numlevels}")
    print(f"Initialization date: {args.YYYYMMDDHH}")
    print(f"Data filename: {args.data_filename}")
    print(f"Weight filename: {args.wgt_filename}")
    print(f"SE Inic output file: {args.se_inic}")
    print(f"Dycore: {args.dycore}")
    print(f"RDADIR: {args.RDADIR}")
    print(f"MPAS as CAM: {args.mpas_as_cam}")
    print(f"Compress file: {args.compress_file}")
    print(f"Write floats: {args.write_floats}")
    print(f"Add cloud vars: {args.add_cloud_vars}")
    print(f"Adjust config: {args.adjust_config}")
    print(f"Model topo file: {args.model_topo_file}")
    print(f"MOD in topo: {args.mod_in_topo}")
    print(f"MOD remap file: {args.mod_remap_file}")

    # Assume these are parsed from command-line arguments
    dycore = args.dycore
    numlevels = args.numlevels
    model_topo_file = args.model_topo_file
    YYYYMMDDHH = args.YYYYMMDDHH
    data_filename = args.data_filename
    wgt_filename = args.wgt_filename
    mpasfile = args.mpasfile
    datasource = args.datasource
    model_topo_file = args.model_topo_file
    mod_in_topo = args.mod_in_topo
    mod_remap_file = args.mod_remap_file
    write_floats = args.write_floats
    RDADIR = args.RDADIR

    # Check if we are using MPAS and the model topo file is provided and exists
    if dycore == "mpas" and model_topo_file:
        if os.path.isfile(model_topo_file):
            print(f"Setting mpasfile to {model_topo_file}")
            mpasfile = model_topo_file
        else:
            print(f"Model topo file {model_topo_file} not found, exiting...")
            #sys.exit(1)

    # Check if mpasfile is set if dycore is "mpas"
    if dycore == "mpas" and mpasfile is None:
        print("No mpasfile passed in but mpas dycore selected.")
        print("MPAS init file is needed to mimic + get zcoord from. Exiting...")
        sys.exit(1)

    # Error check if CAM data is used as input
    if datasource == "CAM":
        if not mod_in_topo:
            print("Using CAM data, but no mod_in_topo passed in, exiting...")
            sys.exit(1)
        if not mod_remap_file:
            print("Using CAM data, but no mod_remap_file passed in, exiting...")
            sys.exit(1)

    # Toggle whether the output streams will be floats or doubles
    write_type = "float" if write_floats else "double"
    print(f"Output type set to: {write_type}")

    # ===== Getting date from YYYYMMDDHH
    dtime_map = [4, 2, 2, 2]
    yearstr, monthstr, daystr, cyclestr = split_by_lengths(str(YYYYMMDDHH), dtime_map)

    print(f"Regridding analysis from: {yearstr} {monthstr} {daystr} {cyclestr}Z")

    print("---------------------------------------------------------")
    print(f"Using this file: {data_filename}")
    print(f"Using this remap: {wgt_filename}")

    hya, hyb, hyai, hybi, lev, ilev = load_cam_levels(PATHTOHERE, numlevels)

    # Load the data file
    grb_file = load_grb_file(data_filename)

    # Process the grb_file to extract lat/lon/lev
    grblat, grblon, grblev = process_grb_file(grb_file, datasource, RDADIR, yearstr, monthstr, daystr, mod_remap_file)

    print("Level information")
    print(f"Number: {len(grblev)}")
    print(f"Max: {np.max(grblev)}")
    print(f"Min: {np.min(grblev)}")

    nlat = len(grblat)
    nlon = len(grblon)
    print(f"nlat: {nlat}, nlon: {nlon}")

    data_vars = load_era5rda_data(RDADIR, yearstr, monthstr, daystr, cyclestr, dycore)

    # Access the specific DataArray for 'SP' from the Dataset
    ps = data_vars['ps']
    t_gfs = data_vars['t_gfs']
    u_gfs = data_vars['u_gfs']
    v_gfs = data_vars['v_gfs']
    q_gfs = data_vars['q_gfs']
    cldliq_gfs = data_vars['cldliq_gfs']
    cldice_gfs = data_vars['cldice_gfs']

    # For MPAS dycore, access W and Z from their respective datasets
    if dycore == "mpas":
        w_gfs = data_vars['w_gfs']
        z_gfs = data_vars['z_gfs']

    # Now print the shapes as you intended
    print(f"Surface pressure shape: {ps.shape}")
    print(f"Temperature shape: {t_gfs.shape}")
    print(f"U wind shape: {u_gfs.shape}")
    print(f"V wind shape: {v_gfs.shape}")
    print(f"Specific humidity shape: {q_gfs.shape}")
    print(f"Cloud liquid water content shape: {cldliq_gfs.shape}")
    print(f"Cloud ice water content shape: {cldice_gfs.shape}")

    if dycore == "mpas":
        print(f"W shape: {w_gfs.shape}")
        print(f"Z shape: {z_gfs.shape}")

    print("="*65)
    print("************ NATIVE DATA")
    print(f"Max T: {t_gfs.max()}   min T: {t_gfs.min()}")
    print(f"Max U: {u_gfs.max()}   min U: {u_gfs.min()}")
    print(f"Max V: {v_gfs.max()}   min V: {v_gfs.min()}")
    print(f"Max Q: {q_gfs.max()}   min Q: {q_gfs.min()}")
    print(f"Max PS: {ps.max()}   min PS: {ps.min()}")
    print(f"Max CLDICE: {cldice_gfs.max()}   min CLDICE: {cldice_gfs.min()}")
    print(f"Max CLDLIQ: {cldliq_gfs.max()}   min CLDLIQ: {cldliq_gfs.min()}")
    print("="*65)

    t_cam = vertremap.interpolate_to_hybrid_levels(grblev, t_gfs, ps, hya, hyb)
    u_cam = vertremap.interpolate_to_hybrid_levels(grblev, u_gfs, ps, hya, hyb)
    v_cam = vertremap.interpolate_to_hybrid_levels(grblev, v_gfs, ps, hya, hyb)
    q_cam = vertremap.interpolate_to_hybrid_levels(grblev, q_gfs, ps, hya, hyb)
    cldice_cam = vertremap.interpolate_to_hybrid_levels(grblev, cldice_gfs, ps, hya, hyb)
    cldliq_cam = vertremap.interpolate_to_hybrid_levels(grblev, cldliq_gfs, ps, hya, hyb)

    # Create an xarray.Dataset
    ds = xr.Dataset(
        {
            "t_cam": (["level", "latitude", "longitude"], t_cam.astype(np.float32)),
            "u_cam": (["level", "latitude", "longitude"], u_cam.astype(np.float32)),
            "v_cam": (["level", "latitude", "longitude"], v_cam.astype(np.float32)),
            "q_cam": (["level", "latitude", "longitude"], q_cam.astype(np.float32)),
            "cldliq_cam": (["level", "latitude", "longitude"], cldliq_cam.astype(np.float32)),
            "cldice_cam": (["level", "latitude", "longitude"], cldice_cam.astype(np.float32))
        },
        coords={
            "latitude": (["latitude"], grblat),
            "longitude": (["longitude"], grblon)
        }
    )
    output_filename = "py_era5_on_hybrid.nc"
    ds.to_netcdf(output_filename)

    print("=================================================================")
    print("************ AFTER VERTICAL INTERP")
    print(f"Max T: {np.max(t_cam)}   min T: {np.min(t_cam)}")
    print(f"Max U: {np.max(u_cam)}   min U: {np.min(u_cam)}")
    print(f"Max V: {np.max(v_cam)}   min V: {np.min(v_cam)}")
    print(f"Max Q: {np.max(q_cam)}   min Q: {np.min(q_cam)}")
    print(f"Max PS: {np.max(ps)}   min PS: {np.min(ps)}")
    print(f"Max CLDICE: {np.max(cldice_cam)}   min CLDICE: {np.min(cldice_cam)}")
    print(f"Max CLDLIQ: {np.max(cldliq_cam)}   min CLDLIQ: {np.min(cldliq_cam)}")
    print("=================================================================")

    print("==CLEAN after vert interp")
    del u_gfs, v_gfs, t_gfs, q_gfs, cldice_gfs, cldliq_gfs

    ps_fv, selat, selon = horizremap.remap_with_weights_wrapper(ps, wgt_filename)
    t_fv, selat, selon = horizremap.remap_with_weights_wrapper(t_cam, wgt_filename)
    u_fv, selat, selon = horizremap.remap_with_weights_wrapper(u_cam, wgt_filename)
    v_fv, selat, selon = horizremap.remap_with_weights_wrapper(v_cam, wgt_filename)
    q_fv, selat, selon = horizremap.remap_with_weights_wrapper(q_cam, wgt_filename)
    cldice_fv, selat, selon = horizremap.remap_with_weights_wrapper(cldice_cam, wgt_filename)
    cldliq_fv, selat, selon = horizremap.remap_with_weights_wrapper(cldliq_cam, wgt_filename)

    # Create an xarray.Dataset
    ds = xr.Dataset(
        {
            "lat": (["ncol"], selat.astype(np.float32)),
            "lon": (["ncol"], selon.astype(np.float32)),
            "ps_fv": (["ncol"], ps_fv.astype(np.float32)),
            "t_fv": (["level", "ncol"], t_fv.astype(np.float32)),
            "u_fv": (["level", "ncol"], u_fv.astype(np.float32)),
            "v_fv": (["level", "ncol"], v_fv.astype(np.float32)),
            "q_fv": (["level", "ncol"], q_fv.astype(np.float32)),
            "cldliq_fv": (["level", "ncol"], cldliq_fv.astype(np.float32)),
            "cldice_fv": (["level", "ncol"], cldice_fv.astype(np.float32))
        },
    )
    output_filename = "py_era5_regrid.nc"
    ds.to_netcdf(output_filename)

if __name__ == "__main__":
    main()

