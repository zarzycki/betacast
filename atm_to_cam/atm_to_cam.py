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
from scipy.interpolate import RegularGridInterpolator

from pyfuncs import *
import vertremap
import horizremap
import topoadjust
from constants import p0, NC_FLOAT_FILL, dtime_map, QMINTHRESH, QMAXTHRESH, CLDMINTHRESH

def main():

    ps_wet_to_dry=False
    output_diag=True

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
    adjust_config = args.adjust_config
    add_cloud_vars = args.add_cloud_vars
    compress_file = args.compress_file
    se_inic = args.se_inic

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

    # Cleanup
    del data_vars

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

    t_cam = vertremap.pressure_to_hybrid(grblev, t_gfs, ps, hya, hyb)
    u_cam = vertremap.pressure_to_hybrid(grblev, u_gfs, ps, hya, hyb)
    v_cam = vertremap.pressure_to_hybrid(grblev, v_gfs, ps, hya, hyb)
    q_cam = vertremap.pressure_to_hybrid(grblev, q_gfs, ps, hya, hyb)
    cldice_cam = vertremap.pressure_to_hybrid(grblev, cldice_gfs, ps, hya, hyb)
    cldliq_cam = vertremap.pressure_to_hybrid(grblev, cldliq_gfs, ps, hya, hyb)

    # Use the print_debug_file function to create and save the xarray.Dataset
    print_debug_file(
          "py_era5_on_hybrid.nc",
          t_cam=(["level", "latitude", "longitude"], t_cam),
          u_cam=(["level", "latitude", "longitude"], u_cam),
          v_cam=(["level", "latitude", "longitude"], v_cam),
          q_cam=(["level", "latitude", "longitude"], q_cam),
          cldliq_cam=(["level", "latitude", "longitude"], cldliq_cam),
          cldice_cam=(["level", "latitude", "longitude"], cldice_cam),
          latitude=(["latitude"], grblat),
          longitude=(["longitude"], grblon)
    )

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
    t_fv, _, _ = horizremap.remap_with_weights_wrapper(t_cam, wgt_filename)
    u_fv, _, _ = horizremap.remap_with_weights_wrapper(u_cam, wgt_filename)
    v_fv, _, _ = horizremap.remap_with_weights_wrapper(v_cam, wgt_filename)
    q_fv, _, _ = horizremap.remap_with_weights_wrapper(q_cam, wgt_filename)
    cldice_fv, _, _ = horizremap.remap_with_weights_wrapper(cldice_cam, wgt_filename)
    cldliq_fv, _, _ = horizremap.remap_with_weights_wrapper(cldliq_cam, wgt_filename)

    if dycore == "se":
        print_debug_file("py_era5_regrid.nc",
                     lat=(["ncol"], selat),
                     lon=(["ncol"], selon),
                     ps_fv=(["ncol"], ps_fv),
                     t_fv=(["level", "ncol"], t_fv),
                     u_fv=(["level", "ncol"], u_fv),
                     v_fv=(["level", "ncol"], v_fv),
                     q_fv=(["level", "ncol"], q_fv),
                     cldliq_fv=(["level", "ncol"], cldliq_fv),
                     cldice_fv=(["level", "ncol"], cldice_fv))
    elif dycore == "fv":
        print_debug_file("py_era5_regrid.nc",
                     lat=(["lat"], selat),
                     lon=(["lon"], selon),
                     ps_fv=(["lat", "lon"], ps_fv),
                     t_fv=(["level", "lat", "lon"], t_fv),
                     u_fv=(["level", "lat", "lon"], u_fv),
                     v_fv=(["level", "lat", "lon"], v_fv),
                     q_fv=(["level", "lat", "lon"], q_fv),
                     cldliq_fv=(["level", "lat", "lon"], cldliq_fv),
                     cldice_fv=(["level", "lat", "lon"], cldice_fv))

    print("=" * 65)
    print("************ AFTER HORIZONTAL INTERP")
    print(f"Max T: {np.max(t_fv):.6f}   min T: {np.nanmin(t_fv):.6f}")
    print(f"Max U: {np.max(u_fv):.6f}   min U: {np.nanmin(u_fv):.6f}")
    print(f"Max V: {np.max(v_fv):.6f}   min V: {np.nanmin(v_fv):.6f}")
    print(f"Max Q: {np.max(q_fv):.6f}   min Q: {np.nanmin(q_fv):.6f}")
    print(f"Max PS: {np.max(ps):.6f}   min PS: {np.nanmin(ps):.6f}")
    print(f"Max CLDICE: {np.max(cldice_fv):.6f}   min CLDICE: {np.nanmin(cldice_fv):.6f}")
    print(f"Max CLDLIQ: {np.max(cldliq_fv):.6f}   min CLDLIQ: {np.nanmin(cldliq_fv):.6f}")
    print("=" * 65)

    print("==CLEAN after horizontal interp")
    del t_cam, u_cam, v_cam, q_cam, cldice_cam, cldliq_cam, ps

    grid_dims = ps_fv.shape
    if dycore == "se":
        ncol = grid_dims
    elif dycore == "fv":
        nfvlat, nfvlon = grid_dims

    if dycore == "fv":
        print("TOPOADJUST_FV: unpacking fv vars")
        ps_fv = latlon_to_ncol(ps_fv)
        t_fv = latlon_to_ncol(t_fv)
        q_fv = latlon_to_ncol(q_fv)
        u_fv = latlon_to_ncol(u_fv)
        v_fv = latlon_to_ncol(v_fv)
        cldice_fv = latlon_to_ncol(cldice_fv)
        cldliq_fv = latlon_to_ncol(cldliq_fv)

    correct_or_not = topoadjust.topo_adjustment(ps_fv, t_fv, q_fv, u_fv, v_fv, cldliq_fv, cldice_fv, hya, hyb, dycore, model_topo_file, datasource, grb_file, lev, yearstr, monthstr, daystr, cyclestr, wgt_filename, adjust_config, RDADIR, add_cloud_vars)

    if ps_wet_to_dry:
        if output_diag:
            ps_fv_before = np.copy(ps_fv)
        ps_fv, pw_fv = ps_wet_to_dry_conversion(ps_fv, q_fv, hyai, hybi, p0, verbose=True)
#         print_debug_file(
#                         "py_era5_tpw.nc",
#                         lat=(["ncol"], selat),
#                         lon=(["ncol"], selon),
#                         pw_fv=(["ncol"], pw_fv),
#                         ps_fv_after=(["ncol"], ps_fv),
#                         ps_fv=(["ncol"], ps_fv_before)
#                         )

    # Repack FV
    if dycore == "fv":
        ps_fv, t_fv, q_fv, u_fv, v_fv, cldice_fv, cldliq_fv, correct_or_not = repack_fv(ps_fv, t_fv, q_fv, u_fv, v_fv, cldice_fv, cldliq_fv, grid_dims, correct_or_not=correct_or_not)

    if dycore == "se":
        print_debug_file("py_era5_topoadjust.nc",
                         lat=(["ncol"], selat),
                         lon=(["ncol"], selon),
                         ps_fv=(["ncol"], ps_fv),
                         correct_or_not=(["ncol"], correct_or_not),
                         t_fv=(["level", "ncol"], t_fv),
                         u_fv=(["level", "ncol"], u_fv),
                         v_fv=(["level", "ncol"], v_fv),
                         q_fv=(["level", "ncol"], q_fv),
                         cldliq_fv=(["level", "ncol"], cldliq_fv),
                         cldice_fv=(["level", "ncol"], cldice_fv))
    elif dycore == "fv":
        print_debug_file("py_era5_topoadjust.nc",
                         lat=(["lat"], selat),
                         lon=(["lon"], selon),
                         ps_fv=(["lat", "lon"], ps_fv),
                         correct_or_not=(["lat", "lon"], correct_or_not),
                         t_fv=(["level", "lat", "lon"], t_fv),
                         u_fv=(["level", "lat", "lon"], u_fv),
                         v_fv=(["level", "lat", "lon"], v_fv),
                         q_fv=(["level", "lat", "lon"], q_fv),
                         cldliq_fv=(["level", "lat", "lon"], cldliq_fv),
                         cldice_fv=(["level", "lat", "lon"], cldice_fv))


    q_fv = clip_and_count(q_fv, min_thresh=QMINTHRESH, max_thresh=QMAXTHRESH, var_name="Q")
    cldliq_fv = clip_and_count(cldliq_fv, min_thresh=CLDMINTHRESH, var_name="CLDLIQ")
    cldice_fv = clip_and_count(cldice_fv, min_thresh=CLDMINTHRESH, var_name="CLDICE")

    if dycore == "se":
        ps_fv = numpy_to_dataarray(ps_fv, dims=['ncol'], attrs={'units': 'Pa', "_FillValue": np.float32(NC_FLOAT_FILL)})
        u_fv = numpy_to_dataarray(u_fv, dims=['lev', 'ncol'], attrs={'units': 'm/s', "_FillValue": np.float32(NC_FLOAT_FILL)})
        v_fv = numpy_to_dataarray(v_fv, dims=['lev', 'ncol'], attrs={'units': 'm/s', "_FillValue": np.float32(NC_FLOAT_FILL)})
        t_fv = numpy_to_dataarray(t_fv, dims=['lev', 'ncol'], attrs={'units': 'K', "_FillValue": np.float32(NC_FLOAT_FILL)})
        q_fv = numpy_to_dataarray(q_fv, dims=['lev', 'ncol'], attrs={'units': 'kg/kg', "_FillValue": np.float32(NC_FLOAT_FILL)})
        cldliq_fv = numpy_to_dataarray(cldliq_fv, dims=['lev', 'ncol'], attrs={'units': 'kg/kg', "_FillValue": np.float32(NC_FLOAT_FILL)})
        cldice_fv = numpy_to_dataarray(cldice_fv, dims=['lev', 'ncol'], attrs={'units': 'kg/kg', "_FillValue": np.float32(NC_FLOAT_FILL)})
        selat = numpy_to_dataarray(selat, dims=['ncol'], attrs={"_FillValue": -900., "long_name": "latitude", "units": "degrees_north"})
        selon = numpy_to_dataarray(selon, dims=['ncol'], attrs={"_FillValue": -900., "long_name": "longitude", "units": "degrees_east"})
        if 'correct_or_not' in locals():
            correct_or_not = numpy_to_dataarray(correct_or_not, dims=['ncol'], attrs={"_FillValue": -1.0}, name='correct_or_not')

        ps_fv = add_time_define_precision(ps_fv, write_type, True)
        u_fv = add_time_define_precision(u_fv, write_type, True)
        v_fv = add_time_define_precision(v_fv, write_type, True)
        t_fv = add_time_define_precision(t_fv, write_type, True)
        q_fv = add_time_define_precision(q_fv, write_type, True)
        cldliq_fv = add_time_define_precision(cldliq_fv, write_type, True)
        cldice_fv = add_time_define_precision(cldice_fv, write_type, True)

    elif dycore == "fv":

        print("FV: need to interpolate u/v to slat/slon")
        nfvslat = nfvlat - 1
        nfvslon = nfvlon

        fvlat = np.linspace(-90, 90, nfvlat)
        delfvlon = 360.0 / nfvlon
        fvlon = np.linspace(0, 360.0 - delfvlon, nfvlon)

        fvslat = (fvlat[:-1] + fvlat[1:]) / 2.0
        fvslon = fvlon - (delfvlon / 2.0)

        print("FV: getting weights")
        w_stag = latRegWgt(fvslat, "double", 0)

        # Interpolate u and v to slat and slon
        us_fv = np.zeros((numlevels, nfvslat, nfvlon))
        vs_fv = np.zeros((numlevels, nfvlat, nfvslon))

        slat_lon_mesh = np.array(np.meshgrid(fvslat, fvlon)).T.reshape(-1, 2)
        lat_slon_mesh = np.array(np.meshgrid(fvlat, fvslon)).T.reshape(-1, 2)

        # Create a cyclic ghost point for slon interpolation
        ghost_lon = fvlon[-1] - 360.0
        tmp_fvlon = np.hstack(([ghost_lon], fvlon))

        for level in range(numlevels):
            # Extend v to match ghost point since this var is interpolated in slon
            v_fv_extended = np.hstack((v_fv[level, :, -1][:, np.newaxis], v_fv[level, :, :]))
            u_interpolator = RegularGridInterpolator((fvlat, fvlon), u_fv[level, :, :], method='linear')
            v_interpolator = RegularGridInterpolator((fvlat, tmp_fvlon), v_fv_extended, method='linear')
            us_fv[level, :, :] = u_interpolator(slat_lon_mesh).reshape(nfvslat, nfvlon)
            vs_fv[level, :, :] = v_interpolator(lat_slon_mesh).reshape(nfvlat, nfvslon)

        # Clean up so we aren't confused later
        del ghost_lon, tmp_fvlon, v_interpolator, u_interpolator, v_fv_extended

        ps_fv = numpy_to_dataarray(ps_fv, dims=['lat', 'lon'], attrs={'units': 'Pa', "_FillValue": np.float32(NC_FLOAT_FILL)})
        u_fv = numpy_to_dataarray(u_fv, dims=['lev', 'lat', 'lon'], attrs={'units': 'm/s', "_FillValue": np.float32(NC_FLOAT_FILL)})
        v_fv = numpy_to_dataarray(v_fv, dims=['lev', 'lat', 'lon'], attrs={'units': 'm/s', "_FillValue": np.float32(NC_FLOAT_FILL)})
        us_fv = numpy_to_dataarray(us_fv, dims=['lev', 'slat', 'lon'], attrs={'units': 'm/s', "_FillValue": np.float32(NC_FLOAT_FILL)})
        vs_fv = numpy_to_dataarray(vs_fv, dims=['lev', 'lat', 'slon'], attrs={'units': 'm/s', "_FillValue": np.float32(NC_FLOAT_FILL)})
        t_fv = numpy_to_dataarray(t_fv, dims=['lev', 'lat', 'lon'], attrs={'units': 'K', "_FillValue": np.float32(NC_FLOAT_FILL)})
        q_fv = numpy_to_dataarray(q_fv, dims=['lev', 'lat', 'lon'], attrs={'units': 'kg/kg', "_FillValue": np.float32(NC_FLOAT_FILL)})
        cldliq_fv = numpy_to_dataarray(cldliq_fv, dims=['lev', 'lat', 'lon'], attrs={'units': 'kg/kg', "_FillValue": np.float32(NC_FLOAT_FILL)})
        cldice_fv = numpy_to_dataarray(cldice_fv, dims=['lev', 'lat', 'lon'], attrs={'units': 'kg/kg', "_FillValue": np.float32(NC_FLOAT_FILL)})
        selat = numpy_to_dataarray(fvlat, dims=['lat'], attrs={"_FillValue": -900., "long_name": "latitude", "units": "degrees_north"})
        selon = numpy_to_dataarray(fvlon, dims=['lon'], attrs={"_FillValue": -900., "long_name": "longitude", "units": "degrees_east"})
        seslat = numpy_to_dataarray(fvslat, dims=['slat'], attrs={"_FillValue": -900., "long_name": "latitude", "units": "degrees_north"})
        seslon = numpy_to_dataarray(fvslon, dims=['slon'], attrs={"_FillValue": -900., "long_name": "longitude", "units": "degrees_east"})
        if 'correct_or_not' in locals():
            correct_or_not = numpy_to_dataarray(correct_or_not, dims=['lat', 'lon'], attrs={"_FillValue": -1.0}, name='correct_or_not')

        ps_fv = add_time_define_precision(ps_fv, write_type, False)
        u_fv = add_time_define_precision(u_fv, write_type, False)
        v_fv = add_time_define_precision(v_fv, write_type, False)
        us_fv = add_time_define_precision(us_fv, write_type, False, lat_dim="slat")
        vs_fv = add_time_define_precision(vs_fv, write_type, False, lon_dim="slon")
        t_fv = add_time_define_precision(t_fv, write_type, False)
        q_fv = add_time_define_precision(q_fv, write_type, False)
        cldliq_fv = add_time_define_precision(cldliq_fv, write_type, False)
        cldice_fv = add_time_define_precision(cldice_fv, write_type, False)

    # Create CF-compliant time
    time, time_atts = create_cf_time(int(yearstr), int(monthstr), int(daystr), int(cyclestr))
    time = numpy_to_dataarray(time, dims=['time'], attrs=time_atts)

    # Reload from the template xarray for metadata purposes
    hya, hyb, hyai, hybi, lev, ilev = load_cam_levels(PATHTOHERE, numlevels, load_xarray = True)

    # Data set to be written out
    ds = xr.Dataset(
        {
            "PS": ps_fv,
            "U": u_fv,
            "V": v_fv,
            "T": t_fv,
            "Q": q_fv,
            "CLDLIQ": cldliq_fv,
            "CLDICE": cldice_fv,
            "hyam": hya,
            "hybm": hyb,
            "hyai": hyai,
            "hybi": hybi,
            "lev": lev,
            "ilev": ilev,
            "time": time
        }
    )

    ds["lat"] = selat
    ds["lon"] = selon
    if dycore == "fv":
        ds["US"] = us_fv
        ds["VS"] = vs_fv
        ds["slat"] = seslat
        ds["slon"] = seslon

    if 'correct_or_not' in locals():
        ds["correct_or_not"] = correct_or_not

    # Add global attributes
    ds.attrs.update({
        "title": "Betacast-generated ncdata file",
        "source_file": data_filename,
        "wgt_file": wgt_filename,
        "init_date": YYYYMMDDHH,
        "creation_date": np.datetime_as_string(np.datetime64('now')),
        "dycore": dycore,
        "datasource": datasource
    })

    # Turn off fill value in relevant coordinate variables
    encoding = {
        var: {'_FillValue': None} for var in ["hyam", "hybm", "hyai", "hybi", "lev", "ilev", "time"]
    }

    # Determine file format and compression based on variable size
    if compress_file:
        for var in list(ds.data_vars) + list(ds.coords):
            if var not in encoding:
                if var in encoding:
                    encoding[var].update({"zlib": True, "complevel": 1})
                else:
                    encoding[var] = {"zlib": True, "complevel": 1}
        ds.to_netcdf(se_inic, format="NETCDF4_CLASSIC", unlimited_dims=["time"], encoding=encoding)
    else:
        netcdf_format = "NETCDF4" if var_max_size >= 4e9 else "NETCDF3_64BIT"
        ds.to_netcdf(se_inic, format=netcdf_format, unlimited_dims=["time"], encoding=encoding)

if __name__ == "__main__":
    main()

