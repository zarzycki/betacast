# Core modules
import os
import sys
import logging

# Third party modules
import numpy as np
import xarray as xr

# Betacast modules
module_path = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', 'py_functions'))
if module_path not in sys.path:
    sys.path.append(module_path)
import pyfuncs
import vertremap
import horizremap
import topoadjust
import packing
import loaddata
import meteo
from constants import (
    p0, NC_FLOAT_FILL, dtime_map, QMINTHRESH, QMAXTHRESH, CLDMINTHRESH,
    ps_wet_to_dry, output_diag, w_smooth_iter, damp_upper_winds_mpas, grav
)

args = pyfuncs.parse_args()

pyfuncs.configure_logging(args.verbose)


def main():

    logger = logging.getLogger(__name__)
    logger.info("Main function started")

    # Is the Betacast path available to us?
    local_only, PATHTOHERE = pyfuncs.get_betacast_path("atm_to_cam")

    args_dict = {
        "Datasource": args.datasource,
        "Number of levels": args.numlevels,
        "Initialization date": args.YYYYMMDDHH,
        "Data filename": args.data_filename,
        "Weight filename": args.wgt_filename,
        "SE Inic output file": args.se_inic,
        "Dycore": args.dycore,
        "RDADIR": args.RDADIR,
        "MPAS as CAM": args.mpas_as_cam,
        "Compress file": args.compress_file,
        "Write floats": args.write_floats,
        "Add cloud vars": args.add_cloud_vars,
        "Adjust config": args.adjust_config,
        "Model topo file": args.model_topo_file,
        "MOD in topo": args.mod_in_topo,
        "MOD remap file": args.mod_remap_file,
        "write_debug_files": args.write_debug_files,
        "write_debug_dir": args.write_debug_dir,
    }

    for key, value in args_dict.items():
        logging.info(f"{key}: {value}")

    # Assume these are parsed from command-line arguments
    dycore = args.dycore
    numlevels = args.numlevels
    model_topo_file = args.model_topo_file
    YYYYMMDDHH = args.YYYYMMDDHH
    data_filename = args.data_filename
    wgt_filename = args.wgt_filename
    mpasfile = args.mpasfile
    datasource = args.datasource
    mod_in_topo = args.mod_in_topo
    mod_remap_file = args.mod_remap_file
    write_floats = args.write_floats
    RDADIR = args.RDADIR
    adjust_config = args.adjust_config
    add_cloud_vars = args.add_cloud_vars
    compress_file = args.compress_file
    se_inic = args.se_inic
    mpas_as_cam = args.mpas_as_cam
    write_debug_files = args.write_debug_files
    DEBUGDIR = args.write_debug_dir

    if write_debug_files:
        pyfuncs.create_folder(DEBUGDIR)
        pyfuncs.delete_files(DEBUGDIR, "py_*.nc")

    # Check if we are using MPAS and the model topo file is provided and exists
    if dycore == "mpas" and model_topo_file:
        if os.path.isfile(model_topo_file):
            logging.info(f"Setting mpasfile to {model_topo_file}")
            mpasfile = model_topo_file
        else:
            logging.error(f"MPAS set and model topo file {model_topo_file} not found, exiting...")
            sys.exit(1)

    # Check if mpasfile is set if dycore is "mpas"
    if dycore == "mpas" and mpasfile is None:
        logging.error("No mpasfile passed in but mpas dycore selected.")
        logging.error("MPAS init file is needed to mimic + get zcoord from. Exiting...")
        sys.exit(1)

    # Error check if CAM data is used as input
    if datasource == "CAM":
        if not mod_in_topo:
            logging.error("Using CAM data, but no mod_in_topo passed in, exiting...")
            sys.exit(1)
        if not mod_remap_file:
            logging.error("Using CAM data, but no mod_remap_file passed in, exiting...")
            sys.exit(1)

    # Toggle whether the output streams will be floats or doubles
    write_type = "float" if write_floats else "double"
    logging.info(f"Output type set to: {write_type}")

    # ===== Getting date from YYYYMMDDHH
    yearstr, monthstr, daystr, cyclestr = pyfuncs.split_by_lengths(str(YYYYMMDDHH), dtime_map)

    logging.info(f"Regridding analysis from: {yearstr} {monthstr} {daystr} {cyclestr}Z")

    logging.info("---------------------------------------------------------")
    logging.info(f"Using this file: {data_filename}")
    logging.info(f"Using this remap: {wgt_filename}")

    hya, hyb, hyai, hybi, lev, ilev = loaddata.load_cam_levels(PATHTOHERE, numlevels)

    # IMPORTANT! data_vars should be organized top-to-bottom when loaddata returns
    # (i.e., lowest pressure/highest z at 0 index of lev)
    # I attempt to account for this elsewhere in the code with flips, but make no promises
    if datasource == 'GFS' or datasource == 'HWRF':
        data_vars = loaddata.load_CFSR_data(data_filename, dycore)
    elif datasource == 'ERA5RDA':
        data_vars = loaddata.load_ERA5RDA_data(RDADIR, data_filename, yearstr, monthstr, daystr, cyclestr, dycore)
    elif datasource == 'CAM':
        data_vars = loaddata.load_cam_data(data_filename, YYYYMMDDHH, mod_in_topo, mod_remap_file, dycore, write_debug_files=write_debug_files,write_debug_dir=DEBUGDIR)

    logging.info("Input Data Level information")
    logging.info(f"Number: {len(data_vars['lev'])}")
    logging.info(f"Max: {np.max(data_vars['lev'])}")
    logging.info(f"Min: {np.min(data_vars['lev'])}")
    logging.info(f"Input Coords: nlat: {len(data_vars['lat'])}, nlon: {len(data_vars['lon'])}")

    if dycore == "mpas":
        # Initialize variables
        data_vars['theta'] = data_vars['t'].copy()
        data_vars['rho'] = data_vars['t'].copy()

        # Create a 3-D pressure field from constant pressure surfaces
        data_vars['pres'] = np.broadcast_to(data_vars['lev'][:, np.newaxis, np.newaxis], data_vars['t'].shape)

        # If w is omega (Pa/s), convert to vertical velocity in m/s
        if data_vars['w_is_omega']:
            data_vars['w'] = meteo.omega_to_w(data_vars['w'], data_vars['pres'], data_vars['t'])

        # Smooth w if needed
        if w_smooth_iter > 0:
            #data_vars['w'] = pyfuncs.smooth_with_gaussian(data_vars['w'], w_smooth_iter)
            data_vars['w'] = pyfuncs.smooth_with_smth9(data_vars['w'], w_smooth_iter)

        # If z is reported in geopotential, convert to geometric height
        if data_vars['z_is_phi']:
            data_vars['z'] = data_vars['z'] / grav

        # Calculate potential temperature using full pressure and actual T
        data_vars['theta'] = meteo.pot_temp(data_vars['pres'], data_vars['t'])

        # Calculate density from pressure, moisture, and temperature
        data_vars['rho'] = meteo.calculate_rho_gfs(data_vars['pres'], data_vars['q'], data_vars['t'])

    # Print diagnostics
    pyfuncs.print_min_max_dict(data_vars)

    if dycore == 'fv' or dycore == 'se':

        if write_debug_files:
            pyfuncs.print_debug_file(
                DEBUGDIR+"/"+"py_era5_before_interp.nc",
                ps_cam=(["lat", "lon"], data_vars['ps']),
                t_cam=(["lev_p", "lat", "lon"], data_vars['t']),
                u_cam=(["lev_p", "lat", "lon"], data_vars['u']),
                v_cam=(["lev_p", "lat", "lon"], data_vars['v']),
                q_cam=(["lev_p", "lat", "lon"], data_vars['q']),
                cldliq_cam=(["lev_p", "lat", "lon"], data_vars['cldliq']),
                cldice_cam=(["lev_p", "lat", "lon"], data_vars['cldice']),
                lat=(["lat"], data_vars['lat']),
                lon=(["lon"], data_vars['lon'])
            )

        data_vint = vertremap.pres2hyb_all(data_vars, data_vars['ps'], hya, hyb)

        if write_debug_files:
            pyfuncs.print_debug_file(
                DEBUGDIR+"/"+"py_era5_on_hybrid.nc",
                ps_cam=(["latitude", "longitude"], data_vint['ps']),
                t_cam=(["level", "latitude", "longitude"], data_vint['t']),
                u_cam=(["level", "latitude", "longitude"], data_vint['u']),
                v_cam=(["level", "latitude", "longitude"], data_vint['v']),
                q_cam=(["level", "latitude", "longitude"], data_vint['q']),
                cldliq_cam=(["level", "latitude", "longitude"], data_vint['cldliq']),
                cldice_cam=(["level", "latitude", "longitude"], data_vint['cldice']),
                latitude=(["latitude"], data_vint['lat']),
                longitude=(["longitude"], data_vint['lon'])
            )

        pyfuncs.print_min_max_dict(data_vint)

    if dycore == 'mpas':
        # Load the MPAS init file and get the zgrid
        logging.info(f"getting grid from {mpasfile}")
        mpas_file = xr.open_dataset(mpasfile)
        mpas_data = {}
        mpas_data['z'] = mpas_file['zgrid'].values
        mpas_data['ncell'], mpas_data['nlevi'] = mpas_data['z'].shape
        mpas_data['nlev'] = mpas_data['nlevi'] - 1

        logging.info(f"MPAS GRID: nlev: {mpas_data['nlev']}    nlevi: {mpas_data['nlevi']}     ncell: {mpas_data['ncell']}")
        pyfuncs.print_min_max_dict(mpas_data)

        # Set up some dummy arrays since we aren't going to vert interp first
        data_vint = data_vars.copy()

        pyfuncs.print_min_max_dict(data_vint)

    data_horiz = horizremap.remap_all(data_vint, wgt_filename, dycore=dycore)

    if write_debug_files:
        if dycore == "se":
            pyfuncs.print_debug_file(DEBUGDIR+"/"+"py_era5_regrid.nc",
                         lat=(["ncol"], data_horiz['lat']),
                         lon=(["ncol"], data_horiz['lon']),
                         ps_fv=(["ncol"], data_horiz['ps']),
                         t_fv=(["level", "ncol"], data_horiz['t']),
                         u_fv=(["level", "ncol"], data_horiz['u']),
                         v_fv=(["level", "ncol"], data_horiz['v']),
                         q_fv=(["level", "ncol"], data_horiz['q']),
                         cldliq_fv=(["level", "ncol"], data_horiz['cldliq']),
                         cldice_fv=(["level", "ncol"], data_horiz['cldice']))
        elif dycore == "fv":
            pyfuncs.print_debug_file(DEBUGDIR+"/"+"py_era5_regrid.nc",
                         lat=(["lat"], data_horiz['lat']),
                         lon=(["lon"], data_horiz['lon']),
                         ps_fv=(["lat", "lon"], data_horiz['ps']),
                         t_fv=(["level", "lat", "lon"], data_horiz['t']),
                         u_fv=(["level", "lat", "lon"], data_horiz['u']),
                         v_fv=(["level", "lat", "lon"], data_horiz['v']),
                         q_fv=(["level", "lat", "lon"], data_horiz['q']),
                         cldliq_fv=(["level", "lat", "lon"], data_horiz['cldliq']),
                         cldice_fv=(["level", "lat", "lon"], data_horiz['cldice']))
        elif dycore == "mpas":
            pyfuncs.print_debug_file(DEBUGDIR+"/"+"py_era5_regrid.nc",
                         lat=(["ncol"], data_horiz['lat']),
                         lon=(["ncol"], data_horiz['lon']),
                         ps_fv=(["ncol"], data_horiz['ps']),
                         t_fv=(["level", "ncol"], data_horiz['t']),
                         u_fv=(["level", "ncol"], data_horiz['u']),
                         v_fv=(["level", "ncol"], data_horiz['v']),
                         q_fv=(["level", "ncol"], data_horiz['q']),
                         cldliq_fv=(["level", "ncol"], data_horiz['cldliq']),
                         cldice_fv=(["level", "ncol"], data_horiz['cldice']),
                         z_fv=(["level", "ncol"], data_horiz['z']),
                         theta_fv=(["level", "ncol"], data_horiz['theta']),
                         rho_fv=(["level", "ncol"], data_horiz['rho']),
                         w_fv=(["level", "ncol"], data_horiz['w']))

    pyfuncs.print_min_max_dict(data_horiz)

    if dycore == 'mpas':
        logging.info("Performing vertical interpolation at each MPAS column...")
        # Call the processing function with your data
        data_horiz = vertremap.interpolate_mpas_columns_wrapper(
            mpas_data, data_horiz
        )
        logging.info("... done performing vertical interpolation at each MPAS column!")

        # Set it so flow doesn't go through the lower boundary condition
        data_horiz['w'] = meteo.noflux_boundary_condition(data_horiz['w'], mpas_data['nlev'])

        # Damp the top few layers of MPAS winds
        if damp_upper_winds_mpas:
            data_horiz['u'], data_horiz['v'] = meteo.mpas.damp_upper_level_winds(data_horiz['u'], data_horiz['v'], mpas_data['nlev'])

        if not mpas_as_cam:
            # put u + v on cell edges...
            logging.info("Projecting u + v to velocity normal to edge...")
            data_horiz['uNorm'] = horizremap.uv_cell_to_edge(data_horiz['u'], data_horiz['v'], mpas_data['nlev'], mpas_file['lonEdge'].values, mpas_file['latEdge'].values,
                                        mpas_file['lonCell'].values, mpas_file['latCell'].values, mpas_file['edgeNormalVectors'].values, mpas_file['cellsOnEdge'].values)

            logging.info("... done projecting u + v to velocity normal to edge!")

            # Clip relevant variables
            data_horiz['q'] = pyfuncs.clip_and_count(data_horiz['q'], min_thresh=QMINTHRESH, max_thresh=QMAXTHRESH, var_name="Q")

            # If not MPAS as CAM, we can just end here.
            logging.info("Writing MPAS file...")
            mpas_file['u'].values[0, :, :] = data_horiz['uNorm'].T
            mpas_file['qv'].values[0, :, :] = data_horiz['q'].T
            mpas_file['rho'].values[0, :, :] = data_horiz['rho'].T
            mpas_file['theta'].values[0, :, :] = data_horiz['theta'].T
            mpas_file['w'].values[0, :, :] = data_horiz['w'].T
            mpas_file.to_netcdf(se_inic, format='NETCDF4')
            mpas_file.close()
            logging.info("Done generating MPAS initial condition file, exiting...")
            sys.exit(0)

    grid_dims = data_horiz['ps'].shape
    if dycore == "se" or dycore == "mpas":
        ncol = grid_dims
    elif dycore == "fv":
        nfvlat, nfvlon = grid_dims

    if dycore == "fv":
        data_horiz = packing.unpack_fv(data_horiz)

    if dycore == "se" or dycore == "fv":
        data_horiz['correct_or_not'] = topoadjust.topo_adjustment(data_horiz, dycore, model_topo_file, adjust_config)

    if ps_wet_to_dry:
        if output_diag:
            ps_fv_before = np.copy(data_horiz['ps'])
        data_horiz['ps'], data_horiz['pw'] = meteo.ps_wet_to_dry_conversion(data_horiz['ps'], data_horiz['q'], hyai, hybi, p0, verbose=True)
#         pyfuncs.print_debug_file(
#                         "py_era5_tpw.nc",
#                         lat=(["ncol"], selat),
#                         lon=(["ncol"], selon),
#                         pw_fv=(["ncol"], pw_fv),
#                         ps_fv_after=(["ncol"], ps_fv),
#                         ps_fv=(["ncol"], ps_fv_before)
#                         )

    # Repack FV
    if dycore == "fv":
        data_horiz = packing.repack_fv(data_horiz, grid_dims)

    if write_debug_files:
        if dycore == "se":
            pyfuncs.print_debug_file(DEBUGDIR+"/"+"py_era5_topoadjust.nc",
                            lat=(["ncol"], data_horiz['lat']),
                            lon=(["ncol"], data_horiz['lon']),
                            ps_fv=(["ncol"], data_horiz['ps']),
                            correct_or_not=(["ncol"], data_horiz['correct_or_not']),
                            t_fv=(["level", "ncol"], data_horiz['t']),
                            u_fv=(["level", "ncol"], data_horiz['u']),
                            v_fv=(["level", "ncol"], data_horiz['v']),
                            q_fv=(["level", "ncol"], data_horiz['q']),
                            cldliq_fv=(["level", "ncol"], data_horiz['cldliq']),
                            cldice_fv=(["level", "ncol"], data_horiz['cldice']))
        elif dycore == "fv":
            pyfuncs.print_debug_file(DEBUGDIR+"/"+"py_era5_topoadjust.nc",
                            lat=(["lat"], data_horiz['lat']),
                            lon=(["lon"], data_horiz['lon']),
                            ps_fv=(["lat", "lon"], data_horiz['ps']),
                            correct_or_not=(["lat", "lon"], data_horiz['correct_or_not']),
                            t_fv=(["level", "lat", "lon"], data_horiz['t']),
                            u_fv=(["level", "lat", "lon"], data_horiz['u']),
                            v_fv=(["level", "lat", "lon"], data_horiz['v']),
                            q_fv=(["level", "lat", "lon"], data_horiz['q']),
                            cldliq_fv=(["level", "lat", "lon"], data_horiz['cldliq']),
                            cldice_fv=(["level", "lat", "lon"], data_horiz['cldice']))

    if dycore == "fv":
        logging.info("FV: need to interpolate u/v to slat/slon")

        # Initialize FV grid
        data_horiz['fvlat'], data_horiz['fvlon'], data_horiz['fvslat'], data_horiz['fvslon'] = packing.initialize_fv_grid(nfvlat, nfvlon)

        # Interpolate u and v to slat/slon
        data_horiz['us'], data_horiz['vs'] = packing.interpolate_uv_to_slat_slon(
            data_horiz, numlevels, data_horiz['fvlat'], data_horiz['fvlon'],
            data_horiz['fvslat'], data_horiz['fvslon']
        )

    data_horiz['q'] = pyfuncs.clip_and_count(data_horiz['q'], min_thresh=QMINTHRESH, max_thresh=QMAXTHRESH, var_name="Q")
    data_horiz['cldliq'] = pyfuncs.clip_and_count(data_horiz['cldliq'], min_thresh=CLDMINTHRESH, var_name="CLDLIQ")
    data_horiz['cldice'] = pyfuncs.clip_and_count(data_horiz['cldice'], min_thresh=CLDMINTHRESH, var_name="CLDICE")

    if datasource == "HWRF":
        logging.info(f"{datasource} replacing nan with _FillValue {NC_FLOAT_FILL} since regional")
        data_horiz = pyfuncs.replace_nans_with_fill_value(data_horiz, ['ps', 'u', 'v', 't', 'q', 'cldliq', 'cldice'], NC_FLOAT_FILL)

    logging.info(f"Begin writing output file: {se_inic}")
    out_data = {}

    if dycore == "se":
        out_data['ps'] = pyfuncs.numpy_to_dataarray(data_horiz['ps'], dims=['ncol'], attrs={'units': 'Pa', "_FillValue": NC_FLOAT_FILL})
        out_data['u'] = pyfuncs.numpy_to_dataarray(data_horiz['u'], dims=['lev', 'ncol'], attrs={'units': 'm/s', "_FillValue": NC_FLOAT_FILL})
        out_data['v'] = pyfuncs.numpy_to_dataarray(data_horiz['v'], dims=['lev', 'ncol'], attrs={'units': 'm/s', "_FillValue": NC_FLOAT_FILL})
        out_data['t'] = pyfuncs.numpy_to_dataarray(data_horiz['t'], dims=['lev', 'ncol'], attrs={'units': 'K', "_FillValue": NC_FLOAT_FILL})
        out_data['q'] = pyfuncs.numpy_to_dataarray(data_horiz['q'], dims=['lev', 'ncol'], attrs={'units': 'kg/kg', "_FillValue": NC_FLOAT_FILL})
        out_data['cldliq'] = pyfuncs.numpy_to_dataarray(data_horiz['cldliq'], dims=['lev', 'ncol'], attrs={'units': 'kg/kg', "_FillValue": NC_FLOAT_FILL})
        out_data['cldice'] = pyfuncs.numpy_to_dataarray(data_horiz['cldice'], dims=['lev', 'ncol'], attrs={'units': 'kg/kg', "_FillValue": NC_FLOAT_FILL})
        out_data['lat'] = pyfuncs.numpy_to_dataarray(data_horiz['lat'], dims=['ncol'], attrs={"_FillValue": -900., "long_name": "latitude", "units": "degrees_north"})
        out_data['lon'] = pyfuncs.numpy_to_dataarray(data_horiz['lon'], dims=['ncol'], attrs={"_FillValue": -900., "long_name": "longitude", "units": "degrees_east"})
        out_data['correct_or_not'] = pyfuncs.numpy_to_dataarray(data_horiz['correct_or_not'], dims=['ncol'], attrs={"_FillValue": -1.0}, name='correct_or_not')

        out_data['ps'] = pyfuncs.add_time_define_precision(out_data['ps'], write_type, True)
        out_data['u'] = pyfuncs.add_time_define_precision(out_data['u'], write_type, True)
        out_data['v'] = pyfuncs.add_time_define_precision(out_data['v'], write_type, True)
        out_data['t'] = pyfuncs.add_time_define_precision(out_data['t'], write_type, True)
        out_data['q'] = pyfuncs.add_time_define_precision(out_data['q'], write_type, True)
        out_data['cldliq'] = pyfuncs.add_time_define_precision(out_data['cldliq'], write_type, True)
        out_data['cldice'] = pyfuncs.add_time_define_precision(out_data['cldice'], write_type, True)

    elif dycore == "fv":
        out_data['ps'] = pyfuncs.numpy_to_dataarray(data_horiz['ps'], dims=['lat', 'lon'], attrs={'units': 'Pa', "_FillValue": NC_FLOAT_FILL})
        out_data['u'] = pyfuncs.numpy_to_dataarray(data_horiz['u'], dims=['lev', 'lat', 'lon'], attrs={'units': 'm/s', "_FillValue": NC_FLOAT_FILL})
        out_data['v'] = pyfuncs.numpy_to_dataarray(data_horiz['v'], dims=['lev', 'lat', 'lon'], attrs={'units': 'm/s', "_FillValue": NC_FLOAT_FILL})
        out_data['us'] = pyfuncs.numpy_to_dataarray(data_horiz['us'], dims=['lev', 'slat', 'lon'], attrs={'units': 'm/s', "_FillValue": NC_FLOAT_FILL})
        out_data['vs'] = pyfuncs.numpy_to_dataarray(data_horiz['vs'], dims=['lev', 'lat', 'slon'], attrs={'units': 'm/s', "_FillValue": NC_FLOAT_FILL})
        out_data['t'] = pyfuncs.numpy_to_dataarray(data_horiz['t'], dims=['lev', 'lat', 'lon'], attrs={'units': 'K', "_FillValue": NC_FLOAT_FILL})
        out_data['q'] = pyfuncs.numpy_to_dataarray(data_horiz['q'], dims=['lev', 'lat', 'lon'], attrs={'units': 'kg/kg', "å": NC_FLOAT_FILL})
        out_data['cldice'] = pyfuncs.numpy_to_dataarray(data_horiz['cldice'], dims=['lev', 'lat', 'lon'], attrs={'units': 'kg/kg', "_FillValue": NC_FLOAT_FILL})
        out_data['cldliq'] = pyfuncs.numpy_to_dataarray(data_horiz['cldliq'], dims=['lev', 'lat', 'lon'], attrs={'units': 'kg/kg', "_FillValue": NC_FLOAT_FILL})
        out_data['lat'] = pyfuncs.numpy_to_dataarray(data_horiz['lat'], dims=['lat'], attrs={"_FillValue": -900., "long_name": "latitude", "units": "degrees_north"})
        out_data['lon'] = pyfuncs.numpy_to_dataarray(data_horiz['lon'], dims=['lon'], attrs={"_FillValue": -900., "long_name": "longitude", "units": "degrees_east"})
        out_data['slat'] = pyfuncs.numpy_to_dataarray(data_horiz['fvslat'], dims=['slat'], attrs={"_FillValue": -900., "long_name": "latitude", "units": "degrees_north"})
        out_data['slon'] = pyfuncs.numpy_to_dataarray(data_horiz['fvslon'], dims=['slon'], attrs={"_FillValue": -900., "long_name": "longitude", "units": "degrees_east"})
        out_data['correct_or_not'] = pyfuncs.numpy_to_dataarray(data_horiz['correct_or_not'], dims=['lat', 'lon'], attrs={"_FillValue": -1.0}, name='correct_or_not')

        out_data['ps'] = pyfuncs.add_time_define_precision(out_data['ps'], write_type, False)
        out_data['u'] = pyfuncs.add_time_define_precision(out_data['u'], write_type, False)
        out_data['v'] = pyfuncs.add_time_define_precision(out_data['v'], write_type, False)
        out_data['us'] = pyfuncs.add_time_define_precision(out_data['us'], write_type, False, lat_dim="slat")
        out_data['vs'] = pyfuncs.add_time_define_precision(out_data['vs'], write_type, False, lon_dim="slon")
        out_data['t'] = pyfuncs.add_time_define_precision(out_data['t'], write_type, False)
        out_data['q'] = pyfuncs.add_time_define_precision(out_data['q'], write_type, False)
        out_data['cldice'] = pyfuncs.add_time_define_precision(out_data['cldice'], write_type, False)
        out_data['cldliq'] = pyfuncs.add_time_define_precision(out_data['cldliq'], write_type, False)

    elif dycore == "mpas":
        # Need to flip the 3-D arrays!
        out_data['ps'] = pyfuncs.numpy_to_dataarray(data_horiz['ps'], dims=['ncol'], attrs={'units': 'Pa', "_FillValue": NC_FLOAT_FILL})
        out_data['u'] = pyfuncs.numpy_to_dataarray(data_horiz['u'][::-1, :], dims=['lev', 'ncol'], attrs={'units': 'm/s', "_FillValue": NC_FLOAT_FILL})
        out_data['v'] = pyfuncs.numpy_to_dataarray(data_horiz['v'][::-1, :], dims=['lev', 'ncol'], attrs={'units': 'm/s', "_FillValue": NC_FLOAT_FILL})
        out_data['t'] = pyfuncs.numpy_to_dataarray(data_horiz['t'][::-1, :], dims=['lev', 'ncol'], attrs={'units': 'K', "_FillValue": NC_FLOAT_FILL})
        out_data['q'] = pyfuncs.numpy_to_dataarray(data_horiz['q'][::-1, :], dims=['lev', 'ncol'], attrs={'units': 'kg/kg', "_FillValue": NC_FLOAT_FILL})
        out_data['lat'] = pyfuncs.numpy_to_dataarray(data_horiz['lat'], dims=['ncol'], attrs={"_FillValue": -900., "long_name": "latitude", "units": "degrees_north"})
        out_data['lon'] = pyfuncs.numpy_to_dataarray(data_horiz['lon'], dims=['ncol'], attrs={"_FillValue": -900., "long_name": "longitude", "units": "degrees_east"})

        out_data['ps'] = pyfuncs.add_time_define_precision(out_data['ps'], write_type, True)
        out_data['u'] = pyfuncs.add_time_define_precision(out_data['u'], write_type, True)
        out_data['v'] = pyfuncs.add_time_define_precision(out_data['v'], write_type, True)
        out_data['t'] = pyfuncs.add_time_define_precision(out_data['t'], write_type, True)
        out_data['q'] = pyfuncs.add_time_define_precision(out_data['q'], write_type, True)

    # Create CF-compliant time
    time, time_atts = pyfuncs.create_cf_time(int(yearstr), int(monthstr), int(daystr), int(cyclestr))
    out_data['time'] = pyfuncs.numpy_to_dataarray(time, dims=['time'], attrs=time_atts)

    # Data set to be written out, base variables for all models
    ds = xr.Dataset(
        {
            "PS": out_data['ps'],
            "U": out_data['u'],
            "V": out_data['v'],
            "T": out_data['t'],
            "Q": out_data['q'],
            "lat": out_data['lat'],
            "lon": out_data['lon'],
        }
    )

    if dycore == "se" or dycore == "fv":
        # Reload from the template xarray for metadata purposes
        hya, hyb, hyai, hybi, lev, ilev = loaddata.load_cam_levels(PATHTOHERE, numlevels, load_xarray = True)
        ds["CLDLIQ"] = out_data['cldliq']
        ds["CLDICE"] = out_data['cldice']
        ds["hyam"] = hya
        ds["hybm"] = hyb
        ds["hyai"] = hyai
        ds["hybi"] = hybi
        ds["lev"] = lev
        ds["ilev"] = ilev
        ds["time"] = time

    if dycore == "fv":
        # Add staggered variables for fv
        ds["US"] = out_data['us']
        ds["VS"] = out_data['vs']
        ds["slat"] = out_data['slat']
        ds["slon"] = out_data['slon']

    # Output correct_or_not if available
    if 'correct_or_not' in locals():
        ds["correct_or_not"] = out_data['correct_or_not']

    # Add global attributes
    ds.attrs.update({
        "title": "Betacast-generated ncdata file",
        "source_file": data_filename,
        "wgt_file": wgt_filename,
        "init_date": YYYYMMDDHH,
        "creation_date": str(np.datetime64('now')),
        "dycore": dycore,
        "datasource": datasource
    })

    # Turn off fill value in relevant coordinate variables
    vars_to_check = ["hyam", "hybm", "hyai", "hybi", "lev", "ilev", "time"]
    existing_vars = [var for var in vars_to_check if var in ds.variables]
    encoding = {var: {'_FillValue': None} for var in existing_vars}

    # Determine file format and compression based on variable size
    if compress_file:
        for var in list(ds.data_vars) + list(ds.coords):
            if var in encoding:  # update it
                encoding[var].update({"zlib": True, "complevel": 1})
            else:   # add it
                encoding[var] = {"zlib": True, "complevel": 1}
        ds.to_netcdf(se_inic, format="NETCDF4_CLASSIC", unlimited_dims=["time"], encoding=encoding)
    else:
        var_max_size = pyfuncs.print_and_return_varsize(ds["PS"], ds["U"], ds["V"], ds["T"], ds["Q"])
        netcdf_format = "NETCDF4" if var_max_size >= 4e9 else "NETCDF3_64BIT"
        ds.to_netcdf(se_inic, format=netcdf_format, unlimited_dims=["time"], encoding=encoding)


if __name__ == "__main__":

    main()