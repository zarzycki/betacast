# Core modules
import os
import sys
import logging

# Third party modules
import numpy as np
import xarray as xr
import netCDF4 as nc

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
import py_seedfuncs
from constants import (
    p0, NC_FLOAT_FILL, dtime_map, QMINTHRESH, QMAXTHRESH, CLDMINTHRESH,
    ps_wet_to_dry, output_diag, w_smooth_iter, damp_upper_winds_mpas, grav
)

args = pyfuncs.parse_args()

pyfuncs.configure_logging(args.verbose)

# Set nc fill values
nc.default_fillvals['f4'] = NC_FLOAT_FILL
nc.default_fillvals['f8'] = float(NC_FLOAT_FILL)

def main():

    logger = logging.getLogger(__name__)
    logger.info("Main function started")

    # Is the Betacast path available to us?
    BETACAST, PATHTOHERE = pyfuncs.get_betacast_path()
    TEMPLATESPATH = os.path.join(BETACAST, 'grids/templates/')
    VERTCOORDSPATH = os.path.join(BETACAST, 'grids/vert-coords/')

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
    augment_tcs = args.augment_tcs
    vortex_namelist = args.vortex_namelist

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

    if (dycore == 'mpas'):
        load_vert_templates = False
    else:
        load_vert_templates = True
    logging.info(f"load_vert_templates: {load_vert_templates}")

    # Toggle whether the output streams will be floats or doubles
    write_type = "float" if write_floats else "double"
    logging.info(f"Output type set to: {write_type}")

    # ===== Getting date from YYYYMMDDHH
    yearstr, monthstr, daystr, cyclestr = pyfuncs.split_by_lengths(str(YYYYMMDDHH), dtime_map)

    logging.info(f"Regridding analysis from: {yearstr} {monthstr} {daystr} {cyclestr}Z")

    logging.info("---------------------------------------------------------")
    logging.info(f"Using this file: {data_filename}")
    logging.info(f"Using this remap: {wgt_filename}")

    if load_vert_templates:
        hya, hyb, hyai, hybi, lev, ilev = loaddata.load_cam_levels(TEMPLATESPATH, numlevels)

    # IMPORTANT! data_vars should be organized top-to-bottom when loaddata returns
    # (i.e., lowest pressure/highest z at 0 index of lev)
    # I attempt to account for this elsewhere in the code with flips, but make no promises
    if datasource == 'GFS' or datasource == 'HWRF' or datasource == 'HRRR':   # NCEP pressure level data
        data_vars = loaddata.load_CFSR_data(data_filename, dycore)
    elif datasource == 'HRRRml':   # HRRR model level data
        data_vars = loaddata.load_HRRRml_data(data_filename, dycore)
    elif datasource == 'SAMPLE':
        data_vars = loaddata.load_SAMPLE_data(data_filename, dycore)
    elif datasource == 'ERA5RDA':
        data_vars = loaddata.load_ERA5RDA_data(RDADIR, data_filename, yearstr, monthstr, daystr, cyclestr, dycore)
    elif datasource == 'ERA5mlRDA':
        data_vars = loaddata.load_ERA5mlRDA_data(RDADIR, data_filename, VERTCOORDSPATH, yearstr, monthstr, daystr, cyclestr, dycore)
        pyfuncs.log_resource_usage()
    elif datasource == 'CAM':
        data_vars = loaddata.load_cam_data(data_filename, YYYYMMDDHH, mod_in_topo, mod_remap_file, dycore, write_debug_files=write_debug_files,write_debug_dir=DEBUGDIR)

    logging.info("Input Data Level information")
    logging.info(f"Number: {len(data_vars['lev'])}")
    logging.info(f"Max: {np.max(data_vars['lev'])}")
    logging.info(f"Min: {np.min(data_vars['lev'])}")
    logging.info(f"Input Coords: nlat: {len(data_vars['lat'])}, nlon: {len(data_vars['lon'])}")

    if dycore == "mpas":
        # Create a 3-D pressure field from constant pressure surfaces
        data_vars['pres'] = np.broadcast_to(data_vars['lev'][:, np.newaxis, np.newaxis], data_vars['t'].shape)

        # If z is reported in geopotential, convert to geometric height
        if data_vars['z_is_phi']:
            data_vars['z'] = data_vars['z'] / grav

        # Calculate pressure interfaces as midpoints between levels
        # NOTE: pint does NOT include any surface pressure, so it's the interface between bottom, 2nd-from-bottom plev
        interfaces = (data_vars['lev'][:-1] + data_vars['lev'][1:]) / 2
        interfaces = np.insert(interfaces, 0, data_vars['lev'][0] / 2)  # Add top interface (half of first level)
        data_vars['pint'] = np.broadcast_to(interfaces[:, np.newaxis, np.newaxis],data_vars['pres'].shape)

    # Print diagnostics
    pyfuncs.log_resource_usage()
    pyfuncs.print_min_max_dict(data_vars)

    # TC augmentation code
    # How logic works...
    # If augment_tcs is true, use internal TCvitals file to take existing vortex and deepen (or weaken)
    # if vortex_namelist is provided and is valid, read a vortex_namelist file and do seeding/unseeding from scratch
    # augment_tcs takes priority
    if augment_tcs or (vortex_namelist and os.path.isfile(vortex_namelist)):
        logging.info(f"Starting TC operation...")

        if augment_tcs and (vortex_namelist and os.path.isfile(vortex_namelist)):
            logging.info(f"WARNING: Both TC augmentation and vortex_namelist are TRUE")
            logging.info(f"WARNING: TC augmentation prioritized, remove augment_tcs if this isn't what you want")

        if augment_tcs:
            # From the TC vitals file, figure out what we need
            tcs_at_this_time = py_seedfuncs.parse_tcvitals(
                os.path.join(BETACAST, 'cyclone-tracking/fin-tcvitals/combined/ALL_combined_tcvitals.dat.gz'),
                YYYYMMDDHH
                )
            num_tcs = len(tcs_at_this_time['lats'])
            logging.info(f"augment_tcs: {augment_tcs}; dealing with {num_tcs} TCs.")
        elif vortex_namelist and os.path.isfile(vortex_namelist):
            num_tcs = 1
            logging.info(f"We have a valid vortex_namelist: {vortex_namelist}; dealing with {num_tcs} TCs.")
        else:
            logging.info("Cannot seed/unseed without {augment_tcs} or {vortex_namelist}")
            num_tcs = 0

        # Check if we found any valid storms
        if num_tcs > 0:
            # Loop over each storm
            for storm_idx in range(num_tcs):
                if augment_tcs:
                    storm_name = tcs_at_this_time['storm_names'][storm_idx]
                    storm_lat = tcs_at_this_time['lats'][storm_idx]
                    storm_lon = tcs_at_this_time['lons'][storm_idx]
                    storm_pressure = tcs_at_this_time['pressures'][storm_idx]
                    storm_rmw = tcs_at_this_time['rmw'][storm_idx]

                    logging.info(f"Processing storm {storm_idx + 1}/{len(tcs_at_this_time['lats'])}: {storm_name}")
                    logging.info(f"  Location: {storm_lat}, {storm_lon}")
                    logging.info(f"  Pressure: {storm_pressure} mb, RMW: {storm_rmw} km")

                    # Skip storms with missing critical data
                    if storm_pressure == -999 or storm_rmw == -999:
                        logging.warning(f"Skipping {storm_name} due to missing pressure or RMW data")
                        continue

                    # Update input_dict with storm-specific values
                    input_dict = {
                        'deltaMax': 0.25,
                        'target_rmw': float(storm_rmw)*1000.,        # Use RMW from TCVitals
                        'minp': float(storm_pressure),         # Use pressure from TCVitals
                        'psminlat': float(storm_lat),          # Use latitude from TCVitals
                        'psminlon': float(storm_lon)           # Use longitude from TCVitals
                    }
                else:
                    # We are going to just straight up read the namelist
                    storm_name = "Synth"
                    input_dict = vortex_namelist

                # Get TC seeding parameters for this storm
                tc_seed_params = py_seedfuncs.read_tcseed_settings(input_dict)
                logging.info(f"Initial parameters for {storm_name}: {tc_seed_params}")

                if tc_seed_params['invert_vortex'] or augment_tcs:
                    # Find/optimize parameters if needed
                    logging.info(f"Optimizing missing TC parameters for {storm_name}...")
                    tc_seed_params = py_seedfuncs.find_fill_parameters(
                        data_vars['ps'],
                        data_vars['t'],
                        data_vars['lat'],
                        data_vars['lon'],
                        data_vars['lev'],
                        tc_seed_params
                    )
                    logging.info(f"Updated parameters for {storm_name}: {tc_seed_params}")

                # Apply TC seeding/unseeding for this storm
                logging.info(f"Applying TC seeding/unseeding for {storm_name}...")
                if dycore == "mpas":
                    data_vars = py_seedfuncs.apply_tc_seeding(
                        data_vars,
                        tc_seed_params,
                        update_z=True,
                        add_deltas_to_data=write_debug_files
                    )
                else:
                    data_vars = py_seedfuncs.apply_tc_seeding(
                        data_vars,
                        tc_seed_params,
                        update_z=False,
                        add_deltas_to_data=write_debug_files
                    )

                logging.info(f"Completed processing for {storm_name}")

            logging.info(f"Finished processing all {num_tcs} storms")

            if write_debug_files:
                pyfuncs.print_debug_file(
                    DEBUGDIR + "/py_era5_tcseed.nc",
                    ps_cam=(["lat", "lon"], data_vars['ps']),
                    ps_vx_cam=(["lat", "lon"], data_vars['ps_vx']),
                    t_cam=(["lev_p", "lat", "lon"], data_vars['t']),
                    t_vx_cam=(["lev_p", "lat", "lon"], data_vars['t_vx']),
                    u_cam=(["lev_p", "lat", "lon"], data_vars['u']),
                    u_vx_cam=(["lev_p", "lat", "lon"], data_vars['u_vx']),
                    v_cam=(["lev_p", "lat", "lon"], data_vars['v']),
                    v_vx_cam=(["lev_p", "lat", "lon"], data_vars['v_vx']),
                    q_cam=(["lev_p", "lat", "lon"], data_vars['q']),
                    q_vx_cam=(["lev_p", "lat", "lon"], data_vars['q_vx']),
                    lat=(["lat"], data_vars['lat']),
                    lon=(["lon"], data_vars['lon']),
                    **({
                        'z_cam': (["lev_p", "lat", "lon"], data_vars['z']),
                        'z_vx_cam': (["lev_p", "lat", "lon"], data_vars['z_vx'])
                    } if dycore == 'mpas' else {})
                )
                del data_vars['ps_vx']
                del data_vars['t_vx']
                del data_vars['u_vx']
                del data_vars['v_vx']
                del data_vars['q_vx']
                if dycore == 'mpas':
                    del data_vars['z_vx']

        else:
            logging.info("No valid storms found in TCVitals data - skipping TC augmentation")
            logging.info("Proceeding without TC modifications")

    else:
        logging.info("No TC modification requested")

    # Do vertical interpolation here for pressure, sigma, and hybrid targets
    if dycore == 'fv' or dycore == 'se':

        if write_debug_files:
            debug_vars = {
                "ps_cam": (["lat", "lon"], data_vars["ps"]),
                "t_cam": (["lev_p", "lat", "lon"], data_vars["t"]),
                "u_cam": (["lev_p", "lat", "lon"], data_vars["u"]),
                "v_cam": (["lev_p", "lat", "lon"], data_vars["v"]),
                "q_cam": (["lev_p", "lat", "lon"], data_vars["q"]),
                "cldliq_cam": (["lev_p", "lat", "lon"], data_vars["cldliq"]),
                "cldice_cam": (["lev_p", "lat", "lon"], data_vars["cldice"]),
            }
            if datasource not in ('HRRRml', 'HRRR', 'HWRF', 'RAP'):
                debug_vars["lat"] = (["lat"], data_vars["lat"])
                debug_vars["lon"] = (["lon"], data_vars["lon"])
            pyfuncs.print_debug_file(DEBUGDIR + "/py_era5_before_interp.nc", **debug_vars)

        data_vint = vertremap.pres2hyb_all(data_vars, data_vars['ps'], hya, hyb)

        if write_debug_files:
            debug_vars = {
                "ps_cam": (["latitude", "longitude"], data_vint["ps"]),
                "t_cam": (["level", "latitude", "longitude"], data_vint["t"]),
                "u_cam": (["level", "latitude", "longitude"], data_vint["u"]),
                "v_cam": (["level", "latitude", "longitude"], data_vint["v"]),
                "q_cam": (["level", "latitude", "longitude"], data_vint["q"]),
                "cldliq_cam": (["level", "latitude", "longitude"], data_vint["cldliq"]),
                "cldice_cam": (["level", "latitude", "longitude"], data_vint["cldice"]),
            }
            if datasource not in ('HRRRml', 'HRRR', 'HWRF', 'RAP'):
                debug_vars["latitude"] = (["latitude"], data_vint["lat"])
                debug_vars["longitude"] = (["longitude"], data_vint["lon"])
            pyfuncs.print_debug_file(DEBUGDIR + "/py_era5_on_hybrid.nc", **debug_vars)

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

        # These things we want to do as close to interpolation as possible
        # Initialize variables
        data_vars['theta'] = data_vars['t'].copy()
        data_vars['rho'] = data_vars['t'].copy()

        # Calculate potential temperature using full pressure and actual T
        data_vars['theta'] = meteo.pot_temp(data_vars['pres'], data_vars['t'])

        # Calculate dry air density from pressure, moisture, and temperature
        data_vars['rho'] = meteo.calculate_rho_dry(data_vars['pres'], data_vars['q'], data_vars['t'])

        # If w is omega (Pa/s), convert to vertical velocity in m/s
        if data_vars['w_is_omega']:
            data_vars['w'] = meteo.omega_to_w(data_vars['w'], data_vars['pres'], data_vars['t'])

        # Smooth w if needed
        if w_smooth_iter > 0:
            #data_vars['w'] = pyfuncs.smooth_with_gaussian(data_vars['w'], w_smooth_iter)
            data_vars['w'] = pyfuncs.smooth_with_smth9(data_vars['w'], w_smooth_iter)

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
            logging.info(f"Writing MPAS file: {se_inic}...")
            write_mpas_xarray=False
            if write_mpas_xarray:  # Use the open xarray dataset to overwrite
                mpas_file['u'].values[0, :, :] = data_horiz['uNorm'].T
                mpas_file['qv'].values[0, :, :] = data_horiz['q'].T
                mpas_file['rho'].values[0, :, :] = data_horiz['rho'].T
                mpas_file['theta'].values[0, :, :] = data_horiz['theta'].T
                mpas_file['w'].values[0, :, :] = 0.06 * data_horiz['w'].T
                logging.info(f"Beginning actual file write")
                mpas_file.to_netcdf(se_inic, format='NETCDF4')
                logging.info(f"Closing the xr file")
                mpas_file.close()
            else:   # Create a new NetCDF file with the same structure

                # Get dimensions from the xarray dataset
                dims = mpas_file.sizes

                # Close the xr file now that we have everything we need
                logging.info(f"Closing the xr file")
                mpas_file.close()

                # Open src + target MPAS files simultaneously using NetCDF4
                with nc.Dataset(mpasfile, mode='r') as src_file, nc.Dataset(se_inic, mode='w', format='NETCDF4') as new_file:

                    # Copy global file attributes from src -> target
                    for attr_name in src_file.ncattrs():
                        logging.debug(f"... copying global attr {attr_name}")
                        new_file.setncattr(attr_name, src_file.getncattr(attr_name))

                    # Create dimensions - all with fixed size (no unlimited for init conditions)
                    for dim_name, dim in src_file.dimensions.items():
                        dim_size = len(dim)
                        logging.debug(f"... creating dim {dim_name}")
                        new_file.createDimension(dim_name, dim_size)

                    # Create and copy variables
                    total_mpas_vars = len(src_file.variables)
                    logging.info(f"Total variables to process: {total_mpas_vars}")

                    # First create all the variables (metadata only)
                    for var_name, var in src_file.variables.items():
                        logging.info(f"Creating variable {var_name}")

                        # Special handling for Time variable
                        if var_name == 'Time':
                            var_type = 'f8'  # Use double for Time
                        # Handle string variables correctly
                        elif hasattr(var, 'dtype') and (var.dtype.kind == 'S' or var.dtype.kind == 'U'):
                            if var.dimensions and 'StrLen' in var.dimensions:
                                var_type = 'S1'  # Char type
                            else:
                                var_type = str  # String type
                        else:
                            var_type = var.dtype

                        # Create the variable with contiguous storage (no chunking, which slows down PIO)
                        new_var = new_file.createVariable(
                            var_name,
                            var_type,
                            var.dimensions,
                            fill_value=getattr(var, '_FillValue', None),
                            contiguous=True
                        )

                        # Copy variable attributes except storage-related ones
                        for attr_name in var.ncattrs():
                            if attr_name not in ['_FillValue', '_Storage', '_ChunkSizes']:
                                try:
                                    setattr(new_var, attr_name, var.getncattr(attr_name))
                                except Exception as e:
                                    logging.warning(f"Error setting attribute {attr_name} for {var_name}: {e}")

                    # Define a list of variables to modify later
                    skip_vars = ['u', 'qv', 'rho', 'theta', 'w']

                    # Now copy data, one variable at a time
                    for idx, var_name in enumerate(src_file.variables, 1):
                        if var_name in skip_vars:
                            logging.info(f"Skipping data copy for {var_name} ({idx}/{total_mpas_vars}) - will set manually later")
                            continue

                        src_var = src_file.variables[var_name]
                        dst_var = new_file.variables[var_name]

                        logging.info(f"Copying data for {var_name} ({idx}/{total_mpas_vars})")

                        # Handle string variables specially
                        if hasattr(src_var, 'dtype') and (src_var.dtype.kind == 'S' or src_var.dtype.kind == 'U'):
                            if len(src_var.dimensions) == 1 and isinstance(dst_var.dtype, str):
                                # This is a single string
                                dst_var[0] = str(src_var[0])
                            elif 'StrLen' in src_var.dimensions:
                                # This is a char array
                                dst_var[:] = src_var[:]
                            else:
                                # Handle 1D array of strings
                                for i in range(len(src_var)):
                                    dst_var[i] = str(src_var[i])
                        else:
                            # Handle large variables in chunks to reduce memory usage
                            if src_var.size > 10_000_000 and len(src_var.shape) > 1:
                                shape = src_var.shape
                                chunk_size = max(1, min(100, shape[0] // 10))

                                for i in range(0, shape[0], chunk_size):
                                    end = min(i + chunk_size, shape[0])
                                    logging.debug(f"... chunk {i}:{end}")
                                    dst_var[i:end] = src_var[i:end]
                            else:
                                # Small enough variable to copy directly
                                dst_var[:] = src_var[:]

                        logging.info(f"Finished copying {var_name}")

                    # Write out Betacast modified vars last
                    logging.info(f"Writing new MPAS fields to target netcdf file")
                    new_file.variables['u'][0, :, :] = data_horiz['uNorm'].T
                    new_file.variables['qv'][0, :, :] = data_horiz['q'].T
                    new_file.variables['rho'][0, :, :] = data_horiz['rho'].T
                    new_file.variables['theta'][0, :, :] = data_horiz['theta'].T
                    new_file.variables['w'][0, :, :] = 0.06 * data_horiz['w'].T

                logging.info(f"Done generating MPAS initial condition file: {se_inic}, exiting...")

            sys.exit(0)

    grid_dims = data_horiz['ps'].shape
    if dycore == "se" or dycore == "mpas":
        nncol = grid_dims
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

    add_pmid = False
    if add_pmid:
        data_horiz['pmid'] = meteo.compute_pmid(data_horiz['t'],data_horiz['rho'])
        pyfuncs.print_min_max_dict(data_horiz)

    if datasource == "HWRF":
        logging.info(f"{datasource} replacing nan with _FillValue {NC_FLOAT_FILL} since regional")
        data_horiz = pyfuncs.replace_nans_with_fill_value(data_horiz, ['ps', 'u', 'v', 't', 'q', 'cldliq', 'cldice'], NC_FLOAT_FILL)

    logging.info(f"Begin writing output file: {se_inic}")

    write_xarray = False

    if write_xarray:
        out_data = {}

        if dycore == "se":
            out_data['ps'] = pyfuncs.numpy_to_dataarray(data_horiz['ps'], dims=['ncol'], attrs={'units': 'Pa', "_FillValue": NC_FLOAT_FILL})
            out_data['u'] = pyfuncs.numpy_to_dataarray(data_horiz['u'], dims=['lev', 'ncol'], attrs={'units': 'm/s', "_FillValue": NC_FLOAT_FILL})
            out_data['v'] = pyfuncs.numpy_to_dataarray(data_horiz['v'], dims=['lev', 'ncol'], attrs={'units': 'm/s', "_FillValue": NC_FLOAT_FILL})
            out_data['t'] = pyfuncs.numpy_to_dataarray(data_horiz['t'], dims=['lev', 'ncol'], attrs={'units': 'K', "_FillValue": NC_FLOAT_FILL})
            out_data['q'] = pyfuncs.numpy_to_dataarray(data_horiz['q'], dims=['lev', 'ncol'], attrs={'units': 'kg/kg', "_FillValue": NC_FLOAT_FILL})
            if add_cloud_vars:
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
            if add_cloud_vars:
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
            if add_cloud_vars:
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
            if add_cloud_vars:
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
            hya, hyb, hyai, hybi, lev, ilev = loaddata.load_cam_levels(TEMPLATESPATH, numlevels, load_xarray = True)
            if add_cloud_vars:
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
        logging.debug("Turning off _FillValue in relevant vars")
        vars_to_check = ["hyam", "hybm", "hyai", "hybi", "lev", "ilev", "time"]
        existing_vars = [var for var in vars_to_check if var in ds.variables]
        encoding = {var: {'_FillValue': None} for var in existing_vars}

        logging.info(f"Writing {se_inic}")
        # Determine file format and compression based on variable size
        if compress_file:
            for var in list(ds.data_vars) + list(ds.coords):
                if var in encoding:  # update it
                    encoding[var].update({"zlib": True, "complevel": 1})
                else:   # add it
                    encoding[var] = {"zlib": True, "complevel": 1}
            ds.load().to_netcdf(se_inic, format="NETCDF4_CLASSIC", unlimited_dims=["time"], encoding=encoding)
        else:
            var_max_size = pyfuncs.print_and_return_varsize(ds["PS"], ds["U"], ds["V"], ds["T"], ds["Q"])
            netcdf_format = "NETCDF4" if var_max_size >= 4e9 else "NETCDF3_64BIT"
            ds.load().to_netcdf(se_inic, format=netcdf_format, unlimited_dims=["time"], encoding=encoding)

    else:

        def replace_nans_with_fill(data, fill_value=-9999.0):
            return np.where(np.isnan(data), fill_value, data)

        if compress_file:
            netcdf_format = "NETCDF4_CLASSIC"
            compression_opts = {'zlib': True, 'complevel': 1}
        else:
            var_max_size = pyfuncs.print_and_return_varsize(data_horiz["ps"], data_horiz["u"], data_horiz["v"], data_horiz["t"], data_horiz["q"])
            netcdf_format = "NETCDF4" if var_max_size >= 4e9 else "NETCDF3_64BIT"
            compression_opts = {}

        nc_file = nc.Dataset(se_inic, 'w', format=netcdf_format)

        nc_file.createDimension('time', None)  # None makes the dimension unlimited

        time_nc = nc_file.createVariable('time', 'f8', ('time',), **compression_opts)
        time, time_atts = pyfuncs.create_cf_time(int(yearstr), int(monthstr), int(daystr), int(cyclestr))
        for attr, value in time_atts.items():
            setattr(time_nc, attr, value)
        time_nc[:] = time

        nc_file.createDimension('lev', numlevels)
        nc_file.createDimension('ilev', numlevels + 1)

        # Horiz coordinates
        if dycore == "se" or dycore == "mpas":
            nc_file.createDimension('ncol', nncol[0])
            lat_nc = nc_file.createVariable('lat', 'f8', ('ncol',), fill_value=-900., **compression_opts)
            lon_nc = nc_file.createVariable('lon', 'f8', ('ncol',), fill_value=-900., **compression_opts)
            lat_nc[:] = data_horiz['lat']
            lon_nc[:] = data_horiz['lon']
        elif dycore == "fv":
            nc_file.createDimension('lat', nfvlat)
            nc_file.createDimension('lon', nfvlon)
            lat_nc = nc_file.createVariable('lat', 'f8', ('lat',), fill_value=-900., **compression_opts)
            lon_nc = nc_file.createVariable('lon', 'f8', ('lon',), fill_value=-900., **compression_opts)
            lat_nc[:] = data_horiz['lat']
            lon_nc[:] = data_horiz['lon']

        # Create nc variables for state fields
        ps_nc = nc_file.createVariable('PS', 'f4', ('time', 'ncol') if dycore == "se" or dycore == "mpas" else ('time', 'lat', 'lon'), fill_value=NC_FLOAT_FILL, **compression_opts)
        u_nc = nc_file.createVariable('U', 'f4', ('time', 'lev', 'ncol') if dycore == "se" or dycore == "mpas" else ('time', 'lev', 'lat', 'lon'), fill_value=NC_FLOAT_FILL, **compression_opts)
        v_nc = nc_file.createVariable('V', 'f4', ('time', 'lev', 'ncol') if dycore == "se" or dycore == "mpas" else ('time', 'lev', 'lat', 'lon'), fill_value=NC_FLOAT_FILL, **compression_opts)
        t_nc = nc_file.createVariable('T', 'f4', ('time', 'lev', 'ncol') if dycore == "se" or dycore == "mpas" else ('time', 'lev', 'lat', 'lon'), fill_value=NC_FLOAT_FILL, **compression_opts)
        q_nc = nc_file.createVariable('Q', 'f4', ('time', 'lev', 'ncol') if dycore == "se" or dycore == "mpas" else ('time', 'lev', 'lat', 'lon'), fill_value=NC_FLOAT_FILL, **compression_opts)
        ps_nc.units = 'Pa'
        u_nc.units = 'm/s'
        v_nc.units = 'm/s'
        t_nc.units = 'K'
        q_nc.units = 'kg/kg'
        if dycore != "mpas" and add_cloud_vars:
            cldliq_nc = nc_file.createVariable('CLDLIQ', 'f4', ('time', 'lev', 'ncol') if dycore == "se" or dycore == "mpas" else ('time', 'lev', 'lat', 'lon'), fill_value=NC_FLOAT_FILL, **compression_opts)
            cldice_nc = nc_file.createVariable('CLDICE', 'f4', ('time', 'lev', 'ncol') if dycore == "se" or dycore == "mpas" else ('time', 'lev', 'lat', 'lon'), fill_value=NC_FLOAT_FILL, **compression_opts)
            cldliq_nc.units = "kg/kg"
            cldice_nc.units = "kg/kg"
        if 'correct_or_not' in locals():
            correct_or_not_nc = nc_file.createVariable('correct_or_not', 'f4', ('time', 'ncol') if dycore == "se" or dycore == "mpas" else ('time', 'lat', 'lon'), fill_value=-1.0, **compression_opts)
        if add_pmid:
            pmid_nc = nc_file.createVariable('PMID', 'f4', ('time', 'lev', 'ncol') if dycore == "se" or dycore == "mpas" else ('time', 'lev', 'lat', 'lon'), fill_value=NC_FLOAT_FILL, **compression_opts)
        # Place data_horiz data into var fields
        if dycore == "fv":
            ps_nc[0, :, :] = replace_nans_with_fill(data_horiz['ps'],fill_value=NC_FLOAT_FILL)
            u_nc[0, :, :, :] = replace_nans_with_fill(data_horiz['u'],fill_value=NC_FLOAT_FILL)
            v_nc[0, :, :, :] = replace_nans_with_fill(data_horiz['v'],fill_value=NC_FLOAT_FILL)
            t_nc[0, :, :, :] = replace_nans_with_fill(data_horiz['t'],fill_value=NC_FLOAT_FILL)
            q_nc[0, :, :, :] = replace_nans_with_fill(data_horiz['q'],fill_value=NC_FLOAT_FILL)
            if add_cloud_vars:
                cldliq_nc[0, :, :, :] = replace_nans_with_fill(data_horiz['cldliq'],fill_value=NC_FLOAT_FILL)
                cldice_nc[0, :, :, :] = replace_nans_with_fill(data_horiz['cldice'],fill_value=NC_FLOAT_FILL)
            if 'correct_or_not' in locals():
                correct_or_not_nc[0, :, :] = replace_nans_with_fill(data_horiz['correct_or_not'],fill_value=NC_FLOAT_FILL)
        elif dycore == "se":
            ps_nc[0, :] = replace_nans_with_fill(data_horiz['ps'],fill_value=NC_FLOAT_FILL)
            u_nc[0, :, :] = replace_nans_with_fill(data_horiz['u'],fill_value=NC_FLOAT_FILL)
            v_nc[0, :, :] = replace_nans_with_fill(data_horiz['v'],fill_value=NC_FLOAT_FILL)
            t_nc[0, :, :] = replace_nans_with_fill(data_horiz['t'],fill_value=NC_FLOAT_FILL)
            q_nc[0, :, :] = replace_nans_with_fill(data_horiz['q'],fill_value=NC_FLOAT_FILL)
            if add_cloud_vars:
                cldliq_nc[0, :, :] = replace_nans_with_fill(data_horiz['cldliq'],fill_value=NC_FLOAT_FILL)
                cldice_nc[0, :, :] = replace_nans_with_fill(data_horiz['cldice'],fill_value=NC_FLOAT_FILL)
            if 'correct_or_not' in locals():
                correct_or_not_nc[0, :] = replace_nans_with_fill(data_horiz['correct_or_not'],fill_value=NC_FLOAT_FILL)
        elif dycore == "mpas":
            ps_nc[0, :] = replace_nans_with_fill(data_horiz['ps'],fill_value=NC_FLOAT_FILL)
            u_nc[0, :, :] = replace_nans_with_fill(data_horiz['u'][::-1, :],fill_value=NC_FLOAT_FILL)
            v_nc[0, :, :] = replace_nans_with_fill(data_horiz['v'][::-1, :],fill_value=NC_FLOAT_FILL)
            t_nc[0, :, :] = replace_nans_with_fill(data_horiz['t'][::-1, :],fill_value=NC_FLOAT_FILL)
            q_nc[0, :, :] = replace_nans_with_fill(data_horiz['q'][::-1, :],fill_value=NC_FLOAT_FILL)
            if 'correct_or_not' in locals():
                correct_or_not_nc[0, :] = replace_nans_with_fill(data_horiz['correct_or_not'],fill_value=NC_FLOAT_FILL)
            if add_pmid:
                pmid_nc[0, :, :] = replace_nans_with_fill(data_horiz['pmid'][::-1, :],fill_value=NC_FLOAT_FILL)

        # If the model has a hybrid coordinate, do that here
        if dycore == "se" or dycore == "fv":
            hya, hyb, hyai, hybi, lev, ilev = loaddata.load_cam_levels(TEMPLATESPATH, numlevels, load_xarray=False)

            hyam_nc = nc_file.createVariable('hyam', 'f8', ('lev',), fill_value=NC_FLOAT_FILL, **compression_opts)
            hybm_nc = nc_file.createVariable('hybm', 'f8', ('lev',), fill_value=NC_FLOAT_FILL, **compression_opts)
            hyai_nc = nc_file.createVariable('hyai', 'f8', ('ilev',), fill_value=NC_FLOAT_FILL, **compression_opts)
            hybi_nc = nc_file.createVariable('hybi', 'f8', ('ilev',), fill_value=NC_FLOAT_FILL, **compression_opts)
            lev_nc = nc_file.createVariable('lev', 'f8', ('lev',), fill_value=NC_FLOAT_FILL, **compression_opts)
            ilev_nc = nc_file.createVariable('ilev', 'f8', ('ilev',), fill_value=NC_FLOAT_FILL, **compression_opts)

            # Let's not allow nans here, so no filling
            hyam_nc[:] = hya
            hybm_nc[:] = hyb
            hyai_nc[:] = hyai
            hybi_nc[:] = hybi
            lev_nc[:] = lev
            ilev_nc[:] = ilev

        # If FV, add staggered information
        if dycore == "fv":
            nc_file.createDimension('slat', nfvlat - 1)
            nc_file.createDimension('slon', nfvlon)

            slat_nc = nc_file.createVariable('slat', 'f8', ('slat',), fill_value=-900., **compression_opts)
            slon_nc = nc_file.createVariable('slon', 'f8', ('slon',), fill_value=-900., **compression_opts)
            us_nc = nc_file.createVariable('US', 'f4', ('time', 'lev', 'slat', 'lon'), fill_value=NC_FLOAT_FILL, **compression_opts)
            vs_nc = nc_file.createVariable('VS', 'f4', ('time', 'lev', 'lat', 'slon'), fill_value=NC_FLOAT_FILL, **compression_opts)

            slat_nc[:] = data_horiz['fvslat']
            slon_nc[:] = data_horiz['fvslon']
            us_nc[0, :, :, :] = replace_nans_with_fill(data_horiz['us'],fill_value=NC_FLOAT_FILL)
            vs_nc[0, :, :, :] = replace_nans_with_fill(data_horiz['vs'],fill_value=NC_FLOAT_FILL)

        # Add global attributes
        nc_file.title = "Betacast-generated ncdata file"
        nc_file.source_file = data_filename
        nc_file.wgt_file = wgt_filename
        nc_file.init_date = YYYYMMDDHH
        nc_file.creation_date = str(np.datetime64('now'))
        nc_file.dycore = dycore
        nc_file.datasource = datasource

        # Close the file
        nc_file.close()


if __name__ == "__main__":

    main()