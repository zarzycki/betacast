import numpy as np
import xarray as xr
import datetime
import argparse
import sys
import glob
import cftime
from scipy.ndimage import gaussian_filter
from timing import start_time, print_elapsed_time

from vertremap import *
from cam2cam import *
import horizremap
from hy2pres import interp_hybrid_to_pressure

from constants import (
    grav, kappa, p0, Rd, Rv_over_Rd, t_freeze_K, dtime_map, OMEGA_LAT_THRESH
)

def parse_args():
    parser = argparse.ArgumentParser(description="Process command-line arguments for climate data processing.")

    # Required arguments
    parser.add_argument('--datasource', type=str, required=True,
                        help='Data source: CFSR, ERA5, ERAI, GFS, CAM, RAP')
    parser.add_argument('--numlevels', type=int, required=True,
                        help='Number of vertical levels')
    parser.add_argument('--YYYYMMDDHH', type=int, required=True,
                        help='Initialization date in YYYYMMDDHH format')
    parser.add_argument('--data_filename', type=str, required=True,
                        help='Full path to file containing initial information')
    parser.add_argument('--wgt_filename', type=str, required=True,
                        help='Full path to ESMF weight file from ANL -> MOD')
    parser.add_argument('--se_inic', type=str, required=True,
                        help='Full path of file to write as output')

    # Optional arguments with defaults
    parser.add_argument('--dycore', type=str, default='se',
                        help='Dycore option: fv, se, mpas (default: se)')
    parser.add_argument('--RDADIR', type=str, default='',
                        help='Base path to RDA folder (default: "")')
    parser.add_argument('--mpas_as_cam', action='store_true',
                        help='If set, write MPAS in CAM physics output (default: False)')
    parser.add_argument('--compress_file', action='store_true',
                        help='If set, will attempt NetCDF "chunking" compression (default: False)')
    parser.add_argument('--write_floats', action='store_true',
                        help='If set, write outputs as single precision (default: False)')
    parser.add_argument('--add_cloud_vars', action='store_true', default=True,
                        help='If set, add CLDICE and CLDLIQ to output file (default: True)')
    parser.add_argument('--adjust_config', type=str, default=' ',
                        help='String defining how to perform hydro adjustment (default: "")')
    parser.add_argument('--model_topo_file', type=str, default='',
                        help='File containing PHIS for FV or SE, or an MPAS inic file (default: "")')
    parser.add_argument('--mod_in_topo', type=str, default='',
                        help='Full path to PHIS field from downscaling MOD for pressure surface calculation (default: "")')
    parser.add_argument('--mod_remap_file', type=str, default='',
                        help='Full path to ESMF weight file that goes downscaling MOD -> ANL (default: "")')
    parser.add_argument('--mpasfile', type=str, default='',
                        help='XXXXXXX')
    parser.add_argument('-v', '--verbose', action='store_true', help='Increase output verbosity for debugging')

    return parser.parse_args()

def split_by_lengths(s, lengths):
    """
    Splits a string `s` into parts defined by `lengths`.
    """
    parts = []
    index = 0
    for length in lengths:
        parts.append(s[index:index + length])
        index += length
    return parts

def process_grb_file(grb_file, datasource, RDADIR, yearstr, monthstr, daystr, mod_remap_file):
    print("---------------------------------------------------------")
    print("Loading lat/lon/lev")

    if datasource in ["GFS", "CFSR", "HWRF"]:
        grblat = grb_file.coords["latitude"].values
        grblon = grb_file.coords["longitude"].values
        grblev = grb_file.coords['isobaricInhPa'].values
        grblev *= 100.0

    elif datasource == "RAP":
        rhlevName = grb_file["RH_P0_L100_GLC0"].dims[0]
        print(f"rhlev varname: {rhlevName}")
        grblat = grb_file["gridlat_0"].values
        grblon = grb_file["gridlon_0"].values
        grblev = grb_file["lv_ISBL0"].values
        rhlev = grb_file[rhlevName].values

    elif datasource in ["ERAI", "ERA5"]:
        grblat = grb_file["latitude"].values.astype(float)
        grblon = grb_file["longitude"].values.astype(float)
        grblev = grb_file["level"].values.astype(float)
        if datasource == "ERA5":
            grblev *= 100.0

    elif datasource == "ERA5RDA":
        pl_dir = f"{RDADIR}/e5.oper.an.pl/{yearstr}{monthstr}"
        sf_dir = f"{RDADIR}/e5.oper.an.sfc/{yearstr}{monthstr}"
        rda_find_pattern = f"{pl_dir}/e5.oper.an.pl.128_130_t.ll025sc.{yearstr}{monthstr}{daystr}00_*.nc"

        # Use glob to find all matching files
        rda_files = glob.glob(rda_find_pattern)

        if not rda_files:
            raise FileNotFoundError(f"No files found matching pattern: {rda_find_pattern}")

        # Open all matching files as a single dataset
        rda_file = xr.open_mfdataset(rda_files, combine='by_coords')

        grblat = rda_file["latitude"].values.astype(float)
        grblon = rda_file["longitude"].values.astype(float)
        grblev = rda_file["level"].values.astype(float) * 100.0

    elif datasource == "CAM":
        grblev = np.array([
            20, 30, 50, 100, 200, 300, 500, 700, 1000, 2000, 3000, 5000, 7000, 10000,
            15000, 20000, 25000, 30000, 35000, 40000, 45000, 50000, 55000, 60000,
            65000, 70000, 75000, 80000, 82500, 85000, 87500, 90000, 92000, 93000,
            94000, 95000, 96000, 97000, 98000, 99000, 99500, 100000, 100500, 101000, 102000, 103000
        ]).astype(float)
        print(f"CAM: Using mod_remap_file: {mod_remap_file}")
        tmp = ESMF_regrid_with_weights(grb_file["PS"][0, :, :], mod_remap_file, False)
        grblat = tmp.coords['lat'].values
        grblon = tmp.coords['lon'].values

    else:
        print("Exiting when Loading lat/lon/lev")
        sys.exit(1)

    return grblat, grblon, grblev




# def load_era5rda_data(RDADIR, yearstr, monthstr, daystr, cyclestr, dycore):
#     p0 = 100000.0  # Pa
#     opt = 0  # Placeholder for options, adjust as needed
#
#     # Initialize storage for variables
#     ps = t_gfs = u_gfs = v_gfs = q_gfs = cldliq_gfs = cldice_gfs = None
#     w_gfs = z_gfs = None
#     w_is_omega = z_is_phi = False
#
#     # Define directories
#     pl_dir = f"{RDADIR}/e5.oper.an.pl/{yearstr}{monthstr}"
#     sf_dir = f"{RDADIR}/e5.oper.an.sfc/{yearstr}{monthstr}"
#
#     # Load surface pressure
#     rda_find = glob.glob(f"{sf_dir}/e5.oper.an.sfc.128_134_sp.ll025sc.{yearstr}{monthstr}0100_*.nc")
#     if not rda_find:
#         raise FileNotFoundError("No matching files found for surface pressure")
#     rda_file = xr.open_dataset(rda_find[0])
#     rda_time = rda_file["time"]
#     rda_thistime = find_closest_time(rda_time, yearstr, monthstr, daystr, cyclestr)
#     ps = rda_file["SP"].sel(time=rda_thistime).values
#
#     # Load temperature
#     rda_find = glob.glob(f"{pl_dir}/e5.oper.an.pl.128_130_t.ll025sc.{yearstr}{monthstr}{daystr}00_*.nc")
#     if not rda_find:
#         raise FileNotFoundError("No matching files found for temperature")
#     rda_file = xr.open_dataset(rda_find[0])
#     rda_time = rda_file["time"]
#     rda_thistime = find_closest_time(rda_time, yearstr, monthstr, daystr, cyclestr)
#     t_gfs = rda_file["T"].sel(time=rda_thistime).values
#
#     # Load U wind
#     rda_find = glob.glob(f"{pl_dir}/e5.oper.an.pl.128_131_u.ll025uv.{yearstr}{monthstr}{daystr}00_*.nc")
#     if not rda_find:
#         raise FileNotFoundError("No matching files found for U wind")
#     rda_file = xr.open_dataset(rda_find[0])
#     u_gfs = rda_file["U"].sel(time=rda_thistime).values
#
#     # Load V wind
#     rda_find = glob.glob(f"{pl_dir}/e5.oper.an.pl.128_132_v.ll025uv.{yearstr}{monthstr}{daystr}00_*.nc")
#     if not rda_find:
#         raise FileNotFoundError("No matching files found for V wind")
#     rda_file = xr.open_dataset(rda_find[0])
#     v_gfs = rda_file["V"].sel(time=rda_thistime).values
#
#     # Load specific humidity
#     rda_find = glob.glob(f"{pl_dir}/e5.oper.an.pl.128_133_q.ll025sc.{yearstr}{monthstr}{daystr}00_*.nc")
#     if not rda_find:
#         raise FileNotFoundError("No matching files found for specific humidity")
#     rda_file = xr.open_dataset(rda_find[0])
#     q_gfs = rda_file["Q"].sel(time=rda_thistime).values
#
#     # Load cloud liquid water content
#     rda_find = glob.glob(f"{pl_dir}/e5.oper.an.pl.128_246_clwc.ll025sc.{yearstr}{monthstr}{daystr}00_*.nc")
#     if not rda_find:
#         raise FileNotFoundError("No matching files found for cloud liquid water content")
#     rda_file = xr.open_dataset(rda_find[0])
#     cldliq_gfs = rda_file["CLWC"].sel(time=rda_thistime).values
#
#     # Load cloud ice water content
#     rda_find = glob.glob(f"{pl_dir}/e5.oper.an.pl.128_247_ciwc.ll025sc.{yearstr}{monthstr}{daystr}00_*.nc")
#     if not rda_find:
#         raise FileNotFoundError("No matching files found for cloud ice water content")
#     rda_file = xr.open_dataset(rda_find[0])
#     cldice_gfs = rda_file["CIWC"].sel(time=rda_thistime).values
#
#     # Load W and Z if dycore is MPAS
#     if dycore == "mpas":
#         # Load vertical velocity
#         rda_find = glob.glob(f"{pl_dir}/e5.oper.an.pl.128_135_w.ll025sc.{yearstr}{monthstr}{daystr}00_*.nc")
#         if not rda_find:
#             raise FileNotFoundError("No matching files found for vertical velocity")
#         rda_file = xr.open_dataset(rda_find[0])
#         w_gfs = rda_file["W"].sel(time=rda_thistime).values
#         w_is_omega = True
#
#         # Load geopotential height
#         rda_find = glob.glob(f"{pl_dir}/e5.oper.an.pl.128_129_z.ll025sc.{yearstr}{monthstr}{daystr}00_*.nc")
#         if not rda_find:
#             raise FileNotFoundError("No matching files found for geopotential height")
#         rda_file = xr.open_dataset(rda_find[0])
#         z_gfs = rda_file["Z"].sel(time=rda_thistime).values
#         z_is_phi = True
#
#     return ps, t_gfs, u_gfs, v_gfs, q_gfs, cldliq_gfs, cldice_gfs, w_gfs, z_gfs, w_is_omega, z_is_phi

def load_ERA5RDA_variable(varname, the_dir, var_code, yearstr, monthstr, daystr, cyclestr, return_coords=False):
    """Helper function to load a variable from a NetCDF file."""
    rda_find = glob.glob(f"{the_dir}/{var_code}.{yearstr}{monthstr}{daystr}00_*.nc") or \
               glob.glob(f"{the_dir}/{var_code}.{yearstr}{monthstr}0100_*.nc")
    if not rda_find:
        raise FileNotFoundError(f"No matching files found for {var_code}")
    rda_file = xr.open_dataset(rda_find[0])
    rda_time = rda_file["time"]
    rda_thistime = find_closest_time(rda_time, yearstr, monthstr, daystr, cyclestr)
    rda_data = rda_file[varname].sel(time=rda_thistime, method='nearest').values
    if return_coords:
        return rda_data, rda_file["latitude"].values.astype(float), rda_file["longitude"].values.astype(float), rda_file["level"].values.astype(float) * 100.0
    else:
        return rda_data

def load_ERA5RDA_data(RDADIR, data_filename, yearstr, monthstr, daystr, cyclestr, dycore):
    # Define directories
    pl_dir = f"{RDADIR}/e5.oper.an.pl/{yearstr}{monthstr}"
    sf_dir = f"{RDADIR}/e5.oper.an.sfc/{yearstr}{monthstr}"

#         grblat = rda_file["latitude"].values.astype(float)
#         grblon = rda_file["longitude"].values.astype(float)
#         grblev = rda_file["level"].values.astype(float) * 100.0

    # Dictionary to store the variables
    data_vars = {}

    _, data_vars['lat'], data_vars['lon'], data_vars['lev'] = load_ERA5RDA_variable('T', pl_dir, "e5.oper.an.pl.128_130_t.ll025sc", yearstr, monthstr, daystr, cyclestr, return_coords=True)

    # Load required variables
    data_vars['ps'] = load_ERA5RDA_variable('SP', sf_dir, "e5.oper.an.sfc.128_134_sp.ll025sc", yearstr, monthstr, daystr, cyclestr)
    data_vars['t'] = load_ERA5RDA_variable('T', pl_dir, "e5.oper.an.pl.128_130_t.ll025sc", yearstr, monthstr, daystr, cyclestr)
    data_vars['u'] = load_ERA5RDA_variable('U', pl_dir, "e5.oper.an.pl.128_131_u.ll025uv", yearstr, monthstr, daystr, cyclestr)
    data_vars['v'] = load_ERA5RDA_variable('V', pl_dir, "e5.oper.an.pl.128_132_v.ll025uv", yearstr, monthstr, daystr, cyclestr)
    data_vars['q'] = load_ERA5RDA_variable('Q', pl_dir, "e5.oper.an.pl.128_133_q.ll025sc", yearstr, monthstr, daystr, cyclestr)
    data_vars['cldliq'] = load_ERA5RDA_variable('CLWC', pl_dir, "e5.oper.an.pl.128_246_clwc.ll025sc", yearstr, monthstr, daystr, cyclestr)
    data_vars['cldice'] = load_ERA5RDA_variable('CIWC', pl_dir, "e5.oper.an.pl.128_247_ciwc.ll025sc", yearstr, monthstr, daystr, cyclestr)

    if dycore == 'mpas':
        data_vars['w'] = load_ERA5RDA_variable('W', pl_dir, "e5.oper.an.pl.128_135_w.ll025sc", yearstr, monthstr, daystr, cyclestr)
        data_vars['w_is_omega'] = True
        data_vars['z'] = load_ERA5RDA_variable('Z', pl_dir, "e5.oper.an.pl.128_129_z.ll025sc", yearstr, monthstr, daystr, cyclestr)
        data_vars['z_is_phi'] = True

    data_vars['ts'] = load_ERA5RDA_variable('VAR_2T', sf_dir, "e5.oper.an.sfc.128_167_2t.ll025sc", yearstr, monthstr, daystr, cyclestr)

    ds = load_ERA5_file(data_filename)
    data_vars['phis'] = ds["Z"].isel(time=0).values
    ds.close

    return data_vars

def load_CFSR_file(data_filename,TYPE_OF_LEVEL,VAR_SHORTNAME):
    return xr.open_dataset(
        data_filename,
        engine="cfgrib",
        backend_kwargs={
            'filter_by_keys': {
                'typeOfLevel': TYPE_OF_LEVEL,
                'shortName':  VAR_SHORTNAME
            }
        }
    )

def load_ERA5_file(data_filename):
    # Load the data file using xarray
    return xr.open_dataset(
        data_filename
    )

def load_CAM_file(data_filename):
    # Load the data file using xarray
    return xr.open_dataset(
        data_filename
    )

def load_CFSR_variable(grb_file, varname):
    """Helper function to load a variable from a CFSR NetCDF file."""
    if varname in grb_file.variables:
        return grb_file[varname].values
    else:
        raise KeyError(f"Variable {varname} not found in the GRIB file.")

def get_CFSR_levels_for_var(grb_file_name, variable):
    grb_file = load_CFSR_file(grb_file_name, 'isobaricInhPa', variable)
    levels = grb_file.coords['isobaricInhPa'].values * 100.0
    grb_file.close()
    return levels

def get_CFSR_coords_for_var(grb_file_name, variable, coord, filter_arg='isobaricInhPa'):
    grb_file = load_CFSR_file(grb_file_name, filter_arg, variable)
    coords = grb_file.coords[coord].values
    grb_file.close()
    return coords

def load_and_extract_CFSR_variable(grb_file_name, level_type, variable):
    """
    Load and extract a variable from a CFSR GRIB file.

    Parameters:
    -----------
    grb_file_name : str
        The name of the GRIB file.

    level_type : str
        The type of level to filter by (e.g., 'surface', 'isobaricInhPa').

    variable : str
        The variable name to load from the GRIB file (e.g., 'sp', 't', 'u', 'v').

    Returns:
    --------
    np.ndarray
        The extracted data variable.
    """
    grb_file = load_CFSR_file(grb_file_name, level_type, variable)
    data = load_CFSR_variable(grb_file, variable)
    grb_file.close()
    return data

def load_CFSR_data(grb_file_name, dycore):

    # Dictionary to store the variables
    data_vars = {}

    # Use temperature to get base three coordinates
    data_vars['lat'] = get_CFSR_coords_for_var(grb_file_name,'t','latitude')
    data_vars['lon'] = get_CFSR_coords_for_var(grb_file_name,'t','longitude')
    data_vars['lev'] = get_CFSR_levels_for_var(grb_file_name, 't')

    data_vars['cldlev'] = get_CFSR_levels_for_var(grb_file_name, 'clwmr')
    data_vars['rhlev'] = get_CFSR_levels_for_var(grb_file_name, 'r')
    data_vars['windlev'] = get_CFSR_levels_for_var(grb_file_name, 'u')

    data_vars['ps'] = load_and_extract_CFSR_variable(grb_file_name, 'surface', 'sp')
    data_vars['t'] = load_and_extract_CFSR_variable(grb_file_name, 'isobaricInhPa', 't')
    data_vars['u'] = load_and_extract_CFSR_variable(grb_file_name, 'isobaricInhPa', 'u')
    data_vars['v'] = load_and_extract_CFSR_variable(grb_file_name, 'isobaricInhPa', 'v')

    if dycore == "mpas":
        data_vars['w'] = load_and_extract_CFSR_variable(grb_file_name, 'isobaricInhPa', 'w')
        data_vars['w_is_omega'] = True
        data_vars['z'] = load_and_extract_CFSR_variable(grb_file_name, 'isobaricInhPa', 'gh')
        data_vars['z_is_phi'] = False

    if data_vars['windlev'].shape != data_vars['lev'].shape:
        print(f"Interpolating u + v wind from {data_vars['windlev'].shape} to {data_vars['lev'].shape} levels")
        data_vars['u'] = int2p_n(windlev, data_vars['u'], grblev, linlog=2, dim=0)
        data_vars['v'] = int2p_n(windlev, data_vars['v'], grblev, linlog=2, dim=0)

    # Load and interpolate additional variables
    data_vars['rh_native'] = load_and_extract_CFSR_variable(grb_file_name, 'isobaricInhPa', 'r')
    data_vars['cldmix_native'] = load_and_extract_CFSR_variable(grb_file_name, 'isobaricInhPa', 'clwmr')

    # Interpolate
    data_vars['rh'] = int2p_n(data_vars['rhlev'],data_vars['rh_native'], data_vars['lev'], linlog=2, dim=0)
    data_vars['cldmix'] = int2p_n(data_vars['cldlev'],data_vars['cldmix_native'], data_vars['lev'], linlog=2, dim=0)

    # Calculate specific humidity from RH
    # Load grblev (Pa), use numpy broadcasting to go from 1-D to 3-D, multiple by 0.01 to convert to mb
    data_vars['q'] = mixhum_ptrh((data_vars['lev'][:, np.newaxis, np.newaxis] * 0.01), data_vars['t'], data_vars['rh'])

    # Sort bad values
    data_vars['cldmix'] = np.where(np.isnan(data_vars['cldmix']), 0, data_vars['cldmix'])
    data_vars['cldmix'] = clip_and_count(data_vars['cldmix'],max_thresh=0.01,var_name="cldmix")

    # Separate cloud ice and water
    data_vars['cldice'] = np.where(data_vars['t'] > t_freeze_K, 0, data_vars['cldmix'])
    data_vars['cldliq'] = np.where(data_vars['t'] > t_freeze_K, data_vars['cldmix'], 0)

    data_vars['phis']  = load_and_extract_CFSR_variable(grb_file_name, 'surface', 'orog')
    data_vars['phis'] = data_vars['phis'] * grav
    data_vars['ts'] = load_and_extract_CFSR_variable(grb_file_name, 'sigma', 't')

    # Cleanup
    del data_vars['cldmix']
    del data_vars['cldlev']
    del data_vars['rhlev']
    del data_vars['windlev']
    del data_vars['rh']
    del data_vars['cldmix_native']
    del data_vars['rh_native']

    return data_vars






def load_cam_data(grb_file_name, YYYYMMDDHH, mod_in_topo, mod_remap_file, dycore, debug=False):

    # Initialize variables
    w_smooth_iter = 1
    data_vars = {}

    grb_file = load_CAM_file(grb_file_name)

    data_vars['lev'] = np.array([
        20, 30, 50, 100, 200, 300, 500, 700, 1000, 2000, 3000, 5000, 7000, 10000,
        15000, 20000, 25000, 30000, 35000, 40000, 45000, 50000, 55000, 60000,
        65000, 70000, 75000, 80000, 82500, 85000, 87500, 90000, 92000, 93000, 94000,
        95000, 96000, 97000, 98000, 99000,
        99500, 100000, 100500, 101000, 102000, 103000
    ], dtype=np.float32)

    yearstr, monthstr, daystr, cyclestr = split_by_lengths(str(YYYYMMDDHH), dtime_map)

    # Access the time coordinate as a NumPy array
    cam_time = grb_file["time"]

#     # Convert the strings to an ISO 8601 datetime format
#     datetime_str = f"{yearstr}-{monthstr.zfill(2)}-{daystr.zfill(2)}T00:00:00"
#     datetime_base = np.datetime64(datetime_str)
#
#     # Add the number of seconds from cyclestr to get the exact datetime
#     datetime_target = datetime_base + np.timedelta64(int(cyclestr), 's')
#     print(f"Target datetime: {datetime_target}")
#
#     # Find the index of the closest time using NumPy operations
#     cam_thistime_ix = np.argmin(np.abs(cam_time - datetime_target))

    cam_thistime_ix = find_closest_time(cam_time.values, yearstr, monthstr, daystr, cyclestr, return_isel=True)
    print(f"Closest time: {cam_thistime_ix}")
    cam_thistime = find_closest_time(cam_time, yearstr, monthstr, daystr, cyclestr)
    print(f"Closest time: {cam_thistime}")

    # Select the closest time
    data_vars['ps'] = grb_file["PS"].sel(time=cam_thistime, method='nearest').values
    data_vars['t'] = grb_file["T"].sel(time=cam_thistime, method='nearest').values
    data_vars['u'] = grb_file["U"].sel(time=cam_thistime, method='nearest').values
    data_vars['v'] = grb_file["V"].sel(time=cam_thistime, method='nearest').values
    data_vars['q'] = grb_file["Q"].sel(time=cam_thistime, method='nearest').values

    print(f"CAM: Using mod_in_topo: {mod_in_topo}")

    # Load topo data
    mod_in_topo_f = xr.open_dataset(mod_in_topo)
    data_vars['phis'] = mod_in_topo_f["PHIS"].values
    data_vars['ts'] = data_vars['t'][-1,:]

    # Horizontal remaps
    data_vars = horizremap.remap_all(data_vars, mod_remap_file, dycore=dycore)

    print_elapsed_time(msg="Horizontal remap")

    # Get relevant vertical coefficients from CAM file
    hyam = grb_file["hyam"].values
    hybm = grb_file["hybm"].values
    hyai = grb_file["hyai"].values
    hybi = grb_file["hybi"].values

    print_min_max_dict(data_vars)

    if dycore == 'mpas':
        # Calculate omega

        print_elapsed_time(msg="Start MPAS")

        data_vars['w_is_omega'] = True
        if grb_file["lat"][0] > grb_file["lat"][1]:
            print("flipping omega since calc needs to be S->N")
            data_vars['w'] = omega_ccm_driver(p0, data_vars['ps'][::-1, :], data_vars['u'][:, ::-1, :], data_vars['v'][:, ::-1, :], data_vars['lat'][::-1], data_vars['lon'], hyam, hybm, hyai, hybi)
            data_vars['w'] = data_vars['w'][:, ::-1, :]  # flip back
        else:
            data_vars['w'] = omega_ccm_driver(p0, data_vars['ps'], data_vars['u'], data_vars['v'], data_vars['lat'], data_vars['lon'], hyam, hybm, hyai, hybi)
        # Remove omega poleward of OMEGA_LAT_THRESH to deal with singularity
        data_vars['w'] = np.where(np.abs(data_vars['lat'])[np.newaxis, :, np.newaxis] > OMEGA_LAT_THRESH, 0.0, data_vars['w'])

        print_elapsed_time(msg="Done omega")

        #data_vars['div'] = ddvfidf_wrapper(data_vars['u'], data_vars['v'], data_vars['lat'], data_vars['lon'], 3)

        tkv = data_vars['t'] * (1. + 0.61 * data_vars['q'])
        data_vars['z'] = cz2ccm(data_vars['ps'], data_vars['phis'], tkv, p0, hyam[::-1], hybm[::-1], hyai[::-1], hybi[::-1])
        data_vars['z_is_phi'] = False
        #data_vars['dpsl'], data_vars['dpsm'] = calculate_gradients(data_vars['ps'], data_vars['lat'], data_vars['lon'])
        #data_vars['div'], data_vars['vort'] = calculate_div_vort(data_vars['lat'], data_vars['lon'], data_vars['u'], data_vars['v'])
        #data_vars['pdel'] = dpres_hybrid_ccm(data_vars['ps'], p0, hyai, hybi)
        #data_vars['pmid'] = pres_hybrid_ccm(data_vars['ps'], p0, hyam, hybm)

        print_elapsed_time(msg="Done Z")

        # Use the print_debug_file function to create and save the xarray.Dataset
        print_debug_file(
              "py_cam_raw.nc",
              ps_cam=(["lat", "lon"], data_vars['ps']),
              phis_cam=(["lat", "lon"], data_vars['phis']),
              ts_cam=(["lat", "lon"], data_vars['ts']),
              t_cam=(["lev", "lat", "lon"], data_vars['t']),
              u_cam=(["lev", "lat", "lon"], data_vars['u']),
              v_cam=(["lev", "lat", "lon"], data_vars['v']),
              q_cam=(["lev", "lat", "lon"], data_vars['q']),
              #tkv_cam=(["lev", "lat", "lon"], data_vars['tkv']),
              z_cam=(["lev", "lat", "lon"], data_vars['z']),
              #dpsl_cam=(["lat", "lon"], data_vars['dpsl']),
              #dpsm_cam=(["lat", "lon"], data_vars['dpsm']),
              #div_cam=(["lev", "lat", "lon"], data_vars['div']),
              #pdel_cam=(["lev", "lat", "lon"], data_vars['pdel']),
              #pmid_cam=(["lev", "lat", "lon"], data_vars['pmid']),
              w_cam=(["lev", "lat", "lon"], data_vars['w']),
              lat=(["lat"], data_vars['lat']),
              lon=(["lon"], data_vars['lon'])
        )

    print_min_max_dict(data_vars)

    # Vertically interpolate the CAM hybrid levels to constant pressure surfaces
    data_vars = interp_hybrid_to_pressure_wrapper(
        data_vars=data_vars,
        ps=data_vars['ps'],
        hyam=hyam,
        hybm=hybm,
        new_levels=data_vars['lev']
        )

    print_elapsed_time(msg="Done hyb-pres in CAM")

    data_vars['cldice'] = np.zeros_like(data_vars['t'])
    data_vars['cldliq'] = np.zeros_like(data_vars['t'])

    return data_vars

def interp_hybrid_to_pressure_wrapper(data_vars, ps, hyam, hybm, new_levels, lev_dim=0, method='log', extrapolate=True):

    allowable_interp_vars = ['t', 'u', 'v', 'q', 'cldice', 'cldliq', 'z', 'theta', 'rho', 'w']

    for var_name, data in data_vars.items():
        if isinstance(data, np.ndarray) and data.ndim == 3:
            variable_type = 'temperature' if var_name == 't' else 'geopotential' if var_name == 'z' else 'other'
            print(f"interp_hybrid_to_pressure for variable {var_name}, using {variable_type}")
            data_vars[var_name] = interp_hybrid_to_pressure(
                data=data,
                ps=ps,
                hyam=hyam,
                hybm=hybm,
                new_levels=new_levels,
                lev_dim=lev_dim,
                method=method,
                extrapolate=extrapolate,
                variable=variable_type,
                t_bot=data_vars['ts'],
                phi_sfc=data_vars['phis']
            )
    return data_vars

def interpolate_to_uniform_levels(data_var):
    ##### int2p_n

    # Placeholder for actual interpolation logic
    # Implement int2p_n or equivalent interpolation in Python
    return data_var  # Modify with actual interpolation

def calculate_specific_humidity(t_gfs, rh_gfs):
    # Placeholder for the actual calculation of specific humidity
    # Implement mixhum_ptrh or equivalent in Python
    return rh_gfs  # Modify with actual calculation

def mixhum_ptrh(p, tk, rh, iswit=2):
    """
    Computes the specific humidity or mixing ratio from pressure, temperature, and relative humidity.

    Parameters:
    -----------
    p : float or numpy.ndarray
        Atmospheric pressure (hPa or mb).
    tk : float or numpy.ndarray
        Temperature in Kelvin.
    rh : float or numpy.ndarray
        Relative humidity in percent.
    iswit : int
        - If iswit=1, the output will be the mixing ratio.
        - If iswit=2, the output will be the specific humidity.
        - If iswit is negative, the units will be kg/kg; otherwise, g/kg.

    Returns:
    --------
    qw : float or numpy.ndarray
        Mixing ratio or specific humidity in the units specified by iswit.
    """

    # Constants
    T0 = 273.15  # Reference temperature in Kelvin
    EP = 0.622  # Ratio of molecular weights of water vapor to dry air
    ONEMEP = 1.0 - EP  # 1 - EP
    ES0 = 6.11  # Reference saturation vapor pressure in hPa or mb
    A = 17.269  # Coefficient for saturation vapor pressure calculation
    B = 35.86   # Coefficient for saturation vapor pressure calculation

    # Calculate the saturation vapor pressure (EST)
    est = ES0 * np.exp((A * (tk - T0)) / (tk - B))

    # Calculate the saturation mixing ratio (QST)
    qst = (EP * est) / (p - ONEMEP * est)

    # Calculate mixing ratio or specific humidity
    qw = qst * (rh * 0.01)

    # Convert to specific humidity if iswit=2
    if abs(iswit) == 2:
        qw = qw / (1.0 + qw)

    # Convert to kg/kg if iswit is negative
    if iswit < 0:
        qw = qw * 1000.0

    return qw







def interpolate_to_levels(ppin, xxin, ppout, linlog=1, xmsg=np.nan):
    """
    Interpolate or extrapolate data from one set of pressure levels to another.

    Parameters:
    -----------
    ppin : numpy.ndarray
        Input pressure levels (1D array).
    xxin : numpy.ndarray
        Input variable values at ppin levels (1D array).
    ppout : numpy.ndarray
        Output pressure levels (1D array).
    linlog : int, optional
        Flag indicating interpolation method:
        - 1: linear interpolation (default)
        - 0: logarithmic interpolation
        - Negative: extrapolation allowed using the nearest valid slope.
    xmsg : float, optional
        Missing data marker. Default is np.nan.

    Returns:
    --------
    xxout : numpy.ndarray
        Interpolated variable values at ppout levels.
    """

    # Initialize output array with missing values
    xxout = np.full_like(ppout, xmsg)

    # Error check: make sure we have enough points
    if len(ppin) < 2 or len(ppout) < 1:
        return xxout  # Return all missing values

    # Determine if input and output pressures need to be reordered
    npin = len(ppin)
    nout = len(ppout)

    ppin_reordered = ppin[::-1] if ppin[0] < ppin[1] else ppin
    xxin_reordered = xxin[::-1] if ppin[0] < ppin[1] else xxin
    ppout_reordered = ppout[::-1] if ppout[0] < ppout[1] else ppout

    # Filter out missing data in xxin
    valid_mask = (xxin_reordered != xmsg) & (ppin_reordered != xmsg)
    pin_valid = ppin_reordered[valid_mask]
    xin_valid = xxin_reordered[valid_mask]

    if len(pin_valid) < 2:
        return xxout  # Not enough valid data points

    # Perform interpolation
    if linlog == 1:
        # Linear interpolation
        xxout = np.interp(ppout_reordered, pin_valid, xin_valid)
    else:
        # Logarithmic interpolation
        xxout = np.interp(np.log(ppout_reordered), np.log(pin_valid), xin_valid)

    # Extrapolation if linlog < 0
    if linlog < 0:
        if ppout_reordered[0] > pin_valid[0]:  # Extrapolate above range
            slope = (xin_valid[1] - xin_valid[0]) / (pin_valid[1] - pin_valid[0])
            xxout[ppout_reordered > pin_valid[0]] = xin_valid[0] + slope * (ppout_reordered[ppout_reordered > pin_valid[0]] - pin_valid[0])
        if ppout_reordered[-1] < pin_valid[-1]:  # Extrapolate below range
            slope = (xin_valid[-1] - xin_valid[-2]) / (pin_valid[-1] - pin_valid[-2])
            xxout[ppout_reordered < pin_valid[-1]] = xin_valid[-1] + slope * (ppout_reordered[ppout_reordered < pin_valid[-1]] - pin_valid[-1])

    # Reorder output back to original if necessary
    if ppout[0] < ppout[1]:
        xxout = xxout[::-1]

    return xxout



def find_closest_time(times, yearstr, monthstr, daystr, cyclestr, return_isel=False):
    if isinstance(times, np.ndarray):
        # Convert the strings to an ISO 8601 datetime format
        datetime_str = f"{yearstr}-{monthstr.zfill(2)}-{daystr.zfill(2)}T00:00:00"
        datetime_base = np.datetime64(datetime_str)

        # Add the number of seconds from cyclestr to get the exact datetime
        datetime_target = datetime_base + np.timedelta64(int(cyclestr), 's')
        print(f"Target datetime: {datetime_target}")

        closest_index = np.argmin(np.abs(times - datetime_target))
        if return_isel:
            closest_time = closest_index
        else:
            closest_time = times[closest_index]
    else:
        target_time = np.datetime64(f"{yearstr}-{monthstr}-{daystr}T{cyclestr[:2]}:00")
        closest_time = min(times.values, key=lambda x: abs(x - target_time))

    return closest_time





def load_cam_levels(PATHTOHERE, numlevels, load_xarray=False):
    """
    Loads the CAM levels data from a specified template file.

    Parameters:
    PATHTOHERE: Path to the directory containing the template file
    numlevels: Number of vertical levels to load

    Returns:
    hya, hyb, hyai, hybi, lev, ilev: The hybrid coefficients and levels
    """
    # Load the CAM levels file
    template_path = f"{PATHTOHERE}/templates/L{numlevels}template.nc"
    fC = xr.open_dataset(template_path)

    # Extract the hybrid coefficients and levels
    if load_xarray:
        hya = fC["hyam"]
        hyb = fC["hybm"]
        hyai = fC["hyai"]
        hybi = fC["hybi"]
        lev = fC["lev"]
        ilev = fC["ilev"]
    else:
        hya = fC["hyam"].values
        hyb = fC["hybm"].values
        hyai = fC["hyai"].values
        hybi = fC["hybi"].values
        lev = fC["lev"].values
        ilev = fC["ilev"].values

    print("---------------------------------------------------------")
    print("Loading CAM levels")
    print(f"Loading {numlevels} level data")

    return hya, hyb, hyai, hybi, lev, ilev





def prcwater_dp(Q, DP, QMSG=np.nan, DPMSG=np.nan):
    """
    Calculate precipitable water given specific humidity and layer thickness.

    Parameters:
    - Q: Specific humidity array [kg/kg; dimensionless].
    - DP: Layer thickness array [Pa].
    - QMSG: Missing value indicator for Q (optional, default: NaN).
    - DPMSG: Missing value indicator for DP (optional, default: NaN).

    Returns:
    - prcwat: Precipitable water [kg/m2].
    """

    # Initialize precipitable water
    prcwat = 0.0

    # Count valid layers
    valid_layers = 0

    # Loop over each layer
    for q, dp in zip(Q, DP):
        if not np.isnan(q) and not np.isnan(dp):
            if q != QMSG and dp != DPMSG:
                valid_layers += 1
                prcwat += q * abs(dp)

    # Final precipitable water calculation
    if valid_layers > 0:
        prcwat /= grav
    else:
        prcwat = QMSG

    return prcwat


def ps_wet_to_dry_conversion(ps_fv, q_fv, hyai, hybi, p0, verbose=False):
    """
    Converts wet surface pressure to dry surface pressure by subtracting
    the total column precipitable water.

    Parameters:
    - ps_fv: Surface pressure array (1D).
    - q_fv: Specific humidity array (2D, with shape [levels, columns]).
    - hyai: Hybrid A interface coefficients (1D).
    - hybi: Hybrid B interface coefficients (1D).
    - p0: Reference pressure (scalar).
    - verbose: If True, print intermediate results for every 10000th column.

    Returns:
    - ps_fv: Updated dry surface pressure array.
    """

    ncol = ps_fv.shape[0]
    pw_fv = np.zeros_like(ps_fv)

    for kk in range(ncol):
        pi_orig = hyai * p0 + hybi * ps_fv[kk]
        nlevp1 = pi_orig.size
        dp = pi_orig[1:nlevp1] - pi_orig[0:nlevp1-1]
        pw_fv[kk] = prcwater_dp(q_fv[:, kk], dp)  # prcwater_dp is assumed to be a pre-defined function
        ps_fv[kk] = ps_fv[kk] - pw_fv[kk] * grav

        if verbose and kk % 10000 == 0:
            print(dp)
            print(f"Correcting PS: {kk} of {ncol-1} from {ps_fv[kk] + pw_fv[kk] * 9.81} to {ps_fv[kk]} since TPW: {pw_fv[kk]}")

    print("Done!")

    return ps_fv, pw_fv

def numpy_to_dataarray(numpy_array, dims, coords=None, attrs=None, name=None):
    """
    Converts a numpy array to an xarray DataArray.

    Parameters:
    - numpy_array: The numpy array to convert.
    - dims: A list of dimension names for the DataArray.
    - coords: Optional dictionary of coordinate arrays.
    - attrs: Optional dictionary of attributes to add to the DataArray.
    - name: Optional name for the DataArray.

    Returns:
    - data_array: The resulting xarray DataArray.
    """
    data_array = xr.DataArray(data=numpy_array, dims=dims, coords=coords, attrs=attrs, name=name)
    return data_array

def add_time_define_precision(var_in, precision, isncol, lat_dim="lat", lon_dim="lon"):
    """
    Adds a time dimension and converts the precision of the input variable.

    Parameters:
    var_in (xarray.DataArray or numpy.ndarray): Input variable.
    precision (str): Desired precision ("float", "single", "double").
    isncol (bool): True if the variable is ncol; False otherwise.
    lat_dim (str): Name of the latitude dimension (default is "lat").
    lon_dim (str): Name of the longitude dimension (default is "lon").

    Returns:
    xarray.DataArray: Output variable with time dimension and specified precision.
    """
    var_dims = var_in.shape

    if isncol:
        if len(var_dims) == 1:  # ncol -> time, ncol
            print(f"Converting a ncol -> ncol {precision} with time attached!")
            ncol = var_dims[0]
            if precision in ["float", "single"]:
                var_out = xr.DataArray(np.zeros((1, ncol), dtype=np.float32), dims=["time", "ncol"])
            elif precision == "double":
                var_out = xr.DataArray(np.zeros((1, ncol), dtype=np.float64), dims=["time", "ncol"])
            else:
                raise ValueError("Invalid precision specified")
            var_out[0, :] = var_in.astype(var_out.dtype)

        elif len(var_dims) == 2:  # lev, ncol -> time, lev, ncol
            print(f"Converting a lev, ncol -> lev, ncol {precision} with time attached!")
            nlev, ncol = var_dims
            if precision in ["float", "single"]:
                var_out = xr.DataArray(np.zeros((1, nlev, ncol), dtype=np.float32), dims=["time", "lev", "ncol"])
            elif precision == "double":
                var_out = xr.DataArray(np.zeros((1, nlev, ncol), dtype=np.float64), dims=["time", "lev", "ncol"])
            else:
                raise ValueError("Invalid precision specified")
            var_out[0, :, :] = var_in.astype(var_out.dtype)

    else:  # Not ncol
        if len(var_dims) == 2:  # lat_dim, lon_dim -> time, lat_dim, lon_dim
            print(f"Converting a {lat_dim}, {lon_dim} -> {lat_dim}, {lon_dim} {precision} with time attached!")
            nlat, nlon = var_dims
            if precision in ["float", "single"]:
                var_out = xr.DataArray(np.zeros((1, nlat, nlon), dtype=np.float32), dims=["time", lat_dim, lon_dim])
            elif precision == "double":
                var_out = xr.DataArray(np.zeros((1, nlat, nlon), dtype=np.float64), dims=["time", lat_dim, lon_dim])
            else:
                raise ValueError("Invalid precision specified")
            var_out[0, :, :] = var_in.astype(var_out.dtype)

        elif len(var_dims) == 3:  # lev, lat_dim, lon_dim -> time, lev, lat_dim, lon_dim
            print(f"Converting a lev, {lat_dim}, {lon_dim} -> lev, {lat_dim}, {lon_dim} {precision} with time attached!")
            nlev, nlat, nlon = var_dims
            if precision in ["float", "single"]:
                var_out = xr.DataArray(np.zeros((1, nlev, nlat, nlon), dtype=np.float32), dims=["time", "lev", lat_dim, lon_dim])
            elif precision == "double":
                var_out = xr.DataArray(np.zeros((1, nlev, nlat, nlon), dtype=np.float64), dims=["time", "lev", lat_dim, lon_dim])
            else:
                raise ValueError("Invalid precision specified")
            var_out[0, :, :, :] = var_in.astype(var_out.dtype)

    # Copy any metadata over to the output variable
    var_out.attrs = var_in.attrs

    return var_out


def calculate_rho_gfs(pres_gfs, q_gfs, t_gfs, rho_d_algo=1):
    """
    Calculates the density of air (rho) using either a simple or advanced method.

    Parameters:
    -----------
    pres_gfs : numpy.ndarray
        Pressure field, typically with shape (nlev, ncol).
    q_gfs : numpy.ndarray
        Specific humidity field, typically with shape (nlev, ncol).
    t_gfs : numpy.ndarray
        Temperature field, typically with shape (nlev, ncol).
    rho_d_algo : int, optional
        Algorithm choice for calculating dry air density.
        If 1, uses a simplified formula; otherwise, uses a more complex CAM/MPAS formula.
        Default is 1.

    Returns:
    --------
    rho_gfs : numpy.ndarray
        Density of air (rho) with the same shape as pres_gfs.
    """
    if rho_d_algo == 1:
        presdry_gfs = pres_gfs / (1. + q_gfs)
        rho_gfs = presdry_gfs / (Rd * t_gfs)
    else:
        # rho_d calculation from internal CAM/MPAS code")
        rho_gfs = pres_gfs / (Rd * t_gfs * (1. + Rv_over_Rd * q_gfs))

    return rho_gfs



def pot_temp(p, t):
    """
    Compute potential temperature.

    Parameters:
    -----------
    p : numpy.ndarray
        Array containing pressure levels (Pa).
    t : numpy.ndarray
        Array containing temperatures (K).

    Returns:
    --------
    theta : numpy.ndarray
        A multi-dimensional array of the same size and shape as t, containing potential temperature values.
    """

    # Calculate potential temperature using the formula theta = t * (p0 / p) ** kappa
    theta = t * (p0 / p) ** kappa

    return theta


def omega_to_w(omega, p, t):
    """
    Converts OMEGA (Pa/s) to W (m/s) using a first-order approximation.

    Parameters:
    -----------
    omega : numpy.ndarray
        Vertical velocity in pressure coordinates (Pa/s).
    p : numpy.ndarray
        Pressure levels (Pa).
    t : numpy.ndarray
        Temperature (K).

    Returns:
    --------
    w : numpy.ndarray
        Vertical velocity in height coordinates (m/s).
    """

    # Ensure omega, p, and t have the same shape
    if not (omega.shape == p.shape == t.shape):
        raise ValueError("omega, p, and t must have the same shape")

    # Calculate density: rho = p / (RGAS * t)
    rho = p / (Rd * t)

    # Convert omega to w using the first-order approximation: w = -omega / (rho * GRAV)
    w = -omega / (rho * grav)

    return w


def clip_and_count(arr, min_thresh=None, max_thresh=None, var_name="Variable"):
    """
    Clips the values in the array to the specified minimum and maximum thresholds,
    and returns the clipped array along with the number of adjustments made.

    Parameters:
    - arr: numpy.ndarray
        The array to be clipped.
    - min_thresh: float or None
        The minimum threshold value. Values below this will be set to min_thresh.
    - max_thresh: float or None
        The maximum threshold value. Values above this will be set to max_thresh.
    - var_name: str
        The name of the variable for diagnostic output.

    Returns:
    - clipped_arr: numpy.ndarray
        The clipped array.
    """

    # Apply the maximum threshold if specified
    if max_thresh is not None:
        max_adjustments = np.sum(arr > max_thresh)
        arr = np.where(arr > max_thresh, max_thresh, arr)
        print(f"Number of {var_name} adjustments for values above {max_thresh}: {max_adjustments}")

    # Apply the minimum threshold if specified
    if min_thresh is not None:
        min_adjustments = np.sum(arr < min_thresh)
        arr = np.where(arr < min_thresh, min_thresh, arr)
        print(f"Number of {var_name} adjustments for values below {min_thresh}: {min_adjustments}")

    if np.any(np.isnan(arr)):
        print(f"WARNING: {var_name} is missing data...")
        print("This is expected for a regional dataset, but not desirable for a global one")

    return arr





def create_cf_time(year, month, day, hour, base_time="1850-01-01 00:00:00", time_units="days since", calendar="standard"):
    """
    Create a CF-compliant time array from year, month, day, and hour.

    Parameters:
    - year, month, day, hour: Arrays or scalars of year, month, day, and hour
    - base_time: The reference time in the format "YYYY-MM-DD HH:MM:SS"
    - time_units: The time units for the output (e.g., "days since" or "hours since")

    Returns:
    - time: A numpy array of CF-compliant time values
    - time_attrs: A dictionary with attributes to assign to the time variable in xarray
    """

    base_year, base_month, base_day, base_hour = [int(val) for val in base_time.split()[0].split('-')] + [int(base_time.split()[1].split(':')[0])]

    base_time_cftime = cftime.DatetimeNoLeap(base_year, base_month, base_day, base_hour)

    nptime = cftime.date2num(cftime.DatetimeNoLeap(year, month, day, hour), units=f"{time_units} {base_time}")

    nptime = np.array([nptime], dtype=np.float64)

    time_attrs = {
        "units": f"{time_units} {base_time}",
        "calendar": f"{calendar}"
    }

    return nptime, time_attrs

def print_all_variables_info(name_width=22, type_width=40, shape_width=18, value_width=50, size_width=15):
    # Get the current frame
    frame = sys._getframe(1)

    # Get all local variables in the current frame
    local_vars = frame.f_locals

    print(f"{'Variable Name':<{name_width}} {'Type':<{type_width}} {'Shape':<{shape_width}} {'Value':<{value_width}} {'Size (bytes)':<{size_width}}")
    print("=" * (name_width + type_width + shape_width + value_width + size_width))

    for var_name, var_value in local_vars.items():
        var_type = str(type(var_value))[:type_width]
        var_shape = ''
        var_size = sys.getsizeof(var_value)

        # Determine the shape and size for numpy arrays and xarray DataArrays
        if isinstance(var_value, (np.ndarray, xr.DataArray)):
            var_shape = str(var_value.shape)[:shape_width]
            var_size = var_value.nbytes if isinstance(var_value, np.ndarray) else var_value.nbytes
        elif isinstance(var_value, xr.Dataset):
            var_shape = str({dim: size for dim, size in var_value.sizes.items()})[:shape_width]
            var_size = sum([var.nbytes for var in var_value.data_vars.values()])

        # Truncate the value for xarray.Dataset to avoid multi-line printing
        if isinstance(var_value, xr.Dataset):
            display_value = "<xarray.Dataset>"[:value_width]
        elif isinstance(var_value, xr.DataArray):
            flattened_value = var_value.values.flatten()
            display_value = f"{flattened_value[0]}..." if flattened_value.size > 0 else "Empty"
        elif isinstance(var_value, np.ndarray):
            display_value = f"{var_value.flatten()[0]}..." if var_value.size > 0 else "Empty"
        else:
            display_value = str(var_value)[:value_width]

        print(f"{var_name:<{name_width}} {var_type:<{type_width}} {var_shape:<{shape_width}} {display_value:<{value_width}} {var_size:<{size_width}}")


def print_debug_file(output_filename, **kwargs):
    """
    Create an xarray Dataset from the given variables and save it to a NetCDF file.

    Parameters:
    - output_filename: Name of the output NetCDF file
    - **kwargs: Arbitrary number of variables to include in the Dataset.
                Each key should be the variable name.
                Each value should be a tuple where the first element is a list of dimension names,
                and the second element is the variable data (numpy array).
    """
    # Prepare the data dictionary for the Dataset
    data_dict = {}
    for var_name, (dims, var_data) in kwargs.items():
        if isinstance(var_data, np.ndarray):
            # Add the variable to the Dataset dictionary
            data_dict[var_name] = (dims, var_data.astype(np.float32))
        else:
            raise ValueError(f"Unsupported data type for variable '{var_name}': {type(var_data)}")

    # Create an xarray.Dataset
    ds = xr.Dataset(data_dict)

    # Save the dataset to a NetCDF file
    ds.to_netcdf(output_filename)








def latRegWgt(lat, nType="float", opt=0):
    """
    Generates [sin(lat+dlat/2) - sin(lat-dlat/2)] weights for equally spaced (regular) global grids that will sum to 2.0.

    Parameters:
    - lat: 1D array-like of latitudes of the global grid (in degrees).
    - nType: The type of variable to be returned ("float" or "double").
    - opt: Not used. Set to 0.

    Returns:
    - weights: A 1D array of weights of the same size as lat.
    """

    # Convert lat to radians
    lat_rad = np.radians(lat)
    dlat = np.abs(np.diff(lat_rad).mean())  # Assume uniform spacing and get the mean difference

    # Calculate weights
    weights = np.sin(lat_rad + dlat / 2) - np.sin(lat_rad - dlat / 2)

    # Ensure the weights sum to 2.0
    weights *= 2.0 / np.sum(weights)

    # Convert to the specified type
    if nType == "float":
        weights = weights.astype(np.float32)
    elif nType == "double":
        weights = weights.astype(np.float64)
    else:
        raise ValueError("nType must be 'float' or 'double'.")

    return weights



def smooth_with_gaussian(var, numiter, sigma=1, truncate=4.0):

    for ii in range(numiter):
        print(f"SMOOTH ITER: {ii+1}")
        var = gaussian_filter(var, sigma=sigma, truncate=truncate)

    return var


def print_and_return_varsize(ps_fv, u_fv, v_fv, t_fv, q_fv):
    """
    Function to print and return the maximum size of xarray DataArray variables in bytes.

    Parameters:
    -----------
    ps_fv, u_fv, v_fv, t_fv, q_fv : xarray.DataArray
        Variables whose sizes are to be determined.

    Returns:
    --------
    max_size : int
        The maximum size among the provided variables in bytes.
    """

    ps_sz = ps_fv.nbytes
    u_sz = u_fv.nbytes
    v_sz = v_fv.nbytes
    t_sz = t_fv.nbytes
    q_sz = q_fv.nbytes

    print("~~~ VAR SIZES (bytes)")
    print(f"ps: {ps_sz}")
    print(f"u: {u_sz}")
    print(f"v: {v_sz}")
    print(f"t: {t_sz}")
    print(f"q: {q_sz}")

    max_size = max(ps_sz, u_sz, v_sz, t_sz, q_sz)
    print(f"Var max size: {max_size}")

    return max_size

def print_min_max_dict(data_vars):
    """
    Prints the minimum and maximum values for each key in the data_vars dictionary.

    Parameters:
    -----------
    data_vars : dict
        A dictionary where each key corresponds to a variable name and its value is a NumPy array or xarray DataArray.
    """
    print("="*65)
    print("************")
    for key, data in data_vars.items():
        if isinstance(data, (np.ndarray, xr.DataArray)):
            print(f"{key.upper()} shape: {data.shape} | Max: {data.max()} | Min: {data.min()}")
    print("="*65)



















