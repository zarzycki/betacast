import numpy as np
import xarray as xr
import datetime
import argparse
import sys
import glob

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
    parser.add_argument('--adjust_config', type=str, default='',
                        help='String defining how to perform hydro adjustment (default: "")')
    parser.add_argument('--model_topo_file', type=str, default='',
                        help='File containing PHIS for FV or SE, or an MPAS inic file (default: "")')
    parser.add_argument('--mod_in_topo', type=str, default='',
                        help='Full path to PHIS field from downscaling MOD for pressure surface calculation (default: "")')
    parser.add_argument('--mod_remap_file', type=str, default='',
                        help='Full path to ESMF weight file that goes downscaling MOD -> ANL (default: "")')
    parser.add_argument('--mpasfile', type=str, default='',
                        help='XXXXXXX')

    return parser.parse_args()

def define_constants():
    global gamma_s, gamma_d, grav, Rd, Rv, Rv_over_Rd
    gamma_s = 5.0 / 1000.0  # Environmental lapse rate in K/m
    gamma_d = 9.8 / 1000.0  # Adiabatic lapse rate in K/m
    grav = 9.80616
    Rd = 287.058
    Rv = 461.6
    Rv_over_Rd = Rv / Rd

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

def load_grb_file(data_filename):
    # Load the data file using xarray
    return xr.open_dataset(data_filename)

def process_grb_file(grb_file, datasource, RDADIR, yearstr, monthstr, daystr, mod_remap_file):
    print("---------------------------------------------------------")
    print("Loading lat/lon/lev")

    if datasource in ["GFS", "CFSR", "HWRF"]:
        cldlevName = grb_file["CLWMR_P0_L100_GLL0"].dims[0]
        print(f"cldlev varname: {cldlevName}")
        rhlevName = grb_file["RH_P0_L100_GLL0"].dims[0]
        print(f"rhlev varname: {rhlevName}")
        grblat = grb_file["lat_0"].values
        grblon = grb_file["lon_0"].values
        grblev = grb_file["lv_ISBL0"].values
        rhlev = grb_file[rhlevName].values
        cldlev = grb_file[cldlevName].values

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

def load_variable(varname, sf_dir, pl_dir, var_code, yearstr, monthstr, daystr, cyclestr):
    """Helper function to load a variable from a NetCDF file."""
    rda_find = glob.glob(f"{pl_dir}/{var_code}.{yearstr}{monthstr}{daystr}00_*.nc") or \
               glob.glob(f"{sf_dir}/{var_code}.{yearstr}{monthstr}0100_*.nc")
    if not rda_find:
        raise FileNotFoundError(f"No matching files found for {var_code}")
    rda_file = xr.open_dataset(rda_find[0])
    rda_time = rda_file["time"]
    rda_thistime = find_closest_time(rda_time, yearstr, monthstr, daystr, cyclestr)
    rda_data = rda_file[varname].sel(time=rda_thistime, method='nearest').values
    return rda_data

def load_era5rda_data(RDADIR, yearstr, monthstr, daystr, cyclestr, dycore):
    # Define directories
    pl_dir = f"{RDADIR}/e5.oper.an.pl/{yearstr}{monthstr}"
    sf_dir = f"{RDADIR}/e5.oper.an.sfc/{yearstr}{monthstr}"

    # Dictionary to store the variables
    data_vars = {}

    # Load required variables
    data_vars['ps'] = load_variable('SP', sf_dir, pl_dir, "e5.oper.an.sfc.128_134_sp.ll025sc", yearstr, monthstr, daystr, cyclestr)
    data_vars['t_gfs'] = load_variable('T', pl_dir, pl_dir, "e5.oper.an.pl.128_130_t.ll025sc", yearstr, monthstr, daystr, cyclestr)
    data_vars['u_gfs'] = load_variable('U', pl_dir, pl_dir, "e5.oper.an.pl.128_131_u.ll025uv", yearstr, monthstr, daystr, cyclestr)
    data_vars['v_gfs'] = load_variable('V', pl_dir, pl_dir, "e5.oper.an.pl.128_132_v.ll025uv", yearstr, monthstr, daystr, cyclestr)
    data_vars['q_gfs'] = load_variable('Q', pl_dir, pl_dir, "e5.oper.an.pl.128_133_q.ll025sc", yearstr, monthstr, daystr, cyclestr)
    data_vars['cldliq_gfs'] = load_variable('CLWC', pl_dir, pl_dir, "e5.oper.an.pl.128_246_clwc.ll025sc", yearstr, monthstr, daystr, cyclestr)
    data_vars['cldice_gfs'] = load_variable('CIWC', pl_dir, pl_dir, "e5.oper.an.pl.128_247_ciwc.ll025sc", yearstr, monthstr, daystr, cyclestr)

    if dycore == "mpas":
        data_vars['w_gfs'] = load_variable('W', pl_dir, pl_dir, "e5.oper.an.pl.128_135_w.ll025sc", yearstr, monthstr, daystr, cyclestr)
        data_vars['w_is_omega'] = True
        data_vars['z_gfs'] = load_variable('Z', pl_dir, pl_dir, "e5.oper.an.pl.128_129_z.ll025sc", yearstr, monthstr, daystr, cyclestr)
        data_vars['z_is_phi'] = True

    return data_vars

def find_closest_time(times, yearstr, monthstr, daystr, cyclestr):
    target_time = np.datetime64(f"{yearstr}-{monthstr}-{daystr}T{cyclestr[:2]}:00")
    closest_time = min(times.values, key=lambda x: abs(x - target_time))
    return closest_time






def load_cam_levels(PATHTOHERE, numlevels):
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

