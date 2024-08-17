import numpy as np
import xarray as xr
import datetime
import argparse
import sys
import glob
from constants import grav
import cftime

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

def add_time_define_precision(var_in, precision, isncol):
    """
    Adds a time dimension and converts the precision of the input variable.

    Parameters:
    var_in (xarray.DataArray or numpy.ndarray): Input variable.
    precision (str): Desired precision ("float", "single", "double").
    isncol (bool): True if the variable is ncol; False otherwise.

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
        if len(var_dims) == 2:  # nlat, nlon -> time, nlat, nlon
            print(f"Converting a nlat, nlon -> nlat, nlon {precision} with time attached!")
            nlat, nlon = var_dims
            if precision in ["float", "single"]:
                var_out = xr.DataArray(np.zeros((1, nlat, nlon), dtype=np.float32), dims=["time", "lat", "lon"])
            elif precision == "double":
                var_out = xr.DataArray(np.zeros((1, nlat, nlon), dtype=np.float64), dims=["time", "lat", "lon"])
            else:
                raise ValueError("Invalid precision specified")
            var_out[0, :, :] = var_in.astype(var_out.dtype)

        elif len(var_dims) == 3:  # lev, nlat, nlon -> time, lev, nlat, nlon
            print(f"Converting a lev, nlat, nlon -> lev, nlat, nlon {precision} with time attached!")
            nlev, nlat, nlon = var_dims
            if precision in ["float", "single"]:
                var_out = xr.DataArray(np.zeros((1, nlev, nlat, nlon), dtype=np.float32), dims=["time", "lev", "lat", "lon"])
            elif precision == "double":
                var_out = xr.DataArray(np.zeros((1, nlev, nlat, nlon), dtype=np.float64), dims=["time", "lev", "lat", "lon"])
            else:
                raise ValueError("Invalid precision specified")
            var_out[0, :, :, :] = var_in.astype(var_out.dtype)

    # Copy any metadata over to the output variable
    var_out.attrs = var_in.attrs

    return var_out

import cftime
import numpy as np
import xarray as xr



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








def ncol_to_latlon(var_out, nlat, nlon):
    vardims = var_out.shape

    if len(vardims) == 1:
        var_in = var_out.reshape((nlat, nlon), order='F')
    elif len(vardims) == 2:
        nlev = vardims[0]
        ncol = vardims[1]
        print(f"repacking -> nlev: {nlev}    nlat: {nlat}    nlon: {nlon}    ncol: {ncol}")
        var_in = np.empty((nlev, nlat, nlon), dtype=var_out.dtype)
        for ii in range(nlev):
            var_in[ii, :, :] = var_out[ii, :].reshape((nlat, nlon), order='F')
    else:
        raise ValueError(f"{vardims} dims not supported")

    return var_in

def repack_fv(ps_fv, t_fv, q_fv, u_fv, v_fv, cldice_fv, cldliq_fv, dim_sePS, correct_or_not=None):
    """
    Repack FV variables back to lat-lon dimensions.

    Parameters:
    - ps_fv: Surface pressure field.
    - t_fv: Temperature field.
    - q_fv: Specific humidity field.
    - u_fv: U wind component field.
    - v_fv: V wind component field.
    - cldice_fv: Cloud ice field.
    - cldliq_fv: Cloud liquid field.
    - dim_sePS: Tuple containing (nlat, nlon) dimensions.
    - correct_or_not: Optional field to indicate corrections (if available).

    Returns:
    - Tuple of repacked variables in the same order as the input.
    """
    nlat, nlon = dim_sePS

    print("Repacking FV variables...")

    ps_fv = ncol_to_latlon(ps_fv, nlat, nlon)
    t_fv = ncol_to_latlon(t_fv, nlat, nlon)
    q_fv = ncol_to_latlon(q_fv, nlat, nlon)
    u_fv = ncol_to_latlon(u_fv, nlat, nlon)
    v_fv = ncol_to_latlon(v_fv, nlat, nlon)
    cldice_fv = ncol_to_latlon(cldice_fv, nlat, nlon)
    cldliq_fv = ncol_to_latlon(cldliq_fv, nlat, nlon)

    if correct_or_not is not None:
        correct_or_not = ncol_to_latlon(correct_or_not, nlat, nlon)
        return ps_fv, t_fv, q_fv, u_fv, v_fv, cldice_fv, cldliq_fv, correct_or_not

    return ps_fv, t_fv, q_fv, u_fv, v_fv, cldice_fv, cldliq_fv