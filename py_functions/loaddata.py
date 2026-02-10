import numpy as np
import xarray as xr
import glob
import sys
import cfgrib

import pyfuncs
import vertremap
import cam2cam
import horizremap
import meteo

from constants import (
    grav, p0, t_freeze_K, dtime_map, OMEGA_LAT_THRESH
)

import logging
logger = logging.getLogger(__name__)


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
    template_path = f"{PATHTOHERE}/L{numlevels}template.nc"
    fC = xr.open_dataset(template_path)

    # Extract the hybrid coefficients and levels
    if load_xarray:
        data = {
            'hya':  fC["hyam"],
            'hyb':  fC["hybm"],
            'hyai': fC["hyai"],
            'hybi': fC["hybi"],
            'lev':  fC["lev"],
            'ilev': fC["ilev"],
        }
    else:
        data = {
            'hya':  fC["hyam"].values,
            'hyb':  fC["hybm"].values,
            'hyai': fC["hyai"].values,
            'hybi': fC["hybi"].values,
            'lev':  fC["lev"].values,
            'ilev': fC["ilev"].values,
        }

    logging.info("---------------------------------------------------------")
    logging.info("Loading CAM levels")
    logging.info(f"Loading {numlevels} level data")

    return data

def load_ERA5RDA_variable(varname, the_dir, var_code, yearstr, monthstr, daystr, cyclestr, return_coords=False, return_hycoef=False):
    """
    Helper function to load a variable from an ERA5 RDA file file.

    Parameters:
    -----------
    return_coords : bool, optional (default False)
        If True, return coordinates along with the data.
    return_hycoef : bool, optional (default False)
        If True, return hybrid coefficients instead of coordinates.
        This option overrides return_coords if both are set to True.

    Returns:
    --------
    If return_hycoef is True:
        rda_data, hyam, hybm, ps
    Elif return_coords is True:
        rda_data, latitude, longitude, level
    Else:
        rda_data
    """
    rda_find = glob.glob(f"{the_dir}/{var_code}.{yearstr}{monthstr}{daystr}00_*.nc") or \
               glob.glob(f"{the_dir}/{var_code}.{yearstr}{monthstr}0100_*.nc")
    if not rda_find:
        raise FileNotFoundError(f"No matching files found for {var_code}")
    rda_file = xr.open_dataset(rda_find[0])
    rda_time = rda_file["time"]
    rda_thistime = pyfuncs.find_closest_time(rda_time, yearstr, monthstr, daystr, cyclestr)
    rda_data = rda_file[varname].sel(time=rda_thistime, method='nearest').values

    logging.debug(f"load_ERA5RDA_variable: Getting {varname} from {rda_find[0]}")

    if return_hycoef:
        # ERA5 "interface" is "half-level"
        # Formula: "p = a + b*ps" ;
        hyam = rda_file["a_model"].values / p0 if "a_model" in rda_file else None
        hybm = rda_file["b_model"].values if "b_model" in rda_file else None
        hyai = rda_file["a_half"].values / p0 if "a_half" in rda_file else None
        hybi = rda_file["b_half"].values if "b_half" in rda_file else None

        if hyam is None or hybm is None or hyai is None or hybi is None:
            raise ValueError("All hybrid coefficients not found in the file")

        return rda_data, rda_file["latitude"].values.astype(float), rda_file["longitude"].values.astype(float), rda_file["level"].values.astype(float) * 100.0, hyam, hybm, hyai, hybi
    elif return_coords:
        return rda_data, rda_file["latitude"].values.astype(float), rda_file["longitude"].values.astype(float), rda_file["level"].values.astype(float) * 100.0
    else:
        return rda_data


def load_SAMPLE_data(data_filename, dycore=None):

    # Dictionary to store the variables
    data_vars = {}

    # Load the GRIB file
    grb_file = load_ERA5_file(data_filename)

    data_vars['ps'] = grb_file.sp[0,:,:].values
    data_vars['t'] = grb_file.t[0,:,:,:].values
    data_vars['u'] = grb_file.u[0,:,:,:].values
    data_vars['v'] = grb_file.v[0,:,:,:].values
    data_vars['q'] = grb_file.q[0,:,:,:].values
    data_vars['cldice'] = grb_file.ciwc[0,:,:,:].values
    data_vars['cldliq'] = grb_file.clwc[0,:,:,:].values

    data_vars['lat'] = grb_file.latitude.values.astype(float)
    data_vars['lon'] = grb_file.longitude.values.astype(float)
    data_vars['lev'] = grb_file.pressure_level.values.astype(float) * 100.

    data_vars['ts'] = grb_file.t2m[0,:,:].values
    data_vars['phis'] = grb_file.phis[0,:,:].values

    # Load additional variables for MPAS dycore
    if dycore == 'mpas':
        data_vars['w'] = grb_file.w[0,:,:,:].values
        data_vars['w_is_omega'] = True
        data_vars['z'] = grb_file.z[0,:,:,:].values
        data_vars['z_is_phi'] = True

    # Flip ordering of levels to go from bottom to top
    data_vars = flip_level_dimension(data_vars)

    # Close the file
    grb_file.close()

    return data_vars


def load_ERA5RDA_data(RDADIR, data_filename, yearstr, monthstr, daystr, cyclestr, dycore, get_chemistry=False):
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

    # Create a 3-D pressure field from constant pressure surfaces
    data_vars['pres'] = np.broadcast_to(data_vars['lev'][:, np.newaxis, np.newaxis], data_vars['t'].shape)

    if dycore == 'mpas':
        data_vars['w'] = load_ERA5RDA_variable('W', pl_dir, "e5.oper.an.pl.128_135_w.ll025sc", yearstr, monthstr, daystr, cyclestr)
        data_vars['w_is_omega'] = True
        data_vars['z'] = load_ERA5RDA_variable('Z', pl_dir, "e5.oper.an.pl.128_129_z.ll025sc", yearstr, monthstr, daystr, cyclestr)
        data_vars['z_is_phi'] = True

    if get_chemistry:
        data_vars['o3'] = load_ERA5RDA_variable('O3', pl_dir, "e5.oper.an.pl.128_203_o3.ll025sc", yearstr, monthstr, daystr, cyclestr)

    data_vars['ts'] = load_ERA5RDA_variable('VAR_2T', sf_dir, "e5.oper.an.sfc.128_167_2t.ll025sc", yearstr, monthstr, daystr, cyclestr)

    ds = load_ERA5_file(data_filename)
    data_vars['phis'] = ds["Z"].isel(time=0).values
    ds.close()

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

def open_nc_dataset(data_filename):
    """Open a NetCDF file as an xarray Dataset."""
    return xr.open_dataset(data_filename)

# Aliases for readability at call sites
load_ERA5_file = open_nc_dataset
load_CAM_file = open_nc_dataset


def inspect_GRIB_file(data_filename):
    """Safely inspect GRIB file contents by opening all datasets separately"""
    datasets = cfgrib.open_datasets(data_filename)

    print(f"Found {len(datasets)} separate datasets in the file:")
    for i, ds in enumerate(datasets):
        print(f"\nDataset {i}:")
        print(f"  Variables: {list(ds.data_vars.keys())}")
        print(f"  Coordinates: {list(ds.coords.keys())}")
        if 'typeOfLevel' in ds.attrs:
            print(f"  typeOfLevel: {ds.attrs['typeOfLevel']}")
        if hasattr(ds, 'level') and 'level' in ds.coords:
            print(f"  Levels: {ds.coords['level'].values}")

    # Close all datasets
    for ds in datasets:
        ds.close()


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

def load_CFSR_data(grb_file_name, dycore, get_chemistry=False):

    # Inspect the file
    #inspect_GRIB_file(grb_file_name)

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

    if get_chemistry:
        data_vars['o3'] = load_and_extract_CFSR_variable(grb_file_name, 'isobaricInhPa', 'o3mr')

    # Create a 3-D pressure field from constant pressure surfaces
    data_vars['pres'] = np.broadcast_to(data_vars['lev'][:, np.newaxis, np.newaxis], data_vars['t'].shape)

    if dycore == "mpas":
        data_vars['w'] = load_and_extract_CFSR_variable(grb_file_name, 'isobaricInhPa', 'w')
        data_vars['w_is_omega'] = True
        data_vars['z'] = load_and_extract_CFSR_variable(grb_file_name, 'isobaricInhPa', 'gh')
        data_vars['z_is_phi'] = False

    if data_vars['windlev'].shape != data_vars['lev'].shape:
        logging.info(f"Interpolating u + v wind from {data_vars['windlev'].shape} to {data_vars['lev'].shape} levels")
        data_vars['u'] = vertremap.int2p_n(data_vars['windlev'], data_vars['u'], data_vars['lev'], linlog=2, dim=0)
        data_vars['v'] = vertremap.int2p_n(data_vars['windlev'], data_vars['v'], data_vars['lev'], linlog=2, dim=0)

    # Load and interpolate additional variables
    data_vars['rh_native'] = load_and_extract_CFSR_variable(grb_file_name, 'isobaricInhPa', 'r')
    data_vars['cldmix_native'] = load_and_extract_CFSR_variable(grb_file_name, 'isobaricInhPa', 'clwmr')

    # Interpolate
    data_vars['rh'] = vertremap.int2p_n(data_vars['rhlev'],data_vars['rh_native'], data_vars['lev'], linlog=2, dim=0)
    data_vars['cldmix'] = vertremap.int2p_n(data_vars['cldlev'],data_vars['cldmix_native'], data_vars['lev'], linlog=2, dim=0)

    # Calculate specific humidity from RH
    # Load grblev (Pa), use numpy broadcasting to go from 1-D to 3-D, multiple by 0.01 to convert to mb
    data_vars['q'] = meteo.mixhum_ptrh((data_vars['lev'][:, np.newaxis, np.newaxis] * 0.01), data_vars['t'], data_vars['rh'])

    # Sort bad values
    data_vars['cldmix'] = np.where(np.isnan(data_vars['cldmix']), 0, data_vars['cldmix'])
    data_vars['cldmix'] = pyfuncs.clip_and_count(data_vars['cldmix'],max_thresh=0.01,var_name="cldmix")

    # Separate cloud ice and water
    data_vars['cldice'] = np.where(data_vars['t'] > t_freeze_K, 0, data_vars['cldmix'])
    data_vars['cldliq'] = np.where(data_vars['t'] > t_freeze_K, data_vars['cldmix'], 0)

    # Surface geopotential
    data_vars['phis']  = load_and_extract_CFSR_variable(grb_file_name, 'surface', 'orog')
    data_vars['phis'] = data_vars['phis'] * grav

    # Surface temperature -- try "surface" first, then 0.995 sigma if surface unavailable.
    try:
        data_vars['ts'] = load_and_extract_CFSR_variable(grb_file_name, 'surface', 't')
    except Exception as e_surface:
        logging.warning(f"Failed to load variable with 'sigma' level: {e_surface}")
        try:
            data_vars['ts'] = load_and_extract_CFSR_variable(grb_file_name, 'sigma', 't')
        except Exception as e_sigma:
            logging.error(f"Failed to load variable with 'surface' level: {e_sigma}")
            logging.error("Both attempts to load the variable failed. Exiting script.")
            sys.exit(1)

    pyfuncs.print_min_max_dict(data_vars)

    # Cleanup
    del data_vars['cldmix']
    del data_vars['cldlev']
    del data_vars['rhlev']
    del data_vars['windlev']
    del data_vars['rh']
    del data_vars['cldmix_native']
    del data_vars['rh_native']

    # Because GFS is bottom -> top, we have to flip to be top-to-bottom
    data_vars = flip_level_dimension(data_vars)

    return data_vars

def load_HRRRml_data(grb_file_name, dycore):
    """
    Load HRRR data and interpolate to constant pressure levels.

    HRRR has variables on hybrid levels that need to be interpolated to
    constant pressure surfaces using the 3D pressure field.
    """

    # Inspect the file
    #inspect_GRIB_file(grb_file_name)

    # Dictionary to store the variables
    data_vars = {}

    # Use surface variables to get ref coordinates
    data_vars['lat'] = get_CFSR_coords_for_var(grb_file_name, 'sp', 'latitude', 'surface')
    data_vars['lon'] = get_CFSR_coords_for_var(grb_file_name, 'sp', 'longitude', 'surface')

    # This really should be edited to get the nominal levels but I can't find the
    # hybrid sigma coords for HRRR
    #data_vars['lev'] = vertremap.hrrr_eta_to_plev()     # HRRRv1 and v2
    data_vars['lev'] = vertremap.__pres_lev_hrrr__       # HRRRv3 and v4

    # Load surface pressure
    data_vars['ps'] = load_and_extract_CFSR_variable(grb_file_name, 'surface', 'sp')

    # Load surface geopotential
    data_vars['phis'] = load_and_extract_CFSR_variable(grb_file_name, 'surface', 'orog')
    data_vars['phis'] = data_vars['phis'] * grav

    # Load 3D pressure field on hybrid levels
    pres_3d_hybrid = load_and_extract_CFSR_variable(grb_file_name, 'hybrid', 'pres')

    # Load atmospheric variables on hybrid levels
    t_hybrid = load_and_extract_CFSR_variable(grb_file_name, 'hybrid', 't')
    u_hybrid = load_and_extract_CFSR_variable(grb_file_name, 'hybrid', 'u')
    v_hybrid = load_and_extract_CFSR_variable(grb_file_name, 'hybrid', 'v')
    q_hybrid = load_and_extract_CFSR_variable(grb_file_name, 'hybrid', 'q')

    # Optional variables for MPAS dycore
    if dycore == "mpas":
        w_hybrid = load_and_extract_CFSR_variable(grb_file_name, 'hybrid', 'w')
        # For geopotential height, try hybrid first, then isobaric
        try:
            z_hybrid = load_and_extract_CFSR_variable(grb_file_name, 'hybrid', 'gh')
            z_is_on_hybrid = True
        except Exception:
            logging.warning("Geopotential height not found on hybrid levels, trying isobaric...")
            z_isobaric = load_and_extract_CFSR_variable(grb_file_name, 'isobaricInhPa', 'gh')
            z_isobaric_levels = get_CFSR_levels_for_var(grb_file_name, 'gh')
            z_is_on_hybrid = False

    # Cloud variables on hybrid levels
    try:
        cldmix_hybrid = load_and_extract_CFSR_variable(grb_file_name, 'hybrid', 'clwmr')
    except Exception:
        logging.warning("Cloud variables not found on hybrid levels, setting to None")
        cldmix_hybrid = None

    # Now interpolate all hybrid-level variables to constant pressure levels
    logging.info("Interpolating variables from hybrid to pressure levels...")
    logging.debug(f"Temperature shape: {t_hybrid.shape}")
    logging.debug(f"Pressure shape: {pres_3d_hybrid.shape}")

    # Determine level dimension (0 or 1)
    level_dim = 1
    if len(t_hybrid.shape) == 4:  # (time, lev, lat, lon)
        level_dim = 1
    elif len(t_hybrid.shape) == 3:  # (lev, lat, lon)
        level_dim = 0

    logging.debug(f"Using level dimension: {level_dim}")

    # Interpolate main atmospheric variables
    data_vars['t'] = vertremap.int2p_n_extrap(
        pin=pres_3d_hybrid, xin=t_hybrid, pout=data_vars['lev'],
        linlog=1, dim=level_dim, flip_p=True,
        extrapolate=True, variable='other', ps=data_vars['ps'], phi_sfc=data_vars['phis']
    )

    data_vars['u'] = vertremap.int2p_n_extrap(
        pin=pres_3d_hybrid, xin=u_hybrid, pout=data_vars['lev'],
        linlog=1, dim=level_dim, flip_p=True,
        extrapolate=True, variable='other', ps=data_vars['ps'], phi_sfc=data_vars['phis']
    )

    data_vars['v'] = vertremap.int2p_n_extrap(
        pin=pres_3d_hybrid, xin=v_hybrid, pout=data_vars['lev'],
        linlog=1, dim=level_dim, flip_p=True,
        extrapolate=True, variable='other', ps=data_vars['ps'], phi_sfc=data_vars['phis']
    )

    data_vars['q'] = vertremap.int2p_n_extrap(
        pin=pres_3d_hybrid, xin=q_hybrid, pout=data_vars['lev'],
        linlog=1, dim=level_dim, flip_p=True,
        extrapolate=True, variable='other', ps=data_vars['ps'], phi_sfc=data_vars['phis']
    )

    # Create a 3-D pressure field from constant pressure surfaces
    data_vars['pres'] = np.broadcast_to(data_vars['lev'][:, np.newaxis, np.newaxis], data_vars['t'].shape)

    # MPAS-specific variables
    if dycore == "mpas":
        data_vars['w'] = vertremap.int2p_n_extrap(
            pin=pres_3d_hybrid, xin=w_hybrid, pout=data_vars['lev'],
            linlog=1, dim=level_dim, flip_p=True,
            extrapolate=True, variable='other', ps=data_vars['ps'], phi_sfc=data_vars['phis']
        )
        data_vars['w_is_omega'] = True

        if z_is_on_hybrid:
            data_vars['z'] = vertremap.int2p_n_extrap(
                pin=pres_3d_hybrid, xin=z_hybrid, pout=data_vars['lev'],
                linlog=1, dim=level_dim, flip_p=True,
                extrapolate=True, variable='other', ps=data_vars['ps'], phi_sfc=data_vars['phis']
            )
        else:
            # Interpolate from isobaric to our target pressure levels
            data_vars['z'] = vertremap.int2p_n_extrap(
                pin=z_isobaric_levels[:, np.newaxis, np.newaxis],
                xin=z_isobaric, pout=data_vars['lev'],
                linlog=1, dim=level_dim, flip_p=True,
                extrapolate=True, variable='other', ps=data_vars['ps'], phi_sfc=data_vars['phis']
            )
        data_vars['z_is_phi'] = False

    if cldmix_hybrid is not None:
        data_vars['cldmix'] = vertremap.int2p_n_extrap(
            pin=pres_3d_hybrid, xin=cldmix_hybrid, pout=data_vars['lev'],
            linlog=1, dim=level_dim, flip_p=True,
            extrapolate=True, variable='other', ps=data_vars['ps'], phi_sfc=data_vars['phis']
        )
    else:
        # Create dummy cloud mixing ratio if not available
        logging.warning("Cloud mixing ratio not available, creating dummy field")
        data_vars['cldmix'] = np.zeros_like(data_vars['t'])

    # Clean up cloud data
    data_vars['cldmix'] = np.where(np.isnan(data_vars['cldmix']), 0, data_vars['cldmix'])
    data_vars['cldmix'] = pyfuncs.clip_and_count(data_vars['cldmix'], max_thresh=0.01, var_name="cldmix")

    # Separate cloud ice and water
    data_vars['cldice'] = np.where(data_vars['t'] > t_freeze_K, 0, data_vars['cldmix'])
    data_vars['cldliq'] = np.where(data_vars['t'] > t_freeze_K, data_vars['cldmix'], 0)

    # Surface temperature - try surface first
    try:
        data_vars['ts'] = load_and_extract_CFSR_variable(grb_file_name, 'surface', 't')
    except Exception as e_surface:
        logging.warning(f"Failed to load surface temperature: {e_surface}")
        try:
            # Try 2m temperature as backup
            data_vars['ts'] = load_and_extract_CFSR_variable(grb_file_name, 'heightAboveGround', 't2m')
            logging.info("Using 2m temperature as surface temperature")
        except Exception as e_2m:
            logging.warning(f"Failed to load 2m temperature: {e_2m}")
            logging.warning("Setting surface temperature from lowest atmospheric level")
            data_vars['ts'] = data_vars['t'][-1, :, :]  # Use lowest pressure level

    # Print diagnostics
    pyfuncs.print_min_max_dict(data_vars)

    # Clean up temporary variables
    if 'cldmix' in data_vars:
        del data_vars['cldmix']

    # Because GFS is bottom -> top, we have to flip to be top-to-bottom
    logging.debug(f"lev before flip: {data_vars['lev']}")
    data_vars = flip_level_dimension(data_vars)
    logging.debug(f"lev after flip: {data_vars['lev']}")

    logging.info(f"Successfully loaded HRRR data interpolated to {len(data_vars['lev'])} pressure levels")

    return data_vars





def flip_level_dimension(data_vars, state_dimensions=3):
    """
    Flips the 'lev' dimension and all relevant 3-D arrays in data_vars if 'lev' is in descending order.

    Parameters:
    data_vars (dict): A dictionary containing arrays, where 'lev' is a key representing the levels and
                      3-D arrays need to be flipped if 'lev' is descending.

    Returns:
    dict: An updated version of data_vars with 'lev' and relevant 3-D arrays flipped if needed.
    """
    if np.all(np.diff(data_vars['lev']) < 0):
        logging.info("'lev' is in descending order, flipping all 3-D variables and 'lev' array.")

        # Flip 'lev' array
        data_vars['lev'] = data_vars['lev'][::-1]

        # Flip all 3-D arrays in data_vars
        for var_name, var_value in data_vars.items():
            if isinstance(var_value, np.ndarray) and var_value.ndim == state_dimensions:
                logging.debug(f"Flipping 3-D array: {var_name}")
                data_vars[var_name] = var_value[::-1, :, :]

    return data_vars


def load_ERA5mlRDA_data(RDADIR, data_filename, VERT_COORD_PATH, yearstr, monthstr, daystr, cyclestr, dycore):
    # Define directories
    ml_dir = f"{RDADIR}/e5.oper.an.ml/{yearstr}{monthstr}"

    # Dictionary to store the variables
    data_vars = {}

    # Get levs data from a CSV file
    # https://confluence.ecmwf.int/display/UDOC/L137+model+level+definitions
    data_vars['lev'] = np.array(pyfuncs.load_csv_column(f"{VERT_COORD_PATH}/era5_levs.csv", "pf [hPa]"), dtype=np.float32)
    data_vars['lev'] = data_vars['lev'] * 100.    # mb to Pa
    logging.debug(f"ERA5 model level nominal plevs: {data_vars['lev'][:3]} ... {data_vars['lev'][-3:]}")

    _, data_vars['lat'], data_vars['lon'], _, hyam, hybm, hyai, hybi = load_ERA5RDA_variable('T', ml_dir, "e5.oper.an.ml.0_5_0_0_0_t.regn320sc", yearstr, monthstr, daystr, cyclestr, return_hycoef=True)

    # Load required variables
    data_vars['ps'] = load_ERA5RDA_variable('SP', ml_dir, "e5.oper.an.ml.128_134_sp.regn320sc", yearstr, monthstr, daystr, cyclestr)
    data_vars['t'] = load_ERA5RDA_variable('T', ml_dir, "e5.oper.an.ml.0_5_0_0_0_t.regn320sc", yearstr, monthstr, daystr, cyclestr)
    data_vars['u'] = load_ERA5RDA_variable('U', ml_dir, "e5.oper.an.ml.0_5_0_2_2_u.regn320uv", yearstr, monthstr, daystr, cyclestr)
    data_vars['v'] = load_ERA5RDA_variable('V', ml_dir, "e5.oper.an.ml.0_5_0_2_3_v.regn320uv", yearstr, monthstr, daystr, cyclestr)
    data_vars['q'] = load_ERA5RDA_variable('Q', ml_dir, "e5.oper.an.ml.0_5_0_1_0_q.regn320sc", yearstr, monthstr, daystr, cyclestr)

    if dycore == 'mpas':
        data_vars['w'] = load_ERA5RDA_variable('W', ml_dir, "e5.oper.an.ml.0_5_0_2_8_w.regn320sc", yearstr, monthstr, daystr, cyclestr)
        data_vars['w_is_omega'] = True

    # Get lowest model level t for ts
    data_vars['ts'] = data_vars['t'][-1, :, :]

    ds = load_ERA5_file(data_filename)
    data_vars['phis'] = ds["Z"].isel(time=0).values
    ds.close()

    if dycore == 'mpas':
        data_vars['w_is_omega'] = False
        tkv = data_vars['t'] * (1. + 0.61 * data_vars['q'])
        data_vars['z'] = cam2cam.cz2ccm(data_vars['ps'], data_vars['phis'], tkv, p0, hyam[::-1], hybm[::-1], hyai[::-1], hybi[::-1])
        data_vars['z_is_phi'] = False

    # Convert dry to wet?
    # data_vars['psdry'], data_vars['pw'] = meteo.ps_wet_to_dry_conversion(data_vars['ps'], data_vars['q'], hyai, hybi, p0, verbose=True)

    # Vertically interpolate the hybrid levels to constant pressure surfaces
    # We are essentially going to set the p levs to the nominal p levs from ECMWF
    data_vars = vertremap.interp_hybrid_to_pressure_wrapper(
        data_vars=data_vars,
        ps=data_vars['ps'],
        hyam=hyam,
        hybm=hybm,
        new_levels=data_vars['lev']
        )

    # Create a 3-D pressure field from constant pressure surfaces
    data_vars['pres'] = np.broadcast_to(data_vars['lev'][:, np.newaxis, np.newaxis], data_vars['t'].shape)

    # No cloud ice/liq, so set these to zero same shape as interpolated t
    data_vars['cldice'] = np.zeros_like(data_vars['t'])
    data_vars['cldliq'] = np.zeros_like(data_vars['t'])

    pyfuncs.print_min_max_dict(data_vars)

    return data_vars










def load_cam_data(grb_file_name, YYYYMMDDHH, mod_in_topo, mod_remap_file, dycore, write_debug_files=False, write_debug_dir="./"):

    data_vars = {}

    grb_file = load_CAM_file(grb_file_name)

    data_vars['lev'] = np.array([
        20, 30, 50, 100, 200, 300, 500, 700, 1000, 2000, 3000, 5000, 7000, 10000,
        15000, 20000, 25000, 30000, 35000, 40000, 45000, 50000, 55000, 60000,
        65000, 70000, 75000, 80000, 82500, 85000, 87500, 90000, 92000, 93000, 94000,
        95000, 96000, 97000, 98000, 99000,
        99500, 100000, 100500, 101000, 102000, 103000
    ], dtype=np.float32)

    yearstr, monthstr, daystr, cyclestr = pyfuncs.split_by_lengths(str(YYYYMMDDHH), dtime_map)

    cam_time = grb_file["time"]

    cam_thistime_ix = pyfuncs.find_closest_time(cam_time.values, yearstr, monthstr, daystr, cyclestr, return_isel=True)
    logging.info(f"Closest time: {cam_thistime_ix}")
    cam_thistime = pyfuncs.find_closest_time(cam_time, yearstr, monthstr, daystr, cyclestr)
    logging.info(f"Closest time: {cam_thistime}")

    # Select the closest time
    data_vars['ps'] = grb_file["PS"].sel(time=cam_thistime, method='nearest').values
    data_vars['t'] = grb_file["T"].sel(time=cam_thistime, method='nearest').values
    data_vars['u'] = grb_file["U"].sel(time=cam_thistime, method='nearest').values
    data_vars['v'] = grb_file["V"].sel(time=cam_thistime, method='nearest').values
    data_vars['q'] = grb_file["Q"].sel(time=cam_thistime, method='nearest').values

    logging.info(f"CAM: Using mod_in_topo: {mod_in_topo}")

    # Load topo data
    mod_in_topo_f = xr.open_dataset(mod_in_topo)
    data_vars['phis'] = mod_in_topo_f["PHIS"].values
    data_vars['ts'] = data_vars['t'][-1,:]

    # Horizontal remaps
    data_vars = horizremap.remap_all(data_vars, mod_remap_file, dycore=dycore)

    # Get relevant vertical coefficients from CAM file
    hyam = grb_file["hyam"].values
    hybm = grb_file["hybm"].values
    hyai = grb_file["hyai"].values
    hybi = grb_file["hybi"].values

    pyfuncs.print_min_max_dict(data_vars)

    if dycore == 'mpas':
        # Calculate omega

        data_vars['w_is_omega'] = True
        if grb_file["lat"][0] > grb_file["lat"][1]:
            logging.info("flipping omega since calc needs to be S->N")
            data_vars['w'] = cam2cam.omega_ccm_driver(p0, data_vars['ps'][::-1, :], data_vars['u'][:, ::-1, :], data_vars['v'][:, ::-1, :], data_vars['lat'][::-1], data_vars['lon'], hyam, hybm, hyai, hybi)
            data_vars['w'] = data_vars['w'][:, ::-1, :]  # flip back
        else:
            data_vars['w'] = cam2cam.omega_ccm_driver(p0, data_vars['ps'], data_vars['u'], data_vars['v'], data_vars['lat'], data_vars['lon'], hyam, hybm, hyai, hybi)
        # Remove omega poleward of OMEGA_LAT_THRESH to deal with singularity
        data_vars['w'] = np.where(np.abs(data_vars['lat'])[np.newaxis, :, np.newaxis] > OMEGA_LAT_THRESH, 0.0, data_vars['w'])

        #data_vars['div'] = ddvfidf_wrapper(data_vars['u'], data_vars['v'], data_vars['lat'], data_vars['lon'], 3)

        tkv = data_vars['t'] * (1. + 0.61 * data_vars['q'])
        data_vars['z'] = cam2cam.cz2ccm(data_vars['ps'], data_vars['phis'], tkv, p0, hyam[::-1], hybm[::-1], hyai[::-1], hybi[::-1])
        data_vars['z_is_phi'] = False
        #data_vars['dpsl'], data_vars['dpsm'] = calculate_gradients(data_vars['ps'], data_vars['lat'], data_vars['lon'])
        #data_vars['div'], data_vars['vort'] = calculate_div_vort(data_vars['lat'], data_vars['lon'], data_vars['u'], data_vars['v'])
        #data_vars['pdel'] = dpres_hybrid_ccm(data_vars['ps'], p0, hyai, hybi)
        #data_vars['pmid'] = pres_hybrid_ccm(data_vars['ps'], p0, hyam, hybm)

        if write_debug_files:
            pyfuncs.print_debug_file(
                  write_debug_dir+"/"+"py_cam_raw.nc",
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

    pyfuncs.print_min_max_dict(data_vars)

    # Vertically interpolate the CAM hybrid levels to constant pressure surfaces
    data_vars = vertremap.interp_hybrid_to_pressure_wrapper(
        data_vars=data_vars,
        ps=data_vars['ps'],
        hyam=hyam,
        hybm=hybm,
        new_levels=data_vars['lev']
        )

    # Create a 3-D pressure field from constant pressure surfaces
    data_vars['pres'] = np.broadcast_to(data_vars['lev'][:, np.newaxis, np.newaxis], data_vars['t'].shape)

    data_vars['cldice'] = np.zeros_like(data_vars['t'])
    data_vars['cldliq'] = np.zeros_like(data_vars['t'])

    return data_vars
