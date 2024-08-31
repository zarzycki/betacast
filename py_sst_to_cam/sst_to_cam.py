import argparse
import numpy as np
import xarray as xr
import logging
from datetime import datetime
from scipy.interpolate import RegularGridInterpolator
import sys
import os

module_path = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', 'py_functions'))
if module_path not in sys.path:
    sys.path.append(module_path)
import pyfuncs

def linint2(xi, yi, fi, fiCyclicX, xo, yo):
    """
    Bilinear interpolation from one rectilinear grid to another.

    Parameters:
    - xi: 1D or 2D array of X coordinates of the input data.
    - yi: 1D or 2D array of Y coordinates of the input data.
    - fi: 2D or 3D array of input data to be interpolated.
    - fiCyclicX: Boolean indicating if the X dimension is cyclic.
    - xo: 1D array of X coordinates of the output data.
    - yo: 1D array of Y coordinates of the output data.

    Returns:
    - Interpolated array of the same dimensions as fi, except the last two dimensions.
    """

    # Handle cyclic conditions
    if fiCyclicX:
        # Extend the fi array and xi array for cyclic boundary
        xi_extended = np.concatenate((xi, xi[:1] + 360.0), axis=-1)
        fi_extended = np.concatenate((fi, fi[..., :1]), axis=-1)
        xi = xi_extended
        fi = fi_extended

    # Create the interpolator object
    interpolator = RegularGridInterpolator((yi, xi), fi, method='linear', bounds_error=False, fill_value=np.nan)

    # Create a meshgrid for target coordinates
    xo_mesh, yo_mesh = np.meshgrid(xo, yo, indexing='xy')

    # Prepare the points for interpolation
    target_points = np.vstack((yo_mesh.ravel(), xo_mesh.ravel())).T

    # Perform the interpolation
    interpolated_data = interpolator(target_points)

    # Reshape the interpolated data to match the target grid shape
    interpolated_data = interpolated_data.reshape(yo_mesh.shape)

    return interpolated_data

def parse_initdate(initdate):
    yyyy = int(initdate[:4])
    mm = int(initdate[4:6])
    dd = int(initdate[6:8])
    hh = int(initdate[8:])
    return yyyy, mm, dd, hh

def process_sst_ice(initdate, predict_docn, inputres, datasource, sstDataFile, iceDataFile, SST_write_file):
    # Constants
    TTHRESH = 271.9  # Initial cut for ice vs. open ocean
    KtoC = 273.15    # Kelvin to Celsius conversion
    smooth_ice = True
    smooth_iter = 3

    logging.info(f"inputres set to: {inputres}")

    readinFile = f"./test_files/domains/domain.ocn.{inputres}.nc"

    do_anom = predict_docn == 1
    logging.info(f"do_anom (i.e., predict_docn) set to: {do_anom}")

    yyyy, mm, dd, hh = parse_initdate(initdate)

    logging.info(f"running date {initdate}")

    # Load SST files that we want to overwrite
    logging.info(f"Writing {inputres} SSTs")

    in_ds = xr.open_dataset(readinFile)

    if datasource == "GDAS":
        sst_file = xr.open_dataset(sstDataFile)
        print(sst_file)
        sstlat = sst_file['latitude'].values
        sstlon = sst_file['longitude'].values
        sst_gfs = sst_file['t'].values

        logging.info(f"For GDAS, set all native cells less than {TTHRESH} to ice covered")
        ice_gfs = np.where(sst_gfs > TTHRESH, 0.0, 1.0)

        if smooth_ice:
            logging.info(f"smoothing derived ice field {smooth_iter} times!")
            ice_gfs = pyfuncs.smooth_with_smth9(ice_gfs, smooth_iter, p=0.50, q=0.25)

        sst_gfs -= KtoC  # Convert SST from K to degC

    elif datasource == "NOAAOI":
        sst_file = xr.open_dataset(sstDataFile)
        ice_file = xr.open_dataset(iceDataFile)
        time = sst_file['time']

        date = datetime(yyyy, mm, dd)
        sst_gfs = sst_file.sel(time=date)['sst'].values
        ice_gfs = ice_file.sel(time=date)['icec'].values
        sstlat = sst_file['lat']
        sstlon = sst_file['lon']

    fvlat = in_ds['yc'].values
    fvlon = in_ds['xc'].values

    if datasource == "GDAS":
        sst_fv = linint2(sstlon, sstlat[::-1], sst_gfs[::-1], True, fvlon, fvlat)
        ice_fv = linint2(sstlon, sstlat[::-1], ice_gfs[::-1], True, fvlon, fvlat)
    elif datasource == "NOAAOI":
        sst_fv = linint2(sstlon, sstlat, sst_gfs, True, fvlon, fvlat)
        ice_fv = linint2(sstlon, sstlat, ice_gfs, True, fvlon, fvlat)

    # Handling missing values
    sst_fv = np.where(np.isnan(sst_fv), TTHRESH, sst_fv)
    sst_fv = np.where(sst_fv > 500, TTHRESH, sst_fv)

    # SST to double precision
    sst_fv_dbl = sst_fv.astype(np.float64)

    # Correct SST time records
    sst_fv_dbl_time = np.empty((12, fvlat.shape[0], fvlon.shape[0]), dtype=np.float64)

    if do_anom:
        in_climo = xr.open_dataset("sst_1x1.nc")
        climo_SST = in_climo['SST_cpl'].values

        doy = (datetime(yyyy, mm, dd) - datetime(yyyy, 1, 1)).days + 1
        sst_anom = sst_fv_dbl - climo_SST[doy-1]

        for i in range(12):
            sst_fv_dbl_time[i] = climo_SST[i] + sst_anom
    else:
        for i in range(12):
            sst_fv_dbl_time[i] = sst_fv_dbl

    # Handling ice data
    ice_fv = np.where(np.isnan(ice_fv), 0.0, ice_fv)
    ice_fv = np.where(ice_fv > 500, 1.0, ice_fv)

    ice_fv_dbl = ice_fv.astype(np.float64)
    ice_fv_dbl_time = np.empty_like(sst_fv_dbl_time)

    if do_anom:
        climo_ice = in_climo['ice_cov'].values
        ice_anom = ice_fv_dbl - climo_ice[doy-1]

        for i in range(12):
            ice_fv_dbl_time[i] = climo_ice[i] + ice_anom
    else:
        for i in range(12):
            ice_fv_dbl_time[i] = ice_fv_dbl

    ice_fv_dbl_time = np.clip(ice_fv_dbl_time, 0.0, 1.0)

    # Time and date setup
    nc_time = np.array([15.5, 45, 74.5, 105, 135.5, 166, 196.5, 227.5, 258, 288.5, 319, 349.5], dtype=np.float64)
    nc_date = np.array([116, 215, 316, 416, 516, 616, 716, 816, 916, 1016, 1116, 1216], dtype=np.int32)
    nc_datesec = np.array([43200, 0, 43200, 0, 43200, 0, 43200, 43200, 0, 43200, 0, 43200], dtype=np.int32)

    pyfuncs.print_all_variables_info()

    # Writing to NetCDF
    out_ds = xr.Dataset(
        {
            'SST_cpl': (['time', 'lat', 'lon'], sst_fv_dbl_time),
            'ice_cov': (['time', 'lat', 'lon'], ice_fv_dbl_time),
            'sst_snapshot': (['lat', 'lon'], sst_fv_dbl),
            'ice_snapshot': (['lat', 'lon'], ice_fv_dbl),
            'date': (['time'], nc_date),
            'datesec': (['time'], nc_datesec)
        },
        coords={
            'time': nc_time,
            'lat': fvlat,
            'lon': fvlon
        }
    )

    out_ds.to_netcdf(SST_write_file)
    logging.info(f"Data written to {SST_write_file}")

def parse_args_sst():
    parser = argparse.ArgumentParser(description="Process SST and Ice data for CESM F compsets.")

    # Required arguments
    parser.add_argument('--initdate', type=str, required=True,
                        help='Initialization date in format YYYYMMDDHH')
    parser.add_argument('--predict_docn', type=int, choices=[0, 1], required=True,
                        help='Flag for predict_docn: 0 (false) or 1 (true)')
    parser.add_argument('--inputres', type=str, required=True,
                        help='Input resolution (e.g., "180x360")')
    parser.add_argument('--datasource', type=str, required=True,
                        help='Data source (e.g., "GDAS", "NOAAOI")')
    parser.add_argument('--sstDataFile', type=str, required=True,
                        help='Full path to the SST data file')
    parser.add_argument('--iceDataFile', type=str, required=True,
                        help='Full path to the Ice data file')
    parser.add_argument('--SST_write_file', type=str, required=True,
                        help='Full path to the output SST file')

    # Optional arguments with defaults
    parser.add_argument('--smooth_ice', action='store_true', default=False,
                        help='If set, smooth the ice field (default: False)')
    parser.add_argument('--smooth_iter', type=int, default=3,
                        help='Number of iterations for smoothing the ice field (default: 3)')
    parser.add_argument('--TTHRESH', type=float, default=271.9,
                        help='Initial threshold for ice vs. open ocean (default: 271.9)')
    parser.add_argument('--KtoC', type=float, default=273.15,
                        help='Kelvin to Celsius conversion constant (default: 273.15)')
    parser.add_argument('-v', '--verbose', action='store_true',
                        help='Increase output verbosity for debugging')

    return parser.parse_args()

def main():
    args = parse_args_sst()

    pyfuncs.configure_logging(args.verbose)

    # Pass the parsed arguments to your processing function
    process_sst_ice(
        args.initdate,
        args.predict_docn,
        args.inputres,
        args.datasource,
        args.sstDataFile,
        args.iceDataFile,
        args.SST_write_file
    )

if __name__ == "__main__":
    main()