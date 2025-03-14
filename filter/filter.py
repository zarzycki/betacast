#!/usr/bin/env python3
"""
Filter script for Betacast init.
"""

import sys
import numpy as np
import xarray as xr
import argparse
import logging

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    datefmt='%Y-%m-%d %H:%M:%S',
)
logger = logging.getLogger(__name__)

# Define fill value constant to match NCL
NC_FLOAT_FILL = 9.96921e+36

def parse_args():
    """
    Parse command line arguments.

    Returns:
        argparse.Namespace: Parsed command line arguments
    """
    parser = argparse.ArgumentParser(description='Filter netCDF variables')
    parser.add_argument('--endhour', type=float, required=True,
                        help='Simulation length in hours')
    parser.add_argument('--tcut', type=float, required=True,
                        help='Filter cutoff time in hours')
    parser.add_argument('--filtfile_name', type=str, required=True,
                        help='Input file name')
    parser.add_argument('--writefile_name', type=str, required=True,
                        help='Output file name')
    return parser.parse_args()

def calculate_filter_coefficients(numtime, endhour, tcut):
    """
    Calculate filter coefficients h_k based on the given parameters.

    Args:
        numtime: Number of time steps
        endhour: Simulation length in hours
        tcut: Filter cutoff time in hours

    Returns:
        np.ndarray: Filter coefficients h_k
    """
    pi = np.pi
    h_k = np.zeros(numtime, dtype=np.float64)
    omega = np.zeros(numtime, dtype=np.float64)
    M = numtime - 1
    N = M // 2
    deltat = float(endhour) / float(M)  # Timestep (hours)
    theta_c = 2 * pi * deltat / tcut
    eps = 0

    logger.info("-----------------")
    logger.info(f"Using a t_cut of: {tcut}")
    logger.info(f"Simulation is {endhour} hours long")
    logger.info(f"timestep is: {deltat}")

    # Get range of k values
    k = np.arange(-N, N+1)

    # Calculate filter coefficients
    for i in range(len(k)):
        if k[i] != 0:
            # Formula for non-zero k values
            omega[i] = np.sin(((eps + k[i]) * pi) / (N + 1)) / (((eps + k[i]) * pi) / (N + 1) + eps)
            h_k[i] = omega[i] * np.sin((eps + k[i]) * theta_c) / ((eps + k[i]) * pi)
        else:
            # Special case for k=0, using limit approximation
            k_plus = 0.0001
            k_minus = -0.0001
            omega_plus = np.sin((k_plus * pi) / (N + 1)) / ((k_plus * pi) / (N + 1))
            h_k_plus = omega_plus * np.sin(k_plus * theta_c) / (k_plus * pi)
            omega_minus = np.sin((k_minus * pi) / (N + 1)) / ((k_minus * pi) / (N + 1))
            h_k_minus = omega_minus * np.sin(k_minus * theta_c) / (k_minus * pi)
            h_k[i] = (h_k_plus + h_k_minus) / 2

    # Normalize h_k so sum equals 1
    h_k = h_k * (1 / np.sum(h_k))
    logger.info(f"Filter coefficients: min={h_k.min():.6f}, max={h_k.max():.6f}, sum={np.sum(h_k):.6f}")

    return h_k

def apply_filter(x, h_k, thisVar, nlev, ncol):
    """
    Apply filter to variable data.

    Args:
        x: Variable data to filter
        h_k: Filter coefficients
        thisVar: Variable name
        nlev: Number of vertical levels
        ncol: Number of columns

    Returns:
        np.ndarray: Filtered data
    """
    # Get the input data type
    dtype = x.dtype

    if thisVar == "PS":
        # PS is a 2D variable (time, ncol)
        x_filt = np.zeros((1, ncol), dtype=dtype)

        # Apply filter to each column using vectorized operations
        for i in range(ncol):
            x_filt[0, i] = np.dot(x[:, i], h_k)

    else:
        # Other variables are 3D (time, lev, ncol)
        x_filt = np.zeros((1, nlev, ncol), dtype=dtype)

        # Apply filter to each level and column
        for i in range(ncol):
            for j in range(nlev):
                x_filt[0, j, i] = np.dot(x[:, j, i], h_k)

    return x_filt

def get_dimensions(filtfile):
    """
    Extract dimensions from the input dataset.

    Args:
        filtfile: Input dataset

    Returns:
        tuple: Number of levels, columns, time values, and time count

    Raises:
        ValueError: If 'T' variable is not found in the input file
    """
    if 'T' in filtfile:
        dimT = filtfile['T'].shape
        nlev = dimT[1]
        ncol = dimT[2]
    else:
        raise ValueError("T variable not found in input file")

    time = filtfile['time'].values
    numtime = len(time)

    logger.debug(f"Dimensions - nlev: {nlev}, ncol: {ncol}, numtime: {numtime}")

    return nlev, ncol, time, numtime

# Parse command line arguments
args = parse_args()

# Extract parameters
endhour = args.endhour
tcut = args.tcut
filtfile_name = args.filtfile_name
writefile_name = args.writefile_name

# Print filter information
logger.info("FILTER PROCESSING")
logger.info(f"Using filter file: {filtfile_name}")
logger.info(f"Will output post-filter file here: {writefile_name}")

# Open input file
logger.info("Opening input file...")
filtfile = xr.open_dataset(filtfile_name)

# Read the output file - this file should already exist (created by cp command)
output_ds = xr.open_dataset(writefile_name)

# Fix lat/lon arrays if they exist in the filtfile but are zeros in output_ds
if 'lat' in filtfile and 'lat' in output_ds:
    if np.all(output_ds.lat == 0) and not np.all(filtfile.lat == 0):
        logger.info("Updating lat values from input file")
        output_ds['lat'] = filtfile['lat']

if 'lon' in filtfile and 'lon' in output_ds:
    if np.all(output_ds.lon == 0) and not np.all(filtfile.lon == 0):
        logger.info("Updating lon values from input file")
        output_ds['lon'] = filtfile['lon']

# Define variables to filter
vars = ["PS", "T", "U", "V", "Q", "CLDLIQ", "CLDICE"]

# Get dimensions
nlev, ncol, time, numtime = get_dimensions(filtfile)
logger.info(f"Input dimensions: nlev={nlev}, ncol={ncol}, numtime={numtime}")

# Calculate filter coefficients
h_k = calculate_filter_coefficients(numtime, endhour, tcut)

# Prepare encoding dictionary for each variable
encoding = {}

# Process each variable
for z in range(len(vars)):
    thisVar = vars[z]
    logger.info(f"Processing {thisVar}")

    if thisVar in filtfile:
        # Get variable data
        x = filtfile[thisVar].values

        # Apply filter
        x_filt = apply_filter(x, h_k, thisVar, nlev, ncol)

        logger.info(f"Writing {thisVar}...")

        # Create the appropriate dimensions for this variable
        if thisVar == "PS":
            # For PS, which is 2D
            output_ds[thisVar].values = x_filt
        else:
            # For 3D variables
            output_ds[thisVar].values = x_filt

        # Set the encoding for this variable
        var_dtype = output_ds[thisVar].dtype
        encoding[thisVar] = {'_FillValue': np.array(NC_FLOAT_FILL, dtype=var_dtype)}

# Add encoding for lat/lon if they exist
if 'lat' in output_ds:
    encoding['lat'] = {'_FillValue': -900.0}
if 'lon' in output_ds:
    encoding['lon'] = {'_FillValue': -900.0}

# Write the modified dataset to disk
logger.info(f"Writing output to {writefile_name}")
output_ds.to_netcdf(path=writefile_name, mode='w', encoding=encoding)
logger.info("done")