import numpy as np
import xarray as xr
from datetime import datetime
import os
import subprocess
import logging
import shutil
from scipy.spatial import SphericalVoronoi
import matplotlib.pyplot as plt
import warnings
try:
    from numba import jit, prange
    HAS_NUMBA = True
    print("hooray numba!")
except ImportError:
    HAS_NUMBA = False
    print("Warning: numba not found, falling back to slower implementation")


def mirrorP2P(p1, po):
    """Mirrors point p1 with respect to po."""
    dVec = p1 - po
    return po - dVec


def get_att_value(opt, attname, default_val):
    """
    This function checks if the given attribute (attname) exists in the 'opt' dictionary.
    If it exists and is not missing, it returns its value. Otherwise, it returns the
    provided default value.

    Parameters:
    - opt: A dictionary containing various options and their values.
    - attname: The name of the attribute to check for.
    - default_val: The default value to return if the attribute is not found or is missing.

    Returns:
    - The value of the attribute if it exists, otherwise the default value.
    """
    if isinstance(opt, dict) and attname in opt:
        return opt.get(attname, default_val)
    else:
        return default_val


def isatt_logical_true(opt, att):
    return opt.get(att, False) is True


def isatt_logical_false(opt, att):
    return opt.get(att, False) is False


def check_grid_type(grid_type):
    # Check for rectilinear, curvilinear, or unstructured grid types
    if grid_type in ["rectilinear", "curvilinear", "unstructured"]:
        return grid_type

    # Check for "NxM" format like "1x1", "2x3", "0.25x0.25"
    if "x" in grid_type:
        str_split = grid_type.split("x")
        if len(str_split) != 2:
            print("check_grid_type: invalid format for NxM type of grid.")
            return "none"
        try:
            dlat = float(str_split[0])
            dlon = float(str_split[1])
        except ValueError:
            print("check_grid_type: invalid format for NxM type of grid.")
            return "none"
        return {"type": "degree", "dlat": dlat, "dlon": dlon}

    # Check for "Ndeg" format like "1deg", "0.25deg"
    if "deg" in grid_type:
        str_split = grid_type.split("deg")
        if len(str_split) != 1:
            print("check_grid_type: invalid format for Ndeg type of grid.")
            return "none"
        try:
            dlat = float(str_split[0])
        except ValueError:
            print("check_grid_type: invalid format for Ndeg type of grid.")
            return "none"
        return {"type": "degree", "dlat": dlat, "dlon": dlat}

    # Check for "G64" format for Gaussian grids
    if "G" in grid_type:
        str_split = grid_type.split("G")
        if len(str_split) != 1:
            print("check_grid_type: invalid format for gaussian grid.")
            return "none"
        try:
            nlon = int(str_split[0])
        except ValueError:
            print("check_grid_type: invalid format for gaussian grid.")
            return "none"
        if nlon % 2 != 0:
            print("check_grid_type: invalid number for gaussian grid.")
            return "none"
        return {"type": "gaussian", "nlon": nlon, "nlat": nlon // 2}

    print("check_grid_type: invalid grid type")
    return "none"


def isLatorLon(fid, var_names):
    """
    This function checks if the 'units' attribute of a set of variable names
    contains 'north' or 'east'. It returns 'latitude' if 'north' is found,
    'longitude' if 'east' is found, and 'unknown' otherwise.

    Parameters:
    - fid: xarray.Dataset or netCDF4.Dataset
    - var_names: List of variable names (strings) to check

    Returns:
    - List of strings indicating 'latitude', 'longitude', or 'unknown'
    """
    output = ["unknown"] * len(var_names)

    for i, var_name in enumerate(var_names):
        if var_name in fid.variables:
            var = fid[var_name]
            if "units" in var.attrs:
                units = var.attrs["units"].lower()
                if "north" in units:
                    output[i] = "latitude"
                elif "east" in units:
                    output[i] = "longitude"

    return output


def check_both_atts(opt, att1, att2, attval):
    """
    This function checks if two attributes are set in the 'opt' dictionary and
    if they are set to different values. If they are different, the function sets
    both attributes to the provided value 'attval' and prints a warning.

    Parameters:
    - opt: A dictionary containing various options and their values.
    - att1: The first attribute to check.
    - att2: The second attribute to check.
    - attval: The value to set both attributes to if they differ.

    Returns:
    - None
    """
    if att1 in opt and att2 in opt and opt[att1] != opt[att2]:
        print(f"check_both_atts: you have set {att1} and {att2} to different values.")
        print(f"      Setting them both to {attval} to be safe.")
        opt[att1] = attval
        opt[att2] = attval


def totypeof(inVar, outType):
    """
    Coerces the type of a variable to the given type.

    Parameters:
    - inVar: The input variable to be converted.
    - outType: The target data type as a string.

    Returns:
    - outVar: The variable converted to the specified type.
    """
    outVar = None

    if outType == "double":
        outVar = np.array(inVar, dtype=np.float64)
    elif outType == "float":
        outVar = np.array(inVar, dtype=np.float32)
    elif outType == "integer":
        outVar = np.array(inVar, dtype=np.int32)
    elif outType == "int64":
        outVar = np.array(inVar, dtype=np.int64)
    elif outType == "uint64":
        outVar = np.array(inVar, dtype=np.uint64)
    elif outType == "long":  # In Python, long and int are the same on most systems
        outVar = np.array(inVar, dtype=np.int64)
    elif outType == "ulong":
        outVar = np.array(inVar, dtype=np.uint64)
    elif outType == "uint":
        outVar = np.array(inVar, dtype=np.uint32)
    elif outType == "short":
        outVar = np.array(inVar, dtype=np.int16)
    elif outType == "ushort":
        outVar = np.array(inVar, dtype=np.uint16)
    elif outType == "byte":
        outVar = np.array(inVar, dtype=np.int8)
    elif outType == "ubyte":
        outVar = np.array(inVar, dtype=np.uint8)
    elif outType == "string":
        outVar = np.array(inVar, dtype=str)
    elif outType == "character":  # This would typically convert to a numpy char array
        outVar = np.array(inVar, dtype=np.string_)
    else:
        raise ValueError("totypeof: Error: conversion to the provided type is not supported.")

    return outVar


def get_mask_name(opt):
    """
    This function returns the grid mask name, if any, being provided by the user.
    It handles both the old case of "Mask2D" being used (pre NCL V6.3.0),
    and the newer "GridMask" being used.

    Parameters:
    - opt: A dictionary containing options, potentially including mask names.

    Returns:
    - A string representing the name of the grid mask or an empty string if none is found.
    """
    if not opt:
        return ""

    valid_mask_names = ["GridMask", "Mask2D"]

    for mask_name in valid_mask_names:
        if mask_name in opt:
            return mask_name

    return ""


def rotate_latlon(lat, lon, Rot):
    """
    Rotates the given lat/lon grid by the specified rotation angle.
    This operation modifies the input latitude and longitude arrays directly.

    Parameters:
    - lat: 2D numpy array of latitudes.
    - lon: 2D numpy array of longitudes.
    - Rot: Rotation angle in degrees.

    Returns:
    - (lat_rotated, lon_rotated): Tuple of rotated latitudes and longitudes.
    """
    if lat.shape != lon.shape:
        raise ValueError("rotate_latlon: lat and lon must have the same dimensions.")

    # Convert degrees to radians
    d2r = np.pi / 180.0

    # Create rotation matrix
    RotateMat = np.array([[np.cos(Rot * d2r), -np.sin(Rot * d2r)],
                          [np.sin(Rot * d2r),  np.cos(Rot * d2r)]])

    # Flatten lat and lon arrays for matrix multiplication
    tmpPoints = np.vstack((lat.flatten(), lon.flatten()))

    # Apply rotation
    rotated_points = RotateMat @ tmpPoints

    # Reshape back to original dimensions
    lat_rotated = rotated_points[0, :].reshape(lat.shape)
    lon_rotated = rotated_points[1, :].reshape(lon.shape)

    return lat_rotated, lon_rotated


def gaus(nlat):
    """
    Computes Gaussian latitudes and weights.

    Parameters:
    nlat : int
        Number of latitude points PER hemisphere.

    Returns:
    np.ndarray
        A 2D array of shape (2*nlat, 2) containing Gaussian latitudes (column 0)
        and Gaussian weights (column 1).
    """
    # Ensure nlat is an integer
    nlat = int(nlat)

    # Generate the points and weights using the Gauss-Legendre quadrature
    x, w = np.polynomial.legendre.leggauss(nlat)

    # Convert the x values to latitudes (in degrees)
    latitudes = np.arcsin(x) * (180.0 / np.pi)

    # Weights are already normalized, but we need to double them
    # for both hemispheres
    weights = w * 2.0

    # Combine the northern and southern hemispheres
    latitudes = np.concatenate([latitudes[::-1], -latitudes])
    weights = np.concatenate([weights[::-1], weights])

    # Return a 2D array where the first column is latitudes
    # and the second column is weights
    return np.column_stack((latitudes, weights))


def calc_SCRIP_corners_boundaries(lat2d, lon2d, grid_corner_lat, grid_corner_lon, Opt):
    DEBUG = isatt_logical_true(Opt, "Debug")
    nlat, nlon = lat2d.shape
    nlatp1 = nlat + 1
    nlonp1 = nlon + 1

    # Create arrays to hold 2D grid corner coordinates
    grid_corner_lat2d = np.zeros((nlatp1, nlonp1), dtype=lat2d.dtype)
    grid_corner_lon2d = np.zeros((nlatp1, nlonp1), dtype=lon2d.dtype)

    # Calculate interior locations of the corner grid
    for i in range(1, nlat):
        for j in range(1, nlon):
            grid_corner_lat2d[i, j] = np.mean([lat2d[i-1, j-1], lat2d[i, j-1], lat2d[i-1, j], lat2d[i, j]])
            grid_corner_lon2d[i, j] = np.mean([lon2d[i-1, j-1], lon2d[i, j-1], lon2d[i-1, j], lon2d[i, j]])

    # Calculate bottom locations of the corner grid
    for j in range(1, nlon):
        grid_corner_lat2d[0, j] = np.mean([lat2d[0, j], lat2d[0, j-1]])
        grid_corner_lon2d[0, j] = np.mean([lon2d[0, j], lon2d[0, j-1]])

    # Calculate top locations of the corner grid
    for j in range(1, nlon):
        grid_corner_lat2d[nlat, j] = np.mean([lat2d[nlat-1, j], lat2d[nlat-1, j-1]])
        grid_corner_lon2d[nlat, j] = np.mean([lon2d[nlat-1, j], lon2d[nlat-1, j-1]])

    # Calculate left locations of the corner grid
    for i in range(1, nlat):
        grid_corner_lat2d[i, 0] = np.mean([lat2d[i, 0], lat2d[i-1, 0]])
        grid_corner_lon2d[i, 0] = np.mean([lon2d[i, 0], lon2d[i-1, 0]])

    # Calculate right locations of the corner grid
    for i in range(1, nlat):
        grid_corner_lat2d[i, nlon] = np.mean([lat2d[i, nlon-1], lat2d[i-1, nlon-1]])
        grid_corner_lon2d[i, nlon] = np.mean([lon2d[i, nlon-1], lon2d[i-1, nlon-1]])

    # Set the four corners of the corner grid
    grid_corner_lat2d[0, 0] = lat2d[0, 0]
    grid_corner_lon2d[0, 0] = lon2d[0, 0]
    grid_corner_lat2d[0, nlon] = lat2d[0, nlon-1]
    grid_corner_lon2d[0, nlon] = lon2d[0, nlon-1]
    grid_corner_lat2d[nlat, nlon] = lat2d[nlat-1, nlon-1]
    grid_corner_lon2d[nlat, nlon] = lon2d[nlat-1, nlon-1]
    grid_corner_lat2d[nlat, 0] = lat2d[nlat-1, 0]
    grid_corner_lon2d[nlat, 0] = lon2d[nlat-1, 0]

    # Map the 2D grid corner lat/lon arrays to the 1D grid corner lat/lon arrays
    n = 0
    for i in range(nlat):
        for j in range(nlon):
            grid_corner_lat[n, 0] = grid_corner_lat2d[i, j]
            grid_corner_lat[n, 1] = grid_corner_lat2d[i, j+1]
            grid_corner_lat[n, 2] = grid_corner_lat2d[i+1, j+1]
            grid_corner_lat[n, 3] = grid_corner_lat2d[i+1, j]
            grid_corner_lon[n, 0] = grid_corner_lon2d[i, j]
            grid_corner_lon[n, 1] = grid_corner_lon2d[i, j+1]
            grid_corner_lon[n, 2] = grid_corner_lon2d[i+1, j+1]
            grid_corner_lon[n, 3] = grid_corner_lon2d[i+1, j]
            n += 1

    if DEBUG:
        print("calc_SCRIP_corners_boundaries")
        print(f"min/max grid_corner_lat: {grid_corner_lat.min()} / {grid_corner_lat.max()}")
        print(f"min/max grid_corner_lon: {grid_corner_lon.min()} / {grid_corner_lon.max()}")


def calc_SCRIP_corners_noboundaries(lat2d, lon2d, grid_corner_lat, grid_corner_lon, Opt):
    DEBUG = isatt_logical_true(Opt, "Debug")
    delta = Opt.get("GridCornerDelta", 0.25)

    if lat2d.shape != lon2d.shape:
        raise ValueError("Latitude and longitude must have the same dimensions.")

    nlat, nlon = lat2d.shape
    nlatp1 = nlat + 1
    nlonp1 = nlon + 1
    nlatp2 = nlat + 2
    nlonp2 = nlon + 2

    if DEBUG:
        print("calc_SCRIP_corners_noboundaries")
        print(f"     min/max original lat: {lat2d.min()} / {lat2d.max()}")
        print(f"     min/max original lon: {lon2d.min()} / {lon2d.max()}")

    # Extend the lat/lon grid (needed to calculate the corners at the boundaries).
    Extlat2d = np.zeros((nlatp2, nlonp2), dtype=lat2d.dtype)
    Extlon2d = np.zeros((nlatp2, nlonp2), dtype=lat2d.dtype)

    # The middle grid is exactly the original lat2d/lon2d arrays
    Extlat2d[1:nlat+1, 1:nlon+1] = lat2d
    Extlon2d[1:nlat+1, 1:nlon+1] = lon2d

    # Bottom row, minus corners
    Extlat2d[0, 1:nlon+1] = mirrorP2P(lat2d[1, :], lat2d[0, :])
    Extlon2d[0, 1:nlon+1] = mirrorP2P(lon2d[1, :], lon2d[0, :])

    # Top row, minus corners
    Extlat2d[nlatp1, 1:nlon+1] = mirrorP2P(lat2d[nlat-2, :], lat2d[nlat-1, :])
    Extlon2d[nlatp1, 1:nlon+1] = mirrorP2P(lon2d[nlat-2, :], lon2d[nlat-1, :])

    # Left column, minus corners
    Extlat2d[1:nlat+1, 0] = mirrorP2P(lat2d[:, 1], lat2d[:, 0])
    Extlon2d[1:nlat+1, 0] = mirrorP2P(lon2d[:, 1], lon2d[:, 0])

    # Right column, minus corners
    Extlat2d[1:nlat+1, nlonp1] = mirrorP2P(lat2d[:, nlon-2], lat2d[:, nlon-1])
    Extlon2d[1:nlat+1, nlonp1] = mirrorP2P(lon2d[:, nlon-2], lon2d[:, nlon-1])

    # Lower left corner
    Extlat2d[0, 0] = mirrorP2P(lat2d[1, 1], lat2d[0, 0])
    Extlon2d[0, 0] = mirrorP2P(lon2d[1, 1], lon2d[0, 0])

    # Upper right corner
    Extlat2d[nlatp1, nlonp1] = mirrorP2P(lat2d[nlat-2, nlon-2], lat2d[nlat-1, nlon-1])
    Extlon2d[nlatp1, nlonp1] = mirrorP2P(lon2d[nlat-2, nlon-2], lon2d[nlat-1, nlon-1])

    # Lower right corner
    Extlat2d[0, nlonp1] = mirrorP2P(lat2d[1, nlon-2], lat2d[0, nlon-1])
    Extlon2d[0, nlonp1] = mirrorP2P(lon2d[1, nlon-2], lon2d[0, nlon-1])

    # Upper left corner
    Extlat2d[nlatp1, 0] = mirrorP2P(lat2d[nlat-2, 1], lat2d[nlat-1, 0])
    Extlon2d[nlatp1, 0] = mirrorP2P(lon2d[nlat-2, 1], lon2d[nlat-1, 0])

    if DEBUG:
        print("calc_SCRIP_corners_noboundaries")
        print(f"     min/max Extlat2d: {Extlat2d.min()} / {Extlat2d.max()}")
        print(f"     min/max Extlon2d: {Extlon2d.min()} / {Extlon2d.max()}")

    # Calculate the cell center of the extended grid, which would be the corner coordinates for the original grid.
    tmp_lat = Extlat2d[:, 1:nlonp1+1] + Extlat2d[:, 0:nlon+1]
    ExtGridCenter_lat = ((tmp_lat[1:nlatp1+1, :] + tmp_lat[0:nlat+1, :]) * delta).flatten(order="C")
    del tmp_lat

    tmp_lon = Extlon2d[:, 1:nlonp1+1] + Extlon2d[:, 0:nlon+1]
    ExtGridCenter_lon = ((tmp_lon[1:nlatp1+1, :] + tmp_lon[0:nlat+1, :]) * delta).flatten(order="C")
    del tmp_lon

    if DEBUG:
        print("calc_SCRIP_corners_noboundaries")
        print(f"     min/max ExtGridCenter_lat: {ExtGridCenter_lat.min()} / {ExtGridCenter_lat.max()}")
        print(f"     min/max ExtGridCenter_lon: {ExtGridCenter_lon.min()} / {ExtGridCenter_lon.max()}")

    # Extract the grid cell corners
    ii, jj = np.meshgrid(np.arange(nlon), np.arange(nlat))

    ii = ii.flatten(order="C")
    jj = jj.flatten(order="C")

    grid_corner_lat[:, 0] = ExtGridCenter_lat[jj * nlonp1 + ii]
    grid_corner_lat[:, 1] = ExtGridCenter_lat[jj * nlonp1 + (ii + 1)]
    grid_corner_lat[:, 2] = ExtGridCenter_lat[(jj + 1) * nlonp1 + (ii + 1)]
    grid_corner_lat[:, 3] = ExtGridCenter_lat[(jj + 1) * nlonp1 + ii]

    grid_corner_lon[:, 0] = ExtGridCenter_lon[jj * nlonp1 + ii]
    grid_corner_lon[:, 1] = ExtGridCenter_lon[jj * nlonp1 + (ii + 1)]
    grid_corner_lon[:, 2] = ExtGridCenter_lon[(jj + 1) * nlonp1 + (ii + 1)]
    grid_corner_lon[:, 3] = ExtGridCenter_lon[(jj + 1) * nlonp1 + ii]

    if DEBUG:
        print("calc_SCRIP_corners_noboundaries")
        print(f"     min/max grid_corner_lat: {grid_corner_lat.min()} / {grid_corner_lat.max()}")
        print(f"     min/max grid_corner_lon: {grid_corner_lon.min()} / {grid_corner_lon.max()}")


def curvilinear_to_SCRIP(FName, lat2d, lon2d, Opt):
    # Check for options
    PrintTimings = True
    if PrintTimings:
        start_time = datetime.now()

    # Determine the file format
    if Opt.get("NetCDFType", "").lower() == "netcdf4":
        nc_file_type = "NETCDF4"
    else:
        nc_file_type = "NETCDF3_CLASSIC"

    if lat2d.shape != lon2d.shape:
        raise ValueError("Latitude and longitude must have the same dimensions.")

    nlat, nlon = lat2d.shape
    grid_size = nlat * nlon
    grid_corners = 4

    FTitle = Opt.get("Title", f"curvilinear_to_SCRIP ({nlat},{nlon})")

    # Handle grid mask
    grid_mask_name = Opt.get("MaskName", "")
    if grid_mask_name:
        grid_mask = Opt.get(grid_mask_name)
        if grid_mask.shape != lat2d.shape:
            raise ValueError(f"Opt@{grid_mask_name} is not the correct dimensionality.")
    else:
        grid_mask = np.ones(lat2d.shape, dtype=np.int32)  # Ensure int32

    # Flatten latitude and longitude arrays
    lat_flat = lat2d.flatten()
    lon_flat = lon2d.flatten()
    grid_size = len(lat_flat)

    # Initialize grid corner arrays
    grid_corner_lat = np.zeros((grid_size, grid_corners), dtype=lat2d.dtype)
    grid_corner_lon = np.zeros((grid_size, grid_corners), dtype=lon2d.dtype)

    # Get the grid lat/lon corners
    if "GridCornerLat" in Opt and "GridCornerLon" in Opt:
        print("curvilinear_to_SCRIP: using grid corners provided by user...")
        GridCornerLat = np.array(Opt["GridCornerLat"], dtype=lat2d.dtype)
        GridCornerLon = np.array(Opt["GridCornerLon"], dtype=lon2d.dtype)
        grid_corner_lat = GridCornerLat.reshape((grid_size, grid_corners))
        grid_corner_lon = GridCornerLon.reshape((grid_size, grid_corners))
    else:
        # Estimate the grid cell corners; it's better if the user provides them.
        print("curvilinear_to_SCRIP: calculating grid corners...")

        if np.any(np.abs(lat2d) == 90):
            print("curvilinear_to_SCRIP: one or more lat values are at the poles, so calculating grid corners using calc_SCRIP_corners_boundaries...")
            calc_SCRIP_corners_boundaries(lat2d, lon2d, grid_corner_lat, grid_corner_lon, Opt)
        else:
            print("curvilinear_to_SCRIP: no lat values are at the poles, so calculating grid corners using calc_SCRIP_corners_noboundaries...")
            calc_SCRIP_corners_noboundaries(lat2d, lon2d, grid_corner_lat, grid_corner_lon, Opt)
            if np.any(np.abs(grid_corner_lat) > 90):
                print("curvilinear_to_SCRIP: calc_SCRIP_corners_noboundaries produced out-of-range latitude values. Trying calc_SCRIP_corners_boundaries...")
                calc_SCRIP_corners_boundaries(lat2d, lon2d, grid_corner_lat, grid_corner_lon, Opt)

    # Create NetCDF file
    with xr.Dataset() as ds:
        # Define dimensions
        ds = ds.assign_coords({
            'grid_size': np.arange(grid_size, dtype=np.int32),  # Ensure int32
            'grid_corners': np.arange(4, dtype=np.int32),       # Ensure int32
            'grid_rank': np.arange(2, dtype=np.int32)           # Ensure int32
        })

        # Define variables
        ds['grid_dims'] = xr.DataArray(np.array([nlon, nlat], dtype=np.int32), dims=["grid_rank"])
        ds['grid_center_lat'] = xr.DataArray(lat_flat, dims=["grid_size"])
        ds['grid_center_lon'] = xr.DataArray(lon_flat, dims=["grid_size"])
        ds['grid_imask'] = xr.DataArray(grid_mask.flatten().astype(np.int32), dims=["grid_size"])  # Ensure int32
        ds['grid_corner_lat'] = xr.DataArray(grid_corner_lat, dims=["grid_size", "grid_corners"])
        ds['grid_corner_lon'] = xr.DataArray(grid_corner_lon, dims=["grid_size", "grid_corners"])

        # Define the variable attributes
        ds['grid_center_lat'].attrs['units'] = "degrees"
        ds['grid_center_lon'].attrs['units'] = "degrees"
        ds['grid_imask'].attrs['units'] = "unitless"
        ds['grid_corner_lat'].attrs['units'] = "degrees"
        ds['grid_corner_lon'].attrs['units'] = "degrees"

        # Define the global attributes
        ds.attrs['date_created'] = datetime.now().strftime("%a %b %d %H:%M:%S %Z %Y")
        ds.attrs['Createdby'] = "ESMF_regridding.ncl"
        ds.attrs['Conventions'] = "SCRIP"
        ds.attrs['title'] = FTitle

        existing_vars = [var for var in ds.variables]
        encoding = {var: {'_FillValue': None} for var in existing_vars}

        # Save to NetCDF
        ds.to_netcdf(FName, format=nc_file_type, encoding=encoding)

    if PrintTimings:
        print(f"curvilinear_to_SCRIP took {datetime.now() - start_time}")


def rectilinear_to_SCRIP(FName, lat, lon, Opt):

    # Check for options
    PrintTimings = isatt_logical_true(Opt, "PrintTimings")
    if PrintTimings:
        start_time = datetime.now()

    nlat = lat.size
    nlon = lon.size

    # Conform the lat/lon to 2D arrays using meshgrid
    grid_center_lon, grid_center_lat = np.meshgrid(lon, lat)

    print(grid_center_lon.shape)
    print(grid_center_lat.shape)

    # Convert to SCRIP
    curvilinear_to_SCRIP(FName, grid_center_lat, grid_center_lon, Opt)

    if PrintTimings:
        print(f"rectilinear_to_SCRIP took {datetime.now() - start_time}")


def latlon_to_SCRIP(FName, GridType, Opt):
    PrintTimings = Opt.get("PrintTimings", False)
    if PrintTimings:
        start_time = datetime.now()

    # Check for "1x1", "2x3", "0.25x0.25" format
    grid_type = check_grid_type(GridType)
    if grid_type['type'] == "none":
        raise ValueError("latlon_to_SCRIP: invalid format for grid type")

    if grid_type['type'] == "degree":
        dlat = grid_type['dlat']
        dlon = grid_type['dlon']
    elif grid_type['type'] == "gaussian":
        nlat = grid_type['nlat']
        nlon = grid_type['nlon']
    else:
        raise ValueError("latlon_to_SCRIP: unsupported grid type")

    # Check for Opt attributes
    FTitle = Opt.get("Title", f"{GridType} grid")
    Rot = Opt.get("Rotation", 0.0)

    LLCorner = np.array(Opt.get("LLCorner", [-90, -180]))
    URCorner = np.array(Opt.get("URCorner", [90, 180]))

    # Check the LLCorner and URCorner variables
    if abs(LLCorner[0]) > 90 or abs(URCorner[0]) > 90:
        raise ValueError("latlon_to_SCRIP: lat corners must be within -90 to 90")

    BoxDiag = URCorner - LLCorner

    if Rot != 0.0:
        print("latlon_to_SCRIP: The box is rotated. Spacing may be different than what you meant.")
        rotate_latlon(BoxDiag[0], BoxDiag[1], -1.0 * Rot)

    if BoxDiag[0] <= 0.0:
        raise ValueError("latlon_to_SCRIP: The corner latitudes do not meet the criteria of a bounding box.")
    if BoxDiag[1] <= 0.0:
        raise ValueError("latlon_to_SCRIP: The corner longitudes do not meet the criteria of a bounding box.")

    if grid_type['type'] == "degree":
        nlat = int(BoxDiag[0] / dlat + 1)
        nlon = int(BoxDiag[1] / dlon + 1)
        lat = np.linspace(0,BoxDiag[0], nlat)
        lon = np.linspace(0,BoxDiag[1], nlon)
    elif grid_type['type'] == "gaussian":
        gau_info = gaus(nlat // 2)
        lat = gau_info[:, 0]  # Gaussian latitudes
        lon = np.linspace(-180, 180, nlon, endpoint=False)

    # Generate the cell center Lat/Lon
    grid_center_lat = np.tile(lat[:, np.newaxis], (1, nlon))
    grid_center_lon = np.tile(lon[np.newaxis, :], (nlat, 1))

    if get_mask_name(Opt) == "":
        Opt["GridMask"] = np.ones((nlat, nlon))

    # Rotate back the box, if needed
    if Rot != 0.0:
        rotate_latlon(grid_center_lat, grid_center_lon, Rot)

    # Translate the grid
    if grid_type['type'] != "gaussian":
        grid_center_lat += LLCorner[0]
        grid_center_lon += LLCorner[1]

    Opt["Title"] = FTitle

    # Make sure curvilinear_to_SCRIP doesn't print timings
    set_pt = False
    if "PrintTimings" in Opt:
        set_pt = True
        orig_pt = Opt["PrintTimings"]
        Opt["PrintTimings"] = False

    curvilinear_to_SCRIP(FName, grid_center_lat, grid_center_lon, Opt)

    if set_pt:
        Opt["PrintTimings"] = orig_pt

    if PrintTimings:
        print(f"latlon_to_SCRIP took {datetime.now() - start_time}")


def check_SCRIP_dims(ds, verbose=False):
    """
    Checks if the required dimensions in the SCRIP file are present.

    Parameters:
    ds (xarray.Dataset): The dataset to check.
    verbose (bool): If True, print which dimension is missing.

    Returns:
    bool: True if all required dimensions are present, False otherwise.
    """
    scrip_dim_names = ["grid_rank", "grid_size", "grid_corners"]

    for dim in scrip_dim_names:
        if dim not in ds.dims:
            if verbose:
                print(f"check_SCRIP_dims: {dim} not defined.")
            return False
    return True


def check_SCRIP_vars(ds, verbose=False):
    """
    Checks if the required variables in the SCRIP file are present.

    Parameters:
    ds (xarray.Dataset): The dataset to check.
    verbose (bool): If True, print which variable is missing.

    Returns:
    bool: True if all required variables are present, False otherwise.
    """
    scrip_vars = [
        "grid_dims", "grid_center_lat", "grid_center_lon",
        "grid_imask", "grid_corner_lat", "grid_corner_lon"
    ]

    for var in scrip_vars:
        if var not in ds.variables:
            if verbose:
                print(f"check_SCRIP_vars: {var} not defined.")
            return False
    return True


def is_SCRIP(filename, verbose=False):
    """
    Checks if the input filename follows the SCRIP standard.

    Parameters:
    filename (str): The path to the file.
    verbose (bool): If True, print which checks failed.

    Returns:
    bool: True if the file follows the SCRIP standard, False otherwise.
    """
    try:
        ds = xr.open_dataset(filename)
    except FileNotFoundError:
        if verbose:
            print(f"is_SCRIP: '{filename}' doesn't exist.")
        return False

    return check_SCRIP_dims(ds, verbose) and check_SCRIP_vars(ds, verbose)


def esmf_regrid_gen_weights(srcGridFile, dstGridFile, wgtFile, opt):
    """
    Generate regrid weights using ESMF_RegridWeightGen.

    Parameters:
    -----------
    srcGridFile : str
        The name of the NetCDF source grid file (ESMF or SCRIP).

    dstGridFile : str
        The name of the NetCDF destination grid file (ESMF or SCRIP).

    wgtFile : str
        The name of the NetCDF weight file to create.

    opt : dict
        Dictionary containing optional attributes:
            - 'Debug' (bool): Turns on debug print statements (default=False).
            - 'SrcESMF' (bool): Whether the source grid file is an ESMF file (default=False).
            - 'DstESMF' (bool): Whether the destination grid file is an ESMF file (default=False).
            - 'SrcRegional' (bool): Whether the source grid file is regional (default=False).
            - 'DstRegional' (bool): Whether the destination grid file is regional (default=False).
            - 'IgnoreUnmappedPoints' (bool): Ignore unmapped points (default=True).
            - 'InterpMethod' (str): Interpolation method, e.g., "bilinear", "patch", "conserve" (default="bilinear").
            - 'Pole' (str): Handle poles, e.g., "all", "none" (default depends on 'InterpMethod').
            - 'RemoveSrcFile' (bool): Remove source grid file after weight generation (default=False).
            - 'RemoveDstFile' (bool): Remove destination grid file after weight generation (default=False).
            - 'RemoveFiles' (bool): Remove both source and destination grid files after weight generation (default=False).
            - 'LargeFile' (bool): Use --64bit_offset for writing large files.
            - 'NetCDFType' (str): NetCDF format type, e.g., "netcdf3", "netcdf4".
            - 'NormalizeWeights' (bool): Normalize weights on the destination grid (default=False).

    Returns:
    --------
    None
    """

    # Handle mutually exclusive attributes
    check_both_atts(opt, "RemoveFiles", "RemoveDstFile", False)
    check_both_atts(opt, "RemoveFiles", "RemoveSrcFile", False)
    check_both_atts(opt, "WgtForceOverwrite", "ForceOverwrite", False)
    check_both_atts(opt, "WgtOverwrite", "Overwrite", False)
    check_both_atts(opt, "WgtNetCDFType", "NetCDFType", "netcdf3")

    # Define the ESMF executable path
    esmf_exec = "ESMF_RegridWeightGen"
    esmf_bindir = os.getenv("ESMFBINDIR", "")
    path_to_esmf = os.path.join(esmf_bindir, esmf_exec) if esmf_bindir else esmf_exec

    # Check if ESMF_RegridWeightGen exists
    if not shutil.which(path_to_esmf):
        logging.error(f"ESMF_regrid_gen_weights: could not find {esmf_exec} executable.")
        return

    # Construct the ESMF command
    esmf_cmd = [path_to_esmf, "--source", srcGridFile, "--destination", dstGridFile, "--weight", wgtFile]

    # Handle optional attributes using get_att_value
    if get_att_value(opt, "SrcESMF", False):
        esmf_cmd.append("--src_loc corner")
    if get_att_value(opt, "DstESMF", False):
        esmf_cmd.append("--dst_loc corner")
    interp_method = get_att_value(opt, "InterpMethod", None)
    if interp_method:
        esmf_cmd.extend(["--method", interp_method])
    pole = get_att_value(opt, "Pole", None)
    if pole:
        esmf_cmd.extend(["--pole", pole])
    if get_att_value(opt, "SrcRegional", False):
        esmf_cmd.append("--src_regional")
    if get_att_value(opt, "DstRegional", False):
        esmf_cmd.append("--dst_regional")
    if get_att_value(opt, "SrcESMF", False):
        esmf_cmd.append("--src_type ESMF")
    if get_att_value(opt, "DstESMF", False):
        esmf_cmd.append("--dst_type ESMF")
    if get_att_value(opt, "IgnoreUnmappedPoints", True):
        esmf_cmd.append("--ignore_unmapped")
    if get_att_value(opt, "LargeFile", False):
        esmf_cmd.append("--64bit_offset")
    if get_att_value(opt, "NetCDFType", "").lower() == "netcdf4":
        esmf_cmd.append("--netcdf4")
    if get_att_value(opt, "NormalizeWeights", False):
        esmf_cmd.append("--norm_type fracarea")

    # Execute the command
    result = subprocess.run(esmf_cmd, capture_output=True, text=True)
    if result.returncode != 0:
        logging.error(f"ESMF_regrid_gen_weights: failed with output:\n{result.stderr}")
        return

    # Remove files if requested
    if opt.get("RemoveFiles", False) or opt.get("RemoveSrcFile", False):
        os.remove(srcGridFile)
    if opt.get("RemoveFiles", False) or opt.get("RemoveDstFile", False):
        os.remove(dstGridFile)

    # Log the success message
    logging.info("ESMF_regrid_gen_weights: Completed weight generation successfully.")


def lat_lon_to_xyz(lat, lon):
    """Convert latitude/longitude to 3D Cartesian coordinates on unit sphere."""
    lat_rad, lon_rad = np.radians(lat), np.radians(lon)
    cos_lat = np.cos(lat_rad)
    return np.column_stack((
        cos_lat * np.cos(lon_rad),
        cos_lat * np.sin(lon_rad),
        np.sin(lat_rad)
    ))

def xyz_to_lat_lon(xyz):
    """Convert 3D Cartesian coordinates to latitude/longitude."""
    xyz_normalized = xyz / np.linalg.norm(xyz, axis=1)[:, np.newaxis]
    lat = np.degrees(np.arcsin(xyz_normalized[:, 2]))
    lon = np.degrees(np.arctan2(xyz_normalized[:, 1], xyz_normalized[:, 0]))
    return lat, lon

def sort_vertices_spherical(vertices, center_point):
    """
    Sort vertices counterclockwise around a center point on a sphere.

    Parameters:
    -----------
    vertices : ndarray
        Array of shape (N, 3) containing vertex coordinates in 3D
    center_point : ndarray
        Array of shape (3,) containing the center point coordinates

    Returns:
    --------
    ndarray
        Indices that will sort vertices counterclockwise
    """
    # Normalize all points to ensure they're on unit sphere
    center_norm = center_point / np.linalg.norm(center_point)
    vertices_norm = vertices / np.linalg.norm(vertices, axis=1)[:, np.newaxis]

    # Get North pole vector (or different reference if too close to poles)
    if abs(np.dot(center_norm, [0, 0, 1])) > 0.99:
        ref = np.array([1, 0, 0])
    else:
        ref = np.array([0, 0, 1])

    # Get reference direction (tangent to sphere at center_point)
    ref_dir = ref - center_norm * np.dot(ref, center_norm)
    ref_dir = ref_dir / np.linalg.norm(ref_dir)

    # Get vector perpendicular to both center_norm and ref_dir
    perp_dir = np.cross(center_norm, ref_dir)

    # For each vertex, compute the tangent vector from center to vertex
    # (removing the radial component)
    tangent_vectors = vertices_norm - center_norm * \
                     np.dot(vertices_norm, center_norm)[:, np.newaxis]

    # Normalize tangent vectors
    norms = np.linalg.norm(tangent_vectors, axis=1)
    valid = norms > 1e-10
    tangent_vectors[valid] /= norms[valid, np.newaxis]

    # Compute angles using atan2
    # Project onto our reference directions
    x_proj = np.dot(tangent_vectors, ref_dir)
    y_proj = np.dot(tangent_vectors, perp_dir)
    angles = np.arctan2(y_proj, x_proj)

    # Handle numerical instability
    angles[~valid] = 0

    return np.argsort(angles)

def normalize_lon_0_360(lon):
    """
    Normalize longitude to [0,360) range.
    """
    return lon % 360

def find_spherical_corners(lat, lon, debug=False):
    """
    Find corners for each point using spherical Voronoi tessellation with
    guaranteed counterclockwise ordering.
    """
    # Normalize input longitudes first
    lon = normalize_lon_0_360(lon)

    # Convert to 3D coordinates
    points_3d = lat_lon_to_xyz(lat, lon)

    # Create spherical Voronoi tessellation
    sv = SphericalVoronoi(points_3d)

    grid_size = len(lat)

    # Determine maximum number of corners across all regions
    max_corners = max(len(region) for region in sv.regions)
    if debug:
        print(f"Maximum number of corners in any cell: {max_corners}")

    corner_lats = np.full((grid_size, max_corners), np.nan)
    corner_lons = np.full((grid_size, max_corners), np.nan)

    # Get regions and vertices
    regions = sv.regions
    vertices = sv.vertices

    if debug:
        print(f"Found {len(regions)} regions")
        print(f"Found {len(vertices)} vertices")

    # Process each point and its corresponding region
    for i in range(grid_size):
        if debug and i % 1000 == 0:
            print(f"Processing point {i}/{grid_size}")

        region = regions[i]
        if len(region) > 0:
            # Get vertices for this region
            region_vertices = vertices[region]

            # Sort vertices counterclockwise
            sorted_idx = sort_vertices_spherical(region_vertices, points_3d[i])
            sorted_vertices = region_vertices[sorted_idx]

            # Convert to lat/lon
            rlat, rlon = xyz_to_lat_lon(sorted_vertices)

            # Normalize longitudes to [0,360)
            rlon = normalize_lon_0_360(rlon)

            # Handle cases where corners cross the prime meridian
            if lon[i] < 30 or lon[i] > 330:
                mean_lon = np.mean(rlon)
                if abs(mean_lon - lon[i]) > 180:
                    if lon[i] < 30:
                        rlon[rlon > 180] -= 360
                    else:
                        rlon[rlon < 180] += 360
                rlon = normalize_lon_0_360(rlon)

            # Store all corners
            n_corners = len(rlat)
            corner_lats[i, :n_corners] = rlat
            corner_lons[i, :n_corners] = rlon

            # Fill remaining slots with the last corner
            if n_corners < max_corners:
                corner_lats[i, n_corners:] = rlat[-1]
                corner_lons[i, n_corners:] = rlon[-1]

    return corner_lats, corner_lons

if HAS_NUMBA:
    @jit(nopython=True)
    def calculate_single_area(points):
        """Calculate area for a single cell using numba acceleration."""
        if len(points) < 3:
            return 0.0

        area = 0.0
        point0 = points[0]

        for j in range(len(points)-1):
            point1 = points[j]
            point2 = points[j+1]

            # Calculate cross product
            cross_x = point1[1]*point2[2] - point1[2]*point2[1]
            cross_y = point1[2]*point2[0] - point1[0]*point2[2]
            cross_z = point1[0]*point2[1] - point1[1]*point2[0]

            numerator = abs(point0[0]*cross_x + point0[1]*cross_y + point0[2]*cross_z)
            denominator = 1.0 + (point0[0]*point1[0] + point0[1]*point1[1] + point0[2]*point1[2]) + \
                              (point0[0]*point2[0] + point0[1]*point2[1] + point0[2]*point2[2]) + \
                              (point1[0]*point2[0] + point1[1]*point2[1] + point1[2]*point2[2])

            area += np.arctan2(numerator, denominator)

        return abs(area)

    @jit(nopython=True, parallel=True)
    def calculate_areas_numba(xyz, valid_counts):
        """Calculate areas for all cells using parallel numba acceleration."""
        grid_size = len(xyz)
        areas = np.zeros(grid_size)

        for i in prange(grid_size):
            if valid_counts[i] >= 3:
                valid_points = xyz[i, :valid_counts[i]]
                areas[i] = calculate_single_area(valid_points)

        return areas

def calculate_spherical_areas(corner_lats, corner_lons):
    """Calculate areas of spherical polygons efficiently using numba if available."""
    grid_size = len(corner_lats)

    # Convert all corners to xyz at once
    xyz = lat_lon_to_xyz(corner_lats.ravel(), corner_lons.ravel())
    xyz = xyz.reshape(grid_size, -1, 3)

    # Count valid points for each cell
    valid_mask = ~np.isnan(corner_lats)
    valid_counts = np.sum(valid_mask, axis=1)

    if HAS_NUMBA:
        # Use numba-accelerated version
        areas = calculate_areas_numba(xyz, valid_counts)
    else:
        # Fall back to non-numba version
        areas = np.zeros(grid_size)
        for i in range(grid_size):
            if valid_counts[i] >= 3:
                points = xyz[i, :valid_counts[i]]
                area = 0
                point0 = points[0]
                for j in range(len(points)-1):
                    point1 = points[j]
                    point2 = points[j+1]

                    cross_prod = np.cross(point1, point2)
                    numerator = abs(np.dot(point0, cross_prod))
                    denominator = 1.0 + np.dot(point0, point1) + \
                                np.dot(point0, point2) + np.dot(point1, point2)

                    area += np.arctan2(numerator, denominator)
                areas[i] = abs(area)

    return areas

def unstructured_to_scrip(filename, lat, lon, Opt):
    """
    Convert latitude and longitude arrays to SCRIP unstructured grid format.

    Parameters:
    filename : str
        Path to output SCRIP file
    lat : array_like
        Latitude values
    lon : array_like
        Longitude values
    Opt : dict
        Dictionary of options containing:
        - Title (str): Title for the grid file
        - Debug (bool): Enable debug output
        - GridCornerLat (array, optional): Corner latitudes if pre-calculated
        - GridCornerLon (array, optional): Corner longitudes if pre-calculated
        - GridMask (array, optional): Grid mask array
        - GridArea (array, optional): Grid area array
    """
    # Get options from Opt dictionary
    title = Opt.get("Title", f"SCRIP Grid ({len(lat)} points)")
    debug = Opt.get("Debug", False)
    grid_corners = None
    if "GridCornerLat" in Opt and "GridCornerLon" in Opt:
        grid_corners = (Opt["GridCornerLat"], Opt["GridCornerLon"])
    grid_mask = Opt.get("GridMask", None)
    grid_area = Opt.get("GridArea", None)

    # Convert inputs to numpy arrays and check dimensions
    lat = np.asarray(lat).ravel()
    lon = np.asarray(lon).ravel()

    if lat.size != lon.size:
        raise ValueError("Latitude and longitude must have the same number of elements")

    grid_size = lat.size
    grid_rank = 1  # Unstructured grid

    # Create corners if not provided
    if grid_corners is None:
        if debug:
            print("Creating spherical Voronoi tessellation...")
        corner_lats, corner_lons = find_spherical_corners(lat, lon, debug)
        grid_corners = (corner_lats, corner_lons)

    corner_lats, corner_lons = grid_corners
    num_corners = corner_lats.shape[1]

    # Calculate areas if not provided
    if grid_area is None:
        if debug:
            print("Calculating spherical areas...")
        grid_area = calculate_spherical_areas(corner_lats, corner_lons)

    # Create default mask if not provided
    if grid_mask is None:
        grid_mask = np.ones(grid_size, dtype=np.int32)

    if debug:
        print("Writing SCRIP file...")

    # Create xarray Dataset
    ds = xr.Dataset(
        coords={
            'grid_size': np.arange(grid_size),
            'grid_corners': np.arange(num_corners),
            'grid_rank': np.arange(grid_rank)
        },
        data_vars={
            'grid_dims': (['grid_rank'], [grid_size]),
            'grid_center_lat': ('grid_size', lat, {'units': 'degrees'}),
            'grid_center_lon': ('grid_size', lon, {'units': 'degrees'}),
            'grid_corner_lat': (
                ['grid_size', 'grid_corners'],
                corner_lats,
                {'units': 'degrees', '_FillValue': -9999.}
            ),
            'grid_corner_lon': (
                ['grid_size', 'grid_corners'],
                corner_lons,
                {'units': 'degrees', '_FillValue': -9999.}
            ),
            'grid_imask': ('grid_size', grid_mask, {
                'units': 'unitless',
                '_FillValue': -9999
            }),
            'grid_area': ('grid_size', grid_area, {
                'units': 'radians^2',
                'long_name': 'area weights'
            })
        },
        attrs={
            'Conventions': 'SCRIP',
            'title': title,
            'created_by': 'Python unstructured_to_scrip function'
        }
    )

    # Save to netCDF file
    ds.to_netcdf(filename)

    if debug:
        print(f"Created SCRIP file with {grid_size} points and {num_corners} corners per cell")


def spherical_polygon_area(lats, lons):
    """
    Calculate the area of a spherical polygon using L'Huilier's formula.

    Parameters:
    -----------
    lats, lons : array-like
        Arrays of latitude and longitude in radians

    Returns:
    -------
    float
        Area in steradians
    """
    if len(lats) < 3:
        return 0.0

    # Convert to 3D coordinates
    points = np.column_stack([
        np.cos(lats) * np.cos(lons),
        np.cos(lats) * np.sin(lons),
        np.sin(lats)
    ])

    # Calculate the area using spherical excess
    total = 0
    n = len(points)
    for i in range(n):
        j = (i + 1) % n
        k = (i + 2) % n

        # Get vectors for triangle
        a = points[i]
        b = points[j]
        c = points[k]

        # Calculate the spherical excess for this triangle
        numerator = abs(np.dot(a, np.cross(b, c)))
        denominator = (1 + np.dot(a, b) + np.dot(b, c) + np.dot(c, a))
        angle = 2 * np.arctan2(numerator, denominator)
        total += angle

    # Subtract out extra triangles
    area = total - (n - 2) * np.pi
    return abs(area)

def check_sphere_coverage(lat, lon, corner_lats, corner_lons, debug=False, create_plot=False):
    """
    Check if the cells properly cover the sphere.
    """
    # Convert to radians
    lat_rad = np.radians(lat)
    lon_rad = np.radians(lon)
    corner_lats_rad = np.radians(corner_lats)
    corner_lons_rad = np.radians(corner_lons)

    # Calculate areas
    areas = np.zeros(len(lat))
    for i in range(len(lat)):
        valid = ~np.isnan(corner_lats_rad[i])
        if np.sum(valid) >= 3:
            areas[i] = spherical_polygon_area(
                corner_lats_rad[i, valid],
                corner_lons_rad[i, valid]
            )

    # Rest of the function remains the same...
    total_area = np.sum(areas)
    sphere_area = 4 * np.pi
    coverage_fraction = total_area / sphere_area

    print(f"DEBUG: Total computed area: {total_area:.4f}")
    print(f"DEBUG: Sphere area: {sphere_area:.4f}")
    print(f"DEBUG: Coverage fraction: {coverage_fraction:.4f}")
    print(f"DEBUG: Number of non-zero areas: {np.sum(areas > 0)}")
    print(f"DEBUG: Area range: [{np.min(areas):.6f}, {np.max(areas):.6f}]")

    # Check for gaps between cells
    has_gaps = coverage_fraction < 0.99  # Allow 1% tolerance

    # Check for cell validity
    valid_cells = np.all(~np.isnan(corner_lats), axis=1) & np.all(~np.isnan(corner_lons), axis=1)
    invalid_cells = np.where(~valid_cells)[0]

    # Check for suspicious areas
    mean_area = np.mean(areas[areas > 0])
    std_area = np.std(areas[areas > 0])
    suspicious_areas = np.where((areas > mean_area + 3*std_area) |
                              (areas < mean_area - 3*std_area))[0]

    if debug and len(suspicious_areas) > 0:
        print("\nSuspicious Cells Found:")
        print("----------------------")
        print(f"Mean area: {mean_area:.6f}")
        print(f"Std dev:  {std_area:.6f}")
        print(f"Expected range: [{mean_area - 3*std_area:.6f}, {mean_area + 3*std_area:.6f}]")
        print("\nDetailed cell information:")
        print("  Cell ID    Latitude    Longitude        Area      Deviation")
        print("--------------------------------------------------------")
        for cell_id in suspicious_areas:
            deviation = (areas[cell_id] - mean_area) / std_area
            print(f"  {cell_id:6d}    {lat[cell_id]:9.4f}    {lon[cell_id]:9.4f}    {areas[cell_id]:9.6f}    {deviation:9.2f}σ")

    if create_plot:
        # Create visualization
        plt.figure(figsize=(15, 5))

        # Plot cell areas
        plt.subplot(131)
        plt.hist(areas, bins=50)
        plt.title('Cell Areas Distribution')
        plt.xlabel('Area')
        plt.ylabel('Count')

        # Plot points on sphere
        plt.subplot(132)
        plt.scatter(lon, lat, c=areas, cmap='viridis')
        plt.colorbar(label='Area')
        plt.title('Cell Areas on Sphere')
        plt.xlabel('Longitude')
        plt.ylabel('Latitude')

        # Plot coverage gaps if any
        plt.subplot(133)
        plt.scatter(lon, lat, c='blue', alpha=0.5, label='Valid cells')
        if len(invalid_cells) > 0:
            plt.scatter(lon[invalid_cells], lat[invalid_cells],
                       c='red', label='Invalid cells')
        if len(suspicious_areas) > 0:
            plt.scatter(lon[suspicious_areas], lat[suspicious_areas],
                       c='yellow', label='Suspicious areas')
        plt.legend()
        plt.title('Coverage Issues')
        plt.xlabel('Longitude')
        plt.ylabel('Latitude')

        plt.tight_layout()
        plt.show()

    return {
        'total_area': total_area,
        'coverage_fraction': coverage_fraction,
        'has_gaps': has_gaps,
        'num_invalid_cells': len(invalid_cells),
        'invalid_cell_indices': invalid_cells,
        'suspicious_areas': suspicious_areas,
        'mean_cell_area': mean_area,
        'std_cell_area': std_area,
        'min_area': np.min(areas),
        'max_area': np.max(areas)
    }


def test_grid_coverage(filename):
    """
    Test the coverage of a SCRIP grid file using xarray.
    """
    ds = xr.open_dataset(filename)

    lat = ds.grid_center_lat.values
    lon = ds.grid_center_lon.values
    corner_lats = ds.grid_corner_lat.values
    corner_lons = ds.grid_corner_lon.values

    results = check_sphere_coverage(lat, lon, corner_lats, corner_lons, debug=True)

    print("\nGrid Coverage Analysis:")
    print(f"Total area covered: {results['total_area']:.2f} (should be close to 4π ≈ 12.57)")
    print(f"Coverage fraction: {results['coverage_fraction']:.2%}")
    print(f"Has gaps: {results['has_gaps']}")
    print(f"Number of invalid cells: {results['num_invalid_cells']}")
    if results['num_invalid_cells'] > 0:
        print(f"Invalid cell indices: {results['invalid_cell_indices']}")
    print(f"Number of suspicious areas: {len(results['suspicious_areas'])}")
    print(f"\nCell area statistics:")
    print(f"Mean area: {results['mean_cell_area']:.6f}")
    print(f"Std dev: {results['std_cell_area']:.6f}")
    print(f"Min area: {results['min_area']:.6f}")
    print(f"Max area: {results['max_area']:.6f}")

    return results

def debug_cell(cell_id, corner_lats, corner_lons, center_lat, center_lon):
    """
    Debug a specific cell's geometry and ordering.

    Parameters:
    -----------
    cell_id : int
        The ID of the cell to debug
    corner_lats, corner_lons : ndarray
        Arrays containing corner coordinates
    center_lat, center_lon : float
        Center coordinates of the cell
    """
    import matplotlib.pyplot as plt

    # Get the corners for this cell
    clat = corner_lats[cell_id]
    clon = corner_lons[cell_id]

    # Remove any NaN values
    valid = ~np.isnan(clat)
    clat = clat[valid]
    clon = clon[valid]

    print(f"\nDEBUG INFO FOR CELL {cell_id}")
    print("----------------------------------------")
    print(f"Center coordinates: lat={center_lat:.6f}, lon={center_lon:.6f}")
    print("\nCorner coordinates:")
    for i, (lat, lon) in enumerate(zip(clat, clon)):
        print(f"Corner {i}: lat={lat:.6f}, lon={lon:.6f}")

    # Calculate angles between successive corners
    center_xyz = lat_lon_to_xyz(np.array([center_lat]), np.array([center_lon]))[0]
    corner_xyz = lat_lon_to_xyz(clat, clon)

    angles = []
    for i in range(len(corner_xyz)):
        v1 = corner_xyz[i]
        v2 = corner_xyz[(i+1) % len(corner_xyz)]
        # Compute signed angle between vectors in the tangent plane
        angle = np.arctan2(
            np.dot(np.cross(v1-center_xyz, v2-center_xyz), center_xyz),
            np.dot(v1-center_xyz, v2-center_xyz)
        )
        angles.append(np.degrees(angle))

    print("\nAngles between successive corners (degrees):")
    angle_sum = 0
    for i, angle in enumerate(angles):
        angle_sum += angle
        print(f"Angle {i}->{i+1}: {angle:.2f}")
    print(f"Sum of angles: {angle_sum:.2f} (should be close to 360 for CCW)")

#     # Plot the cell
#     plt.figure(figsize=(10, 10))
#
#     # Plot corners
#     plt.scatter(clon, clat, c='blue', label='Corners')
#
#     # Plot center
#     plt.scatter(center_lon, center_lat, c='red', label='Center')
#
#     # Connect corners with lines
#     for i in range(len(clat)):
#         next_i = (i + 1) % len(clat)
#         plt.plot([clon[i], clon[next_i]], [clat[i], clat[next_i]], 'b-')
#         # Add corner numbers
#         plt.annotate(f'C{i}', (clon[i], clat[i]))
#
#     plt.title(f'Cell {cell_id} Geometry')
#     plt.xlabel('Longitude')
#     plt.ylabel('Latitude')
#     plt.legend()
#     plt.grid(True)
#     plt.show()

    # Check for potential issues
    print("\nPotential issues:")

    # Check for very close corners
    for i in range(len(corner_xyz)):
        for j in range(i+1, len(corner_xyz)):
            dist = np.linalg.norm(corner_xyz[i] - corner_xyz[j])
            if dist < 1e-6:
                print(f"WARNING: Corners {i} and {j} are very close together")

    # Check for colinear points
    for i in range(len(corner_xyz)):
        v1 = corner_xyz[i] - center_xyz
        v2 = corner_xyz[(i+1) % len(corner_xyz)] - center_xyz
        cross_prod = np.linalg.norm(np.cross(v1, v2))
        if cross_prod < 1e-6:
            print(f"WARNING: Corners {i} and {i+1} are nearly colinear with center")

    # Check for clockwise ordering
    if angle_sum < 0:
        print("WARNING: Corners appear to be ordered clockwise")

    return angles

def test_specific_cell(filename, cell_id):
    """
    Test a specific cell in a SCRIP grid file.
    """
    ds = xr.open_dataset(filename)

    corner_lats = ds.grid_corner_lat.values
    corner_lons = ds.grid_corner_lon.values
    center_lat = ds.grid_center_lat.values[cell_id]
    center_lon = ds.grid_center_lon.values[cell_id]

    debug_cell(cell_id, corner_lats, corner_lons, center_lat, center_lon)


def validate_scrip_file(filename):
    """
    Validate a SCRIP grid file and print results.

    Parameters:
    -----------
    filename : str
        Path to SCRIP file to validate
    """
    print("\nValidating SCRIP file:", filename)
    print("============================")

    # Get coverage results
    results = test_grid_coverage(filename)

    # Check critical issues
    has_issues = False

    if results['coverage_fraction'] < 0.99:
        print(f"ERROR: Incomplete coverage: {results['coverage_fraction']:.1%}")
        has_issues = True

    if results['num_invalid_cells'] > 0:
        print(f"ERROR: Found {results['num_invalid_cells']} invalid cells")
        print(f"       First few invalid cells: {results['invalid_cell_indices'][:5]}")
        has_issues = True

    if len(results['suspicious_areas']) > 0:
        print(f"WARNING: Found {len(results['suspicious_areas'])} cells with suspicious areas")

    # Print basic statistics
    print(f"\nGrid size: {results['total_area']:.2f}")
    print(f"Mean cell area: {results['mean_cell_area']:.6f}")
    print(f"Area range: [{results['min_area']:.6f}, {results['max_area']:.6f}]")

    if not has_issues:
        print("\nNo critical issues found.")

    return not has_issues  # Return True if valid, False if issues found