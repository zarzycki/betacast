from math import sin, cos, sqrt, pi, atan, exp, radians, degrees
import re
import os
import gzip
import logging
import shutil
import sys

# Third party
import numpy as np
import matplotlib.pyplot as plt

# Betacast
import meteo
import packing
from constants import (
    grav
)

logger = logging.getLogger(__name__)

def calc_mass_weighted_integral(var3d, pdel, area, gravit=grav):
    """
    Calculate mass-weighted global integral of a 3D field
    Parameters:
    -----------
    var3d : numpy.ndarray
        3D array of variable (time, lev, ncol) or (lev, ncol)
    pdel : numpy.ndarray
        Pressure thickness (delta-p) array matching var3d dimensions
    area : numpy.ndarray
        Grid cell areas (ncol)
    gravit : float, optional
        Gravitational constant, default 9.81 m/s^2
    Returns:
    --------
    float
        Mass-weighted global integral
    """
    # Calculate mass-weighted column integral
    column_integral = np.sum(var3d * pdel, axis=0) / gravit
    # Calculate global sum (not mean)
    global_integral = np.sum(column_integral * area)
    return global_integral


def calc_inverse_mass_weighted_field(target_integral, pdel, area, gravit=grav):
    """
    Calculate a uniform 3D field that would produce the target mass-weighted global integral
    Parameters:
    -----------
    target_integral : float
        Target value for the mass-weighted global integral
    pdel : numpy.ndarray
        Pressure thickness (delta-p) array (time, lev, ncol) or (lev, ncol)
    area : numpy.ndarray
        Grid cell areas (ncol)
    gravit : float, optional
        Gravitational constant, default 9.81 m/s^2
    Returns:
    --------
    numpy.ndarray
        3D field that when integrated gives target_integral
    """
    # Get total mass weighting (without area normalization)
    mass_weight = np.sum((pdel/gravit) * area[None,:])
    # Calculate uniform value needed
    uniform_value = target_integral / mass_weight
    # Create ND field with this uniform value
    if len(pdel.shape) == 3:
        result = np.full_like(pdel[0], uniform_value)
    else:
        result = np.full_like(pdel, uniform_value)
    return result


def convert_lon(tmplon, jcode):
    """
    Convert longitude based on jcode.
    jcode >= 0 --> lon: 0 to 360
    jcode <  0 --> lon: -180 to 180
    """
    tmplon = np.asarray(tmplon)
    if jcode >= 0:
        tmplon[tmplon < 0] += 360
    else:
        tmplon[tmplon < -180] += 360
        tmplon[tmplon > 180] -= 360
    return tmplon


def replace_or_add_variable(file_path, variable, value):
    """
    Replace or add a variable assignment in a text file.

    Removes any existing assignments to the variable and appends a new
    assignment at the end of the file.

    Args:
        file_path (str): Path to the file to modify.
        variable (str): Variable name to replace or add.
        value (str): Value to assign to the variable.

    Example:
        >>> replace_or_add_variable('config.txt', 'DEBUG', 'True')
        # If config.txt contained "DEBUG = False", it will be replaced with
        # "DEBUG = True" at the end of the file.
    """
    with open(file_path, 'r') as file:
        lines = file.readlines()

    # Eliminate any leading or trailing spaces in the delete command
    re_pattern = re.compile(rf"^\s*{variable}\s*=\s*")
    lines = [line for line in lines if not re_pattern.match(line)]

    # Add the new variable at the end of the file
    lines.append(f"{variable} = {value}\n")

    # Write the modified content back to the file
    with open(file_path, 'w') as file:
        file.writelines(lines)


def gc_latlon(lat1, lon1, lat2, lon2, npts, iu):
    """
    Calculate the great-circle distance between two points on the surface of the Earth, and interpolate points.

    Parameters:
    lat1, lon1 : float or array-like
        Latitude and longitude of the first point(s). Can be scalars or arrays of the same size as lat2/lon2.
    lat2, lon2 : float or array-like
        Latitude and longitude of the second point(s). Must have the same size.
    npts : int
        Number of equally-spaced points to interpolate along the great circle.
    iu : int
        Unit and longitude conversion code. Positive means 0-360, negative means -180 to 180.
        abs(iu) = 1: radians, 2: degrees, 3: meters, 4: kilometers.

    Returns:
    dist : float or array-like
        Great-circle distance between the points.
    Attributes:
        gclat : array-like
            Latitudes of the interpolated points, including endpoints.
        gclon : array-like
            Longitudes of the interpolated points, including endpoints.
        spacing : float
            The distance between interpolated points.
        units : str
            Units of the distance ('radians', 'degrees', 'meters', 'kilometers').
    """
    # Convert lat/lon to radians for calculations
    lat1_rad = np.radians(lat1)
    lon1_rad = np.radians(lon1)
    lat2_rad = np.radians(lat2)
    lon2_rad = np.radians(lon2)

    # Calculate great circle distance (Haversine formula)
    dlat = lat2_rad - lat1_rad
    dlon = lon2_rad - lon1_rad
    a = np.sin(dlat / 2)**2 + np.cos(lat1_rad) * np.cos(lat2_rad) * np.sin(dlon / 2)**2
    c = 2 * np.arctan2(np.sqrt(a), np.sqrt(1 - a))

    # Select appropriate radius based on 'iu'
    if abs(iu) == 1:
        radius = 1.0  # Radians
        unit_str = "radians"
    elif abs(iu) == 2:
        radius = 180 / np.pi  # Degrees
        unit_str = "degrees"
    elif abs(iu) == 3:
        radius = 6371220.0  # Meters
        unit_str = "meters"
    elif abs(iu) == 4:
        radius = 6371.220  # Kilometers
        unit_str = "kilometers"
    else:
        raise ValueError("Invalid value for 'iu'. Must be 1, 2, 3, or 4.")

    dist = radius * c  # Great circle distance in the specified units

    # If npts <= 2, return just the distance and endpoints
    if npts <= 2:
        return dist, {
            "gclat": np.concatenate([np.atleast_1d(lat1), np.atleast_1d(lat2)]),
            "gclon": np.concatenate([np.atleast_1d(lon1), np.atleast_1d(lon2)]),
            "spacing": dist,
            "units": unit_str
        }

    # Interpolate points along the great circle path
    f = np.linspace(0, 1, npts)  # Fraction along the path
    sin_lat1 = np.sin(lat1_rad)
    sin_lat2 = np.sin(lat2_rad)
    cos_lat1 = np.cos(lat1_rad)
    cos_lat2 = np.cos(lat2_rad)

    # Calculate interpolated latitudes and longitudes
    A = np.sin((1 - f) * c) / np.sin(c)
    B = np.sin(f * c) / np.sin(c)
    x = A * cos_lat1 * np.cos(lon1_rad) + B * cos_lat2 * np.cos(lon2_rad)
    y = A * cos_lat1 * np.sin(lon1_rad) + B * cos_lat2 * np.sin(lon2_rad)
    z = A * sin_lat1 + B * sin_lat2

    lat_interp = np.arctan2(z, np.sqrt(x**2 + y**2))
    lon_interp = np.arctan2(y, x)

    # Convert lat/lon back to degrees for output
    lat_interp_deg = np.degrees(lat_interp)
    lon_interp_deg = np.degrees(lon_interp)

    # Convert longitudes to the correct range based on 'iu'
    lon_interp_deg = convert_lon(lon_interp_deg, iu)

    # Calculate the average spacing between the interpolated points
    spacing = dist / (npts - 1)

    # Return distance and interpolated points
    return dist, {
        "gclat": lat_interp_deg,
        "gclon": lon_interp_deg,
        "spacing": spacing,
        "units": unit_str
    }


def convert_lon(lon, iu):
    """
    Convert longitudes to either 0-360 or -180 to 180 based on the iu flag.

    Parameters:
    lon : array-like
        Longitude values to be converted.
    iu : int
        Positive: convert to 0-360 range, Negative: convert to -180 to 180 range.

    Returns:
    Converted longitude values.
    """
    lon = np.asarray(lon)
    if iu > 0:
        lon = np.where(lon < 0, lon + 360, lon)
    else:
        lon = np.where(lon > 180, lon - 360, lon)
    return lon


def tctestcase(cen_lon, cen_lat, dp, rp, zp, exppr, gamma_, lon, lat, p, z, zcoords, psin, uin, vin, Tin, qin, invert_vortex, modify_q, modify_q_mult, return_anoms=False, debug=False):
    if debug:
        print(f"cen_lon: {cen_lon}, cen_lat: {cen_lat}, dp: {dp}, rp: {rp}, zp: {zp}, exppr: {exppr}, gamma_: {gamma_}, lon: {lon}, lat: {lat}, p: {p}, z: {z}, zcoords: {zcoords}, psin: {psin}, uin: {uin}, vin: {vin}, Tin: {Tin}, qin: {qin}, invert_vortex: {invert_vortex}, modify_q: {modify_q}, modify_q_mult: {modify_q_mult}")

    if rp <= 0.0:
        rp = 273000.  # Default value
    if dp <= 0.0:
        dp = 2280.  # Default value
    if zp <= 0.0:
        zp = 12000.  # Default value
    if gamma_ <= 0.0:
        gamma_ = 0.006  # Default value
    if exppr <= 0.0:
        exppr = 1.5  # Default value
    if cen_lat <= -900.:
        cen_lat = 20.  # Default value
    if cen_lon <= -900.:
        cen_lon = -40.  # Default value
    if modify_q_mult < 0.:
        modify_q_mult = 1.0  # Default value

    # Constants
    a = 6371220.  # Earth's Radius (m)
    Rd = 287.0  # Ideal gas const dry air (J kg^-1 K^1)
    g = 9.80616  # Gravity (m s^2)
    omega = 7.292115e-5  # angular velocity 1/s
    convert = 180. / pi  # conversion factor: radians to degrees
    q0 = 0.021  # q at surface from Jordan
    Ts0 = 302.0  # Surface temperature (SST)
    p00 = psin #101500.  # global mean surface pressure
    p0 = 100000.  # p for model level calculation
    zq1 = 3000.  # Height 1 for q calculation
    zq2 = 8000.  # Height 2 for q calculation
    exppz = 2.0  # Exponent for z dependence of p
    ztrop = 20000.  # Tropopause Height
    qtrop = 1.e-11  # Tropopause specific humidity
    constTv = 0.608  # Constant for Virtual Temp Conversion
    deltaz = 2.e-5  # Small number to ensure convergence in FPI
    epsilon = 1.e-25  # Small number to avoid dividing by zero in wind calc
    exponent = Rd * gamma_ / g  # exponent
    T0 = Ts0 * (1. + constTv * q0)  # Surface temp
    Ttrop = T0 - gamma_ * ztrop  # Tropopause temp
    ptrop = p00 * (Ttrop / T0) ** (1. / exponent)  # Tropopause pressure

    if debug:
        print(f"cen_lon: {cen_lon}, cen_lat: {cen_lat}, dp: {dp}, rp: {rp}, zp: {zp}, exppr: {exppr}, gamma_: {gamma_}, lon: {lon}, lat: {lat}, p: {p}, z: {z}, zcoords: {zcoords}, psin: {psin}, uin: {uin}, vin: {vin}, Tin: {Tin}, qin: {qin}, invert_vortex: {invert_vortex}, modify_q: {modify_q}, modify_q_mult: {modify_q_mult}")

    # Coriolis parameter
    f = 2. * omega * sin(radians(cen_lat))

    # Great circle distance calculation (assuming `gc_latlon` is provided)
    gr, _ = gc_latlon(cen_lat, cen_lon, lat, lon, 2, 3)

    if debug:
        print(f"gr: {gr}")
        print(f"dp: {dp}")
        print(f"rp: {rp}")
        print(f"exppr: {exppr}")
        print(f"p00: {p00}")
        print(f"T0: {T0}")
        print(f"exponent: {exponent}")

    if zcoords == 1:
        height = z
        if height > ztrop:
            p = ptrop * safe_exp(-g * (height - ztrop) / (Rd * Ttrop))
        else:
            p = (p00 - dp * safe_exp(-(gr / rp) ** exppr) * safe_exp(-(height / zp) ** exppz)) * \
                ((T0 - gamma_ * height) / T0) ** (1 / exponent)
    else:
        ps = -dp * safe_exp(-(gr / rp) ** exppr) + p00
        height = (T0 / gamma_) * (1. - (p / ps) ** exponent)

    # Initialize U and V (wind components)
    d1 = sin(radians(cen_lat)) * cos(radians(lat)) - cos(radians(cen_lat)) * sin(radians(lat)) * cos(radians(lon) - radians(cen_lon))
    d2 = cos(radians(cen_lat)) * sin(radians(lon) - radians(cen_lon))
    d = max(epsilon, sqrt(d1 ** 2. + d2 ** 2.))

    ufac = d1 / d
    vfac = d2 / d

    # Flip ufac and vfac in SH
    if cen_lat < 0:
        ufac = -ufac
        vfac = -vfac

    if height > ztrop:
        u = uin
        v = vin
    else:
        if debug:
            print(f"height: {height}, zp: {zp}, exppz: {exppz}, dp: {dp}, rp: {rp}, gr: {gr}")
        vt = (-abs(f) * gr / 2 + sqrt((f * gr / 2) ** 2 - (exppr * (gr / rp) ** exppr) *
                                 (Rd * (T0 - gamma_ * height)) / (exppz * height * Rd * (T0 - gamma_ * height) / (g * zp ** exppz) + 1. - p00 / dp * safe_exp((gr / rp) ** exppr) * safe_exp((height / zp) ** exppz))))
        v = vin + vfac * vt
        u = uin + ufac * vt

    # Calculate temperature and specific humidity
    Tvin = Tin * (1 + constTv * qin)
    t = (Tvin) / (1 + constTv * qin) / (1 + exppz * Rd * Tvin * height / (g * zp ** exppz * (1 - p00 / dp * safe_exp((gr / rp) ** exppr) * safe_exp((height / zp) ** exppz))))

    if modify_q:
        rh = 0.263 * p * qin / safe_exp((17.67 * (Tin - 273.15)) / (Tin - 29.65))
        if height > ztrop:
            q = qin
        else:
            q = rh / (0.263 * p) * safe_exp((17.67 * (t - 273.15)) / (t - 29.65))
    else:
        q = qin

    # Assemble outputs
    output = np.zeros(5, dtype=float)
    delv = v - vin
    delu = u - uin
    delt = t - Tin
    delq = q - qin
    delps = ps - psin

    if modify_q:
        delq *= modify_q_mult
        q = qin + delq

    if return_anoms:
        output[0] = delv
        output[1] = delu
        output[2] = delq
        output[3] = delt
        output[4] = delps
        if invert_vortex:
            output = -output
    else:
        if invert_vortex:
            output[0] = vin - delv
            output[1] = uin - delu
            output[2] = qin - delq
            output[3] = Tin - delt
            output[4] = psin - delps
        else:
            output[0] = v
            output[1] = u
            output[2] = q
            output[3] = t
            output[4] = ps

    return output


def get_rp_from_dp_rmw(cen_lat, dp, target_rmw, debug=False):
    # Constants
    a = 6371220.  # Earth's Radius (m)
    Rd = 287.0  # Ideal gas const dry air (J kg^-1 K^1)
    g = 9.80616  # Gravity (m s^2)
    omega = 7.292115e-5  # angular velocity 1/s
    convert = 180. / pi  # conversion factor: radians to degrees
    Ts0 = 302.0  # Surface temperature (SST)
    p00 = 101500.  # global mean surface pressure
    q0 = 0.021  # q at surface from Jordan
    gamma_ = 0.007  # lapse rate
    exponent = Rd * gamma_ / g  # exponent

    # Initial parameters
    rpi = 200000.  # First guess
    step = 5000.  # Step size in meters
    maxiter = 300  # Maximum iterations
    err_stop = 111.  # Stop when error is smaller than this

    T0 = Ts0 * (1. + 0.608 * q0)
    f = 2. * omega * sin(radians(cen_lat))  # Coriolis parameter

    for jj in range(maxiter):
        rp = rpi
        rr = np.linspace(0, 1000000., 10001)
        vt = np.zeros_like(rr)

        for ii, r in enumerate(rr):
            T1 = -(abs(f) * r) / 2.
            T2 = (f ** 2 * r ** 2) / 4.
            NUM = (3. / 2.) * ((r / rp) ** (3. / 2.)) * T0 * Rd
            DEN = 1. - (p00 / dp) * safe_exp((r / rp) ** (3. / 2.))
            diff = T2 - (NUM / DEN)
            if diff < 0:
                if debug:
                    print(f"Skipping r={r:.1f}: sqrt arg < 0 (T2={T2:.2e}, NUM/DEN={NUM / DEN:.2e}, diff={diff:.2e})")
                vt[ii] = 0.0  # Or np.nan if you prefer to catch it later
            else:
                vt[ii] = T1 + sqrt(diff)


        vmax = np.max(vt)
        rmw = rr[np.argmax(vt)]

        err_here = rmw - target_rmw
        if debug:
            print(f"iter: {jj}  rpi: {rpi}  rmw: {rmw}  target: {target_rmw}  err: {abs(err_here)}")

        if abs(err_here) < err_stop:
            break

        if jj > 0 and erri * err_here < 0:
            step /= 2.

        if err_here >= 0:
            rpi -= step
        else:
            rpi += step

        erri = err_here

        if jj == maxiter - 1:
            print("Did not converge!")
            print(f"Setting rpi to target_rmw ({target_rmw})!")
            rpi = target_rmw

    print(f"dp={dp}, rp={rpi}")
    return rpi


def keyword_values(namelist_file, key, return_type, default=None, verbose=False):
    """
    Reads a namelist file and extracts the value for the specified key.

    Parameters:
    -----------
    namelist_file (str): Path to the namelist file.
    key (str): The key whose value is to be extracted.
    return_type (str): The type to which the value should be converted. Options: "int", "float", "bool", "str".
    default: Optional default value to return if the key is not found (default is None).

    Returns:
    --------
    The value associated with the key in the specified return type, or the default value if the key is not found.
    """
    types = ["int", "float", "bool", "str"]
    # Validate return_type
    if return_type not in types:
        raise ValueError(f"Unsupported return_type '{return_type}'. Must be one of {types}.")

    def strip_unescaped_quotes(value):
        """
        Strips unescaped single and double quotes from a string.
        Preserves escaped quotes (\' and \").
        """
        if not isinstance(value, str):
            return value

        # Process the string character by character
        result = []
        i = 0
        while i < len(value):
            # Check for escaped quotes
            if i < len(value) - 1 and value[i] == '\\' and value[i + 1] in ['"', "'"]:
                result.append(value[i + 1])  # Keep the quote but drop the escape character
                i += 2
            # Skip unescaped quotes
            elif value[i] in ['"', "'"]:
                i += 1
            else:
                result.append(value[i])
                i += 1

        return ''.join(result).strip()

    try:
        # Read the namelist file
        with open(namelist_file, 'r') as f:
            lines = f.readlines()

        # Remove whitespace and split key-value pairs
        for line in lines:
            # Ignore empty lines and comments
            line = line.strip()
            if line and not line.startswith('#'):
                # Split key and value
                if '=' in line:
                    k, v = [x.strip() for x in line.split('=', 1)]  # Strip whitespace around key and value
                    # If we found the key
                    if k == key:
                        # Strip quotes before type conversion
                        v = strip_unescaped_quotes(v)

                        # Handle conversion to requested type
                        if return_type == "int":
                            value = int(v)
                        elif return_type == "float":
                            value = float(v)
                        elif return_type == "bool":
                            value = v.lower() in ['true', 't', '1']
                        elif return_type == "str":
                            value = os.path.expandvars(v)
                        if verbose:
                            logging.info(f"[keyword_values] Key '{key}' found with value '{value}' (type: {return_type})")
                        return value

        # If key is not found, return the default
        if verbose:
            logging.info(f"[keyword_values] Key '{key}' not found. Returning default: {default}")
        return default

    except Exception as e:
        # If any exception occurs (e.g., file not found), return the default
        return default


def radialAvg2D_unstruc(data, lat, lon, deltaMax, psminlat, psminlon, outerRad, mergeInnerBins):
    # km per degree at the equator
    kmInDeg = 111.32
    kmGrid = kmInDeg * deltaMax
    logging.info(f"The max lat/lon km grid spacing at equator is {kmGrid} km")

    ncol = len(lat)
    pi = np.pi
    d2r = pi / 180.0

    # Prepare radial bins
    timesGrid = 1.1  # We want each radius bin to be timesGrid * kmGrid
    nx = int(outerRad / (timesGrid * kmGrid))
    logging.info(f"Number of bins is equal to {nx}")

    if mergeInnerBins:
        numMerge = 2  # Number of innermost radial bins to merge
        logging.info(f"Merging innermost {numMerge} bins because radial average is so small.")
        numMergeMinusOne = numMerge - 1
        origRadiusArr = np.linspace(0, outerRad, nx + numMergeMinusOne)
        radiusArr = np.zeros(len(origRadiusArr) - numMergeMinusOne)
        radiusArr[0] = origRadiusArr[0]
        radiusArr[1:] = origRadiusArr[numMerge:]
    else:
        logging.info("Not merging any innermost bins -- be careful that your inner bins have > 1 pt.")
        radiusArr = np.linspace(0, outerRad, nx)

    numRadBins = len(radiusArr)

    rad_thevar_hit = np.zeros(numRadBins, dtype=int)
    rad_thevar_cum = np.zeros(numRadBins, dtype=float)

    # Iterate over each point in the dataset
    logging.info("Starting loop")
    for i in range(ncol):
        # Use the gc_latlon function to calculate the great circle distance
        gcdist, _ = gc_latlon(psminlat, psminlon, lat[i], lon[i], 2, 4)

        if gcdist <= outerRad:
            bin_idx = np.argmin(np.abs(radiusArr - gcdist))
            rad_thevar_hit[bin_idx] += 1
            rad_thevar_cum[bin_idx] += data[i]

    logging.info(f"Minimum number of hits per gridbox: {np.min(rad_thevar_hit)}")
    if np.min(rad_thevar_hit) == 0:
        logging.info("WARNING: Some radial bins have zero hits - radial average will be unreliable")
    logging.info(f"Maximum number of hits per gridbox: {np.max(rad_thevar_hit)}")

    # Calculate the radial average
    rad_thevar = np.divide(rad_thevar_cum, rad_thevar_hit, out=np.zeros_like(rad_thevar_cum), where=rad_thevar_hit != 0)

    # Returning as dictionary for simplicity, with 'radius' key and the radial averages
    return {
        'radius': radiusArr,
        'radial_average': rad_thevar,
        'hit_count': rad_thevar_hit
    }


def radialAvg3D_unstruc(data, lat, lon, lev, deltaMax, psminlat, psminlon, outerRad, mergeInnerBins):
    # km per degree at the equator
    kmInDeg = 111.32
    kmGrid = kmInDeg * deltaMax
    logging.info(f"The max lat/lon km grid spacing at equator is {kmGrid} km")

    ncol = len(lat)
    nlev = len(lev)
    pi = np.pi
    d2r = pi / 180.0

    # Prepare radial bins
    timesGrid = 1.1  # We want each radius bin to be timesGrid * kmGrid
    nx = int(outerRad / (timesGrid * kmGrid))
    logging.info(f"Number of bins is equal to {nx}")

    if mergeInnerBins:
        numMerge = 2  # Number of innermost radial bins to merge
        logging.info(f"Merging innermost {numMerge} bins because radial average is so small.")
        numMergeMinusOne = numMerge - 1
        origRadiusArr = np.linspace(0, outerRad, nx + numMergeMinusOne)
        radiusArr = np.zeros(len(origRadiusArr) - numMergeMinusOne)
        radiusArr[0] = origRadiusArr[0]
        radiusArr[1:] = origRadiusArr[numMerge:]
    else:
        logging.info("Not merging any innermost bins -- be careful that your inner bins have > 1 pt.")
        radiusArr = np.linspace(0, outerRad, nx)

    numRadBins = len(radiusArr)

    rad_thevar_hit = np.zeros((nlev, numRadBins), dtype=int)
    rad_thevar_cum = np.zeros((nlev, numRadBins), dtype=float)

    total_hits = 0
    total_loops = 0

    # Iterate over each column in the dataset
    logging.info("Starting loop")
    for i in range(ncol):
        total_loops += 1
        # Use the gc_latlon function to calculate the great circle distance
        gcdist, _ = gc_latlon(psminlat, psminlon, lat[i], lon[i], 2, 4)

        if gcdist <= outerRad:
            bin_idx = np.argmin(np.abs(radiusArr - gcdist))
            rad_thevar_hit[:, bin_idx] += 1
            rad_thevar_cum[:, bin_idx] += data[:, i]
            total_hits += 1

    # Handle grid boxes with no hits
    rad_thevar_hit = np.where(rad_thevar_hit == 0, np.nan, rad_thevar_hit)

    logging.info(f"Minimum number of hits per gridbox: {np.nanmin(rad_thevar_hit)}")
    logging.info(f"Maximum number of hits per gridbox: {np.nanmax(rad_thevar_hit)}")
    logging.info(f"Total hits found: {total_hits}")
    logging.info(f"Total loops conducted: {total_loops}")

    # Calculate the radial average
    rad_thevar = np.divide(rad_thevar_cum, rad_thevar_hit, out=np.zeros_like(rad_thevar_cum), where=rad_thevar_hit != 0)

    # Returning the result as a dictionary for simplicity
    return {
        'radius': radiusArr,
        'radial_average': rad_thevar,
        'hit_count': rad_thevar_hit
    }


def q_for_optim(p, alpha=0.02):
    """
    Calculate specific humidity (kg/kg) as a function of pressure (mb)
    with a shape parameter to control the rate of decrease

    Parameters:
    -----------
    p : float or array-like
        Pressure in millibars (mb)
    alpha : float, optional
        Shape parameter controlling how rapidly q falls off with height
        alpha > 1: faster decrease
        alpha < 1: slower decrease

    Returns:
    --------
    q : float or array-like
        Specific humidity in kg/kg
    """
    # Calculate constants based on constraints and shape parameter
    # q(1000 mb) = 0.21 kg/kg and q(50 mb) = 1e-8 kg/kg

    # For numerical stability when using different alpha values
    p1 = 1000.0  # surface pressure
    p2 = 50.0    # upper atmosphere pressure
    q1 = 0.21    # surface specific humidity
    q2 = 1e-08   # upper atmosphere specific humidity

    # Calculate parameter B
    B = np.log(q1/q2) / (p1**alpha - p2**alpha)

    # Calculate parameter A
    A = q1 / safe_exp(B * p1**alpha)

    # Calculate specific humidity
    return A * safe_exp(B * p**alpha)


def read_tcseed_settings(vortex_input):
    """
    Read TC settings from namelist file or dictionary and set defaults.

    Parameters:
    -----------
    vortex_input : str or dict
        Either path to the namelist file or dictionary containing TC settings

    Returns:
    --------
    tc_settings : dict
        Dictionary containing all TC settings and optimization parameters
    """

    # Define all default values
    defaults = {
        # TC location and basic parameters
        'invert_vortex': False,
        'deltaMax': 1.0,
        'psminlat': 20.0,
        'psminlon': -40.0,
        'modify_q': False,
        'modify_q_mult': 1.0,
        'gamma_': 0.0065,
        # Optional: Read existing optimized parameters if available
        'minp': None,
        'target_rmw': None,
        'rp': None,
        'dp': None,
        'zp': None,
        'exppr': None,
        'degree_adjust': 15.0,
        # Optimization settings (can be overridden in namelist)
        'truePS_scan_radius': 0.5,
        'rad_for_corr': 800.0,
        'n_p_steps': 32,
        'n_t_steps': 60,
        # Optimization ranges (can be customized in namelist)
        'rp_min': 80000.0,
        'rp_max': 300000.0,
        'dp_min': 200.0,
        'dp_max': 6000.0,
        'exppr_min': 1.1,
        'exppr_max': 1.9,
        'zp_min': 6000.0,
        'zp_max': 16000.0,
        # Round 2 refinement factor
        'mult_factor': 0.52,
    }

    # Type mapping for file reading
    type_mapping = {
        'invert_vortex': 'bool',
        'modify_q': 'bool',
        'deltaMax': 'float',
        'psminlat': 'float',
        'psminlon': 'float',
        'modify_q_mult': 'float',
        'gamma_': 'float',
        'minp': 'float',
        'target_rmw': 'float',
        'rp': 'float',
        'dp': 'float',
        'zp': 'float',
        'exppr': 'float',
        'degree_adjust': 'float',
        'truePS_scan_radius': 'float',
        'rad_for_corr': 'float',
        'n_p_steps': 'int',
        'n_t_steps': 'int',
        'rp_min': 'float',
        'rp_max': 'float',
        'dp_min': 'float',
        'dp_max': 'float',
        'exppr_min': 'float',
        'exppr_max': 'float',
        'zp_min': 'float',
        'zp_max': 'float',
        'mult_factor': 'float',
    }

    if isinstance(vortex_input, str):
        # Read from namelist file

        tc_settings = {}
        for key, default_value in defaults.items():
            tc_settings[key] = keyword_values(
                vortex_input,
                key if key != 'gamma_' else 'gamma',  # Handle gamma_ special case
                type_mapping[key],
                default=default_value
            )

        # Store namelist path for posterity
        tc_settings['namelist_path'] = vortex_input
        logging.info(f"TC settings loaded from {vortex_input}")

    elif isinstance(vortex_input, dict):
        # Use dictionary input with defaults
        tc_settings = defaults.copy()

        # Update with provided values
        for key, value in vortex_input.items():
            if key in tc_settings:
                tc_settings[key] = value
            else:
                logging.warning(f"Unknown setting '{key}' provided, ignoring")

        # No namelist path since we're using dict input
        tc_settings['namelist_path'] = None
        logging.info("TC settings loaded from dictionary input")

    else:
        raise TypeError("vortex_input must be either a string (file path) or dictionary")

    # Check for required parameters
    required_params = ['psminlat', 'psminlon']
    missing_params = []

    for param in required_params:
        if tc_settings[param] is None:
            missing_params.append(param)

    if missing_params:
        error_msg = f"Missing required parameters: {', '.join(missing_params)}"
        logging.error(error_msg)
        print(f"ERROR: {error_msg}")
        print("Required parameters: target_rmw, minp, psminlat, psminlon")
        sys.exit(1)

    # Validate gamma_ value
    if tc_settings['gamma_'] < 0.0:
        tc_settings['gamma_'] = 0.0065
        logging.info(f"Setting gamma_ to default of {tc_settings['gamma_']}")

    return tc_settings


def find_actual_center(guess_lat, guess_lon, lat, lon, pressure, search_radius=5.0, psunits='Pa'):
    """
    Find the actual TC center by locating minimum pressure within a search radius.

    Parameters
    ----------
    guess_lat, guess_lon : float
        Initial guess for TC center coordinates (degrees)
    lat, lon : array-like (1-D unstructured)
        Grid point latitude and longitude arrays
    pressure : array-like (1-D unstructured)
        Surface pressure values at each grid point
    search_radius : float, optional
        Search radius in great circle degrees, default=15.0
    search_radius : string, optional
        Units of the pressure array, default='Pa'

    Returns
    -------
    actual_lat, actual_lon : float
        Coordinates of the located TC center
    """

    # Calculate great circle distances (GCD) from guessed center to all grid points
    gcdist, _ = gc_latlon(guess_lat, guess_lon, lat, lon, 2, 2)

    logging.info(f"find_actual_center: search_radius: {search_radius}")

    # Create a masked array where points outside search radius are set to NaN
    masked_pressure = np.where(gcdist > search_radius, np.nan, pressure)

    # Report number of cells remaining after masking
    cells_remaining = np.sum(~np.isnan(masked_pressure))
    total_cells = len(pressure)
    percentage_remaining = (cells_remaining / total_cells) * 100
    logging.info(f"find_actual_center: Cells remaining after masking: {cells_remaining} of {total_cells} ({percentage_remaining:.6f}%)")

    # Find the index of the minimum pressure point
    min_index = np.nanargmin(masked_pressure)

    # Get coordinates of the actual center
    actual_lat = lat[min_index]
    actual_lon = lon[min_index]

    # Log the results
    logging.info(f"Initial TC center guess: lat={guess_lat:.2f}, lon={guess_lon:.2f}")
    logging.info(f"Actual TC center found ({pressure[min_index]:.1f} {psunits}): lat={actual_lat:.2f}, lon={actual_lon:.2f}")

    return actual_lat, actual_lon, pressure[min_index]


def find_fill_parameters(mps, T, lat, lon, lev, tc_settings, plot_bestfits=True, debug=False):
    """
    Find and fill missing TC parameters (rp, dp, zp, exppr, cen_lat, cen_lon)
    by optimizing against observed atmospheric data.

    Parameters:
    -----------
    data_vars : dict
        Dictionary containing atmospheric fields with keys: 'ps', 't', 'lat', 'lon', 'lev'
    tc_settings : dict
        Dictionary containing TC settings from read_tcseed_settings()

    Returns:
    --------
    tc_settings : dict
        Dictionary with all 6 optimized TC parameters plus metadata
    """

    # Extract settings
    psminlat = tc_settings['psminlat']
    psminlon = tc_settings['psminlon']
    gamma_ = tc_settings['gamma_']
    modify_q = tc_settings['modify_q']
    modify_q_mult = tc_settings['modify_q_mult']
    deltaMax = tc_settings['deltaMax']
    truePS_scan_radius = tc_settings['truePS_scan_radius']
    rad_for_corr = tc_settings['rad_for_corr']
    n_p_steps = tc_settings['n_p_steps']
    n_t_steps = tc_settings['n_t_steps']

    # Convert pressure levels if needed (assumes mb if max < 2000)
    if np.max(lev) < 2000:
        lev_mb = lev
    else:
        # We have Pa but the TC code needs hPa
        lev_mb = lev / 100.0

    if mps.ndim == 2:
        # PS is a 2D grid, assume lat/lon, need to ncol-ize
        logging.info("Source data is structured")
        is_unstructured = False
        lat_flat, lon_flat = packing.flatten_1D_lat_lon(lat, lon)
        mps_flat = packing.latlon_to_ncol(mps)
    else:
        # Unstructured grid: lat(ncol), lon(ncol), ps(ncol)
        logging.info("Source data is unstructured")
        is_unstructured = True
        lat_flat = lat
        lon_flat = lon
        mps_flat = mps

    # Create radial profiles of pressure and temperature
    # Handle temperature array dimensions for radial averaging
    if not is_unstructured:
        # For radialAvg3D_unstruc, we need 1D lat/lon arrays and T as (lev, ncol)
        T_for_radial = packing.latlon_to_ncol(T)
    else:
        # T is (lev, ncol) - unstructured case
        T_for_radial = T

    # Find the actual minimum pressure location
    cen_lat_actual, cen_lon_actual, _ = find_actual_center(psminlat, psminlon, lat_flat, lon_flat, mps_flat, search_radius=truePS_scan_radius)

    # Calculate radial averages
    # Both variables get keys: radius, radial_average, hit_count
    # 2D vars like psl are nrad
    # 3D vars like Tavg are nlev x nrad
    rad_psl = radialAvg2D_unstruc(mps_flat, lat_flat, lon_flat, deltaMax, cen_lat_actual, cen_lon_actual, rad_for_corr, False)
    Tavg = radialAvg3D_unstruc(T_for_radial, lat_flat, lon_flat, lev_mb, deltaMax, cen_lat_actual, cen_lon_actual, rad_for_corr, False)

    # Restrict to valid pressure levels (1000 > p > 50 mb)
    valid_plev_mask = (lev_mb <= 950.0) & (lev_mb >= 70.0)
    lev_mb = lev_mb[valid_plev_mask]
    Tavg['radial_average'] = Tavg['radial_average'][valid_plev_mask, :]

#     print(lat_flat)
#     print(lon_flat)
#     print(mps_flat)
#     print(T_for_radial[0,:])
#     print("----")
#
#     print(rad_psl['radial_average'])
#     print(rad_psl['radius'])
#     print(rad_psl['hit_count'])
#
#     print(type(Tavg['radial_average']))
#     print(Tavg['radial_average'].shape)

    logging.info(f"Radial pressure profile: min={np.min(rad_psl['radial_average']):.0f} Pa, max={np.max(rad_psl['radial_average']):.0f} Pa")

    # Calculate temperature anomaly (center minus outer radius)
    Tenv = Tavg['radial_average'][:, -1]
    Tanom = Tavg['radial_average'] - Tenv[:, np.newaxis]

#     print(Tanom)
#     print(lev_mb)
#     print(Tenv)
#     print("***")

    # Convert radius to degrees longitude along the latitude line
    mylon = rad_psl['radius'] / (111.0 * np.cos(np.radians(cen_lat_actual)))
    psl_amb = rad_psl['radial_average'][-1]  # Ambient pressure at outer radius
    logging.info(f"Ambient pressure at outer radius is: {psl_amb}")

    # Set up empty arrays to hold analysis estimates for curve fitting
    anlPS = np.full(len(mylon), 100000.0)
    anlT = np.full((len(lev_mb), len(mylon)), 0.0)

    # Set up optimization parameter ranges
    rp_arr = np.linspace(tc_settings['rp_min'], tc_settings['rp_max'], n_p_steps)
    dp_arr = np.linspace(tc_settings['dp_min'], tc_settings['dp_max'], n_p_steps)
    exppr_arr = np.linspace(tc_settings['exppr_min'], tc_settings['exppr_max'], n_p_steps)

    # OPTIMIZATION ROUND 1 - Broad search for pressure parameters
    logging.info("Starting TC parameter optimization round 1 (broad search)...")
    rmsd = np.zeros((4, n_p_steps**3))
    counter = 0

    # These are the "base parameters" that are used for dummy profiles
    # Assume calm atmosphere
    base_zp = 10000.0
    plev_for_optim = lev_mb[-1]
    tin_for_optim = Tenv[-1]
    uin_for_optim = 0.0
    vin_for_optim = 0.0

    for aa, rp in enumerate(rp_arr):
        for bb, dp in enumerate(dp_arr):
            for cc, exppr in enumerate(exppr_arr):

                # Calculate analytic pressure profile along the latitude line
                for ii, lon_pt in enumerate(mylon):
                    tmp = tctestcase(0.0, cen_lat_actual, dp, rp, base_zp, exppr, gamma_,
                                   lon_pt, cen_lat_actual, plev_for_optim * 100., -999, 0,
                                   psl_amb, uin_for_optim, vin_for_optim, tin_for_optim,
                                   q_for_optim(plev_for_optim), False, modify_q, modify_q_mult)
                    anlPS[ii] = tmp[4]  # Get out only sfc pressure

                # Calculate Root Mean Square Deviation
                rmsd[0, counter] = np.sqrt(np.mean((anlPS - rad_psl['radial_average'])**2))
                rmsd[1, counter] = rp
                rmsd[2, counter] = dp
                rmsd[3, counter] = exppr
                counter += 1

    # Find best configuration from round 1
    best_idx = np.argmin(rmsd[0, :])   # best "run" by RMSD
    best_rp, best_dp, best_exppr = rmsd[1, best_idx], rmsd[2, best_idx], rmsd[3, best_idx]
    logging.info(f"Round 1 best: RMSD={rmsd[0, best_idx]:.2f} Pa, rp={best_rp:.0f} m, dp={best_dp:.0f} Pa, exppr={best_exppr:.3f}")

    # OPTIMIZATION ROUND 2 - Refined search
    logging.info("Starting TC parameter optimization round 2 (refined search)...")
    mult_factor = tc_settings['mult_factor']
    rp_range = (rp_arr[1] - rp_arr[0]) * mult_factor
    dp_range = (dp_arr[1] - dp_arr[0]) * mult_factor
    exppr_range = (exppr_arr[1] - exppr_arr[0]) * mult_factor

    # Create new arrays centered on round 1 best with smaller range
    rp_arr = np.linspace(best_rp - rp_range, best_rp + rp_range, n_p_steps)
    dp_arr = np.linspace(best_dp - dp_range, best_dp + dp_range, n_p_steps)
    exppr_arr = np.linspace(best_exppr - exppr_range, best_exppr + exppr_range, n_p_steps)

    rmsd = np.zeros((4, n_p_steps**3))
    counter = 0

    for aa, rp in enumerate(rp_arr):
        for bb, dp in enumerate(dp_arr):
            for cc, exppr in enumerate(exppr_arr):

                for ii, lon_pt in enumerate(mylon):
                    tmp = tctestcase(0.0, cen_lat_actual, dp, rp, base_zp, exppr, gamma_,
                                   lon_pt, cen_lat_actual, plev_for_optim * 100., -999, 0,
                                   psl_amb, uin_for_optim, vin_for_optim, tin_for_optim,
                                   q_for_optim(plev_for_optim/100.), False, modify_q, modify_q_mult)
                    anlPS[ii] = tmp[4]

                rmsd[0, counter] = np.sqrt(np.mean((anlPS - rad_psl['radial_average'])**2))
                rmsd[1, counter] = rp
                rmsd[2, counter] = dp
                rmsd[3, counter] = exppr
                counter += 1

    # Final best pressure parameters from round 2
    best_idx = np.argmin(rmsd[0, :])
    final_rp, final_dp, final_exppr = rmsd[1, best_idx], rmsd[2, best_idx], rmsd[3, best_idx]
    logging.info(f"Round 2 best: RMSD={rmsd[0, best_idx]:.2f} Pa, rp={final_rp:.0f} m, dp={final_dp:.0f} Pa, exppr={final_exppr:.3f}")

    # Now use the final versions of rp, dp, and exppr to estimate zp, which controls warm core height

    # OPTIMIZATION ROUND 3 - Temperature parameter (zp)
    logging.info("Starting temperature parameter optimization (zp)...")
    zp_arr = np.linspace(tc_settings['zp_min'], tc_settings['zp_max'], n_t_steps)
    corr = np.zeros((2, n_t_steps))

    # Loop over the possible zps
    for dd, zp in enumerate(zp_arr):
        # Loop over lons + lats and build analytic T
        for ii, lon_pt in enumerate(mylon):
            for jj, lev_pt in enumerate(lev_mb):
                tmp = tctestcase(0.0, cen_lat_actual, final_dp, final_rp, zp, final_exppr, gamma_,
                               lon_pt, cen_lat_actual, lev_pt * 100.0, -999, 0,
                               psl_amb, uin_for_optim, vin_for_optim, Tenv[jj], q_for_optim(jj),
                               False, modify_q, modify_q_mult)
                anlT[jj, ii] = tmp[3]  # Temperature

        # Calculate temperature anomaly from analytic profile
        anlT_anom = anlT[:, 0] - anlT[:, -1]  # Center minus outer

        # Calculate correlation between analytic and observed temperature anomalies
        corr[0, dd] = np.corrcoef(anlT_anom, Tanom[:,0])[0, 1]
        corr[1, dd] = zp

    # Find best configuration for T anomaly
    if debug:
        nlev = corr.shape[1]  # Get the number of levels
        for i in range(nlev):
            print(f"{corr[0, i]:.6f} {corr[1, i]:.2f}")

    # Find best zp (highest correlation)
    best_idx = np.argmax(corr[0, :])
    final_zp = corr[1, best_idx]
    logging.info(f"Best temperature fit: correlation={corr[0, best_idx]:.3f}, zp={final_zp:.0f} m")

    # Update only the 6 optimized parameters in tc_settings dict
    tc_settings['rp'] = final_rp
    tc_settings['dp'] = final_dp
    tc_settings['zp'] = final_zp
    tc_settings['exppr'] = final_exppr
    tc_settings['psminlat'] = cen_lat_actual
    tc_settings['psminlon'] = cen_lon_actual

    logging.info("FINAL OPTIMIZED TC PARAMETERS:")
    logging.info(f"  rp (radius parameter): {final_rp:.0f} m")
    logging.info(f"  dp (pressure depth): {final_dp:.0f} Pa")
    logging.info(f"  zp (height parameter): {final_zp:.0f} m")
    logging.info(f"  exppr (radial exponent): {final_exppr:.3f}")
    logging.info(f"  Center: lat={cen_lat_actual:.3f}°, lon={cen_lon_actual:.3f}°")

    if plot_bestfits:

        # For final analysis, we take the best values and create final anlT and anlPS arrays
        # for plotting purposes

        # Iterate over radial points and levels
        for ii, lon_pt in enumerate(mylon):
            for jj, lev_pt in enumerate(lev_mb):
                # Call the tctestcase function and unpack the result into anlT and anlPS
                tmp = tctestcase(0.0, cen_lat_actual, final_dp, final_rp, final_zp, final_exppr, gamma_,
                               lon_pt, cen_lat_actual, lev_pt * 100.0, -999, 0,
                               psl_amb, 0.0, 0.0, Tenv[jj], q_for_optim(jj),
                               False, modify_q, modify_q_mult)

                # Assign the temperature at level jj and radius ii
                anlT[jj, ii] = tmp[3]  # Equivalent to tmp(3) in NCL (temperature)

            # Assign the pressure at radius ii
            anlPS[ii] = tmp[4]  # Equivalent to tmp(4) in NCL (pressure)

        # Compute the temperature anomaly
        anlT_anom = anlT[:, 0] - anlT[:, -1]  # Subtract the first column from the last column

        # Plot data
        plot_best_fit_profiles(
            radius=rad_psl['radius'],
            rad_psl_avg=rad_psl['radial_average'],
            anlPS=anlPS,
            Tanom=Tanom[:,0],
            anlT_anom=anlT_anom,
            lev=lev_mb
        )

    return tc_settings


def plot_best_fit_profiles(radius, rad_psl_avg, anlPS, Tanom, anlT_anom, lev, output_path='best_fit_profiles.png', show=False):

    fig, (ax1, ax2) = plt.subplots(nrows=2, figsize=(8, 12))

    # Pressure profile
    ax1.plot(radius, rad_psl_avg, label='Radial Avg Pressure', color='black')
    ax1.plot(radius, anlPS, label='Analytic Pressure', linestyle='--', color='red')
    ax1.set_title('Radial Pressure Profile')
    ax1.set_xlabel('Radius (km)')
    ax1.set_ylabel('Surface pressure')
    ax1.legend()
    ax1.set_xlim(0, 900)
    ax1.set_ylim(93000, 101500)

    # Temperature anomaly
    ax2.plot(Tanom, lev, label='Temperature Anomaly', color='black')
    ax2.plot(anlT_anom, lev, linestyle='--', color='red')
    ax2.invert_yaxis()
    ax2.set_title('Temperature Anomaly Profile')
    ax2.set_xlabel('Temperature Anomaly (K)')
    ax2.set_ylabel('hybrid level at midpoints (1000*(A+B))')
    ax2.set_xlim(-10, 10)
    ax2.set_ylim(1010, 40)

    plt.tight_layout()

    if show:
        plt.show()
    else:
        plt.savefig(output_path)
        print(f"Plot saved to {output_path}")


def apply_tc_seeding(data_vars, tc_params, update_z=False, add_deltas_to_data=False, min_rmw=15000.0):
    """
    Apply TC seeding/unseeding to atmospheric data using optimized parameters.

    Parameters:
    -----------
    data_vars : dict
        Dictionary containing atmospheric fields with keys: 'ps', 't', 'u', 'v', 'q', 'lat', 'lon', 'lev'
    tc_params : dict
        TC parameters from read_tcseed_settings/find_fill_parameters
    add_deltas_to_data : bool
        XXXXX

    Returns:
    --------
    data_vars : dict
        Modified data_vars with TC seeding applied
    """

    # If update_z is True, we need to make sure we have z and pint in data_vars before going further
    if update_z:
        required_keys = ['z', 'pint']
        missing_keys = [key for key in required_keys if key not in data_vars]
        if missing_keys:
            print(f"Error: correct_z=True but missing required keys: {missing_keys}")
            sys.exit(1)

    # Extract TC parameters
    invert_vortex = tc_params.get('invert_vortex')
    modify_q = tc_params.get('modify_q')
    modify_q_mult = tc_params.get('modify_q_mult')
    gamma_ = tc_params.get('gamma_')
    cen_lat = tc_params.get('psminlat')
    cen_lon = tc_params.get('psminlon')
    zp = tc_params.get('zp')
    exppr = tc_params.get('exppr')
    deg_bnd = tc_params.get('degree_adjust')

    if not invert_vortex:
        minp = tc_params.get('minp')
        target_rmw = tc_params.get('target_rmw')

        # If minp is really small, set to 995 mb
        if minp < -19999.0:
            minp = 995.0
        # If minp is negative, but plausible
        elif 0.0 < minp <= -19999.0:
            print("Using negative minp as dp")

        # If target_rmw is less than 0, set to 200km
        if target_rmw < 0.0:
            logging.info("Negative RMW passed in, setting to 200km:")
            target_rmw = 200000.0
        # If target_rmw is less than 0, set to 200km
        if target_rmw < min_rmw:
            logging.info(f"Target RMW requested {target_rmw} is less than min_rmw, capping at {min_rmw}")
            target_rmw = min_rmw
    else:
        rp = tc_params.get('rp')
        dp = tc_params.get('dp')

    # Get coordinate and data arrays and count some things
    lat = data_vars['lat']
    lon = data_vars['lon']
    lev = data_vars['lev']
    nlev = len(lev)

    # Handle different grid structures (structured vs unstructured)
    if lat.ndim == 1 and lon.ndim == 1 and len(lat) != len(lon):
        # Structured grid: lat(nlat), lon(nlon), ps(nlat, nlon), etc.
        nlat, nlon = len(lat), len(lon)
        ncol = nlat * nlon
        is_structured = True

        # Create 2D coordinate arrays and flatten
        lat_work, lon_work = packing.flatten_1D_lat_lon(lat, lon)

        ps = packing.latlon_to_ncol(data_vars['ps'])
        u = packing.latlon_to_ncol(data_vars['u'])
        v = packing.latlon_to_ncol(data_vars['v'])
        t = packing.latlon_to_ncol(data_vars['t'])
        q = packing.latlon_to_ncol(data_vars['q'])
        if update_z:
            logging.info("z updates request in seeding code, loading z and pint")
            z = packing.latlon_to_ncol(data_vars['z'])
            pint = packing.latlon_to_ncol(data_vars['pint'])
    else:
        # Unstructured grid: lat(ncol), lon(ncol), ps(ncol), etc.
        ncol = len(lat)
        is_structured = False

        # Work directly with the data
        lat_work = lat
        lon_work = lon

        ps = data_vars['ps'].copy()
        u = data_vars['u'].copy()
        v = data_vars['v'].copy()
        t = data_vars['t'].copy()
        q = data_vars['q'].copy()
        if update_z:
            z = data_vars['z'].copy()

    # Calculate great circle distances
    gcdist, _ = gc_latlon(cen_lat, cen_lon, lat_work, lon_work, 2, 2)

    # Find relevant parameters for TC seed
    if not invert_vortex:
        if minp > 0.0:
            minix = np.argmin(gcdist)
            logging.info(f'minix: {minix} --> {lat_work[minix]}, {lon_work[minix]}')
            # if the user gave us minp, we subtract from ambient at vortex center
            ambps = ps[minix]
            logging.info(f"ambient ps at minix: {ambps}")
            dp = ambps - (minp * 100.0)
            # if dp here is negative, it means the requested vortex is shallower than existing
            if dp <= 0.0:
                logging.info(f"Negative dp, no adjustment, returning")
                logging.info(f"If you want to fill a vortex, toggle invert_vortex")
                return data_vars
        else:
            # If the user gave us a negative minp to begin with, it's actually dp
            dp = -minp
        rp = get_rp_from_dp_rmw(cen_lat, dp, target_rmw)
        logging.info(f"VORTEX DERIVED: rp: {rp}  -- dp: {dp}")

    # Create copies of unmodified fields to query later if needed
    u_orig = np.copy(u)
    v_orig = np.copy(v)
    ps_orig = np.copy(ps)
    t_orig = np.copy(t)
    q_orig = np.copy(q)
    if update_z:
        z_orig = np.copy(z)

    # Add the vortex seed to the state fields
    for ii in range(ncol):
        if gcdist[ii] <= deg_bnd:
            if update_z:

                # Estimate the vertical column at this ii before
                z_before = meteo.hydro_bottom_to_top(ps[ii], pint[:,ii], t[:, ii], q[:, ii])

            for kk in range(nlev):
                #p = hyam[kk] * P0 + hybm[kk] * ps[ii]
                # Note, -999 is a failsafe for "z" since zcoords = should be pressure.
                theArr = tctestcase(cen_lon, cen_lat, dp, rp, zp, exppr, gamma_, lon_work[ii], lat_work[ii], lev[kk], -999, 0, ps[ii], u[kk, ii], v[kk, ii], t[kk, ii], q[kk, ii], invert_vortex, modify_q, modify_q_mult)
                v[kk, ii] = theArr[0]
                u[kk, ii] = theArr[1]
                q[kk, ii] = theArr[2]
                t[kk, ii] = theArr[3]
                if kk == 0:
                    # Just update ps once
                    ps[ii] = theArr[4]

            if update_z:
                # Estimate the vertical column after using new ps, t, q
                z_after = meteo.hydro_bottom_to_top(ps[ii], pint[:,ii], t[:, ii], q[:, ii])
                # Calculate a delta from the estimates and apply that to the input z
                # Probably safest to do it this way in case the hydro approx isn't perfect
                z[:,ii] = z[:,ii] + (z_after - z_before)

    if is_structured:
        # Convert back to lat/lon structure
        data_vars['ps'] = packing.ncol_to_latlon(ps, nlat, nlon)
        data_vars['u'] = packing.ncol_to_latlon(u, nlat, nlon)
        data_vars['v'] = packing.ncol_to_latlon(v, nlat, nlon)
        data_vars['t'] = packing.ncol_to_latlon(t, nlat, nlon)
        data_vars['q'] = packing.ncol_to_latlon(q, nlat, nlon)
        if update_z:
            data_vars['z'] = packing.ncol_to_latlon(z, nlat, nlon)

        if add_deltas_to_data:
            # Initialize or accumulate delta variables
            # First storm: data_vars.get('ps_vx', 0) returns 0 (default), so ps_vx = 0 + ps_delta (initialization)
            # Subsequent TCs: data_vars.get('ps_vx', 0) returns the existing value, so ps_vx = existing_value + ps_delta (accumulation)
            data_vars['ps_vx'] = data_vars.get('ps_vx', 0) + packing.ncol_to_latlon((ps - ps_orig), nlat, nlon)
            data_vars['t_vx'] = data_vars.get('t_vx', 0) + packing.ncol_to_latlon((t - t_orig), nlat, nlon)
            data_vars['u_vx'] = data_vars.get('u_vx', 0) + packing.ncol_to_latlon((u - u_orig), nlat, nlon)
            data_vars['v_vx'] = data_vars.get('v_vx', 0) + packing.ncol_to_latlon((v - v_orig), nlat, nlon)
            data_vars['q_vx'] = data_vars.get('q_vx', 0) + packing.ncol_to_latlon((q - q_orig), nlat, nlon)
            if update_z:
                data_vars['z_vx'] = data_vars.get('z_vx', 0) + packing.ncol_to_latlon((z - z_orig), nlat, nlon)

    else:
        # Keep unstructured format
        data_vars['ps'] = ps
        data_vars['u'] = u
        data_vars['v'] = v
        data_vars['t'] = t
        data_vars['q'] = q
        if update_z:
            data_vars['z'] = z

        if add_deltas_to_data:
            data_vars['ps_vx'] = data_vars.get('ps_vx', 0) + (ps - ps_orig)
            data_vars['t_vx'] = data_vars.get('t_vx', 0) + (t - t_orig)
            data_vars['u_vx'] = data_vars.get('u_vx', 0) + (u - u_orig)
            data_vars['v_vx'] = data_vars.get('v_vx', 0) + (v - v_orig)
            data_vars['q_vx'] = data_vars.get('q_vx', 0) + (q - q_orig)
            if update_z:
                data_vars['z_vx'] = data_vars.get('z_vx', 0) + (z - z_orig)

    logging.info("TC seeding completed successfully")

    return data_vars


def parse_tcvitals(filename, datetime_input):
    """
    Parse TCVitals file and return dictionary of lists for cyclones matching the given date/time.

    Parameters:
    -----------
    filename : str
        Path to the TCVitals file
    datetime_input : str or int
        Date and time in YYYYMMDDHH format (e.g., "2002090106" or 2002090106 for Sept 1, 2002 06:00)

    Returns:
    --------
    cyclone_data : dict
        Dictionary with keys 'lats', 'lons', 'pressures', 'rmw', 'storm_names'
        Each value is a list of the corresponding data for all matching cyclones
    """

    # Convert to string if input is integer
    if isinstance(datetime_input, int):
        datetime_str = f"{datetime_input:010d}"  # Ensure 10 digits with leading zeros
    elif isinstance(datetime_input, str):
        datetime_str = datetime_input
    else:
        raise TypeError(f"datetime_input must be int or str, got {type(datetime_input)}")

    # Validate and parse the datetime string
    if len(datetime_str) != 10:
        raise ValueError(f"datetime_input must be 10 digits (YYYYMMDDHH), got: {datetime_str}")

    try:
        yyyy = datetime_str[:4]
        mm = datetime_str[4:6]
        dd = datetime_str[6:8]
        hh = datetime_str[8:10]

        # Validate ranges
        if not (1 <= int(mm) <= 12):
            raise ValueError(f"Invalid month: {mm}")
        if not (1 <= int(dd) <= 31):
            raise ValueError(f"Invalid day: {dd}")
        if not (0 <= int(hh) <= 23):
            raise ValueError(f"Invalid hour: {hh}")

    except ValueError as e:
        raise ValueError(f"Invalid datetime format '{datetime_str}': {e}")

    # Format the target date and time for matching
    target_date = f"{yyyy}{mm}{dd}"
    target_time = f"{hh}00"

    # Lists to store extracted data
    lats = []
    lons = []
    pressures = []
    rmw_values = []
    storm_names = []

    # Check if file exists
    if not os.path.exists(filename):
        raise FileNotFoundError(f"TCVitals file not found: {filename}")

    # Parse the file
    open_func = gzip.open if filename.endswith('.gz') else open
    with open_func(filename, 'rt') as f:
        for line_num, line in enumerate(f, 1):
            line = line.strip()
            if not line:
                continue

            # Split the line into fields
            fields = line.split()

            # Check if we have enough fields (minimum expected)
            if len(fields) < 15:
                print(f"Warning: Line {line_num} has insufficient fields, skipping")
                continue

            try:
                # Extract date and time (fields 3 and 4)
                file_date = fields[3]
                file_time = fields[4]

                # Check if this record matches our target date/time
                if file_date == target_date and file_time == target_time:

                    storm_ok = True

                    # Extract storm information
                    storm_id = fields[1]
                    storm_name = fields[2]

                    # Extract latitude (field 5) - format like "249N" or "167N"
                    lat_str = fields[5]
                    if lat_str.endswith('N'):
                        lat = float(lat_str[:-1]) / 10.0  # Convert to decimal degrees
                    elif lat_str.endswith('S'):
                        lat = -float(lat_str[:-1]) / 10.0  # South is negative
                    else:
                        print(f"Warning: Unexpected latitude format '{lat_str}' on line {line_num}")
                        storm_ok = False

                    # Extract longitude (field 6) - format like "1270W" or "1101W"
                    lon_str = fields[6]
                    if lon_str.endswith('W'):
                        lon = -float(lon_str[:-1]) / 10.0  # West is negative
                    elif lon_str.endswith('E'):
                        lon = float(lon_str[:-1]) / 10.0  # East is positive
                    else:
                        print(f"Warning: Unexpected longitude format '{lon_str}' on line {line_num}")
                        storm_ok = False

                    # Extract pressure (field 9)
                    pressure_str = fields[9]
                    if pressure_str == '-999' or pressure_str == '-99':
                        pressure = -999  # Missing data
                        storm_ok = False
                    else:
                        pressure = int(pressure_str)

                    # Extract RMW (field 13)
                    # Estimated radius at which the maximum wind occurs in kilometers from the storm center
                    rmw_str = fields[13]
                    if rmw_str == '-999' or rmw_str == '-99':
                        rmw = -999  # Missing data
                        storm_ok = False
                    else:
                        rmw = int(rmw_str)

                    # Store the extracted data
                    if storm_ok:
                        lats.append(lat)
                        lons.append(lon)
                        pressures.append(pressure)
                        rmw_values.append(rmw)
                        storm_names.append(f"{storm_id}_{storm_name}")

            except (ValueError, IndexError) as e:
                print(f"Warning: Error parsing line {line_num}: {e}")
                continue

    # Handle duplicate storm names - keep only the one with minimum pressure
    if len(storm_names) > 1:
        # Create a dictionary to track storm names and their indices
        storm_indices = {}
        for i, name in enumerate(storm_names):
            if name not in storm_indices:
                storm_indices[name] = []
            storm_indices[name].append(i)

        # Find indices to keep (one per unique storm name with minimum pressure)
        indices_to_keep = []
        for name, indices in storm_indices.items():
            if len(indices) == 1:
                # No duplicates for this storm
                indices_to_keep.append(indices[0])
            else:
                # Find the index with minimum pressure
                min_pressure = float('inf')
                best_index = indices[0]
                for idx in indices:
                    if pressures[idx] < min_pressure:
                        min_pressure = pressures[idx]
                        best_index = idx
                indices_to_keep.append(best_index)
                print(f"Info: Found {len(indices)} duplicates for storm '{name}', keeping the one with pressure {min_pressure} mb")

        # Sort indices to maintain order
        indices_to_keep.sort()

        # Create new lists with only the kept storms
        lats = [lats[i] for i in indices_to_keep]
        lons = [lons[i] for i in indices_to_keep]
        pressures = [pressures[i] for i in indices_to_keep]
        rmw_values = [rmw_values[i] for i in indices_to_keep]
        storm_names = [storm_names[i] for i in indices_to_keep]

    # Return as dictionary of lists
    cyclone_data = {
        'lats': lats,
        'lons': lons,
        'pressures': pressures,
        'rmw': rmw_values,
        'storm_names': storm_names
    }

    return cyclone_data


def print_tcvitals_summary(filename, datetime_input):
    """
    Print a summary of cyclones found for the given date/time in a TCVitals file
    """
    cyclone_data = parse_tcvitals(filename, datetime_input)

    # Convert to string for parsing if needed
    if isinstance(datetime_input, int):
        datetime_str = f"{datetime_input:010d}"
    else:
        datetime_str = datetime_input

    # Parse datetime for display
    yyyy, mm, dd, hh = datetime_str[:4], datetime_str[4:6], datetime_str[6:8], datetime_str[8:10]

    print(f"Found {len(cyclone_data['lats'])} cyclones for {yyyy}-{mm}-{dd} {hh}:00")
    print(f"Storm IDs: {', '.join(cyclone_data['storm_names'])}")
    print()

    for i in range(len(cyclone_data['lats'])):
        print(f"Storm {i+1}: {cyclone_data['storm_names'][i]}")
        print(f"  Lat: {cyclone_data['lats'][i]}°")
        print(f"  Lon: {cyclone_data['lons'][i]}°")
        print(f"  Pressure: {cyclone_data['pressures'][i]} mb")
        print(f"  RMW: {cyclone_data['rmw'][i]} km")
        print()


def safe_exp(arg, max_arg=300):
    """Safely evaluate exp(arg), clamping to avoid overflow."""
    if not np.isfinite(arg):
        #logging.info(f"safe_exp received non-finite input: {arg}")
        return np.exp(max_arg)
    if arg > max_arg:
        #logging.info(f"Clamping exp({arg:.2f}) to exp({max_arg}) to avoid overflow")
        return np.exp(max_arg)
    return np.exp(arg)


def calculate_wind_components_unstructured(U, V, lat, lon, psminlat, psminlon, movespd=0, movedir=0):
    """
    Calculate radial and tangential wind components relative to storm center on unstructured grid.

    Parameters:
    -----------
    U, V : array (ncol, nlev) or (ncol,)
        Zonal and meridional wind components (m/s)
    lat, lon : array (ncol,)
        Latitude and longitude coordinates (degrees)
    psminlat, psminlon : float
        Storm center coordinates (degrees)
    movespd : float, optional
        Storm movement speed in m/s (default: 0)
    movedir : float, optional
        Storm movement direction in meteorological degrees
        (0=N, 90=E, 180=S, 270=W, default: 0)

    Returns:
    --------
    v_rad, v_theta : array
        Radial and tangential wind components (same shape as U, V)
        v_rad: positive = outflow from center
        v_theta: positive = counterclockwise rotation
    """
    import numpy as np

    # Constants
    pi = np.pi
    d2r = pi / 180.0
    r2d = 180.0 / pi

    # Convert coordinates to radians
    lonr = lon * d2r
    latr = lat * d2r
    psminlat_r = psminlat * d2r
    psminlon_r = psminlon * d2r

    # Calculate longitude difference
    deltalong = lonr - psminlon_r

    # Calculate bearing angle from storm center to each grid point
    # This uses the standard bearing calculation formula
    arr1 = np.sin(deltalong) * np.cos(latr)
    arr2 = (np.cos(psminlat_r) * np.sin(latr) -
            np.sin(psminlat_r) * np.cos(latr) * np.cos(deltalong))
    dir_angle_r = np.arctan2(arr1, arr2)

    # Calculate wind direction angle
    phi_r = np.arctan2(U, V)

    # Make copies for manipulation
    U_adj = U.copy()
    V_adj = V.copy()

    # Remove storm motion if specified
    if movespd > 0:
        print(f"Adjusting U and V fields to remove TC forward motion: {movespd} m/s at {movedir}°")
        U_adj = U_adj - movespd * np.sin(movedir * d2r)
        V_adj = V_adj - movespd * np.cos(movedir * d2r)

    # Calculate wind speed magnitude
    WIND = np.sqrt(U_adj**2 + V_adj**2)

    # Handle broadcasting for 3D case (ncol, nlev)
    if U.ndim == 2:
        # Broadcast direction angle to match 3D wind data
        dir_angle_r_bc = dir_angle_r[:, np.newaxis]
        # Transform to radial/tangential components
        v_theta = WIND * np.sin(dir_angle_r_bc - phi_r)
        v_rad = -WIND * np.cos(dir_angle_r_bc - phi_r)
    else:
        # 2D case
        v_theta = WIND * np.sin(dir_angle_r - phi_r)
        v_rad = -WIND * np.cos(dir_angle_r - phi_r)

    return v_rad, v_theta
