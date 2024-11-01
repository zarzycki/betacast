import numpy as np
from math import sin, cos, sqrt, pi, atan, exp, radians, degrees
import re


def calc_mass_weighted_integral(var3d, pdel, area, gravit=9.81):
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

def calc_inverse_mass_weighted_field(target_integral, pdel, area, gravit=9.81):
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
    with open(file_path, 'r') as file:
        lines = file.readlines()

    # Eliminate any leading or trailing spaces in the delete command
    re_pattern = re.compile(rf"^\s*{variable}\s*=")
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



def tctestcase(cen_lon, cen_lat, dp, rp, zp, exppr, gamma_, lon, lat, p, z, zcoords, psin, uin, vin, Tin, qin, invert_vortex, modify_q, modify_q_mult, debug=False):
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
            p = ptrop * exp(-g * (height - ztrop) / (Rd * Ttrop))
        else:
            p = (p00 - dp * exp(-(gr / rp) ** exppr) * exp(-(height / zp) ** exppz)) * \
                ((T0 - gamma_ * height) / T0) ** (1 / exponent)
    else:
        ps = -dp * exp(-(gr / rp) ** exppr) + p00
        height = (T0 / gamma_) * (1. - (p / ps) ** exponent)

    # Initialize U and V (wind components)
    d1 = sin(radians(cen_lat)) * cos(radians(lat)) - cos(radians(cen_lat)) * sin(radians(lat)) * cos(radians(lon) - radians(cen_lon))
    d2 = cos(radians(cen_lat)) * sin(radians(lon) - radians(cen_lon))
    d = max(epsilon, sqrt(d1 ** 2. + d2 ** 2.))

    ufac = d1 / d
    vfac = d2 / d

    if height > ztrop:
        u = uin
        v = vin
    else:
        if debug:
            print(f"height: {height}, zp: {zp}, exppz: {exppz}, dp: {dp}, rp: {rp}, gr: {gr}")
        vt = (-f * gr / 2 + sqrt((f * gr / 2) ** 2 - (exppr * (gr / rp) ** exppr) *
                                 (Rd * (T0 - gamma_ * height)) / (exppz * height * Rd * (T0 - gamma_ * height) / (g * zp ** exppz) + 1. - p00 / dp * exp((gr / rp) ** exppr) * exp((height / zp) ** exppz))))
        v = vin + vfac * vt
        u = uin + ufac * vt

    # Calculate temperature and specific humidity
    Tvin = Tin * (1 + constTv * qin)
    t = (Tvin) / (1 + constTv * qin) / (1 + exppz * Rd * Tvin * height / (g * zp ** exppz * (1 - p00 / dp * exp((gr / rp) ** exppr) * exp((height / zp) ** exppz))))

    if modify_q:
        rh = 0.263 * p * qin / exp((17.67 * (Tin - 273.15)) / (Tin - 29.65))
        if height > ztrop:
            q = qin
        else:
            q = rh / (0.263 * p) * exp((17.67 * (t - 273.15)) / (t - 29.65))
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
            T1 = -(f * r) / 2.
            T2 = (f ** 2 * r ** 2) / 4.
            NUM = (3. / 2.) * ((r / rp) ** (3. / 2.)) * T0 * Rd
            DEN = 1. - (p00 / dp) * exp((r / rp) ** (3. / 2.))
            vt[ii] = T1 + sqrt(T2 - (NUM / DEN))

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
            exit()

    print(f"dp={dp}, rp={rpi}")
    return rpi


def keyword_values(namelist_file, key, return_type):
    """
    Reads a namelist file and extracts the value for the specified key.
    Parameters:
    namelist_file (str): Path to the namelist file.
    key (str): The key whose value is to be extracted.
    return_type (str): The type to which the value should be converted. Options: "int", "float", "bool", "str".
    Returns:
    The value associated with the key in the specified return type.
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
                        return int(v)
                    elif return_type == "float":
                        return float(v)
                    elif return_type == "bool":
                        return v.lower() in ['true', 't', '1']
                    elif return_type == "str":
                        return v

    # If key is not found
    raise KeyError(f"Key '{key}' not found in the namelist file.")


def radialAvg2D_unstruc(data, lat, lon, deltaMax, psminlat, psminlon, outerRad, mergeInnerBins):
    # km per degree at the equator
    kmInDeg = 111.32
    kmGrid = kmInDeg * deltaMax
    print(f"The max lat/lon km grid spacing at equator is {kmGrid} km")

    ncol = len(lat)
    pi = np.pi
    d2r = pi / 180.0

    # Prepare radial bins
    timesGrid = 1.1  # We want each radius bin to be timesGrid * kmGrid
    nx = int(outerRad / (timesGrid * kmGrid))
    print(f"Number of bins is equal to {nx}")

    if mergeInnerBins:
        numMerge = 2  # Number of innermost radial bins to merge
        print(f"Merging innermost {numMerge} bins because radial average is so small.")
        numMergeMinusOne = numMerge - 1
        origRadiusArr = np.linspace(0, outerRad, nx + numMergeMinusOne)
        radiusArr = np.zeros(len(origRadiusArr) - numMergeMinusOne)
        radiusArr[0] = origRadiusArr[0]
        radiusArr[1:] = origRadiusArr[numMerge:]
    else:
        print("Not merging any innermost bins -- be careful that your inner bins have > 1 pt.")
        radiusArr = np.linspace(0, outerRad, nx)

    numRadBins = len(radiusArr)

    rad_thevar_hit = np.zeros(numRadBins, dtype=int)
    rad_thevar_cum = np.zeros(numRadBins, dtype=float)

    # Iterate over each point in the dataset
    print("Starting loop")
    for i in range(ncol):
        # Use the gc_latlon function to calculate the great circle distance
        gcdist, _ = gc_latlon(psminlat, psminlon, lat[i], lon[i], 2, 4)

        if gcdist <= outerRad:
            bin_idx = np.argmin(np.abs(radiusArr - gcdist))
            rad_thevar_hit[bin_idx] += 1
            rad_thevar_cum[bin_idx] += data[i]

    print(f"Minimum number of hits per gridbox: {np.min(rad_thevar_hit)}")
    print(f"Maximum number of hits per gridbox: {np.max(rad_thevar_hit)}")

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
    print(f"The max lat/lon km grid spacing at equator is {kmGrid} km")

    ncol = len(lat)
    nlev = len(lev)
    pi = np.pi
    d2r = pi / 180.0

    # Prepare radial bins
    timesGrid = 1.1  # We want each radius bin to be timesGrid * kmGrid
    nx = int(outerRad / (timesGrid * kmGrid))
    print(f"Number of bins is equal to {nx}")

    if mergeInnerBins:
        numMerge = 2  # Number of innermost radial bins to merge
        print(f"Merging innermost {numMerge} bins because radial average is so small.")
        numMergeMinusOne = numMerge - 1
        origRadiusArr = np.linspace(0, outerRad, nx + numMergeMinusOne)
        radiusArr = np.zeros(len(origRadiusArr) - numMergeMinusOne)
        radiusArr[0] = origRadiusArr[0]
        radiusArr[1:] = origRadiusArr[numMerge:]
    else:
        print("Not merging any innermost bins -- be careful that your inner bins have > 1 pt.")
        radiusArr = np.linspace(0, outerRad, nx)

    numRadBins = len(radiusArr)

    rad_thevar_hit = np.zeros((nlev, numRadBins), dtype=int)
    rad_thevar_cum = np.zeros((nlev, numRadBins), dtype=float)

    # Iterate over each column in the dataset
    print("Starting loop")
    for i in range(ncol):
        # Use the gc_latlon function to calculate the great circle distance
        gcdist, _ = gc_latlon(psminlat, psminlon, lat[i], lon[i], 2, 4)

        if gcdist <= outerRad:
            bin_idx = np.argmin(np.abs(radiusArr - gcdist))
            rad_thevar_hit[:, bin_idx] += 1
            rad_thevar_cum[:, bin_idx] += data[:, i]

    # Handle grid boxes with no hits
    rad_thevar_hit = np.where(rad_thevar_hit == 0, np.nan, rad_thevar_hit)

    print(f"Minimum number of hits per gridbox: {np.nanmin(rad_thevar_hit)}")
    print(f"Maximum number of hits per gridbox: {np.nanmax(rad_thevar_hit)}")

    # Calculate the radial average
    rad_thevar = np.divide(rad_thevar_cum, rad_thevar_hit, out=np.zeros_like(rad_thevar_cum), where=rad_thevar_hit != 0)

    # Returning the result as a dictionary for simplicity
    return {
        'radius': radiusArr,
        'radial_average': rad_thevar,
        'hit_count': rad_thevar_hit
    }