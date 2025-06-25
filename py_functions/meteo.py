import numpy as np
import logging
import numba as nb
# import xarray as xr
# import datetime
# import argparse
# import os
# import sys
# import glob
# import cftime
# from scipy.ndimage import gaussian_filter
from constants import (
    p0, grav, Rd, Rv_over_Rd, kappa, mpas_uv_damping_coeffs
)
import py_seedfuncs

@nb.jit(nopython=True)
def hydro_bottom_to_top(ps, pint, T, q, zsfc=0.0):
    """
    Compute midlayer geopotential height via hydrostatic integration from surface pressure upward.

    Parameters
    ----------
    ps : float
        Surface pressure [Pa].
    pint : ndarray
        Interface pressures [Pa], length nlev, ordered top-to-bottom or bottom-to-top.
    T : ndarray
        Temperature [K] at each layer midpoint, shape (nlev,). Ordered same as pint.
    q : ndarray
        Water vapor mixing ratio [kg/kg] at each layer midpoint, shape (nlev,). Ordered same as pint.
    zsfc : float, optional
        Surface geopotential height [m], default is 0.

    Returns
    -------
    zmid : ndarray
        Geopotential height [m] at layer midpoints, shape (nlev,). Preserves input vertical ordering.

    Notes
    -----
    Interface pressures are clamped such that no level extends below the surface pressure.
    """

    nlev = pint.size

    # Flip inputs if they run top-to-bottom
    flip_vert = False
    if pint[0] < pint[nlev-1]:
        flip_vert = True
        pint2 = np.empty(nlev, dtype=pint.dtype)
        T2    = np.empty(nlev, dtype=T.dtype)
        q2    = np.empty(nlev, dtype=q.dtype)
        for i in range(nlev):
            pint2[i] = pint[nlev-1-i]
            T2[i]    = T[nlev-1-i]
            q2[i]    = q[nlev-1-i]
        pint = pint2
        T    = T2
        q    = q2

    # Build full interface pressures, clamp at ps
    p_int = np.empty(nlev+1, dtype=pint.dtype)
    p_int[0] = ps
    for i in range(nlev):
        p_int[i+1] = pint[i] if pint[i] <= ps else ps

    # Hydrostatic integration bottom-to-top
    zh_int = np.empty(nlev+1, dtype=pint.dtype)
    zh_int[0] = zsfc

    # constants (Rd/g)
    RDAG  = Rd / grav

    for k in range(nlev):
        # virtual temp in layer k
        tkv = T[k] * (1.0 + 0.61 * q[k])
        # ΔΦ = RDAG * tkv * ln(p_lower / p_upper)
        zh_int[k+1] = zh_int[k] + RDAG * tkv * np.log(p_int[k] / p_int[k+1])

    zmid = (zh_int[:-1] + zh_int[1:]) / 2

    # Flip result back if needed
    if flip_vert:
        z2mid = zmid[::-1]
        return z2mid

    return zmid




def compute_pmid(t, rho):
    """
    Compute the mid-level pressure (pmid) from temperature (T) and density (rho).

    Parameters:
    - t: numpy array, temperature (can be 2D, 3D, 4D, etc.).
    - rho: numpy array, density (must have the same shape as t).

    Returns:
    - pmid: numpy array, computed pressure with the same shape as t and rho.
    """
    # Check if input arrays have the same shape
    if t.shape != rho.shape:
        raise ValueError("Temperature and density arrays must have the same shape.")

    # Compute pressure using the ideal gas law: pmid = rho * R_dry_air * T
    pmid = rho * Rd * t

    return pmid


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


@nb.jit(nopython=True)
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
    prcwat = 0.0
    valid_layers = 0
    for i in range(len(Q)):
        q, dp = Q[i], DP[i]
        if not np.isnan(q) and not np.isnan(dp):
            if q != QMSG and dp != DPMSG:
                valid_layers += 1
                prcwat += q * abs(dp)
    if valid_layers > 0:
        prcwat /= grav
    else:
        prcwat = QMSG
    return prcwat


def deswskewt(t):
    """
    Python implementation of DESWSKEWT function from NCL.
    Calculates saturation vapor pressure over water.

    Parameters:
    -----------
    t : float or numpy.ndarray
        Temperature in Celsius

    Returns:
    --------
    es : float or numpy.ndarray
        Saturation vapor pressure in mb (hPa)
    """

    # Standard formula for saturation vapor pressure (mb)
    return 6.112 * np.exp(17.67 * t / (t + 243.5))


def mixhum_ptd(p, tdk, iswit=1):
    """
    Calculates the mixing ratio or specific humidity given pressure and dew point temperature.
    Attempt to match NCL's mixhum_ptd function.

    Parameters:
    -----------
    p : float or numpy.ndarray
        Atmospheric pressure (Pa).
    tdk : float or numpy.ndarray
        Dew point temperature (K).
    iswit : int
        - If iswit=1, the output will be the mixing ratio (kg/kg).
        - If iswit=2, the output will be the specific humidity (kg/kg).
        - If iswit is negative, the units will be g/kg.

    Returns:
    --------
    wmr : float or numpy.ndarray
        Mixing ratio or specific humidity in the units specified by iswit.
    """
    # Constants (stolen from NCL)
    T0 = 273.15  # K to C conversion
    PA2MB = 0.01  # Pa to mb conversion
    EPS = 0.62197  # Ratio of molecular weights (water/dry air)

    # Unit conversion
    p_mb = p * PA2MB  # Convert Pa to mb (hPa)
    td_c = tdk - T0   # Convert K to C

    # Implement NCL's DWMRSKEWT calculation
    x = 0.02 * (td_c - 12.5 + 7500.0 / p_mb)
    wfw = 1.0 + 4.5e-6 * p_mb + 1.4e-3 * x * x
    fwesw = wfw * deswskewt(td_c)
    r = EPS * fwesw / (p_mb - fwesw)

    # Convert r to g/kg
    wmr = 1000.0 * r

    # Convert to kg/kg for standard output
    wmr = wmr * 0.001

    # Convert to specific humidity if iswit=2
    if abs(iswit) == 2:
        wmr = wmr / (wmr + 1.0)

    # Convert to g/kg if iswit<0
    if iswit < 0:
        wmr = wmr * 1000.0

    return wmr


@nb.jit(nopython=True)
def ps_wet_to_dry_conversion_core(ps_fv, q_fv, hyai, hybi, p0):
    """Numba-optimized core of the conversion function"""
    ncol = ps_fv.shape[1]
    pw_fv = np.zeros_like(ps_fv)

    for kk in range(ncol):
        pi_orig = hyai * p0 + hybi * ps_fv[0, kk]
        nlevp1 = pi_orig.size
        dp = pi_orig[1:nlevp1] - pi_orig[0:nlevp1 - 1]
        pw_fv[0, kk] = prcwater_dp(q_fv[:, 0, kk], dp)
        ps_fv[0, kk] = ps_fv[0, kk] - pw_fv[0, kk] * grav

    return ps_fv, pw_fv


def ps_wet_to_dry_conversion(ps_fv, q_fv, hyai, hybi, p0, verbose=False):
    """
    Converts wet surface pressure to dry surface pressure by subtracting
    the total column precipitable water.
    Parameters:
    - ps_fv: Surface pressure array (1D or 2D).
    - q_fv: Specific humidity array (2D or 3D, with shape [levels, ...]).
    - hyai: Hybrid A interface coefficients (1D).
    - hybi: Hybrid B interface coefficients (1D).
    - p0: Reference pressure (scalar).
    - verbose: If True, print intermediate results for every 10000th column.
    Returns:
    - ps_fv: Updated dry surface pressure array.
    - pw_fv: Precipitable water array.
    """
    # Determine input dimensions
    if ps_fv.ndim == 1:
        is_1d = True
        ps_fv = ps_fv.reshape(1, -1)
        q_fv = q_fv.reshape(q_fv.shape[0], 1, -1)
    elif ps_fv.ndim == 2:
        is_1d = False
        nlat, nlon = ps_fv.shape
        ps_fv = ps_fv.reshape(1, -1)
        q_fv = q_fv.reshape(q_fv.shape[0], 1, -1)
    else:
        raise ValueError("ps_fv must be 1D or 2D")

    # Call the Numba-optimized core function
    ps_fv, pw_fv = ps_wet_to_dry_conversion_core(ps_fv, q_fv, hyai, hybi, p0)

    if verbose:
        ncol = ps_fv.shape[1]
        for kk in range(0, ncol, 10000):
            logging.info(f"Correcting PS: {kk} of {ncol-1} from {ps_fv[0, kk] + pw_fv[0, kk] * grav} to {ps_fv[0, kk]} since TPW: {pw_fv[0, kk]}")

    logging.info("Done!")

    # Reshape outputs to original dimensions
    if is_1d:
        ps_fv = ps_fv.reshape(-1)
        pw_fv = pw_fv.reshape(-1)
    else:
        ps_fv = ps_fv.reshape(nlat, nlon)
        pw_fv = pw_fv.reshape(nlat, nlon)

    return ps_fv, pw_fv


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


def calculate_rho_dry(pres_gfs, q_gfs, t_gfs, rho_d_algo=1):
    """
    Calculates the dry air density (rho_dry) using either a simple or advanced method.

    Parameters:
    -----------
    pres_gfs : numpy.ndarray
        Pressure field
    q_gfs : numpy.ndarray
        Specific humidity field
    t_gfs : numpy.ndarray
        Temperature field

    Returns:
    --------
    rho_dry : numpy.ndarray
        Dry air density with the same shape as pres_gfs.
    """
    Rd = 287.05  # Specific gas constant for dry air in J/(kg·K)
    Rv = 461.5  # Specific gas constant for water vapor in J/(kg·K)
    Rv_over_Rd = Rv / Rd

    e = (q_gfs * pres_gfs) / (0.622 + 0.378 * q_gfs)
    pd = pres_gfs - e
    rho_dry = pd / (Rd * t_gfs)

    return rho_dry


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


def damp_upper_level_winds(u, v, nlev, damping_coeffs=None):
    """
    Damps the upper-level MPAS winds using specified damping coefficients.

    Parameters:
    -----------
    u : numpy.ndarray
        U-component of wind, typically with shape (nlev, ncol).
    v : numpy.ndarray
        V-component of wind, typically with shape (nlev, ncol).
    mpas_nlev : int
        Number of vertical levels.
    damping_coeffs : list or numpy.ndarray, optional
        Damping coefficients for the upper levels. Defaults to values in constants.py.

    Returns:
    --------
    u : numpy.ndarray
        U-component of wind after damping.
    v : numpy.ndarray
        V-component of wind after damping.
    """
    if damping_coeffs is None:
        damping_coeffs = mpas_uv_damping_coeffs

    logging.info(f"Damping upper level MPAS winds with coefficients: {damping_coeffs}")
    u[:, nlev - 1] *= damping_coeffs[0]
    u[:, nlev - 2] *= damping_coeffs[1]
    u[:, nlev - 3] *= damping_coeffs[2]
    v[:, nlev - 1] *= damping_coeffs[0]
    v[:, nlev - 2] *= damping_coeffs[1]
    v[:, nlev - 3] *= damping_coeffs[2]
    logging.info("... done damping upper level MPAS winds")

    return u, v


def noflux_boundary_condition(w, nlev):
    logging.info("Setting lower BC for W so flow can't go through surface...")
    w[:, 0] = 0.0
    logging.info("... done setting lower BC for W so flow can't go through surface")

    return w


@nb.jit(nopython=True)
def calculate_geopotential_height_unstructured(PS, T, Q, PRES):
    """
    Calculate geopotential height for unstructured grid using hydrostatic integration.

    Parameters:
    -----------
    PS : array (ncol,)
        Surface pressure [Pa]
    T : array (ncol, nlev)
        Temperature [K] at layer midpoints
    Q : array (ncol, nlev)
        Water vapor mixing ratio [kg/kg]
    PRES : array (ncol, nlev)
        Pressure [Pa] at layer midpoints

    Returns:
    --------
    Z : array (ncol, nlev)
        Geopotential height [m] at layer midpoints
    """
    ncol, nlev = T.shape
    Z = np.zeros_like(T)
    top_to_bottom = False

    for i in range(ncol):

        pint = np.zeros(nlev)  # Interface pressures at TOP of each layer

        if PRES[i, 0] < PRES[i, -1]:  # Top-to-bottom ordering
            for k in range(nlev):
                top_to_bottom = True
                layer_idx = nlev - 1 - k
                if k == nlev - 1:  # Top layer - extrapolate above
                    pint[k] = PRES[i, 0] * 0.5
                else:  # All other layers - midpoint between this and next layer up
                    pint[k] = (PRES[i, layer_idx] + PRES[i, layer_idx-1]) / 2.0
        else:  # Bottom-to-top ordering
            for k in range(nlev):
                if k == nlev - 1:  # Top layer - extrapolate above
                    pint[k] = PRES[i, -1] * 0.5
                else:  # All other layers - midpoint between this and next layer up
                    pint[k] = (PRES[i, k] + PRES[i, k+1]) / 2.0

        # Call hydrostatic integration - surface interface will be PS
        if top_to_bottom:
            Z[i, ::-1] = hydro_bottom_to_top(PS[i], pint, T[i, ::-1], Q[i, ::-1], zsfc=0.0)
        else:
            Z[i, :] = hydro_bottom_to_top(PS[i], pint, T[i, :], Q[i, :], zsfc=0.0)

    return Z
