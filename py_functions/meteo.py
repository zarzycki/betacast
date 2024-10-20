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
