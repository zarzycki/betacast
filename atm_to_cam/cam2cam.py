import numpy as np
import xarray as xr

def omega_ccm(u, v, div, dpsl, dpsm, pmid, pdel, psfc, hybdif, hybm, nprlev):
    """
    Computes the vertical pressure velocity (omega) using model diagnostic methods.

    Parameters:
    -----------
    u, v : ndarray
        Zonal and meridional wind components. The three rightmost dimensions must be (lev, lat, lon).

    div : ndarray
        Divergence with the same shape as u and v.

    dpsl, dpsm : ndarray
        Longitudinal and latitudinal components of grad ln(psfc). Shape is (lat, lon).

    pmid : ndarray
        Mid-level pressure values. Same shape as u and v.

    pdel : ndarray
        Layer pressure thickness values. Same shape as u and v.

    psfc : ndarray
        Surface pressure with dimensions (lat, lon).

    hybdif : ndarray
        The difference between the hybrid interface coefficients [eg, hybi(k+1)-hybi(k)].

    hybm : ndarray
        The hybrid B coefficients. Same size as the level dimension of u and v.

    nprlev : int
        Number of pure pressure levels.

    Returns:
    --------
    omega : ndarray
        Vertical pressure velocity with the same shape as u and v.
    """

    # Print the shapes of the input variables
    print(f"u shape: {u.shape}")
    print(f"v shape: {v.shape}")
    print(f"div shape: {div.shape}")
    print(f"dpsl shape: {dpsl.shape}")
    print(f"dpsm shape: {dpsm.shape}")
    print(f"pmid shape: {pmid.shape}")
    print(f"pdel shape: {pdel.shape}")
    print(f"psfc shape: {psfc.shape}")
    print(f"hybdif shape: {hybdif.shape}")
    print(f"hybm shape: {hybm.shape}")
    print(f"nprlev: {nprlev}")

    klev, jlat, ilon = u.shape[-3], u.shape[-2], u.shape[-1]

    print(f"ilon: {ilon}, jlat: {jlat}, klev: {klev}")

    omega = np.zeros_like(u)

    # Initialize partial sums
    suml = np.zeros((jlat, ilon))  # Correct shape here

    # Inverse of pmid (no need to check for division by zero)
    rpmid = 1.0 / pmid

    # Pure pressure part: top level
    hkk = 0.5 * rpmid[0, :, :]
    omega[0, :, :] = -hkk * div[0, :, :] * pdel[0, :, :]
    suml += div[0, :, :] * pdel[0, :, :]

    if 1 >= nprlev:
        vgpk = (u[0, :, :] * dpsl + v[0, :, :] * dpsm) * psfc
        tmp = vgpk * hybdif[0]
        omega[0, :, :] += hybm[0] * rpmid[0, :, :] * vgpk - hkk * tmp
        suml += tmp

    # Integrals to level above bottom
    for k in range(1, klev - 1):
        hkk = 0.5 * rpmid[k, :, :]
        hlk = rpmid[k, :, :]
        omega[k, :, :] = -hkk * div[k, :, :] * pdel[k, :, :] - hlk * suml
        suml += div[k, :, :] * pdel[k, :, :]

        if k >= nprlev:
            vgpk = (u[k, :, :] * dpsl + v[k, :, :] * dpsm) * psfc
            tmp = vgpk * hybdif[k]
            omega[k, :, :] += hybm[k] * rpmid[k, :, :] * vgpk - hkk * tmp
            suml += tmp

    # Pure pressure part: bottom level
    hkk = 0.5 * rpmid[klev - 1, :, :]
    hlk = rpmid[klev - 1, :, :]
    omega[klev - 1, :, :] = -hkk * div[klev - 1, :, :] * pdel[klev - 1, :, :] - hlk * suml

    if klev >= nprlev:
        vgpk = (u[klev - 1, :, :] * dpsl + v[klev - 1, :, :] * dpsm) * psfc
        omega[klev - 1, :, :] += hybm[klev - 1] * rpmid[klev - 1, :, :] * vgpk - hkk * vgpk * hybdif[klev - 1]

    # Rescale omega/p to omega
    omega *= pmid

    return omega


def calculate_gradients(psfc, lat, lon):
    """
    Calculate the longitudinal and latitudinal components of the gradient of ln(psfc).

    Parameters:
    -----------
    psfc : ndarray
        2D array of surface pressure values (Pa) with dimensions (lat, lon).
    lat : ndarray
        1D array of latitude values (degrees).
    lon : ndarray
        1D array of longitude values (degrees).

    Returns:
    --------
    dpsl : ndarray
        Longitudinal component of grad ln(psfc) with the same shape as psfc.
    dpsm : ndarray
        Latitudinal component of grad ln(psfc) with the same shape as psfc.
    """

    # Convert latitude and longitude to radians for calculations
    lat_rad = np.deg2rad(lat)
    lon_rad = np.deg2rad(lon)

    # Calculate the logarithm of surface pressure
    ln_psfc = np.log(psfc)

    # Earth's radius in meters
    R = 6371000.0

    # Calculate the gradient of ln(psfc)
    dln_psfc_dy, dln_psfc_dx = np.gradient(ln_psfc)

    # Manually calculate one-sided differences at the boundaries for dln_psfc_dx
    #dln_psfc_dx[:, 0] = (ln_psfc[:, 1] - ln_psfc[:, 0]) / (lon_rad[1] - lon_rad[0])
    #dln_psfc_dx[:, -1] = (ln_psfc[:, -1] - ln_psfc[:, -2]) / (lon_rad[-1] - lon_rad[-2])
    #dln_psfc_dy[0, :] = (ln_psfc[1, :] - ln_psfc[0, :]) / (lat_rad[1] - lat_rad[0])
    #dln_psfc_dy[-1, :] = (ln_psfc[-1, :] - ln_psfc[-2, :]) / (lat_rad[-1] - lat_rad[-2])

    # Scale gradients by physical distances
    dpsl = dln_psfc_dx / (R * np.cos(lat_rad[:, np.newaxis]) * np.gradient(lon_rad))
    dpsm = dln_psfc_dy / (R * np.gradient(lat_rad)[:, np.newaxis])

    # Set the top and bottom rows of dpsl to zero
    dpsl[0, :] = 0.0
    dpsl[-1, :] = 0.0

    # Check if lat is descending
    if lat[0] > lat[-1]:
        print("Flipping dpsm because lats are descending")
        dpsm = -dpsm

    return dpsl, dpsm


def calculate_div_vort(lat, lon, u, v):
    """
    Calculate the divergence and vorticity of the wind field using centered finite differences.

    Parameters:
    -----------
    lat : ndarray
        1D array of latitude values (degrees).
    lon : ndarray
        1D array of longitude values (degrees).
    u : ndarray
        2D or 3D array of zonal wind component with dimensions (..., lat, lon).
    v : ndarray
        2D or 3D array of meridional wind component with dimensions (..., lat, lon).

    Returns:
    --------
    div : ndarray
        Divergence of the wind field with the same shape as u and v.
    vort : ndarray
        Vorticity of the wind field with the same shape as u and v.
    """

    # Convert latitude and longitude to radians for calculations
    lat_rad = np.deg2rad(lat)
    lon_rad = np.deg2rad(lon)

    # Earth's radius in meters
    R = 6371000.0

    # Calculate gradients along the last two dimensions (lat, lon)
    du_dy = np.gradient(u, axis=-2)
    du_dx = np.gradient(u, axis=-1)
    dv_dy = np.gradient(v, axis=-2)
    dv_dx = np.gradient(v, axis=-1)

    # Scale gradients by physical distances
    du_dx = du_dx / (R * np.cos(lat_rad)[:, np.newaxis] * np.gradient(lon_rad))
    dv_dy = dv_dy / (R * np.gradient(lat_rad)[:, np.newaxis])

    # Divergence: du/dx + dv/dy
    div = du_dx + dv_dy

    # Vorticity: dv/dx - du/dy
    du_dy_scaled = du_dy / (R * np.gradient(lat_rad)[:, np.newaxis])
    dv_dx_scaled = dv_dx / (R * np.cos(lat_rad)[:, np.newaxis] * np.gradient(lon_rad))
    vort = dv_dx_scaled - du_dy_scaled

    return div, vort







def omega_ccm_driver(p0, psfc, u, v, lat, lon, hyam, hybm, hyai, hybi):
    """
    Driver function to calculate intermediate quantities and then compute omega (vertical pressure velocity).

    Parameters:
    -----------
    p0 : float
        Surface reference pressure in Pa.
    psfc : numpy.ndarray
        2D or 3D array of surface pressures in Pa.
    u, v : numpy.ndarray
        3D or 4D arrays of zonal and meridional wind (m/s).
    hyam, hybm : numpy.ndarray
        1D arrays of hybrid A and B coefficients (mid-level). Must have the same dimension as the level dimension of u and v.
    hyai, hybi : numpy.ndarray
        1D arrays of interface hybrid A and B coefficients. Must have the same dimension as the interface levels.

    Returns:
    --------
    omega : numpy.ndarray
        3D or 4D array of vertical pressure velocity (Pa/s).
    """

    # Check the dimensions of inputs
    rank_ps = psfc.ndim
    rank_u = u.ndim

    if not (rank_u in [3, 4] and rank_ps in [2, 3]):
        raise ValueError(f"omega_ccm_driver: expected rank_u=3 or 4 and rank_ps=2 or 3, but got rank_u={rank_u} and rank_ps={rank_ps}")

    # Get dimensions
    if rank_u == 3:
        klev, nlat, mlon = u.shape
    else:
        ntim, klev, nlat, mlon = u.shape

    # Check if lat is descending
    if lat[0] > lat[-1]:
        print("Flipping lats in omega_ccm")
        flip_lats=True
    else:
        flip_lats=False

    # Initialize omega array
    omega = np.zeros_like(u)

    # Calculate hybrid differences
    hybd = np.diff(hybi)

    # Determine the number of pure pressure levels
    nprlev = np.argmax(hybi != 0) if np.any(hybi != 0) else 0

    # Calculate layer pressure thicknesses (pdel) and mid-level pressures (pmid)
    if flip_lats:
        pdel = dpres_hybrid_ccm(psfc[::-1,:], p0, hyai, hybi)
        pmid = pres_hybrid_ccm(psfc[::-1,:], p0, hyam, hybm)
    else:
        pdel = dpres_hybrid_ccm(psfc, p0, hyai, hybi)
        pmid = pres_hybrid_ccm(psfc, p0, hyam, hybm)

    # Calculate gradients of log(psfc)
    dpsl, dpsm = calculate_gradients(psfc, lat, lon)

    # Calculate divergence on Gaussian grid
    div, vort = calculate_div_vort(lat, lon, u, v)

    # Call omega_ccm to calculate omega
    omega = omega_ccm(u, v, div, dpsl, dpsm, pmid, pdel, psfc, hybd, hybm, nprlev)

    return omega

def pres_hybrid_ccm(psfc, p0, hya, hyb, pmsg=np.nan):
    """
    Calculate pressure at hybrid levels.

    Parameters:
    -----------
    psfc : numpy.ndarray
        2D array of surface pressures in Pa (shape: [lat, lon]).
    p0 : float
        Base pressure in Pa.
    hya : numpy.ndarray
        1D array of "a" or pressure hybrid coefficients (shape: [klev]).
    hyb : numpy.ndarray
        1D array of "b" or sigma coefficients (shape: [klev]).
    pmsg : float, optional
        Missing value indicator, defaults to NaN.

    Returns:
    --------
    phy : numpy.ndarray
        3D array of pressure at hybrid levels (shape: [klev, lat, lon]).
    """
    klev = len(hya)
    nlat, mlon = psfc.shape
    phy = np.full((klev, nlat, mlon), pmsg)

    for kl in range(klev):
        for nl in range(nlat):
            for ml in range(mlon):
                if not np.isnan(psfc[nl, ml]):
                    phy[kl, nl, ml] = hya[kl] * p0 + hyb[kl] * psfc[nl, ml]

    return phy

def dpres_hybrid_ccm(psfc, p0, hyai, hybi, pmsg=np.nan):
    """
    Calculate delta pressure between hybrid levels.

    Parameters:
    -----------
    psfc : numpy.ndarray
        2D array of surface pressures in Pa (shape: [lat, lon]).
    p0 : float
        Base pressure in Pa.
    hyai : numpy.ndarray
        1D array of interface "a" or pressure hybrid coefficients (shape: [klev+1]).
    hybi : numpy.ndarray
        1D array of interface "b" or sigma coefficients (shape: [klev+1]).
    pmsg : float, optional
        Missing value indicator, defaults to NaN.

    Returns:
    --------
    dphy : numpy.ndarray
        3D array of delta pressure between hybrid levels (shape: [klev-1, lat, lon]).
    """
    klev = len(hyai) - 1
    nlat, mlon = psfc.shape
    dphy = np.full((klev, nlat, mlon), pmsg)

    for nl in range(nlat):
        for ml in range(mlon):
            if not np.isnan(psfc[nl, ml]):
                for kl in range(klev):
                    pa = p0 * hyai[kl] + hybi[kl] * psfc[nl, ml]
                    pb = p0 * hyai[kl + 1] + hybi[kl + 1] * psfc[nl, ml]
                    dphy[kl, nl, ml] = abs(pb - pa)

    return dphy











def cz2ccm(ps, phis, tv, p0, hyam, hybm, hyai, hybi, debug=False):
    """
    Calculate geopotential height using the hybrid coordinate system.

    Parameters:
    -----------
    ps : ndarray
        2D array of surface pressures (Pa) with dimensions (lat, lon).
    phis : ndarray
        2D array of surface geopotential with dimensions (lat, lon).
    tv : ndarray
        3D array of virtual temperature (K) with dimensions (lev, lat, lon).
    p0 : float
        Base pressure (Pa).
    hyam : ndarray
        1D array of hybrid A coefficients for mid-levels.
    hybm : ndarray
        1D array of hybrid B coefficients for mid-levels.
    hyai : ndarray
        1D array of hybrid A coefficients for interfaces.
    hybi : ndarray
        1D array of hybrid B coefficients for interfaces.

    Returns:
    --------
    z2 : ndarray
        3D array of geopotential height (m) with dimensions (lev, lat, lon).
    """

    if debug:
        print("cz2ccm function called with the following inputs:")
        print(f"ps.shape: {ps.shape}")
        print(f"phis.shape: {phis.shape}")
        print(f"tv.shape: {tv.shape}")
        print(f"p0: {p0}")
        print(f"hyam.shape: {hyam.shape}")
        print(f"hybm.shape: {hybm.shape}")
        print(f"hyai.shape: {hyai.shape}")
        print(f"hybi.shape: {hybi.shape}")
        print(hyam)
        print(hybm)
        print(hyai)
        print(hybi)

    nlat, mlon = ps.shape
    klev = len(hyam)
    klev1 = len(hyai)

    z2 = np.zeros((klev, nlat, mlon))
    pmln = np.zeros((klev + 1, nlat, mlon))
    pterm = np.zeros((klev, nlat, mlon))

    if debug:
        print(f"nlat: {nlat}, mlon: {mlon}, klev: {klev}, klev1: {klev1}")
        print(f"z2.shape: {z2.shape}")
        print(f"pmln.shape: {pmln.shape}")
        print(f"pterm.shape: {pterm.shape}")

    hyba = np.zeros((2, klev + 1))
    hybb = np.zeros((2, klev + 1))

    # Copy to temporary arrays
    hyba[0, :] = hyai
    hybb[0, :] = hybi

    # Copy midpoint coefficients to the second row (index 1) of HYBA and HYBB
    hyba[1, 1:klev1] = hyam
    hybb[1, 1:klev1] = hybm

    if debug:
        print("Intermediate arrays initialized.")
        print(f"hyba.shape: {hyba.shape}")
        print(f"hybb.shape: {hybb.shape}")

    pmln[0, :, :] = np.log(p0 * hyba[1, klev] + ps * hybb[0, klev])
    pmln[-1, :, :] = np.log(p0 * hyba[1, 0] + ps * hybb[0, 0])
    for k in range(klev-1, 0, -1):
        idx = klev - k
        arg = p0 * hyba[1, idx] + ps * hybb[1, idx]
        pmln[k, :, :] = np.where(arg > 0.0, np.log(arg), 0.0)

    # Calculate geopotential height Z2
    R = 287.04  # Gas constant for dry air (J/(kg*K))
    G0 = 9.80616  # Gravity (m/s^2)
    RBYG = R / G0

    # Eq 3.a.109.2
    for k in range(1, klev-1):
        pterm[k, :, :] = RBYG * tv[k, :, :] * 0.5 * (pmln[k+1, :, :] - pmln[k-1, :, :])

    # Eq 3.a.109.5 and 3.a.109.2
    for k in range(0, klev-1):
        z2[k, :, :] = phis / G0 + RBYG * tv[k, :, :] * 0.5 * (pmln[k+1, :, :] - pmln[k, :, :])

    # Step 5: Special Case for Last Layer (3.a.109.5)
    k = klev-1
    z2[k, :, :] = phis / G0 + RBYG * tv[k, :, :] * (np.log(ps * hybb[0, 0]) - pmln[k, :, :])

    # Eq 3.a.109.4
    for k in range(0, klev-1):
        l = klev-1
        z2[k, :, :] = z2[k, :, :] + RBYG * tv[l, :, :] * (np.log(ps * hybb[0, 0]) - 0.5 * (pmln[l-1, :, :] + pmln[l, :, :]))

    # Add thickness of the remaining full layers (Eq 3.a.109.3)
    for k in range(0,klev-2):
        for l in range(k+1, klev):
            z2[k, :, :] = z2[k, :, :] + pterm[l, :, :]

    return z2





def ddvfidf_wrapper(u, v, glat, glon, iopt, xmsg=np.nan):
    """
    Wrapper function to handle 2D and 3D input arrays for calculating the divergence of the wind field.

    Parameters:
    -----------
    u : ndarray
        2D or 3D array of zonal wind component with dimensions (nlat, mlon) or (nlev, nlat, mlon).
    v : ndarray
        2D or 3D array of meridional wind component with dimensions (nlat, mlon) or (nlev, nlat, mlon).
    glat : ndarray
        1D array of latitude values (degrees) with dimensions (nlat).
    glon : ndarray
        1D array of longitude values (degrees) with dimensions (mlon).
    xmsg : float
        Missing value indicator.
    iopt : int
        Option for handling cyclic grids.
        - 1 or 3: Cyclic in longitude.
        - 2: Non-cyclic in longitude.

    Returns:
    --------
    div : ndarray
        2D or 3D array of divergence with dimensions (nlat, mlon) or (nlev, nlat, mlon).
    ier : int
        Error code, 0 if successful.
    """

    # Check if the input is 2D or 3D
    if u.ndim == 2:
        # 2D case: Directly call the ddvfidf function
        print("2d in ddvfidf_wrapper")
        div = ddvfidf(u, v, glat, glon, iopt, xmsg)

    elif u.ndim == 3:
        # 3D case: Initialize the output divergence array
        print("3d in ddvfidf_wrapper")
        nlev = u.shape[0]
        nlat = u.shape[1]
        mlon = u.shape[2]
        div = np.full((nlev, nlat, mlon), xmsg)

        # Loop over the levels and compute divergence for each level
        for lev in range(nlev):
            div[lev, :, :] = ddvfidf(u[lev, :, :], v[lev, :, :], glat, glon, iopt, xmsg)

    else:
        raise ValueError("Input arrays u and v must be either 2D or 3D.")

    return div



def ddvfidf(u, v, glat, glon, iopt, xmsg=np.nan):
    """
    Calculate the divergence of the wind field using centered finite differences.

    Parameters:
    -----------
    u : ndarray
        2D array of zonal wind component with dimensions (mlon, nlat).
    v : ndarray
        2D array of meridional wind component with dimensions (mlon, nlat).
    glat : ndarray
        1D array of latitude values (degrees) with dimensions (nlat).
    glon : ndarray
        1D array of longitude values (degrees) with dimensions (mlon).
    xmsg : float
        Missing value indicator.
    iopt : int
        Option for handling cyclic grids.
        - 1 or 3: Cyclic in longitude.
        - 2: Non-cyclic in longitude.

    Returns:
    --------
    dv : ndarray
        2D array of divergence with the same dimensions as u and v.
    ier : int
        Error code, 0 if successful.
    """

    # Transpose u and v to ensure they have shape (mlon, nlat)
    transposed = False
    if u.shape != (len(glon), len(glat)):
        u = u.T
        v = v.T
        transposed = True

    # Check if latitude is descending and flip if necessary
    flipped = False
    if glat[0] > glat[-1]:
        glat = glat[::-1]
        u = u[:, ::-1]
        v = v[:, ::-1]
        flipped = True

    # Constants
    RE = 6.37122e6  # Earth's radius in meters
    RAD = np.pi / 180.0
    RCON = RE * RAD
    nlat = len(glat)
    mlon = len(glon)

    print(u.shape)
    print(v.shape)
    print(len(glat))
    print(len(glon))

    # Initialize output array and error code
    dv = np.full((mlon, nlat), xmsg)

    # Pre-compute cos(lat) and tan(lat)/RE
    clat = np.cos(RAD * glat)
#     tlatre = np.zeros(nlat)
#     for nl in range(nlat):
#         if abs(glat[nl]) < 90.0:
#             tlatre[nl] = np.tan(RAD * glat[nl]) / RE
#         else:
#             if glat[nl] == 90.0:
#                 polat = 0.5 * (glat[nl] + glat[nl - 1])
#             elif glat[nl] == -90.0:
#                 polat = 0.5 * (glat[nl] + glat[nl - 1])
#             else:
#                 polat = glat[nl]  # Or any other appropriate value for polar regions
#
#             tlatre[nl] = np.tan(RAD * polat) / RE

    tlatre = np.zeros(nlat)
    for nl in range(nlat):
        if abs(glat[nl]) < 89.9999:
            tlatre[nl] = np.tan(RAD * glat[nl]) / RE
        else:
            if glat[nl] >= 89.9999:
                print(f"NP {nl}")
                polat = 0.5 * (glat[nl] + glat[nl - 1])
                print(polat)
                tlatre[nl] = np.tan(RAD * polat) / RE
            else:  # This corresponds to the case where glat[nl] == -90.0
                print(f"SP {nl}")
                polat = 0.5 * (glat[nl] + glat[nl + 1])
                print(polat)
                tlatre[nl] = np.tan(RAD * polat) / RE

    print(tlatre)

    # Calculate 1/dy and 1/(2*dy)
    dybot = 1.0 / (RCON * (glat[1] - glat[0]))
    dytop = 1.0 / (RCON * (glat[nlat-1] - glat[nlat-2]))
    dy2 = np.zeros(nlat)
    for nl in range(1, nlat - 1):
        dy2[nl] = 1.0 / (RCON * (glat[nl + 1] - glat[nl - 1]))

    print(dybot)
    print(dytop)
    print(dy2)

    # Calculate 1/dx and 1/(2*dx)
    dlon = glon[1] - glon[0]
    dlon2 = glon[2] - glon[0]
    dx = np.zeros(nlat)
    dx2 = np.zeros(nlat)
    for nl in range(nlat):
        if abs(glat[nl]) <= 89.99999:
            dx[nl] = 1.0 / (RCON * dlon * clat[nl])
            dx2[nl] = 1.0 / (RCON * dlon2 * clat[nl])
        else:
            dx[nl] = 0.0
            dx2[nl] = 0.0

    # Set subscript range
    jopt = abs(iopt)
    if jopt == 1 or jopt == 3:
        mlstrt = 0
        mlend = mlon
    else:
        mlstrt = 1
        mlend = mlon - 1

    # Calculate divergence in the grid body
    for ml in range(mlstrt, mlend):
        mlm1 = ml - 1
        mlp1 = ml + 1
        if ml == 0:
            mlm1 = mlon - 1
        if ml == mlon - 1:
            mlp1 = 0

        for nl in range(1, nlat - 1):
            dv[ml, nl] = (v[ml, nl + 1] - v[ml, nl - 1]) * dy2[nl] + \
                         (u[mlp1, nl] - u[mlm1, nl]) * dx2[nl] - \
                         v[ml, nl] * tlatre[nl]

        # Bottom and top boundaries (nl = 1, nlat)
        if jopt >= 2:
            # Bottom boundary (nl = 1)
            dv[ml, 0] = (v[ml, 1] - v[ml, 0]) * dybot + \
                        (u[mlp1, 0] - u[mlm1, 0]) * dx2[0] - \
                        v[ml, 0] * tlatre[0]

            # Top boundary (nl = nlat)
            dv[ml, nlat - 1] = (v[ml, nlat - 1] - v[ml, nlat - 2]) * dytop + \
                               (u[mlp1, nlat - 1] - u[mlm1, nlat - 1]) * dx2[nlat - 1] - \
                               v[ml, nlat - 1] * tlatre[nlat - 1]

    # Left and right boundaries (ml = 1, mlon) for JOPT = 2
    if jopt == 2:
        for nl in range(1, nlat - 1):
            # Left boundary (ml = 1)
            if v[0, nl + 1] != xmsg and v[0, nl - 1] != xmsg and \
               u[1, nl] != xmsg and u[0, nl] != xmsg and \
               v[0, nl] != xmsg:
                dv[0, nl] = (v[0, nl + 1] - v[0, nl - 1]) * dy2[nl] + \
                            (u[1, nl] - u[0, nl]) * dx[nl] - \
                            v[0, nl] * tlatre[nl]

            # Right boundary (ml = mlon)
            if v[mlon - 1, nl + 1] != xmsg and v[mlon - 1, nl - 1] != xmsg and \
               u[mlon - 1, nl] != xmsg and u[mlon - 2, nl] != xmsg and \
               v[mlon - 1, nl] != xmsg:
                dv[mlon - 1, nl] = (v[mlon - 1, nl + 1] - v[mlon - 1, nl - 1]) * dy2[nl] + \
                                   (u[mlon - 1, nl] - u[mlon - 2, nl]) * dx[nl] - \
                                   v[mlon - 1, nl] * tlatre[nl]

    count_nan_or_zero = np.count_nonzero(np.isnan(dv) | (dv == 0))
    print("Number of NaN values or values equal to 0 in dv:", count_nan_or_zero)

    # Special handling for poles (lat = ±90)
    if abs(glat[0]) >= 89.9999:
        dv[:, 0] = np.mean(dv[:, 0])
    if abs(glat[-1]) >= 89.9999:
        dv[:, -1] = np.mean(dv[:, -1])

    count_nan_or_zero = np.count_nonzero(np.isnan(dv) | (dv == 0))
    print("Number of NaN values or values equal to 0 in dv:", count_nan_or_zero)
    #dv[np.isnan(dv)] = 0

    if flipped:
        dv = dv[:,::-1]
    if transposed:
        dv = dv.T

    return dv

