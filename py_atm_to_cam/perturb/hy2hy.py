import numpy as np

def hyi2hyo(p0: float,
            hyai: np.ndarray,
            hybi: np.ndarray,
            ps: np.ndarray,
            xi: np.ndarray,
            hyao: np.ndarray,
            hybo: np.ndarray,
            intflg: int = 0) -> np.ndarray:
    """
    Interpolates data from one set of hybrid levels to another set of hybrid levels.

    Parameters
    ----------
    p0 : float
        Surface reference pressure. Must have same units as ps.
    hyai : ndarray
        1D array of hybrid coefficients A for input levels (top-to-bottom). Unitless.
    hybi : ndarray
        1D array of hybrid coefficients B for input levels (top-to-bottom). Unitless.
    ps : ndarray
        N-dimensional array of pressures. Must be one dimension smaller than xi.
    xi : ndarray
        N-dimensional array to interpolate. Rightmost dims must be (level,lat,lon).
    hyao : ndarray
        1D array of hybrid coefficients A for output levels (top-to-bottom). Unitless.
    hybo : ndarray
        1D array of hybrid coefficients B for output levels (top-to-bottom). Unitless.
    intflg : int, optional
        Extrapolation type. 0=missing values outside range, 1=use nearest level.

    Returns
    -------
    ndarray
        Interpolated data on new hybrid levels
    """
    # Input validation
    if ps.ndim < 2:
        raise ValueError("ps must have at least 2 dimensions")
    if xi.ndim < 3:
        raise ValueError("xi must have at least 3 dimensions")

    # Get dimensions
    klevi = hyai.shape[0]  # number of input levels
    klevo = hyao.shape[0]  # number of output levels

    # Get lat/lon dimensions from rightmost dims of xi
    nlat = xi.shape[-2]
    mlon = xi.shape[-1]

    # Verify dimensions match
    if ps.shape[-2:] != (nlat, mlon):
        raise ValueError("ps must have same rightmost dimensions as xi")
    if hyai.shape[0] != hybi.shape[0]:
        raise ValueError("hyai and hybi must have same length")
    if hyao.shape[0] != hybo.shape[0]:
        raise ValueError("hyao and hybo must have same length")
    if xi.shape[-3] != klevi:
        raise ValueError("xi level dimension must match length of hyai/hybi")

    # Calculate size of leftmost dimensions
    leftmost_shape = xi.shape[:-3]
    size_leftmost = int(np.prod(leftmost_shape))

    # Reshape arrays for processing
    ps_flat = ps.reshape(size_leftmost, nlat * mlon)
    xi_flat = xi.reshape(size_leftmost, klevi, nlat * mlon)

    # Initialize output array
    out_shape = leftmost_shape + (klevo, nlat, mlon)
    xo = np.zeros(out_shape, dtype=xi.dtype)
    xo_flat = xo.reshape(size_leftmost, klevo, nlat * mlon)

    # Get missing value (use default if none provided)
    missing_value = getattr(xi, 'missing_value', np.nan)

    # Loop over each point in leftmost dimensions
    for i in range(size_leftmost):
        # Calculate pressure levels for input grid
        pin = np.zeros((klevi, nlat * mlon))
        for k in range(klevi):
            pin[k] = hyai[k] * p0 + hybi[k] * ps_flat[i]

        # Calculate pressure levels for output grid
        pout = np.zeros((klevo, nlat * mlon))
        for k in range(klevo):
            pout[k] = hyao[k] * p0 + hybo[k] * ps_flat[i]

        # Do log-pressure interpolation for each lat/lon point
        for j in range(nlat * mlon):
            pin_col = pin[:, j]
            pout_col = pout[:, j]
            xi_col = xi_flat[i, :, j]

            # Convert to log pressure
            log_pin = np.log(pin_col)
            log_pout = np.log(pout_col)

            # Do the interpolation
            if intflg == 0:
                # Set values outside range to missing
                xo_flat[i, :, j] = np.interp(
                    log_pout,
                    log_pin,
                    xi_col,
                    left=missing_value,
                    right=missing_value
                )
            else:
                # Use nearest value outside range
                xo_flat[i, :, j] = np.interp(
                    log_pout,
                    log_pin,
                    xi_col,
                    left=xi_col[0],
                    right=xi_col[-1]
                )

    return xo