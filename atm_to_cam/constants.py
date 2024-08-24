gamma_s = 5.0 / 1000.0  # Environmental lapse rate in K/m
gamma_d = 9.8 / 1000.0  # Adiabatic lapse rate in K/m
grav = 9.80616
Rd = 287.058
Rv = 461.6
cp = 1004.
Rv_over_Rd = Rv / Rd
p0 = 100000
kappa = Rd / cp

vert_interp_thresh = 0.1  # ps corr diff (Pa) req. to interp vert profiles
extrap_threshold = 5000.  # maximum ps corr diff (Pa) to allow extrapolation

NC_FLOAT_FILL=9.96921e+36

dtime_map = [4, 2, 2, 2]

QMAXTHRESH=0.05
QMINTHRESH=1.0e-12
CLDMINTHRESH=0.0
OMEGA_LAT_THRESH = 88.0

rho_d_algo = 1

damp_upper_winds_mpas = False
mpas_uv_damping_coeffs = [0.90, 0.95, 0.98]

ps_wet_to_dry=False
output_diag=True

w_smooth_iter=0

t_freeze_K=273.15