gamma_s = 5.0 / 1000.0  # Environmental lapse rate in K/m
gamma_d = 9.8 / 1000.0  # Adiabatic lapse rate in K/m
grav = 9.80616
Rd = 287.058
Rv = 461.6
Rv_over_Rd = Rv / Rd
p0 = 100000

vert_interp_thresh = 0.1  # ps corr diff (Pa) req. to interp vert profiles
extrap_threshold = 5000.  # maximum ps corr diff (Pa) to allow extrapolation

NC_FLOAT_FILL=9.96921e+36

dtime_map = [4, 2, 2, 2]