#!/usr/bin/env python3
# Python port of cut-templates.ncl
#
# Removes the interfaces listed in delLevs from an existing template. Midpoints are
# then recomputed as the average of adjacent interfaces.
#
# python cut-templates.py

import os
import numpy as np
import netCDF4 as nc

top_template = "./L56template.nc"

tfile = nc.Dataset(top_template, "r")
tfile.set_auto_mask(False)
thyam = tfile.variables["hyam"][:]
thybm = tfile.variables["hybm"][:]
thyai = tfile.variables["hyai"][:]
thybi = tfile.variables["hybi"][:]
tlev  = tfile.variables["lev"][:]
tilev = tfile.variables["ilev"][:]

tnlev = len(tlev)
tnint = len(tilev)

keeplev  = np.zeros(tnlev, dtype=bool)   # unused (as in NCL)
keepilev = np.zeros(tnint, dtype=bool)

keepilev[:] = True

#delLevs = [tnint-2, tnint-3, tnint-5, tnint-6, tnint-8]
delLevs = [tnint-2]

for ii in range(0, tnint):
    if any(np.array(delLevs) == ii):
        keepilev[ii] = False

ohyai = thyai[np.where(keepilev)[0]]
ohybi = thybi[np.where(keepilev)[0]]
oilev = tilev[np.where(keepilev)[0]]

onint = len(ohyai)
onlev = onint-1

# NCL new() fills with the default _FillValue of the type
ohyam = np.full(onlev, nc.default_fillvals[thyam.dtype.str[1:]], dtype=thyam.dtype)
ohybm = np.full(onlev, nc.default_fillvals[thybm.dtype.str[1:]], dtype=thybm.dtype)
olev  = np.full(onlev, nc.default_fillvals[tlev.dtype.str[1:]], dtype=tlev.dtype)

# now calculate midpoints, which are just literally the midpoints

for jj in range(0, onlev):
    ohyam[jj] = (ohyai[jj] + ohyai[jj+1]) / 2.
    ohybm[jj] = (ohybi[jj] + ohybi[jj+1]) / 2.
    olev[jj]  = (oilev[jj] + oilev[jj+1]) / 2.

# Metadata, as NCL carries it: olev/ohyam/ohybm only have the _FillValue from new().
# Subscripting with ind() copies the input attributes (and the ilev coordinate, whose
# attributes end up on the output ilev) in reverse order, with no _FillValue added.
olev_atts  = {"_FillValue": nc.default_fillvals[olev.dtype.str[1:]]}
ohyam_atts = {"_FillValue": nc.default_fillvals[ohyam.dtype.str[1:]]}
ohybm_atts = {"_FillValue": nc.default_fillvals[ohybm.dtype.str[1:]]}
ohyai_atts = {k: tfile.variables["hyai"].getncattr(k) for k in reversed(tfile.variables["hyai"].ncattrs())}
ohybi_atts = {k: tfile.variables["hybi"].getncattr(k) for k in reversed(tfile.variables["hybi"].ncattrs())}
oilev_atts = {k: tfile.variables["ilev"].getncattr(k) for k in reversed(tfile.variables["ilev"].ncattrs())}

#===================================================================
# write output template
#===================================================================
out_template = "./Z_L"+str(onlev)+"template.nc"
if os.path.exists(out_template):   # remove any pre-existing file
    os.remove(out_template)
ofile = nc.Dataset(out_template, "w", format="NETCDF3_CLASSIC")  # open output netCDF file

#===================================================================
# output variables; order matches NCL, which writes the lev (ilev)
# coordinate variable immediately after the first variable using it
#===================================================================

ofile.createDimension("lev", onlev)
ofile.createDimension("ilev", onint)
for name, dim, data, atts in [("hyam", "lev",  ohyam, ohyam_atts),
                              ("lev",  "lev",  olev,  olev_atts),
                              ("hybm", "lev",  ohybm, ohybm_atts),
                              ("hyai", "ilev", ohyai, ohyai_atts),
                              ("ilev", "ilev", oilev.astype(np.float32).astype(np.float64), oilev_atts),   # todouble(tofloat(oilev))
                              ("hybi", "ilev", ohybi, ohybi_atts)]:
    var = ofile.createVariable(name, data.dtype, (dim,), fill_value=False)
    var.setncatts(atts)   # setncatts allows _FillValue in NCL's attribute position
    var[:] = data

ofile.close()
tfile.close()
