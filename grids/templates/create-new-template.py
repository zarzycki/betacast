#!/usr/bin/env python3
# Python port of create-new-template.ncl
#
# Doubles the vertical resolution of an existing template by inserting a new interface
# halfway between each pair of input interfaces. Midpoints are then recomputed as the
# average of adjacent interfaces. If keeplowestlevel, the lowest input layer is not split.
#
# python create-new-template.py

import os
import numpy as np
import netCDF4 as nc

inttype = "double"   # double, quadruple
keeplowestlevel = True

in_template = "./L30template.nc"
ifile = nc.Dataset(in_template, "r")
ifile.set_auto_mask(False)

ihyam = ifile.variables["hyam"][:]
ihybm = ifile.variables["hybm"][:]
ihyai = ifile.variables["hyai"][:]
ihybi = ifile.variables["hybi"][:]
ilev  = ifile.variables["lev"][:]
iilev = ifile.variables["ilev"][:]

inlev = len(ilev)
inint = len(iilev)
onlev = 2*inlev
onint = onlev+1

if keeplowestlevel:
    onlev = onlev-1
    onint = onint-1

print("(0)\tNumber of input levels = "+str(inlev))
print("(0)\tNumber of input interfaces = "+str(inint))
print("(0)\tNumber of output levels = "+str(onlev))
print("(0)\tNumber of output interfaces = "+str(onint))

# NCL new() fills with the default _FillValue of the type
ohyam = np.full(onlev, nc.default_fillvals[ihyam.dtype.str[1:]], dtype=ihyam.dtype)
ohybm = np.full(onlev, nc.default_fillvals[ihybm.dtype.str[1:]], dtype=ihybm.dtype)
ohyai = np.full(onint, nc.default_fillvals[ihyai.dtype.str[1:]], dtype=ihyai.dtype)
ohybi = np.full(onint, nc.default_fillvals[ihybi.dtype.str[1:]], dtype=ihybi.dtype)
olev  = np.full(onlev, nc.default_fillvals[ilev.dtype.str[1:]], dtype=ilev.dtype)
oilev = np.full(onint, nc.default_fillvals[iilev.dtype.str[1:]], dtype=iilev.dtype)

for ii in range(0, onint):
    if ii == 0:
        ohyai[0] = ihyai[0]
        ohybi[0] = ihybi[0]
        oilev[0] = iilev[0]
    elif ii == onint-1:
        ohyai[ii] = ihyai[inint-1]
        ohybi[ii] = ihybi[inint-1]
        oilev[ii] = iilev[inint-1]
    else:
        if ii % 2 == 0:
            ohyai[ii] = ihyai[ii//2]
            ohybi[ii] = ihybi[ii//2]
            oilev[ii] = iilev[ii//2]
        else:
            ohyai[ii] = ( ihyai[(ii+1)//2] + ihyai[(ii-1)//2] ) / 2.
            ohybi[ii] = ( ihybi[(ii+1)//2] + ihybi[(ii-1)//2] ) / 2.
            oilev[ii] = ( iilev[(ii+1)//2] + iilev[(ii-1)//2] ) / 2.
            print("(0)\t"+str(ii)+" "+"%.16g" % oilev[ii])

# now calculate midpoints, which are just literally the midpoints

for jj in range(0, onlev):
    ohyam[jj] = (ohyai[jj] + ohyai[jj+1]) / 2.
    ohybm[jj] = (ohybi[jj] + ohybi[jj+1]) / 2.
    olev[jj]  = (oilev[jj] + oilev[jj+1]) / 2.

# Metadata, as NCL carries it: olev/ohyam/ohybm only have the _FillValue from new().
# The ohyai/ohybi element assignments copy the input attributes (and the input ilev
# coordinate, whose attributes end up on the output ilev), plus _FillValue from new().
olev_atts  = {"_FillValue": nc.default_fillvals[olev.dtype.str[1:]]}
ohyam_atts = {"_FillValue": nc.default_fillvals[ohyam.dtype.str[1:]]}
ohybm_atts = {"_FillValue": nc.default_fillvals[ohybm.dtype.str[1:]]}
ohyai_atts = {k: ifile.variables["hyai"].getncattr(k) for k in ifile.variables["hyai"].ncattrs() if k != "_FillValue"}
ohyai_atts["_FillValue"] = nc.default_fillvals[ohyai.dtype.str[1:]]
ohybi_atts = {k: ifile.variables["hybi"].getncattr(k) for k in ifile.variables["hybi"].ncattrs() if k != "_FillValue"}
ohybi_atts["_FillValue"] = nc.default_fillvals[ohybi.dtype.str[1:]]
oilev_atts = {k: ifile.variables["ilev"].getncattr(k) for k in ifile.variables["ilev"].ncattrs() if k != "_FillValue"}
oilev_atts["_FillValue"] = nc.default_fillvals[oilev.dtype.str[1:]]

#===================================================================
# write output template
#===================================================================
out_template = "./_L"+str(onlev)+"template.nc"
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
ifile.close()
