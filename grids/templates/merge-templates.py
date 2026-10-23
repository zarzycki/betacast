#!/usr/bin/env python3
# Python port of merge-templates.ncl
#
# Merges two templates: interfaces from top_template above splitlev and interfaces from
# bot_template at or below splitlev. Midpoints are then recomputed as the average of
# adjacent interfaces.
#
# python merge-templates.py

import os
import numpy as np
import netCDF4 as nc


def ncl_print(name, x):
    # Mimic NCL print() of a 1-D numeric variable with no attributes/coordinates
    ncltype = {"i4": "integer", "f4": "float", "f8": "double"}[x.dtype.str[1:]]
    fmt = {"i4": "%d", "f4": "%2.7g", "f8": "%4.16g"}[x.dtype.str[1:]]
    print("\n\nVariable: "+name)
    print("Type: "+ncltype)
    print("Total Size: "+str(x.nbytes)+" bytes")
    print("            "+str(x.size)+" values")
    print("Number of Dimensions: 1")
    print("Dimensions and sizes:\t["+str(x.size)+"]")
    print("Coordinates: ")
    for i in range(0, x.size):
        print("("+str(i)+")\t"+fmt % x[i])


keeplowestlevel = False

top_template = "./L30template.nc"
bot_template = "./L59template.nc"
splitlev = 600.

tfile = nc.Dataset(top_template, "r")
tfile.set_auto_mask(False)
thyam = tfile.variables["hyam"][:]
thybm = tfile.variables["hybm"][:]
thyai = tfile.variables["hyai"][:]
thybi = tfile.variables["hybi"][:]
tlev  = tfile.variables["lev"][:]
tilev = tfile.variables["ilev"][:]

bfile = nc.Dataset(bot_template, "r")
bfile.set_auto_mask(False)
bhyam = bfile.variables["hyam"][:]
bhybm = bfile.variables["hybm"][:]
bhyai = bfile.variables["hyai"][:]
bhybi = bfile.variables["hybi"][:]
blev  = bfile.variables["lev"][:]
bilev = bfile.variables["ilev"][:]

tnlev = len(tlev)
tnint = len(tilev)
bnlev = len(blev)
bnint = len(bilev)

onlev = max([bnlev, tnlev])
onint = max([bnint, tnint])

print("(0)\tNumber of input levels = "+str(tnlev))
print("(0)\tNumber of input interfaces = "+str(tnint))
print("(0)\tNumber of input levels = "+str(bnlev))
print("(0)\tNumber of input interfaces = "+str(bnint))
print("(0)\tNumber of output levels = "+str(onlev))
print("(0)\tNumber of output interfaces = "+str(onint))

tilev = np.where(tilev < np.float32(splitlev), tilev, -1).astype(tilev.dtype)
bilev = np.where(bilev >= np.float32(splitlev), bilev, -1).astype(bilev.dtype)

indtopi = np.where(tilev > 0)[0].astype(np.int32)
indboti = np.where(bilev > 0)[0].astype(np.int32)

ncl_print("indtopi", indtopi)
ncl_print("indboti", indboti)

numtopi = len(indtopi)
numboti = len(indboti)

onint = numtopi + numboti
onlev = onint - 1

# NCL new() fills with the default _FillValue of the type
ohyam = np.full(onlev, nc.default_fillvals[thyam.dtype.str[1:]], dtype=thyam.dtype)
ohybm = np.full(onlev, nc.default_fillvals[thybm.dtype.str[1:]], dtype=thybm.dtype)
ohyai = np.full(onint, nc.default_fillvals[thyai.dtype.str[1:]], dtype=thyai.dtype)
ohybi = np.full(onint, nc.default_fillvals[thybi.dtype.str[1:]], dtype=thybi.dtype)
olev  = np.full(onlev, nc.default_fillvals[tlev.dtype.str[1:]], dtype=tlev.dtype)
oilev = np.full(onint, nc.default_fillvals[tilev.dtype.str[1:]], dtype=tilev.dtype)

oilev[0:numtopi] = tilev[indtopi]
oilev[numtopi:onint] = bilev[indboti]

ohyai[0:numtopi] = thyai[indtopi]
ohyai[numtopi:onint] = bhyai[indboti]

ohybi[0:numtopi] = thybi[indtopi]
ohybi[numtopi:onint] = bhybi[indboti]

# now calculate midpoints, which are just literally the midpoints

for jj in range(0, onlev):
    ohyam[jj] = (ohyai[jj] + ohyai[jj+1]) / 2.
    ohybm[jj] = (ohybi[jj] + ohybi[jj+1]) / 2.
    olev[jj]  = (oilev[jj] + oilev[jj+1]) / 2.

# Metadata, as NCL carries it: olev/ohyam/ohybm only have the _FillValue from new().
# The ohyai/ohybi subsection assignments copy the top attributes, then the bottom
# attributes overwrite/extend them (likewise for the ilev coordinate, whose attributes
# end up on the output ilev), plus _FillValue from new().
olev_atts  = {"_FillValue": nc.default_fillvals[olev.dtype.str[1:]]}
ohyam_atts = {"_FillValue": nc.default_fillvals[ohyam.dtype.str[1:]]}
ohybm_atts = {"_FillValue": nc.default_fillvals[ohybm.dtype.str[1:]]}
ohyai_atts = {k: tfile.variables["hyai"].getncattr(k) for k in tfile.variables["hyai"].ncattrs() if k != "_FillValue"}
ohyai_atts.update({k: bfile.variables["hyai"].getncattr(k) for k in bfile.variables["hyai"].ncattrs() if k != "_FillValue"})
ohyai_atts["_FillValue"] = nc.default_fillvals[ohyai.dtype.str[1:]]
ohybi_atts = {k: tfile.variables["hybi"].getncattr(k) for k in tfile.variables["hybi"].ncattrs() if k != "_FillValue"}
ohybi_atts.update({k: bfile.variables["hybi"].getncattr(k) for k in bfile.variables["hybi"].ncattrs() if k != "_FillValue"})
ohybi_atts["_FillValue"] = nc.default_fillvals[ohybi.dtype.str[1:]]
oilev_atts = {k: tfile.variables["ilev"].getncattr(k) for k in tfile.variables["ilev"].ncattrs() if k != "_FillValue"}
oilev_atts.update({k: bfile.variables["ilev"].getncattr(k) for k in bfile.variables["ilev"].ncattrs() if k != "_FillValue"})
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
tfile.close()
bfile.close()
