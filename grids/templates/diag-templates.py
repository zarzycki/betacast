#!/usr/bin/env python3
# Python port of diag-templates.ncl
#
# Prints, for each layer of a template, the bottom interface, the midpoint, and the
# layer thickness (all in ilev/lev units).
#
# python diag-templates.py

import numpy as np
import netCDF4 as nc

top_template = "./Z_L55template.nc"

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

print("(0)\tNumber of input levels = "+str(tnlev))
print("(0)\tNumber of input interfaces = "+str(tnint))

ilev_thick = tilev[1:tnint] - tilev[0:tnint-1]

# NCL converts doubles to strings with %.16g
for i in range(0, tnlev):
    print("("+str(i)+")\t"+"%.16g" % tilev[1:tnint][i]+"    "+"%.16g" % tlev[i]+"     "+"%.16g" % ilev_thick[i])

tfile.close()
