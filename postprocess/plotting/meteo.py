#**************************************
#  meteo.py
#  Python port of meteo.ncl (based on NCL meteo_2.ncl)
#
# Concepts illustrated:
#   - Drawing a meteogram
#   - Explicitly setting tickmarks and labels on the bottom X axis
#   - Drawing a bar chart
#   - Changing the width and height of a plot
#   - Attaching plots along the X axis (gsnAttachPlotsXAxis = True)
#
#***********************************************
import os
import glob

import numpy as np
import xarray as xr
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

def gc_latlon(lat1, lon1, lat2, lon2):
    # great circle distance (degrees), equivalent to NCL gc_latlon(...,2,4)
    d2r = np.pi/180.
    lat1, lon1, lat2, lon2 = lat1*d2r, lon1*d2r, lat2*d2r, lon2*d2r
    hav = np.sin((lat2-lat1)/2.)**2 + np.cos(lat1)*np.cos(lat2)*np.sin((lon2-lon1)/2.)**2
    return 2.*np.arcsin(np.sqrt(hav))/d2r

all_files = sorted(glob.glob(os.path.expanduser("~/scratch/output/forecast_nhemitc_30_x4_CAM5_L30/run/2020051600/*h0*nc")))
fall = [xr.open_dataset(f, decode_times=False) for f in all_files]

TESTCIT = ["STATECOLLEGE"]
TESTSTR = ["KUNV"]
TESTLON = [-77.848611]
TESTLAT = [ 40.849167]
#TESTLON = [-96.7898]
#TESTLAT = [ 46.8772]

f1 = fall[0]
lat = f1["lat"].values
lon = f1["lon"].values

for ii in range(len(TESTLAT)):
    gcdist = gc_latlon(TESTLAT[ii], TESTLON[ii], lat, lon)
    ix = np.argmin(gcdist)
    print(ix)
    print(lat[ix])
    print(lon[ix])

tempht = np.concatenate([f["TREFHT"].values[:,ix] for f in fall])
tempht = (tempht-273.15)*9.0/5.0+32.0
rain03 = np.concatenate([f["PRECT"].values[:,ix] for f in fall])
rain03 = rain03*141732*3
flut = np.concatenate([f["FLUT"].values[:,ix] for f in fall])
ptype = np.concatenate([f["PTYPE"].values[:,ix] for f in fall])
#z3 = np.concatenate([f["Z3"].values[:,:,ix] for f in fall]).T   # (lev, time)
taus = np.concatenate([f["time"].values for f in fall])
levels = fall[0]["lev"].values
print(taus)
print(levels)

#***********************************************
#  open a workstation
#***********************************************
fig = plt.figure(figsize=(8., 8.))      # send graphics to PNG file

#***********************************************
# Here are the resources for the flut (base) plot
#***********************************************
# x location, y location (top), width, height
flut_vp = [0.12, 0.90, 0.84, 0.48]
ax_flut = fig.add_axes([flut_vp[0], flut_vp[1]-flut_vp[3], flut_vp[2], flut_vp[3]])
ax_flut.set_xlim(taus.min(), taus.max())          # max value on x-axis
ax_flut.set_ylim(flut.min()-10., flut.max()+10.)  # min/max value on y-axis
ax_flut.set_ylabel("FLUT")                        # set y-axis string
ax_flut.plot(taus, flut, color="black", linewidth=1.0)   # line thickness/color

#***********************************************
# rain histogram only resources
#***********************************************
# attached plots keep their own height, stacked directly below the base plot
rain_vp = [flut_vp[0], flut_vp[1]-flut_vp[3], flut_vp[2], 0.12]
ax_rain = fig.add_axes([rain_vp[0], rain_vp[1]-rain_vp[3], rain_vp[2], rain_vp[3]])
ax_rain.set_ylim(0.0, 0.6)                        # min/max value on y-axis
ax_rain.set_xlim(0.0, taus.max())                 # min/max value on x-axis
ax_rain.set_yticks([0.00, 0.10 , 0.20 , 0.30 , 0.40 , 0.50 , 0.60 ])       # location of explicit labels
ax_rain.set_yticklabels([ "" ,"0.10","0.20","0.30","0.40","0.50",  ""  ])
ax_rain.set_ylabel("3hr rain total")              # Y axis label.

prect_bar_colors = ["none"]*len(taus)
prect_bar_cmap = ["dodgerblue","none","#00FF00","orange","red"]   # NCL "" -> none, NCL green
print(ptype)
for ii in range(len(taus)):
    prect_bar_colors[ii] = prect_bar_cmap[int(ptype[ii])]
# bar chart above ref line 0.0
ax_rain.bar(taus, rain03, width=taus[1]-taus[0], color=prect_bar_colors, edgecolor="black", linewidth=0.5)

#***********************************************
# temperature at surface only resources
#***********************************************
tempsfc_vp = [rain_vp[0], rain_vp[1]-rain_vp[3], rain_vp[2], 0.12]
ax_tsfc = fig.add_axes([tempsfc_vp[0], tempsfc_vp[1]-tempsfc_vp[3], tempsfc_vp[2], tempsfc_vp[3]])
ax_tsfc.set_xlim(taus.min(), taus.max())              # max value on x-axis
ax_tsfc.set_ylim(tempht.min()-2., tempht.max()+2.)    # min/max value on y-axis
ax_tsfc.set_ylabel("Temp at 2m")                      # set y-axis string
ax_tsfc.plot(taus, tempht, color="red", linewidth=1.0)   # line thickness/color

#********************************************************
# attach the plots: only the bottom plot gets x-axis labels
#********************************************************
for ax in (ax_flut, ax_rain, ax_tsfc):
    ax.tick_params(direction="in", top=False, right=True, labelsize=11)   # turn off the top tickmarks
    ax.minorticks_on()
    ax.tick_params(axis="x", which="minor", bottom=False)                # No minor tick marks.
for ax in (ax_flut, ax_rain):
    ax.tick_params(axis="x", labelbottom=False)
ax_tsfc.xaxis.set_major_formatter(matplotlib.ticker.FormatStrFormatter("%.1f"))
ax_tsfc.tick_params(axis="x", labelrotation=90)        # change label angle
fig.align_ylabels([ax_flut, ax_rain, ax_tsfc])

fig.savefig("meteo.png", dpi=100)
plt.close(fig)
