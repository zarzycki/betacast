#========================================================
# This script plots ATCF formatting trajectories for NWP
# output
# Colin Zarzycki (czarzycki@psu.edu)
# Python port of plot_ATCF.ncl
#========================================================
import os
import re

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter

def ncl_tofloat(s):
    # NCL tofloat parses the leading number of a string (e.g., " 226N" -> 226.)
    return float(re.match(r"\s*[-+]?\d*\.?\d*", s).group())

yyyymmddhh = os.environ["YYYYMMDDHH"]
filename = os.environ["ATCFFILEMERGE"]
out_type = "png"

# DEBUG CONFIG
#yyyymmddhh="2021082000"
#filename="./fin-atcf/atcf.tempest."+yyyymmddhh

color_pw = False
cutwind = 0.0

#---Read in file as array of strings so we can parse each line
with open(filename) as f:
    lines = [line for line in f.read().splitlines() if line.strip() != ""]

delim = ","

#---Read fields (NCL str_get_field treats consecutive delimiters as one)
fields = [[x for x in line.split(delim) if x != ""] for line in lines]
_basin = [fld[0] for fld in fields]
_number_only = [fld[1] for fld in fields]
_number = [_basin[ii]+_number_only[ii] for ii in range(len(lines))]
_modelid = [fld[4] for fld in fields]
_fctime = [int(fld[5]) for fld in fields]
_lat = np.array([ncl_tofloat(fld[6]) for fld in fields], dtype=np.float32)/10.
_lon = np.array([ncl_tofloat(fld[7]) for fld in fields], dtype=np.float32)/10.
_wind = np.array([ncl_tofloat(fld[8]) for fld in fields], dtype=np.float32)*0.514444 # kts to m/s
_pres = np.array([ncl_tofloat(fld[9]) for fld in fields], dtype=np.float32)

# See if "S" or "W" are a substring of the lat/lon string and if so, correct
_isLatS = np.array(["S" in fld[6].strip() for fld in fields])
_isLonW = np.array(["W" in fld[7].strip() for fld in fields])

_lat = np.where(_isLatS, -_lat, _lat)
_lon = np.where(_isLonW, 360.-_lon, _lon)

nstorms = len(set(_number))
nmodels = len(set(_modelid))
ntimes = len(set(_fctime))

print("nstorms "+str(nstorms)+"    nmodels "+str(nmodels)+"     ntimes "+str(ntimes))

lat = np.full((nmodels*nstorms, ntimes), np.nan, dtype=np.float32)
lon = lat.copy()
wind = lat.copy()
pres = lat.copy()

uniqNums = sorted(set(_number))
labelNum = [None]*len(uniqNums)
labelLat = np.full(len(uniqNums), np.nan, dtype=np.float32)
labelLon = np.full(len(uniqNums), np.nan, dtype=np.float32)
timeCounter = 0
stormCounter = 0
THISMODEL = _modelid[0]
THISSTORM = _number[0]
for ii in range(len(_lat)):
    if labelNum[uniqNums.index(THISSTORM)] is None:
        IX = uniqNums.index(THISSTORM)
        labelNum[IX] = _number_only[ii]
        labelLat[IX] = _lat[ii]
        labelLon[IX] = _lon[ii]
    if THISMODEL == _modelid[ii] and THISSTORM == _number[ii]:
        lat[stormCounter,timeCounter] = _lat[ii]
        lon[stormCounter,timeCounter] = _lon[ii]
        wind[stormCounter,timeCounter] = _wind[ii]
        pres[stormCounter,timeCounter] = _pres[ii]
        timeCounter = timeCounter + 1
    else:
        # new storm or model, start a new trajectory
        timeCounter = 0
        stormCounter = stormCounter+1
        lat[stormCounter,timeCounter] = _lat[ii]
        lon[stormCounter,timeCounter] = _lon[ii]
        wind[stormCounter,timeCounter] = _wind[ii]
        pres[stormCounter,timeCounter] = _pres[ii]
        # NCL version did not increment here, so 2nd point overwrote genesis point of every traj after the 1st
        timeCounter = timeCounter + 1
        THISMODEL = _modelid[ii]
        THISSTORM = _number[ii]

nstorms = nstorms*nmodels

print("Beginning plot sequence")

numpanels = 4
#                 minlat, maxlat, minlon, maxlon
basins = [ ( 10, 45, 260, 350),   # atlantic
           ( 10, 35,-165, -95),   # epac
           ( 10, 45, 105, 195),   # wpac
           (-60, 60,  30, 350) ]  # global

# panel heights proportional to domain aspect so all panels share a width (as gsn_panel does)
fig = plt.figure(figsize=(7., 10.))
gs = fig.add_gridspec(numpanels, 1, height_ratios=[(b[1]-b[0])/(b[3]-b[2]) for b in basins])

for zz in range(numpanels):
    mpMinLatF, mpMaxLatF, mpMinLonF, mpMaxLonF = basins[zz]

    # center projection on panel so domains crossing the dateline work
    clon = (mpMinLonF+mpMaxLonF)/2.
    proj = ccrs.PlateCarree(central_longitude=clon)
    ax = fig.add_subplot(gs[zz], projection=proj)
    ax.add_feature(cfeature.OCEAN.with_scale("110m"), facecolor="black")
    ax.add_feature(cfeature.LAND.with_scale("110m"), facecolor="#D2B48C")    # Tan
    ax.add_feature(cfeature.LAKES.with_scale("110m"), facecolor="black")
    ax.coastlines(resolution="110m", linewidth=0.5)
    if zz == 0:
        ax.add_feature(cfeature.STATES.with_scale("110m"), edgecolor="black", linewidth=0.3)  # USStates
    if mpMaxLatF-mpMinLatF <= 45:
        latspc = 10.
    else:
        latspc = 30.
    ax.set_xticks(np.arange(-180., 361., 30.), crs=ccrs.PlateCarree())
    ax.set_yticks(np.arange(-90., 91., latspc), crs=ccrs.PlateCarree())
    ax.set_extent([mpMinLonF-clon, mpMaxLonF-clon, mpMinLatF, mpMaxLatF], crs=proj)
    ax.xaxis.set_major_formatter(LongitudeFormatter(degree_symbol=""))
    ax.yaxis.set_major_formatter(LatitudeFormatter(degree_symbol=""))
    ax.tick_params(labelsize=7)

    SSthresh_wind = np.array([17.,33.,43.,49.,58.])

    if color_pw:
        # NOTE, these are negative pressures
        SSthresh_wind = ((SSthresh_wind/4.4/0.514444)**(1/0.76) - 1010)

    linethicknesses = np.array([0.5,1.0,1.0,1.0,1.0,1.0])
    linethicknesses = linethicknesses*4.5

    allcolors = ["red", "#CD6600", "#FFA500", "#FFFF00", "#00FF00", "blue"]   # red, darkorange3, orange, yellow, green, blue

    for i in range(nstorms):
        for j in range(ntimes):

            if color_pw:
                thisWind = -pres[i,j]
            else:
                thisWind = wind[i,j]

            if not np.isnan(thisWind):
                if thisWind < SSthresh_wind[0]:
                    gsLineColor = allcolors[5]
                    gsLineThicknessF = linethicknesses[0]
                elif thisWind >= SSthresh_wind[0] and thisWind < SSthresh_wind[1]:
                    gsLineColor = allcolors[4]
                    gsLineThicknessF = linethicknesses[1]
                elif thisWind >= SSthresh_wind[1] and thisWind < SSthresh_wind[2]:
                    gsLineColor = allcolors[3]
                    gsLineThicknessF = linethicknesses[2]
                elif thisWind >= SSthresh_wind[2] and thisWind < SSthresh_wind[3]:
                    gsLineColor = allcolors[2]
                    gsLineThicknessF = linethicknesses[3]
                elif thisWind >= SSthresh_wind[3] and thisWind < SSthresh_wind[4]:
                    gsLineColor = allcolors[1]
                    gsLineThicknessF = linethicknesses[4]
                else:
                    gsLineColor = allcolors[0]
                    gsLineThicknessF = linethicknesses[5]

                if j != 0:
                    # Geodetic == mpGreatCircleLinesOn
                    ax.plot([lon[i,j],lon[i,j-1]], [lat[i,j],lat[i,j-1]], color=gsLineColor,
                            linewidth=0.35*gsLineThicknessF, solid_capstyle="round", transform=ccrs.Geodetic())

    for jj in range(len(uniqNums)):
        ax.text(labelLon[jj]+0.5, labelLat[jj]-0.5, labelNum[jj], color="white", ha="left", va="top", fontsize=9,
                transform=ccrs.PlateCarree(), clip_on=True)

#******************
# Plot a legend   *
#******************

markerSize = 6

if cutwind == 0:
    lgMarkerColors = ["red", "#CD6600", "#FFA500", "#FFFF00", "#00FF00", "blue"]   # red, darkorange3, orange, yellow, green, blue
    legend_labels = ["Cat. 4/5","Cat. 3","Cat. 2","Cat. 1","TS","TD"]
elif cutwind == 17:
    lgMarkerColors = ["red", "#CD6600", "#FFA500", "#EEEE00", "#00CD00"]           # red, darkorange3, orange, yellow2, green3
    legend_labels = ["Category 4/5","Category 3","Category 2","Category 1","Tropical Storm"]

handles = [Line2D([], [], linestyle="none", marker="o", markersize=markerSize, color=c) for c in lgMarkerColors]

# add legend to last panel, NCL draws the first legend item at the bottom
ax.legend(handles[::-1], legend_labels[::-1], loc="upper left", fontsize=6, framealpha=1.0, facecolor="white",
          edgecolor="black", fancybox=False, borderaxespad=0.3, handletextpad=0.2, labelspacing=0.3)

print("Drawing plot")

fig.tight_layout()
fig.savefig("enstraj."+yyyymmddhh+"."+out_type, dpi=150)
