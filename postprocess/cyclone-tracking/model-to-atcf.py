#========================================================
# This script reads in a TempestExtremes trajectory file
# and produces an ATCF-formatted output given matching tcvitals
# Python port of model-to-atcf.ncl
#========================================================
import os
import sys
from datetime import datetime

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "functions"))
from getTrajectories import getTrajectories

#;;; testing (assuming TE has already generated a track file and TCvitals exist at this time
#env TRAJFILE="trajectories.txt.forecast_nhemitc_30_x4_CAM5_L30" \
#    TCVITFILE="./fin-tcvitals/tcvitals.2024052600" \
#    ATCFFILE="atcf.forecast_nhemitc_30_x4_CAM5_L30.2024052600" \
#    ATCFTECH="CAM5" \
#    python model-to-atcf.py

def gc_latlon(lat1, lon1, lat2, lon2):
    # great circle distance (degrees), equivalent to NCL gc_latlon(...,2,2)
    d2r = np.pi/180.
    lat1, lon1, lat2, lon2 = [float(x)*d2r for x in (lat1, lon1, lat2, lon2)]
    hav = np.sin((lat2-lat1)/2.)**2 + np.cos(lat1)*np.cos(lat2)*np.sin((lon2-lon1)/2.)**2
    return 2.*np.arcsin(np.sqrt(hav))/d2r

def ncl_round(x):
    # NCL round(x,3) rounds half away from zero and returns an integer
    return int(np.sign(x)*np.floor(np.abs(x)+0.5))

thefile = os.environ["TRAJFILE"]
vitalsfile = os.environ["TCVITFILE"]
outfile = os.environ["ATCFFILE"]
TECH = os.environ["ATCFTECH"]
print("thefile: "+thefile)
print("vitalsfile: "+vitalsfile)
print("outfile: "+outfile)
print("TECH: "+TECH)

#--------------------------------------------------------

TIMETHRESH = 12.0 # max number of hours init in tempest can be off from vitals
DISTTHRESH = 5.1  # max horiz distance (GC deg) init in tempest can be off from vitals (slightly > than mergedist)

output_type = "png"
cutwind = 0.0

ncols = 9
headerDelimStr = "start"
isUnstruc = False
nstorms, ntimes, ncols, traj_data = getTrajectories(thefile, ncols, headerDelimStr, isUnstruc)
# NCL reads trajectories as single precision, keep float32 so ATCF ints match
traj_data = traj_data.astype(np.float32)

lon = traj_data[1,:,:]
lat = traj_data[2,:,:]
wind = traj_data[4,:,:]
pres = traj_data[3,:,:]/100.
yyyy = traj_data[5,:,:]
mm = traj_data[6,:,:]
dd = traj_data[7,:,:]
hh = traj_data[8,:,:]

print("Sorted (S->N) genesis points from TempestExtremes tracking")
sorted_lats = np.argsort(lat[:,0], kind="stable")
for ix in sorted_lats:
    print(f"{lat[ix,0]:g} {lon[ix,0]:g} {pres[ix,0]:g}")

tempesthit = [""]*nstorms
#--------------------------------------------------------

with open(vitalsfile) as f:
    data = [line for line in f.read().splitlines() if line.strip() != ""]
v_id = [line[5:8] for line in data]
v_name = [line[9:19] for line in data]
v_yyyy = [int(line[19:23]) for line in data]
v_mm = [int(line[23:25]) for line in data]
v_dd = [int(line[25:27]) for line in data]
v_hh = [int(line[28:30]) for line in data]
v_lat = np.array([float(line[33:36]) for line in data], dtype=np.float32)/10.
v_latdir = [line[36:37] for line in data]
v_lon = np.array([float(line[38:42]) for line in data], dtype=np.float32)/10.
v_londir = [line[42:43] for line in data]

print("Parsed ATCF data")
for ii in range(len(data)):
    print(f"{v_id[ii]} {v_name[ii]} {v_yyyy[ii]} {v_mm[ii]} {v_dd[ii]} {v_hh[ii]} {v_lat[ii]:g} {v_latdir[ii]} {v_lon[ii]:g} {v_londir[ii]}")

#NHC  12L JOSE      20170920 1200 388N 0702W 035 036 0976 1012 0556 28 093 0334 0334 0297 0297 D

num_vit = len(v_id)

for ii in range(num_vit):
    if v_latdir[ii] == "S":
        v_lat[ii] = -v_lat[ii]
    if v_londir[ii] == "W":
        v_lon[ii] = 360.-v_lon[ii]

print("Corrected ATCF lat/lons")
for ii in range(num_vit):
    print(f"{v_lat[ii]:g} {v_lon[ii]:g}")

print("******* Looking for storms!")
for ii in range(num_vit):
    print("    "+str(ii)+" "+v_name[ii])
    v_single_time = datetime(v_yyyy[ii], v_mm[ii], v_dd[ii], v_hh[ii])
    for jj in range(nstorms):
        init_single_time = datetime(int(yyyy[jj,0]), int(mm[jj,0]), int(dd[jj,0]), int(hh[jj,0]))
        diff_single_time = (init_single_time-v_single_time).total_seconds()/3600.   # hours b/w tempest init and tcvitals
        diff_single_dist = gc_latlon(v_lat[ii], v_lon[ii], lat[jj,0], lon[jj,0])     # dist b/w tempest init and tcvitals
        if abs(diff_single_time) <= TIMETHRESH and diff_single_dist < DISTTHRESH:
            print("Storm "+str(jj)+"  in tempest file is "+v_name[ii])
            print(f"---> diff_single_time: {diff_single_time:g}    diff_single_dist: {diff_single_dist:g}")
            tempesthit[jj] = v_id[ii]
            for zz in range(ntimes):
                if not np.isnan(lat[jj,zz]):
                    if zz != 0:
                        disttrav = gc_latlon(lat[jj,zz-1], lon[jj,zz-1], lat[jj,zz], lon[jj,zz])
                    else:
                        disttrav = 0.
                    print(f"{v_id[ii]} {v_name[ii]} {lat[jj,zz]:g} {lon[jj,zz]:g} {pres[jj,zz]:g} {int(yyyy[jj,zz])} {int(mm[jj,zz])} {int(dd[jj,zz])} {int(hh[jj,zz])} dist: {disttrav:g}")
            break

# WRITE ATCF FILE
foundstorms = [jj for jj in range(nstorms) if tempesthit[jj] != ""]
if len(foundstorms) > 0:
    INITDATE = v_yyyy[0]*1000000+v_mm[0]*10000+v_dd[0]*100+v_hh[0]
    inittime = datetime(v_yyyy[0], v_mm[0], v_dd[0], v_hh[0])
    ms_to_kts = 1.94384
    delim = ", "
    TECHNUM = "03"
    TY = "XX"
    FILLER = "34, NEQ,    0,    0,    0,    0, "
    with open(outfile, "a") as fout:
        for ix in foundstorms:
            for zz in range(ntimes):
                if not np.isnan(hh[ix,zz]) and int(hh[ix,zz]) % 6 == 0:

                    temptime = datetime(int(yyyy[ix,zz]), int(mm[ix,zz]), int(dd[ix,zz]), int(hh[ix,zz]))
                    LEADTIME = int((temptime-inittime).total_seconds()/3600.)

                    STORMID = int(tempesthit[ix][0:2])
                    tmpbasin = tempesthit[ix][2:3]
                    if tmpbasin == "L":
                        STORMBASIN = "AL"
                    elif tmpbasin == "E":
                        STORMBASIN = "EP"
                    elif tmpbasin == "W":
                        STORMBASIN = "WP"
                    elif tmpbasin == "C":
                        STORMBASIN = "CP"
                    elif tmpbasin == "B" or tmpbasin == "A":
                        STORMBASIN = "IO"
                    elif tmpbasin == "Q" or tmpbasin == "P" or tmpbasin == "S":
                        STORMBASIN = "SH"
                    else:
                        STORMBASIN = "XX"

                    templat = lat[ix,zz]
                    if templat >= 0.0:
                        LATDIR = "N"
                    else:
                        LATDIR = "S"
                    LAT = int(abs(templat)*np.float32(10.))

                    templon = lon[ix,zz]
                    if templon >= 180.0:
                        LONDIR = "W"
                        templon = templon-np.float32(360.)
                    else:
                        LONDIR = "E"
                    LON = int(abs(templon)*np.float32(10.))

                    WIND = ncl_round(wind[ix,zz]*np.float32(ms_to_kts))
                    PRES = ncl_round(pres[ix,zz])

                    #AL, 01, 2012052106, 03, AP06,  60, 375N,  706W,  18, 1013, XX,  34, NEQ,    0,    0,    0,    0,
                    ATCFSTRING = STORMBASIN+delim+"%0.2i" % STORMID+delim+"%0.10i" % INITDATE+delim+ \
                        TECHNUM+delim+TECH+delim+"%0.3i" % LEADTIME+delim+ \
                        "%0.3i" % LAT+LATDIR+delim+"%0.4i" % LON+LONDIR+delim+ \
                        "%0.3i" % WIND+delim+"%0.4i" % PRES+delim+TY+delim+FILLER

                    fout.write(ATCFSTRING+"\n")
else:
    print("Found no matches to tempest + vitals, cannot write ATCF")

#--------------------------------------------------------
print("Beginning plot sequence")

numpanels = 4
#                 minlat, maxlat, minlon, maxlon
basins = [ (  5, 45, 260, 355),   # atlantic
           (  5, 35,-165, -90),   # epac
           (  5, 45, 100, 195),   # wpac
           (  0, 75,  50, 350) ]  # northern hemisphere

# panel heights proportional to domain aspect so all panels share a width (as gsn_panel does)
fig = plt.figure(figsize=(7., 10.))
gs = fig.add_gridspec(numpanels, 1, height_ratios=[(b[1]-b[0])/(b[3]-b[2]) for b in basins])

for zz in range(numpanels): # different basins
    mpMinLatF, mpMaxLatF, mpMinLonF, mpMaxLonF = basins[zz]

    # center projection on panel so domains crossing the dateline work
    clon = (mpMinLonF+mpMaxLonF)/2.
    proj = ccrs.PlateCarree(central_longitude=clon)
    ax = fig.add_subplot(gs[zz], projection=proj)
    ax.add_feature(cfeature.OCEAN.with_scale("110m"), facecolor="#BFEFFF")   # LightBlue1
    ax.add_feature(cfeature.LAND.with_scale("110m"), facecolor="#D2B48C")    # Tan
    ax.add_feature(cfeature.LAKES.with_scale("110m"), facecolor="#BFEFFF")
    ax.coastlines(resolution="110m", linewidth=0.5)
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

    # We need to use thinner lines for the global map
    linethicknesses = np.array([0.5,1.0,1.3,1.9,2.3,2.6])
    linethicknesses = 4.0*linethicknesses

    for i in range(nstorms):
        if tempesthit[i] != "":
            for j in range(ntimes):

                thisWind = wind[i,j]

                if not np.isnan(thisWind):

                    if thisWind < 17:
                        gsLineColor = "blue"
                        gsLineThicknessF = linethicknesses[0]
                    elif thisWind >= 17 and thisWind < 33:
                        gsLineColor = "#00CD00"   # green3
                        gsLineThicknessF = linethicknesses[1]
                    elif thisWind >= 33 and thisWind < 43:
                        gsLineColor = "#EEEE00"   # yellow2
                        gsLineThicknessF = linethicknesses[2]
                    elif thisWind >= 43 and thisWind < 49:
                        gsLineColor = "#FFA500"   # orange
                        gsLineThicknessF = linethicknesses[3]
                    elif thisWind >= 49 and thisWind < 58:
                        gsLineColor = "#CD6600"   # darkorange3
                        gsLineThicknessF = linethicknesses[4]
                    else:
                        gsLineColor = "red"
                        gsLineThicknessF = linethicknesses[5]

                    if j != 0:
                        if thisWind >= cutwind:
                            # Geodetic == mpGreatCircleLinesOn
                            ax.plot([lon[i,j],lon[i,j-1]], [lat[i,j],lat[i,j-1]], color=gsLineColor,
                                    linewidth=0.35*gsLineThicknessF, solid_capstyle="round", transform=ccrs.Geodetic())

            label = tempesthit[i]
            ax.text(lon[i,0]+0.5, lat[i,0]-0.5, label, ha="left", va="top", fontsize=6,
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
fig.savefig("trajs."+thefile+"."+output_type, dpi=150)
