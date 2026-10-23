#!/usr/bin/env python
# Port of gridplot.ncl -- plot an Exodus (.g) mesh on an orthographic projection
# Usage: python gridplot.py <gridfile>
import os
import sys
import numpy as np
import netCDF4
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
import cartopy.crs as ccrs
import cartopy.feature as cfeature

out_type = "png"

#gridfile = "test_NewOrleans_TC_grid_ne32x8.g"
if len(sys.argv) != 2:
    print("Usage: python gridplot.py <gridfile>")
    sys.exit(1)
gridfile = sys.argv[1]

if not os.path.isfile(gridfile):
    print("Can not find " + gridfile + "!")
    sys.exit(1)

# name for writing output file
gridname = os.path.basename(gridfile)

print("Plotting mesh from " + gridfile)
# (NCL symlinks gridfile -> gridfile.nc; netCDF4 opens the .g directly)

# Set constants
pi4 = np.arctan(1.0)
pi2 = np.arccos(0.0)
pi = pi2 * 2.0

f = netCDF4.Dataset(gridfile, "r")
num_el_blk = len(f.dimensions["num_el_blk"])

if (num_el_blk != 6) and (num_el_blk != 1):
    print("Neither cube nor sphere, can't plot!")
    sys.exit(1)

os.makedirs("./pdfs", exist_ok=True)

if num_el_blk == 6:
    IsCube = True
else:
    IsCube = False

# cartesian coords:
coordx = np.asarray(f.variables["coord"][0, :], dtype="f8")
coordy = np.asarray(f.variables["coord"][1, :], dtype="f8")
coordz = np.asarray(f.variables["coord"][2, :], dtype="f8")

# convert to lat/lon (project from cube to sphere!)
if IsCube:
    tol = 1e-8
    print("Projecting coordinates onto sphere...")
    n_nodes = len(coordx)
    r = np.sqrt(coordx * coordx + coordy * coordy + coordz * coordz)
    lat = np.zeros(n_nodes, dtype="f8")
    lon = np.zeros(n_nodes, dtype="f8")
    face = np.zeros(n_nodes, dtype="f8")
    for i in range(n_nodes):
        face[i] = -1
        if (coordx[i] <= -pi4 + tol) and (face[i] == -1):
            face[i] = 2
            alpha = -coordy[i]
            beta = coordz[i]
        if (coordx[i] >= pi4 - tol) and (face[i] == -1):
            face[i] = 0
            alpha = coordy[i]
            beta = coordz[i]
        if (coordy[i] <= -pi4 + tol) and (face[i] == -1):
            face[i] = 3
            alpha = coordx[i]
            beta = coordz[i]
        if (coordy[i] >= pi4 - tol) and (face[i] == -1):
            face[i] = 1
            alpha = -coordx[i]
            beta = coordz[i]
        if (coordz[i] <= -pi4 + tol) and (face[i] == -1):
            face[i] = 5
            alpha = coordy[i]
            beta = coordx[i]
        if (coordz[i] >= pi4 - tol) and (face[i] == -1):
            face[i] = 4
            alpha = coordy[i]
            beta = -coordx[i]
        if face[i] == -1:
            print("Error: face not assigned to a node!")
            print(str(i) + " " + str(coordx[i]) + " " + str(coordy[i]) + " " + str(coordz[i]))
            sys.exit(1)
        if face[i] < 4:
            lon[i] = alpha + pi2 * face[i]
            lat[i] = np.arctan(np.tan(beta) * np.cos(alpha))
        if face[i] == 4:
            lon[i] = np.arctan2(np.tan(alpha), -np.tan(beta))
            if np.tan(alpha) * np.tan(alpha) + np.tan(beta) * np.tan(beta) < tol:
                lat[i] = pi2
            else:
                lat[i] = np.arctan(1.0 / np.sqrt(np.tan(alpha) * np.tan(alpha) + np.tan(beta) * np.tan(beta)))
        if face[i] == 5:
            lon[i] = np.arctan2(np.tan(alpha), np.tan(beta))  # -pi
            if np.tan(alpha) * np.tan(alpha) + np.tan(beta) * np.tan(beta) < tol:
                lat[i] = -pi2
            else:
                lat[i] = -np.arctan(1.0 / np.sqrt(np.tan(alpha) * np.tan(alpha) + np.tan(beta) * np.tan(beta)))
    print("... done!")
else:
    lat = np.arcsin(coordz)
    lon = np.arctan2(coordy, coordx)

lat = lat * 180.0 / pi
lon = lon * 180.0 / pi

print("number of verticies = " + str(len(lon)))
print("lat min/max = " + str(lat.min()) + " " + str(lat.max()))

for dimname, dim in f.dimensions.items():
    if dimname == "num_elem":
        print("dimname[i]=" + dimname + " " + str(len(dim)))
        nelem = len(dim)

map = np.full((nelem, 4), -99, dtype="i4")

count = 0
fcnt = f.variables["connect1"].shape
print("dimsizes(map1)=" + str(fcnt))
n = fcnt[0]
map[count:count + n, :] = f.variables["connect1"][:]
count = count + n

if IsCube:
    fcnt = f.variables["connect2"].shape
    print("dimsizes(map2)=" + str(fcnt))
    n = fcnt[0]
    map[count:count + n, :] = f.variables["connect2"][:]
    count = count + n

    fcnt = f.variables["connect3"].shape
    print("dimsizes(map2)=" + str(fcnt))
    n = fcnt[0]
    map[count:count + n, :] = f.variables["connect3"][:]
    count = count + n

    fcnt = f.variables["connect4"].shape
    print("dimsizes(map2)=" + str(fcnt))
    n = fcnt[0]
    map[count:count + n, :] = f.variables["connect4"][:]
    count = count + n

    fcnt = f.variables["connect5"].shape
    print("dimsizes(map2)=" + str(fcnt))
    n = fcnt[0]
    map[count:count + n, :] = f.variables["connect5"][:]
    count = count + n

    fcnt = f.variables["connect6"].shape
    print("dimsizes(map2)=" + str(fcnt))
    n = fcnt[0]
    map[count:count + n, :] = f.variables["connect6"][:]
    count = count + n

nelem = count
print(map.shape)
print("number of elements = " + str(nelem))

# Orthographic Projection
# NEUS
#40.7128 N, 74.0060 W
mpCenterLatF = 41.
mpCenterLonF = 286.

# USA
# mpCenterLatF = 45.
# mpCenterLonF = 260.

# ATLANTIC
# mpCenterLatF = 10.
# mpCenterLonF = 330.

# EUROPE
# mpCenterLatF = 45.
# mpCenterLonF = 20.

# Africa
# mpCenterLatF = 15.
# mpCenterLonF = 120.

# Pacfic
# mpCenterLatF = 15.
# mpCenterLonF = 100.

# mpCenterLatF = 15.
# mpCenterLonF = 200.

# RICH WESTPAC 30x4
# mpCenterLatF = 10.
# mpCenterLonF = 150.

# RICH INDIAN 30x4
# mpCenterLatF = 0.
# mpCenterLonF = 70.

# NORTH POLE
# mpCenterLatF = 90.
# mpCenterLonF = 0.

# SOUTH POLE
# mpCenterLatF = -90.
# mpCenterLonF = 0.

# Australasia
# mpCenterLatF = -27.
# mpCenterLonF = 147.

# New Zealand
# mpCenterLatF = -41.
# mpCenterLonF = 175.

fig = plt.figure(figsize=(8, 8))
# vpXF=0.05, vpYF=0.9, vpWidthF=0.9, vpHeightF=0.8 (NDC; YF is the top edge)
ax = fig.add_axes([0.05, 0.9 - 0.8, 0.9, 0.8],
                  projection=ccrs.Orthographic(central_longitude=mpCenterLonF,
                                               central_latitude=mpCenterLatF))
ax.set_global()

plot_outlines = True
if plot_outlines:
    ax.add_feature(cfeature.LAND, facecolor="tan", edgecolor="none")
    ax.add_feature(cfeature.OCEAN, facecolor="lightblue", edgecolor="none")
    ax.add_feature(cfeature.LAKES, facecolor="blue", edgecolor="none")
    ax.add_feature(cfeature.COASTLINE, edgecolor="black", linewidth=0.5)
    ax.add_feature(cfeature.BORDERS, edgecolor="black", linewidth=0.5)  # National
    ax.spines["geo"].set_visible(False)  # mpPerimOn = False
else:
    ax.set_facecolor("white")
    ax.spines["geo"].set_visible(False)

# polygon resources
gsLineThicknessF = 1.0
gsLineColor = "black"

j1 = np.arange(0, 4)
j2 = (j1 + 1) % 4

xlon = np.zeros(5, dtype="f8")
xlat = np.zeros(5, dtype="f8")

segs = []  # collect all polygons and draw once (much faster than one plot call per element)

for i in range(nelem):
    if i % 500 == 0:
        print("i = " + str(i) + "/" + str(nelem - 1))
    xlon[0:4] = lon[map[i, :] - 1]
    xlat[0:4] = lat[map[i, :] - 1]

    for j in range(3):
        if abs(xlon[j + 1] - xlon[j]) > 180.0:
            if xlon[j + 1] > xlon[j]:
                xlon[j] = xlon[j] + 360.
            else:
                xlon[j + 1] = xlon[j + 1] + 360.

    xlon[4] = xlon[0]
    xlat[4] = xlat[0]
    segs.append(np.column_stack((xlon, xlat)))

# True = draw edges as great circles (mpGreatCircleLinesOn); very slow
# True recovers old NCL workflow, though, and doesn't draw lines as straight edges
great_circle_lines = False

if great_circle_lines:
    ax.add_collection(LineCollection(segs, colors=gsLineColor, linewidths=gsLineThicknessF * 0.25,
                                     transform=ccrs.Geodetic()))
else:
    # project every polygon to map coordinates in one vectorized call (cartopy's per-path
    # Geodetic projection is very slow); elements are small so straight segments ~ great circles
    segs = np.array(segs)
    pts = ax.projection.transform_points(ccrs.PlateCarree(), segs[:, :, 0].ravel(), segs[:, :, 1].ravel())
    xy = pts[:, :2].reshape(segs.shape)
    xy = xy[np.isfinite(xy).all(axis=(1, 2))]  # drop polygons on the far side of the globe
    ax.add_collection(LineCollection(xy, colors=gsLineColor, linewidths=gsLineThicknessF * 0.25))

#*****************************************
# plot min extrema over original data
#*****************************************

# polymarkers for cordex corners
# polyres: circle marker, size 8, red
#ax.plot(-79.373333, 43.741667, "o", color="red", markersize=8, transform=ccrs.PlateCarree())
#ax.plot(-79.373333, 47.741667, "o", color="red", markersize=8, transform=ccrs.PlateCarree())

#ax.plot(1.92, 71.84, "o", color="red", markersize=8, transform=ccrs.PlateCarree())
#ax.plot(64.4, 66.65, "o", color="red", markersize=8, transform=ccrs.PlateCarree())
#ax.plot(338.23, 42.36, "o", color="red", markersize=8, transform=ccrs.PlateCarree())
#ax.plot(9.75, 49.68, "o", color="red", markersize=8, transform=ccrs.PlateCarree())
#ax.plot(44.77, 46.72, "o", color="red", markersize=8, transform=ccrs.PlateCarree())
#ax.plot(350.01, 22.20, "o", color="red", markersize=8, transform=ccrs.PlateCarree())
#ax.plot(12.48, 27.34, "o", color="red", markersize=8, transform=ccrs.PlateCarree())
#ax.plot(36.30, 25.36, "o", color="red", markersize=8, transform=ccrs.PlateCarree())

if 0 == 1:
    # gres: red fill, 0.5 opacity
    #polylat = [40, 24, 24, 40]   # SE USA
    #polylon = [265, 265, 287, 287]  # SE USA
    #polylat = [-5., 25., 25., -5.]
    #polylon = [179., 179., 120., 120.]
    # AUSTRALIA CORDEX (9)
    #polylat = [8.76, 3.87, 12.21, -18.03, -24.26, -27.90, -44.28, -52.36, -39.25]
    #polylon = [110.19, 146.16, 182.02, 101.41, 147.63, 199.57, 89.25, 150.03, 206.57]
    # RICO / ATOMIC
    polylat = [0., 35., 35., 0.]
    polylon = [-110., -110., 30., 30.]
    ax.fill(polylon, polylat, facecolor="red", edgecolor="red", alpha=0.5,
            transform=ccrs.Geodetic())

    #TLC (110.19; 8.76)
    #CNB (146.16; 3.87)
    #TRC (182.02; 12.21)
    #CWB (101.41; -18.03)
    #CPD (147.63; -24.26)
    #CEB (199.57; -27.90)
    #BLC (89.25; -44.28)
    #CSB (150.03; -52.36)
    #BRC (206.57; -39.25)

if 1 == 0:
    # AUSTRALIA CORDEX (9)
    polylat = [-25., 35., 35., -25.]
    polylon = [250., 250., 30., 30.]
    ax.plot(polylon, polylat, "o", color="red", markersize=10, transform=ccrs.PlateCarree())

    #TLC (110.19; 8.76)
    #CNB (146.16; 3.87)
    #TRC (182.02; 12.21)
    #CWB (101.41; -18.03)
    #CPD (147.63; -24.26)
    #CEB (199.57; -27.90)
    #BLC (89.25; -44.28)
    #CSB (150.03; -52.36)
    #BRC (206.57; -39.25)

fig.savefig("./pdfs/" + gridname + "_stereo." + out_type, dpi=300)
plt.close(fig)
f.close()
