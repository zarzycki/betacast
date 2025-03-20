#!/usr/bin/env python3
import os
import sys
import argparse
import numpy as np
import xarray as xr
import netCDF4 as nc
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib import cm
import matplotlib.gridspec as gridspec
import matplotlib.ticker as mticker
import cartopy
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from cartopy.util import add_cyclic_point

module_path = os.path.abspath(os.path.join(os.path.dirname(__file__), '../..', 'py_functions'))
if module_path not in sys.path:
    sys.path.append(module_path)
import py_seedfuncs

# Parse command line arg
parser = argparse.ArgumentParser(description='Create nudging window functions')
parser.add_argument('--NLEV', type=int, help='Number of vertical levels')
args = parser.parse_args()

# Check for NLEV as commandline arg
if args.NLEV is not None:
    print(f"Number of Vertical Levels: nlev:")
    nlev = args.NLEV
    print(f"nlev={nlev}")
    print("=====")
else:
    print("Using Default number of Vertical Levels: nlev = 58")
    nlev = 58

print(f"nlev type: {type(nlev)}, value: {nlev}")

use_nclev = True
lev_template = f"../../grids/templates/L{nlev}template.nc"

try:
    with nc.Dataset(lev_template, 'r') as lt:
        nclev = lt.variables['lev'][:]
except Exception as e:
    print(f"Error reading {lev_template}: {e}")
    print("Continuing without nclev data")
    use_nclev = False
    nclev = None

pthi = "user_nl_cam"

# Read parameters from namelist
# Default for the "inverts" is missing --> if not missing, it means we do stuff below
Nudge_Hwin_Invert = py_seedfuncs.keyword_values(pthi, "Nudge_Hwin_Invert", "str", default="missing")
Nudge_Vwin_Invert = py_seedfuncs.keyword_values(pthi, "Nudge_Vwin_Invert", "str", default="missing")

Nudge_Hwin_lat0 = py_seedfuncs.keyword_values(pthi, "Nudge_Hwin_lat0", "float", default=0.0)
Nudge_Hwin_lon0 = py_seedfuncs.keyword_values(pthi, "Nudge_Hwin_lon0", "float", default=0.0)

Nudge_Vwin_lo = py_seedfuncs.keyword_values(pthi, "Nudge_Vwin_lo", "float", default=0.0)
Nudge_Vwin_hi = py_seedfuncs.keyword_values(pthi, "Nudge_Vwin_hi", "float", default=0.0)

Nudge_Hwin_lo = py_seedfuncs.keyword_values(pthi, "Nudge_Hwin_lo", "float", default=0.0)
Nudge_Hwin_hi = py_seedfuncs.keyword_values(pthi, "Nudge_Hwin_hi", "float", default=0.0)

Nudge_Hwin_latWidth = py_seedfuncs.keyword_values(pthi, "Nudge_Hwin_latWidth", "float", default=0.0)
Nudge_Hwin_lonWidth = py_seedfuncs.keyword_values(pthi, "Nudge_Hwin_lonWidth", "float", default=0.0)

Nudge_Hwin_latDelta = py_seedfuncs.keyword_values(pthi, "Nudge_Hwin_latDelta", "float", default=0.0)
Nudge_Hwin_lonDelta = py_seedfuncs.keyword_values(pthi, "Nudge_Hwin_lonDelta", "float", default=0.0)

Nudge_Vwin_Lindex = py_seedfuncs.keyword_values(pthi, "Nudge_Vwin_Lindex", "float", default=0.0)
Nudge_Vwin_Hindex = py_seedfuncs.keyword_values(pthi, "Nudge_Vwin_Hindex", "float", default=0.0)

Nudge_Vwin_Ldelta = py_seedfuncs.keyword_values(pthi, "Nudge_Vwin_Ldelta", "float", default=0.0)
Nudge_Vwin_Hdelta = py_seedfuncs.keyword_values(pthi, "Nudge_Vwin_Hdelta", "float", default=0.0)

# Process invert parameters
if Nudge_Hwin_Invert != "missing":
    val = Nudge_Hwin_Invert.strip()
    if val in [".true.", ".TRUE.", "TRUE", ".T.", ".t.", "T", "True", ".True."]:
        Nudge_Hwin_Invert = True
    else:
        Nudge_Hwin_Invert = False

    print(Nudge_Hwin_Invert)

    if Nudge_Hwin_Invert:
        Nudge_Hwin_lo = 1.0
        Nudge_Hwin_hi = 0.0
    else:
        Nudge_Hwin_lo = 0.0
        Nudge_Hwin_hi = 1.0

if Nudge_Vwin_Invert != "missing":
    val = Nudge_Vwin_Invert.strip()
    if val in [".true.", ".TRUE.", "TRUE", ".T.", ".t.", "T", "True", ".True."]:
        Nudge_Vwin_Invert = True
    else:
        Nudge_Vwin_Invert = False

    print(Nudge_Vwin_Invert)

    if Nudge_Vwin_Invert:
        Nudge_Vwin_lo = 1.0
        Nudge_Vwin_hi = 0.0
    else:
        Nudge_Vwin_lo = 0.0
        Nudge_Vwin_hi = 1.0

# Print parameter values for sanity
print(f"Nudge_Hwin_lo={Nudge_Hwin_lo}")
print(f"Nudge_Hwin_hi={Nudge_Hwin_hi}")
print(f"Nudge_Vwin_lo={Nudge_Vwin_lo}")
print(f"Nudge_Vwin_hi={Nudge_Vwin_hi}")
print(f"Nudge_Hwin_lat0={Nudge_Hwin_lat0}")
print(f"Nudge_Hwin_lon0={Nudge_Hwin_lon0}")
print(f"Nudge_Hwin_latWidth={Nudge_Hwin_latWidth}")
print(f"Nudge_Hwin_lonWidth={Nudge_Hwin_lonWidth}")
print(f"Nudge_Hwin_latDelta={Nudge_Hwin_latDelta}")
print(f"Nudge_Hwin_lonDelta={Nudge_Hwin_lonDelta}")
print(f"Nudge_Vwin_Lindex={Nudge_Vwin_Lindex}")
print(f"Nudge_Vwin_Hindex={Nudge_Vwin_Hindex}")
print(f"Nudge_Vwin_Ldelta={Nudge_Vwin_Ldelta}")
print(f"Nudge_Vwin_Hdelta={Nudge_Vwin_Hdelta}")

# Assign shorter variable names for plotting
Hlo = Nudge_Hwin_lo
Hhi = Nudge_Hwin_hi
Vlo = Nudge_Vwin_lo
Vhi = Nudge_Vwin_hi
lat0 = Nudge_Hwin_lat0
lat_width = Nudge_Hwin_latWidth
lat_delta = Nudge_Hwin_latDelta
lon0 = Nudge_Hwin_lon0
lon_width = Nudge_Hwin_lonWidth
lon_delta = Nudge_Hwin_lonDelta
levH = Nudge_Vwin_Hindex
levH_delta = Nudge_Vwin_Hdelta
levL = Nudge_Vwin_Lindex
levL_delta = Nudge_Vwin_Ldelta

# Create a Horizontal test array
nlon = 360
nlat = 180
lon = np.linspace(0, (nlon*360./(nlon+1)), nlon)
lat = np.linspace(-90., 90., nlat)

# Create 2D arrays using meshgrid
lat2d, lon2d = np.meshgrid(lat, lon, indexing='ij')
Hcoef = np.zeros((nlat, nlon), dtype=np.float32)

# Set lat/lon profiles for window function
lonx = lon2d - lon0
lonx = np.where(lonx <= -180., lonx + 360., lonx)
lonx = np.where(lonx > 180., lonx - 360., lonx)
lon0_min = -(lon_width/2.)
lon0_max = (lon_width/2.)
lon_lo = (lonx - lon0_min)/lon_delta
lon_hi = (lon0_max - lonx)/lon_delta

lat0_min = lat0 - (lat_width/2.)
lat0_max = lat0 + (lat_width/2.)
lat_lo = (lat2d - lat0_min)/lat_delta
lat_hi = (lat0_max - lat2d)/lat_delta

# Calculate min/max of RAW window function
Val1_p = ((1. + np.tanh((180. - lon0_min)/lon_delta))/2.)
Val1_0 = ((1. + np.tanh((0. - lon0_min)/lon_delta))/2.)
Val1_n = ((1. + np.tanh((-179. - lon0_min)/lon_delta))/2.)
Val2_p = ((1. + np.tanh((lon0_max - 180.)/lon_delta))/2.)
Val2_0 = ((1. + np.tanh((lon0_max - 0.)/lon_delta))/2.)
Val2_n = ((1. + np.tanh((lon0_max + 179.)/lon_delta))/2.)

Val3_p = ((1. + np.tanh((90. - lat0_min)/lat_delta))/2.)
Val3_0 = ((1. + np.tanh((lat0 - lat0_min)/lat_delta))/2.)
Val3_n = ((1. + np.tanh((-90. - lat0_min)/lat_delta))/2.)
Val4_p = ((1. + np.tanh((lat0_max - 90.)/lat_delta))/2.)
Val4_0 = ((1. + np.tanh((lat0_max - lat0)/lat_delta))/2.)
Val4_n = ((1. + np.tanh((lat0_max + 90.)/lat_delta))/2.)

Hmax = Val1_0*Val2_0*Val3_0*Val4_0
Htest = np.array([
    Val1_p*Val2_p*Val3_n*Val4_n,
    Val1_p*Val2_p*Val3_p*Val4_p,
    Val1_n*Val2_n*Val3_n*Val4_n,
    Val1_n*Val2_n*Val3_p*Val4_p
])
Hmin = np.min(Htest)

# Compute the RAW window function
for ilat in range(nlat):
    for ilon in range(nlon):
        # Note the different indexing compared to NCL because of meshgrid usage
        Hcoef[ilat, ilon] = ((1. + np.tanh(lon_lo[ilat, ilon]))/2.) * \
                           ((1. + np.tanh(lon_hi[ilat, ilon]))/2.) * \
                           ((1. + np.tanh(lat_lo[ilat, ilon]))/2.) * \
                           ((1. + np.tanh(lat_hi[ilat, ilon]))/2.)

# Scale the Window function to span the values between Hlo and Hhi
if Hmax <= Hmin:
    Hcoef = np.ones_like(Hcoef)
else:
    Hcoef = (Hcoef - Hmin) / (Hmax - Hmin)

Hcoef = Hlo + Hcoef * (Hhi - Hlo)

print("")
print(f"lon0={lon0} lon_width={lon_width} lon_delta={lon_delta}")
print(f"lat0={lat0} lat_width={lat_width} lat_delta={lat_delta}")
print(f"Hlo={Hlo} Hhi={Hhi}")
print("")
print(f" RAW(Hmin)={Hmin} RAW(Hmax)={Hmax}")
print(f" min(Hcoef)={np.min(Hcoef)} max(Hcoef)={np.max(Hcoef)}")
print("")

#----------------------
# VERTICAL WINDOW
#----------------------
# Create a Vertical test array
print(f" nlev={nlev}")
lev = np.linspace(1., nlev, nlev)
Vcoef = np.zeros(nlev, dtype=np.float32)

# Set level profiles for window function
lev0 = (levL + levH) / 2.
ilev = np.argmin(np.abs(lev - lev0))
lev_lo = (lev - levL) / levL_delta
lev_hi = (levH - lev) / levH_delta

# Compute the RAW window function
for i in range(nlev):
    Vcoef[i] = ((1. + np.tanh(lev_lo[i])) / 2.) * ((1. + np.tanh(lev_hi[i])) / 2.)

# Scale the Window function to span the values between Vlo and Vhi
Vmax = np.max(Vcoef)
Vmin = np.min(Vcoef)

if Vmax <= Vmin:
    Vcoef = np.ones_like(Vcoef)
else:
    Vcoef = (Vcoef - Vmin) / (Vmax - Vmin)

Vcoef = Vlo + Vcoef * (Vhi - Vlo)

print("")
print(f"levH={levH} levH_delta={levH_delta}")
print(f"levL={levL} levL_delta={levL_delta}")
print(f"Vlo={Vlo} Vhi={Vhi}")
print("")
print(f" RAW(Vmin)={Vmin} RAW(Vmax)={Vmax}")
print(f" min(Vcoef)={np.min(Vcoef)} max(Vcoef)={np.max(Vcoef)}")
print("")

# Create figure for plotting
fig = plt.figure(figsize=(10, 10))

# Find closest indices for lat0 and lon0
ilat = np.argmin(np.abs(lat - lat0))
ilon = np.argmin(np.abs(lon - lon0))
Hprof_lon = Hcoef[ilat, :]
Hprof_lat = Hcoef[:, ilon]

# Create panel layout
gs = gridspec.GridSpec(3, 2, height_ratios=[1.1, 2.1, 2.1], hspace=0.3, wspace=0.3)

# Plot 1: Lat0 Zonal Window Profile
ax1 = plt.subplot(gs[0, 0])
ax1.plot(lon, Hprof_lon, color='blue')
ax1.set_title("Lat$_0$ Zonal Window Profile")
ax1.set_ylim(0, 1.1)
ax1.set_xlim(np.min(lon), np.max(lon))

# Simplified longitude formatter
def lon_formatter(x, pos):
    if x == 0:
        return "0"
    elif x > 180:
        return f"{360-x:.0f}W"
    else:
        return f"{x:.0f}E"

ax1.set_xticks(np.arange(0, 361, 60))
ax1.xaxis.set_major_formatter(plt.FuncFormatter(lon_formatter))

# Grid lines
ax1.minorticks_on()
ax1.grid(which='major', linestyle='-', linewidth=0.5, alpha=1.0, color='gray')
ax1.grid(which='minor', linestyle=':', linewidth=0.5, alpha=1.0, color='gray')

# Plot 2: Lon0 Meridional Window Profile
ax2 = plt.subplot(gs[0, 1])
ax2.plot(lat, Hprof_lat, color='blue')
ax2.set_title("Lon$_0$ Meridional Window Profile")
ax2.set_ylim(0, 1.1)
ax2.set_xlim(np.min(lat), np.max(lat))

# Simplified latitude formatter
def lat_formatter(x, pos):
    if x == 0:
        return "0"
    elif x < 0:
        return f"{abs(x):.0f}S"
    else:
        return f"{x:.0f}N"

ax2.set_xticks(np.arange(-90, 91, 30))
ax2.xaxis.set_major_formatter(plt.FuncFormatter(lat_formatter))

# Grid lines
ax2.minorticks_on()
ax2.grid(which='major', linestyle='-', linewidth=0.5, alpha=1.0, color='gray')
ax2.grid(which='minor', linestyle=':', linewidth=0.5, alpha=1.0, color='gray')

# Plot 3: Horizontal Window Contour Map
ax3 = plt.subplot(gs[1, :], projection=ccrs.PlateCarree())
cmap = plt.colormaps['tab20b']
levels = np.linspace(0.0, 1.0, 11)

# Add cyclic point to avoid gap at the dateline
Hcoef_cyclic, lon_cyclic = add_cyclic_point(Hcoef, coord=lon, axis=1)
lon_mesh, lat_mesh = np.meshgrid(lon_cyclic, lat)

# Plot the contour map
cs = ax3.contourf(lon_mesh, lat_mesh, Hcoef_cyclic, levels=levels,
                  cmap=cmap, transform=ccrs.PlateCarree(), extend='both')

# Add map features
ax3.coastlines(linewidth=0.5)
ax3.add_feature(cfeature.BORDERS, linewidth=0.5)
ax3.add_feature(cfeature.LAKES, linewidth=0.5, facecolor='none')

# Add grid lines
# Add grid lines
gl = ax3.gridlines(draw_labels=True, linewidth=0.25, color='black', alpha=1.0)
gl.top_labels = False
gl.right_labels = False
gl.xlocator = mticker.FixedLocator(np.arange(-180, 181, 10))  # Every 10 deg lon
gl.ylocator = mticker.FixedLocator(np.arange(-90, 91, 10))    # Every 10 deg lat
gl.xlabel_style = {'size': 8, 'rotation': 0}  # Adjust size as needed
gl.ylabel_style = {'size': 8}

# Set labels only every 30 degrees
gl.xformatter = mticker.FuncFormatter(lambda x, _: f"{int(x)}" if x % 30 == 0 else "")
gl.yformatter = mticker.FuncFormatter(lambda y, _: f"{int(y)}" if y % 30 == 0 else "")

ax3.set_extent([-180, 180, -90, 90], crs=ccrs.PlateCarree())


# Plot 4: Vertical Window Profile
ax4 = plt.subplot(gs[2, :])

ax4.set_box_aspect(0.8)  # Make panel taller relative to width (X height-to-width ratio)

if use_nclev and nclev is not None:
    ax4.plot(Vcoef, nclev, 'ro-', markersize=4)
    ax4.set_ylim(nclev[0], nclev[-1])
    ax4.set_ylabel("Approx model plev")
    ax4.set_yticks(np.arange(0, nclev[-1] + 100, 100))
else:
    ax4.plot(Vcoef, lev, 'ro-', markersize=4)
    ax4.set_ylim(lev[0], lev[-1])
    ax4.set_ylabel("Model Level Index")

ax4.set_title("Vertical Window Profile")
ax4.set_xlabel("Nudging strength")
ax4.set_xlim(-0.05, 1.05)
ax4.invert_yaxis()

# Grid lines
ax4.minorticks_on()
ax4.grid(which='major', linestyle='-', linewidth=0.5, alpha=1.0, color='gray')
ax4.grid(which='minor', linestyle=':', linewidth=0.5, alpha=1.0, color='gray')

# Add colorbar
plt.colorbar(cs, ax=[ax3], orientation='vertical',
             ticks=np.linspace(0, 1.0, 6),
             label="Nudging strength")

# Save figure
plt.savefig("Wcoef.png", dpi=300, bbox_inches='tight')
