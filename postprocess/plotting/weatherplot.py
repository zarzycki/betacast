#========================================================
# This script regrids a single CAM history file to a lat/lon
# grid and plots a set of forecast maps for each region
# Colin Zarzycki
# Python port of weatherplot.ncl
#========================================================
# python weatherplot.py --inisec 43200 --iniday 25 --inimon 05 --iniyear 2018 \
#   --filename /glade/scratch/zarzycki/forecast_natlantic_30_x4_CAM5_L30/run/2018052512/_forecast_natlantic_30_x4_CAM5_L30.cam.h0.2018-05-27-00000.nc \
#   --wgt_file /glade/p/work/zarzycki/maps/forecast_plot_maps/map_natlantic_30_x4_to_0.25x0.25glob_bilinear.nc
import os
import sys
import time
import argparse
from datetime import datetime

import numpy as np
import xarray as xr
import cftime
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap, BoundaryNorm, to_rgba
from scipy.ndimage import minimum_filter
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import cartopy.io.shapereader as shpreader
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "py_functions"))
from horizremap import remap_with_weights_wrapper, linint2
from pyfuncs import poisson_grid_fill

PATH_TO_NCL = os.environ.get("PATH_TO_NCL", os.path.dirname(os.path.abspath(__file__)))

def rainratetodbz(rainrate, units):
    #========================================================
    # This function takes in rainfall rate (here in mm/hr)
    # and converts it to DbZ for plotting as simulated
    # radar reflectivity
    #
    # rainrate = 2D (or 3D in time) array of integrated
    # rainfall rate
    # units, 0 = mm/hr, 1 = m/s, 2 cm/s
    #
    # Colin Zarzycki, zarzycki@umich.edu, 11/14/2012
    #========================================================

    # If units of m/s, convert to mm/hr
    if units == 1:
        rainrate = rainrate*3600000
    elif units == 2:
        rainrate = rainrate*36000

    eps = 10**-25        # very small number to prevent log from = 0
    # Uses common NWS NEXRAD published conversion formula
    #dbz = 10*np.log10(eps + (200 * rainrate**(1.6) ) )

    # Z = 250 R^1.2
    # any NaN that crop up stay NaN (missing)
    with np.errstate(invalid="ignore"):
        dbz = 10*np.log10(eps + (250*rainrate**(1.2)))

    return dbz

def wallClockElapseTime(wcStrt, label):
    print(f"{label}: {time.time()-wcStrt:.1f} s")

def read_colormap_file(filename):
    # NCL .rgb colormap, skip header/comment lines, returns (ncolors, 3) RGB in 0-1
    rgb = []
    with open(filename) as f:
        for line in f:
            line = line.split("#")[0].strip()
            if line == "" or line.startswith(";") or line.lower().startswith("ncolors"):
                continue
            rgb.append([float(x) for x in line.split()[:3]])
    rgb = np.array(rgb)
    if rgb.max() > 1.0:
        rgb = rgb/255.
    return rgb

def manual_levels(minval, maxval, spacing):
    # cnLevelSelectionMode = "ManualLevels"
    return np.arange(minval, maxval+spacing/2., spacing)

def ncl_round(x):
    # NCL toint truncates toward zero
    return int(np.trunc(x))

def uv2vr_cfd(u, v, lat, lon):
    # relative vorticity via centered finite differences, one-sided at the edges (NCL bound_opt=2)
    re = 6.37122e6
    rad = np.pi/180.
    dlon = np.gradient(lon*rad)
    dlat = np.gradient(lat*rad)
    coslat = np.cos(lat*rad)[:, np.newaxis]
    dvdx = np.gradient(v, axis=-1)/(re*coslat*dlon[np.newaxis, :])
    dudy = np.gradient(u, axis=-2)/(re*dlat[:, np.newaxis])
    return dvdx - dudy + (u/re)*np.tan(lat*rad)[:, np.newaxis]

def area_hi2lores(xi, yi, fi, xo, yo):
    # area (overlap) weighted average of fi onto coarser grid xo/yo, missing values ignored
    def overlap_wgts(ci, co):
        ei = np.concatenate(([ci[0]-(ci[1]-ci[0])/2.], (ci[1:]+ci[:-1])/2., [ci[-1]+(ci[-1]-ci[-2])/2.]))
        eo = np.concatenate(([co[0]-(co[1]-co[0])/2.], (co[1:]+co[:-1])/2., [co[-1]+(co[-1]-co[-2])/2.]))
        lo = np.maximum(eo[:-1, np.newaxis], ei[np.newaxis, :-1])
        hi = np.minimum(eo[1:, np.newaxis], ei[np.newaxis, 1:])
        return np.maximum(hi-lo, 0.)
    wx = overlap_wgts(xi, xo)
    wy = overlap_wgts(yi, yo)
    valid = ~np.isnan(fi)
    num = wy @ np.where(valid, fi, 0.) @ wx.T
    den = wy @ valid.astype(float) @ wx.T
    with np.errstate(invalid="ignore", divide="ignore"):
        return np.where(den > 0, num/den, np.nan)

def fill_cmap(cmap, levels, opacity=1.0):
    # cnFillPalette: NCL spans the palette, picking nlevels+1 evenly spaced colors
    if isinstance(cmap, np.ndarray):
        ncolors = cmap.shape[0]
        ix = np.round(np.linspace(0, ncolors-1, len(levels)+1)).astype(int)
        colors = [tuple(cmap[ii]) + (opacity,) for ii in ix]
    else:
        colors = [to_rgba(c, opacity) for c in cmap]
    cm = ListedColormap(colors[1:-1])
    cm.set_under(colors[0])
    cm.set_over(colors[-1])
    norm = BoundaryNorm(levels, cm.N)
    return cm, norm

def map_ticks(minval, maxval):
    for spacing in [1., 2., 5., 10., 20., 30., 60., 90.]:
        if (maxval-minval)/spacing <= 4.:
            break
    return np.arange(np.ceil(minval/spacing)*spacing, maxval+0.001, spacing)

# US states only (mpOutlineBoundarySets = "geophysicalandusstates")
def us_states(resolution):
    shp = shpreader.Reader(shpreader.natural_earth(resolution=resolution, category="cultural", name="admin_1_states_provinces_lakes"))
    return cfeature.ShapelyFeature([rec.geometry for rec in shp.records() if rec.attributes["adm0_a3"] == "USA"], ccrs.PlateCarree())

def raster_fill(ax, x, y, data, cmap, norm, clon):
    # cnFillMode = "RasterFill" + cnRasterSmoothingOn, data are bilinearly interpolated to the image pixels
    # extent given in the (shifted) map projection coordinates so cartopy does not warp the image
    dx = x[1]-x[0]
    dy = y[1]-y[0]
    extent = [x[0]-dx/2.-clon, x[-1]+dx/2.-clon, y[0]-dy/2., y[-1]+dy/2.]
    return ax.imshow(np.ma.masked_invalid(data), cmap=cmap, norm=norm, origin="lower", extent=extent,
                     interpolation="bilinear", transform=ax.projection)

def gsn_csm_contour_map(data, res):
    # filled (raster) contour map with labelbar, returns fig/ax for overlays
    fig = plt.figure(figsize=(8.25, 8.25))
    proj = ccrs.PlateCarree(central_longitude=res["mpCenterLonF"])
    ax = fig.add_axes([0.05, 0.25, 0.9, 0.6], projection=proj)
    ax.set_extent([res["mpMinLonF"], res["mpMaxLonF"], res["mpMinLatF"], res["mpMaxLatF"]], crs=ccrs.PlateCarree())

    if res.get("mpFillOn", False):
        ax.add_feature(cfeature.OCEAN.with_scale(res["mpDataBaseVersion"]), facecolor=res["mpOceanFillColor"], zorder=0)
        ax.add_feature(cfeature.LAND.with_scale(res["mpDataBaseVersion"]), facecolor=res["mpLandFillColor"], zorder=0)
        ax.add_feature(cfeature.LAKES.with_scale(res["mpDataBaseVersion"]), facecolor=res["mpInlandWaterFillColor"], zorder=0)

    # only draw points inside (or just outside) the plotted region
    iy = np.nonzero((lat >= res["mpMinLatF"]-1.) & (lat <= res["mpMaxLatF"]+1.))[0]
    ix = np.nonzero((lon >= res["mpMinLonF"]-1.) & (lon <= res["mpMaxLonF"]+1.))[0]
    cmap, norm = fill_cmap(res["cnFillPalette"], res["cnLevels"], res.get("cnFillOpacityF", 1.0))
    cf = raster_fill(ax, lon[ix], lat[iy], data[np.ix_(iy, ix)], cmap, norm, res["mpCenterLonF"])
    ax.set_extent([res["mpMinLonF"], res["mpMaxLonF"], res["mpMinLatF"], res["mpMaxLatF"]], crs=ccrs.PlateCarree())

    ax.coastlines(resolution=res["mpDataBaseVersion"], color=res["mpGeophysicalLineColor"], linewidth=0.35*res["mpGeophysicalLineThicknessF"])
    ax.add_feature(cfeature.LAKES.with_scale(res["mpDataBaseVersion"]), facecolor="none", edgecolor=res["mpGeophysicalLineColor"], linewidth=0.35*res["mpGeophysicalLineThicknessF"])
    ax.add_feature(us_states(res["mpDataBaseVersion"]), facecolor="none", edgecolor="black", linewidth=0.3*res["mpUSStateLineThicknessF"])

    ax.set_xticks(map_ticks(res["mpMinLonF"], res["mpMaxLonF"]), crs=ccrs.PlateCarree())
    ax.set_yticks(map_ticks(res["mpMinLatF"], res["mpMaxLatF"]), crs=ccrs.PlateCarree())
    ax.xaxis.set_major_formatter(LongitudeFormatter(degree_symbol=""))
    ax.yaxis.set_major_formatter(LatitudeFormatter(degree_symbol=""))
    ax.tick_params(labelsize=9, direction="out", length=3, top=True, right=True)

    ax.set_title(res["gsnLeftString"], loc="left", fontsize=11)
    ax.set_title(res["gsnRightString"], loc="right", fontsize=11)

    pos = ax.get_position()
    cax = fig.add_axes([0.5-res["pmLabelBarWidthF"]/2., pos.y0-0.065, res["pmLabelBarWidthF"], 0.025])
    cb = fig.colorbar(cf, cax=cax, orientation="horizontal", extend="both", extendrect=True, extendfrac="auto", drawedges=True)
    levels = res["cnLevels"]
    # NCL thins labelbar labels when they would overlap
    stride = res.get("lbLabelStride", max(1, int(np.ceil(len(levels)*max(len(f"{x:g}") for x in levels)/100.))))
    cb.set_ticks(levels[::stride])
    cb.set_ticklabels([f"{x:g}" for x in levels[::stride]])
    cb.ax.tick_params(labelsize=10, length=0, labelrotation=res.get("lbLabelAngleF", 0))

    return fig, ax

def gsn_csm_contour(ax, data, res, x=None, y=None):
    # line contours overlaid on an existing map
    if x is None:
        x, y = lon, lat
    cs = ax.contour(x, y, data, levels=res["cnLevels"], colors=res.get("cnLineColor", "black"),
                    linewidths=0.5*res.get("cnLineThicknessF", 1.0), transform=ccrs.PlateCarree())
    if res.get("cnLineLabelsOn", True):
        ax.clabel(cs, fmt="%g", fontsize=11, inline=True)

    if res.get("cnLowLabelsOn", False):
        # mark local minima inside the plotted region with L and value
        win = max(3, int((res["mpMaxLonF"]-res["mpMinLonF"])/12./dlon))
        isLow = (data == minimum_filter(data, size=win, mode="nearest")) & ~np.isnan(data)
        jj, ii = np.nonzero(isLow)
        for j, i in zip(jj, ii):
            if res["mpMinLatF"] < y[j] < res["mpMaxLatF"] and res["mpMinLonF"] < x[i] < res["mpMaxLonF"]:
                ax.text(x[i], y[j], "L", color=res["cnLowLabelFontColor"], fontsize=14, ha="center", va="center", transform=ccrs.PlateCarree())
                ax.text(x[i]+0.01*(res["mpMaxLonF"]-res["mpMinLonF"]), y[j], f"{ncl_round(data[j,i])}", color=res["cnLowLabelFontColor"],
                        fontsize=8, ha="left", va="top", transform=ccrs.PlateCarree())
    return cs

def gsn_add_text_grid(ax, data, stride, minval, txres):
    # write integer values of data on the plot every stride points
    for ii in range(0, nlat, stride):
        for jj in range(0, nlon, stride):
            if lat[ii] >= minlat and lat[ii] <= maxlat and lon[jj] >= minlon and lon[jj] <= maxlon:
                if not np.isnan(data[ii,jj]) and data[ii,jj] >= minval:
                    ax.text(lon[jj], lat[ii], str(ncl_round(data[ii,jj])), fontsize=txres["fontsize"], fontweight="bold",
                            ha="center", va="center", transform=ccrs.PlateCarree(), clip_on=True)

def ESMF_regrid_with_weights(var, wgt_file):
    # output keeps the input type (float), like NCL
    data, dstlat, dstlon = remap_with_weights_wrapper(var, wgt_file)
    return data.astype(var.dtype), dstlat, dstlon

def frame(fig, fname):
    fig.savefig(fname+"."+out_type, dpi=100)
    plt.close(fig)

parser = argparse.ArgumentParser(description="Plot forecast maps from a CAM history file")
parser.add_argument("--iniyear", type=int, required=True, help="Init year (1999)")
parser.add_argument("--inimon", type=int, required=True, help="Init month (12)")
parser.add_argument("--iniday", type=int, required=True, help="Init day (25)")
parser.add_argument("--inisec", type=int, required=True, help="Init seconds (24000)")
parser.add_argument("--filename", type=str, required=True, help="CAM history file")
parser.add_argument("--wgt_file", type=str, required=True, help="ESMF weight file to lat/lon grid")
args = parser.parse_args()

wcStrt = time.time()

iniyear = args.iniyear
inimon = args.inimon
iniday = args.iniday
inisec = args.inisec
filename = args.filename
wgt_file = args.wgt_file

iniyearstr = f"{iniyear:04d}"
inimonstr = f"{inimon:02d}"
inidaystr = f"{iniday:02d}"
inisecstr = f"{inisec:05d}"
inihourstr = f"{inisec//3600:02d}"
inihoursince = datetime(iniyear, inimon, iniday, inisec//3600)

do_interp = True
out_type = "png"

datafile = xr.open_dataset(filename, decode_times=False)

strs = os.path.basename(filename).split(".")
# strs(0) = case
# strs(1) = model
# strs(2) = file (h1, etc)
# strs(3) = date
# strs(4) = extension

datestrs = strs[-2].split("-")
yearstr = datestrs[0]
monstr = datestrs[1]
daystr = datestrs[2]
secstr = datestrs[3]
hourstr = f"{int(secstr)//3600:02d}"

hoursince = datetime(int(yearstr), int(monstr), int(daystr), int(hourstr))
forecastHour = int((hoursince - inihoursince).total_seconds()/3600.)
forecastHourStr = f"{forecastHour:03d}"
print(forecastHourStr)

if forecastHour % 6 != 0:
    print("For now, we are only printing data at 6 hr increments")
    sys.exit(0)

print("Forecast for: "+hourstr+"Z "+monstr+" "+daystr+" "+yearstr)
print("Initialized at: "+inihourstr+"Z "+inimonstr+" "+inidaystr+" "+iniyearstr)

#---Get variables to regrid and the source lat/lon grid
prect_se = datafile["PRECT"]
u500_se = datafile["U500"]
v500_se = datafile["V500"]
z500_se = datafile["Z500"]
tmq_se = datafile["TMQ"]
flut_se = datafile["FLUT"]
u850_se = datafile["U850"]
v850_se = datafile["V850"]
mslp_se = datafile["PSL"]/100
preclav_se = datafile["PRECLav"]
preccav_se = datafile["PRECCav"]
precbsn_se = datafile["PRECBSN"]
ts_se = datafile["TREFHT"]
ptype_se = datafile["PTYPE"]

prectav_se = preclav_se+preccav_se
del preclav_se, preccav_se

print(float(mslp_se.max()))
print(float(mslp_se.min()))

wallClockElapseTime(wcStrt, "Loading and initial processing")

print("Begin interp")

if do_interp:
    prect, lat, lon = ESMF_regrid_with_weights(prect_se, wgt_file)
    u500, _, _ = ESMF_regrid_with_weights(u500_se, wgt_file)
    v500, _, _ = ESMF_regrid_with_weights(v500_se, wgt_file)
    z500, _, _ = ESMF_regrid_with_weights(z500_se, wgt_file)
    tmq, _, _ = ESMF_regrid_with_weights(tmq_se, wgt_file)
    flut, _, _ = ESMF_regrid_with_weights(flut_se, wgt_file)
    u850, _, _ = ESMF_regrid_with_weights(u850_se, wgt_file)
    v850, _, _ = ESMF_regrid_with_weights(v850_se, wgt_file)
    mslp, _, _ = ESMF_regrid_with_weights(mslp_se, wgt_file)
    prectav, _, _ = ESMF_regrid_with_weights(prectav_se, wgt_file)
    precbsn, _, _ = ESMF_regrid_with_weights(precbsn_se, wgt_file)
    ts, _, _ = ESMF_regrid_with_weights(ts_se, wgt_file)
    ptype, _, _ = ESMF_regrid_with_weights(ptype_se, wgt_file)

    guess     = 1                # use zonal means
    is_cyclic = True             # cyclic [global]
    nscan     = 500              # usually much less than this
    eps       = 1.e-2            # variable dependent
    relc      = 0.5              # relaxation coefficient
    opt       = 0                # not used

    poisson_grid_fill( u500, is_cyclic, guess, nscan, eps, relc, opt)
    poisson_grid_fill( v500, is_cyclic, guess, nscan, eps, relc, opt)
    poisson_grid_fill( prect, is_cyclic, guess, nscan, eps, relc, opt)
    poisson_grid_fill( z500, is_cyclic, guess, nscan, eps, relc, opt)
    poisson_grid_fill( tmq, is_cyclic, guess, nscan, eps, relc, opt)
    poisson_grid_fill( flut, is_cyclic, guess, nscan, eps, relc, opt)
    poisson_grid_fill( u850, is_cyclic, guess, nscan, eps, relc, opt)
    poisson_grid_fill( v850, is_cyclic, guess, nscan, eps, relc, opt)
    poisson_grid_fill( mslp, is_cyclic, guess, nscan, eps, relc, opt)
    poisson_grid_fill( prectav, is_cyclic, guess, nscan, eps, relc, opt)
    poisson_grid_fill( precbsn, is_cyclic, guess, nscan, eps, relc, opt)
    poisson_grid_fill( ts, is_cyclic, guess, nscan, eps, relc, opt)
    poisson_grid_fill( ptype, is_cyclic, guess, nscan, eps, relc, opt)

    print("end interp")
else:
    # NOTE: plotting below assumes data on a regular lat/lon grid
    prect = prect_se.values
    u500 = u500_se.values
    v500 = v500_se.values
    z500 = z500_se.values
    tmq = tmq_se.values
    flut = flut_se.values
    u850 = u850_se.values
    v850 = v850_se.values
    mslp = mslp_se.values
    prectav = prectav_se.values
    precbsn = precbsn_se.values
    ts = ts_se.values
    ptype = ptype_se.values
    lat = datafile["lat"].values
    lon = datafile["lon"].values

datafile.close()
del prect_se, u500_se, v500_se, z500_se, tmq_se, flut_se, u850_se, v850_se, mslp_se, prectav_se, precbsn_se, ts_se, ptype_se

wallClockElapseTime(wcStrt, "Interpolating")

ntime, nlat, nlon = mslp.shape

lat2d, lon2d = np.meshgrid(lat, lon, indexing="ij")

dlat = lat[1]-lat[0]
dlon = lon[1]-lon[0]

print(f"dlat: {dlat:g} dlon: {dlon:g}")

#RAINDBZ = rainratetodbz(prect,1)

prect_mmday = prect*3600000          # mm/hr

prectav_in = prectav*10800           # seconds in 3 hours, now gives us meters per timestep
prectav_in = prectav_in*39.3701      # now inches per timestep
prectav_in = np.where(prectav_in < 0.1, np.nan, prectav_in)

precbsn_in = precbsn*10800           # seconds in 3 hours, now gives us meters per timestep
precbsn_in = precbsn_in*39.3701      # now inches per timestep
precbsn_in = precbsn_in*10           # now snow water ratio
precbsn_in = np.where(precbsn_in < 0.1, np.nan, precbsn_in)

vort500 = np.zeros_like(u500)
for tt in range(ntime):
    vort500[tt,:,:] = uv2vr_cfd(u500[tt,:,:], v500[tt,:,:], lat, lon)
vort500 = vort500/1e-5               # 1E-5 1/s

wind = np.sqrt(u850**2 + v850**2)

tsf = (ts-273.15)*9.0/5.0+32.0      # degrees F

ptype = np.trunc(ptype)
ptype = np.where(prect_mmday < 0.1, np.nan, ptype)

#*****************************************************
# ANOMS
#*****************************************************

print("doing anoms")
f_anom = xr.open_dataset(os.path.join(PATH_TO_NCL, "climo-files", "TS-1980-2009-NCEP.nc"), decode_times=False)
ANOMTIME = f_anom["time"]
ANOMDATE = cftime.date2num(cftime.datetime(2000, int(monstr), int(daystr), int(hourstr), calendar="standard"), ANOMTIME.attrs["units"])
TSANOM = f_anom["TS"].sel(time=ANOMDATE, method="nearest")
TSANOM = TSANOM.assign_coords(lon=((TSANOM.lon+180) % 360)-180).sortby("lon").sortby("lat")   # lonFlip
TSANOM_CAM = linint2(TSANOM.lon.values, TSANOM.lat.values, TSANOM.values, True, lon, lat)
TSANOM_CAM_F = (TSANOM_CAM-273.15)*9.0/5.0+32.0
tsf_anom = tsf.copy()
tsf_anom[0,:,:] = tsf[0,:,:]-TSANOM_CAM_F[:,:]
f_anom.close()
del f_anom, ANOMTIME, ANOMDATE, TSANOM, TSANOM_CAM, TSANOM_CAM_F

wallClockElapseTime(wcStrt, "Var. preparation")

#*****************************************************
# PLOT
#*****************************************************

for plotindex in range(3):

    snowstride = int(360.0/dlat)
    tsstride = int(360.0/dlat)
    windstride = int(360.0/dlat)

    # NORTH ATLANTIC
    if plotindex == 0:
        plotname = "natl"
        maxlat = 55.           # NATL
        minlat = 5.
        maxlon = -5.
        minlon = -100.
        snowstride = int(4.5/dlat)
        tsstride = int(6.0/dlat)
        windstride = int(4.0/dlat)
    #if plotindex == 0:
    #    plotname = "natl"
    #    maxlat = 53.           # CONUS
    #    minlat = 21.
    #    maxlon = 295.-360.
    #    minlon = 231.-360.
    #    snowstride = int(2.25/dlat)
    #    tsstride = int(3.125/dlat)
    #    windstride = int(2.0/dlat)
    #elif plotindex == 1:    # ORIG EPAC
    #    plotname = "epac"
    #    maxlat = 40.
    #    minlat = 5.
    #    maxlon = -90.
    #    minlon = -165.
    #    snowstride = int(4.5/dlat)
    #    tsstride = int(6.0/dlat)
    #    windstride = int(4.0/dlat)
    elif plotindex == 1:     # NE USA
        plotname = "epac"
        maxlat = 49.
        minlat = 37.
        maxlon = 295.-360.
        minlon = 275.-360.
        snowstride = int(0.75/dlat)
        tsstride = int(1.0/dlat)
        windstride = int(2.0/dlat)
    #elif plotindex == 1:     # Modif EPAC for Caribbean
    #    plotname = "epac"
    #    maxlat = 30.
    #    minlat = 10.
    #    maxlon = -55.
    #    minlon = -85.
    #elif plotindex == 1:     # Modif EPAC for Gulf
    #    plotname = "epac"
    #    maxlat = 35.
    #    minlat = 10.
    #    maxlon = -70.
    #    minlon = -105.
    #elif plotindex == 1:     # Modif EPAC for SE USA
    #    plotname = "epac"
    #    maxlat = 37.
    #    minlat = 20.
    #    maxlon = 289.
    #    minlon = 269.
    #elif plotindex == 2:     # for Gulf + Carrib
    #    plotname = "float1"
    #    maxlat = 33.
    #    minlat = 8.
    #    maxlon = -58.
    #    minlon = -100.
    #elif plotindex == 1:     # Modif EPAC for CONUS
    #    plotname = "epac"
    #    maxlat = 53.
    #    minlat = 25.
    #    maxlon = 300.
    #    minlon = 230.
    #elif plotindex == 2:     # Modif EPAC for Colorado
    #    plotname = "float1"
    #    maxlat = 41.5
    #    minlat = 36.5
    #    maxlon = 258.5-360.
    #    minlon = 250.5-360.
    #elif plotindex == 2:     # Modif EPAC for WPAC
    #    plotname = "float1"
    #    maxlat = 40.0
    #    minlat = 5.0
    #    maxlon = 170.0
    #    minlon = 75.0
    #    snowstride = int(4.5/dlat)
    #    tsstride = int(6.0/dlat)
    #    windstride = int(4.0/dlat)
    elif plotindex == 2:     # Modif GLOB
        plotname = "float1"
        maxlat = 80.0
        minlat = -80.0
        maxlon = 180.
        minlon = -180.
        snowstride = int(7.5/dlat)
        tsstride = int(15.0/dlat)
        windstride = int(5.0/dlat)
    else:
        plotname = "glob"
        maxlat = 90.           # choose subregion
        minlat = -90.
        maxlon = 180.
        minlon = -180.

    print("********* DOING: "+plotname)
    print("Using strides...")
    print("snowstride: "+str(snowstride))
    print("tsstride: "+str(tsstride))
    print("windstride: "+str(windstride))
    print("------------------------------")
    print("Plotting: "+plotname)

    #----------- Global resources -----------------------------------------

    res = {}
    res["gsnLeftString"] = ""
    res["gsnRightString"] = ""

    res["mpFillOn"] = False                       # turn off gray fill
    res["mpGeophysicalLineColor"] = "Black"       # color of cont. outlines
    res["mpGeophysicalLineThicknessF"] = 3        # thickness of outlines

    res["mpMaxLatF"] = maxlat                     # choose subregion
    res["mpMinLatF"] = minlat
    res["mpMinLonF"] = minlon
    res["mpMaxLonF"] = maxlon
    res["mpCenterLonF"] = (res["mpMinLonF"]+res["mpMaxLonF"])/2.

    # this adds county level data (no county database in cartopy, use hi-res states)
    if plotindex == 9:
        res["mpDataBaseVersion"] = "10m"
        res["mpUSStateLineThicknessF"] = 4.0
    else:
        res["mpDataBaseVersion"] = "50m"
        res["mpUSStateLineThicknessF"] = 1.0

    res["gsnLeftString"] = "VALID "+hourstr+"Z "+monstr+" "+daystr+" "+yearstr+" [+"+forecastHourStr+"]"

    res["pmLabelBarWidthF"] = 0.7

    res_contour = {}
    res_contour["cnLineLabelsOn"] = True
    res_contour["mpMinLatF"] = minlat
    res_contour["mpMaxLatF"] = maxlat
    res_contour["mpMinLonF"] = minlon
    res_contour["mpMaxLonF"] = maxlon

    res_mslp = dict(res_contour)
    if plotindex == 0 or plotindex == 1:
        res_mslp["cnLowLabelsOn"] = True
    else:
        res_mslp["cnLowLabelsOn"] = False
    res_mslp["cnLowLabelFontColor"] = "red"
    res_mslp["cnLineLabelsOn"] = False
    res_mslp["cnLevels"] = manual_levels(900.0, 1060.0, 4.0)

    res_z500 = dict(res_contour)
    res_z500["cnLevels"] = manual_levels(4800.0, 6200.0, 60.0)

    # =============================================================================================================================================

    print("Plotting instant precip...")

    #;;; RAIN, SIMULATED RADAR
    # res_rain = dict(res)
    # res_rain["cnFillPalette"] = read_colormap_file(os.path.join(PATH_TO_NCL, "colormaps", "rainsnowradar.rgb"))
    # res_rain["gsnRightString"] = "Simulated Reflectivity (dBZ)"
    # res_rain["cnLevels"] = np.array([0,5,10,15,20,25,30,35,40,45,50,55,60,65,70])
    # res_rain["lbLabelStride"] = 2
    # fig, ax = gsn_csm_contour_map(RAINDBZ[0,:,:], res_rain)
    # frame(fig, "dbz_"+plotname+"_"+yearstr+monstr+daystr+hourstr)

    #;;; RAIN, MM/DAY
    res_rain = dict(res)
    cmap = read_colormap_file(os.path.join(PATH_TO_NCL, "colormaps", "prcp_1.rgb"))
    ncolors, npairs = cmap.shape
    nnewcolors = 200
    cmap_interp = np.zeros((nnewcolors, npairs))
    for ii in range(npairs):
        cmap_interp[:,ii] = np.interp(np.linspace(0, 1, nnewcolors), np.linspace(0, 1, ncolors), cmap[:,ii])
    cmap_interp[0,:] = cmap[0,:]
    cmap_interp[1,:] = cmap[1,:]

    res_rain["cnFillPalette"] = cmap

    res_rain["cnLevels"] = np.array([0.1,0.25,0.5,1.0,2.0,3.0,4.0,5.0,7.5,10.0,15.0,20.0,25.0,35.0])   # set levels

    res_rain["gsnRightString"] = "Instant. precip. rate (color, mm/hr), MSLP (black contour, hPa)"

    fig, ax = gsn_csm_contour_map(prect_mmday[0,:,:], res_rain)

    # overlay
    gsn_csm_contour(ax, mslp[0,:,:], res_mslp)

    frame(fig, "prect_"+plotname+"_"+yearstr+monstr+daystr+hourstr)

    del cmap, cmap_interp

    wallClockElapseTime(wcStrt, "Instant precip "+plotname)

    # =============================================================================================================================================

    print("Plotting ptype...")

    #;;; ptype, MM/DAY
    res_rainptype = dict(res)
    cmap = read_colormap_file(os.path.join(PATH_TO_NCL, "colormaps", "MPL_Greys.rgb"))
    res_rainptype["cnFillPalette"] = cmap

    res_rainptype["cnLevels"] = np.array([0.1,0.15,0.25,0.5,0.75,1.0,2.0,3.0,4.0,5.0,6.0,8.0,10.0,15.0])   # set levels

    res_rainptype["gsnRightString"] = "Instant. precip. rate (shading, mm/hr), precip. type (color)"

    fig, ax = gsn_csm_contour_map(prect_mmday[0,:,:], res_rainptype)

    cmap = ["white","dodgerblue","orange","green","red","white"]
    ptype_levels = np.array([0,1,2,3,4])
    ptype_cmap, ptype_norm = fill_cmap(cmap, ptype_levels, 0.2)
    iy = np.nonzero((lat >= minlat-1.) & (lat <= maxlat+1.))[0]
    ix = np.nonzero((lon >= minlon-1.) & (lon <= maxlon+1.))[0]
    raster_fill(ax, lon[ix], lat[iy], ptype[0][np.ix_(iy, ix)], ptype_cmap, ptype_norm, res["mpCenterLonF"])
    ax.set_extent([minlon, maxlon, minlat, maxlat], crs=ccrs.PlateCarree())

    # overlay
    gsn_csm_contour(ax, mslp[0,:,:], res_mslp)

    frame(fig, "ptype_"+plotname+"_"+yearstr+monstr+daystr+hourstr)

    del cmap

    wallClockElapseTime(wcStrt, "ptype "+plotname)

    # =============================================================================================================================================

    print("Plotting accumulated precip...")
    res_sumprecip = dict(res)

    cmap = read_colormap_file(os.path.join(PATH_TO_NCL, "colormaps", "prcp_1.rgb"))
    ncolors, npairs = cmap.shape
    nnewcolors = 200
    cmap_interp = np.zeros((nnewcolors, npairs))
    for ii in range(npairs):
        cmap_interp[:,ii] = np.interp(np.linspace(0, 1, nnewcolors), np.linspace(0, 1, ncolors), cmap[:,ii])
    cmap_interp[0,:] = cmap[0,:]
    cmap_interp[1,:] = cmap[1,:]

    res_sumprecip["cnFillPalette"] = cmap_interp

    res_sumprecip["cnLevels"] = np.array([0.05,0.1,0.15,0.2,0.3,0.4,0.5,0.75,1,1.25,1.5,1.75,2,2.5,3,3.5,4,4.5,5,5.5,6,7,8,9,10,12,14,16,20,24])   # set levels
    res_sumprecip["lbLabelStride"] = 1
    res_sumprecip["lbLabelAngleF"] = 50                 # angle labels

    res_sumprecip["gsnRightString"] = "Accum. precipitation (in.)"

    fig, ax = gsn_csm_contour_map(prectav_in[0,:,:], res_sumprecip)
    frame(fig, "sumprect_"+plotname+"_"+yearstr+monstr+daystr+hourstr)

    del cmap, cmap_interp

    wallClockElapseTime(wcStrt, "accumulated precip "+plotname)

    # =============================================================================================================================================

    print("Plotting accumulated snow...")

    res_sumsnow = dict(res)

    cmap = read_colormap_file(os.path.join(PATH_TO_NCL, "colormaps", "wxbell_snow.rgb"))
    res_sumsnow["cnFillPalette"] = cmap[::-1,:]       # reverse color map
    res_sumsnow["cnFillOpacityF"] = 0.75

    res_sumsnow["cnLevels"] = np.array([ 0.25, 0.5, 1., 1.5, 2., 3., 4., 5., 6., 7., 8., 9., 10., 11., 12., 13., 14., 15., 16., 17., 18., 19., 20., 21., 22., 23., 24., 26., 28., 30., 32., 34., 36., 38., 40., 42., 44., 46., 48., 52., 56., 60., 64., 68., 72.])

    res_sumsnow["gsnRightString"] = "Accum. Bourgouin snow (10:1 in.)"

    fig, ax = gsn_csm_contour_map(precbsn_in[0,:,:], res_sumsnow)

    txres = {}
    txres["fontsize"] = 8             # Set the font height

    stride = snowstride

    gsn_add_text_grid(ax, precbsn_in[0,:,:], stride, 1.0, txres)

    frame(fig, "sumsnow_"+plotname+"_"+yearstr+monstr+daystr+hourstr)

    del cmap, stride, txres

    wallClockElapseTime(wcStrt, "accumulated snowfall "+plotname)

    # =============================================================================================================================================

    print("Plotting 500 mb vort + z500 + U/V 500...")

    res_500vort = dict(res)

    cmap = read_colormap_file(os.path.join(PATH_TO_NCL, "colormaps", "WhiteYellowOrangeRed.rgb"))
    res_500vort["cnFillPalette"] = cmap

    res_500vort["cnLevels"] = manual_levels(0.0, 40.0, 2.0)

    res_500vort["gsnRightString"] = "500mb vort/Z/wind"

    fig, ax = gsn_csm_contour_map(vort500[0,:,:], res_500vort)
    gsn_csm_contour(ax, z500[0,:,:], res_z500)

    vecres = {}                        # vector only resources
    vecres["vcRefMagnitudeF"] = 20     # define vector ref mag
    vecres["vcRefLengthF"] = 0.045     # define length of vec ref (fraction of figure width)

    # this is 8 for 0.25deg data, 16 for 0.125deg data
    stride = windstride

    # NCL used CurlyVector, matplotlib only has straight arrows
    Q = ax.quiver(lon[::stride], lat[::stride], u500[0,::stride,::stride], v500[0,::stride,::stride],
                  scale=vecres["vcRefMagnitudeF"]/vecres["vcRefLengthF"], scale_units="width", width=0.001,
                  headwidth=5, headlength=6, color="#333333", transform=ccrs.PlateCarree(), regrid_shape=None)
    ax.quiverkey(Q, 0.93, -0.16, vecres["vcRefMagnitudeF"], str(vecres["vcRefMagnitudeF"])+"\nReference Vector",
                 labelpos="S", fontproperties={"size": 9})

    frame(fig, "500vort_"+plotname+"_"+yearstr+monstr+daystr+hourstr)

    del cmap, stride

    wallClockElapseTime(wcStrt, "500 mb fields "+plotname)

    # ====== Plot TMQ ============

    print("Plotting TMQ...")

    res_tmq = dict(res)

    cmap = read_colormap_file(os.path.join(PATH_TO_NCL, "colormaps", "WhiteBlueGreenYellowRed.rgb"))
    res_tmq["cnFillPalette"] = cmap

    res_tmq["cnLevels"] = manual_levels(20.0, 66.0, 2.0)

    res_tmq["gsnRightString"] = "Total precipitable water (mm)"

    fig, ax = gsn_csm_contour_map(tmq[0,:,:], res_tmq)
    frame(fig, "tmq_"+plotname+"_"+yearstr+monstr+daystr+hourstr)

    del cmap

    wallClockElapseTime(wcStrt, "TMQ "+plotname)

    # ====== Plot TS ============

    print("Plotting TS...")

    res_ts = dict(res)

    cmap = read_colormap_file(os.path.join(PATH_TO_NCL, "colormaps", "amwg256.rgb"))
    res_ts["cnFillPalette"] = cmap

    res_ts["cnLevels"] = manual_levels(-5.0, 110.0, 2.5)

    res_ts["gsnRightString"] = "Surface temperature (degF)"

    fig, ax = gsn_csm_contour_map(tsf[0,:,:], res_ts)

    res_fzTS = dict(res_contour)
    res_fzTS["cnLevels"] = [ 32. ]      # set levels
    res_fzTS["cnLineLabelsOn"] = False
    res_fzTS["cnLineColor"] = "#CD0000"   # red3
    res_fzTS["cnLineThicknessF"] = 4.     # doubles thickness of zero contour
    gsn_csm_contour(ax, tsf[0,::2,::2], res_fzTS, x=lon[::2], y=lat[::2])

    txres = {}
    txres["fontsize"] = 8             # Set the font height

    stride = tsstride

    gsn_add_text_grid(ax, tsf[0,:,:], stride, -np.inf, txres)

    frame(fig, "ts_"+plotname+"_"+yearstr+monstr+daystr+hourstr)

    del cmap, stride, txres

    wallClockElapseTime(wcStrt, "TS "+plotname)

    # ====== Plot TS ANOM ============

    print("Plotting TS ANOM...")

    res_tsanom = dict(res)

    cmap = read_colormap_file(os.path.join(PATH_TO_NCL, "colormaps", "MPL_RdBu.rgb"))
    res_tsanom["cnFillPalette"] = cmap[::-1,:]

    res_tsanom["cnLevels"] = manual_levels(-20.0, 20.0, 2.0)

    res_tsanom["gsnRightString"] = "Sfc. temp. anom. (1980-2009 NCEP base) (degF)"

    # area average to every 4th point, then plot on that coarser grid
    tsf_tmp = area_hi2lores(lon, lat, tsf_anom[0,:,:], lon[::4], lat[::4])
    lon_full, lat_full = lon, lat
    lon, lat = lon[::4], lat[::4]
    fig, ax = gsn_csm_contour_map(tsf_tmp, res_tsanom)
    lon, lat = lon_full, lat_full

    frame(fig, "tsanom_"+plotname+"_"+yearstr+monstr+daystr+hourstr)

    del tsf_tmp, cmap

    wallClockElapseTime(wcStrt, "TS anom "+plotname)

    # ====== Plot FLUT ============

    print("Plotting FLUT...")

    res_flut = dict(res)

    cmap = read_colormap_file(os.path.join(PATH_TO_NCL, "colormaps", "MPL_gist_yarg.rgb"))
    res_flut["cnFillPalette"] = cmap

    res_flut["cnLevels"] = manual_levels(120.0, 250.0, 10.)

    res_flut["mpFillOn"] = True
    res_flut["mpOceanFillColor"] = "black"
    res_flut["mpLandFillColor"] = "darkgreen"
    res_flut["mpInlandWaterFillColor"] = "black"

    res_flut["gsnRightString"] = "Outgoing longwave rad. (W/m2)"

    res_flut["mpGeophysicalLineColor"] = "#00FF00"   # NCL green

    #flut = np.where(flut >= 245, np.nan, flut)
    fig, ax = gsn_csm_contour_map(flut[0,:,:], res_flut)
    frame(fig, "flut_"+plotname+"_"+yearstr+monstr+daystr+hourstr)

    del res_flut, cmap

    wallClockElapseTime(wcStrt, "FLUT "+plotname)

    # ====== Plot 850 WIND ============

    print("Plotting 850mb wind + MSLP...")

    res_wind = dict(res)

    cmap = read_colormap_file(os.path.join(PATH_TO_NCL, "colormaps", "WhBlGrYeRe.rgb"))
    res_wind["cnFillPalette"] = cmap

    res_wind["cnLevels"] = manual_levels(4.0, 46.0, 2.0)

    res_wind["gsnRightString"] = "MSLP/850mb wind"

    fig, ax = gsn_csm_contour_map(wind[0,:,:], res_wind)
    gsn_csm_contour(ax, mslp[0,:,:], res_mslp)

    frame(fig, "wind_"+plotname+"_"+yearstr+monstr+daystr+hourstr)

    del res_wind, cmap

    wallClockElapseTime(wcStrt, "850mb/MSLP "+plotname)

    # ========================= Plot 250-850 SHEAR ===============================
    # (shear850250 not computed, U250/V250 not read in)

    # delete res at end
    del res
