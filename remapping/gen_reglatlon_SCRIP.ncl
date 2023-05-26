load "$NCARG_ROOT/lib/ncarg/nclscripts/esmf/ESMF_regridding.ncl"

begin

srcfilename="/glade/work/zarzycki/2023_mawar/mawar02w.2023052300.hwrfprs.storm.0p015.f000.grb2"
dstDir = "./anl_scrip/"
dstGridName = "hwrf_storm_scrip.nc"
metatitle = "HWRF"
appendLL = False      ; append corner lat lon to SCRIP filename? (useful for regional nests that may move)

; ====================================================================================

; check if things are brought in
if isvar("SRCFILENAME")
  srcfilename=SRCFILENAME
end if
if isvar("DSTDIR")
  dstDir = DSTDIR
end if
if isvar("DSTGRIDNAME")
  dstGridName = DSTGRIDNAME
end if
if isvar("METATITLE")
  metatitle = METATITLE
end if
if isvar("APPENDLL")
  appendLL=True
end if

; ====================================================================================

; load data
src_file = addfile(srcfilename,"r")

; load lat/lon
lat = src_file->lat_0
lon = src_file->lon_0

; figure out top lat/lon corners
top_left_lat = lat(0)
top_left_lon = lon(0)

if appendLL
  dstGridName = dstGridName+"_"+top_left_lat+"_"+top_left_lon+".nc"
end if

print("Creating SCRIP file from: "+srcfilename)
print("Writing SCRIP file to: "+dstDir+"/"+dstGridName)
print("Metaname for SCRIP file is: "+metatitle)

Opt                = True
Opt@ForceOverwrite = True
Opt@PrintTimings   = True
Opt@Title          = metatitle

rectilinear_to_SCRIP(dstDir+"/"+dstGridName,lat,lon,Opt)

end
