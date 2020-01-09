load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/esmf/ESMF_regridding.ncl"

begin

src_file = addfile("/glade/scratch/zarzycki/diana-etc/1980/atl_30_x4_refine_nochem.cam.h3.1980-06-16-00000.nc","r")  

dstDir = "/glade/p/work/zarzycki/grids/scrip/"
dstGridName = "caminterp_0.5_NHonly_SCRIP.nc"

lat = src_file->lat(::2)
lon = src_file->lon(::2)

Opt                = True
Opt@ForceOverwrite = True
Opt@PrintTimings   = True
Opt@Title          = "Rectilinear 0.5deg NH grid"

rectilinear_to_SCRIP(dstDir+dstGridName,lat,lon,Opt)



end
