load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/esmf/ESMF_regridding.ncl"

begin

;src_file = addfile("/glade/scratch/zarzycki/diana-etc/1980/atl_30_x4_refine_nochem.cam.h3.1980-06-16-00000.nc","r")  
src_file = addfile("/global/homes/w/wrliu/scratch/sewx_newld/ECMWF/ERA-Int_2012102500.nc","r")  

dstDir = "./"
dstGridName = "ERA5_SCRIP.nc"

lat = src_file->latitude
lon = src_file->longitude

Opt                = True
Opt@ForceOverwrite = True
Opt@PrintTimings   = True
Opt@Title          = "ERA5 regular 0.25deg grid"

rectilinear_to_SCRIP(dstDir+dstGridName,lat,lon,Opt)



end
