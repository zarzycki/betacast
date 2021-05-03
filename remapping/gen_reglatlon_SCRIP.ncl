load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/esmf/ESMF_regridding.ncl"

begin

;src_file = addfile("/glade/scratch/zarzycki/diana-etc/1980/atl_30_x4_refine_nochem.cam.h3.1980-06-16-00000.nc","r")  
src_file = addfile("/glade/scratch/zarzycki/ERA5-DATM/DATM//TPHWL6Hrly/CMZERA5.v0.c2021.0.5d.TPQWL.1995-09.nc","r")  

dstDir = "./"
dstGridName = "360x720_cruncep_SCRIP.nc"

lat = src_file->lat
lon = src_file->lon

Opt                = True
Opt@ForceOverwrite = True
Opt@PrintTimings   = True
Opt@Title          = "360x720_cruncep for DATM"

rectilinear_to_SCRIP(dstDir+dstGridName,lat,lon,Opt)



end
