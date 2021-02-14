;========================================================
; TC_cross_and_latlon.ncl:   plots the PS, U, V, and T.
;========================================================
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/esmf/ESMF_regridding.ncl"
;========================================================

begin

anltype="gfs"   ; era5 or gfs
anlres="0.50" ;0.50 or 0.25

dstGridFile="conus_30_x8.g_scrip.nc"
gridName="conus_30_x8"

InterpMethod="patch"

if (InterpMethod .eq. "patch") then
  shortInterpName = "patc"
else
  shortInterpName=InterpMethod
end if

dstGridDir="/glade/work/zarzycki/grids/scrip/"

srcGridDir="./scrip/"
srcGridFile=anltype+"_"+anlres+"x"+anlres+"_scrip.nc"

wgtFileDir="/glade/work/zarzycki/maps/gfsmaps/"
wgtFileName="map_"+anltype+"-"+anlres+"_TO_"+gridName+"_"+shortInterpName+".nc"

;----------------------------------------------------------------------

srcGridName=srcGridDir+"/"+srcGridFile
dstGridName=dstGridDir+"/"+dstGridFile

;----------------------------------------------------------------------
; Now generate weights file
;----------------------------------------------------------------------

Opt                      = True
Opt@InterpMethod         = InterpMethod
Opt@ForceOverwrite       = True
Opt@PrintTimings         = True
Opt@NoPETLog             = True
Opt@RemovePETLog         = True
Opt@Debug = True

ESMF_regrid_gen_weights(srcGridName,dstGridName,wgtFileDir+"/"+wgtFileName,Opt)

delete(Opt)

end
