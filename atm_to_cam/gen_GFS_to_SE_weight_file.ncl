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

gridName="ecsnow_30_x0"
InterpMethod="patch"

if (InterpMethod .eq. "patch") then
  shortInterpName = "patc"
else
  shortInterpName=InterpMethod
end if

dstGridDir="/glade/p/work/zarzycki/grids/scrip/"
dstGridFile=gridName+".g_scrip.nc"

srcGridDir="/glade/p/work/zarzycki/grids/scrip/"
srcGridFile="gfs_0.50x0.50_scrip.nc"

wgtFileDir="/glade/p/work/zarzycki/maps/gfsmaps/"
wgtFileName="map_gfs0.50_TO_"+gridName+"_"+shortInterpName+".nc"

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

ESMF_regrid_gen_weights(srcGridName,dstGridName,wgtFileDir+"/"+wgtFileName,Opt)

delete(Opt)

end
