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

;dstGridFile="ne30np4_091226_pentagons.nc"
dstGridFile="ne120np4_SCRIP.nc"
GFSres="0.25" ;0.50 or 0.25
gridName="ne120np4"
InterpMethod="patch"

if (InterpMethod .eq. "patch") then
  shortInterpName = "patc"
else
  shortInterpName=InterpMethod
end if

dstGridDir="/glade/work/zarzycki/grids/scrip/"

srcGridDir="/glade/work/zarzycki/grids/scrip/"
srcGridFile="ERA5_SCRIP.nc"

wgtFileDir="/glade/work/zarzycki/maps/gfsmaps/"
wgtFileName="map_era"+GFSres+"_TO_"+gridName+"_"+shortInterpName+".nc"

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
