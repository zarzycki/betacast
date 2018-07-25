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

dstGridFile="ne30np4_091226_pentagons.nc"
GFSres="0.50" ;0.50 or 0.25
gridName="ne30np4"
InterpMethod="patch"

if (InterpMethod .eq. "patch") then
  shortInterpName = "patc"
else
  shortInterpName=InterpMethod
end if

dstGridDir="/glade/p/work/zarzycki/grids/scrip/"

srcGridDir="/glade/p/work/zarzycki/grids/scrip/"
srcGridFile="gfs_"+GFSres+"x"+GFSres+"_scrip.nc"

wgtFileDir="/glade/p/work/zarzycki/maps/gfsmaps/"
wgtFileName="map_gfs"+GFSres+"_TO_"+gridName+"_"+shortInterpName+".nc"

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
