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

;gridName="fv0.23x0.31"
gridName="fv0.9x1.25"
InterpMethod="patch"

if (InterpMethod .eq. "patch") then
  shortInterpName = "patc"
else
  shortInterpName=InterpMethod
end if

dstGridDir="./fv_scrip/"
dstGridFile="SCRIPgrid_0.9x1.25_nomask_c191014.nc"

srcGridDir="/glade/u/home/zarzycki/betacast/remapping/scrip/"
srcGridFile="era5_0.25x0.25_scrip.nc"

wgtFileDir="./"
wgtFileName="map_era5_0.25_TO_"+gridName+"_"+shortInterpName+".nc"

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

print("successfully generated... "+wgtFileDir+"/"+wgtFileName)

end
