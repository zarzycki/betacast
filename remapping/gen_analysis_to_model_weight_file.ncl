load "$NCARG_ROOT/lib/ncarg/nclscripts/esmf/ESMF_regridding.ncl"

begin

; =*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*
; ============= User settings

anltype="era5"   ; era5 or gfs
anlres="0.25" ;0.50 or 0.25

; shortname for grid (for naming purposes)
dstGridName="conus_30_x8"
; full path to model SCRIP file
dstGridFile="/glade/work/zarzycki/grids/scrip/conus_30_x8.g_scrip.nc"

; where to save the map file?
wgtFileDir="/glade/work/zarzycki/maps/gfsmaps/"

; =*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*
; ============= do not edit below this line

srcGridFile="./anl_scrip/"+anltype+"_"+anlres+"x"+anlres+"_scrip.nc"

InterpMethod="patch"
if (InterpMethod .eq. "patch") then
  shortInterpName = "patc"
else
  shortInterpName=InterpMethod
end if

wgtFileName="map_"+anltype+"-"+anlres+"_TO_"+dstGridName+"_"+shortInterpName+".nc"

;----------------------------------------------------------------------

srcGridName=srcGridFile
dstGridName=dstGridFile

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

print("successfully generated... "+wgtFileDir+"/"+wgtFileName)

end
