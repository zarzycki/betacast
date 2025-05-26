load "$NCARG_ROOT/lib/ncarg/nclscripts/esmf/ESMF_regridding.ncl"

begin

INDEX="009"
CENLON=-80.1
CENLAT=30.8

gridName="mpasa3-60-tclf"+INDEX
InterpMethod="patch"  ;bilinear, patch, conserve, nearestdtos, neareststod
regional=True
regdomain="tclf"+INDEX

; RLL GRID
outres="0.025x0.025"
srcGridDir="/glade/u/home/zarzycki/work/grids/scrip/"
srcGridFile="mpasa3-60-tclf"+INDEX+"_scrip.nc"
dstGridDir="/glade/work/zarzycki/grids/scrip/"

if (regional) then
  dstGridFile=outres+"_reg_SCRIP.nc"
  wgtFileName="map_"+gridName+"_to_"+outres+"reg_"+InterpMethod+".nc"
else
  dstGridFile=outres+"_SCRIP.nc"
  wgtFileName="map_"+gridName+"_to_"+outres+"glob_"+InterpMethod+".nc"
end if
wgtFileDir="/glade/derecho/scratch/zarzycki/"

;----------------------------------------------------------------------

srcGridName=srcGridDir+"/"+srcGridFile
dstGridName=dstGridDir+"/"+dstGridFile

;----------------------------------------------------------------------
; First generate a regular SCRIP file
;----------------------------------------------------------------------

Opt                = True
Opt@ForceOverwrite = True
Opt@PrintTimings   = True
Opt@Debug = True

centerLat = CENLAT
if (CENLON .lt. 0.0) then
  centerLon = 360.0 + CENLON
else
  centerLon = CENLON
end if
width_div_2 = 5.0
Opt@LLCorner      = (/centerLat-width_div_2,centerLon-width_div_2 /)
Opt@URCorner      = (/centerLat+width_div_2,centerLon+width_div_2 /)

latlon_to_SCRIP(dstGridName,outres,Opt)

;---Clean up
delete(Opt)

;----------------------------------------------------------------------
; Now generate weights file
;----------------------------------------------------------------------

Opt                      = True
Opt@InterpMethod         = InterpMethod
Opt@ForceOverwrite       = True
Opt@PrintTimings         = True
Opt@SrcRegional=True
ESMF_regrid_gen_weights(srcGridName,dstGridName,wgtFileDir+"/"+wgtFileName,Opt)
print("Done writing: "+wgtFileDir+"/"+wgtFileName)
delete(Opt)

system("rm -fv "+dstGridName)

end
