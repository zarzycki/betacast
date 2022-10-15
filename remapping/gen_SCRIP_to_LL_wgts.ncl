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

gridName="mp15-120-atlantic"
InterpMethod="patch"  ;bilinear, patch, conserve, nearestdtos, neareststod
regional=True
regdomain="florida" ;asd, atlantic, tctest, florida

; RLL GRID
;outres="0.25x0.25"
;outres="0.2x0.2"
outres="0.125x0.125"
;outres="0.0625x0.0625"
;outres="0.025x0.025"
;outres="0.5x0.5"
;outres="1x1"
srcGridDir="/glade/u/home/zarzycki/work/grids/scrip/"
;srcGridFile="NARR_SCRIP.nc"
;srcGridFile=gridName+".g_scrip.nc"
;srcGridDir="/glade/work/zarzycki/ASD2017_files/grids/scrip/"
;srcGridFile="mp15a-120a-US.grid.170118.nc"
;srcGridFile="mp120a_grid_140708.nc"
;srcGridFile="ne30pg3_scrip_170611.nc"
;srcGridFile="ne30np4_091226_pentagons.nc"
;srcGridFile="mpasa3-60-florida_scrip.nc"
srcGridFile="mpasa15natl_scrip.nc"
;srcGridFile="uniform_240.g_scrip.nc"
;srcGridDir="/glade/p/cesmdata/cseg/mapping/grids/"
;srcGridFile="ne120np4_pentagons_100310.nc"

; RLL GRID
;dstGridDir="."
dstGridDir="/glade/work/zarzycki/grids/scrip/"
if (regional) then
  dstGridFile=outres+"_reg_SCRIP.nc"
  wgtFileName="map_"+gridName+"_to_"+outres+"reg_"+InterpMethod+".nc"
else
  dstGridFile=outres+"_SCRIP.nc"
  wgtFileName="map_"+gridName+"_to_"+outres+"glob_"+InterpMethod+".nc"
end if
wgtFileDir="/glade/scratch/zarzycki/"

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

; ASD regional
if (regional) then
  if (regdomain .eq. "asd") then
    Opt@LLCorner      = (/10.0d, 230.0d /) 
    Opt@URCorner      = (/55.0d, 299.0d /)
  elseif (regdomain .eq. "atlantic") then
    Opt@LLCorner      = (/5.0d, 250.0d /)
    Opt@URCorner      = (/55.0d,355.0d /)
  elseif (regdomain .eq. "tctest") then
    Opt@LLCorner      = (/5.0d, 250.0d /)
    Opt@URCorner      = (/55.0d,355.0d /)
  elseif (regdomain .eq. "florida") then
    Opt@LLCorner      = (/15.0d,271.0d /)
    Opt@URCorner      = (/35.0d,296.0d /)
  elseif (regdomain .eq. "narr") then
    Opt@LLCorner      = (/12.0d,215.0d /)
    Opt@URCorner      = (/65.0d,298.0d /)
  else
  end if
end if

;if (fileexists(dstGridName)) then
;  print("Latlon scrip already exists...")
;else
  latlon_to_SCRIP(dstGridName,outres,Opt)
;end if

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

delete(Opt)

end
