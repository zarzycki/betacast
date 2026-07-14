load "$NCARG_ROOT/lib/ncarg/nclscripts/esmf/ESMF_regridding.ncl"

begin

;STRING="013"

shortgrid="mp3a"
InterpMethod="patch"  ;bilinear, patch, conserve, nearestdtos, neareststod
regional=True
regdomain="tclf"
gridName=regdomain+STRING+"-mp3a"

; RLL GRID
;outres="0.25x0.25"
;outres="0.2x0.2"
;outres="0.125x0.125"
;outres="0.0625x0.0625"
outres="0.025x0.025"
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
if (shortgrid .eq. "mp3a") then
  srcGridFile="mpasa3-60-tclf"+STRING+"_scrip.nc"
end if
;srcGridFile="mpasa15natl_scrip.nc"
;srcGridFile="uniform_240.g_scrip.nc"
;srcGridDir="/glade/p/cesmdata/cseg/mapping/grids/"
;srcGridFile="ne120np4_pentagons_100310.nc"

; RLL GRID
;dstGridDir="."
dstGridDir="/glade/derecho/scratch/zarzycki/"
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

LATHALFRANGE=10.0
LONHALFRANGE=12.0

; ASD regional
if (regional) then
  if (regdomain .eq. "asd") then
    Opt@LLCorner      = (/10.0d, 230.0d /)
    Opt@URCorner      = (/55.0d, 299.0d /)
  elseif (regdomain .eq. "tclf") then
    if (STRING .eq. "001") then
      CENLAT=27.1
      CENLON=-95.7
    elseif (STRING .eq. "002") then
      CENLAT=28.4
      CENLON=-92.9
    elseif (STRING .eq. "003") then
      CENLAT=28.4
      CENLON=-89.4
    elseif (STRING .eq. "004") then
      CENLAT=29.3
      CENLON=-86.7
    elseif (STRING .eq. "005") then
      CENLAT=28.5
      CENLON=-84.2
    elseif (STRING .eq. "006") then
      CENLAT=26.1
      CENLON=-83.0
    elseif (STRING .eq. "007") then
      CENLAT=24.8
      CENLON=-79.8
    elseif (STRING .eq. "008") then
      CENLAT=28.0
      CENLON=-79.5
    elseif (STRING .eq. "009") then
      CENLAT=30.8
      CENLON=-80.1
    elseif (STRING .eq. "010") then
      CENLAT=32.7
      CENLON=-78.1
    elseif (STRING .eq. "011") then
      CENLAT=34.3
      CENLON=-75.5
    elseif (STRING .eq. "012") then
      CENLAT=37.3
      CENLON=-74.6
    elseif (STRING .eq. "013") then
      CENLAT=39.9
      CENLON=-73.1
    elseif (STRING .eq. "014") then
      CENLAT=41.9
      CENLON=-69.2
    else
    end if
    Opt@LLCorner      = (/(CENLAT-LATHALFRANGE),(CENLON+360.0)-LONHALFRANGE /)
    Opt@URCorner      = (/(CENLAT+LATHALFRANGE),(CENLON+360.0)+LONHALFRANGE /)

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

system("ncks -O -4 -L 1 "+wgtFileDir+"/"+wgtFileName+" "+wgtFileDir+"/"+wgtFileName)

end
