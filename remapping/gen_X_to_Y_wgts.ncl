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

InterpMethod="patch"   ;bilinear, patch, conserve
srcName="era5_0.25"
srcGridName="./anl_scrip/era5_0.25x0.25_scrip.nc"
dstName="360x720_cruncep"
dstGridName="/glade/u/home/zarzycki/betacast/remapping/360x720_cruncep_SCRIP.nc"
wgtFileDir="/glade/scratch/zarzycki/"


if (InterpMethod .eq. "patch") then
  interpString="patc"
end if
if (InterpMethod .eq. "bilinear") then
  interpString="blin"
end if
if (InterpMethod .eq. "conserve") then
  interpString="aave"
end if
cdate=systemfunc("date +%y%m%d")
wgtFileName="map_"+srcName+"_TO_"+dstName+"_"+interpString+"."+cdate+".nc"


;----------------------------------------------------------------------
; Now generate weights file
;----------------------------------------------------------------------

Opt                      = True
Opt@InterpMethod         = InterpMethod
Opt@ForceOverwrite       = True
Opt@PrintTimings         = True

; check for regional domains
if (srcName .eq. "gland4km") then
  Opt@SrcRegional=True
end if
if (dstName .eq. "gland4km") then
  Opt@DstRegional=True
end if

ESMF_regrid_gen_weights(srcGridName,dstGridName,wgtFileDir+"/"+wgtFileName,Opt)

delete(Opt)

end
