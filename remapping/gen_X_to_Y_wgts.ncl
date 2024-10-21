load "$NCARG_ROOT/lib/ncarg/nclscripts/esmf/ESMF_regridding.ncl"
;========================================================

begin

InterpMethod="patch"   ;bilinear, patch, conserve
dstName="era5_0.25"
dstGridName="../grids/anl_scrip/era5_0.25x0.25_scrip.nc"
srcName="rrmx4"
srcGridName="./rrm.x4_ESMF.nc"
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

; figure out which file type is input by which variables are on the file...
f = addfile(srcGridName,"r")
if ( isfilevar(f, "nodeCoords") .and. isfilevar(f, "elementConn")) then
  srcType="esmf"
else if ( isfilevar(f, "grid_center_lat") .and. isfilevar(f, "grid_corner_lat")) then
  srcType="scrip"
else if ( isfilevar(f, "eb_prop1") ) then
  srcType="exodus"
  print("WAIT: I think you have an exodus file since "+srcGridName+" has eb_prop1 as a var on it.")
  print("Do not support Exodus, generate the physics grid (i.e., w/ ncol) and pass that in, please!")
  exit
else
  srcType="model"
end if
end if
end if
print("determined source grid "+srcGridName+"   type to be: "+srcType)

; figure out which file type is input by which variables are on the file...
f = addfile(dstGridName,"r")
if ( isfilevar(f, "nodeCoords") .and. isfilevar(f, "elementConn")) then
  dstType="esmf"
else if ( isfilevar(f, "grid_center_lat") .and. isfilevar(f, "grid_corner_lat")) then
  dstType="scrip"
else if ( isfilevar(f, "eb_prop1") ) then
  dstType="exodus"
  print("WAIT: I think you have an exodus file since "+dstGridName+" has eb_prop1 as a var on it.")
  print("Do not support Exodus, generate the physics grid (i.e., w/ ncol) and pass that in, please!")
  exit
else
  dstType="model"
end if
end if
end if
print("determined destination grid "+dstGridName+"   type to be: "+dstType)
;----------------------------------------------------------------------
; Now generate weights file
;----------------------------------------------------------------------

Opt                      = True
Opt@InterpMethod         = InterpMethod
Opt@ForceOverwrite       = True
Opt@PrintTimings         = True

if (srcType .eq. "esmf") then
  Opt@SrcESMF = True
end if
if (dstType .eq. "esmf") then
  Opt@DstESMF = True
end if

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
