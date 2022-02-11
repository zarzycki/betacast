load "$NCARG_ROOT/lib/ncarg/nclscripts/esmf/ESMF_regridding.ncl"

begin

; =*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*
; ============= User settings

; hrrr_3km_scrip.nc, era5_0.25x0.25, gfs_0.50x0.50, gfs_0.25x0.25, rap_13km, hrrr_3km
anlgrid="hrrr_3km" 

; shortname for grid (for naming purposes)
dstGridName="rrm.x4_ESMF.nc"
; full path to model SCRIP file
;dstGridFile="/glade/u/home/zarzycki/work/grids/scrip/ne30np4_091226_pentagons.nc"
dstGridFile="./rrm.x4_ESMF.nc"

; where to save the map file?
wgtFileDir="./"

; =*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*
; ============= do not edit below this line

if (anlgrid .ne. "era5_0.25x0.25" .and. anlgrid .ne. "gfs_0.25x0.25" .and. anlgrid .ne. "gfs_0.50x0.50" .and. anlgrid .ne. "rap_13km" .and. anlgrid .ne. "hrrr_3km") then
  print("unsupported analysis grid!   : "+anlgrid)
  exit
end if
srcGridFile="./anl_scrip/"+anlgrid+"_scrip.nc"

; if using a regional grid, ensure that ESMF knows that...
regionalGrid=False
if (anlgrid .eq. "rap_13km" .or. anlgrid .eq. "hrrr_3km") then
  regionalGrid=True
end if

; figure out which file type is input by which variables are on the file...
f = addfile(dstGridFile,"r")
if ( isfilevar(f, "nodeCoords") .and. isfilevar(f, "elementConn")) then
  dstType="esmf"
else if ( isfilevar(f, "grid_center_lat") .and. isfilevar(f, "grid_corner_lat")) then
  dstType="scrip"
else if ( isfilevar(f, "eb_prop1") ) then
  dstType="exodus"
  print("WAIT: I think you have an exodus file since "+dstGridFile+" has eb_prop1 as a var on it.")
  print("Do not support Exodus, generate the physics grid (i.e., w/ ncol) and pass that in, please!")
  exit
else
  dstType="model"
end if
end if
end if
print("determined input grid type to be: "+dstType)

if (dstType .eq. "model") then
  print("Generating grid file for model!")
  
  if ( .not. isfilevar(f, "lat") .or. .not. isfilevar(f, "lon")) then
    print("dst file not ESMF or SCRIP nor does it have lat/lon coords, existing!")
  end if
  
  lat=f->lat
  lon=f->lon

  ; gen SE grid
  Opt_se               = True
  Opt_se@ForceOverwrite = True
  Opt_se@PrintTimings   = True
  Opt_se@Title          = "SE Grid"
  seGridName = "grid_se.nc"
  unstructured_to_ESMF(seGridName,lat,lon,Opt_se)
  
  delete(Opt_se)
  delete(lat)
  delete(lon)
  delete(f)
end if

InterpMethod="patch"
if (InterpMethod .eq. "patch") then
  shortInterpName = "patc"
else
  shortInterpName=InterpMethod
end if

wgtFileName="map_"+anlgrid+"_TO_"+dstGridName+"_"+shortInterpName+".nc"

;----------------------------------------------------------------------
; Now generate weights file
;----------------------------------------------------------------------

print("Generating weights!")

srcGridName=srcGridFile
dstGridName=dstGridFile

Opt                      = True
Opt@InterpMethod         = InterpMethod
Opt@ForceOverwrite       = True
Opt@PrintTimings         = True
Opt@NoPETLog             = True
Opt@RemovePETLog         = True
Opt@Debug = True
if (regionalGrid) then
  Opt@SrcRegional=True
end if
if (dstType .eq. "model" .or. dstType .eq. "esmf") then
  Opt@DstESMF = True
  if (dstType .eq. "model") then
    dstGridName = seGridName
  end if
end if

ESMF_regrid_gen_weights(srcGridName,dstGridName,wgtFileDir+"/"+wgtFileName,Opt)

delete(Opt)
if (dstType .eq. "model") then
  print("removing online generated grid descriptor")
  system("rm "+seGridName)
end if

print("successfully generated...   "+wgtFileDir+"/"+wgtFileName)

end
