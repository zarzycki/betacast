load "$NCARG_ROOT/lib/ncarg/nclscripts/esmf/ESMF_regridding.ncl"

begin

; =*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*
; ============= User settings

; era5_0.25x0.25, gfs_0.50x0.50, gfs_0.25x0.25
anlgrid="gfs_0.25x0.25" 

; shortname for grid (for naming purposes)
dstGridName="ne30"
; full path to model SCRIP file
dstGridFile="/glade/u/home/zarzycki/scratch/F-betacast-F2000climo/run/scrip/F-betacast-F2000climo.cam.h0.2017-09-07-00000.nc"

; where to save the map file?
wgtFileDir="./"

; =*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*
; ============= do not edit below this line

if (anlgrid .ne. "era5_0.25x0.25" .and. anlgrid .ne. "gfs_0.25x0.25" .and. anlgrid .ne. "gfs_0.50x0.50") then
  print("unsupported analysis grid!   : "+anlgrid)
  exit
end if
srcGridFile="./anl_scrip/"+anlgrid+"_scrip.nc"

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
