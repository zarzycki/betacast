load "$NCARG_ROOT/lib/ncarg/nclscripts/esmf/ESMF_regridding.ncl"

begin

; =*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*
; ============= User settings

; hrrr_3km_scrip.nc, era5_0.25x0.25, era5_0.3gaus, gfs_0.50x0.50, gfs_0.25x0.25, rap_13km, hrrr_3km, hwrf_storm
anlgrid="era5_0.25x0.25"
anlgridpath="../grids/anl_scrip/"

; shortname for grid (for naming purposes)
dstGridName="ne16np4"
; full path to model SCRIP file
;dstGridFile="/glade/u/home/zarzycki/work/ASD2017_files/grids/scrip/mp120a_grid_140708.nc"
dstGridFile="../grids/model_scrip/ne16np4_scrip_171002.nc"
;dstGridFile="/glade/work/zarzycki/grids/scrip/ne0np4natlanticwat.ne30x4.g_scrip.nc"
;dstGridFile="./rrm.x4_ESMF.nc"

; where to save the map file?
wgtFileDir="./"

; set this flag to True if you need to generated a CAM->ANL grid (mainly for model downscaling)
flip_model_and_analysis=False

; =*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*
; ============= possible CLI

; check if things are brought in
if isvar("ANLGRID")
  anlgrid=ANLGRID
end if
if isvar("ANLGRIDPATH")
  anlgridpath=ANLGRIDPATH
end if
if isvar("DSTGRIDNAME")
  dstGridName=DSTGRIDNAME
end if
if isvar("DSTGRIDFILE")
  dstGridFile=DSTGRIDFILE
end if
if isvar("WGTFILEDIR")
  wgtFileDir=WGTFILEDIR
end if
if isvar("FLIP_MODEL_AND_ANALYSIS")
  flip_model_and_analysis=FLIP_MODEL_AND_ANALYSIS
end if

; =*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*
; ============= do not edit below this line

if (anlgrid .ne. "era5_0.25x0.25" .and. anlgrid .ne. "era5_0.3gaus" .and. anlgrid .ne. "gfs_0.25x0.25" .and. anlgrid .ne. "gfs_0.50x0.50" .and. anlgrid .ne. "rap_13km" .and. anlgrid .ne. "hrrr_3km" .and. anlgrid .ne. "hwrf_storm") then
  print("unsupported analysis grid!   : "+anlgrid)
  exit
end if
srcGridFile=anlgridpath+"/"+anlgrid+"_scrip.nc"

; if using a regional grid, ensure that ESMF knows that...
regionalGrid=False
if (anlgrid .eq. "rap_13km" .or. anlgrid .eq. "hrrr_3km" .or. anlgrid .eq. "hwrf_storm") then
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


;----------------------------------------------------------------------
; Now generate weights file
;----------------------------------------------------------------------

print("Generating weights!")

if (.not.flip_model_and_analysis) then
  wgtFileName="map_"+anlgrid+"_TO_"+dstGridName+"_"+shortInterpName+".nc"
  srcGridName=srcGridFile
  dstGridName=dstGridFile
else
  wgtFileName="map_"+dstGridName+"_TO_"+anlgrid+"_"+shortInterpName+".nc"
  srcGridName=dstGridFile
  dstGridName=srcGridFile
end if

Opt                      = True
Opt@InterpMethod         = InterpMethod
Opt@ForceOverwrite       = True
Opt@PrintTimings         = True
Opt@NoPETLog             = True
Opt@RemovePETLog         = True
Opt@Debug = True
Opt@NetCDFType="netcdf4"

if (regionalGrid) then
  if (.not.flip_model_and_analysis) then
    Opt@SrcRegional=True
  else
    Opt@DstRegional=True
  end if
end if

if (dstType .eq. "model" .or. dstType .eq. "esmf") then
  if (.not.flip_model_and_analysis) then
    Opt@DstESMF = True
    if (dstType .eq. "model") then
      dstGridName = seGridName
    end if
  else
    Opt@SrcESMF = True
    if (dstType .eq. "model") then
      srcGridName = seGridName
    end if
  end if
end if

print(Opt)

  ESMF_regrid_gen_weights(srcGridName,dstGridName,wgtFileDir+"/"+wgtFileName,Opt)


delete(Opt)
if (dstType .eq. "model") then
  print("removing online generated grid descriptor")
  system("rm "+seGridName)
end if

if (.not.ismissing(systemfunc("type ncks"))) then
  print("NCO compressing: "+wgtFileDir+"/"+wgtFileName)
  beforefilesize=tofloat(systemfunc("stat -c %s "+wgtFileDir+"/"+wgtFileName))
  system("ncks -O -4 -L 1 "+wgtFileDir+"/"+wgtFileName+" "+wgtFileDir+"/"+wgtFileName)
  afterfilesize=tofloat(systemfunc("stat -c %s "+wgtFileDir+"/"+wgtFileName))
  compressionpercentage = 100.0 * (1 - afterfilesize/beforefilesize)
  print("Compression percentage: "+compressionpercentage+" %")
else
  print("Cannot compress "+wgtFileDir+"/"+wgtFileName+", no ncks")
end if

print("successfully generated...   "+wgtFileDir+"/"+wgtFileName)

print("done")
status_exit(9)
end
