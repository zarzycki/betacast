if ismissing(getenv("BETACAST")) then
  print("We are running local-only atm_to_cam. Set export BETACAST env to run elsewhere")
  ; if we don't have BETACAST, check to see where are located...
  ; we cannot do relative paths
  this_script_name=get_script_name()
  slash_index=str_index_of_substr(this_script_name, "/", -1)
  if (.not. ismissing(slash_index)) then
    print("The script is called from another directory")
    print("I found: "+this_script_name+" which contains a /")
    print("This isn't supported! Either set BETACAST as env param or run from local dir!")
    exit
  else
    print("The script is called from the local directory")
    PATHTOHERE=systemfunc("pwd")
  end if
  ; get BETACAST based on current dir assuming one up
  BETACAST = systemfunc("echo $(dirname $(pwd))")
else
  print("Not local only!")
  PATHTOHERE=getenv("BETACAST")+"/land-spinup/gen_datm/gen-datm/"
  BETACAST=getenv("BETACAST")
end if

begin

;d=addfile("/glade/p/cgd/tss/CTSM_datm_forcing_data/atm_forcing.datm7.cruncep_qianFill.0.5d.v7.c160715/domain.lnd.360x720_cruncep.130305.nc","r")
;d=addfile("/glade/p/cgd/tss/CTSM_datm_forcing_data/atm_forcing.datm7.GSWP3.0.5d.v1.c170516/domain.lnd.360x720_gswp3.0v1.c170606.nc","r")
dfilename=PATHTOHERE+"/era5-domain.nc"
print("Loading domain file: "+dfilename)
d=addfile(dfilename,"r")
timestride=1
STIX=0
do_q=True
do_flds=True
;outdirbase="/glade/derecho/scratch/zarzycki/ERA5-DATM/DATM_FORCING/ERA5/"
;pathtorawera5="/glade/u/home/zarzycki/scratch/ERA5-DATM/"
datafilename="CMZERA5.v0.c2021.0.5d"
;/glade/p/cgd/tss/CTSM_datm_forcing_data/atm_forcing.datm7.cruncep_qianFill.0.5d.v7.c160715/TPHWL6Hrly/clmforc.cruncep.V7.c2016.0.5d.TPQWL.1986-08.nc

; CLI args!! (uncomment if working just on this script)
;RAWERA5FILE="/glade/u/home/zarzycki/scratch/ERA5-DATM/out.1992.01.nc"
;YYYY="1992"
;MM="01"

greg_to_noleap=True

newgrid=False
wgt_filename="/glade/u/home/zarzycki/scratch/map_era5_0.25_TO_360x720_cruncep_patc.210315.nc"

; ----- internal settings
filenames = (/"TPQWL","Solar","Prec"/)
outfolders = (/"TPQW","Solar","Precip"/)

;==============================================================

is_leap_month=False
if ( isleapyear(toint(YYYY)) .and. toint(MM) .eq. 2) then
  print("We have a February leap year for "+YYYY+" "+MM)
  is_leap_month=True
end if

;==============================================================

print("Adding "+RAWERA5FILE+" ...")
;f=addfile(pathtorawera5+"out."+YYYY+"."+MM+".nc","r")
f=addfile(RAWERA5FILE,"r")

;mixhum_ptd
;"/glade/campaign/collections/rda/data/ds633.0/e5.oper.fc.sfc.meanflux/198608/e5.oper.fc.sfc.meanflux.235_058_msdrswrf.ll025sc.1986080106_1986081606.nc"
;msdwswrf
;msdwlwrf
;https://confluence.ecmwf.int/pages/viewpage.action?pageId=82870405#ERA5:datadocumentation-Parameterlistings

print("Getting coordinates from domain file...")
lon_datm_2D=d->xc
lat_datm_2D=d->yc

; assuming domain is regular lat-lon grid, pull first lat/lon to get 1D lat/lon arrays
lon_datm = lon_datm_2D(0,:)
lat_datm = lat_datm_2D(:,0)

print("Reading relevant coordinates from reanalysis...")
lon_era5=f->longitude
lat_era5=f->latitude(::-1)

print("Figuring out ENIX")
ENIX=dimsizes(f->time)-1
if (is_leap_month .and. greg_to_noleap) then
  ; if this is Feb/leap and we want no leap, find all days that have DD "29" and cut off
  ; ENIX before we pull those from the raw data
  print("Correcting ENIX from "+ENIX+" to...")
  tmp_times = cd_calendar(f->time,0)
  tmp_days = toint(tmp_times(:,2))
  ENIX = num(tmp_days .le. 28) - 1
  print("... to "+ENIX)
end if

; get full time
time_era5=f->time(STIX:ENIX:timestride)
print(cd_calendar(time_era5,-3))

print("Pull required variables off files...")
; get ERA5 data
u10_era5 = short2flt( f->u10(STIX:ENIX:timestride,::-1,:) )
v10_era5 = short2flt( f->v10(STIX:ENIX:timestride,::-1,:) )
tbot_era5 = short2flt( f->t2m(STIX:ENIX:timestride,::-1,:) )
tdew_era5 = short2flt( f->d2m(STIX:ENIX:timestride,::-1,:) )
psrf_era5 = short2flt( f->sp(STIX:ENIX:timestride,::-1,:) )
prec_era5 = short2flt( f->mtpr(STIX:ENIX:timestride,::-1,:) )
fsds_era5 = short2flt( f->ssrd(STIX:ENIX:timestride,::-1,:) )
if (do_flds) then
  flds_era5 = short2flt( f->strd(STIX:ENIX:timestride,::-1,:) )
end if

print("READ FROM ERA5 STATS:")
print("u10_era5 max: "+max(u10_era5)+"   min: "+min(u10_era5))
print("v10_era5 max: "+max(v10_era5)+"   min: "+min(v10_era5))
print("tbot_era5 max: "+max(tbot_era5)+"   min: "+min(tbot_era5))
print("tdew_era5 max: "+max(tdew_era5)+"   min: "+min(tdew_era5))
print("psrf_era5 max: "+max(psrf_era5)+"   min: "+min(psrf_era5))
print("prec_era5 max: "+max(prec_era5)+"   min: "+min(prec_era5))
print("fsds_era5 max: "+max(fsds_era5)+"   min: "+min(fsds_era5))
if (do_flds) then
  print("flds_era5 max: "+max(flds_era5)+"   min: "+min(flds_era5))
end if

print("Processing wind...")
wind_era5 = u10_era5
wind_era5 = sqrt(u10_era5^2. + v10_era5^2.)

if (do_q) then
  print("Processing Q...")
  qbot_era5 = mixhum_ptd(psrf_era5, tdew_era5, 2)
  copy_VarCoords(tdew_era5,qbot_era5)
end if

if (newgrid) then
  if (wgt_filename .ne. "") then
    print("Using ESMF mapping")
    Opt = True
    Opt@CopyVarCoords = False  ; turn this off to allow overwrite of coord types
    tbot_datm  = ESMF_regrid_with_weights(tbot_era5(:,::-1,:),wgt_filename,Opt)
    if (do_q) then
      qbot_datm = ESMF_regrid_with_weights(qbot_era5(:,::-1,:),wgt_filename,Opt)
    else
      qbot_datm = ESMF_regrid_with_weights(tdew_era5(:,::-1,:),wgt_filename,Opt)
    end if
    wind_datm  = ESMF_regrid_with_weights(wind_era5(:,::-1,:),wgt_filename,Opt)
    psrf_datm  = ESMF_regrid_with_weights(psrf_era5(:,::-1,:),wgt_filename,Opt)
    prec_datm  = ESMF_regrid_with_weights(prec_era5(:,::-1,:),wgt_filename,Opt)
    fsds_datm  = ESMF_regrid_with_weights(fsds_era5(:,::-1,:),wgt_filename,Opt)
    if (do_flds) then
      flds_datm  = ESMF_regrid_with_weights(flds_era5(:,::-1,:),wgt_filename,Opt)
    end if
  else
    print("Bilinear interpolation")
    tbot_datm = linint2(lon_era5,lat_era5,tbot_era5,True,lon_datm,lat_datm, 0)
    if (do_q) then
      qbot_datm = linint2(lon_era5,lat_era5,qbot_era5,True,lon_datm,lat_datm, 0)
    else
      qbot_datm = linint2(lon_era5,lat_era5,tdew_era5,True,lon_datm,lat_datm, 0)
    end if
    wind_datm = linint2(lon_era5,lat_era5,wind_era5,True,lon_datm,lat_datm, 0)
    psrf_datm = linint2(lon_era5,lat_era5,psrf_era5,True,lon_datm,lat_datm, 0)
    prec_datm = linint2(lon_era5,lat_era5,prec_era5,True,lon_datm,lat_datm, 0)
    fsds_datm = linint2(lon_era5,lat_era5,fsds_era5,True,lon_datm,lat_datm, 0)
    if (do_flds) then
      flds_datm = linint2(lon_era5,lat_era5,flds_era5,True,lon_datm,lat_datm, 0)
    end if
  end if
else
  print("Copying vars")
  tbot_datm  = tbot_era5
  delete(tbot_era5)
  if (do_q) then
    qbot_datm = qbot_era5
  else
    qbot_datm = tdew_era5
  end if
  wind_datm  = wind_era5
  delete(wind_era5)
  psrf_datm  = psrf_era5
  delete(psrf_era5)
  prec_datm  = prec_era5
  delete(prec_era5)
  fsds_datm  = fsds_era5
  delete(fsds_era5)
  if (do_flds) then
    flds_datm  = flds_era5
    delete(flds_era5)
  end if
end if

print("Enforce zero checks")
eps = 1.0e-8
prec_datm = where(prec_datm .lt. 0.0, 0.0, prec_datm)
fsds_datm = where(fsds_datm .le. 0.0, eps, fsds_datm)
if (do_flds) then
  flds_datm = where(flds_datm .le. 0.0, eps, flds_datm)
end if
if (do_q) then
  qbot_datm = where(qbot_datm .le. 0.0, eps, qbot_datm)
end if

print("Enforce specified capping checks")

maxwind=45.0
minwind=0.0
print("wind_datm max: "+max(wind_datm)+"   min: "+min(wind_datm))
print("min/max set by user: "+minwind+" "+maxwind)
print("Number of wind_datm over max to be corrected: "+num(wind_datm .gt. maxwind))
print("Number of wind_datm under min to be corrected: "+num(wind_datm .lt. minwind))
wind_datm = where(wind_datm .gt. maxwind, maxwind, wind_datm)
wind_datm = where(wind_datm .lt. minwind, minwind, wind_datm)

maxpsrf=120000.0
minpsrf=30000.0
print("psrf_datm max: "+max(psrf_datm)+"   min: "+min(psrf_datm))
print("min/max set by user: "+minpsrf+" "+maxpsrf)
print("Number of psrf_datm over max to be corrected: "+num(psrf_datm .gt. maxpsrf))
print("Number of psrf_datm under min to be corrected: "+num(psrf_datm .lt. minpsrf))
psrf_datm = where(psrf_datm .gt. maxpsrf, maxpsrf, psrf_datm)
psrf_datm = where(psrf_datm .lt. minpsrf, minpsrf, psrf_datm)

maxtbot=340.0
mintbot=175.0
print("tbot_datm max: "+max(tbot_datm)+"   min: "+min(tbot_datm))
print("min/max set by user: "+mintbot+" "+maxtbot)
print("Number of tbot_datm over max to be corrected: "+num(tbot_datm .gt. maxtbot))
print("Number of tbot_datm under min to be corrected: "+num(tbot_datm .lt. mintbot))
tbot_datm = where(tbot_datm .gt. maxtbot, maxtbot, tbot_datm)
tbot_datm = where(tbot_datm .lt. mintbot, mintbot, tbot_datm)

if (do_q) then
  maxqbot=0.2
  print("qbot_datm max: "+max(qbot_datm)+"   min: "+min(qbot_datm))
  print("max set by user: "+maxqbot)
  print("Number of qbot_datm over max to be corrected: "+num(qbot_datm .gt. maxqbot))
  qbot_datm = where(qbot_datm .gt. maxqbot, maxqbot, qbot_datm)
end if

print("Dummy variables")
zbot_datm = tbot_datm
zbot_datm = 10.
zbot_datm@long_name = "reference height"

print("Convert units")
; convert from J/s to W/m2 (over 1 hour)
fsds_datm = fsds_datm / 3600.
if (do_flds) then
  flds_datm = flds_datm / 3600.
end if

print("Set metadata")
tbot_datm!0="time"
tbot_datm!1="lat"
tbot_datm!2="lon"

print("Set time units")
if (greg_to_noleap) then
  new_time_units = "days since 1900-01-01 00:00:00"
  utc_date = cd_calendar(time_era5,0)
  time_opt=0
  time_opt@calendar="noleap"
  time_opt@return_type="float"
  time_datm=cd_inv_calendar(tointeger(utc_date(:,0)),tointeger(utc_date(:,1)),tointeger(utc_date(:,2)),tointeger(utc_date(:,3)),tointeger(utc_date(:,4)),utc_date(:,5),new_time_units,time_opt)
  delete(utc_date)
  delete(time_opt)
else
  new_time_units = "days since 1986-08-01 00:00:00"
  time_old_float = tofloat(time_era5)
  copy_VarMeta(time_era5,time_old_float)
  time_datm = time_to_newtime(time_old_float, new_time_units)
end if
time_datm!0="time"
time_datm&time = time_datm
time_datm@long_name = "observation time"

print("Convert types")
; 1-D lats and lons
lat_datm_float = tofloat(lat_datm)
copy_VarMeta(lat_datm,lat_datm_float)
lon_datm_float = tofloat(lon_datm)
copy_VarMeta(lon_datm,lon_datm_float)
lat_datm_float@long_name = "latitude"
lat_datm_float@units = "degrees_north"
lat_datm_float@mode = "time-invariant"
lon_datm_float@long_name = "longitude"
lon_datm_float@units = "degrees_east"
lon_datm_float@mode = "time-invariant"

; 2-D lats and lons
lat_datm_2D_float = tofloat(lat_datm_2D)
copy_VarMeta(lat_datm_2D,lat_datm_2D_float)
lon_datm_2D_float = tofloat(lon_datm_2D)
copy_VarMeta(lon_datm_2D,lon_datm_2D_float)
lat_datm_2D_float!0="lat"
lat_datm_2D_float!1="lon"
lon_datm_2D_float!0="lat"
lon_datm_2D_float!1="lon"
lat_datm_2D_float@long_name = "latitude"
lat_datm_2D_float@units = "degrees_north"
lat_datm_2D_float@mode = "time-invariant"
lon_datm_2D_float@long_name = "longitude"
lon_datm_2D_float@units = "degrees_east"
lon_datm_2D_float@mode = "time-invariant"

print("Assigning coordinates to time-varying fields")
tbot_datm&lat  := lat_datm_float
tbot_datm&lon  := lon_datm_float
tbot_datm&time := time_datm

print("Copying var coords")
delete(wind_datm&time)
copy_VarCoords(tbot_datm,wind_datm)
delete(psrf_datm&time)
copy_VarCoords(tbot_datm,psrf_datm)
delete(qbot_datm&time)
copy_VarCoords(tbot_datm,qbot_datm)
delete(zbot_datm&time)
copy_VarCoords(tbot_datm,zbot_datm)
delete(prec_datm&time)
copy_VarCoords(tbot_datm,prec_datm)
delete(fsds_datm&time)
copy_VarCoords(tbot_datm,fsds_datm)
if (do_flds) then
  delete(flds_datm&time)
  copy_VarCoords(tbot_datm,flds_datm)
end if

print("Assigning units and other metadata for time-dependent variables")
tbot_datm@mode = "time-dependent"
tbot_datm@units = "K"
wind_datm@mode = "time-dependent"
wind_datm@units = "m/s"
psrf_datm@mode = "time-dependent"
psrf_datm@units = "Pa"
qbot_datm@mode = "time-dependent"
if (do_q) then
  qbot_datm@units = "kg/kg"
else
  qbot_datm@units = "K"
end if
zbot_datm@mode = "time-dependent"
zbot_datm@units = "m"
prec_datm@mode = "time-dependent"
prec_datm@units = "mm H2O / sec"
prec_datm@long_name = "PRECTmms total precipitation"
fsds_datm@mode = "time-dependent"
fsds_datm@units = "W/m**2"
if (do_flds) then
  flds_datm@mode = "time-dependent"
  flds_datm@units = "W/m**2"
end if

print("Set edge vars...")
EDGEW =   0.
EDGEE = 360.
EDGES = -90.
EDGEN =  90.
EDGEW!0 = "scalar"
EDGEE!0 = "scalar"
EDGES!0 = "scalar"
EDGEN!0 = "scalar"
EDGEW@long_name = "western edge in atmospheric data"
EDGEE@long_name = "eastern edge in atmospheric data"
EDGES@long_name = "southern edge in atmospheric data"
EDGEN@long_name = "northern edge in atmospheric data"
EDGEW@mode = "time-invariant"
EDGEE@mode = "time-invariant"
EDGES@mode = "time-invariant"
EDGEN@mode = "time-invariant"


;======= OUTPUT FILE


do ii = 0,dimsizes(filenames)-1

  fullfilename=outdirbase+"/"+outfolders(ii)+"/"+datafilename+"."+filenames(ii)+"."+tostring(YYYY)+"-"+tostring(MM)+".nc"

  print("Writing "+fullfilename+" output!")

  system("mkdir -p "+outdirbase+"/"+outfolders(ii))

  system("/bin/rm -f "+fullfilename)   ; remove any pre-existing file

  ;setfileoption("nc","format","netcdf4")

;!  ncdf = addfile(fullfilename ,"c")  ; open output netCDF file
;!
;!  ;===================================================================
;!  ; create global attributes of the file (optional)
;!  ;===================================================================
;!  fAtt               = True            ; assign file attributes
;!  fAtt@title         = "NCL Simple Approach to netCDF Creation"
;!  fAtt@source_file   =  "original-file.nc"
;!  fAtt@Conventions   = "None"
;!  fAtt@creation_date = systemfunc ("date")
;!  fileattdef( ncdf, fAtt )            ; copy file attributes
;!
;!  ;===================================================================
;!  ; make time an UNLIMITED dimension; recommended  for most applications
;!  ;===================================================================
;!  filedimdef(ncdf,"time",-1,True)
;!
;!  ;===================================================================
;!  ; output variables directly; NCL will call appropriate functions
;!  ; to write the meta data associated with each variable
;!  ;===================================================================
;!  ncdf->LONGXY = lon_datm_2D_float
;!  ncdf->LATIXY = lat_datm_2D_float
;!  ncdf->EDGEW  = EDGEW
;!  ncdf->EDGEE  = EDGEE
;!  ncdf->EDGES  = EDGES
;!  ncdf->EDGEN  = EDGEN
;!
;!  if (ii .eq. 0) then
;!    ncdf->TBOT  = tbot_datm
;!    ncdf->WIND  = wind_datm
;!    ncdf->PSRF  = psrf_datm
;!    if (do_q) then
;!     ncdf->QBOT  = qbot_datm
;!    else
;!     ncdf->TDEW  = qbot_datm
;!    end if
;!    ncdf->ZBOT  = zbot_datm
;!    ncdf->FLDS  = flds_datm
;!  else if (ii .eq. 2) then
;!    ncdf->PRECTmms  = prec_datm
;!  else
;!    ncdf->FSDS  = fsds_datm
;!  end if
;!  end if

;===================================================================
    ; Assume variables T, PS and ORO exist and that they have
    ; associated meta data: (a) coordinate variables time, lev, lat, lon
    ; and (b) attributes
    ;===================================================================

  thesedims = dimsizes(tbot_datm)
  ntim  = thesedims(0)
	nlat  = thesedims(1)
	nlon  = thesedims(2)

  ncdf = addfile(fullfilename ,"c")  ; open output netCDF file

  ;===================================================================
  ; explicitly declare file definition mode. Improve efficiency.
  ;===================================================================
  setfileoption(ncdf,"DefineMode",True)

  ;===================================================================
  ; create global attributes of the file
  ;===================================================================
  fAtt               = True
	fAtt@title         = "DATM forcing using ERA5: "+YYYY+" "+MM
	fAtt@source_file   = RAWERA5FILE
	fAtt@Conventions   = "None"
	fAtt@creation_date = systemfunc ("date")
	fAtt@notes   = "Generated with Betacast toolkit"
	fileattdef( ncdf, fAtt )            ; copy file attributes

    ;===================================================================
    ; predefine the coordinate variables and their dimensionality
    ; Note: to get an UNLIMITED record dimension, we set the dimensionality
    ; to -1 (or the actual size) and set the dimension name to True.
    ;===================================================================
  dimNames = (/"time", "lat", "lon", "scalar"/)
	dimSizes = (/ -1   ,  nlat,  nlon , 1 /)
	dimUnlim = (/ True , False, False, False/)
	filedimdef(ncdf,dimNames,dimSizes,dimUnlim)


       filevardef(ncdf, "time" ,typeof(time_datm),getvardims(time_datm))
       filevardef(ncdf, "LONGXY"  ,typeof(lon_datm_2D_float),getvardims(lon_datm_2D_float))
       filevardef(ncdf, "LATIXY"  ,typeof(lat_datm_2D_float),getvardims(lat_datm_2D_float))
       filevardef(ncdf, "EDGEW"  ,typeof(EDGEW),getvardims(EDGEW))
       filevardef(ncdf, "EDGEE"  ,typeof(EDGEE),getvardims(EDGEE))
       filevardef(ncdf, "EDGES"  ,typeof(EDGES),getvardims(EDGES))
       filevardef(ncdf, "EDGEN"  ,typeof(EDGEN),getvardims(EDGEN))

       if (ii .eq. 0) then
         filevardef(ncdf, "TBOT"  ,typeof(tbot_datm),getvardims(tbot_datm))
         filevardef(ncdf, "WIND"  ,typeof(wind_datm),getvardims(wind_datm))
         filevardef(ncdf, "PSRF"  ,typeof(psrf_datm),getvardims(psrf_datm))
         filevardef(ncdf, "QBOT"  ,typeof(qbot_datm),getvardims(qbot_datm))
         filevardef(ncdf, "ZBOT"  ,typeof(zbot_datm),getvardims(zbot_datm))
         if (do_flds) then
           filevardef(ncdf, "FLDS"  ,typeof(flds_datm),getvardims(flds_datm))
         end if
       else if (ii .eq. 2) then
         filevardef(ncdf, "PRECTmms"  ,typeof(prec_datm),getvardims(prec_datm))
       else
         filevardef(ncdf, "FSDS"  ,typeof(fsds_datm),getvardims(fsds_datm))
       end if
       end if


    ;===================================================================
    ; Copy attributes associated with each variable to the file
    ; All attributes associated with each variable will be copied.
    ;====================================================================

       filevarattdef(ncdf,"time",time_datm)
       filevarattdef(ncdf,"LONGXY",lon_datm_2D_float)
       filevarattdef(ncdf,"LATIXY",lat_datm_2D_float)
       filevarattdef(ncdf,"EDGEW",EDGEW)
       filevarattdef(ncdf,"EDGEE",EDGEE)
       filevarattdef(ncdf,"EDGES",EDGES)
       filevarattdef(ncdf,"EDGEN",EDGEN)

       if (ii .eq. 0) then
         filevarattdef(ncdf, "TBOT"  , tbot_datm)
         filevarattdef(ncdf, "WIND"  , wind_datm)
         filevarattdef(ncdf, "PSRF"  ,psrf_datm)
         filevarattdef(ncdf, "QBOT"  ,qbot_datm)
         filevarattdef(ncdf, "ZBOT"  ,zbot_datm)
         if (do_flds) then
           filevarattdef(ncdf, "FLDS"  ,flds_datm)
         end if
       else if (ii .eq. 2) then
         filevarattdef(ncdf, "PRECTmms"  ,prec_datm)
       else
         filevarattdef(ncdf, "FSDS"  ,fsds_datm)
       end if
       end if

    ;===================================================================
    ; explicitly exit file definition mode. **NOT REQUIRED**
    ;===================================================================
        setfileoption(ncdf,"DefineMode",False)

    ;===================================================================
    ; output only the data values since the dimensionality and such have
    ; been predefined. The "(/", "/)" syntax tells NCL to only output the
    ; data values to the predefined locations on the file.
    ;====================================================================
       ncdf->time   = (/time_datm/)
       ncdf->LONGXY    = (/lon_datm_2D_float/)
       ncdf->LATIXY    = (/lat_datm_2D_float/)
       ncdf->EDGEW    = (/EDGEW/)
       ncdf->EDGEE      = (/EDGEE/)
       ncdf->EDGES      = (/EDGES/)
       ncdf->EDGEN  = (/EDGEN/)

       if (ii .eq. 0) then
         ncdf->TBOT   = (/tbot_datm/)
         ncdf->WIND   = (/wind_datm/)
         ncdf->PSRF   = (/psrf_datm/)
         ncdf->QBOT   = (/qbot_datm/)
         ncdf->ZBOT   = (/zbot_datm/)
         if (do_flds) then
           ncdf->FLDS   = (/flds_datm/)
         end if
       else if (ii .eq. 2) then
         ncdf->PRECTmms   = (/prec_datm/)
       else
         ncdf->FSDS   = (/fsds_datm/)
       end if
       end if

end do

status_exit(9)

end