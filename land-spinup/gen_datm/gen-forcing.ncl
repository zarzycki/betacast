d=addfile("/glade/p/cgd/tss/CTSM_datm_forcing_data/atm_forcing.datm7.cruncep_qianFill.0.5d.v7.c160715/domain.lnd.360x720_cruncep.130305.nc","r")
timestride=6
do_q=True
;/glade/p/cgd/tss/CTSM_datm_forcing_data/atm_forcing.datm7.cruncep_qianFill.0.5d.v7.c160715/TPHWL6Hrly/clmforc.cruncep.V7.c2016.0.5d.TPQWL.1986-08.nc

print("Adding files...")
ubot_era5_f=addfile("/glade/u/home/zarzycki/rda/ds633.0/e5.oper.an.sfc/198608/e5.oper.an.sfc.128_165_10u.ll025sc.1986080100_1986083123.nc","r")
vbot_era5_f=addfile("/glade/u/home/zarzycki/rda/ds633.0/e5.oper.an.sfc/198608/e5.oper.an.sfc.128_166_10v.ll025sc.1986080100_1986083123.nc","r")
ps_era5_f=addfile("/glade/u/home/zarzycki/rda/ds633.0/e5.oper.an.sfc/198608/e5.oper.an.sfc.128_134_sp.ll025sc.1986080100_1986083123.nc","r")
tbot_era5_f=addfile("/glade/u/home/zarzycki/rda/ds633.0/e5.oper.an.sfc/198608/e5.oper.an.sfc.128_167_2t.ll025sc.1986080100_1986083123.nc","r")
tdew_era5_f=addfile("/glade/u/home/zarzycki/rda/ds633.0/e5.oper.an.sfc/198608/e5.oper.an.sfc.128_168_2d.ll025sc.1986080100_1986083123.nc","r")

;mixhum_ptd
;"/glade/u/home/zarzycki/rda/ds633.0/e5.oper.fc.sfc.meanflux/198608/e5.oper.fc.sfc.meanflux.235_058_msdrswrf.ll025sc.1986080106_1986081606.nc"
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
lon_era5=tbot_era5_f->longitude
lat_era5=tbot_era5_f->latitude(::-1)
time_era5=tbot_era5_f->time(3::timestride)

print("Pull required variables off files...")
; get ERA5 data
u10_era5 = ubot_era5_f->VAR_10U(3::timestride,::-1,:)
v10_era5 = vbot_era5_f->VAR_10V(3::timestride,::-1,:)
tbot_era5 = tbot_era5_f->VAR_2T(3::timestride,::-1,:)
tdew_era5 = tdew_era5_f->VAR_2D(3::timestride,::-1,:)
psrf_era5 = ps_era5_f->SP(3::timestride,::-1,:)

print("Processing wind...")
wind_era5 = u10_era5
wind_era5 = sqrt(u10_era5^2. + v10_era5^2.)

if (do_q) then
  print("Processing Q...")
  qbot_era5 = mixhum_ptd(psrf_era5, tdew_era5, 2)
end if

print("Bilinear interpolation")
tbot_datm = linint2(lon_era5,lat_era5,tbot_era5,True,lon_datm,lat_datm, 0)
if (do_q) then
  qbot_datm = linint2(lon_era5,lat_era5,qbot_era5,True,lon_datm,lat_datm, 0)
else
  qbot_datm = linint2(lon_era5,lat_era5,tdew_era5,True,lon_datm,lat_datm, 0)
end if
wind_datm = linint2(lon_era5,lat_era5,wind_era5,True,lon_datm,lat_datm, 0)
psrf_datm = linint2(lon_era5,lat_era5,psrf_era5,True,lon_datm,lat_datm, 0)

print("Dummy variables")
zbot_datm = tbot_datm
zbot_datm = 10.

print("Set metadata")
tbot_datm!0="time"
tbot_datm!1="lat"
tbot_datm!2="lon"

print("Set time units")
new_time_units = "days since 1986-08-01 00:00:00"
time_old_float = tofloat(time_era5)
copy_VarMeta(time_era5,time_old_float)
time_datm = time_to_newtime(time_old_float, new_time_units) 
time_datm!0="time"
time_datm&time = time_datm
time_datm@long_name = "observation time"
;time_datm@units = "noleap"

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
tbot_datm&lat = lat_datm_float
tbot_datm&lon = lon_datm_float
tbot_datm&time = time_datm

print("Copying var coords")
copy_VarCoords(tbot_datm,wind_datm)
copy_VarCoords(tbot_datm,psrf_datm)
copy_VarCoords(tbot_datm,qbot_datm)
copy_VarCoords(tbot_datm,zbot_datm)

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
print("Writing output!")
system("/bin/rm -f simple.nc")   ; remove any pre-existing file
ncdf = addfile("simple.nc" ,"c")  ; open output netCDF file

;===================================================================
; create global attributes of the file (optional)
;===================================================================
fAtt               = True            ; assign file attributes
fAtt@title         = "NCL Simple Approach to netCDF Creation"
fAtt@source_file   =  "original-file.nc"
fAtt@Conventions   = "None"
fAtt@creation_date = systemfunc ("date")
fileattdef( ncdf, fAtt )            ; copy file attributes

;===================================================================
; make time an UNLIMITED dimension; recommended  for most applications
;===================================================================
filedimdef(ncdf,"time",-1,True) 

;===================================================================
; output variables directly; NCL will call appropriate functions
; to write the meta data associated with each variable
;===================================================================
ncdf->LONGXY = lon_datm_2D_float
ncdf->LATIXY = lat_datm_2D_float
ncdf->TBOT  = tbot_datm
ncdf->WIND  = wind_datm
ncdf->PSRF  = psrf_datm
if (do_q) then
  ncdf->QBOT  = qbot_datm
else
  ncdf->TDEW  = qbot_datm
end if
ncdf->ZBOT  = zbot_datm

ncdf->EDGEW  = EDGEW
ncdf->EDGEE  = EDGEE
ncdf->EDGES  = EDGES
ncdf->EDGEN  = EDGEN
