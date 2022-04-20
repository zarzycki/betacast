load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"

begin

; Last updated 9/7/2017 by Colin Zarzycki

; This file takes either GDAS SST or NOAAOI SST and interpolates to 1x1 SST datafile for
; use in CESM F compsets as a datastream.
; NOTE: Given initdate it populates the entire datastream with that days SST/ice fields
; therefore this method is only suitable for short-term forecasts (<14 days) due to the assumption
; of persistence. Future work may use initial anomalies on top of a daily-varying SST field
; particularly for S2S work.

; NOAAOI EXAMPLE
;ncl sst_interp.ncl 'initdate="2015033000"' \
;      'datasource="NOAAOI"' \
;      'sstDataFile = "/glade/p/work/zarzycki/NOAAOI/sst.day.mean.2015.v2.nc"' \
;      'iceDataFile = "/glade/p/work/zarzycki/NOAAOI/icec.day.mean.2015.v2.nc"' \
;      'SST_write_file = "/glade/p/work/zarzycki/NOAAOI/sst_1x1.nc"'

;initdate=2012082400
;SST_write_dir = "/glade/p/work/zarzycki/sewx/SST/"
;datasource = "GDAS"

; Note: copy_orig = True is the old way where a file was copied and then overwritten
; Much "smarter" to just do this from scratch...
copy_orig = False
if (copy_orig) then
  print("copying SST file from repo to "+SST_write_file)
  system("cp ./sst_1x1.nc "+SST_write_file)
  readinFile=SST_write_file
else
  readinFile="./sst_1x1.nc"
end if

if (.not. isvar("predict_docn")) then
  do_anom=False  
else
  if (toint(predict_docn) .eq. 1) then
    do_anom=True ; this superimposes the observed SST + ice onto the climatological seasonal cycle
  else if (toint(predict_docn) .eq. 0) then
    do_anom=False
  else
    print("unknown predict_docn option passed in, must be 0 (false) or 1 (true).")
    exit
  end if
  end if
end if
print("do_anom (i.e., predict_docn) set to: "+do_anom)

print("running date "+initdate)
KtoC=273.15

; GDAS ice settings
TTHRESH=271.95   ; initial cut for ice vs. open ocean
smooth_ice=True  ; smooth the ice field?
smooth_iter=3    ; if yes, number of iters?

; Get YYYYMMDDHH from the initdate string
dtime_map = (/4,2,2,2/)
tst_dtime = str_split_by_length(initdate, dtime_map)
yyyy=toint(tst_dtime(0))
mm=toint(tst_dtime(1))
dd=toint(tst_dtime(2))
hh=toint(tst_dtime(3))

; Load SST files that we want to overwrite
print("Writing 1x1 SSTs")
if (copy_orig) then
  in = addfile(SST_write_file,"w")
else
  in = addfile(readinFile,"r")
end if
;print("Writing 0.23x0.31 SSTs")
;in = addfile(SST_write_dir+"sst_0.23x0.31.nc","w")

if (datasource .eq. "GDAS") then
  sst_file = addfile(sstDataFile,"r")
  sstlat = sst_file->lat_0
  sstlon = sst_file->lon_0
  sst_gfs = sst_file->TMP_P0_L1_GLL0(:,:)
  ; handle ice
  print("For GDAS, set all native cells less than "+TTHRESH+" to ice covered")
  ice_gfs = sst_gfs
  ice_gfs = where(sst_gfs .gt. TTHRESH , 0.0, 1.0)
  if (smooth_ice) then
    print("smoothing derived ice field "+smooth_iter+" times!")
    do SM = 0,smooth_iter-1
      ice_gfs = smth9_Wrap(ice_gfs, 0.50,  0.25, True)
    end do
  end if
  ; Convert SST from K -> degC
  sst_gfs = sst_gfs - KtoC
else if (datasource .eq. "NOAAOI") then
  sst_file = addfile(sstDataFile,"r")
  ice_file = addfile(iceDataFile,"r")
  time = sst_file->time
  date = cd_inv_calendar(yyyy,mm,dd,0,0,0,time@units,0)
  sst_gfs = sst_file->sst({date},:,:)
  ice_gfs = ice_file->icec({date},:,:)
  sstlat = sst_file->lat
  sstlon = sst_file->lon
else
end if
end if

fvlat  = in->lat
fvlon  = in->lon
numfvlat = dimsizes(fvlat)
numfvlon = dimsizes(fvlon)

; wks   = gsn_open_wks ("x11", "coneff")       ; open ps file
; gsn_define_colormap(wks,"gui_default")        ; choose color map
; res                 = True                    ; plot mods desired
; res@gsnDraw         = False                   ; don't draw yet
; res@gsnFrame        = False                   ; don't advance frame yet
; res@cnFillOn        = True                    ; turn on color
; res@gsnSpreadColors = True                    ; use full colormap
; plot = gsn_csm_hov(wks, sst_gfs, res) 
; draw(plot)
; frame(wks)

if (datasource .eq. "GDAS") then
  ; Need to flip raw GDAS lats/lons to match SST file
  sst_fv = linint2_Wrap (sstlon,sstlat(::-1),sst_gfs(::-1,:),True,fvlon,fvlat,0)
  ice_fv = linint2_Wrap (sstlon,sstlat(::-1),ice_gfs(::-1,:),True,fvlon,fvlat,0)
else if (datasource .eq. "NOAAOI") then
  sst_fv = linint2_Wrap(sstlon,sstlat,sst_gfs,True,fvlon,fvlat,0)
  ice_fv = linint2_Wrap(sstlon,sstlat,ice_gfs,True,fvlon,fvlat,0)
else
end if
end if

sst_fv = linmsg(sst_fv,-1)
sst_fv = linmsg_n(sst_fv,-1,0)

print("Sorting bad SST")
do i = 0,dimsizes(fvlat)-1
    do j = 0,dimsizes(fvlon)-1
            if (ismissing(sst_fv(i,j))) then
                print("Found missing SST")
                sst_fv(i,j) = 271.5
            else if (sst_fv(i,j) .gt. 500) then
                print("Found bad value")
                sst_fv(i,j) = 271.5
            end if
            end if
        sst_fv(i,j) = sst_fv(i,j); - 273.15
    end do
end do

print("Converting SST floats to doubles")
sst_fv_dbl = tofloat(sst_fv)
print("Copying SST metadata")
copy_VarMeta(sst_fv,sst_fv_dbl)
delete(sst_fv)
print("Correcting SST time records")
; Add time record
sst_fv_dbl_time = new((/12,numfvlat,numfvlon/),float)


if (do_anom) then

  print("SST_TO_CAM: Creating observationally-consistent seasonal cycle")
  ; create daily time array
  daily_time = fspan (0.0, 364.0, 365)
  daily_time@units = "days since 0000-01-01 00:00:00"
  daily_time@calendar = "365_day"
  
  ; read monthly climo SST and interp to daily SST
  climo_SST = in->SST_cpl
  climo_SST_daily = linint1_n_Wrap (climo_SST&time, climo_SST, True, daily_time, 0, 0)

  ; get DoY for the forecast day
  yyyy@calendar = "noleap"           ; also "365"
  doy = day_of_year(yyyy,mm,dd)  ; doy = (/61,60,60/)
  doyint = toint(doy) - 1

  ; calculate sst_anom as observed SST - that day's climo
  sst_anom = sst_fv_dbl
  sst_anom = sst_fv_dbl - climo_SST_daily(doyint,:,:)
  
  ; project that anomaly onto the monthly values
  do i = 0,11
    sst_fv_dbl_time(i,:,:) = climo_SST(i,:,:) + sst_anom(:,:)
  end do
  
else

  print("SST_TO_CAM: Creating time-invariant SST forcing")

  do i = 0,11
    sst_fv_dbl_time(i,:,:) = sst_fv_dbl
  end do
  
end if
  
print("Sorting bad ice")
do i = 0,dimsizes(fvlat)-1
  do j = 0,dimsizes(fvlon)-1
          if (ismissing(ice_fv(i,j))) then
              ;print("Found missing ice")
              ice_fv(i,j) = 0
          else if (ice_fv(i,j) .gt. 500) then
              print("Found bad ice")
              ice_fv(i,j) = 1 
          end if
          end if
  end do
end do

print("Converting ICE floats to doubles")
ice_fv_dbl = tofloat(ice_fv)
print("Copying ICE metadata")
copy_VarMeta(ice_fv,ice_fv_dbl)
delete(ice_fv)
print("Correcting ICE time records")
ice_fv_dbl_time = new((/12,numfvlat,numfvlon/),float)

if (do_anom) then

  climo_ice = in->ice_cov
  climo_ice_daily = linint1_n_Wrap (climo_ice&time, climo_ice, True, daily_time, 0, 0)
  
  ice_anom = ice_fv_dbl
  ice_anom = ice_fv_dbl - climo_ice_daily(doyint,:,:)
  
  do i = 0,11
    ice_fv_dbl_time(i,:,:) = climo_ice(i,:,:) + ice_anom(:,:)
  end do

else

  do i = 0,11
    ice_fv_dbl_time(i,:,:) = ice_fv_dbl
  end do

end if ; if do_anom

; clip ice cover frac to values between 0 and 1
ice_fv_dbl_time = where(ice_fv_dbl_time .gt. 1.0,1.0,ice_fv_dbl_time)
ice_fv_dbl_time = where(ice_fv_dbl_time .lt. 0.0,0.0,ice_fv_dbl_time)

; change units from yr 0000 to 0001
nc_time = in->time
nc_time@units="days since 0001-01-01 00:00:00"

date = new(1,integer)
datasec = new(1,integer)
timecoord = new(1,double)

date = 116
datasec = 43200
timecoord = 15.5

date!0 = "time"
datasec!0 = "time"
timecoord!0 = "time"

; Adding relevant coordinates
sst_fv_dbl_time!0 = "time"
sst_fv_dbl_time!1 = "lat"
sst_fv_dbl_time!2 = "lon"
sst_fv_dbl_time&time = nc_time
sst_fv_dbl_time&lat = fvlat
sst_fv_dbl_time&lon = fvlon
ice_fv_dbl_time!0 = "time"
ice_fv_dbl_time!1 = "lat"
ice_fv_dbl_time!2 = "lon"
ice_fv_dbl_time&time = nc_time
ice_fv_dbl_time&lat = fvlat
ice_fv_dbl_time&lon = fvlon

; Write to NetCDF
print("Writing file...")

if (copy_orig) then
  ; write the sst and ice streams to be read in CIME
  in->SST_cpl = (/sst_fv_dbl_time/)
  in->ice_cov = (/ice_fv_dbl_time/)
  ; print snapshot of interpolated data from obs
  in->sst_snapshot = sst_fv_dbl
  in->ice_snapshot = ice_fv_dbl
else
  system("/bin/rm -f -v "+SST_write_file)
  out = addfile(SST_write_file,"c")
       
  fAtt               = True
  fAtt@creation_date = systemfunc ("date")
  fileattdef(out, fAtt)
  filedimdef(out,"time",-1,True) 
  
  ; copy over some vars. Not sure top 2 are needed.
  out->date = in->date
  out->datesec = in->datesec
  out->lat = in->lat
  out->lon = in->lon
  
  ; updated time coord with new units
  out->time = nc_time

  ; output SST and ice fields
  out->SST_cpl = sst_fv_dbl_time
  out->ice_cov = ice_fv_dbl_time
  ; print snapshot of interpolated data from obs
  out->sst_snapshot = sst_fv_dbl
  out->ice_snapshot = ice_fv_dbl

end if

; for debugging 
;system("/bin/rm -f simple.nc")   ; remove any pre-existing file
;ncdf = addfile("simple.nc" ,"c")  ; open output netCDF file
;filedimdef(ncdf,"time",-1,True) 
;ncdf->climo_SST_daily = climo_SST_daily
;ncdf->sst_snapshot = sst_fv_dbl
;ncdf->ice_snapshot = ice_fv_dbl

print("done")

status_exit(9)

end
