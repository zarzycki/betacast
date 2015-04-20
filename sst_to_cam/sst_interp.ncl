load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"

begin

;initdate=2012082400
;SST_write_dir = "/glade/p/work/zarzycki/sewx/SST/"
;datasource = "GDAS"

print(initdate)
KtoC=273.15

if (datasource .eq. "GDAS") then
  do_ice = False
else if (datasource .eq. "NOAAOI") then
  do_ice = True
else
end if
end if

; Get YYYYMMDDHH from the initdate string
dtime_map = (/4,2,2,2/)
tst_dtime = str_split_by_length(initdate, dtime_map)
yyyy=toint(tst_dtime(0))
mm=toint(tst_dtime(1))
dd=toint(tst_dtime(2))
hh=toint(tst_dtime(3))

; Load SST files that we want to overwrite
do filenum = 0,0 ;Begin multiple file loop
if (filenum .eq. 0) then
    print("Writing 1x1 SSTs")
    in = addfile(SST_write_file,"w")
else if (filenum .eq. 1) then
    print("Writing 0.23x0.31 SSTs")
    in = addfile(SST_write_dir+"sst_0.23x0.31.nc","w")
else
end if
end if

if (datasource .eq. "GDAS") then
  sst_file = addfile(sstDataFile,"r")
  sstlat = sst_file->lat_0
  sstlon = sst_file->lon_0
  sst_gfs = sst_file->TMP_P0_L1_GLL0(:,:)
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

do i = 0,11
    sst_fv_dbl_time(i,:,:) = sst_fv_dbl
end do


if (do_ice) then
  ;ice_fv = linmsg(ice_fv,-1)
  ;ice_fv = linmsg_n(ice_fv,-1,0)
  
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
  
  do i = 0,11
    ice_fv_dbl_time(i,:,:) = ice_fv_dbl
  end do

end if

print("Convert to float")

  date = new(1,integer)
  datasec = new(1,integer)
  timecoord = new(1,double)
  
  date = 116
  datasec = 43200
  timecoord = 15.5
  
  date!0 = "time"
  datasec!0 = "time"
  timecoord!0 = "time"
  
; Write to NetCDF
  print("Writing file...")
  
  in->SST_cpl = (/sst_fv_dbl_time/)
  if (do_ice) then
    in->ice_cov = (/ice_fv_dbl_time/)
  end if

  print("done")
  
delete(fvlat)
delete(fvlon)
delete(sst_fv_dbl)
if (do_ice) then
  delete(ice_fv_dbl)
  delete(ice_fv_dbl_time)
end if
delete(sst_fv_dbl_time)
delete(numfvlat)
delete(numfvlon)

end do ; End multiple file loop

status_exit(9)

end
