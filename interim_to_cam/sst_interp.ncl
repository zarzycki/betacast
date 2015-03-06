; INPUTS
; dataSource (string) - valid values are GDAS, NOAAOI, ERAI
; YYYYMMDDHH (string or integer) - date string (only currently used for NOAAOI data)
; sstDataFile (string) - full path to data file with SST
; iceDataFile (string) - full path to data file with ice


load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"

begin

;ncl sst_interp.ncl YYYYMMDDHH=2005082412 'dataSource="ERAI"' 'sstDataFile = "/home/zarzycki/getECMWFdata/ERA-Int_sfc_2005082412.nc"' 'iceDataFile = "/home/zarzycki/getECMWFdata/ERA-Int_sfc_2005082412.nc"'
;ncl sst_interp.ncl YYYYMMDDHH=2012082400 'dataSource="NOAAOI"'

if (.not. isvar("dataSource"))
  print("Need to specify data source... exiting...")
  exit()
end if

if (dataSource .eq. "GDAS" ) then
  do_ice = False
else if (dataSource .eq. "NOAAOI" .or. dataSource .eq. "ERAI") then
  do_ice = True
else
  print("Unsupported data source")
  exit()
end if
end if

if (dataSource .eq. "NOAAOI" .and. .not. isvar("sstDataFile"))
  sstDataFile = "/home/zarzycki/sst.day.mean.2012.v2.nc"
  print("No user sstDataFile so specified a NOAAOI sst file: "+sstDataFile)
end if
if (dataSource .eq. "NOAAOI" .and. .not. isvar("iceDataFile"))
  iceDataFile = "/home/zarzycki/icec.day.mean.2012.v2.nc"
  print("No user iceDataFile so specified a NOAAOI ice file: "+iceDataFile)
end if

if (.not. isvar("sstDataFile")) then
  print("Need to provide sst data file... exiting...")
  exit()
end if

if (.not. isvar("YYYYMMDDHH"))
  if (dataSource .eq. "NOAAOI") then
    print("We need YYYYMMDDHH if using "+dataSource+" data")
    exit()
  else
    print("No YYYYMMDDHH specified")
    print("We are adding today's date with the 00Z cycle")
    YYYYMMDDHH=systemfunc("date -u +%Y%m%d"+"00")
  end if
end if

if (.not. isvar("iceDataFile"))
  if (do_ice) then
    print("No ice file specified but do_ice = true... exiting...")
    exit()
  else
    print("No ice file specified. do_ice = false, so empty stringing it")
    iceDataFile=""
  end if
end if

CAMSSTdir = "/home/zarzycki/sewx/SST/" ; directory where 1x1 CAM file is stored
;dataSource = "ERAI" ; GDAS, NOAAOI, ERAI

;==========================================================================

; Get YYYYMMDDHH from the initdate string
dtime_map = (/4,2,2,2/)
tst_dtime = str_split_by_length(tostring(YYYYMMDDHH), dtime_map)
yyyy=toint(tst_dtime(0))
mm=toint(tst_dtime(1))
dd=toint(tst_dtime(2))
hh=toint(tst_dtime(3))

; Load SST files that we want to overwrite
do filenum = 0,0 ;Begin multiple file loop
if (filenum .eq. 0) then
    print("Writing 1x1 SSTs")
    in = addfile(CAMSSTdir+"sst_1x1.nc","w")
else if (filenum .eq. 1) then
    print("Writing 0.23x0.31 SSTs")
    in = addfile(CAMSSTdir+"sst_0.23x0.31.nc","w")
else
end if
end if

if (dataSource .eq. "GDAS") then
  sst_file = addfile(sstDataFile,"r")
  sstlat = sst_file->lat_0
  sstlon = sst_file->lon_0
  sstFromData = sst_file->TMP_P0_L1_GLL0(:,:)
else if (dataSource .eq. "NOAAOI") then
  sst_file = addfile(sstDataFile,"r")
  ice_file = addfile(iceDataFile,"r")
  time = sst_file->time
  date = cd_inv_calendar(yyyy,mm,dd,0,0,0,time@units,0)
  sstFromData = sst_file->sst({date},:,:)*.01
  iceFromData = ice_file->icec({date},:,:)*.01
  sstlat = sst_file->lat
  sstlon = sst_file->lon
else if (dataSource .eq. "ERAI") then
  sst_file = addfile(sstDataFile,"r")
  ice_file = addfile(sstDataFile,"r")
  sstFromData = short2flt(sst_file->sst(0,:,:))
  iceFromData = short2flt(ice_file->ci(0,:,:))
  sstlat = sst_file->latitude
  sstlon = sst_file->longitude
else
end if
end if
end if

destLat  = in->lat
destLon  = in->lon
numdestLat = dimsizes(destLat)
numdestLon = dimsizes(destLon)

; wks   = gsn_open_wks ("x11", "coneff")       ; open ps file
; gsn_define_colormap(wks,"gui_default")        ; choose color map
; res                 = True                    ; plot mods desired
; res@gsnDraw         = False                   ; don't draw yet
; res@gsnFrame        = False                   ; don't advance frame yet
; res@cnFillOn        = True                    ; turn on color
; res@gsnSpreadColors = True                    ; use full colormap
; plot = gsn_csm_hov(wks, sstFromData, res) 
; draw(plot)
; frame(wks)

if (dataSource .eq. "GDAS" .or. dataSource .eq. "ERAI") then
  ; Need to flip raw GDAS lats/lons to match SST file
  sstDest = linint2_Wrap (sstlon,sstlat(::-1),sstFromData(::-1,:),True,destLon,destLat,0)
  if (do_ice) then  
    iceDest = linint2_Wrap (sstlon,sstlat(::-1),iceFromData(::-1,:),True,destLon,destLat,0)
  end if
else if (dataSource .eq. "NOAAOI") then
  sstDest = linint2_Wrap(sstlon,sstlat,sstFromData,True,destLon,destLat,0)
  if (do_ice) then
    iceDest = linint2_Wrap(sstlon,sstlat,iceFromData,True,destLon,destLat,0)
  end if
else
end if
end if

sstDest = linmsg(sstDest,-1)
sstDest = linmsg_n(sstDest,-1,0)

print("Sorting bad SST")
do i = 0,dimsizes(destLat)-1
    do j = 0,dimsizes(destLon)-1
            if (ismissing(sstDest(i,j))) then
                print("Found missing SST")
                sstDest(i,j) = 271.5
            else if (sstDest(i,j) .gt. 500) then
                print("Found bad value")
                sstDest(i,j) = 271.5
            end if
            end if
        sstDest(i,j) = sstDest(i,j) - 273.15
    end do
end do

print("Converting SST floats to doubles")
sstDest_dbl = tofloat(sstDest)
print("Copying SST metadata")
copy_VarMeta(sstDest,sstDest_dbl)
delete(sstDest)
print("Correcting SST time records")
; Add time record
sstDest_dbl_time = new((/12,numdestLat,numdestLon/),float)

do i = 0,11
    sstDest_dbl_time(i,:,:) = sstDest_dbl
end do


if (do_ice) then
  ;iceDest = linmsg(iceDest,-1)
  ;iceDest = linmsg_n(iceDest,-1,0)
  
  print("Sorting bad ice")
  do i = 0,dimsizes(destLat)-1
      do j = 0,dimsizes(destLon)-1
              if (ismissing(iceDest(i,j))) then
                  ;print("Found missing ice")
                  iceDest(i,j) = 0
              else if (iceDest(i,j) .gt. 500) then
                  print("Found bad ice")
                  iceDest(i,j) = 1 
              end if
              end if
      end do
  end do
  
  print("Converting ICE floats to doubles")
  iceDest_dbl = tofloat(iceDest)
  print("Copying ICE metadata")
  copy_VarMeta(iceDest,iceDest_dbl)
  delete(iceDest)
  print("Correcting ICE time records")
  iceDest_dbl_time = new((/12,numdestLat,numdestLon/),float)
  
  do i = 0,11
    iceDest_dbl_time(i,:,:) = iceDest_dbl
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
  
  in->SST_cpl = (/sstDest_dbl_time/)
  if (do_ice) then
    in->ice_cov = (/iceDest_dbl_time/)
  end if

  print("done")
  
delete(destLat)
delete(destLon)
delete(sstDest_dbl)
if (do_ice) then
  delete(iceDest_dbl)
  delete(iceDest_dbl_time)
end if
delete(sstDest_dbl_time)
delete(numdestLat)
delete(numdestLon)

end do ; End multiple file loop

end
