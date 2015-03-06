load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"

begin

do filenum = 0,0 ;Begin multiple file loop

SST_write_dir = "/home/zarzycki/SE_weather/SST/"

;filenum = 0
; Select SST file (native grid for FV, 1x1 grid for CESM
if (filenum .eq. 0) then
    print("Writing 1x1 SSTs")
;    fV = addfile ("/glade/proj3/cseg/inputdata/atm/cam/sst/sst_HadOIBl_bc_1x1_clim_c101029.nc", "r")
    in = addfile(SST_write_dir+"sst_1x1.nc","w")
else if (filenum .eq. 1) then
    print("Writing 0.23x0.31 SSTs")
;    fV = addfile ("/glade/proj3/cseg/inputdata/atm/cam/sst/sst_HadOIBl_bc_0.23x0.31_clim_c110526.nc", "r")
    in = addfile(SST_write_dir+"sst_0.23x0.31.nc","w")
else
end if
end if

datestr = systemfunc("date +'%Y%m%d'")		; Get date string for parsing

rawhour = systemfunc("date -u +'%H'")

sst_file_dir="/home/zarzycki/SE_weather/GFS/"
sst_file_name="gdas1.t00z.sstgrb.grib2"
sst_file = addfile(sst_file_dir+sst_file_name,"r")



sstlat = sst_file->lat_0
sstlon = sst_file->lon_0
sst_gfs = sst_file->TMP_P0_L1_GLL0(:,:)

fvlat  = in->lat
fvlon  = in->lon

sst_fv = linint2_Wrap (sstlon,sstlat(::-1),sst_gfs(::-1,:),True,fvlon,fvlat,0)
;ice_fv = linint2_Wrap (sstlon,sstlat(::-1),ice_gfs(::-1,:),True,fvlon,fvlat,0)

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
        sst_fv(i,j) = sst_fv(i,j) - 273.15
    end do
end do

;ice_fv = linmsg(ice_fv,-1)
;ice_fv = linmsg_n(ice_fv,-1,0)

;print("Sorting bad ice")
;do i = 0,dimsizes(fvlat)-1
;    do j = 0,dimsizes(fvlon)-1
;            if (ismissing(ice_fv(i,j))) then
;                print("Found missing ice")
;                ice_fv(i,j) = 1
;            else if (ice_fv(i,j) .gt. 500) then
;                print("Found bad ice")
;                ice_fv(i,j) = 1 
;            end if
;            end if
;    end do
;end do

print("Converting floats to doubles")
sst_fv_dbl = tofloat(sst_fv)
;ice_fv_dbl = tofloat(ice_fv)

print("Copying metadata")
copy_VarMeta(sst_fv,sst_fv_dbl)
;copy_VarMeta(ice_fv,ice_fv_dbl)

; clean up some stuff
delete(sst_fv)
;delete(ice_fv)

numfvlat = dimsizes(fvlat)
numfvlon = dimsizes(fvlon)

print("Correcting time records")
; Add time record

sst_fv_dbl_time = new((/12,numfvlat,numfvlon/),float)
;ice_fv_dbl_time = new((/12,numfvlat,numfvlon/),float)

do i = 0,11
    ;ice_fv_dbl_time(i,:,:) = ice_fv_dbl
    sst_fv_dbl_time(i,:,:) = sst_fv_dbl
end do

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
  ;in->ice_cov = (/ice_fv_dbl_time/)

  
  
  print("done")
  
delete(fvlat)
delete(fvlon)
delete(sst_fv_dbl)
;delete(ice_fv_dbl)
;delete(ice_fv_dbl_time)
delete(sst_fv_dbl_time)
delete(numfvlat)
delete(numfvlon)

end do ; End multiple file loop

end
