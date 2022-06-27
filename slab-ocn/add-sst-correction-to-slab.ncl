load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"

begin

SST_write_name="climo_1x1_slab_exp009.nc"
SST_write_dir = "/glade/u/home/zarzycki/work/slab/docn_files/"
in = addfile(SST_write_dir+"/"+SST_write_name,"w")

corr = in->qdp
lat = in->yc
lon = in->xc

printVarSummary(corr)

SST_bias_dir="/glade/u/home/zarzycki/work/slab/bias_files/slab_exp009/"
monthstr=(/"01","02","03","04","05","06","07","08","09","10","11","12"/)

do i = 0,11
corrFile=SST_bias_dir+"/diff"+monthstr(i)+".nc"
if (isfilepresent(corrFile))
  print("Found correction from month: "+monthstr(i))
  f = addfile(corrFile,"r")
  corrLat = f->lat
  corrLon = f->lon
  corrVar = f->SST(0,:,:)
  ;print(corrLon)
  corrOut = linint2_Wrap (corrLon,corrLat,corrVar,True,lon,lat,0)
  printVarSummary(corrOut)
  corr(i,:,:) = corrOut
else
  corr(i,:,:) = 0.0
end if
  
end do



regrid =False
;-- if (regrid) then
;--   fvlat  = in->yc
;--   fvlon  = in->xc
;--   numfvlat = dimsizes(fvlat)
;--   numfvlon = dimsizes(fvlon)
;-- 
;--     sst_fv = linint2_Wrap(sstlon,sstlat,sst_gfs,True,fvlon,fvlat,0)
;--     ice_fv = linint2_Wrap(sstlon,sstlat,ice_gfs,True,fvlon,fvlat,0)
;-- 
;-- sst_fv = linmsg(sst_fv,-1)
;-- sst_fv = linmsg_n(sst_fv,-1,0)
;-- 
;-- print("Sorting bad SST")
;-- do i = 0,dimsizes(fvlat)-1
;--     do j = 0,dimsizes(fvlon)-1
;--             if (ismissing(sst_fv(i,j))) then
;--                 print("Found missing SST")
;--                 sst_fv(i,j) = 271.5
;--             else if (sst_fv(i,j) .gt. 500) then
;--                 print("Found bad value")
;--                 sst_fv(i,j) = 271.5
;--             end if
;--             end if
;--         sst_fv(i,j) = sst_fv(i,j); - 273.15
;--     end do
;-- end do
;-- 
;-- print("Converting SST floats to doubles")
;-- sst_fv_dbl = tofloat(sst_fv)
;-- print("Copying SST metadata")
;-- copy_VarMeta(sst_fv,sst_fv_dbl)
;-- delete(sst_fv)
;-- print("Correcting SST time records")
;-- ; Add time record
;-- sst_fv_dbl_time = new((/12,numfvlat,numfvlon/),float)
;-- 
;-- do i = 0,11
;--     sst_fv_dbl_time(i,:,:) = sst_fv_dbl
;-- end do
;-- 
;-- 
;-- end if


; Write to NetCDF
  print("Writing file...")
  in->qdp = (/corr/)
  

end
