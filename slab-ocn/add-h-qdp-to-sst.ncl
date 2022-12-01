begin

;sfile="/glade/u/home/zarzycki/work/cesm2_2_0/components/cam/tools/icesst/bcgen/final_sst_cice.mon.mean.COBE2.EOFnn.nc"
WHICH_H = "GODAS"   ; Z16 or GODAS
WHICH_Q = "Z16"   ; Z16 or NONE
copy_ice_cov_if_avail=True   ; do we want to copy ice if it's on the input?

; get fields from SST stream
s = addfile(sfile,"r")
SST_cpl = s->SST_cpl
SST_cpl_dims = dimsizes(SST_cpl)
ntim = SST_cpl_dims(0)
nlat = SST_cpl_dims(1)
nlon = SST_cpl_dims(2)
lat = SST_cpl&lat
lon = SST_cpl&lon

; if we want to copy ice and it's available, let's load and then populate to new file
if (copy_ice_cov_if_avail .and. isfilevar(s, "ice_cov")) then
  copy_ice = True
else
  copy_ice = False
end if

if (copy_ice) then
  ice_cov = s->ice_cov
end if

; get boundary layer heights
if (WHICH_H .eq. "Z16") then
  afile="./nc_files/hq_z16.nc"
  a = addfile(afile,"r")
  hblt = a->hblt
else if (WHICH_H .eq. "GODAS") then
  afile="./nc_files/dbss_obml.mon.ltm.1991-2020.nc"
  a = addfile(afile,"r")
  hblt_o = a->dbss_obml
  hblt = linint2_Wrap (hblt_o&lon,hblt_o&lat,hblt_o, True, lon,lat, 0)
else
  print("unsupported WHICH_H")
  exit
end if
end if

; use existing correction or set to zero
if (WHICH_Q .eq. "Z16") then
  qfile="./nc_files/hq_z16.nc"
  q = addfile(qfile,"r")
  qdp = q->qdp
else if (WHICH_Q .eq. "NONE") then
  qdp = hblt
  qdp = 0.0
else
  print("unsupported WHICH_H")
  exit
end if
end if

nyears = ntim / 12

hblt_full = SST_cpl
qdp_full  = SST_cpl

do ii = 0,nyears-1
  ix = ii*12
  hblt_full(ix:ix+11,:,:) = (/hblt/)
  qdp_full(ix:ix+11,:,:)  = (/qdp/)
end do

;SST_cpl(1638,{-30:30},:) = SST_cpl(1638,{-30:30},:) - 10.
;SST_cpl(1639,{-30:30},:) = SST_cpl(1639,{-30:30},:) - 10.
;SST_cpl(1640,{-30:30},:) = SST_cpl(1640,{-30:30},:) - 10.

;hblt_full(1638,{0:90},:) = 10000000.
;hblt_full(1638,{-90:0},:) = 0.1

;qdp_full(1638,{0:80},:) = 15.
;qdp_full(1639,{0:80},:) = 15.
;qdp_full(1640,{0:80},:) = 15.

; first interpolate along longitude and then latitude to fill in missing values.
hblt_full = linmsg_n(hblt_full,-1,2)
hblt_full = linmsg_n(hblt_full,-1,1)
hblt_full = where(hblt_full .lt. 5.0,5.0,hblt_full)

times=SST_cpl&time
print(cd_calendar(times,3))

system("/bin/rm -f simple.nc")   ; remove any pre-existing file
ncdf = addfile("simple.nc" ,"c")  ; open output netCDF file

fAtt               = True            ; assign file attributes
fAtt@title         = "NCL Simple Approach to netCDF Creation"
fAtt@source_file   =  "original-file.nc"
fAtt@Conventions   = "None"
fAtt@creation_date = systemfunc ("date")
fileattdef( ncdf, fAtt )            ; copy file attributes

filedimdef(ncdf,"time",-1,True) 

if (copy_ice) then
  ncdf->ice_cov=ice_cov
end if
ncdf->SST_cpl=SST_cpl
ncdf->hblt=hblt_full
ncdf->qdp=qdp_full

end
