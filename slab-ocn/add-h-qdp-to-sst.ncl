begin

sfile="/global/homes/c/czarzyck/m2637/betacast/som/docn_files/sst_HadOIBl_bc_1x1_1850_2021_c120422.nc"
afile="./hq_z16.nc"

s = addfile(sfile,"r")
a = addfile(afile,"r")

; get fields from SST stream
SST_cpl = s->SST_cpl

; get 
hblt = a->hblt
qdp = a->qdp

SST_cpl_dims = dimsizes(SST_cpl)
ntim = SST_cpl_dims(0)
nlat = SST_cpl_dims(1)
nlon = SST_cpl_dims(2)

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

;hblt_full = where(ismissing(hblt_full),1000000.,hblt_full)
hblt_full = linmsg_n(hblt_full,-1,2)
hblt_full = where(hblt_full .lt. 10.0,10.0,hblt_full)

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

ncdf->SST_cpl=SST_cpl
ncdf->hblt=hblt_full
ncdf->qdp=qdp_full

end
