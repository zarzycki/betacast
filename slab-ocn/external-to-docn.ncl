begin


f=addfile(infile,"r")

r=addfile("../sst_to_cam/sst_1x1.nc","r")

lat_in = f->lat(::-1)
lon_in = f->lon
time_in = f->time
lat_out = r->lat
lon_out = r->lon

var_in = f->sst(:,::-1,:)
ice_in = f->icec(:,::-1,:)

var_in = linmsg(var_in,-1)
var_in = linmsg_n(var_in,-1,1)

ice_in = linmsg(ice_in,-1)
ice_in = where(ismissing(ice_in),0.0,ice_in)

alldates = cd_calendar(time_in, 0)
yyyy = toint(alldates(:,0))
mm = toint(alldates(:,1))
dd = toint(alldates(:,2))
sssss = toint(alldates(:,5))

do ii = 0,dimsizes(time_in)-1
  if (mm(ii) .eq. 2) then
    dd(ii) = 15
  else
    dd(ii) = 16
  end if
  
  if (mm(ii) .eq. 2 .or. mm(ii) .eq. 4  .or. mm(ii) .eq. 6  .or. mm(ii) .eq. 9  .or. mm(ii) .eq. 11) then
    sssss(ii) = 0
  else
    sssss(ii) = 43200
  end if

end do

date = yyyy*10000 + mm*100 + dd
date!0 = "time"

datesec = sssss
datesec!0 = "time"

print(date+" "+datesec)

SST_write_file=infile+"_CIMEified.nc"
system("/bin/rm -f -v "+SST_write_file)
out = addfile(SST_write_file,"c")
     
fAtt               = True
fAtt@creation_date = systemfunc ("date")
fileattdef(out, fAtt)
filedimdef(out,"time",-1,True) 

; copy over some vars. Not sure top 2 are needed.
;out->date = in->date
;out->datesec = in->datesec
out->lat = lat_in
out->lon = lon_in

; updated time coord with new units
out->time = time_in
out->date = (/date/)
out->datesec = (/datesec/)


; output SST and ice fields
out->SST = var_in
out->ICEFRAC = ice_in







end