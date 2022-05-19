load "./fcns_seeding.ncl"

begin

;2004 8 31 0   19.772 281.867     1002.52  104.73
;2004 9 1 0   20.673 273.8     994.01  76.27
;(1)	1993 8 2 0   24.988 294.532     1004.89  121.38
; 11, 325, 991

minp = 991.00 
target_rmw = 121380
cen_lat    = 11.0
cen_lon    = 325.0

lat_bnd = 12.
lon_bnd = 12.
doplot=False

system("rm ./testing.nc")
system("cp ./FHIST-ne30-betacast-ERA5-x001_INIC.nc ./testing.nc")

inputFile = addfile("./testing.nc","w")

lat = inputFile->lat
lon = inputFile->lon
lev = inputFile->lev
ncol = dimsizes(lat)

hyai = inputFile->hyai
hybi = inputFile->hyai
hyam = inputFile->hyam
hybm = inputFile->hybm
;P0   = inputFile->P0
P0   = 100000.

u = inputFile->U
v = inputFile->V
ps = inputFile->PS
t = inputFile->T
q = inputFile->Q

;to_delete=(/"_FillValue_original","missing_value_original"/)
;to_delete=(/"_FillValue","missing_value"/)
;to_delete=(/"maximum_value","minimum_value"/)
to_delete=(/"time"/)

delete_VarAtts(u,to_delete)
delete_VarAtts(v,to_delete)
delete_VarAtts(t,to_delete)
delete_VarAtts(q,to_delete)
delete_VarAtts(ps,to_delete)

; find model point of cen_lat, cen_lon
gcdist = gc_latlon(cen_lat,cen_lon,lat,lon,2,4)
minix = minind(gcdist)

;======= get dp

ambps = ps(0,minix)
dp = ambps - minp*100.

;======= get rp

rp = get_rp_from_dp_rmw(cen_lat,dp,target_rmw)

;======= seed TC

do ii = 0,ncol-1

  if (mod(ii,1000) .eq. 0) then
    print("At ncol: "+ii+" of "+ncol)
  end if

  if (lat(ii) .lt. (cen_lat + lat_bnd) .and. lat(ii) .gt. (cen_lat - lat_bnd) .and. \
    lon(ii) .lt. (cen_lon + lon_bnd) .and. lon(ii) .gt. (cen_lon - lon_bnd)) then

    do kk = 0,dimsizes(lev)-1

      ; calc pressure level
      p = hyam(kk)*P0 + hybm(kk)*ps(0,ii)

      theArr = tctestcase(cen_lon,cen_lat,dp,rp,lon(ii),lat(ii),p,-999,0,ps(0,ii),u(0,kk,ii),v(0,kk,ii),t(0,kk,ii),q(0,kk,ii))

      v(0,kk,ii) = totype(theArr(0), typeof(v))
      u(0,kk,ii) = totype(theArr(1), typeof(u))
      q(0,kk,ii) = totype(theArr(2), typeof(q))
      t(0,kk,ii) = totype(theArr(3), typeof(t))

    end do
    
    ; we only want to update ps once per column
    ps(0,ii) = totype(theArr(4), typeof(ps))
    if (ps(0,ii) .lt. 80000) then
      print(" "+ps(0,ii))
    end if

  end if
end do

; Write file!

inputFile->PS=ps
inputFile->U=u
inputFile->V=v
inputFile->T=t
inputFile->Q=q

if (doplot) then
  wks = gsn_open_wks("x11","tc_mask")
  
  res                      = True
  res@gsnMaximize          = True
  res@gsnSpreadColors     = True 	 
  res@cnFillOn             = True              ; turn on color
  res@cnLinesOn            = False             ; turn off contour lines
  res@cnLineLabelsOn       = False
  res@lbLabelAutoStride   = True         ; Clean up labelbar labels.

  res@cnLevelSelectionMode = "ManualLevels"	; manually set the contour levels with the following 3 resources
  res@cnMinLevelValF  = 96000.			; set the minimum contour level
  res@cnMaxLevelValF  = 102000.			; set the maximum contour level
  res@cnLevelSpacingF = 500.			; set the interval between contours
  
  res@mpMinLatF = cen_lat - 15.
  res@mpMaxLatF = cen_lat + 15.
  res@mpMinLonF = cen_lon - 15.
  res@mpMaxLonF = cen_lon + 15.

  res@sfXArray            = lon        ; Required to tell NCL where to
  res@sfYArray            = lat        ; overlay data on globe.
  map = gsn_csm_contour_map_ce(wks,ps(0,:),res)
end if

end
