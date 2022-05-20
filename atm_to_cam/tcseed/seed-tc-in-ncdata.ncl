load "./fcns_seeding.ncl"

begin

;2004 8 31 0   19.772 281.867     1002.52  104.73
;2004 9 1 0   20.673 273.8     994.01  76.27
;(1)	1993 8 2 0   24.988 294.532     1004.89  121.38
; 11, 325, 991

minp = 994.01
target_rmw = 76270
cen_lat    = 24.75
cen_lon    = 360.-85.75
invert_vortex = False

deg_bnd = 12.
doplot=False

origfile="sample_ne120.nc"
seedfile="sample_ne120_unseed.nc"
system("rm -v "+seedfile)
system("cp -v "+origfile+" "+seedfile)

inputFile = addfile(seedfile,"w")

lat = inputFile->lat
lon = inputFile->lon
lev = inputFile->lev
ncol = dimsizes(lat)

hyai = inputFile->hyai
hybi = inputFile->hyai
hyam = inputFile->hyam
hybm = inputFile->hybm
P0   = 100000.

u = inputFile->U
v = inputFile->V
ps = inputFile->PS
t = inputFile->T
q = inputFile->Q

; make a copy of orig vars on file
inputFile->U_orig = u
inputFile->V_orig = v
inputFile->PS_orig = ps
inputFile->T_orig = t
inputFile->Q_orig = ps

; find model point of cen_lat, cen_lon
gcdist = gc_latlon(cen_lat,cen_lon,lat,lon,2,2)
minix = minind(gcdist)

;======= get dp and rp

if (.not. invert_vortex) then
  ambps = ps(0,minix)
  print("ambient ps: "+ambps)
  dp = ambps - minp*100.
  rp = get_rp_from_dp_rmw(cen_lat,dp,target_rmw)
else
  ; for inversion, we need to estimate...
  ; dp and target_rmw -> get_stuff(ps,windbot,initiallat,initiallon)
  ; dp should be difference between vortex min and environment at ring away from storm
  ; rmw should be radius of maximum wind scaled?
  dp = 3500.
  target_rmw=200000.
  cen_lat    = 24.75
  cen_lon    = 360.-85.75
  rp = get_rp_from_dp_rmw(cen_lat,dp,target_rmw)
  ;rp = 82031.2
end if

print("dp: "+dp)
print("rp: "+rp)

;======= seed TC

do ii = 0,ncol-1

  if (mod(ii,1000) .eq. 0) then
    print("At ncol: "+ii+" of "+ncol)
  end if

  if ( gcdist(ii) .le. deg_bnd ) then
  
    do kk = 0,dimsizes(lev)-1
      ; calc pressure level
      p = hyam(kk)*P0 + hybm(kk)*ps(0,ii)
      ; get perturbed TC field
      theArr = tctestcase(cen_lon,cen_lat,dp,rp,lon(ii),lat(ii),p,-999,0,ps(0,ii),u(0,kk,ii),v(0,kk,ii),t(0,kk,ii),q(0,kk,ii),invert_vortex)
      ; convert and write to state arrays
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
    
  end if  ; if this cell is an edited cell
  
end do

; Write file!
inputFile->PS=ps
inputFile->U=u
inputFile->V=v
inputFile->T=t
inputFile->Q=q

; add metadata
inputFile@dp = dp
inputFile@target_rmw = target_rmw
inputFile@rp = rp

if (doplot) then
  wks = gsn_open_wks("x11","tc_mask")
  
  res                      = True
  res@gsnMaximize          = True
  res@gsnSpreadColors     = True 	 
  res@cnFillOn             = True              ; turn on color
  res@cnLinesOn            = False             ; turn off contour lines
  res@cnLineLabelsOn       = False
  res@lbLabelAutoStride   = True         ; Clean up labelbar labels.

  res@mpMinLatF = cen_lat - 15.
  res@mpMaxLatF = cen_lat + 15.
  res@mpMinLonF = cen_lon - 15.
  res@mpMaxLonF = cen_lon + 15.

  res@sfXArray            = lon        ; Required to tell NCL where to
  res@sfYArray            = lat        ; overlay data on globe.

  res@cnLevelSelectionMode = "ManualLevels"	; manually set the contour levels with the following 3 resources
  res@cnMinLevelValF  = 96000.			; set the minimum contour level
  res@cnMaxLevelValF  = 104000.			; set the maximum contour level
  res@cnLevelSpacingF = 500.			; set the interval between contours
  map = gsn_csm_contour_map_ce(wks,ps(0,:),res)
  
  res@cnLevelSelectionMode = "ManualLevels"	; manually set the contour levels with the following 3 resources
  res@cnMinLevelValF  = -20.
  res@cnMaxLevelValF  = 20.
  res@cnLevelSpacingF = 2.			; set the interval between contours
  map = gsn_csm_contour_map_ce(wks,u(0,26,:),res)
end if

end
