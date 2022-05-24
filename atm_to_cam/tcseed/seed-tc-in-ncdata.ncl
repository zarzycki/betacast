load "./fcns_seeding.ncl"

begin

invert_vortex = False    ; do we want to ADD a vortex (False) or REMOVE a vortex (True)
deg_bnd = 15.            ; what is the radius where the code will influence the analysis (GCD)
doplot=False             ; should we plot a diagnostic map
modify_q = True          ; do we also modify Q based on RH and the T anomaly?
origfile="/glade/u/home/zarzycki/work/sewx/INIC/sample_ne120_unseed.nc"
seedfile="/glade/u/home/zarzycki/work/sewx/INIC/sample_ne120_unseed2.nc"

if (.not. invert_vortex) then
  ; here we are seeding the vortex
  ; rp and dp will come from minp and target_rmw
  ; minp -> what is the SLP min we want to inject?
  ; target_rmw -> what is the RMW we want?
  ; setting zp, gamma, exppr, modify_q_mult to -1 uses defaults
  minp = 990.0
  target_rmw = 150000.
  cen_lat    = 22.75
  cen_lon    = 360.-92.0
  zp = -1.
  gamma_ = -1.
  exppr = -1.
  modify_q_mult = -1.
else
  ; here we are removing the vortex, these settings come from find-tc-fill-params
  cen_lat = 21.87653196344101
  cen_lon = 306.2441153488365
  rp       = 179761 ; Radius for calculation of PS
  dp       = 1641.97   ; Delta P for calculation of PS
  zp       = 14305.1   ; Height for calculation of P
  gamma_   = 0.0065    ; lapse rate
  exppr    = 1.12118     ; Exponent for r dependence of p
  modify_q_mult = 2.5
end if

wcSeedStrt     = systemfunc("date")

print("copying files...")
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
inputFile->Q_orig = q

; get gcdist from projected center lat/lon to full lat/lon arrays
gcdist = gc_latlon(cen_lat,cen_lon,lat,lon,2,2)

;======= get dp and rp
if (.not. invert_vortex) then
  ; figure out true model point of center lat
  ;tmps = ps(0,:)
  ;tmps = where(gcdist .gt. 5.0 , 10000000. , tmps)
  ;minix = minind(tmps)
  ;print(lat(minix)+" "+lon(minix))
  ;delete(tmps)
  
  minix = minind(gcdist)
  print(lat(minix)+" "+lon(minix))
  
  ambps = ps(0,minix)
  print("ambient ps: "+ambps)
  dp = ambps - minp*100.
  rp = get_rp_from_dp_rmw(cen_lat,dp,target_rmw)
  print("derived dp: "+dp)
  print("derived rp: "+rp)
end if

;======= seed TC

do ii = 0,ncol-1

  if (mod(ii,1000) .eq. 0) then
    print("At ncol: "+ii+" of "+ncol)
  end if

  ; only edit cells that are within our defined radius of TC seed
  if ( gcdist(ii) .le. deg_bnd ) then
  
    do kk = 0,dimsizes(lev)-1
      ; calc pressure level
      p = hyam(kk)*P0 + hybm(kk)*ps(0,ii)
      ; get perturbed TC field
      theArr = tctestcase(todouble(cen_lon),todouble(cen_lat),dp,rp,zp,exppr,gamma_,lon(ii),lat(ii),p,-999,0,ps(0,ii),u(0,kk,ii),v(0,kk,ii),t(0,kk,ii),q(0,kk,ii),invert_vortex,modify_q,modify_q_mult)
      ; convert and write to state arrays
      v(0,kk,ii) = totype(theArr(0), typeof(v))
      u(0,kk,ii) = totype(theArr(1), typeof(u))
      q(0,kk,ii) = totype(theArr(2), typeof(q))
      t(0,kk,ii) = totype(theArr(3), typeof(t))
    end do
    
    ; we only need to update ps once per column
    ps(0,ii) = totype(theArr(4), typeof(ps))
    if (ps(0,ii) .lt. 80000) then
      print("For diagnostic purposes, PS below 800mb --> "+ps(0,ii))
    end if
    
  end if  ; if this cell is an edited cell
  
end do  ; end column loop

; Write file!
inputFile->PS=ps
inputFile->U=u
inputFile->V=v
inputFile->T=t
inputFile->Q=q

; add metadata
;inputFile@invert_vortex = invert_vortex
;inputFile@dp = dp
;inputFile@target_rmw = target_rmw
;inputFile@rp = rp
;inputFile@zp = zp
;inputFile@gamma_ = gamma_
;inputFile@exppr = exppr
;inputFile@cen_lat = cen_lat
;inputFile@cen_lon = cen_lon
;inputFile@modify_q = modify_q
;inputFile@modify_q_mult = modify_q_mult

wallClockElapseTime(wcSeedStrt, "Time to seed", 0)

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
