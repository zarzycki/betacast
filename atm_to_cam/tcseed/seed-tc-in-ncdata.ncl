if ismissing(getenv("BETACAST")) then
  print("We are running local-only atm_to_cam. Set export BETACAST env to run elsewhere")
  ; if we don't have BETACAST, check to see where are located...
  ; we cannot do relative paths
  this_script_name=get_script_name()
  slash_index=str_index_of_substr(this_script_name, "/", -1)
  if (.not. ismissing(slash_index)) then
    print("The script is called from another directory")
    print("I found: "+this_script_name+" which contains a /")
    print("This isn't supported! Either set BETACAST as env param or run from local dir!")
    exit
  else
    print("The script is called from the local directory")
    PATHTOHERE=systemfunc("pwd")
  end if
  ; get BETACAST based on current dir assuming one up
  BETACAST = systemfunc("echo $(dirname $(pwd))")
else
  print("Not local only!")
  PATHTOHERE=getenv("BETACAST")+"/atm_to_cam/tcseed/"
  BETACAST=getenv("BETACAST")
end if

print(""+BETACAST)
print(""+PATHTOHERE)

loadscript(BETACAST+"/functions/fcns_seeding.ncl")
loadscript(BETACAST+"/functions/namelist.ncl")

begin

;debug restart file vortex add
;system("cp /glade/u/home/zarzycki/scratch/APE3-test/run/ORIG.nc "+seedfile)

;pthi="../../namelists/unseed.default.nl"
invert_vortex  = keyword_values( pthi, "invert_vortex", "logical")
;seedfile = keyword_values( pthi, "new_file", "string")
modify_q = keyword_values( pthi, "modify_q", "logical")
if (modify_q) then
  modify_q_mult  = keyword_values( pthi, "modify_q_mult", "float")
else
  modify_q_mult=1.0
end if
gamma_  = keyword_values( pthi, "gamma", "float")
cen_lat  = keyword_values( pthi, "psminlat", "float")
cen_lon  = keyword_values( pthi, "psminlon", "float")
zp  = keyword_values( pthi, "zp", "float")
exppr  = keyword_values( pthi, "exppr", "float")
restart_file  = keyword_values( pthi, "restart_file", "logical")

deg_bnd=15.      ; what is the radius where the code will influence the analysis (GCD)
doplot=False     ; should we plot a diagnostic map
;restart_file=True

if (.not. invert_vortex) then
  ; here we are seeding the vortex
  ; rp and dp will come from minp and target_rmw
  ; minp -> what is the SLP min we want to inject?
  ; target_rmw -> what is the RMW we want?
  ; setting zp, gamma, exppr, modify_q_mult to -1 uses defaults
  minp       = keyword_values( pthi, "minp", "float")
  target_rmw = keyword_values( pthi, "target_rmw", "float")

  print(minp)
  if (minp .lt. -19999.0) then
    minp = 995.0
    print("SEED TC: setting minp to default of "+minp)
  else if (minp .ge. -19999.0 .and. minp .le. 0.0) then
    print("SEED TC: Using the negative value of minp as dp: "+minp)
  else
    print("SEED TC: using minp from input file as: "+minp)
  end if
  end if

  if (target_rmw .lt. 0.0) then
    target_rmw = 200000.
    print("SEED TC: setting target_rmw to default of "+target_rmw)
  end if
else
  ; here we are removing the vortex, these settings come from find-tc-fill-params
  rp  = keyword_values( pthi, "rp", "float")
  dp  = keyword_values( pthi, "dp", "float")
end if

wcSeedStrt     = systemfunc("date")

;print("copying files...")
;system("rm -v "+seedfile)
;system("cp -v "+origfile+" "+seedfile)

inputFile = addfile(seedfile,"w")

if (restart_file) then
  if (isfilevar(inputFile, "lat_d")) then
    lat = inputFile->lat_d
    lon = inputFile->lon_d
    area = inputFile->area_d
  else
    if (isfilevar(inputFile, "lat")) then
      lat = inputFile->lat
      lon = inputFile->lon
      area = inputFile->area
    else
      print("Cannot find lat or lat_d on input file")
      exit
    end if
  end if
else
  lat = inputFile->lat
  lon = inputFile->lon
end if

hyai = inputFile->hyai
hybi = inputFile->hybi
hyam = inputFile->hyam
hybm = inputFile->hybm
P0   = 100000.

ncol = dimsizes(lat)
nlev = dimsizes(hybm)
nlevi = dimsizes(hybi)

if (restart_file) then
  ; based on CESM2 dry mass coordinate *but* not using physgrid
  ; physgrid requires dp to be updated on the fvm grid as well...
  u = inputFile->U
  v = inputFile->V
  ps = inputFile->PSDRY
  t = inputFile->T
  dpq = inputFile->dpQ   ; this is the pressure contribution from q to PS

  print("getting dp3d")
  ; uses same formulation as in dp_coupling
  dp3d = t
  dp3d = 0.0
  do ii = 0,ncol-1
    dp3d(0,:,ii) = ((hyai(1:nlevi-1) - hyai(0:nlevi-2)) * P0) + ((hybi(1:nlevi-1) - hybi(0:nlevi-2)) * ps(0,ii))
  end do
  print("min/max dp3d: "+min(dp3d)+" "+max(dp3d))

  print("getting q from dpq + dp3d")
  q = dpq
  q = dpq / dp3d
  print("min/max q: "+min(q)+" "+max(q))

  ; calculate before dry mass
  avg_ps_in = sum(ps(0,:) * area) / sum(area)
  print("Average PSDRY in: "+avg_ps_in)

  ; dry to total PS
  corrMass = ps
  corrMass = 0.0
  corrMass = corrMass + dim_sum_n(inputFile->dpQ,1)
  corrMass = corrMass + dim_sum_n(inputFile->dpCLDICE,1)
  corrMass = corrMass + dim_sum_n(inputFile->dpCLDLIQ,1)
  print("Initial wet pressure for mass correction: "+avg(corrMass))
  ps = ps + corrMass
  print("Total PS after adding condensates: "+sum(ps(0,:) * area) / sum(area))
  ; NOTE: now ps is total pressure!
else
  u = inputFile->U
  v = inputFile->V
  ps = inputFile->PS
  t = inputFile->T
  q = inputFile->Q
end if

u_orig = u
v_orig = v
ps_orig = ps
t_orig = t
q_orig = q

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
  if (minp .gt. 0.0)
    dp = ambps - minp*100.
  else
    dp = -minp
  end if
  rp = get_rp_from_dp_rmw(cen_lat,dp,target_rmw)
  print("derived dp: "+dp)
  print("derived rp: "+rp)
end if

;======= seed TC
; print some diagnostics
print("**************")
print("invert_vortex: "+invert_vortex)
print("cen_lat: "+cen_lat)
print("cen_lon: "+cen_lon)
print("dp: "+dp)
print("rp: "+rp)
print("zp: "+zp)
print("gamma: "+gamma_)
print("exppr: "+exppr)
print("modify_q: "+modify_q)
print("modify_q_mult: "+modify_q_mult)
print("**************")

do ii = 0,ncol-1

  if (mod(ii,1000) .eq. 0) then
    print("At ncol: "+ii+" of "+ncol)
  end if

  ; only edit cells that are within our defined radius of TC seed
  if ( gcdist(ii) .le. deg_bnd ) then

    do kk = 0,nlev-1
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
if (restart_file) then

  ; update dpq by assuming that the percent change in q due to vortex adding equals
  ; the percent change in dpq
  q_perc = q / q_orig
  print("qperc min/max: "+min(q_perc)+" "+max(q_perc))
  dpq = q_perc * dpq

  ; recalculate the dry PS from total PS based on new dpq and unmod'ed condensates
  print("Total PS after adding vortex: "+sum(ps(0,:) * area) / sum(area))
  corrMass = 0.0
  corrMass = corrMass + dim_sum_n(dpq,1)
  corrMass = corrMass + dim_sum_n(inputFile->dpCLDICE,1)
  corrMass = corrMass + dim_sum_n(inputFile->dpCLDLIQ,1)
  print("New wet pressure for mass correction: "+avg(corrMass))
  ps = ps - corrMass

  ; calculate new global dry mass and add/substract PS to ensure PSDRY before = PSDRY after
  avg_ps_out = sum(ps(0,:) * area) / sum(area)
  print("Average PSDRY out: "+avg_ps_out)
  ps = ps - (avg_ps_out - avg_ps_in)
  avg_ps_final = sum(ps(0,:) * area) / sum(area)
  print("Final PSDRY out: "+avg_ps_final)

  ; write to file
  inputFile->PSDRY=ps
  inputFile->U=u
  inputFile->V=v
  inputFile->T=t
  inputFile->dpQ=dpq
else
  inputFile->PS=ps
  inputFile->U=u
  inputFile->V=v
  inputFile->T=t
  inputFile->Q=q
end if

if (.not. restart_file) then
  ; make a copy of orig vars on file
  if (.not. isfilevar(inputFile, "U_orig")) then
    print("copying u -> u_orig")
    inputFile->U_orig = u_orig
    inputFile->V_orig = v_orig
    inputFile->PS_orig = ps_orig
    inputFile->T_orig = t_orig
    inputFile->Q_orig = q_orig
  else
    if (.not. isfilevar(inputFile, "U_orig2")) then
      print("copying u -> u_orig2")
      inputFile->U_orig2 = u_orig
      inputFile->V_orig2 = v_orig
      inputFile->PS_orig2 = ps_orig
      inputFile->T_orig2 = t_orig
      inputFile->Q_orig2 = q_orig
    else
      print("We already have U_orig + U_orig2 on file, not adding any more!")
    end if
  end if
end if

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

  res@cnLevelSelectionMode = "ManualLevels"  ; manually set the contour levels with the following 3 resources
  res@cnMinLevelValF  = 96000.      ; set the minimum contour level
  res@cnMaxLevelValF  = 104000.      ; set the maximum contour level
  res@cnLevelSpacingF = 500.      ; set the interval between contours
  map = gsn_csm_contour_map_ce(wks,ps(0,:),res)

  res@cnLevelSelectionMode = "ManualLevels"  ; manually set the contour levels with the following 3 resources
  res@cnMinLevelValF  = -20.
  res@cnMaxLevelValF  = 20.
  res@cnLevelSpacingF = 2.      ; set the interval between contours
  map = gsn_csm_contour_map_ce(wks,u(0,26,:),res)
end if

status_exit(9)

end
