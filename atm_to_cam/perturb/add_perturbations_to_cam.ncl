load "../../functions/namelist.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/esmf/ESMF_regridding.ncl"

begin

case  = keyword_values( pthi, "case", "string")
basedir  = keyword_values( pthi, "basedir", "string")
start_month = keyword_values( pthi, "start_month", "integer") 
end_month = keyword_values( pthi, "end_month", "integer") 
current_year = keyword_values( pthi, "current_year", "integer") 
comp_year = keyword_values( pthi, "comp_year", "integer") 

correct_sfc  = keyword_values( pthi, "correct_sfc", "logical")
plevs  = keyword_values( pthi, "plevs", "logical")
update_pressure  = keyword_values( pthi, "update_pressure", "logical")
update_winds  = keyword_values( pthi, "update_winds", "logical")
do_ps_corr  = keyword_values( pthi, "do_ps_corr", "logical")
esmf_remap  = keyword_values( pthi, "esmf_remap", "logical")
smooth_deltas  = keyword_values( pthi, "smooth_deltas", "logical")
if (smooth_deltas) then   ; only read iter if smooth_deltas is true!
  smooth_delta_iter  = keyword_values( pthi, "smooth_delta_iter", "integer")
end if
output_atm_diag  = keyword_values( pthi, "output_atm_diag", "logical")
keep_esmf  = keyword_values( pthi, "keep_esmf", "logical")
extra_diags_atm  = keyword_values( pthi, "extra_diags_atm", "logical")

; print configurations
print("************* Running ATM perturbation code *************")
print("Case: "+case)
print("start_month: "+start_month)
print("end_month: "+end_month)
print("current_year: "+current_year)
print("comp_year: "+comp_year)
print("BEFOREPERTFILE: "+BEFOREPERTFILE)
print("AFTERPERTFILE: "+AFTERPERTFILE)

print("plevs: "+plevs)
print("correct_sfc: "+correct_sfc)
print("update_pressure: "+update_pressure)
print("update_winds: "+update_winds)
print("do_ps_corr: "+do_ps_corr)

print("esmf_remap: "+esmf_remap)
print("keep_esmf: "+keep_esmf)

print("smooth_deltas: "+smooth_deltas)
if (smooth_deltas) then
  print("smooth_delta_iter: "+smooth_delta_iter)
end if

print("output_atm_diag: "+output_atm_diag)
print("extra_diags_atm: "+extra_diags_atm)

print("****************************************************")

system("cp" + " " + BEFOREPERTFILE + " " + AFTERPERTFILE)
fileName=AFTERPERTFILE
camFile=addfile(fileName,"w")
T=camFile->T
Q=camFile->Q
U=camFile->U
V=camFile->V
lat=camFile->lat
lon=camFile->lon
camlev=camFile->lev
hyam=camFile->hyam
hybm=camFile->hybm
hyai=camFile->hyai
hybi=camFile->hybi
P0 = 100000.
P0@units="Pa"
PS=camFile->PS
ncol = dimsizes(lat)
nlev = dimsizes(hyam)

if (case .eq. "CAMC20C") ; Can only by used for attribution
  deltaFilePS = addfile(basedir+"/"+case+"_plev/delta_ps_CAM5-1-1degree_All-Hist_est1_v2.0_1996-2016.nc_Climatology_2016-1996.nc","r")
  deltaFileT  = addfile(basedir+"/"+case+"_plev/delta_ta_CAM5-1-1degree_All-Hist_est1_v2.0_1996-2016.nc_Climatology_2016-1996.nc","r")
  deltaFileQ  = addfile(basedir+"/"+case+"_plev/delta_hus_CAM5-1-1degree_All-Hist_est1_v2.0_1996-2016.nc_Climatology_2016-1996.nc","r")
  deltaPS_in  = deltaFilePS->delta_ps_Climatology_Monthly(start_month-1:end_month-1,0,:,:)
  deltaT_in   = deltaFileT->delta_ta_Climatology_Monthly(start_month-1:end_month-1,::-1,:,:)
  deltaQ_in   = deltaFileQ->delta_hus_Climatology_Monthly(start_month-1:end_month-1,::-1,:,:)
  deltaPS     = dim_avg_n_Wrap(deltaPS_in,0)
  deltaT      = dim_avg_n_Wrap(deltaT_in,0)
  deltaQ      = dim_avg_n_Wrap(deltaQ_in,0)
else if (case .eq. "CESMLENS") ; Can be used for attribution and future projections
  if (plevs) then
    deltaFilePS     = addfile(basedir+"/"+case+"_plev/ens_PS_anom.nc","r")
    deltaFileT      = addfile(basedir+"/"+case+"_plev/ens_T_anom.nc","r")
    deltaFileQ      = addfile(basedir+"/"+case+"_plev/ens_Q_anom.nc","r") 
  else
    deltaFilePS     = addfile(basedir+"/"+case+"_mlev/PS/ens_PS_anom.nc","r")
    deltaFileT      = addfile(basedir+"/"+case+"_mlev/T/ens_T_anom.nc","r")
    deltaFileQ      = addfile(basedir+"/"+case+"_mlev/Q/ens_Q_anom.nc","r")
    hyam_in=deltaFileT->hyam
    hybm_in=deltaFileT->hybm
  end if
  deltaPS_in      = flt2dble(deltaFilePS->PS(current_year*12-1920*12+start_month-1:current_year*12-1920*12+end_month-1,:,:))
  deltaT_in       = flt2dble(deltaFileT->T(current_year*12-1920*12+start_month-1:current_year*12-1920*12+end_month-1,:,:,:))
  deltaQ_in       = flt2dble(deltaFileQ->Q(current_year*12-1920*12+start_month-1:current_year*12-1920*12+end_month-1,:,:,:))
  deltaPS_current = dim_avg_n_Wrap(deltaPS_in,0)
  deltaT_current  = dim_avg_n_Wrap(deltaT_in,0)
  deltaQ_current  = dim_avg_n_Wrap(deltaQ_in,0)
  if (comp_year .lt. 1920) then
    deltaPS_comp = deltaPS_current
    deltaT_comp  = deltaT_current
    deltaQ_comp  = deltaQ_current
    deltaPS_comp = 0.0
    deltaT_comp = 0.0
    deltaQ_comp = 0.0
  else
    deltaPS_in   = flt2dble(deltaFilePS->PS(comp_year*12-1920*12+start_month-1:comp_year*12-1920*12+end_month-1,:,:))
    deltaT_in    = flt2dble(deltaFileT->T(comp_year*12-1920*12+start_month-1:comp_year*12-1920*12+end_month-1,:,:,:))
    deltaQ_in    = flt2dble(deltaFileQ->Q(comp_year*12-1920*12+start_month-1:comp_year*12-1920*12+end_month-1,:,:,:))
    deltaPS_comp = dim_avg_n_Wrap(deltaPS_in,0)
    deltaT_comp  = dim_avg_n_Wrap(deltaT_in,0)
    deltaQ_comp  = dim_avg_n_Wrap(deltaQ_in,0)
  end if
  deltaPS = deltaPS_comp - deltaPS_current
  deltaT  = deltaT_comp  - deltaT_current
  deltaQ  = deltaQ_comp  - deltaQ_current
  copyatt(deltaPS,deltaPS_current)
  copyatt(deltaT,deltaT_current)
  copyatt(deltaQ,deltaQ_current)
end if
end if

print(deltaT_in&time)

print("filling in deltas...")
deltaT  = where(ismissing(deltaT),0,deltaT)
deltaQ  = where(ismissing(deltaQ),0,deltaQ)
deltaPS = where(ismissing(deltaPS),0,deltaPS)

if (esmf_remap) then
  print("Generating weights!")
  ; gen SE grid
  Opt_se               = True
  Opt_se@ForceOverwrite = True
  Opt_se@PrintTimings   = True
  Opt_se@Title          = "SE Grid"
  seGridName = "grid_se_"+ncol+".nc"
  if (.not. keep_esmf .or. .not.isfilepresent(seGridName)) then unstructured_to_ESMF(seGridName,lat,lon,Opt_se) end if

  ; gen RLL grid
  Opt_ll = True
  Opt_ll@ForceOverwrite = True
  Opt_ll@PrintTimings   = True
  Opt_ll@Title          = "Deltas grid"
  llGridName="grid_deltas.nc"
  if (.not. keep_esmf .or. .not.isfilepresent(llGridName)) then rectilinear_to_SCRIP(llGridName,deltaT&lat,deltaT&lon,Opt_ll) end if

  ;system ("echo $PATH")
  ;system ("echo $ESMFBINDIR")
  ;system ("which ESMF_RegridWeightGen")

  Opt                      = True
  Opt@InterpMethod         = "patch"
  Opt@ForceOverwrite       = True
  Opt@PrintTimings         = True
  Opt@SrcESMF        = True    
  Opt@DstESMF        = False
  wgtFileName1 = "map_se_"+ncol+"_to_ll.nc"
  if (.not. keep_esmf .or. .not.isfilepresent(wgtFileName1)) then ESMF_regrid_gen_weights(seGridName,llGridName,wgtFileName1,Opt) end if
  
  Opt@SrcESMF        = False    
  Opt@DstESMF        = True
  wgtFileName2 = "map_ll_to_se_"+ncol+".nc"
  if (.not. keep_esmf .or. .not.isfilepresent(wgtFileName2)) then ESMF_regrid_gen_weights(llGridName,seGridName,wgtFileName2,Opt) end if
  
  delete(Opt)
end if

print("Beginning PS to lat-lon interp...")
if (esmf_remap) then
  Opt                = True
  PS_deltaGrid = ESMF_regrid_with_weights(PS(0,:),wgtFileName1,Opt)
  delete(Opt)
else
  PS_deltaGrid=dsgrid2(lat,lon,PS(0,:),deltaT&lat,deltaT&lon)
end if
PS_deltaGrid@units="Pa"

print("Beginning vertical interp interp")
if (case .eq. "CAMC20C")
  deltaTCAM  = pres2hybrid_Wrap(deltaT&plev,PS_deltaGrid,P0,deltaT,hyam,hybm,1)
  deltaQCAM  = pres2hybrid_Wrap(deltaQ&plev,PS_deltaGrid,P0,deltaQ,hyam,hybm,1)
else if (case .eq. "CESMLENS") 
  if (plevs) then
    deltaTCAM  = pres2hybrid_Wrap(deltaT&plev,PS_deltaGrid,P0,deltaT,hyam,hybm,1)
    deltaQCAM  = pres2hybrid_Wrap(deltaQ&plev,PS_deltaGrid,P0,deltaQ,hyam,hybm,1)
  else
    deltaTCAM  = hyi2hyo_Wrap(P0,hyam_in,hybm_in,PS_deltaGrid,deltaT,hyam,hybm,1)
    deltaQCAM  = hyi2hyo_Wrap(P0,hyam_in,hybm_in,PS_deltaGrid,deltaQ,hyam,hybm,1)
  end if
  ; Starting to think that pres2hybrid_Wrap can't handle missing values
  ;deltaTCAM  = where(abs(deltaTCAM) .gt. 1e+10,0,deltaTCAM)
  ;deltaQCAM  = where(abs(deltaQCAM) .gt. 1e+10,0,deltaQCAM)
  ;deltaTCAM  = where(ismissing(deltaTCAM),0,deltaTCAM)
  ;deltaQCAM  = where(ismissing(deltaQCAM),0,deltaQCAM)
end if
end if
deltaPSCAM = deltaPS ; no need to vert interp PS
delete([/deltaT,deltaQ,deltaPS/])

if (smooth_deltas) then
  do ii = 0,smooth_delta_iter-1
    deltaQCAM  := smth9_Wrap(deltaQCAM, 0.50, 0.25, True)
    deltaTCAM  := smth9_Wrap(deltaTCAM, 0.50, 0.25, True)
    deltaPSCAM := smth9_Wrap(deltaPSCAM, 0.50, 0.25, True)
  end do
end if
  
if (extra_diags_atm) then
  dp = dpres_hybrid_ccm(PS_deltaGrid,P0,hyai,hybi)   ; dp(ntime,klevi-1,nlat,nlon)
  copy_VarCoords(deltaQCAM,dp)    ; dp(time,lev,lat,lon)
  pw = prcwater_dp(deltaQCAM(lat|:,lon|:,lev|:),dp(lat|:,lon|:,lev|:))    
  pw@long_name = "total column precipitable water"
  pw@units     = "kg/m2"
end if

print("Beginning lat-lon to SE interp")
if (esmf_remap) then
  Opt                = True
  deltaTCAM_interp  = ESMF_regrid_with_weights(deltaTCAM,wgtFileName2,Opt)
  deltaQCAM_interp  = ESMF_regrid_with_weights(deltaQCAM,wgtFileName2,Opt)
  deltaPSCAM_interp = ESMF_regrid_with_weights(deltaPSCAM,wgtFileName2,Opt)
  delete(Opt)
else
  deltaTCAM_interp  = linint2_points_Wrap(deltaTCAM&lon,deltaTCAM&lat,deltaTCAM,True,lon,lat,0)
  deltaQCAM_interp  = linint2_points_Wrap(deltaQCAM&lon,deltaQCAM&lat,deltaQCAM,True,lon,lat,0)
  deltaPSCAM_interp = linint2_points_Wrap(deltaPSCAM&lon,deltaPSCAM&lat,deltaPSCAM,True,lon,lat,0)
end if

deltaTCAM_interp  = where(ismissing(deltaTCAM_interp),0,deltaTCAM_interp)
deltaQCAM_interp  = where(ismissing(deltaQCAM_interp),0,deltaQCAM_interp)
deltaPSCAM_interp = where(ismissing(deltaPSCAM_interp),0,deltaPSCAM_interp)

if (do_ps_corr) then
  if (.not. esmf_remap) then
    print("do_ps_corr not supported without ESMF right now...")
    break
  end if
  ;---> perform an empirical correction to PS with emphasis over low PS areas
  anom_scaling=3.0   ; vertical average reference Tanom for scaling (Tanom/anom_scaling)
  print("Doing empirical ps_corr with anom_scaling: "+anom_scaling)
  
  print("Calculate pint and dpint")
  ; do weighted integral of deltaTCAM_interp
  dpint = deltaTCAM_interp
  dpint = 0.
  pint = new((/nlev+1,ncol/),typeof(hyai))
  pint = 0.
  
  pint = conform(pint,hyai,0)*P0 + conform(pint,hybi,0)*conform(pint,PS(0,:),1)
  ;pint = hyai*P0 + hybi*PS(0,ii)
  dpint = pint(1:nlev,:) - pint(0:nlev-1,:)
  ;do ii = 0,ncol-1
  ;  pint = hyai*P0 + hybi*PS(0,ii)
  ;  dpint(:,ii) = pint(1:nlev) - pint(0:nlev-1)
  ;end do
  
  print("Calculate weighted T anomaly and sign")
  Tanom = wgt_vertical_n(deltaTCAM_interp, dpint, 0, 0)
  ; Find where column is warmer and where is colder
  Tsign = where(Tanom .ge. 0.0, 1.0, -1.0)
  
  ; Map native -> RLL to smooth
  printVarSummary(Tsign)
  printVarSummary(deltaPSCAM)
  
  print("Remap to RLL")
  if (esmf_remap) then
    Tsign_deltaGrid=ESMF_regrid_with_weights(Tsign,wgtFileName1,False)
  else
    Tsign_deltaGrid=dsgrid2(lat,lon,Tsign,deltaPSCAM&lat,deltaPSCAM&lon)
  end if
  print("Smooth")
  ; smooth to not have step function flips from + to -
  smoothiter=50
  do ii = 0,smoothiter-1
    Tsign_deltaGrid := smth9(Tsign_deltaGrid, 0.50, 0.25, True)
  end do
  
  print(num(ismissing(Tsign_deltaGrid)))
  
  print("Remap to SE")
  ; back to native
  if (esmf_remap) then
    Tsign = ESMF_regrid_with_weights(Tsign_deltaGrid,wgtFileName2,False)
  else
    Tsign = linint2_points_Wrap(deltaPSCAM&lon,deltaPSCAM&lat,Tsign_deltaGrid,True,lon,lat,0)
  end if

  print("before "+num(ismissing(Tsign)))


  if (any(ismissing(Tsign))) then
    print("Found some missing Tsign data")
    print(num(ismissing(Tsign)))
    Tsignorig=Tsign
    do ii = 0,ncol-1
      if (ismissing(Tsign(ii))) then
        thisLat=lat(ii)
        thisLon=lon(ii)
        gcdist = gc_latlon(thisLat,thisLon,lat,lon,2,2)
        ; set the point where gc dist = 0 (i.e., this point) to a huge number
        ; so it isn't flagged as minix (which doesn't help us!)
        gcdist = where(ismissing(Tsignorig), 999999., gcdist)
        Tsign(ii) = Tsign(minind(gcdist))
        print("Replacing Tsign at "+ii+" with "+minind(gcdist))
        print("this lat/lon "+lat(ii)+" "+lon(ii)+"   nearest: "+lat(minind(gcdist))+" "+lon(minind(gcdist)))
      end if
    end do
    delete(Tsignorig)
  end if
  
  print("after "+num(ismissing(Tsign)))
  
  ; scale the multiplier based on the magnitude of the column anomaly.
  Tsign = Tsign * (abs(tofloat(Tanom))/anom_scaling)
  
  print(num(ismissing(Tsign)))

  ; correction coefficients derived from ne30 run
  rc=-0.007590353
  rc_intercept=771.9941
  print("CORR: using dPSL = "+rc+"*PS + "+rc_intercept)
  ; generate PS corr based on PS and constants
  PScorr = PS(0,:)*rc+rc_intercept
  
  print(num(ismissing(PScorr)))

  ; Scale correction by Tanom + heating/cooling
  PScorr = PScorr * Tsign
  
  print(num(ismissing(PScorr)))

  PS(0,:)  = PS(0,:)  + PScorr
end if

; Update arrays
T(0,:,:) = T(0,:,:) + deltaTCAM_interp(:,:)
Q(0,:,:) = Q(0,:,:) + deltaQCAM_interp(:,:)
if (update_pressure) then
  PS(0,:)  = PS(0,:)  + deltaPSCAM_interp(:)
end if
if (update_winds) then
  U(0,:,:) = U(0,:,:) + deltaUCAM_interp(:,:)
  V(0,:,:) = V(0,:,:) + deltaVCAM_interp(:,:)
end if

if (correct_sfc) then
  do ii = 0,nlev-1
    T(0,ii,:) = T(0,ii,:) - deltaTCAM_interp(nlev-1,:)
    Q(0,ii,:) = Q(0,ii,:) - deltaQCAM_interp(nlev-1,:)
  end do
end if

; where moisture is negative due to deltas, reset
print("Reset "+num(Q.le.0)+" Q elements for being negative")
Q=where(Q.le.0,1.0e-9,Q)


if(any(ismissing(Q))) then
  print("Q data contains some missing values. Beware.")
  print(num(ismissing(Q)))
end if

if(any(ismissing(T))) then
  print("T data contains some missing values. Beware.")
  print(num(ismissing(T)))
end if

if(any(ismissing(PS))) then
  print("PS data contains some missing values. Beware.")
  print(num(ismissing(PS)))
end if

print("Writing T and Q")
camFile->T = T
camFile->Q = Q
if (update_winds) then
  print("Writing U and V")
  camFile->U = U
  camFile->V = V
end if
if (update_pressure .or. do_ps_corr) then
  print("Writing PS")
  camFile->PS = PS
end if

if (output_atm_diag) then
  diag_filename="deltas_atm.nc"
  print("outputting diags to "+diag_filename)
  
  system("/bin/rm -f "+diag_filename)   ; remove any pre-existing file
  ncdf = addfile(diag_filename ,"c")  ; open output netCDF file

  fAtt               = True            ; assign file attributes
  fAtt@creation_date = systemfunc ("date")
  fileattdef( ncdf, fAtt )            ; copy file attributes

  filedimdef(ncdf,"time",-1,True) 

  ncdf->PS_deltaGrid=PS_deltaGrid
  ncdf->deltaPSCAM=deltaPSCAM
  ncdf->deltaTCAM=deltaTCAM
  ncdf->deltaQCAM=deltaQCAM

  if (extra_diags_atm) then
    ncdf->pw=pw
  end if
  
  if (do_ps_corr) then
    Opt=True
    if (esmf_remap) then
      ncdf->Tanom  = ESMF_regrid_with_weights(Tanom,wgtFileName1,Opt)
      ncdf->Tsign  = ESMF_regrid_with_weights(Tsign,wgtFileName1,Opt)
      ncdf->PScorr = ESMF_regrid_with_weights(PScorr,wgtFileName1,Opt)
      delete(Opt)
    else
      ncdf->Tanom   = dsgrid2(lat,lon,Tanom, deltaPSCAM&lat,deltaPSCAM&lon)
      ncdf->Tsign   = dsgrid2(lat,lon,Tsign, deltaPSCAM&lat,deltaPSCAM&lon)
      ncdf->PScorr  = dsgrid2(lat,lon,PScorr,deltaPSCAM&lat,deltaPSCAM&lon)
    end if
  end if

end if

status_exit(9)

end


