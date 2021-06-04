load "$NCARG_ROOT/lib/ncarg/nclscripts/esmf/ESMF_regridding.ncl"

begin

case         ="CESMLENS" ; CAMC20C or CESMLENS
current_year = 1996
comp_year    = 2019    ;2019 (+1), 2044 (+2), 2064 (+3), 2081/2082 (+4)   ; valid 1921-onwards
plevs        = False
output_diag  = False 
do_sw_fluxes = False
do_lw_fluxes = True

;------ internal settings
;keep_esmf=False
;esmf_remap = True
start_month  = 1
end_month    = 12

if (case .eq. "CAMC20C") ; Can only by used for attribution
  deltaFilePS = addfile("/glade/work/alyssas/deltas/delta_ps_CAM5-1-1degree_All-Hist_est1_v2.0_1996-2016.nc_Climatology_2016-1996.nc","r")
  deltaFileT  = addfile("/glade/work/alyssas/deltas/delta_ta_CAM5-1-1degree_All-Hist_est1_v2.0_1996-2016.nc_Climatology_2016-1996.nc","r")
  deltaFileQ  = addfile("/glade/work/alyssas/deltas/delta_hus_CAM5-1-1degree_All-Hist_est1_v2.0_1996-2016.nc_Climatology_2016-1996.nc","r")
  deltaPS_in  = deltaFilePS->delta_ps_Climatology_Monthly(start_month-1:end_month-1,0,:,:)
  deltaT_in   = deltaFileT->delta_ta_Climatology_Monthly(start_month-1:end_month-1,::-1,:,:)
  deltaQ_in   = deltaFileQ->delta_hus_Climatology_Monthly(start_month-1:end_month-1,::-1,:,:)
  deltaPS     = dim_avg_n_Wrap(deltaPS_in,0)
  deltaT      = dim_avg_n_Wrap(deltaT_in,0)
  deltaQ      = dim_avg_n_Wrap(deltaQ_in,0)
else if (case .eq. "CESMLENS") ; Can be used for attribution and future projections
  if (plevs) then
    deltaFilePS     = addfile("/glade/work/kareed/CESM_LENS_deltas/ens_PS_anom.nc","r")
    deltaFileT      = addfile("/glade/work/kareed/CESM_LENS_deltas/ens_T_anom.nc","r")
    deltaFileQ      = addfile("/glade/work/kareed/CESM_LENS_deltas/ens_Q_anom.nc","r")
  else
    deltaFilePS     = addfile("/glade/work/kareed/CESM_LENS_deltas/ens_PS_anom.nc","r")
    deltaFileT      = addfile("/glade/u/home/zarzycki/scratch/CESM_LENS_temp/T/ens_T_anom.nc","r")
    deltaFileQ      = addfile("/glade/u/home/zarzycki/scratch/CESM_LENS_temp/Q/ens_Q_anom.nc","r")
    if (do_sw_fluxes) then
      deltaFileFSDS      = addfile("/glade/u/home/zarzycki/scratch/CESM_LENS_temp/FSDS/ens_FSDS_anom.nc","r")
    end if
    if (do_lw_fluxes) then
      deltaFileFLDS      = addfile("/glade/u/home/zarzycki/scratch/CESM_LENS_temp/FLDS/ens_FLDS_anom.nc","r")
    end if
    deltaFileT_2    = addfile("/glade/scratch/zarzycki/CESM_LENS_temp/T/ens_012.nc","r")
    hyam_in=deltaFileT_2->hyam
    hybm_in=deltaFileT_2->hybm
  end if
  
  shouldermon=1
  deltaPS_in      = flt2dble(deltaFilePS->PS(current_year*12-1920*12+start_month-1-shouldermon:current_year*12-1920*12+end_month-1,:,:))
  deltaT_in       = flt2dble(deltaFileT->T(current_year*12-1920*12+start_month-1-shouldermon:current_year*12-1920*12+end_month-1,:,:,:))
  deltaQ_in       = flt2dble(deltaFileQ->Q(current_year*12-1920*12+start_month-1-shouldermon:current_year*12-1920*12+end_month-1,:,:,:))
  deltaPS_current = dim_avg_n_Wrap(deltaPS_in,0)
  deltaT_current  = dim_avg_n_Wrap(deltaT_in,0)
  deltaQ_current  = dim_avg_n_Wrap(deltaQ_in,0)
  if (do_sw_fluxes) then
    deltaFSDS_in       = flt2dble(deltaFileFSDS->FSDS(current_year*12-1920*12+start_month-1-shouldermon:current_year*12-1920*12+end_month-1,:,:))
    deltaFSDS_current  = dim_avg_n_Wrap(deltaFSDS_in,0)
  end if
  if (do_lw_fluxes) then
    deltaFLDS_in       = flt2dble(deltaFileFLDS->FLDS(current_year*12-1920*12+start_month-1-shouldermon:current_year*12-1920*12+end_month-1,:,:))
    deltaFLDS_current  = dim_avg_n_Wrap(deltaFLDS_in,0)
  end if
  
  if (comp_year .lt. 1920) then
    deltaPS_comp = deltaPS_current
    deltaT_comp  = deltaT_current
    deltaQ_comp  = deltaQ_current
    deltaPS_comp = 0.0
    deltaT_comp = 0.0
    deltaQ_comp = 0.0
    if (do_sw_fluxes) then
      deltaFSDS_comp  = deltaFSDS_current
      deltaFSDS_comp = 0.0
    end if
    if (do_lw_fluxes) then
      deltaFLDS_comp  = deltaFLDS_current
      deltaFLDS_comp = 0.0
    end if
  else
    deltaPS_in   = flt2dble(deltaFilePS->PS(comp_year*12-1920*12+start_month-1-shouldermon:comp_year*12-1920*12+end_month-1,:,:))
    deltaT_in    = flt2dble(deltaFileT->T(comp_year*12-1920*12+start_month-1-shouldermon:comp_year*12-1920*12+end_month-1,:,:,:))
    deltaQ_in    = flt2dble(deltaFileQ->Q(comp_year*12-1920*12+start_month-1-shouldermon:comp_year*12-1920*12+end_month-1,:,:,:))
    deltaPS_comp = dim_avg_n_Wrap(deltaPS_in,0)
    deltaT_comp  = dim_avg_n_Wrap(deltaT_in,0)
    deltaQ_comp  = dim_avg_n_Wrap(deltaQ_in,0)
    if (do_sw_fluxes) then
      deltaFSDS_in    = flt2dble(deltaFileFSDS->FSDS(comp_year*12-1920*12+start_month-1-shouldermon:comp_year*12-1920*12+end_month-1,:,:))
      deltaFSDS_comp  = dim_avg_n_Wrap(deltaFSDS_in,0)
    end if
    if (do_lw_fluxes) then
      deltaFLDS_in    = flt2dble(deltaFileFLDS->FLDS(comp_year*12-1920*12+start_month-1-shouldermon:comp_year*12-1920*12+end_month-1,:,:))
      deltaFLDS_comp  = dim_avg_n_Wrap(deltaFLDS_in,0)
    end if
  end if
  deltaPS = deltaPS_comp - deltaPS_current
  deltaT  = deltaT_comp  - deltaT_current
  deltaQ  = deltaQ_comp  - deltaQ_current
  copyatt(deltaPS,deltaPS_current)
  copyatt(deltaT,deltaT_current)
  copyatt(deltaQ,deltaQ_current)
  if (do_sw_fluxes) then
    deltaFSDS  = deltaFSDS_comp  - deltaFSDS_current
    copyatt(deltaFSDS,deltaFSDS_current)
  end if
  if (do_lw_fluxes) then
    deltaFLDS  = deltaFLDS_comp  - deltaFLDS_current
    copyatt(deltaFLDS,deltaFLDS_current)
  end if
  
end if
end if

print("filling in deltas...")
deltaT  = where(ismissing(deltaT),0,deltaT)
deltaQ  = where(ismissing(deltaQ),0,deltaQ)
deltaPS = where(ismissing(deltaPS),0,deltaPS)
if (do_sw_fluxes) then
  deltaFSDS  = where(ismissing(deltaFSDS),0,deltaFSDS)
end if
if (do_lw_fluxes) then
  deltaFLDS  = where(ismissing(deltaFLDS),0,deltaFLDS)
end if

print(deltaT_in&time@units)
time_deltas = cd_calendar(deltaT_in&time,0)
time_deltas_doy = day_of_year(toint(time_deltas(:,0)),toint(time_deltas(:,1)),toint(time_deltas(:,2)))
time_deltas_doy(12) = 366
print(time_deltas_doy)

; what is our datm file?
;datm_file_name="/glade/u/home/zarzycki/scratch/ERA5-DATM/DATM_perturb/TPHWL6Hrly/CMZERA5.v0.c2021.0.5d.TPQWL.1995-01.nc"
;datm2_file_name="/glade/u/home/zarzycki/scratch/ERA5-DATM/DATM_perturb/Precip6Hrly/CMZERA5.v0.c2021.0.5d.Prec.1995-01.nc"

print(datm_file_name+"")
print(datm2_file_name+"")
print(datm3_file_name+"")
datm_file=addfile(datm_file_name,"w")
datm2_file=addfile(datm2_file_name,"w")
datm3_file=addfile(datm3_file_name,"w")

;print(datm_file->time@units)
datm_deltas = cd_calendar(datm_file->time,0)
datm_deltas_doy = day_of_year(toint(datm_deltas(:,0)),toint(datm_deltas(:,1)),toint(datm_deltas(:,2)))

datm_deltas_doy := tofloat(datm_deltas_doy) + tofloat(datm_deltas(:,3))/24.
;print(datm_deltas_doy)

print("Time interpolation")
; interp LENS to time deltas in DATM
interp_T     = linint1_n_Wrap (time_deltas_doy, deltaT_in(:,29,:,:), False, datm_deltas_doy, 0, 0)
interp_Q     = linint1_n_Wrap (time_deltas_doy, deltaQ_in(:,29,:,:), False, datm_deltas_doy, 0, 0)
if (do_sw_fluxes) then
  interp_FSDS     = linint1_n_Wrap (time_deltas_doy, deltaFSDS_in(:,:,:), False, datm_deltas_doy, 0, 0)
end if
if (do_lw_fluxes) then
  interp_FLDS     = linint1_n_Wrap (time_deltas_doy, deltaFLDS_in(:,:,:), False, datm_deltas_doy, 0, 0)
end if
interp_PRECT = interp_T
interp_PRECT = 1.0 + (0.03 * interp_T)

print("Spatial interpolation output")
; now interp the time interpolated data to the DATM grid
datm_lat = datm_file->LATIXY(:,0)
datm_lon = datm_file->LONGXY(0,:)
interp_T_DATM     = linint2_Wrap (interp_T&lon,interp_T&lat,interp_T,True,datm_lon,datm_lat, 0)
interp_Q_DATM     = linint2_Wrap (interp_Q&lon,interp_Q&lat,interp_Q,True,datm_lon,datm_lat, 0)
if (do_sw_fluxes) then
  interp_FSDS_DATM     = linint2_Wrap (interp_FSDS&lon,interp_FSDS&lat,interp_FSDS,True,datm_lon,datm_lat, 0)
end if
if (do_lw_fluxes) then
  interp_FLDS_DATM     = linint2_Wrap (interp_FLDS&lon,interp_FLDS&lat,interp_FLDS,True,datm_lon,datm_lat, 0)
end if
interp_PRECT_DATM = linint2_Wrap (interp_PRECT&lon,interp_PRECT&lat,interp_PRECT,True,datm_lon,datm_lat, 0)

print("Read and update DATM data")
TBOT = datm_file->TBOT
TBOT_orig = TBOT
TBOT = TBOT + totype(interp_T_DATM, typeof(TBOT))

QBOT = datm_file->QBOT
QBOT_orig = QBOT
QBOT = QBOT + totype(interp_Q_DATM, typeof(QBOT))

if (do_sw_fluxes) then
  FSDS = datm3_file->FSDS
  FSDS_orig = FSDS
  FSDS = FSDS * totype(interp_FSDS_DATM, typeof(FSDS)) 
end if
if (do_lw_fluxes) then
  FLDS = datm_file->FLDS
  FLDS_orig = FLDS
  FLDS = FLDS * totype(interp_FLDS_DATM, typeof(FLDS)) 
end if

PRECTmms = datm2_file->PRECTmms
PRECTmms_orig = PRECTmms
PRECTmms = PRECTmms * totype(interp_PRECT_DATM, typeof(PRECTmms)) 

; write new fields to perturbation file
print("Writing output")
datm_file->TBOT      = TBOT
datm_file->QBOT      = QBOT
datm2_file->PRECTmms = PRECTmms
if (do_sw_fluxes) then
  datm3_file->FSDS = FSDS
end if
if (do_lw_fluxes) then
  datm_file->FLDS = FLDS
end if

if (output_diag) then
  diag_filename="datm.nc"
  print("outputting diags to "+diag_filename)
  
  system("/bin/rm -f "+diag_filename)   ; remove any pre-existing file
  ncdf = addfile(diag_filename ,"c")  ; open output netCDF file

  fAtt               = True            ; assign file attributes
  fAtt@creation_date = systemfunc ("date")
  fileattdef( ncdf, fAtt )            ; copy file attributes

  filedimdef(ncdf,"time",-1,True) 

  ncdf->TBOT_orig= TBOT_orig(0:4,:,:)
  ncdf->TBOT     = TBOT(0:4,:,:)
  ncdf->QBOT_orig= QBOT_orig(0:4,:,:)
  ncdf->QBOT     = QBOT(0:4,:,:)
  ncdf->PRECTmms_orig = PRECTmms_orig(0:4,:,:)
  ncdf->PRECTmms = PRECTmms(0:4,:,:)
  if (do_sw_fluxes) then
    ncdf->FSDS_orig= FSDS_orig(0:4,:,:)
    ncdf->FSDS     = FSDS(0:4,:,:)
  end if
  if (do_lw_fluxes) then
    ncdf->FLDS_orig= FLDS_orig(0:4,:,:)
    ncdf->FLDS     = FLDS(0:4,:,:)
  end if
  
  ncdf->interp_T_DATM = interp_T_DATM(0:4,:,:)
  ncdf->interp_Q_DATM = interp_Q_DATM(0:4,:,:)
  ncdf->interp_PRECT_DATM  = interp_PRECT_DATM(0:4,:,:)
  if (do_sw_fluxes) then
    ncdf->interp_FSDS_DATM = interp_FSDS_DATM(0:4,:,:)
  end if
  if (do_lw_fluxes) then
    ncdf->interp_FLDS_DATM = interp_FLDS_DATM(0:4,:,:)
  end if
  
end if

status_exit(9)

end


