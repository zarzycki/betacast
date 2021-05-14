load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/esmf/ESMF_regridding.ncl"

begin

case         ="CESMLENS" ; CAMC20C or CESMLENS
start_month  = 1 ;6
end_month    = 12 ;11
current_year = 1996 ; CESMLENS only
comp_year    = 2081 
correct_sfc=False
plevs=False
update_pressure=False
update_winds=False
do_ps_corr=True
esmf_remap = True
smooth_deltas = False
smooth_delta_iter = 10
output_diag=False 
keep_esmf=False
extra_diags=True

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
  if (comp_year .lt. 1920) then
    deltaPS_comp = deltaPS_current
    deltaT_comp  = deltaT_current
    deltaQ_comp  = deltaQ_current
    deltaPS_comp = 0.0
    deltaT_comp = 0.0
    deltaQ_comp = 0.0
  else
    deltaPS_in   = flt2dble(deltaFilePS->PS(comp_year*12-1920*12+start_month-1-shouldermon:comp_year*12-1920*12+end_month-1,:,:))
    deltaT_in    = flt2dble(deltaFileT->T(comp_year*12-1920*12+start_month-1-shouldermon:comp_year*12-1920*12+end_month-1,:,:,:))
    deltaQ_in    = flt2dble(deltaFileQ->Q(comp_year*12-1920*12+start_month-1-shouldermon:comp_year*12-1920*12+end_month-1,:,:,:))
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


print("filling in deltas...")
deltaT  = where(ismissing(deltaT),0,deltaT)
deltaQ  = where(ismissing(deltaQ),0,deltaQ)
deltaPS = where(ismissing(deltaPS),0,deltaPS)

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
datm_file=addfile(datm_file_name,"w")
datm2_file=addfile(datm2_file_name,"w")

;print(datm_file->time@units)
datm_deltas = cd_calendar(datm_file->time,0)
datm_deltas_doy = day_of_year(toint(datm_deltas(:,0)),toint(datm_deltas(:,1)),toint(datm_deltas(:,2)))

datm_deltas_doy := tofloat(datm_deltas_doy) + tofloat(datm_deltas(:,3))/24.
;print(datm_deltas_doy)

print("Time interpolation")
; interp LENS to time deltas in DATM
interp_T     = linint1_n_Wrap (time_deltas_doy, deltaT_in(:,29,:,:), False, datm_deltas_doy, 0, 0)
interp_Q     = linint1_n_Wrap (time_deltas_doy, deltaQ_in(:,29,:,:), False, datm_deltas_doy, 0, 0)
interp_PRECT = interp_T
interp_PRECT = 1.0 + (0.03 * interp_T)

print("Spatial interpolation output")
; now interp the time interpolated data to the DATM grid
datm_lat = datm_file->LATIXY(:,0)
datm_lon = datm_file->LONGXY(0,:)
interp_T_DATM     = linint2_Wrap (interp_T&lon,interp_T&lat,interp_T,True,datm_lon,datm_lat, 0)
interp_Q_DATM     = linint2_Wrap (interp_Q&lon,interp_Q&lat,interp_Q,True,datm_lon,datm_lat, 0)
interp_PRECT_DATM = linint2_Wrap (interp_PRECT&lon,interp_PRECT&lat,interp_PRECT,True,datm_lon,datm_lat, 0)

print("Read and update DATM data")
TBOT = datm_file->TBOT
TBOT_orig = TBOT
TBOT = TBOT + totype(interp_T_DATM, typeof(TBOT))

QBOT = datm_file->QBOT
QBOT_orig = QBOT
QBOT = QBOT + totype(interp_Q_DATM, typeof(QBOT)) 

PRECTmms = datm2_file->PRECTmms
PRECTmms_orig = PRECTmms
PRECTmms = PRECTmms * totype(interp_PRECT_DATM, typeof(PRECTmms)) 

; write new fields to perturbation file
print("Writing output")
datm_file->TBOT      = TBOT
datm_file->QBOT      = QBOT
datm2_file->PRECTmms = PRECTmms







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
  
end if

status_exit(9)

end


