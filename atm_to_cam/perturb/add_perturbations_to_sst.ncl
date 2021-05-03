load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/esmf/ESMF_regridding.ncl"

begin

case         ="CESMLENS" ; CAMC20C or CESMLENS
start_month  = 9 ;6 
end_month    = 9 ;11 
current_year = 1996 ; CESMLENS only
comp_year    = 2081 ; CESMLENS only enter 1850 or 1920-2100

system("cp "+BEFOREPERTFILE+" "+AFTERPERTFILE)
fileName=AFTERPERTFILE
camFile=addfile(fileName,"w")
lat=camFile->lat
lon=camFile->lon
ntimes=dimsizes(camFile->date)
SST=camFile->SST_cpl

if (case .eq. "CAMC20C") ; Can only by used for attribution
  deltaFileSST = addfile("/glade/work/alyssas/deltas/delta_ts_CAM5-1-1degree_All-Hist_est1_v2.0_1996-2016.nc_Climatology_2016-1996.nc","r")
  deltaSST_in  = deltaFileSST->delta_ts_Climatology_Monthly(start_month-1:end_month-1,0,:,:)
  deltaSST     = dim_avg_n_Wrap(deltaSST_in,0)
else if (case .eq. "CESMLENS") ; Can be used for attribution and future projections
  deltaFileSST     = addfile("/glade/work/kareed/CESM_LENS_deltas/ens_SST_anom.nc","r")
  deltaSST_in      = deltaFileSST->TS(current_year*12-1920*12+start_month-1:current_year*12-1920*12+end_month-1,:,:)
  deltaSST_current = dim_avg_n_Wrap(deltaSST_in,0) 
  if (comp_year .lt. 1920) then
    deltaSST_comp = deltaSST_current
    deltaSST_comp = 0.0     
  else
    deltaSST_in   = deltaFileSST->TS(comp_year*12-1920*12+start_month-1:comp_year*12-1920*12+end_month-1,:,:)
    deltaSST_comp = dim_avg_n_Wrap(deltaSST_in,0)
  end if  
  deltaSST = deltaSST_comp - deltaSST_current 
  copyatt(deltaSST,deltaSST_current)
end if
end if

deltaSST=where(ismissing(deltaSST),0,deltaSST)
deltaSST_interp=linint2_Wrap(deltaSST&lon,deltaSST&lat,deltaSST, True, lon,lat, 0)

print("max SST: "+max(deltaSST_interp)+"    min SST: "+min(deltaSST_interp))

do ii = 0,ntimes-1
  SST(ii,:,:) = SST(ii,:,:) + tofloat(deltaSST_interp(:,:))
end do

camFile->SST_cpl = SST

status_exit(9)

end


