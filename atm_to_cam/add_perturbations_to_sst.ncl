load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/esmf/ESMF_regridding.ncl"

begin

case="CAM_RCP85"   ; CAM, MIROC, CAM_RCP85

BEFOREPERTFILE="/glade/p/work/zarzycki/apply-haiyan-perturb/sst_1x1.nc"
AFTERPERTFILE ="/glade/p/work/zarzycki/apply-haiyan-perturb/sst_1x1.nc_DELTA_"+case+".nc"

system("cp "+BEFOREPERTFILE+" "+AFTERPERTFILE)
fileName=AFTERPERTFILE
camFile=addfile(fileName,"w")

if (case .eq. "CAM_RCP85")
  deltaFileSST=addfile("/glade/p/work/zarzycki/apply-haiyan-perturb/RCP85/delta_TS_CESM-CAM5-BGC-LE_209101-210012_minus_199601-200512.nc","r")
  deltaSST = deltaFileSST->delta_TS(10,:,:)
end if

lat=camFile->lat
lon=camFile->lon
ntimes=dimsizes(camFile->date)

SST=camFile->SST_cpl

deltaSST=where(ismissing(deltaSST),0,deltaSST)

deltaSST_interp=linint2_Wrap (deltaSST&lon,deltaSST&lat,deltaSST, True, lon,lat, 0)

print("max SST: "+max(deltaSST_interp)+"    min SST: "+min(deltaSST_interp))

do ii = 0,ntimes-1
  SST(ii,:,:) = SST(ii,:,:) + tofloat(deltaSST_interp(:,:))
end do

camFile->SST_cpl = SST

status_exit(9)

end


