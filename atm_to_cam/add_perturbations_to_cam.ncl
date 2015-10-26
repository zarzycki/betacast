load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/esmf/ESMF_regridding.ncl"

begin

;BEFOREPERTFILE=/glade/scratch/zarzycki/apply-haiyan-perturb/haiyan_48_x8_INIC_filter.nc
;AFTERPERTFILE=/glade/scratch/zarzycki/apply-haiyan-perturb/haiyan_48_x8_INIC_filter.nc_NOTREAL.nc

;system("cp /glade/p/work/zarzycki/sewx/INIC/ecsnow_30_x4_INIC.nc /glade/p/work/zarzycki/sewx/INIC/ecsnow_30_x4_TMP.nc")
;fileName="/glade/p/work/zarzycki/sewx/INIC/ecsnow_30_x4_TMP.nc"
system("cp "+BEFOREPERTFILE+" "+AFTERPERTFILE)
fileName=AFTERPERTFILE
camFile=addfile(fileName,"w")

deltaFileU=addfile("/glade/scratch/zarzycki/apply-haiyan-perturb/delta_ua_cam5.1_ACE_v2.0_2013.nc","r")
deltaFileV=addfile("/glade/scratch/zarzycki/apply-haiyan-perturb/delta_va_cam5.1_ACE_v2.0_2013.nc","r")
deltaFilePS=addfile("/glade/scratch/zarzycki/apply-haiyan-perturb/delta_ps_cam5.1_ACE_v2.0_2013.nc","r")
deltaFileT=addfile("/glade/scratch/zarzycki/apply-haiyan-perturb/delta_ta_cam5.1_ACE_v2.0_2013.nc","r")
deltaFileQ=addfile("/glade/scratch/zarzycki/apply-haiyan-perturb/delta_hus_cam5.1_ACE_v2.0_2013.nc","r")

deltaU = deltaFileU->delta_ua(10,::-1,:,:)
deltaV = deltaFileV->delta_va(10,::-1,:,:)
deltaPS = deltaFilePS->delta_ps(10,::-1,:,:)
deltaT = deltaFileT->delta_ta(10,::-1,:,:)
deltaQ = deltaFileQ->delta_hus(10,::-1,:,:)

hyam=camFile->hyam
hybm=camFile->hybm
P0=camFile->P0
PS=camFile->PS
U=camFile->U
V=camFile->V
T=camFile->T
Q=camFile->Q
lat=camFile->lat
lon=camFile->lon

printVarSummary(PS)
print(deltaU&plev)
print(P0)

print("Beginning PS to lat-lon interp...")
PS_deltaGrid=dsgrid2(lat,lon,PS(0,:),deltaU&lat,deltaU&lon)

printVarSummary(deltaU&plev)
printVarSummary(PS_deltaGrid)
printVarSummary(deltaU)

deltaU=where(ismissing(deltaU),0,deltaU)
deltaV=where(ismissing(deltaV),0,deltaV)
deltaT=where(ismissing(deltaT),0,deltaT)
deltaQ=where(ismissing(deltaQ),0,deltaQ)
deltaPS=where(ismissing(deltaPS),0,deltaPS)
;deltaU = linmsg_n(deltaU,(/-1,2/),1)
;deltaU = linmsg_n(deltaU,(/-1,2/),2)

print("Beginning vertical interp interp")
deltaUCAM=pres2hybrid_Wrap(deltaU&plev,PS_deltaGrid,P0,deltaU,hyam,hybm,1)
deltaVCAM=pres2hybrid_Wrap(deltaV&plev,PS_deltaGrid,P0,deltaV,hyam,hybm,1)
deltaTCAM=pres2hybrid_Wrap(deltaT&plev,PS_deltaGrid,P0,deltaT,hyam,hybm,1)
deltaQCAM=pres2hybrid_Wrap(deltaQ&plev,PS_deltaGrid,P0,deltaQ,hyam,hybm,1)
deltaPSCAM=deltaPS ; no need to vert interp PS

print("Beginning lat-lon to SE interp")
deltaUCAM_interp=linint2_points_Wrap(deltaUCAM&lon,deltaUCAM&lat,deltaUCAM,True,lon,lat,0)
deltaVCAM_interp=linint2_points_Wrap(deltaVCAM&lon,deltaVCAM&lat,deltaVCAM,True,lon,lat,0)
deltaTCAM_interp=linint2_points_Wrap(deltaTCAM&lon,deltaTCAM&lat,deltaTCAM,True,lon,lat,0)
deltaQCAM_interp=linint2_points_Wrap(deltaQCAM&lon,deltaQCAM&lat,deltaQCAM,True,lon,lat,0)
deltaPSCAM_interp=linint2_points_Wrap(deltaPSCAM&lon,deltaPSCAM&lat,deltaPSCAM,True,lon,lat,0)

deltaUCAM_interp=where(ismissing(deltaUCAM_interp),0,deltaUCAM_interp)
deltaVCAM_interp=where(ismissing(deltaVCAM_interp),0,deltaVCAM_interp)
deltaTCAM_interp=where(ismissing(deltaTCAM_interp),0,deltaTCAM_interp)
deltaQCAM_interp=where(ismissing(deltaQCAM_interp),0,deltaQCAM_interp)
deltaPSCAM_interp=where(ismissing(deltaPSCAM_interp),0,deltaPSCAM_interp)

printVarSummary(PS)
printVarSummary(deltaPSCAM_interp)
PS(0,:) = PS(0,:) + deltaPSCAM_interp(0,:)
T(0,:,:) = T(0,:,:) + deltaTCAM_interp(:,:)
V(0,:,:) = V(0,:,:) + deltaVCAM_interp(:,:)
U(0,:,:) = U(0,:,:) + deltaUCAM_interp(:,:)
Q(0,:,:) = Q(0,:,:) + deltaQCAM_interp(:,:)

; where moisture is negative due to deltas, reset
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

camFile->U = U
camFile->V = V
camFile->T = T
camFile->Q = Q
camFile->PS = PS

if (0 .eq. 1) then

  ; ===== OUTPUT
  system("/bin/rm -f simple.nc")   ; remove any pre-existing file
  ncdf = addfile("simple.nc" ,"c")  ; open output netCDF file

  ;===================================================================
  ; create global attributes of the file (optional)
  ;===================================================================
  fAtt               = True            ; assign file attributes
  fAtt@title         = "NCL Simple Approach to netCDF Creation"
  fAtt@source_file   =  "original-file.nc"
  fAtt@Conventions   = "None"
  fAtt@creation_date = systemfunc ("date")
  fileattdef( ncdf, fAtt )            ; copy file attributes

  ;===================================================================
  ; make time an UNLIMITED dimension; recommended  for most applications
  ;===================================================================
  filedimdef(ncdf,"time",-1,True) 

  ;===================================================================
  ; output variables directly; NCL will call appropriate functions
  ; to write the meta data associated with each variable
  ;===================================================================
  ;ncdf->deltaU=deltaU
  ncdf->PS_deltaGrid=PS_deltaGrid
  ncdf->deltaUCAM=deltaUCAM
  ncdf->deltaVCAM=deltaVCAM
  ncdf->deltaTCAM=deltaTCAM
  ncdf->deltaQCAM=deltaQCAM
  ncdf->deltaPSCAM=deltaPSCAM

end if

status_exit(9)

end


