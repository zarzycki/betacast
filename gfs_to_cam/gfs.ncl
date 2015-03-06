load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/esmf/ESMF_regridding.ncl"

begin

; Inputs needed
; machineid (integer)
; numlevels (integer)
; YYYYMMDDHH (string)
; gridname (string)
; ERA_dir (string)
; 

wcStrt = systemfunc("date")

; ===== Getting date from YYYYMMDDHH
dtime_map = (/4,2,2,2/)
splitDate = str_split_by_length(tostring(YYYYMMDDHH),dtime_map)
yearstr=splitDate(0)
monthstr=splitDate(1)
daystr=splitDate(2)
cyclestr=splitDate(3)
print("Regridding analysis from: "+yearstr+" "+monthstr+" "+daystr+" "+cyclestr+"Z")

;eraVertFile = addfile(ERA_dir+"/ERA-Interim_coordvars.nc","r")

print("---------------------------------------------------------")
print("Using this file: "+data_filename)
print("Using this remap: "+wgt_filename)

in = addfile(se_inic,"w")

; ====== load datafile
grb_file = addfile(data_filename,"r")

print("---------------------------------------------------------")
print("Loading lat/lon/lev")
if (datasource .eq. "GFS") then
  grblat = grb_file->lat_0
  grblon = grb_file->lon_0
  grblev = grb_file->lv_ISBL0
  rhlev = grb_file->lv_ISBL6
  cldlev = grb_file->lv_ISBL8
end if
if (datasource .eq. "CFSR") then
  grblat = grb_file->lat_0
  grblon = grb_file->lon_0
  grblev = grb_file->lv_ISBL0
  cldlev = grb_file->lv_ISBL4
end if
if (datasource .eq. "ERAI") then
  grblat = grb_file->latitude
  grblon = grb_file->longitude
  grblev = grb_file->level
end if

nlat = dimsizes(grblat)
nlon = dimsizes(grblon)


print("---------------------------------------------------------")
print("Loading variables")
if (datasource .eq. "GFS") then
  ps = grb_file->PRES_P0_L1_GLL0(:,:)
  p0   = 100000.
  pblh = grb_file->HPBL_P0_L1_GLL0(:,:)
  t_gfs = grb_file->TMP_P0_L100_GLL0(:,:,:)
  u_gfs = grb_file->UGRD_P0_L100_GLL0(:,:,:)
  v_gfs = grb_file->VGRD_P0_L100_GLL0(:,:,:)
  rh_gfs_native = grb_file->RH_P0_L100_GLL0(:,:,:)
  cldmix_gfs_native = grb_file->CLWMR_P0_L100_GLL0(:,:,:)
  print("Interpolating GRIB to uniform vertical levels")
  rh_gfs = int2p_n(rhlev,rh_gfs_native,grblev,2,0)
  cldmix_gfs = int2p_n(cldlev,cldmix_gfs_native,grblev,2,0)
  delete(rh_gfs_native)
  delete(cldmix_gfs_native)
  print("Calculating q from RH")
  q_gfs = mixhum_ptrh(conform(t_gfs,grblev,0)*0.01, t_gfs , rh_gfs, 2 ) ; specific humidity (g/kg)  
  print("Sorting bad values")
  cldmix_gfs = where(ismissing(cldmix_gfs),0,cldmix_gfs) 
  cldmix_gfs = where(cldmix_gfs.gt.0.01,0,cldmix_gfs)    
  print("seperating cloud ice and water")
   ; If T > 0 C, water, if less than < 0, ice
  cldice_gfs = cldmix_gfs
  cldliq_gfs = cldmix_gfs
  cldice_gfs = where(t_gfs.gt.273.15,0,cldice_gfs) 
  cldliq_gfs = where(t_gfs.gt.273.15,cldliq_gfs,0)    
  delete(cldmix_gfs)
  delete(rh_gfs)
end if
if (datasource .eq. "CFSR") then
  ps = grb_file->PRES_P0_L1_GLL0(:,:)
  p0   = 100000.
  pblh = grb_file->PRES_P0_L1_GLL0(:,:)
  t_gfs = grb_file->TMP_P0_L100_GLL0(:,:,:)
  u_gfs = grb_file->UGRD_P0_L100_GLL0(:,:,:)
  v_gfs = grb_file->VGRD_P0_L100_GLL0(:,:,:)
  rh_gfs = grb_file->RH_P0_L100_GLL0(:,:,:)
  cldmix_gfs_native = grb_file->CLWMR_P0_L100_GLL0(:,:,:)
  print("Interpolating GRIB to uniform vertical levels")
  cldmix_gfs = int2p_n(cldlev,cldmix_gfs_native,grblev,2,0)
  delete(cldmix_gfs_native)
  print("Calculating q from RH")
  q_gfs = mixhum_ptrh(conform(t_gfs,grblev,0)*0.01, t_gfs , rh_gfs, 2 ) ; specific humidity (g/kg)  
  print("Sorting bad values")
  cldmix_gfs = where(ismissing(cldmix_gfs),0,cldmix_gfs) 
  cldmix_gfs = where(cldmix_gfs.gt.0.01,0,cldmix_gfs)    
  print("seperating cloud ice and water")
   ; If T > 0 C, water, if less than < 0, ice
  cldice_gfs = cldmix_gfs
  cldliq_gfs = cldmix_gfs
  cldice_gfs = where(t_gfs.gt.273.15,0,cldice_gfs) 
  cldliq_gfs = where(t_gfs.gt.273.15,cldliq_gfs,0)    
  delete(cldmix_gfs)
  delete(rh_gfs)
end if
if (datasource .eq. "ERAI") then
  p0   = 100000.
  ; We need to use short2flt because ERA-Interim data is stored
  ; in short format and requires scale/offset
  ps = short2flt(grb_file->sp(0,:,:))
  pblh       = short2flt(grb_file->sp(0,:,:))
  t_gfs      = short2flt(grb_file->t(0,:,:,:))
  u_gfs      = short2flt(grb_file->u(0,:,:,:))
  v_gfs      = short2flt(grb_file->v(0,:,:,:))
  q_gfs      = short2flt(grb_file->q(0,:,:,:))
  cldice_gfs = short2flt(grb_file->ciwc(0,:,:,:))
  cldliq_gfs = short2flt(grb_file->clwc(0,:,:,:))
end if

print("=================================================================")
print("Max T: "+max(t_gfs)+"   min T: "+min(t_gfs))
print("Max U: "+max(u_gfs)+"   min U: "+min(u_gfs))
print("Max V: "+max(v_gfs)+"   min V: "+min(v_gfs))
print("Max Q: "+max(q_gfs)+"   min Q: "+min(q_gfs))
print("Max PS: "+max(ps)+"   min PS: "+min(ps))
print("Max CLDICE: "+max(cldice_gfs)+"   min CLDICE: "+min(cldice_gfs))
print("Max CLDLIQ: "+max(cldliq_gfs)+"   min CLDLIQ: "+min(cldliq_gfs))
print("Max PBLH: "+max(pblh)+"   min PBLH: "+min(pblh))
print("=================================================================")


print("---------------------------------------------------------")
print("Loading CAM levels")
if (numlevels .eq. 30)
print("Loading 30 level data")
fC   = addfile ("./L30template.nc", "r")
else
print("Loading 26 level data")
fC   = addfile ("./L26template.nc", "r")
end if
hya  = fC->hyam
hyb  = fC->hybm


print("Interpolating to CAM hybrid coordinates")
if (datasource .eq. "GFS" .or. datasource .eq. "CFSR") then
  print("Using pressure to hybrid interpolation")
  t_cam = pres2hybrid_Wrap(grblev,ps,p0,t_gfs,hya,hyb,4)
  u_cam = pres2hybrid_Wrap(grblev,ps,p0,u_gfs,hya,hyb,4)
  v_cam = pres2hybrid_Wrap(grblev,ps,p0,v_gfs,hya,hyb,4)
  q_cam = pres2hybrid_Wrap(grblev,ps,p0,q_gfs,hya,hyb,4)
  cldice_cam = pres2hybrid_Wrap(grblev,ps,p0,cldice_gfs,hya,hyb,4)
  cldliq_cam = pres2hybrid_Wrap(grblev,ps,p0,cldliq_gfs,hya,hyb,4)
end if
if (datasource .eq. "ERAI") then
  print("Using hybrid to hybrid interpolation")
  eraVertFile = addfile("/glade/p/work/zarzycki/getECMWFdata/ERA-Interim_coordvars.nc","r")
  hya_era = eraVertFile->a_model_ave
  hyb_era = eraVertFile->b_model_ave
  hya_era = hya_era/p0
  t_cam = hyi2hyo_Wrap(p0,hya_era,hyb_era,ps,t_gfs,hya,hyb,1)
  u_cam = hyi2hyo_Wrap(p0,hya_era,hyb_era,ps,u_gfs,hya,hyb,1)
  v_cam = hyi2hyo_Wrap(p0,hya_era,hyb_era,ps,v_gfs,hya,hyb,1)
  q_cam = hyi2hyo_Wrap(p0,hya_era,hyb_era,ps,q_gfs,hya,hyb,1)
  cldice_cam = hyi2hyo_Wrap(p0,hya_era,hyb_era,ps,cldice_gfs,hya,hyb,1)
  cldliq_cam = hyi2hyo_Wrap(p0,hya_era,hyb_era,ps,cldliq_gfs,hya,hyb,1)
end if

copy_VarCoords(t_cam,q_cam)
copy_VarCoords(t_cam,cldice_cam)
copy_VarCoords(t_cam,cldliq_cam)

t_cam!0    = "lev"
t_cam!1    = "lat"
t_cam!2    = "lon"
ps!0 = "lat"
ps!1 = "lon"
pblh!0 = "lat"
pblh!1 = "lon"

copy_VarCoords(t_cam,q_cam)
copy_VarCoords(t_cam,cldice_cam)
copy_VarCoords(t_cam,cldliq_cam)
copy_VarCoords(t_cam,u_cam)
copy_VarCoords(t_cam,v_cam)

print("Interpolating horizontal to CAM")
;fV   = addfile ("./L26template.nc", "r")
;fvlat  = fV->lat
;fvlon  = fV->lon
;fvslat = fV->slat
;fvslon = fV->slon

dstFileName = se_inic
dfile       = addfile(dstFileName,"r")
lonCell     = dfile->lon
latCell     = dfile->lat
Opt         = True
Opt@CopyVarCoords = False
ps_fv  = ESMF_regrid_with_weights(ps,wgt_filename,Opt)
pblh_fv  = ESMF_regrid_with_weights(pblh,wgt_filename,Opt)
t_fv  = ESMF_regrid_with_weights(t_cam,wgt_filename,Opt)
u_fv  = ESMF_regrid_with_weights(u_cam,wgt_filename,Opt)
v_fv  = ESMF_regrid_with_weights(v_cam,wgt_filename,Opt)
q_fv  = ESMF_regrid_with_weights(q_cam,wgt_filename,Opt)
cldice_fv  = ESMF_regrid_with_weights(cldice_cam,wgt_filename,Opt)
cldliq_fv  = ESMF_regrid_with_weights(cldliq_cam,wgt_filename,Opt)

sePS = in->PS
dim_sePS=dimsizes(sePS)
ncol=dim_sePS(1)

 print("Converting floats to doubles")
 ; Convert floats to doubles
 ps_fv_dbl = todouble(ps_fv)
 pblh_fv_dbl = todouble(pblh_fv)
 t_fv_dbl = todouble(t_fv)
 u_fv_dbl = todouble(u_fv)
 v_fv_dbl = todouble(v_fv)
 q_fv_dbl = todouble(q_fv)
 cldice_fv_dbl = todouble(cldice_fv)
 cldliq_fv_dbl = todouble(cldliq_fv)
 
 copy_VarMeta(ps_fv,ps_fv_dbl)
 copy_VarMeta(pblh_fv,pblh_fv_dbl)
 copy_VarMeta(t_fv,t_fv_dbl)
 copy_VarMeta(u_fv,u_fv_dbl)
 copy_VarMeta(v_fv,v_fv_dbl)
 copy_VarMeta(q_fv,q_fv_dbl)
 copy_VarMeta(cldliq_fv,cldliq_fv_dbl)
 copy_VarMeta(cldice_fv,cldice_fv_dbl)
; 
; ; clean up some stuff
 delete(ps_fv)
 delete(pblh_fv)
 delete(t_fv)
 delete(u_fv)
 delete(v_fv)
 delete(q_fv)
 delete(cldliq_fv)
 delete(cldice_fv)
; 
 print("Correcting time records")
; ; Add time record
 ps_fv_dbl_time = new((/1,ncol/),double)
 ps_fv_dbl_time(0,:) = ps_fv_dbl
 
 pblh_fv_dbl_time = new((/1,ncol/),double)
 pblh_fv_dbl_time(0,:) = pblh_fv_dbl

 u_fv_dbl_time = new((/1,numlevels,ncol/),double)
 u_fv_dbl_time(0,:,:) = u_fv_dbl
 
 v_fv_dbl_time = new((/1,numlevels,ncol/),double)
 v_fv_dbl_time(0,:,:) = v_fv_dbl

t_fv_dbl_time = new((/1,numlevels,ncol/),double)
t_fv_dbl_time(0,:,:) = t_fv_dbl

q_fv_dbl_time = new((/1,numlevels,ncol/),double)
q_fv_dbl_time(0,:,:) = q_fv_dbl

cldliq_fv_dbl_time = new((/1,numlevels,ncol/),double)
cldliq_fv_dbl_time(0,:,:) = cldliq_fv_dbl

cldice_fv_dbl_time = new((/1,numlevels,ncol/),double)
cldice_fv_dbl_time(0,:,:) = cldice_fv_dbl

ps_fv_dbl_time!0 = "time"
ps_fv_dbl_time!1 = "ncol"
u_fv_dbl_time!0 = "time"
u_fv_dbl_time!1 = "lev"
u_fv_dbl_time!2 = "ncol"
v_fv_dbl_time!0 = "time"
v_fv_dbl_time!1 = "lev"
v_fv_dbl_time!2 = "ncol"
t_fv_dbl_time!0 = "time"
t_fv_dbl_time!1 = "lev"
t_fv_dbl_time!2 = "ncol"
q_fv_dbl_time!0 = "time"
q_fv_dbl_time!1 = "lev"
q_fv_dbl_time!2 = "ncol"
cldliq_fv_dbl_time!0 = "time"
cldliq_fv_dbl_time!1 = "lev"
cldliq_fv_dbl_time!2 = "ncol"
cldice_fv_dbl_time!0 = "time"
cldice_fv_dbl_time!1 = "lev"
cldice_fv_dbl_time!2 = "ncol"

print("Writing output file")
in->PS=(/ps_fv_dbl_time/)
in->U=(/u_fv_dbl_time/)
in->V=(/v_fv_dbl_time/)
in->T=(/t_fv_dbl_time/)
in->Q=(/q_fv_dbl_time/)
in->CLDLIQ=(/cldliq_fv_dbl_time/)
in->CLDICE=(/cldice_fv_dbl_time/)
;  in->PBLH=(/pblh_fv_dbl_time/)


wallClockElapseTime(wcStrt, "Total runtime", 0)
print("done")
  
end


