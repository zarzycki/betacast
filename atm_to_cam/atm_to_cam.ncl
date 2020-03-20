begin

wcStrt = systemfunc("date")

; Inputs needed
; machineid (integer)
; numlevels (integer)
; YYYYMMDDHH (string)
; gridname (string)
; ERA_dir (string)
;model_topo_file="/project/projectdirs/acme/inputdata/atm/cam/topo/USGS-gtopo30_ne30np4_16xdel2-PFC-consistentSGH.nc"

; Constants
gamma_s = 6.5/1000.   ; lapse rate in K/m
grav = 9.80616
Rd = 287.058

dycore="se"   ; fv or se

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

; ====== load datafile
grb_file = addfile(data_filename,"r")

print("---------------------------------------------------------")
print("Loading lat/lon/lev")
if (datasource .eq. "GFS" .or. datasource .eq. "CFSR") then
  dNames=getfilevardims(grb_file,"CLWMR_P0_L100_GLL0")
  cldlevName=dNames(0)
  delete(dNames)
  print("cldlev varname: "+cldlevName)
  dNames=getfilevardims(grb_file,"RH_P0_L100_GLL0")
  rhlevName=dNames(0)
  delete(dNames)
  print("rhlev varname: "+cldlevName)
  grblat = grb_file->lat_0
  grblon = grb_file->lon_0
  grblev = grb_file->lv_ISBL0
  rhlev  = grb_file->$rhlevName$
  cldlev = grb_file->$cldlevName$
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
if (datasource .eq. "GFS" .or. datasource .eq. "CFSR") then
  ps = grb_file->PRES_P0_L1_GLL0(:,:)
  p0   = 100000.
  ;pblh = grb_file->HPBL_P0_L1_GLL0(:,:)
  t_gfs = grb_file->TMP_P0_L100_GLL0(:,:,:)
  u_gfs = grb_file->UGRD_P0_L100_GLL0(:,:,:)

  v_gfs = grb_file->VGRD_P0_L100_GLL0(:,:,:)

  dNames:=getfilevardims(grb_file,"TMP_P0_L100_GLL0")
  templevName=dNames(0)
  dNames:=getfilevardims(grb_file,"UGRD_P0_L100_GLL0")
  windlevName=dNames(0)
  if ( "templevName" .ne. "windlevName" ) then
    u_gfs := int2p_n(grb_file->$windlevName$,u_gfs,grblev,2,0)
    v_gfs := int2p_n(grb_file->$windlevName$,v_gfs,grblev,2,0)
  end if
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
if (datasource .eq. "ERAI") then
  p0   = 100000.
  ; We need to use short2flt because ERA-Interim data is stored
  ; in short format and requires scale/offset
  ps = short2flt(grb_file->sp(0,:,:))
  ;pblh       = short2flt(grb_file->sp(0,:,:))
  t_gfs      = short2flt(grb_file->t(0,:,:,:))
  u_gfs      = short2flt(grb_file->u(0,:,:,:))
  v_gfs      = short2flt(grb_file->v(0,:,:,:))
  q_gfs      = short2flt(grb_file->q(0,:,:,:))
  cldice_gfs = short2flt(grb_file->ciwc(0,:,:,:))
  cldliq_gfs = short2flt(grb_file->clwc(0,:,:,:))
end if

print("=================================================================")
print("************ NATIVE DATA")
print("Max T: "+max(t_gfs)+"   min T: "+min(t_gfs))
print("Max U: "+max(u_gfs)+"   min U: "+min(u_gfs))
print("Max V: "+max(v_gfs)+"   min V: "+min(v_gfs))
print("Max Q: "+max(q_gfs)+"   min Q: "+min(q_gfs))
print("Max PS: "+max(ps)+"   min PS: "+min(ps))
print("Max CLDICE: "+max(cldice_gfs)+"   min CLDICE: "+min(cldice_gfs))
print("Max CLDLIQ: "+max(cldliq_gfs)+"   min CLDLIQ: "+min(cldliq_gfs))
;print("Max PBLH: "+max(pblh)+"   min PBLH: "+min(pblh))
print("=================================================================")


print("---------------------------------------------------------")
print("Loading CAM levels")
print("Loading "+numlevels+" level data")
fC   = addfile ("./L"+numlevels+"template.nc", "r")
hya  = fC->hyam
hyb  = fC->hybm
hyai  = fC->hyai
hybi  = fC->hybi
lev = fC->lev
ilev = fC->ilev


print("Interpolating to CAM hybrid coordinates")
if (datasource .eq. "GFS" .or. datasource .eq. "CFSR") then
  print("Using pressure to hybrid interpolation")
  printVarSummary(t_gfs)
  printVarSummary(u_gfs)
  t_cam = pres2hybrid_Wrap(grblev,ps,p0,t_gfs,hya,hyb,1)
  u_cam = pres2hybrid_Wrap(grblev,ps,p0,u_gfs,hya,hyb,1)
  v_cam = pres2hybrid_Wrap(grblev,ps,p0,v_gfs,hya,hyb,1)
  q_cam = pres2hybrid_Wrap(grblev,ps,p0,q_gfs,hya,hyb,1)
  cldice_cam = pres2hybrid_Wrap(grblev,ps,p0,cldice_gfs,hya,hyb,1)
  cldliq_cam = pres2hybrid_Wrap(grblev,ps,p0,cldliq_gfs,hya,hyb,1)
end if
if (datasource .eq. "ERAI") then
  print("Using hybrid to hybrid interpolation")
  eraVertFile = addfile("./ERA-Interim_coordvars.nc","r")
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

print("=================================================================")
print("************ AFTER VERTICAL INTERP")
print("Max T: "+max(t_cam)+"   min T: "+min(t_cam))
print("Max U: "+max(u_cam)+"   min U: "+min(u_cam))
print("Max V: "+max(v_cam)+"   min V: "+min(v_cam))
print("Max Q: "+max(q_cam)+"   min Q: "+min(q_cam))
print("Max PS: "+max(ps)+"   min PS: "+min(ps))
print("Max CLDICE: "+max(cldice_cam)+"   min CLDICE: "+min(cldice_cam))
print("Max CLDLIQ: "+max(cldliq_cam)+"   min CLDLIQ: "+min(cldliq_cam))
;print("Max PBLH: "+max(pblh)+"   min PBLH: "+min(pblh))
print("=================================================================")


copy_VarCoords(t_cam,q_cam)
copy_VarCoords(t_cam,cldice_cam)
copy_VarCoords(t_cam,cldliq_cam)

t_cam!0    = "lev"
t_cam!1    = "lat"
t_cam!2    = "lon"
ps!0 = "lat"
ps!1 = "lon"
;pblh!0 = "lat"
;pblh!1 = "lon"

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

;dstFileName = se_inic
;dfile       = addfile(dstFileName,"r")
;lonCell     = dfile->lon
;latCell     = dfile->lat
Opt         = True

; Here we need to keep VarCoords on to extract lat/lon
Opt@CopyVarCoords = True
ps_fv  = ESMF_regrid_with_weights(ps,wgt_filename,Opt)
selat=todouble(ps_fv@lat1d)
selon=todouble(ps_fv@lon1d)
delete_VarAtts(ps_fv,(/"lat1d","lon1d"/))
Opt@CopyVarCoords = False  ; now we can turn off
;pblh_fv  = ESMF_regrid_with_weights(pblh,wgt_filename,Opt)
t_fv  = ESMF_regrid_with_weights(t_cam,wgt_filename,Opt)
u_fv  = ESMF_regrid_with_weights(u_cam,wgt_filename,Opt)
v_fv  = ESMF_regrid_with_weights(v_cam,wgt_filename,Opt)
q_fv  = ESMF_regrid_with_weights(q_cam,wgt_filename,Opt)
cldice_fv  = ESMF_regrid_with_weights(cldice_cam,wgt_filename,Opt)
cldliq_fv  = ESMF_regrid_with_weights(cldliq_cam,wgt_filename,Opt)

print("=================================================================")
print("************ AFTER HORIZONTAL INTERP")
print("Max T: "+max(t_fv)+"   min T: "+min(t_fv))
print("Max U: "+max(u_fv)+"   min U: "+min(u_fv))
print("Max V: "+max(v_fv)+"   min V: "+min(v_fv))
print("Max Q: "+max(q_fv)+"   min Q: "+min(q_fv))
print("Max PS: "+max(ps)+"   min PS: "+min(ps))
print("Max CLDICE: "+max(cldice_fv)+"   min CLDICE: "+min(cldice_fv))
print("Max CLDLIQ: "+max(cldliq_fv)+"   min CLDLIQ: "+min(cldliq_fv))
;print("Max PBLH: "+max(pblh)+"   min PBLH: "+min(pblh))
print("=================================================================")

sePS = ps_fv
dim_sePS=dimsizes(sePS)
;printVarSummary(ps_fv)


if (isvar("model_topo_file") .and. fileexists(model_topo_file)) then
  print("Performing hydrostatic correction for surface pressure using "+model_topo_file)
;;;; Do hydrostatic correction for PS

  ; load additional fields from reanalysis
  topo_data = grb_file->HGT_P0_L1_GLL0(:,:)    ; gpm
  sfct_data = grb_file->TMP_P0_L104_GLL0(:,:)  ; 0.995 sigma T in K

  ; correct topo from gpm to m2/s2
  topo_data = topo_data * grav

  ; regrid to model grid
  topo_data_SE = ESMF_regrid_with_weights(topo_data,wgt_filename,Opt)
  sfct_data_SE = ESMF_regrid_with_weights(sfct_data,wgt_filename,Opt)

  ; load model orography
  ttfile = addfile(model_topo_file,"r")
  topo_model_SE = ttfile->PHIS     ;ncol

  ncol=dim_sePS(0)

  do kk = 0,ncol-1

    ; get difference in geopotential
    deltaPhi = topo_data_SE(kk) - topo_model_SE(kk)
    
    ; estimate "model's" surface temperature
    Tsfc_fv = sfct_data_SE(kk) + gamma_s * (deltaPhi / grav)

    ; calculate layer mean temperature for use in hydrostatic
    Tlayermean = (Tsfc_fv + sfct_data_SE(kk)) / 2.

    ; correct very warm and very cold layer means based on Trenberth 1993
    tcorriter=0
    if (Tlayermean .lt. 255.0) then
      Tlayermean = (255.0 + Tlayermean)/2.
      tcorriter=tcorriter+1
    else if (Tlayermean .gt. 290.5) then
      Tlayermean = (290.5 + Tlayermean)/2.
      tcorriter=tcorriter+1
    else
      ; no correction
    end if
    end if

    ; correct
    ps_orig = ps_fv(kk)
    beta = exp( tofloat(deltaPhi) / tofloat(Rd) / tofloat(Tlayermean) )
    ps_fv(kk) = ps_fv(kk) * beta

    ; correct T by shifting Tbot "down" the same delta
    nlev=dimsizes(lev)
    t_orig = t_fv(nlev-1,kk)
    delT = Tsfc_fv - sfct_data_SE(kk)
;    t_fv(nlev-1,kk) = t_fv(nlev-1,kk) + tofloat(delT)

    if (mod(kk,1000) .eq. 0) then
      print("Correcting: "+kk+" of "+(ncol-1)+"   from "+ps_orig+" to "+ps_fv(kk))
    end if

  end do
  print("needed to correct "+tcorriter+" temps for being too cold or too hot")
  delete([/beta,ps_orig,tcorriter,Tlayermean,Tsfc_fv,deltaPhi/])
  delete([/ttfile,topo_model_SE,topo_data_SE,sfct_data_SE,topo_data,sfct_data/])
else if (isvar("model_topo_file") .and. (model_topo_file .eq. " " .or. model_topo_file .eq. "NULL" .or. model_topo_file .eq. "")) then
  print("Empty model topo file entered, not performing hydro adjustment")
  print("continuing...")
else if (isvar("model_topo_file") .and. .not. fileexists(model_topo_file)) then
  print("model_topo_file passed in but cannot find file on Unix system")
  print("if you do not want adjustment, specify NULL in the namelist")
  print("exiting...")
  exit
else
  print("No model topo file passed into script, not performing hydro adjustment")
  print("continuing...")
end if
end if
end if

 print("Converting floats to doubles")
 ; Convert floats to doubles
 ps_fv_dbl = todouble(ps_fv)
 ;pblh_fv_dbl = todouble(pblh_fv)
 t_fv_dbl = todouble(t_fv)
 u_fv_dbl = todouble(u_fv)
 v_fv_dbl = todouble(v_fv)
 q_fv_dbl = todouble(q_fv)
 cldice_fv_dbl = todouble(cldice_fv)
 cldliq_fv_dbl = todouble(cldliq_fv)
 
 copy_VarMeta(ps_fv,ps_fv_dbl)
 ;copy_VarMeta(pblh_fv,pblh_fv_dbl)
 copy_VarMeta(t_fv,t_fv_dbl)
 copy_VarMeta(u_fv,u_fv_dbl)
 copy_VarMeta(v_fv,v_fv_dbl)
 copy_VarMeta(q_fv,q_fv_dbl)
 copy_VarMeta(cldliq_fv,cldliq_fv_dbl)
 copy_VarMeta(cldice_fv,cldice_fv_dbl)
; 
; ; clean up some stuff
 delete(ps_fv)
 ;delete(pblh_fv)
 delete(t_fv)
 delete(u_fv)
 delete(v_fv)
 delete(q_fv)
 delete(cldliq_fv)
 delete(cldice_fv)


if (dycore .eq. "se") then

ncol=dim_sePS(0)

 print("Correcting time records")
; ; Add time record
 ps_fv_dbl_time = new((/1,ncol/),double)
 ps_fv_dbl_time(0,:) = ps_fv_dbl
 
 ;pblh_fv_dbl_time = new((/1,ncol/),double)
 ;pblh_fv_dbl_time(0,:) = pblh_fv_dbl

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

selat!0 = "ncol"
selon!0 = "ncol"

else if (dycore .eq. "fv") then

nlat=dim_sePS(0)
nlon=dim_sePS(1)

; NEED TO INTERP SLAT/SLON
nslat=nlat-1
nslon=nlon
slat=new(nslat,"double")
slon=new(nslon,"double")
;lat=dim_sePS&lat
;lon=dim_sePS&lon
lat=todouble(fspan(-90.,90.,nlat))
dellon=360./nlon
lon=todouble(fspan(0,360.-dellon,nlon))

slat(:)=todouble((lat(0:nlat-2)+lat(1:nlat-1))/2.0)
slon = todouble(lon-(dellon/2.))

w_stag=latRegWgt(slat,"double",0)
print(w_stag)

 print("interpolating FV staggered")
u_fv_dbl_slat = linint2(lon,lat,u_fv_dbl,True,lon,slat,0) ;(ntim,mlat,mlon)
v_fv_dbl_slon = linint2(lon,lat,v_fv_dbl,True,slon,lat,0)

 print("Correcting time records")
; ; Add time record
 ps_fv_dbl_time = new((/1,nlat,nlon/),double)
 ps_fv_dbl_time(0,:,:) = ps_fv_dbl
 
 ;pblh_fv_dbl_time = new((/1,ncol/),double)
 ;pblh_fv_dbl_time(0,:) = pblh_fv_dbl

 u_fv_dbl_time = new((/1,numlevels,nslat,nlon/),double)
 u_fv_dbl_time(0,:,:,:) = u_fv_dbl_slat
 
 v_fv_dbl_time = new((/1,numlevels,nlat,nslon/),double)
 v_fv_dbl_time(0,:,:,:) = v_fv_dbl_slon

t_fv_dbl_time = new((/1,numlevels,nlat,nlon/),double)
t_fv_dbl_time(0,:,:,:) = t_fv_dbl

q_fv_dbl_time = new((/1,numlevels,nlat,nlon/),double)
q_fv_dbl_time(0,:,:,:) = q_fv_dbl

cldliq_fv_dbl_time = new((/1,numlevels,nlat,nlon/),double)
cldliq_fv_dbl_time(0,:,:,:) = cldliq_fv_dbl

cldice_fv_dbl_time = new((/1,numlevels,nlat,nlon/),double)
cldice_fv_dbl_time(0,:,:,:) = cldice_fv_dbl

ps_fv_dbl_time!0 = "time"
ps_fv_dbl_time!1 = "lat"
ps_fv_dbl_time!2 = "lon"

u_fv_dbl_time!0 = "time"
u_fv_dbl_time!1 = "lev"
u_fv_dbl_time!2 = "slat"
u_fv_dbl_time!3 = "lon"

v_fv_dbl_time!0 = "time"
v_fv_dbl_time!1 = "lev"
v_fv_dbl_time!2 = "lat"
v_fv_dbl_time!3 = "slon"

t_fv_dbl_time!0 = "time"
t_fv_dbl_time!1 = "lev"
t_fv_dbl_time!2 = "lat"
t_fv_dbl_time!3 = "lon"

q_fv_dbl_time!0 = "time"
q_fv_dbl_time!1 = "lev"
q_fv_dbl_time!2 = "lat"
q_fv_dbl_time!3 = "lon"

cldliq_fv_dbl_time!0 = "time"
cldliq_fv_dbl_time!1 = "lev"
cldliq_fv_dbl_time!2 = "lat"
cldliq_fv_dbl_time!3 = "lon"

cldice_fv_dbl_time!0 = "time"
cldice_fv_dbl_time!1 = "lev"
cldice_fv_dbl_time!2 = "lat"
cldice_fv_dbl_time!3 = "lon"

cldice_fv_dbl_time@units="kg/kg"
cldliq_fv_dbl_time@units="kg/kg"
q_fv_dbl_time@units="kg/kg"

lat!0="lat"
lat&lat=lat

lon!0="lon"
lon&lon=lon

slat!0="slat"
slat&slat=slat

slon!0="slon"
slon&slon=slon

else
  print("unsupported dycore")
	  exit
	end if
	end if


	; Fix error where Q goes very large (missing values?)
	qthresh=1.0
	q_fv_dbl_time=where(q_fv_dbl_time.gt.qthresh,0.0,q_fv_dbl_time)
	qthresh=1.0e-12
	q_fv_dbl_time=where(q_fv_dbl_time.le.qthresh,0.0,q_fv_dbl_time)
	if(any(ismissing(q_fv_dbl_time))) then
	  print("Q is missing data...")
	end if
	q_fv_dbl_time = where(ismissing(q_fv_dbl_time),0.0,q_fv_dbl_time)

	;q_fv_dbl_time=0.0

	print("Writing output file")
	;in->PS=(/ps_fv_dbl_time/)
	;in->U=(/u_fv_dbl_time/)
	;in->V=(/v_fv_dbl_time/)
	;in->T=(/t_fv_dbl_time/)
	;in->Q=(/q_fv_dbl_time/)
	;in->CLDLIQ=(/cldliq_fv_dbl_time/)
	;in->CLDICE=(/cldice_fv_dbl_time/)

	;  in->PBLH=(/pblh_fv_dbl_time/)

	other_stuff=True
	if (other_stuff) then
	system("/bin/rm -f "+se_inic)   ; remove any pre-existing file
	in2 = addfile(se_inic,"c")  ; open output netCDF file

	;===================================================================
	; create global attributes of the file (optional)
	;===================================================================
	fAtt               = True            ; assign file attributes
	fAtt@title         = "NCL Simple Approach to netCDF Creation"
	fAtt@source_file   =  "original-file.nc"
	fAtt@Conventions   = "None"
	fAtt@creation_date = systemfunc ("date")
	fileattdef( in2, fAtt )            ; copy file attributes

	;===================================================================
	; make time an UNLIMITED dimension; recommended  for most applications
	;===================================================================
	filedimdef(in2,"time",-1,True) 

	;===================================================================
	; output variables directly; NCL will call appropriate functions
	; to write the meta data associated with each variable
	;===================================================================
	in2->PS=ps_fv_dbl_time
	if(dycore .eq. "fv") then
	  in2->US=u_fv_dbl_time
	  in2->VS=v_fv_dbl_time
	else
	  in2->U=u_fv_dbl_time
	  in2->V=v_fv_dbl_time
	end if
	in2->T=t_fv_dbl_time
	in2->Q=q_fv_dbl_time
	in2->CLDLIQ=cldliq_fv_dbl_time
	in2->CLDICE=cldice_fv_dbl_time
	in2->hyam  = hya
	in2->hybm  = hyb
	in2->hyai  = hyai
	in2->hybi  = hybi
	in2->lev = lev
	in2->ilev = ilev
	if(dycore .eq. "fv") then
	  in2->lat=lat
	  in2->lon=lon
	  in2->slat=slat
	  in2->slon=slon
	  in2->w_stag=w_stag
	else  ; put lat/lon arrays on SE grid needed as of c2b08
	  in2->lat=selat
	  in2->lon=selon
	end if

	NUMLIQ=cldliq_fv_dbl_time
	NUMLIQ=0
	NUMICE=cldice_fv_dbl_time
	NUMICE=0

	;in2->NUMICE=NUMICE
	;in2->NUMLIQ=NUMLIQ

	end if

	wallClockElapseTime(wcStrt, "Total runtime", 0)
	print("done")

	status_exit(9)
  
end


