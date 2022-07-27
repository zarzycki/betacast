;Not needed NCL v6.5+
load "$NCARG_ROOT/lib/ncarg/nclscripts/esmf/ESMF_regridding.ncl"
load "helpers.ncl"

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
gamma_s = 5.0/1000.   ; env. lapse rate in K/m
gamma_d = 9.8/1000.   ; adiabatic lapse rate in K/m
grav = 9.80616
Rd = 287.058
Rv = 461.6
Rv_over_Rd = Rv / Rd

; Internal settings
debug=False
rho_d_algo = 2  ; 1 or 2
damp_upper_winds_mpas=False

; dycore options are FV, MPAS, and SE for now
; reverse compat -- if you don't specify dycore, default to SE
if (.not. isvar("dycore")) then
  dycore="se"
end if
print("USING DYCORE:  "+dycore)

; if we are using mpas, decide whether we are writing an MPAS file or CAM-ized file
if (dycore .eq. "mpas") then
  if (.not. isvar("mpas_as_cam")) then
    mpas_as_cam=False
  end if
  print("USING mpas_as_cam:  "+mpas_as_cam)
end if

; checks for whether an input is provided or not -- if not, set to default
if (.not. isvar("RDADIR")) then
  RDADIR=""
end if
if (.not. isvar("write_floats")) then
  write_floats=False
end if
if (.not. isvar("add_cloud_vars")) then
  add_cloud_vars=True
end if
print("USING RDADIR:  "+RDADIR)
print("USING write_floats:  "+write_floats)
print("USING add_cloud_vars:  "+add_cloud_vars)

; check if we are using MPAS if dycore has been passed in
if (isvar("model_topo_file") .and. isfilepresent(model_topo_file) .and. dycore .eq. "mpas" ) then
  print("Setting mpasfile to "+model_topo_file)
  mpasfile=model_topo_file
end if
if (.not. isvar("mpasfile") .and. dycore .eq. "mpas") then
  print("No mpasfile passed in but mpas dycore selected.")
  print("MPAS init file is needed to mimic + get zcoord from.  exiting...")
  exit
end if

; Toggle whether or not the output streams will be floats or doubles
if (write_floats) then
  write_type="float"
else
  write_type="double"
end if

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
  print("rhlev varname: "+rhlevName)
  grblat = grb_file->lat_0
  grblon = grb_file->lon_0
  grblev = grb_file->lv_ISBL0
  rhlev  = grb_file->$rhlevName$
  cldlev = grb_file->$cldlevName$
else if (datasource .eq. "RAP") then
  dNames=getfilevardims(grb_file,"RH_P0_L100_GLC0")
  rhlevName=dNames(0)
  delete(dNames)
  print("rhlev varname: "+rhlevName)
  grblat = grb_file->gridlat_0
  grblon = grb_file->gridlon_0
  grblev = grb_file->lv_ISBL0
  rhlev  = grb_file->$rhlevName$
else if (datasource .eq. "ERAI" .or. datasource .eq. "ERA5") then
  grblat = tofloat(grb_file->latitude)
  grblon = tofloat(grb_file->longitude)
  grblev = tofloat(grb_file->level)
  if (datasource .eq. "ERA5") then
    grblev = grblev * 100.
  end if
else if (datasource .eq. "ERA5RDA") then
  ;use RDA T file to extract dims
  pl_dir=RDADIR+"/e5.oper.an.pl/"+yearstr+monthstr
  sf_dir=RDADIR+"/e5.oper.an.sfc/"+yearstr+monthstr
  rda_find     := systemfunc("ls "+pl_dir+"/e5.oper.an.pl.128_130_t.ll025sc."+yearstr+monthstr+daystr+"00_*.nc")
  rda_file     := addfile(rda_find,"r")
  grblat        = tofloat(rda_file->latitude)
  grblon        = tofloat(rda_file->longitude)
  grblev        = tofloat(rda_file->level)
  grblev        = grblev * 100.
else if (datasource .eq. "CAM") then
  grblev := (/ 100, 200, 300, 500, 700, 1000, 2000, 3000, 5000, 7000, 10000, \
    15000, 20000, 25000, 30000, 35000, 40000, 45000, 50000, 55000, 60000,  \
    65000, 70000, 75000, 80000, 85000, 90000, 92500, 95000, 97500, 100000/)
  grblev := tofloat(grblev)
  mod_remap_file = "/glade/u/home/zarzycki/betacast/remapping/map_ne0np4natlanticref.ne30x4_TO_era5_0.25x0.25_patc.nc"
  tmp  = ESMF_regrid_with_weights(grb_file->PS(0,:),mod_remap_file,False)
  grblat = tmp&lat
  grblon = tmp&lon
  delete(tmp)
else
  print("Exiting when Loading lat/lon/lev")
  exit
end if
end if
end if
end if
end if

print("Level in information")
print("Number: "+dimsizes(grblev))
print("Max: "+max(grblev))
print("Min: "+min(grblev))

nlat = dimsizes(grblat)
nlon = dimsizes(grblon)

print("---------------------------------------------------------")
print("Loading variables")
if (datasource .eq. "GFS" .or. datasource .eq. "CFSR") then
  ps = grb_file->PRES_P0_L1_GLL0(:,:)
  p0   = 100000.
  t_gfs = grb_file->TMP_P0_L100_GLL0(:,:,:)
  u_gfs = grb_file->UGRD_P0_L100_GLL0(:,:,:)
  v_gfs = grb_file->VGRD_P0_L100_GLL0(:,:,:)  
  if (dycore .eq. "mpas") then
    w_gfs = grb_file->VVEL_P0_L100_GLL0(:,:,:)
    w_is_omega    = True
    z_gfs = grb_file->HGT_P0_L100_GLL0(:,:,:)
    z_is_phi      = False
  end if
  
  ; interpolate u, v, and w to t levels
  dNames:=getfilevardims(grb_file,"TMP_P0_L100_GLL0")
  templevName=dNames(0)
  dNames:=getfilevardims(grb_file,"UGRD_P0_L100_GLL0")
  windlevName=dNames(0)
  if ( templevName .ne. windlevName ) then
    print("interpolating u + v wind from "+dimsizes(grb_file->$windlevName$)+" to "+dimsizes(grblev)+" levels")
    u_gfs := int2p_n(grb_file->$windlevName$,u_gfs,grblev,2,0)
    v_gfs := int2p_n(grb_file->$windlevName$,v_gfs,grblev,2,0)
  end if
  if (dycore .eq. "mpas") then
    dNames:=getfilevardims(grb_file,"VVEL_P0_L100_GLL0")
    windlevName=dNames(0)
    if ( templevName .ne. windlevName ) then
      print("interpolating w wind from "+dimsizes(grb_file->$windlevName$)+" to "+dimsizes(grblev)+" levels")
      w_gfs := int2p_n(grb_file->$windlevName$,w_gfs,grblev,2,0)
    end if
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



if (datasource .eq. "RAP") then
  ps = grb_file->PRES_P0_L1_GLC0(:,:)
  p0   = 100000.
  t_gfs = grb_file->TMP_P0_L100_GLC0(:,:,:)
  u_gfs = grb_file->UGRD_P0_L100_GLC0(:,:,:)
  v_gfs = grb_file->VGRD_P0_L100_GLC0(:,:,:)  
  if (dycore .eq. "mpas") then
    w_gfs = grb_file->VVEL_P0_L100_GLC0(:,:,:)
    w_is_omega    = True
    z_gfs = grb_file->HGT_P0_L100_GLC0(:,:,:)
    z_is_phi      = False
  end if
  
  ; interpolate u, v, and w to t levels
  dNames:=getfilevardims(grb_file,"TMP_P0_L100_GLC0")
  templevName=dNames(0)
  dNames:=getfilevardims(grb_file,"UGRD_P0_L100_GLC0")
  windlevName=dNames(0)
  if ( templevName .ne. windlevName ) then
    print("interpolating u + v wind from "+dimsizes(grb_file->$windlevName$)+" to "+dimsizes(grblev)+" levels")
    u_gfs := int2p_n(grb_file->$windlevName$,u_gfs,grblev,2,0)
    v_gfs := int2p_n(grb_file->$windlevName$,v_gfs,grblev,2,0)
  end if
  dNames:=getfilevardims(grb_file,"VVEL_P0_L100_GLC0")
  windlevName=dNames(0)
  if ( templevName .ne. windlevName ) then
    print("interpolating w wind from "+dimsizes(grb_file->$windlevName$)+" to "+dimsizes(grblev)+" levels")
    w_gfs := int2p_n(grb_file->$windlevName$,w_gfs,grblev,2,0)
  end if
    
  rh_gfs_native = grb_file->RH_P0_L100_GLC0(:,:,:)
  cldmix_gfs_native = grb_file->RH_P0_L100_GLC0(:,:,:)
  cldmix_gfs_native = 0.0
  print("Interpolating GRIB to uniform vertical levels")
  rh_gfs = int2p_n(rhlev,rh_gfs_native,grblev,2,0)
  cldmix_gfs = cldmix_gfs_native
  delete(rh_gfs_native)
  delete(cldmix_gfs_native)
  print("Calculating q from RH")
  q_gfs = mixhum_ptrh(conform(t_gfs,grblev,0)*0.01, t_gfs , rh_gfs, 2 ) ; specific humidity (g/kg)  
  print("Sorting bad values")
  cldmix_gfs = where(ismissing(cldmix_gfs),0,cldmix_gfs) 
  cldmix_gfs = where(cldmix_gfs.gt.0.01,0,cldmix_gfs)    
  print("seperating cloud ice and water")
  cldice_gfs = cldmix_gfs
  cldliq_gfs = cldmix_gfs
  cldice_gfs = 0.0
  cldliq_gfs = 0.0
  delete(cldmix_gfs)
  delete(rh_gfs)
end if



if (datasource .eq. "ERAI" .or. datasource .eq. "ERA5") then
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

if (datasource .eq. "ERA5RDA") then
  ; define some stuff that is static
  p0   = 100000.
  opt=0
  opt@return_type = "int"

  ; go through var by var, get relevant file, extract time
  rda_find     := systemfunc("ls "+sf_dir+"/e5.oper.an.sfc.128_134_sp.ll025sc."+yearstr+monthstr+"0100_*.nc")
  print(rda_find+"")
  rda_file     := addfile(rda_find,"r")
  rda_time     := rda_file->time
  rda_thistime := cd_inv_calendar(toint(yearstr),toint(monthstr),toint(daystr),toint(cyclestr),0,0,rda_time@units,opt)
  ps            = rda_file->SP(closest_val(rda_thistime,rda_time),:,:)

  rda_find     := systemfunc("ls "+pl_dir+"/e5.oper.an.pl.128_130_t.ll025sc."+yearstr+monthstr+daystr+"00_*.nc")
  print(rda_find+"")
  rda_file     := addfile(rda_find,"r")
  rda_time     := rda_file->time
  rda_thistime := cd_inv_calendar(toint(yearstr),toint(monthstr),toint(daystr),toint(cyclestr),0,0,rda_time@units,opt)
  t_gfs         = rda_file->T(closest_val(rda_thistime,rda_time),:,:,:)

  rda_find     := systemfunc("ls "+pl_dir+"/e5.oper.an.pl.128_131_u.ll025uv."+yearstr+monthstr+daystr+"00_*.nc")
  print(rda_find+"")
  rda_file     := addfile(rda_find,"r")
  rda_time     := rda_file->time
  rda_thistime := cd_inv_calendar(toint(yearstr),toint(monthstr),toint(daystr),toint(cyclestr),0,0,rda_time@units,opt)
  u_gfs         = rda_file->U(closest_val(rda_thistime,rda_time),:,:,:)

  rda_find     := systemfunc("ls "+pl_dir+"/e5.oper.an.pl.128_132_v.ll025uv."+yearstr+monthstr+daystr+"00_*.nc")
  print(rda_find+"")
  rda_file     := addfile(rda_find,"r")
  rda_time     := rda_file->time
  rda_thistime := cd_inv_calendar(toint(yearstr),toint(monthstr),toint(daystr),toint(cyclestr),0,0,rda_time@units,opt)
  v_gfs         = rda_file->V(closest_val(rda_thistime,rda_time),:,:,:)

  rda_find     := systemfunc("ls "+pl_dir+"/e5.oper.an.pl.128_133_q.ll025sc."+yearstr+monthstr+daystr+"00_*.nc")
  print(rda_find+"")
  rda_file     := addfile(rda_find,"r")
  rda_time     := rda_file->time
  rda_thistime := cd_inv_calendar(toint(yearstr),toint(monthstr),toint(daystr),toint(cyclestr),0,0,rda_time@units,opt)
  q_gfs         = rda_file->Q(closest_val(rda_thistime,rda_time),:,:,:)

  rda_find     := systemfunc("ls "+pl_dir+"/e5.oper.an.pl.128_246_clwc.ll025sc."+yearstr+monthstr+daystr+"00_*.nc")
  print(rda_find+"")
  rda_file     := addfile(rda_find,"r")
  rda_time     := rda_file->time
  rda_thistime := cd_inv_calendar(toint(yearstr),toint(monthstr),toint(daystr),toint(cyclestr),0,0,rda_time@units,opt)
  cldliq_gfs    = rda_file->CLWC(closest_val(rda_thistime,rda_time),:,:,:)

  rda_find     := systemfunc("ls "+pl_dir+"/e5.oper.an.pl.128_247_ciwc.ll025sc."+yearstr+monthstr+daystr+"00_*.nc")
  print(rda_find+"")
  rda_file     := addfile(rda_find,"r")
  rda_time     := rda_file->time
  rda_thistime := cd_inv_calendar(toint(yearstr),toint(monthstr),toint(daystr),toint(cyclestr),0,0,rda_time@units,opt)
  cldice_gfs    = rda_file->CIWC(closest_val(rda_thistime,rda_time),:,:,:)
  
  if (dycore .eq. "mpas") then
    rda_find     := systemfunc("ls "+pl_dir+"/e5.oper.an.pl.128_135_w.ll025sc."+yearstr+monthstr+daystr+"00_*.nc")
    print(rda_find+"")
    rda_file     := addfile(rda_find,"r")
    rda_time     := rda_file->time
    rda_thistime := cd_inv_calendar(toint(yearstr),toint(monthstr),toint(daystr),toint(cyclestr),0,0,rda_time@units,opt)
    w_gfs         = rda_file->W(closest_val(rda_thistime,rda_time),:,:,:)
    w_is_omega    = True
  
    rda_find     := systemfunc("ls "+pl_dir+"/e5.oper.an.pl.128_129_z.ll025sc."+yearstr+monthstr+daystr+"00_*.nc")
    print(rda_find+"")
    rda_file     := addfile(rda_find,"r")
    rda_time     := rda_file->time
    rda_thistime := cd_inv_calendar(toint(yearstr),toint(monthstr),toint(daystr),toint(cyclestr),0,0,rda_time@units,opt)
    z_gfs         = rda_file->Z(closest_val(rda_thistime,rda_time),:,:,:)
    z_is_phi      = True
  end if
end if


if (datasource .eq. "CAM") then

  mod_remap_file = "/glade/u/home/zarzycki/betacast/remapping/map_ne0np4natlanticref.ne30x4_TO_era5_0.25x0.25_patc.nc"
  mod_in_topo = "/glade/u/home/zarzycki/work/unigridFiles/ne0np4natlanticref.ne30x4/topo/topo_ne0np4natlanticref.ne30x4_smooth.nc"

  ;grb_file is loaded above
  opt=0
  opt@return_type = "int"
  cam_time := grb_file->time
  cam_thistime := cd_inv_calendar(toint(yearstr),toint(monthstr),toint(daystr),toint(cyclestr),0,0,cam_time@units,opt)
  delete(opt)
  
  ps_mod = grb_file->PS(closest_val(cam_thistime,cam_time),:)
  p0   = 100000.
  t_mod = grb_file->T(closest_val(cam_thistime,cam_time),:,:)
  u_mod = grb_file->U(closest_val(cam_thistime,cam_time),:,:)
  v_mod = grb_file->V(closest_val(cam_thistime,cam_time),:,:)
  q_mod = grb_file->Q(closest_val(cam_thistime,cam_time),:,:)
  cldice_mod = t_mod
  cldliq_mod = t_mod
  cldice_mod = 0.0
  cldliq_mod = 0.0
  
  mod_in_topo_f = addfile(mod_in_topo,"r")
  phis_mod=mod_in_topo_f->PHIS
    
  ; remap from CAM to lat-lon
  Opt         = True
  ps  = ESMF_regrid_with_weights(ps_mod,mod_remap_file,Opt)
  t_rll  = ESMF_regrid_with_weights(t_mod,mod_remap_file,Opt)
  u_rll  = ESMF_regrid_with_weights(u_mod,mod_remap_file,Opt)
  v_rll  = ESMF_regrid_with_weights(v_mod,mod_remap_file,Opt)
  q_rll  = ESMF_regrid_with_weights(q_mod,mod_remap_file,Opt)
  cldice_rll  = ESMF_regrid_with_weights(cldice_mod,mod_remap_file,Opt)
  cldliq_rll  = ESMF_regrid_with_weights(cldliq_mod,mod_remap_file,Opt)

  phis_rll  = ESMF_regrid_with_weights(phis_mod,mod_remap_file,Opt)
  delete([/ps_mod,t_mod,u_mod,v_mod,q_mod,phis_mod,cldliq_mod,cldice_mod/])
  delete(Opt)
  
  ; get relevant coeffs from CAM file
  hyam=grb_file->hyam
  hybm=grb_file->hybm
  hyai=grb_file->hyai
  hybi=grb_file->hybi
  mod_lev = dimsizes(hyam)

  if (dycore .eq. "mpas") then
    ; calculate omega
    omega_rll = omega_ccm_driver(p0,ps,u_rll,v_rll,hyam,hybm,hyai,hybi)
  
    ; remove omega poleward of OMEGA_LAT_THRESH to deal with singularity
    OMEGA_LAT_THRESH=86.0
    omega_rll = where(abs(conform(omega_rll,omega_rll&lat,1)).gt.OMEGA_LAT_THRESH,0.0,omega_rll)
  
    ; convert OMEGA to W
    p3d_rll = pres_hybrid_ccm(ps,p0,hyam,hybm)
    w_rll = omega_to_w(omega_rll, p3d_rll, t_rll)         ; wc[*][*][*][*] (m/s)
  
    ; calculate Z
    tkv_rll = t_rll
    tkv_rll = t_rll*(1.+0.61*q_rll)  
    ;z_rll = hydro (p3d_rll,tkv_rll,phis_rll/9.80655)
    z_rll = cz2ccm(ps,phis_rll,tkv_rll,p0,hyam(::-1),hybm(::-1),hyai(::-1),hybi(::-1))
    copy_VarCoords(t_rll,z_rll)  
  end if

  ; now go to constant pressure surfaces
;  grblev := (/ 100, 200, 300, 500, 700, 1000, 2000, 3000, 5000, 7000, 10000, \
;    15000, 20000, 25000, 30000, 35000, 40000, 45000, 50000, 55000, 60000,  \
;    65000, 70000, 75000, 80000, 85000, 90000, 92500, 95000, 97500, 100000/)
;  grblev := tofloat(grblev)
  
  P0mb=p0/100.                          ; p0 in mb
  intyp = 1                             ; 1=linear, 2=log, 3=log-log
  kxtrp = True                          ; True=extrapolate
  tbot_mod = t_rll(mod_lev-1,:,:)             ; bottom level temperature

  ; note, this function needs levels in mb
  ; the bit after kxtrp needs to be 1 for T, -1 for Z, 0 for all else
  t_gfs := vinth2p_ecmwf(t_rll,hyam,hybm,grblev/100,ps,intyp,P0mb,\
                 1,kxtrp,1,tbot_mod,phis_rll)
  u_gfs := vinth2p_ecmwf(u_rll,hyam,hybm,grblev/100,ps,intyp,P0mb,\
                 1,kxtrp,0,tbot_mod,phis_rll)               
  v_gfs := vinth2p_ecmwf(v_rll,hyam,hybm,grblev/100,ps,intyp,P0mb,\
                 1,kxtrp,0,tbot_mod,phis_rll)      
  q_gfs := vinth2p_ecmwf(q_rll,hyam,hybm,grblev/100,ps,intyp,P0mb,\
                 1,kxtrp,0,tbot_mod,phis_rll)
  cldice_gfs := vinth2p_ecmwf(cldice_rll,hyam,hybm,grblev/100,ps,intyp,P0mb,\
                 1,kxtrp,0,tbot_mod,phis_rll)       
  cldliq_gfs := vinth2p_ecmwf(cldliq_rll,hyam,hybm,grblev/100,ps,intyp,P0mb,\
                 1,kxtrp,0,tbot_mod,phis_rll)
  if (dycore .eq. "mpas") then
    w_gfs := vinth2p_ecmwf(w_rll,hyam,hybm,grblev/100,ps,intyp,P0mb,\
                   1,kxtrp,0,tbot_mod,phis_rll)       
    z_gfs := vinth2p_ecmwf(z_rll,hyam,hybm,grblev/100,ps,intyp,P0mb,\
                   1,kxtrp,-1,tbot_mod,phis_rll)
  end if
  
  delete([/kxtrp,intyp,P0mb/])
  delete([/t_rll,u_rll,v_rll,q_rll,cldice_rll,cldliq_rll/])
  delete([/hyam,hybm,hyai,hybi,mod_lev/])       
  delete([/cam_time,cam_thistime/])
  if (dycore .eq. "mpas") then
    delete([/w_rll,z_rll,omega_rll/])
    delete([/p3d_rll,tkv_rll/])
  end if

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

if (dycore .eq. "mpas") then
  print("Max Z: "+max(z_gfs)+"   min Z: "+min(z_gfs))
  print("------------Deriving variables----------------")

  ; initialize variables
  theta_gfs = t_gfs
  rho_gfs = t_gfs

  ; create a 3-D pressure field from constant pressure surfaces
  pres_gfs = conform(t_gfs,grblev,0)

  ; if w is omega (Pa/s), convert to vertical velocity in m/s
  if (w_is_omega) then
    w_gfs := omega_to_w(w_gfs, pres_gfs, t_gfs)
  end if

  ; if z is reported in geopotential, convert to geometric height
  if (z_is_phi) then
    ;phis_sfc = grb_file->Z(0,:,:)
    ;z_gfs = z_gfs - conform(z_gfs,phis_sfc,(/1,2/))
    z_gfs = z_gfs / grav  
  end if

  ; potential temperature is calculated using full pressure and actual T
  theta_gfs = pot_temp(pres_gfs, t_gfs, -1, False)

  ; rho we need to calculate with *dry* air
  if (rho_d_algo .eq. 1) then
    presdry_gfs = pres_gfs
    presdry_gfs = pres_gfs/(1.+q_gfs)
    ;presdry_gfs = presdry_gfs*0.995
    rho_gfs = presdry_gfs / Rd / t_gfs
  else
    print("using rho_d calculation from internal CAM/MPAS code")
    rho_gfs = pres_gfs / Rd / t_gfs / (1. + Rv_over_Rd * q_gfs)
  end if

  print("Max RHO_DRY: "+max(rho_gfs)+"   min RHO_DRY: "+min(rho_gfs))
  print("Max THETA: "+max(theta_gfs)+"   min THETA: "+min(theta_gfs))
  print("Max PRES: "+max(pres_gfs)+"   min PRES: "+min(pres_gfs))
  print("Max W: "+max(w_gfs)+"   min W: "+min(w_gfs))
end if

print("=================================================================")

if (dycore .ne. "mpas") then
  print("---------------------------------------------------------")
  print("Loading CAM levels")
  print("Loading "+numlevels+" level data")
  fC   = addfile ("./templates/L"+numlevels+"template.nc", "r")
  hya  = fC->hyam
  hyb  = fC->hybm
  hyai  = fC->hyai
  hybi  = fC->hybi
  lev = fC->lev
  ilev = fC->ilev

  print("Interpolating to CAM hybrid coordinates")
  if (datasource .eq. "GFS" .or. datasource .eq. "CFSR" .or. datasource .eq. "ERA5" .or. datasource .eq. "ERA5RDA" .or. datasource .eq. "RAP" .or. datasource .eq. "CAM") then
    print("Using pressure to hybrid interpolation")
    t_cam = pres2hybrid_Wrap(grblev,ps,p0,t_gfs,hya,hyb,1)
    u_cam = pres2hybrid_Wrap(grblev,ps,p0,u_gfs,hya,hyb,1)
    v_cam = pres2hybrid_Wrap(grblev,ps,p0,v_gfs,hya,hyb,1)
    q_cam = pres2hybrid_Wrap(grblev,ps,p0,q_gfs,hya,hyb,1)
    cldice_cam = pres2hybrid_Wrap(grblev,ps,p0,cldice_gfs,hya,hyb,1)
    cldliq_cam = pres2hybrid_Wrap(grblev,ps,p0,cldliq_gfs,hya,hyb,1)
  end if
  if (datasource .eq. "ERAI") then
    print("Using hybrid to hybrid interpolation")
    eraVertFile = addfile("./templates/ERA-Interim_coordvars.nc","r")
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
  print("=================================================================")

  print("==CLEAN after vert interp")
  delete([/u_gfs,v_gfs,t_gfs,q_gfs,cldice_gfs,cldliq_gfs/])

  ; set var coords
  t_cam!0    = "lev"
  t_cam!1    = "lat"
  t_cam!2    = "lon"
  ps!0 = "lat"
  ps!1 = "lon"

  ; copy var coords over to q, cldliq, and clidice
  copy_VarCoords(t_cam,q_cam)
  copy_VarCoords(t_cam,cldice_cam)
  copy_VarCoords(t_cam,cldliq_cam)
  copy_VarCoords(t_cam,u_cam)
  copy_VarCoords(t_cam,v_cam)

  print("Interpolating horizontal to CAM")
  Opt         = True

  ; Here we need to keep VarCoords on to extract lat/lon
  Opt@CopyVarCoords = True
  ps_fv  = ESMF_regrid_with_weights(ps,wgt_filename,Opt)
  if (dycore .ne. "fv") then
    selat=todouble(ps_fv@lat1d)
    selon=todouble(ps_fv@lon1d)
    delete_VarAtts(ps_fv,(/"lat1d","lon1d"/))
  end if
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
  print("=================================================================")

  print("clearing after horiz interp")
  delete([/u_cam,v_cam,t_cam,q_cam,cldice_cam,cldliq_cam,ps/])

else   ; if mpas
  print("We are doing MPAS!")
  
  ; First, we have to load the MPAS init file and get the zgrid.
  mpasf=addfile(mpasfile,"r")
  mpas_z = mpasf->zgrid
  mpas_dims = dimsizes(mpas_z)
  mpas_ncell = mpas_dims(0)
  mpas_nlevi = mpas_dims(1)
  mpas_nlev = mpas_nlevi - 1
  print("MPAS GRID: nlev: "+mpas_nlev+"    nlevi: "+mpas_nlevi+"     ncell: "+mpas_ncell)
  ; Let's just set up some dummy arrays since we aren't going to vert interp first
  t_cam = t_gfs
  u_cam = u_gfs
  v_cam = v_gfs
  q_cam = q_gfs
  z_cam = z_gfs
  cldice_cam = cldice_gfs
  cldliq_cam = cldliq_gfs
  theta_cam = theta_gfs
  rho_cam = rho_gfs
  w_cam = w_gfs

  print("Interpolating horizontal to "+dycore)
  Opt         = True
  ; Here we need to keep VarCoords on to extract lat/lon
  Opt@CopyVarCoords = True
  ps_fv  = ESMF_regrid_with_weights(ps,wgt_filename,Opt)
  selat=todouble(ps_fv@lat1d)
  selon=todouble(ps_fv@lon1d)
  delete_VarAtts(ps_fv,(/"lat1d","lon1d"/))
  Opt@CopyVarCoords = False  ; now we can turn off
  t_fv  = ESMF_regrid_with_weights(t_cam,wgt_filename,Opt)
  u_fv  = ESMF_regrid_with_weights(u_cam,wgt_filename,Opt)
  v_fv  = ESMF_regrid_with_weights(v_cam,wgt_filename,Opt)
  q_fv  = ESMF_regrid_with_weights(q_cam,wgt_filename,Opt)
  cldice_fv  = ESMF_regrid_with_weights(cldice_cam,wgt_filename,Opt)
  cldliq_fv  = ESMF_regrid_with_weights(cldliq_cam,wgt_filename,Opt)
  z_fv  = ESMF_regrid_with_weights(z_cam,wgt_filename,Opt)
  theta_fv  = ESMF_regrid_with_weights(theta_cam,wgt_filename,Opt)
  rho_fv  = ESMF_regrid_with_weights(rho_cam,wgt_filename,Opt)
  w_fv  = ESMF_regrid_with_weights(w_cam,wgt_filename,Opt)
  
  ; now vertical interpolation
  t_wrf = new((/mpas_nlev,mpas_ncell/),"double")
  rho_wrf = new((/mpas_nlev,mpas_ncell/),"double")
  theta_wrf = new((/mpas_nlev,mpas_ncell/),"double")
  w_wrf = new((/mpas_nlevi,mpas_ncell/),"double")
  q_wrf = new((/mpas_nlev,mpas_ncell/),"double")
  u_wrf = new((/mpas_nlev,mpas_ncell/),"double")
  v_wrf = new((/mpas_nlev,mpas_ncell/),"double")

  print("Performing vertical interpolation at each MPAS column")
  do ix = 0,mpas_ncell-1
    if (mod(ix,10000) .eq. 0) then
      print("... MPAS_VERT_INTERP: "+(100.*ix/(mpas_ncell-1))+"%")
    end if
    zmid = (mpas_z(ix,1:mpas_nlevi-1) + mpas_z(ix,0:mpas_nlevi-2)) / 2.
    zint = mpas_z(ix,:)
    if (mpas_as_cam) then
      t_wrf(:,ix) = z_to_z_interp(t_fv(::-1,ix),z_fv(::-1,ix),zmid,True,True)
    end if
    theta_wrf(:,ix) = z_to_z_interp(theta_fv(::-1,ix),z_fv(::-1,ix),zmid,True,True)
    rho_wrf(:,ix)   = z_to_z_interp(rho_fv(::-1,ix),z_fv(::-1,ix),zmid,True,True)
    w_wrf(:,ix)     = z_to_z_interp(w_fv(::-1,ix),z_fv(::-1,ix),zint,True,True)
    q_wrf(:,ix)     = z_to_z_interp(q_fv(::-1,ix),z_fv(::-1,ix),zmid,True,True)
    ; u and v we don't extrapolate aloft to prevent wind speeds from getting too high
    u_wrf(:,ix)     = z_to_z_interp(u_fv(::-1,ix),z_fv(::-1,ix),zmid,True,False)
    v_wrf(:,ix)     = z_to_z_interp(v_fv(::-1,ix),z_fv(::-1,ix),zmid,True,False)
  end do

  ; delete the analysis data that is on the MPAS grid since we now have interpolated MPAS + vert
  delete([/theta_fv,rho_fv,w_fv,q_fv,u_fv,v_fv/])
 
  print("Setting lower BC for W so flow can't go through surface")
  w_wrf(:,0) = 0.0

  if (damp_upper_winds_mpas) then
    print("Damping upper level MPAS winds")
    mpas_damp_coefs=(/0.90,0.95,0.98/)
    u_wrf(:,mpas_nlev-1) = u_wrf(:,mpas_nlev-1)*mpas_damp_coefs(0)
    u_wrf(:,mpas_nlev-2) = u_wrf(:,mpas_nlev-2)*mpas_damp_coefs(1)
    u_wrf(:,mpas_nlev-3) = u_wrf(:,mpas_nlev-3)*mpas_damp_coefs(2)
    v_wrf(:,mpas_nlev-1) = v_wrf(:,mpas_nlev-1)*mpas_damp_coefs(0)
    v_wrf(:,mpas_nlev-2) = v_wrf(:,mpas_nlev-2)*mpas_damp_coefs(1)
    v_wrf(:,mpas_nlev-3) = v_wrf(:,mpas_nlev-3)*mpas_damp_coefs(2)
  end if

  if (.not. mpas_as_cam) then
    ; put u + v on cell edges...
    print("Projecting u + v to velocity normal to edge")
    uNorm_wrf = uv_cell_to_edge(u_wrf,v_wrf,mpas_nlev,mpasf->lonEdge,mpasf->latEdge,mpasf->lonCell,mpasf->latCell,mpasf->edgeNormalVectors,mpasf->cellsOnEdge)
    ; delete u and v components since we have mapped to edge normals
    delete([/u_wrf,v_wrf/])
    uNorm_wrf!0 = "lev"
    uNorm_wrf!1 = "nEdges"
  else
    ; we just need to keep U and V at cell centers
    u_wrf!0 = "lev"
    u_wrf!1 = "ncol"
    v_wrf!0 = "lev"
    v_wrf!1 = "ncol"
    ; we also have T which we don't have in the MPAS code
    t_wrf!0 = "lev"
    t_wrf!1 = "ncol"
    ; finally, we keep ps_fv which never needed vert interp
    ps_fv!0 = "ncol"
  end if

  theta_wrf!0 = "lev"
  theta_wrf!1 = "ncol"
  rho_wrf!0 = "lev"
  rho_wrf!1 = "ncol"
  w_wrf!0 = "ilev"
  w_wrf!1 = "ncol"
  q_wrf!0 = "lev"
  q_wrf!1 = "ncol"

  if (debug) then
    system("/bin/rm -f mpas_debug.nc")   ; remove any pre-existing file
    ncdf = addfile("mpas_debug.nc" ,"c")  ; open output netCDF file
    fAtt               = True            ; assign file attributes
    fAtt@creation_date = systemfunc ("date")
    fileattdef( ncdf, fAtt )            ; copy file attributes
    filedimdef(ncdf,"time",-1,True) 
    ncdf->theta_wrf = theta_wrf
    ncdf->rho_wrf = rho_wrf
    ncdf->w_wrf = w_wrf
    ncdf->q_wrf = q_wrf
    ncdf->uNorm_wrf = uNorm_wrf
  end if
  
  if (mpas_as_cam) then
    
    print("Writing MPAS as CAM file, presumably for nudging")
    
    ; take _wrf vars, add time coord for writing -- also flip vert index!
    u_fv_dbl_time  = add_time_define_precision(u_wrf(::-1,:),write_type,True)
    v_fv_dbl_time  = add_time_define_precision(v_wrf(::-1,:),write_type,True)
    t_fv_dbl_time  = add_time_define_precision(t_wrf(::-1,:),write_type,True)
    q_fv_dbl_time  = add_time_define_precision(q_wrf(::-1,:),write_type,True)
    ps_fv_dbl_time = add_time_define_precision(ps_fv,write_type,True)

    ; clean up
    delete([/u_wrf,v_wrf,t_wrf,q_wrf,ps_fv/])

    ; cheap way of doing this, should update to efficient way later
    system("/bin/rm -f "+se_inic)   ; remove any pre-existing file
    ncdf = addfile(se_inic,"c")  ; open output netCDF file
    fAtt               = True            ; assign file attributes
    fAtt@creation_date = systemfunc ("date")
    fileattdef( ncdf, fAtt )            ; copy file attributes
    filedimdef(ncdf,"time",-1,True) 
    ncdf->U  = u_fv_dbl_time
    ncdf->V  = v_fv_dbl_time
    ncdf->T  = t_fv_dbl_time
    ncdf->Q  = q_fv_dbl_time
    ncdf->PS = ps_fv_dbl_time

  else
  
    print("Writing MPAS fields to a copy of existing init file")
    system("/bin/rm -f -v "+se_inic)   ; remove any pre-existing file
    system("cp -v "+mpasfile+" "+se_inic)   ; remove any pre-existing file
    nc = addfile(se_inic ,"w")  ; open output netCDF file
    nc->u(0,:,:)     = (/ uNorm_wrf(nEdges|:,lev|:) /)
    nc->qv(0,:,:)    = (/ q_wrf(ncol|:,lev|:) /)
    nc->rho(0,:,:)   = (/ rho_wrf(ncol|:,lev|:) /)
    nc->theta(0,:,:) = (/ theta_wrf(ncol|:,lev|:) /)
    nc->w(0,:,:)     = (/ w_wrf(ncol|:,ilev|:) /)
    
  end if
  
  wallClockElapseTime(wcStrt, "Total runtime", 0)
  print("done with MPAS, now exiting gracefully...")
  status_exit(9)
  
end if

sePS = ps_fv
dim_sePS=dimsizes(sePS)

if (dycore .ne. "mpas" .and. isvar("model_topo_file") .and. fileexists(model_topo_file)) then
  print("Performing hydrostatic correction for surface pressure using "+model_topo_file)

  ; set flags for config
  tempadjustflag=""
  if (isvar("adjust_config") .and. adjust_config .ne. "") then
    tempadjustflag=str_get_cols(adjust_config,0,0)
  end if
  print("tempadjustflag set to: "+tempadjustflag)

  ; load additional fields from reanalysis
  if (datasource .eq. "GFS" .or. datasource .eq. "CFSR") then
    topo_data = grb_file->HGT_P0_L1_GLL0(:,:)    ; gpm
    sfct_data = grb_file->TMP_P0_L104_GLL0(:,:)  ; 0.995 sigma T in K
    ; correct topo from gpm to m2/s2
    topo_data = topo_data * grav
  else if (datasource .eq. "RAP") then
    topo_data = grb_file->HGT_P0_L1_GLC0(:,:)    ; gpm
    sfct_data = grb_file->TMP_P0_L1_GLC0(:,:)  ; TS
    ; correct topo from gpm to m2/s2
    topo_data = topo_data * grav
  else if (datasource .eq. "ERAI" .or. datasource .eq. "ERA5") then
    topo_data = short2flt(grb_file->phis(0,:,:))    ; gpm
    sfct_data = short2flt(grb_file->t2m(0,:,:))  ; 0.995 sigma T in K
  else if (datasource .eq. "ERA5RDA") then
    topo_data     = grb_file->Z(0,:,:)
    rda_find     := systemfunc("ls "+sf_dir+"/e5.oper.an.sfc.128_167_2t.ll025sc."+yearstr+monthstr+"0100_*.nc")
    rda_file     := addfile(rda_find,"r")
    rda_time     := rda_file->time
    rda_thistime := cd_inv_calendar(toint(yearstr),toint(monthstr),toint(daystr),toint(cyclestr),0,0,rda_time@units,opt)
    sfct_data     = rda_file->VAR_2T(closest_val(rda_thistime,rda_time),:,:)
  else if (datasource .eq. "CAM") then
    topo_data = phis_rll
    sfct_data = tbot_mod
  else
    print("no avail reanalysis surface data for corr... exiting...")
    exit
  end if
  end if
  end if
  end if
  end if
  
  ; regrid to model grid
  topo_data_SE = ESMF_regrid_with_weights(topo_data,wgt_filename,Opt)
  sfct_data_SE = ESMF_regrid_with_weights(sfct_data,wgt_filename,Opt)

  ; load model orography
  ttfile = addfile(model_topo_file,"r")
  
  if (dycore .eq. "se") then
    ; if SE, the PHIS field that matches ncol can live in either PHIS or PHIS_d
    if (isfilevar(ttfile, "PHIS") .and. dimsizes(ttfile->PHIS) .eq. dim_sePS(0)) then
      print("SE: Using PHIS for topo adjustment")
      topo_model_SE = ttfile->PHIS     ;ncol  
    else if (isfilevar(ttfile, "PHIS_d") .and. dimsizes(ttfile->PHIS_d) .eq. dim_sePS(0)) then
      print("SE: Using PHIS_d for topo adjustment")
      topo_model_SE = ttfile->PHIS_d     ;ncol
    else
      print("SE: Cannot find valid PHIS or PHIS_d with dimensions that match ncol: "+ncol)
      exit
    end if
    end if
  else
    ; if not SE, we just assume user is smart and gave PHIS in correct dims
    print("Using PHIS for topo adjustment")
    topo_model_SE = ttfile->PHIS
  end if
  
  ; unpack FV
  if (dycore .eq. "fv") then
    print("unpacking fv vars")
    topo_data_SE := ndtooned(topo_data_SE)
    sfct_data_SE := ndtooned(sfct_data_SE)
    topo_model_SE := ndtooned(topo_model_SE)
    ps_fv := latlon_to_ncol(ps_fv)
    ; Unpack lev,lat,lon
    t_fv := latlon_to_ncol(t_fv)
    q_fv := latlon_to_ncol(q_fv)
    u_fv := latlon_to_ncol(u_fv)
    v_fv := latlon_to_ncol(v_fv)
    cldice_fv := latlon_to_ncol(cldice_fv)
    cldliq_fv := latlon_to_ncol(cldliq_fv)
    ncol=dim_sePS(0) * dim_sePS(1)
  else
    ncol=dim_sePS(0)
  end if
  
  ; check number of columns corrected
  vert_corrs = 0
  tcorriter=0
      
  do kk = 0,ncol-1
    if (.not.ismissing(sfct_data_SE(kk)))  ; check if missing (useful for regional data)
      ; get difference in geopotential
      deltaPhi = topo_data_SE(kk) - topo_model_SE(kk)
    
      ; estimate "model's" surface temperature
      Tsfc_fv = sfct_data_SE(kk) + gamma_s * (deltaPhi / grav)

      ; calculate layer mean temperature for use in hydrostatic
      lwr_coef = 0.5   ; use 0.5 for straight average
      if (Tsfc_fv .ge. sfct_data_SE(kk)) then
        Tlayermean = ( lwr_coef*Tsfc_fv + (1.-lwr_coef)*sfct_data_SE(kk) )
      else
        Tlayermean = ( (1.-lwr_coef)*Tsfc_fv + lwr_coef*sfct_data_SE(kk) )
      end if
    
      ; correct very warm and very cold layer means based on Trenberth 1993
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

      if (mod(kk,10000) .eq. 0) then
        print("Correcting PS: "+kk+" of "+(ncol-1)+"   from "+ps_orig+" to "+ps_fv(kk))
      end if
    
      ; correct T by shifting Tbot "down" the same delta
      if (tempadjustflag .eq. "a")    
        if (dycore .eq. "fv")
          if (mod(kk,10000) .eq. 0) then
            print("WARNING: tempadjust not supported for FV")
            print("WARNING: continuing!")
          end if
        else
          nlev=dimsizes(lev)
          t_orig = t_fv(nlev-1,kk)
          delT = Tsfc_fv - sfct_data_SE(kk)
          t_fv(nlev-1,kk) = t_fv(nlev-1,kk) + tofloat(delT)
          if (mod(kk,10000) .eq. 0) then
            print("Correcting TBOT: "+kk+" of "+(ncol-1)+"   from "+t_orig+" to "+t_fv(nlev-1,kk))
          end if
        end if
      end if
    
      ; correct other state variables!
      vert_interp_thresh = 0.1   ; ps corr diff (Pa) req. to interp vert profiles
      extrap_threshold = 5000.   ; maximum ps corr diff (Pa) to allow xtrp
      if ( abs(ps_orig-ps_fv(kk)) .gt. vert_interp_thresh) then
        pm_orig = hya*p0 + hyb*ps_orig
        pm_corr = hya*p0 + hyb*ps_fv(kk)
      
        if ( abs(ps_orig-ps_fv(kk)) .gt. extrap_threshold) then
          linlog = 2    ; no extrapolation outside of orig levs
        else
          linlog = -2   ; allow extrapolation
        end if
      
        tmp := int2p(pm_orig,t_fv(:,kk),pm_corr,linlog)
        t_fv(:,kk) = where(ismissing(tmp),t_fv(:,kk),tmp)

        tmp := int2p(pm_orig,q_fv(:,kk),pm_corr,linlog)
        q_fv(:,kk) = where(ismissing(tmp),q_fv(:,kk),tmp)

        tmp := int2p(pm_orig,u_fv(:,kk),pm_corr,linlog)
        u_fv(:,kk) = where(ismissing(tmp),u_fv(:,kk),tmp)
      
        tmp := int2p(pm_orig,v_fv(:,kk),pm_corr,linlog)
        v_fv(:,kk) = where(ismissing(tmp),v_fv(:,kk),tmp)
      
        ; we'll try to speed things up a bit here if we aren't writing cloud_vars
        if (add_cloud_vars) then
          tmp := int2p(pm_orig,cldice_fv(:,kk),pm_corr,linlog)
          cldice_fv(:,kk) = where(ismissing(tmp),cldice_fv(:,kk),tmp)
      
          tmp := int2p(pm_orig,cldliq_fv(:,kk),pm_corr,linlog)
          cldliq_fv(:,kk) = where(ismissing(tmp),cldliq_fv(:,kk),tmp)
        end if
        
        vert_corrs = vert_corrs+1
        delete(tmp)
      end if
    end if
  end do
  print("needed to correct "+vert_corrs+" vertical profiles for updated PS")
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

; dry pressure?
ps_wet_to_dry=False
if (ps_wet_to_dry) then
  print("wet to dry!")
  ncol=dim_sePS(0)
  pw_fv = ps_fv
  ps_orig = ps_fv
  pw_fv@long_name = "total column precipitable water"
  pw_fv@units     = "kg/m2"
  do kk = 0,ncol-1
    pi_orig = hyai*p0 + hybi*ps_fv(kk)
    nlevp1 = dimsizes(pi_orig)
    dp = tofloat(pi_orig(1:nlevp1-1) - pi_orig(0:nlevp1-2))
    pw_fv(kk) = prcwater_dp (q_fv(:,kk), dp )
    ps_fv(kk) = ps_orig(kk) - pw_fv(kk)*9.81
    if (mod(kk,10000) .eq. 0) then
      print(dp+"")
      print("Correcting PS: "+kk+" of "+(ncol-1)+"   from "+ps_orig(kk)+" to "+ps_fv(kk)+" since TPW: "+pw_fv(kk))
    end if
  end do

  delete(ps_orig)
  delete(pw_fv)
  delete(dp)
  delete(pi_orig)
  print("done!")
end if




; repack
if (dycore .eq. "fv") then
  print("repacking fv vars")
  ps_fv := ncol_to_latlon(ps_fv,dim_sePS(0),dim_sePS(1))
  t_fv := ncol_to_latlon(t_fv,dim_sePS(0),dim_sePS(1))
  q_fv := ncol_to_latlon(q_fv,dim_sePS(0),dim_sePS(1))
  u_fv := ncol_to_latlon(u_fv,dim_sePS(0),dim_sePS(1))
  v_fv := ncol_to_latlon(v_fv,dim_sePS(0),dim_sePS(1))
  cldice_fv := ncol_to_latlon(cldice_fv,dim_sePS(0),dim_sePS(1))
  cldliq_fv := ncol_to_latlon(cldliq_fv,dim_sePS(0),dim_sePS(1))
end if


if (dycore .eq. "se") then
  
  ncol=dim_sePS(0)
  
  print("Correcting time records and converting to "+write_type+"s")
  ps_fv_dbl_time     = add_time_define_precision(ps_fv,write_type,True)
  u_fv_dbl_time      = add_time_define_precision(u_fv,write_type,True)
  v_fv_dbl_time      = add_time_define_precision(v_fv,write_type,True)
  t_fv_dbl_time      = add_time_define_precision(t_fv,write_type,True)
  q_fv_dbl_time      = add_time_define_precision(q_fv,write_type,True)
  cldliq_fv_dbl_time = add_time_define_precision(cldliq_fv,write_type,True)
  cldice_fv_dbl_time = add_time_define_precision(cldice_fv,write_type,True)
  
  ; set coordinates and attributes for se lat/lon arrays
  selat!0 = "ncol"
  selon!0 = "ncol"
  selat@_FillValue = -900.
  selat@long_name = "latitude"
  selat@units = "degrees_north"
  selon@_FillValue = -900.
  selon@long_name = "longitude"
  selon@units = "degrees_east"

else if (dycore .eq. "fv") then

  nlat=dim_sePS(0)
  nlon=dim_sePS(1)
  
  print("FV: need to interpolate u/v to slat/slon")
  nslat=nlat-1
  nslon=nlon
  slat=new(nslat,"double")
  slon=new(nslon,"double")
  
  lat=todouble(fspan(-90.,90.,nlat))
  dellon=360./nlon
  lon=todouble(fspan(0,360.-dellon,nlon))
  
  slat(:)=todouble((lat(0:nlat-2)+lat(1:nlat-1))/2.0)
  slon = todouble(lon-(dellon/2.))
  
  print("FV: getting weights")
  w_stag=latRegWgt(slat,"double",0)
  
  print("FV: interpolating FV staggered")
  u_fv_slat = linint2(lon,lat,u_fv,True,lon,slat,0) ;(ntim,mlat,mlon)
  v_fv_slon = linint2(lon,lat,v_fv,True,slon,lat,0)
  
  print("Correcting time records and converting to "+write_type+"s")
  ps_fv_dbl_time     = add_time_define_precision(ps_fv,write_type,False)
  u_fv_dbl_time      = add_time_define_precision(u_fv,write_type,False)
  v_fv_dbl_time      = add_time_define_precision(v_fv,write_type,False)
  us_fv_dbl_time     = add_time_define_precision(u_fv_slat,write_type,False)
  vs_fv_dbl_time     = add_time_define_precision(v_fv_slon,write_type,False)
  t_fv_dbl_time      = add_time_define_precision(t_fv,write_type,False)
  q_fv_dbl_time      = add_time_define_precision(q_fv,write_type,False)
  cldliq_fv_dbl_time = add_time_define_precision(cldliq_fv,write_type,False)
  cldice_fv_dbl_time = add_time_define_precision(cldice_fv,write_type,False)
  
  ; Correct dim names for staggered coords
  us_fv_dbl_time!2 = "slat"
  vs_fv_dbl_time!3 = "slon"
  
  ; Adjust units
  cldice_fv_dbl_time@units="kg/kg"
  cldliq_fv_dbl_time@units="kg/kg"
  q_fv_dbl_time@units="kg/kg"
  
  lat!0="lat"
  lat&lat=lat
  lat@units="degrees_north"
  
  lon!0="lon"
  lon&lon=lon
  lon@units="degrees_east"
  
  slat!0="slat"
  slat&slat=slat
  slat@units="degrees_north"
  
  slon!0="slon"
  slon&slon=slon
  slon@units="degrees_east"

else
  print("unsupported dycore")
  exit
end if
end if

; clean up more stuff
delete([/ps_fv,t_fv,u_fv,v_fv,q_fv,cldliq_fv,cldice_fv/])
if (dycore .eq. "fv")
  delete([/v_fv_slon,v_fv_slon/])
end if

; Fix error where Q goes very large (missing values?)
qthresh=1.0
q_fv_dbl_time=where(q_fv_dbl_time.gt.qthresh,0.0,q_fv_dbl_time)
qthresh=1.0e-12
q_fv_dbl_time=where(q_fv_dbl_time.le.qthresh,0.0,q_fv_dbl_time)
if(any(ismissing(q_fv_dbl_time))) then
  print("Q is missing data...")
end if
if (datasource .ne. "RAP") then
  print("Filling q in global dataset")
  q_fv_dbl_time = where(ismissing(q_fv_dbl_time),0.0,q_fv_dbl_time)
end if

; create time array
timeunits  = "days since 1850-01-01 00:00:00"
time   = cd_inv_calendar(toint(yearstr),toint(monthstr),toint(daystr),toint(cyclestr),0,0,timeunits,0)
time!0 = "time"
delete(timeunits)

;===================================================================
; Remove any bad metadata we want here
;===================================================================

print("cleaning up output variables in FV/SE/HOMME/ncol part of atm_to_cam...")
vars_to_clean=[/ ps_fv_dbl_time, u_fv_dbl_time, v_fv_dbl_time, t_fv_dbl_time, q_fv_dbl_time /]
if (dycore .eq. "fv") then
  ListAppend(vars_to_clean,us_fv_dbl_time)
  ListAppend(vars_to_clean,vs_fv_dbl_time)
end if
if (add_cloud_vars) then
  ListAppend(vars_to_clean,cldliq_fv_dbl_time)
  ListAppend(vars_to_clean,cldice_fv_dbl_time)
end if
do cc = 0, ListCount(vars_to_clean)-1
  att_to_clean=(/"time", "lat", "slat", "lon", "slon", "lev", "ilev", "ncol"/)
  do jj = 0,dimsizes(att_to_clean)-1
    if(isatt(vars_to_clean[cc],att_to_clean(jj))) then
      print("var: "+cc+" -> cleaning "+att_to_clean(jj)+" attribute from state variables")
      delete_VarAtts(vars_to_clean[cc],att_to_clean(jj))
    end if
  end do
end do

;===================================================================
; Write output file
;===================================================================

print("Writing output file: "+se_inic)

if (isvar("compress_file") .and. compress_file) then
  compress_level=1
  print("writing netcdf4 classic with compression level "+compress_level)
  setfileoption("nc","Format","NetCDF4Classic")
  setfileoption("nc","CompressionLevel",compress_level)
else
  print("writing largefile with default compression settings")
  setfileoption("nc","Format","LargeFile")
end if

system("/bin/rm -f -v "+se_inic)   ; remove any pre-existing file
fout = addfile(se_inic,"c")  ; open output netCDF file   

;===================================================================
; explicitly declare file definition mode. Improve efficiency.
;===================================================================
setfileoption(fout,"DefineMode",True)

;===================================================================
; create global attributes of the file
;===================================================================
fAtt               = True            ; assign file attributes
fAtt@title         = "Betacast-generated ncdata file"  
fAtt@source_file   = data_filename
fAtt@wgt_file      = wgt_filename
fAtt@init_date     = YYYYMMDDHH
;fAtt@Conventions   = "None"   
fAtt@creation_date = systemfunc ("date")       
fAtt@dycore        = dycore
fAtt@datasource    = datasource
fileattdef( fout, fAtt )            ; copy file attributes    

;===================================================================
; predefine the coordinate variables and their dimensionality
;===================================================================

print("predefine the coordinate variables and their dimensionality")
nlev=dimsizes(lev)
nilev=nlev+1
if (dycore .eq. "fv") then
  dimNames = (/"time", "lat", "slat", "lon", "slon", "lev", "ilev"/)  
  dimSizes = (/ -1   ,  dimsizes(lat),dimsizes(slat),dimsizes(lon),dimsizes(slon), nlev, nilev /) 
  dimUnlim = (/ True , False, False, False, False, False, False/)
else
  dimNames = (/"time", "ncol", "lev", "ilev"/)  
  dimSizes = (/ -1   ,  ncol,  nlev, nilev /) 
  dimUnlim = (/ True , False, False, False/)
end if
filedimdef(fout,dimNames,dimSizes,dimUnlim)

;===================================================================
; predefine the the dimensionality of the variables to be written out
;===================================================================
print("predefine the the dimensionality of the variables to be written out")
filevardef(fout, "lev"   ,typeof(lev),getvardims(lev) )      
filevardef(fout, "ilev"  ,typeof(ilev),getvardims(ilev) )
filevardef(fout, "time"  ,typeof(time),getvardims(time) )    
if (dycore .eq. "fv") then
  filevardef(fout, "lat"  ,typeof(lat),getvardims(lat))                          
  filevardef(fout, "lon"  ,typeof(lon),getvardims(lon))
  filevardef(fout, "slat" ,typeof(slat),getvardims(slat))                          
  filevardef(fout, "slon" ,typeof(slon),getvardims(slon))
else                       
  filevardef(fout, "lat"  ,typeof(selat),getvardims(selat))                          
  filevardef(fout, "lon"  ,typeof(selon),getvardims(selon))
end if
filevardef(fout, "hyam"  ,typeof(hya),getvardims(hya))         
filevardef(fout, "hybm"  ,typeof(hyb),getvardims(hyb))         
filevardef(fout, "hyai"  ,typeof(hyai),getvardims(hyai))         
filevardef(fout, "hybi"  ,typeof(hyb),getvardims(hybi))      
filevardef(fout, "PS"    ,typeof(ps_fv_dbl_time)  ,getvardims(ps_fv_dbl_time))
if (dycore .eq. "fv") then
  filevardef(fout, "US"    ,typeof(us_fv_dbl_time)  ,getvardims(us_fv_dbl_time))    
  filevardef(fout, "VS"   ,typeof(vs_fv_dbl_time) ,getvardims(vs_fv_dbl_time))  
end if
filevardef(fout, "U"    ,typeof(u_fv_dbl_time)  ,getvardims(u_fv_dbl_time))    
filevardef(fout, "V"   ,typeof(v_fv_dbl_time) ,getvardims(v_fv_dbl_time))  
filevardef(fout, "T",typeof(t_fv_dbl_time),getvardims(t_fv_dbl_time))  
filevardef(fout, "Q",typeof(q_fv_dbl_time),getvardims(q_fv_dbl_time))
if (add_cloud_vars) then
  filevardef(fout, "CLDLIQ",typeof(cldliq_fv_dbl_time),getvardims(cldliq_fv_dbl_time))  
  filevardef(fout, "CLDICE",typeof(cldice_fv_dbl_time),getvardims(cldice_fv_dbl_time))  
end if

;====================================================================
; Copy attributes associated with each variable to the file
;====================================================================
print("Copy attributes associated with each variable to the file")
filevarattdef(fout,"lev",lev)
filevarattdef(fout,"ilev",ilev)
filevarattdef(fout,"time",time)
if (dycore .eq. "fv") then
  filevarattdef(fout,"lat",lat)
  filevarattdef(fout,"lon",lon)
  filevarattdef(fout,"slat",slat)
  filevarattdef(fout,"slon",slon)
else                       
  filevarattdef(fout,"lat",selat)
  filevarattdef(fout,"lon",selon)
end if
filevarattdef(fout,"hyam",hya)
filevarattdef(fout,"hybm",hyb)
filevarattdef(fout,"hyai",hyai)
filevarattdef(fout,"hybi",hybi)
filevarattdef(fout,"PS" ,ps_fv_dbl_time)
if (dycore .eq. "fv") then
  filevarattdef(fout,"US"  ,us_fv_dbl_time) 
  filevarattdef(fout,"VS"  ,vs_fv_dbl_time)
end if
filevarattdef(fout,"U"  ,u_fv_dbl_time) 
filevarattdef(fout,"V"  ,v_fv_dbl_time)
filevarattdef(fout,"T"  ,t_fv_dbl_time) 
filevarattdef(fout,"Q"   ,q_fv_dbl_time)
if (add_cloud_vars) then
  filevarattdef(fout,"CLDLIQ",cldliq_fv_dbl_time)
  filevarattdef(fout,"CLDICE",cldice_fv_dbl_time)
end if
;===================================================================
; explicitly exit file definition mode. **NOT REQUIRED**
;===================================================================
setfileoption(fout,"DefineMode",False)

;===================================================================
; output only the data values since the dimensionality and such have
; been predefined.
;====================================================================
print("output the data values")
fout->PS   = (/ps_fv_dbl_time/)  
if (dycore .eq. "fv") then 
  fout->US    = (/us_fv_dbl_time/)
  fout->VS    = (/vs_fv_dbl_time/)
end if
fout->U    = (/u_fv_dbl_time/)
fout->V    = (/v_fv_dbl_time/)
fout->T    = (/t_fv_dbl_time/) 
fout->Q      = (/q_fv_dbl_time/)
if (add_cloud_vars) then
  fout->CLDLIQ  = (/cldliq_fv_dbl_time/)
  fout->CLDICE  = (/cldice_fv_dbl_time/)
end if
fout->hyam   = (/hya/)     
fout->hybm    = (/hyb/)
fout->hyai    = (/hyai/)
fout->hybi    = (/hybi/) 
fout->lev     = (/lev/)
fout->ilev    = (/ilev/)
fout->time    = (/time/)
if (dycore .eq. "fv") then
  fout->lat   = (/lat/)
  fout->lon   = (/lon/)
  fout->slat  = (/slat/)
  fout->slon  = (/slon/)
else                       
  fout->lat   = (/selat/)
  fout->lon   = (/selon/)
end if

wallClockElapseTime(wcStrt, "Total runtime", 0)
print("done")

status_exit(9)

end


