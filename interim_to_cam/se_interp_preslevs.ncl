load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/esmf/ESMF_regridding.ncl"

begin

wcStrt = systemfunc("date")

numlevels=30
yearstr="2012"
monthstr="08"
daystr="24"
cyclestr="12"
YYYYMMDDHH=yearstr+monthstr+daystr+cyclestr
;ERA_dir="XXXXXX"

;ncl se_interp.ncl machineid=3 'gridname = "tcforecast_60_x4"' 'ERA_dir="/home/zarzycki/getECMWFdata/"'

;  echo "Sedding commands in NCL interpolation scripts"
;  cp se_interp.ncl $backupdir/se_interp.BAK.$timestamp
;  sed 's?^numlevels.*?numlevels='$numlevels'?' se_interp.ncl > test.ncl
;  mv test.ncl se_interp.ncl
;  sed 's?^monthstr.*?monthstr="'$monthstr'"?' se_interp.ncl > test.ncl
;  mv test.ncl se_interp.ncl
;  sed 's?^yearstr.*?yearstr="'$yearstr'"?' se_interp.ncl > test.ncl
;  mv test.ncl se_interp.ncl
;  sed 's?^daystr.*?daystr="'$daystr'"?' se_interp.ncl > test.ncl
;  mv test.ncl se_interp.ncl
;  sed 's?^cyclestr.*?cyclestr="'$cyclestr'"?' se_interp.ncl > test.ncl
;  mv test.ncl se_interp.ncl
;
;  ncl se_interp.ncl 'gridname = "'$gridname'"' 'gfs_filename = "'$gfs_files_path'/gfs_atm_'$yearstr$monthstr$daystr$cyclestr'.grib2"'
;  ncl sst_interp.ncl initdate=${yearstr}${monthstr}${daystr}${cyclestr} 'sst_file_full = "'$gfs_files_path'/gfs_sst_'$yearstr$monthstr$daystr$cyclestr'.grib2"'

;machineid = 3
;gridname = "tcforecast_60_x4"
gfs_filename = ERA_dir+"/ERA-Int_ml_"+YYYYMMDDHH+".nc"
sfc_filename = ERA_dir+"/ERA-Int_sfc_"+YYYYMMDDHH+".nc"

if (machineid .eq. 1) then
  gfs_data_dir = "/glade/scratch/zarzycki/sewx/GFS/"
;  sst_filename = "gdas1.t00z.sstgrb.grib2"
;  gfs_filename = "gfsanl_4_"+yearstr+monthstr+daystr+"_"+cyclestr+"00_000.grb2"
  wgt_file_dir = "/glade/p/work/zarzycki/maps/gfsmaps/"
  wgt_filename = "map_gfs0.50_TO_"+gridname+"_patc.nc"
  se_inic_dir = "/glade/scratch/zarzycki/sewx/INIC/"
  se_inic_filename = gridname+"_INIC.nc"
else if (machineid .eq. 2) then
  print("Flux not configured yet")
  exit
else if (machineid .eq. 3) then
  gfs_data_dir = "/home/zarzycki/getECMWFdata/"
;  sst_filename = "gdas1.t00z.sstgrb.grib2"
;  gfs_filename = "gfsanl_4_"+yearstr+monthstr+daystr+"_"+cyclestr+"00_000.grb2"
;  gfs_filename = "gfs.t"+cyclestr+"z.pgrb2bf00.grb2"
;  gfs_filename = "gfs.t12z.pgrb2f00.grb2"
  wgt_file_dir = "/home/zarzycki/getECMWFdata/"
  wgt_filename = "map_ERA_TO_"+gridname+"_patc.140530.nc"
;  wgt_file_dir = "/home/zarzycki/maps/wx_to_cam/"
;  wgt_filename = "map_gfs0.50_TO_"+gridname+"_patc.nc"
  se_inic_dir = "/home/zarzycki/sewx/INIC/"
  se_inic_filename = gridname+"_INIC.nc"
else
  print("Invalid Machine ID specified")
  exit
end if
end if
end if

wgt_file = wgt_file_dir+wgt_filename
se_inic = se_inic_dir+se_inic_filename

print("We are using 3-D file: "+gfs_filename)
print("We are using surface file: "+sfc_filename)
print("Using this remap: "+wgt_filename)

; Get strings
datestr = systemfunc("date +'%Y%m%d'")		; Get date string for parsing
rawhour = systemfunc("date -u +'%H'")

in = addfile(se_inic,"w")

grb_file = addfile(gfs_filename,"r")
sfc_file = addfile(sfc_filename,"r")

grblat = grb_file->latitude
grblon = grb_file->longitude
grblev = grb_file->level

p0   = 100000.

; We need to use short2flt because ERA-Interim data is stored
; in short format and requires scale/offset
ps = short2flt(sfc_file->sp(0,:,:))
pblh       = short2flt(sfc_file->sp(0,:,:))
t_gfs      = short2flt(grb_file->t(0,:,:,:))
u_gfs      = short2flt(grb_file->u(0,:,:,:))
v_gfs      = short2flt(grb_file->v(0,:,:,:))
q_gfs      = short2flt(grb_file->q(0,:,:,:))
cldice_gfs = short2flt(grb_file->ciwc(0,:,:,:))
cldliq_gfs = short2flt(grb_file->clwc(0,:,:,:))

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

; This was put here in the event we needed to regrid from ERA-interim to GFS
; although we use ESMF regridding, I'll just eave this here
; need to change above inputs to have a 2 after all of the var names (see below)
;
;print("interpolating ecmwf to gfs coords")
;gfslatlonFile = addfile("/home/zarzycki/sewx/GFS/gfs_atm_2012102100.grib2","r")
;gfslat = gfslatlonFile->lat_0
;gfslon = gfslatlonFile->lon_0
;ps = linint2_Wrap(grblon,grblat(::-1),ps2(::-1,:),True,gfslon,gfslat(::-1),0)
;t_gfs = linint2_Wrap(grblon,grblat(::-1),t_gfs2(:,::-1,:),True,gfslon,gfslat(::-1),0)
;u_gfs = linint2_Wrap(grblon,grblat(::-1),u_gfs2(:,::-1,:),True,gfslon,gfslat(::-1),0)
;v_gfs = linint2_Wrap(grblon,grblat(::-1),v_gfs2(:,::-1,:),True,gfslon,gfslat(::-1),0)
;q_gfs = linint2_Wrap(grblon,grblat(::-1),q_gfs2(:,::-1,:),True,gfslon,gfslat(::-1),0)
;cldmix_gfs = linint2_Wrap(grblon,grblat(::-1),cldmix_gfs2(:,::-1,:),True,gfslon,gfslat(::-1),0)
;pblh = linint2_Wrap(grblon,grblat(::-1),pblh2(::-1,:),True,gfslon,gfslat(::-1),0)

if (numlevels .eq. 30)
print("Loading 30 level data")
fC   = addfile ("./L30template.nc", "r")
else
print("Loading 26 level data")
fC   = addfile ("./L26template.nc", "r")
end if

hya  = fC->hyam                                 ; levels to which the variable
hyb  = fC->hybm                                 ; will be interpolated

print("Interpolating to CAM hybrid coordinates")
t_cam = pres2hybrid_Wrap(grblev,ps,p0,t_gfs,hya,hyb,4)
u_cam = pres2hybrid_Wrap(grblev,ps,p0,u_gfs,hya,hyb,4)
v_cam = pres2hybrid_Wrap(grblev,ps,p0,v_gfs,hya,hyb,4)
q_cam = pres2hybrid_Wrap(grblev,ps,p0,q_gfs,hya,hyb,4)
cldice_cam = pres2hybrid_Wrap(grblev,ps,p0,cldice_gfs,hya,hyb,4)
cldliq_cam = pres2hybrid_Wrap(grblev,ps,p0,cldliq_gfs,hya,hyb,4)

;print("Sorting bad values")
; do i = 0,dimsizes(grblat)-1
;    do j = 0,dimsizes(grblon)-1
;        do k = 0,25
;            if (ismissing(cldmix_cam(k,i,j))) then
;                cldmix_cam(k,i,j) = 0
;            else if (cldmix_cam(k,i,j) .gt. 0.01) then
;                cldmix_cam(k,i,j) = 0
;            end if
;            end if
;        end do
;    end do
; end do

q_cam!0    = "lev"                      ; assign named dimensions
q_cam!1    = "lat"
q_cam!2    = "lon"
       
;q_cam&lev = t_cam&level
;q_cam&lat = t_cam&latitude
;q_cam&lon = t_cam&longitude

;cldmix_cam!0    = "lev"                      ; assign named dimensions
;cldmix_cam!1    = "lat"
;cldmix_cam!2    = "lon"
       
;cldmix_cam&lev = t_cam&level
;cldmix_cam&lat = t_cam&latitude
;cldmix_cam&lon = t_cam&longitude

t_cam!0    = "lev"                      ; assign named dimensions
t_cam!1    = "lat"
t_cam!2    = "lon"
u_cam!0    = "lev"                      ; assign named dimensions
u_cam!1    = "lat"
u_cam!2    = "lon"
v_cam!0    = "lev"                      ; assign named dimensions
v_cam!1    = "lat"
v_cam!2    = "lon"
ps!0 = "lat"
ps!1 = "lon"
pblh!0 = "lat"
pblh!1 = "lon"

print("Regridding horizontal to CAM-SE")
dstFileName = se_inic             ; Destination grid
dfile       = addfile(dstFileName,"r")
lonCell     = dfile->lon
latCell     = dfile->lat

Opt         = True
Opt@CopyVarCoords = False
;Opt@Debug = True
ps_fv  = ESMF_regrid_with_weights(ps,wgt_file,Opt)
pblh_fv  = ESMF_regrid_with_weights(pblh,wgt_file,Opt)
t_fv  = ESMF_regrid_with_weights(t_cam,wgt_file,Opt)
u_fv  = ESMF_regrid_with_weights(u_cam,wgt_file,Opt)
v_fv  = ESMF_regrid_with_weights(v_cam,wgt_file,Opt)
q_fv  = ESMF_regrid_with_weights(q_cam,wgt_file,Opt)
cldice_fv  = ESMF_regrid_with_weights(cldice_cam,wgt_file,Opt)
cldliq_fv  = ESMF_regrid_with_weights(cldliq_cam,wgt_file,Opt)


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
 
; cldmix_fv_dbl = todouble(cldmix_fv)
; sst_fv_dbl = todouble(sst_fv)
; 
; print("Copying metadata")
; ; Copy meta-deta over to newly created doubles
 copy_VarMeta(ps_fv,ps_fv_dbl)
 copy_VarMeta(pblh_fv,pblh_fv_dbl)
 copy_VarMeta(t_fv,t_fv_dbl)
 copy_VarMeta(u_fv,u_fv_dbl)
 copy_VarMeta(v_fv,v_fv_dbl)
 copy_VarMeta(q_fv,q_fv_dbl)
 copy_VarMeta(cldliq_fv,cldliq_fv_dbl)
 copy_VarMeta(cldice_fv,cldice_fv_dbl)
; copy_VarMeta(sst_fv,sst_fv_dbl)
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
; delete(sst_fv)
; 
 print("Correcting time records")
; ; Add time record
 ps_fv_dbl_time = new((/1,ncol/),double)
 ps_fv_dbl_time(0,:) = ps_fv_dbl
; 
 pblh_fv_dbl_time = new((/1,ncol/),double)
 pblh_fv_dbl_time(0,:) = pblh_fv_dbl

; u_fv_dbl_reord = u_fv_dbl(lev|:,slat|:,lon|:) 
 u_fv_dbl_time = new((/1,numlevels,ncol/),double)
 u_fv_dbl_time(0,:,:) = u_fv_dbl
 
; v_fv_dbl_reord = v_fv_dbl(lev|:,lat|:,slon|:) 
 v_fv_dbl_time = new((/1,numlevels,ncol/),double)
 v_fv_dbl_time(0,:,:) = v_fv_dbl

;t_fv_dbl_reord = t_fv_dbl(lev|:,lat|:,lon|:) 
t_fv_dbl_time = new((/1,numlevels,ncol/),double)
t_fv_dbl_time(0,:,:) = t_fv_dbl

;q_fv_dbl_reord = q_fv_dbl(lev|:,lat|:,lon|:) 
q_fv_dbl_time = new((/1,numlevels,ncol/),double)
q_fv_dbl_time(0,:,:) = q_fv_dbl
;cldmix_fv_dbl_reord = cldmix_fv_dbl(lev|:,lat|:,lon|:) 
;cldmix_fv_dbl_time = new((/1,numlevels,ncol/),double)
;cldmix_fv_dbl_time(0,:,:) = cldmix_fv_dbl
;cldice_fv_dbl_time = cldmix_fv_dbl_time

cldliq_fv_dbl_time = new((/1,numlevels,ncol/),double)
cldliq_fv_dbl_time(0,:,:) = cldliq_fv_dbl

cldice_fv_dbl_time = new((/1,numlevels,ncol/),double)
cldice_fv_dbl_time(0,:,:) = cldice_fv_dbl

; print("Seperate cloud ice and water")
; ; If T > 0 C, water, if less than < 0, ice
; do i = 0,ncol-1
;     ;do j = 0,numfvlon-1
;         do k = 0,numlevels-1
;             if (t_fv_dbl_time(0,k,i) .gt. 273.15 ) then
;                 cldice_fv_dbl_time(0,k,i) = 0
;             else
;                 cldmix_fv_dbl_time(0,k,i) = 0
;             end if
;         end do
;     ;end do
; end do


print("Setting dimension names for new output fields")
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
in->PS=(ps_fv_dbl_time)
in->U=(u_fv_dbl_time)
in->V=(v_fv_dbl_time)
in->T=(t_fv_dbl_time)
in->Q=(q_fv_dbl_time)
in->CLDLIQ=(cldliq_fv_dbl_time)
in->CLDICE=(cldice_fv_dbl_time)
; ;  in->PBLH=(/pblh_fv_dbl_time/)
;   
; 
;   if (numlevels .eq. 30)
; ;    system("ncks -O -x -v CLOUD,ICEFRAC,LCWAT,QCWAT,QPERT,TBOT,TCWAT,TPERT,TSICERAD,TSOCN fv.0.23x0.31_L26_input.nc fv.0.23x0.31_L26_input_clean.nc")
;   else 
; ;    system("ncks -O -x -v CLOUD,ICEFRAC,LCWAT,QCWAT,QPERT,TBOT,TCWAT,TPERT,TSICERAD,TSOCN fv.0.23x0.31_L26_input.nc fv.0.23x0.31_L26_input_clean.nc")
;   end if

wallClockElapseTime(wcStrt, "Total runtime", 0)
print("done")
  
end

