load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"

begin

datestr = systemfunc("date +'%Y%m%d'")		; Get date string for parsing

rawhour = systemfunc("date -u +'%H'")

numlevels=30
yearstr="2012"
monthstr="10"
daystr="26"
cyclestr="12"

; GFS ARCHIVE analysis files are -- http://nomads.ncdc.noaa.gov/data/gfsanl/
;  grb_file = addfile("gfs.t06z.pgrb2f"+fcstTime(i)+".grib2","r")
;  grb_file = addfile("/glade/scratch/zarzycki/ISSAC_GFS/"+yearstr+monthstr+daystr+cyclestr+".grib2","r")
  sst_file = addfile("/glade/scratch/zarzycki/INIC/gdas1.t00z.sstgrb.grib2","r")  
  grb_file = addfile("/glade/scratch/zarzycki/INIC/gfsanl_4_20121026_1200_000.grb2","r")
;  grb_file = addfile("gfs.t"+cycle(1)+"z.pgrb2f"+fcstTime(i)+".grib2","r")
;	grb_file = addfile("ruc2.t00z.bgrb13f14.grib2","r")

sstlat = sst_file->lat_0
sstlon = sst_file->lon_0
sst_gfs = sst_file->TMP_P0_L1_GLL0(:,:)

grblat = grb_file->lat_0
grblon = grb_file->lon_0
grblev = grb_file->lv_ISBL0
rhlev = grb_file->lv_ISBL6
cldlev = grb_file->lv_ISBL8

ps = grb_file->PRES_P0_L1_GLL0(:,:)
p0   = 100000.

pblh = grb_file->HPBL_P0_L1_GLL0(:,:)

t_gfs = grb_file->TMP_P0_L100_GLL0(:,:,:)
u_gfs = grb_file->UGRD_P0_L100_GLL0(:,:,:)
v_gfs = grb_file->VGRD_P0_L100_GLL0(:,:,:)
rh_gfs = grb_file->RH_P0_L100_GLL0(:,:,:)
cldmix_gfs = grb_file->CLWMR_P0_L100_GLL0(:,:,:)

print("Interpolating GRIB to uniform vertical levels")
rh_gfs_interp = int2p_n(rhlev,rh_gfs,grblev,2,0)
cld_gfs_interp = int2p_n(cldlev,cldmix_gfs,grblev,2,0)

print("Calculating q from RH")
q_gfs = mixhum_ptrh(conform(t_gfs,grblev,0)*0.01, t_gfs , rh_gfs_interp, 2 ) ; specific humidity (g/kg)

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
cldmix_cam = pres2hybrid_Wrap(grblev,ps,p0,cld_gfs_interp,hya,hyb,4)

print("Sorting bad values")
 do i = 0,dimsizes(grblat)-1
    do j = 0,dimsizes(grblon)-1
        do k = 0,25
            if (ismissing(cldmix_cam(k,i,j))) then
                cldmix_cam(k,i,j) = 0
            else if (cldmix_cam(k,i,j) .gt. 0.01) then
                cldmix_cam(k,i,j) = 0
            end if
            end if
        end do
    end do
 end do

q_cam!0    = "lev"                      ; assign named dimensions
q_cam!1    = "lat"
q_cam!2    = "lon"
       
q_cam&lev = t_cam&lv_ISBL0
q_cam&lat = t_cam&lat_0
q_cam&lon = t_cam&lon_0

cldmix_cam!0    = "lev"                      ; assign named dimensions
cldmix_cam!1    = "lat"
cldmix_cam!2    = "lon"
       
cldmix_cam&lev = t_cam&lv_ISBL0
cldmix_cam&lat = t_cam&lat_0
cldmix_cam&lon = t_cam&lon_0

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

; Create sst array for CAM
sst_cam = sst_gfs
sst_cam!0    = "lat"
sst_cam!1    = "lon"

fV   = addfile ("./L26template.nc", "r")
fvlat  = fV->lat
fvlon  = fV->lon
fvslat = fV->slat
fvslon = fV->slon

print("Interpolating horizontal to CAM")

; Interpolate from GFS grid to CAM FV grid
ps_fv = linint2_Wrap (grblon,grblat(::-1),ps(::-1,:),True,fvlon,fvlat,0)
pblh_fv = linint2_Wrap (grblon,grblat(::-1),pblh(::-1,:),True,fvlon,fvlat,0)
t_fv = linint2_Wrap (grblon,grblat(::-1),t_cam(:,::-1,:),True,fvlon,fvlat,0)
u_fv = linint2_Wrap (grblon,grblat(::-1),u_cam(:,::-1,:),True,fvlon,fvslat,0)
v_fv = linint2_Wrap (grblon,grblat(::-1),v_cam(:,::-1,:),True,fvslon,fvlat,0)
q_fv = linint2_Wrap (grblon,grblat(::-1),q_cam(:,::-1,:),True,fvlon,fvlat,0)
cldmix_fv = linint2_Wrap (grblon,grblat(::-1),cldmix_cam(:,::-1,:),True,fvlon,fvlat,0)

sst_fv = linint2_Wrap (sstlon,sstlat(::-1),sst_cam(::-1,:),True,fvlon,fvlat,0)


u_fv!0    = "lev"                      ; assign named dimensions
u_fv!1    = "slat"
u_fv!2    = "lon"
v_fv!0    = "lev"                      ; assign named dimensions
v_fv!1    = "lat"
v_fv!2    = "slon"

print("Converting floats to doubles")
; Convert floats to doubles
ps_fv_dbl = todouble(ps_fv)
pblh_fv_dbl = todouble(pblh_fv)
t_fv_dbl = todouble(t_fv)
u_fv_dbl = todouble(u_fv)
v_fv_dbl = todouble(v_fv)
q_fv_dbl = todouble(q_fv)
cldmix_fv_dbl = todouble(cldmix_fv)
sst_fv_dbl = todouble(sst_fv)

print("Copying metadata")
; Copy meta-deta over to newly created doubles
copy_VarMeta(ps_fv,ps_fv_dbl)
copy_VarMeta(pblh_fv,pblh_fv_dbl)
copy_VarMeta(t_fv,t_fv_dbl)
copy_VarMeta(u_fv,u_fv_dbl)
copy_VarMeta(v_fv,v_fv_dbl)
copy_VarMeta(q_fv,q_fv_dbl)
copy_VarMeta(cldmix_fv,cldmix_fv_dbl)
copy_VarMeta(sst_fv,sst_fv_dbl)

; clean up some stuff
delete(ps_fv)
delete(pblh_fv)
delete(t_fv)
delete(u_fv)
delete(v_fv)
delete(q_fv)
delete(cldmix_fv)
delete(sst_fv)

numfvlat = dimsizes(fvlat)
numfvlon = dimsizes(fvlon)
numfvslat = dimsizes(fvslat)
numfvslon = dimsizes(fvslon)

print("Correcting time records")
; Add time record
ps_fv_dbl_time = new((/1,numfvlat,numfvlon/),double)
ps_fv_dbl_time(0,:,:) = ps_fv_dbl

pblh_fv_dbl_time = new((/1,numfvlat,numfvlon/),double)
pblh_fv_dbl_time(0,:,:) = pblh_fv_dbl

u_fv_dbl_reord = u_fv_dbl(lev|:,slat|:,lon|:) 
u_fv_dbl_time = new((/1,numlevels,numfvslat,numfvlon/),double)
u_fv_dbl_time(0,:,:,:) = u_fv_dbl_reord

v_fv_dbl_reord = v_fv_dbl(lev|:,lat|:,slon|:) 
v_fv_dbl_time = new((/1,numlevels,numfvlat,numfvslon/),double)
v_fv_dbl_time(0,:,:,:) = v_fv_dbl_reord

t_fv_dbl_reord = t_fv_dbl(lev|:,lat|:,lon|:) 
t_fv_dbl_time = new((/1,numlevels,numfvlat,numfvlon/),double)
t_fv_dbl_time(0,:,:,:) = t_fv_dbl_reord

q_fv_dbl_reord = q_fv_dbl(lev|:,lat|:,lon|:) 
q_fv_dbl_time = new((/1,numlevels,numfvlat,numfvlon/),double)
q_fv_dbl_time(0,:,:,:) = q_fv_dbl_reord

cldmix_fv_dbl_reord = cldmix_fv_dbl(lev|:,lat|:,lon|:) 
cldmix_fv_dbl_time = new((/1,numlevels,numfvlat,numfvlon/),double)
cldmix_fv_dbl_time(0,:,:,:) = cldmix_fv_dbl_reord

sst_fv_dbl_time = new((/1,numfvlat,numfvlon/),double)
sst_fv_dbl_time(0,:,:) = sst_fv_dbl

; Create clone ice array
cldice_fv_dbl_time = cldmix_fv_dbl_time



print("Seperate cloud ice and water")
; If T > 0 C, water, if less than < 0, ice
do i = 0,numfvlat-1
    do j = 0,numfvlon-1
        do k = 0,numlevels-1
            if (t_fv_dbl_time(0,k,i,j) .gt. 273.15 ) then
                cldice_fv_dbl_time(0,k,i,j) = 0
            else
                cldmix_fv_dbl_time(0,k,i,j) = 0
            end if
        end do
    end do
end do

delete(u_fv_dbl_reord)
delete(t_fv_dbl_reord)
delete(v_fv_dbl_reord)
delete(q_fv_dbl_reord)
delete(cldmix_fv_dbl_reord)

; Write to NetCDF

;  system("/bin/rm -f simple.nc")   ; remove any pre-existing file
;  sst_file     = addfile("sst.nc" ,"c")  ; open output netCDF file

; make time and UNLIMITED dimension      ; recommended  for most applications
;  filedimdef(sst_file,"time",-1,True)
  
;  sst_file->SST=(/sst_fv_dbl_time/)

  print("Writing file...")

  if (numlevels .eq. 30)
  in = addfile("/glade/scratch/zarzycki/gfs_to_cam/L30_clean.nc","w")
  else
  in = addfile("/glade/scratch/zarzycki/gfs_to_cam/L26_clean.nc","w")
  end if
 
  in->PS=(/ps_fv_dbl_time/)
  in->US=(/u_fv_dbl_time/)
  in->VS=(/v_fv_dbl_time/)
  in->T=(/t_fv_dbl_time/)
  in->Q=(/q_fv_dbl_time/)
  in->CLDLIQ=(/cldmix_fv_dbl_time/)
  in->CLDICE=(/cldice_fv_dbl_time/)
;  in->PBLH=(/pblh_fv_dbl_time/)
  

  if (numlevels .eq. 30)
;    system("ncks -O -x -v CLOUD,ICEFRAC,LCWAT,QCWAT,QPERT,TBOT,TCWAT,TPERT,TSICERAD,TSOCN fv.0.23x0.31_L26_input.nc fv.0.23x0.31_L26_input_clean.nc")
  else 
;    system("ncks -O -x -v CLOUD,ICEFRAC,LCWAT,QCWAT,QPERT,TBOT,TCWAT,TPERT,TSICERAD,TSOCN fv.0.23x0.31_L26_input.nc fv.0.23x0.31_L26_input_clean.nc")
  end if

  print("done")
  
end
