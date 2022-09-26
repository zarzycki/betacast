load "../../functions/fcns_seeding.ncl"
load "../../functions/radialAvg.ncl"
load "../../functions/namelist.ncl"

begin

;debug restart file vortex add
;system("cp -v /glade/u/home/zarzycki/scratch/APE3-test/run/ORIG.nc "+inic_file)

;pthi="../../namelists/unseed.default.nl"
invert_vortex  = keyword_values( pthi, "invert_vortex", "logical")
if (.not. invert_vortex)
  print("Not inverting vortex, we can exit gracefully!")
  status_exit(9)
end if
;inic_file  = keyword_values( pthi, "new_file", "string")
deltaMax  = keyword_values( pthi, "deltaMax", "float")
psminlat  = keyword_values( pthi, "psminlat", "float")
psminlon  = keyword_values( pthi, "psminlon", "float")
modify_q  = keyword_values( pthi, "modify_q", "logical")
modify_q_mult  = keyword_values( pthi, "modify_q_mult", "float")
gamma_  = keyword_values( pthi, "gamma", "float")
restart_file  = keyword_values( pthi, "restart_file", "logical")

if (gamma_ .lt. 0.0) then
  gamma_ = 0.0065
  print("setting gamma_ to default of "+gamma_)
end if

; ================== INTERNAL SETTINGS ===================

truePS_scan_radius=300.  ; what is radius to scan for true psmin relative to psminlat/psminlon? In km.
rad_for_corr=800.        ; what radius (in km) are we using for radial averages?

n_p_steps = 18
rp_arr = fspan(080000.,300000.,n_p_steps)
dp_arr = fspan(0200.,6000.,n_p_steps)
exppr_arr = fspan(1.1,1.9,n_p_steps)

n_t_steps = 60
zp_arr = fspan(6000.,16000.,n_t_steps)

n_p_vars = 3  ; rp, dp, exppr
n_t_vars = 1  ; zp

; ================== CONSTANTS ===================

pi      = 4.*atan(1.)       ; pi
convert = 180./pi   ; conversion factor: radians to degrees

; =================================================

; load a regridded state field with the TC in it.
f = addfile(inic_file,"r")
if (restart_file) then
  lat = f->lat_d
  lon = f->lon_d
  hybm = f->hybm
  hyam = f->hyam
  P0   = 100000.
  lev  = hyam*P0 + hybm*avg(f->PSDRY)
  lev  = lev/100.  ; needs to be in mb
else
  lat = f->lat
  lon = f->lon
  lev = f->lev
end if

ncol = dimsizes(lat)
nlev = dimsizes(lev)

if (restart_file) then
  ps = f->PSDRY
  corrMass = ps
  corrMass = 0.0
  corrMass = corrMass + dim_sum_n(f->dpQ,1)
  corrMass = corrMass + dim_sum_n(f->dpCLDICE,1)
  corrMass = corrMass + dim_sum_n(f->dpCLDLIQ,1)
  ps = ps + corrMass
  mps = ps(0,:)
  T = f->T(0,:,:)
  mps@units = "Pa"
  mps@long_name = "Dry surface pressure"
  T@units = "K"
  T@long_name = "Air temperature"
  delete(corrMass)
  delete(ps)
else
  mps = f->PS(0,:)
  T = f->T(0,:,:)
end if

; figure out ACTUAL PSL min location from our "best guess" (i.e., vitals)
gcdist = gc_latlon(psminlat,psminlon,lat,lon,2,4)
tmpps  = where(gcdist .gt. truePS_scan_radius,mps@_FillValue,mps)
minix  = minind(tmpps)
psminlat := lat(minix)
psminlon := lon(minix)
delete(gcdist)
delete(tmpps)

; Now create a radial profile of ncdata PSL given psminlat
rad_psl = radialAvg2D_unstruc(mps,lat,lon,deltaMax,psminlat,psminlon,rad_for_corr,False)
nrad_PS = dimsizes(rad_psl)

; Now create a radial profile of ncdata T given psminlat
Tavg = radialAvg3D_unstruc(T,lat,lon,lev,deltaMax,psminlat,psminlon,rad_for_corr,False)

; Gets the number of bins
Tavg_dims = dimsizes(Tavg)
nrad_T = Tavg_dims(1)

; Calculate the anomaly relative to the outermost bin
Tenv = Tavg(:,nrad_T-1)
Tanom = Tavg - conform(Tavg,Tenv,0)

; Rewrite Tavg as Tanom
Tavg := Tanom

; clean up
delete(T)
delete(mps)

; ===================OPTIMIZATION ROUND 1==============

wcOptimStrt     = systemfunc("date")
     
; Need to convert radius associated with radial averaging to degrees along line of constant lat
mylon = rad_psl&radius / (111.*cos(psminlat/convert))

; Get ambient pressure as outermost P from radial average
psl_amb = rad_psl(nrad_PS-1)

; create arrays for "analytic" PS and T profiles that match the radial profiles above.
anlPS=mylon
anlPS=100000.
anlT=new((/nlev,dimsizes(mylon)/),typeof(anlPS))

; Define cen_lat as where the TC is and cen_lon as just 0.0
cen_lat  = psminlat
cen_lon  = 0.0
zp       = 10000.   ; Initial guess for calculation of P
  
; Here, we don't want to invert vortex, we want to force ideal case to match obs.
invert_vortex=False

rp_ext    := ( rp_arr(n_p_steps-1)-rp_arr(0)       )
dp_ext    := ( dp_arr(n_p_steps-1)-dp_arr(0)       )
exppr_ext := ( exppr_arr(n_p_steps-1)-exppr_arr(0) )
print("orig rp range: "+(rp_ext)+"    orig dp range: "+(dp_ext)+"     orig exppr range: "+(exppr_ext))

; initialize "skill" array
rmsd = new((/(1+n_p_vars),toint(n_p_steps^n_p_vars)/),"float")

; Initialize linear counter
counter = 0

; Loop over PS vars of import
do aa = 0,dimsizes(rp_arr)-1
do bb = 0,dimsizes(dp_arr)-1
do cc = 0,dimsizes(exppr_arr)-1
  ; pull rp, dp and exppr from predefined earlier arrays
  rp = rp_arr(aa)
  dp = dp_arr(bb)
  exppr = exppr_arr(cc)

  ; use the "mylon" strip to calculate the analytic TC profile at each of the mylon longitudes along the cen_lat line. Also use kk index to estimate vertical profile
  do ii=0,dimsizes(mylon)-1
    tmp = tctestcase(cen_lon,cen_lat,dp,rp,zp,exppr,gamma_,mylon(ii),cen_lat,lev(nlev-1)*100.,-999,0,psl_amb,0.,0.,Tenv(nlev-1),0.0,invert_vortex,modify_q,modify_q_mult)
    anlPS(ii) = totype(tmp(4),typeof(anlPS))
  end do

  ;print("rp: "+rp+"   dp: "+dp+"    exppr: "+exppr+"        RMSD: "+dim_rmsd_n(anlPS,rad_psl,0))
  
  ; Load skill array
  rmsd(0,counter) = tofloat(dim_rmsd_n(anlPS,rad_psl,0))
  rmsd(1,counter) = rp
  rmsd(2,counter) = dp
  rmsd(3,counter) = exppr
  ; Increment linear counter
  counter = counter+1
end do
end do
end do

; Find "best" config from round 1
best_config = minind(rmsd(0,:))
print("BEST: "+rmsd(0,best_config)+" "+rmsd(1,best_config)+" "+rmsd(2,best_config)+" "+rmsd(3,best_config))

; ===================OPTIMIZATION ROUND 2==============

; "Shrink" the space by taking the delta between successive loops above and "widening" it a bit (mult_factor) for safety
; Note a bin width here = 2*mult_factor
n_p_steps=18
mult_factor=0.52
rp_ext    := ( rp_arr(1)-rp_arr(0)       )*mult_factor
dp_ext    := ( dp_arr(1)-dp_arr(0)       )*mult_factor
exppr_ext := ( exppr_arr(1)-exppr_arr(0) )*mult_factor
print("new rp range: "+(2.*rp_ext)+"    new dp range: "+(2.*dp_ext)+"     new exppr range: "+(2.*exppr_ext))

; Create updated arrays and reset skill and counter
rp_arr    := fspan(rmsd(1,best_config)-rp_ext,rmsd(1,best_config)+rp_ext,n_p_steps)
dp_arr    := fspan(rmsd(2,best_config)-dp_ext,rmsd(2,best_config)+dp_ext,n_p_steps)
exppr_arr := fspan(rmsd(3,best_config)-exppr_ext,rmsd(3,best_config)+exppr_ext,n_p_steps)
rmsd      := new((/(1+n_p_vars),toint(n_p_steps^n_p_vars)/),"float")
counter   := 0

do aa = 0,dimsizes(rp_arr)-1
do bb = 0,dimsizes(dp_arr)-1
do cc = 0,dimsizes(exppr_arr)-1
  ; pull rp, dp and exppr from predefined earlier arrays
  rp = rp_arr(aa)
  dp = dp_arr(bb)
  exppr = exppr_arr(cc)
  ; use the "mylon" strip to calculate the analytic TC profile at each of the mylon longitudes along the cen_lat line. Also use kk index to estimate vertical profile
  do ii=0,dimsizes(mylon)-1
    tmp = tctestcase(cen_lon,cen_lat,dp,rp,zp,exppr,gamma_,mylon(ii),cen_lat,lev(nlev-1)*100.,-999,0,psl_amb,0.,0.,Tenv(nlev-1),0.0,invert_vortex,modify_q,modify_q_mult)
    anlPS(ii) = totype(tmp(4),typeof(anlPS))
  end do
  ;print("rp: "+rp+"   dp: "+dp+"    exppr: "+exppr+"        RMSD: "+dim_rmsd_n(anlPS,rad_psl,0))
  rmsd(0,counter) = tofloat(dim_rmsd_n(anlPS,rad_psl,0))
  rmsd(1,counter) = rp
  rmsd(2,counter) = dp
  rmsd(3,counter) = exppr
  counter = counter+1
end do
end do
end do

best_config = minind(rmsd(0,:))
print("BEST: "+rmsd(0,best_config)+" "+rmsd(1,best_config)+" "+rmsd(2,best_config)+" "+rmsd(3,best_config))

; Save settings from "best" config in Round #2 to use in the T optimization
rp       = rmsd(1,best_config)
dp       = rmsd(2,best_config)
exppr    = rmsd(3,best_config)

; ===================OPTIMIZATION T ANOM==============

; initialize "skill" array
corr = new((/(1+n_t_vars),toint(n_t_steps^n_t_vars)/),"float")
; Initialize linear counter
counter = 0

; Loop over T vars
do dd = 0,dimsizes(zp_arr)-1

  ; pull zp
  zp = zp_arr(dd)

  do ii=0,dimsizes(mylon)-1
    do jj=0,nlev-1
      tmp = tctestcase(cen_lon,cen_lat,dp,rp,zp,exppr,gamma_,mylon(ii),cen_lat,lev(jj)*100.,-999,0,psl_amb,0.,0.,Tenv(jj),0.0,invert_vortex,modify_q,modify_q_mult)
      anlT(jj,ii) = totype(tmp(3),typeof(anlT))
    end do
  end do

  ; Calculate the anomaly from the TC test case by taking the T at the "center" of the vortex and T at the end of the lat strip
  anlT_anom = anlT(:,0) - anlT(:,dimsizes(mylon)-1)

  ;print("zp: "+zp+"          CORR: "+corr(0,counter))
  corr(0,counter) = tofloat(escorc(anlT_anom,Tavg(:,0)))
  corr(1,counter) = zp
  counter = counter+1
end do

best_config = maxind(corr(0,:))
print("BEST: "+corr(0,best_config)+" "+corr(1,best_config))
zp = corr(1,best_config)

; Print the wallclock time
wallClockElapseTime(wcOptimStrt, "Optimization Time", 0)

print("BEST SETTINGS --------------------------")
print("rp: "+rp)
print("dp: "+dp)
print("zp: "+zp)
print("gamma_: "+gamma_)
print("exppr: "+exppr)
print("cen_lat: "+psminlat)
print("cen_lat: "+psminlon)
print("modify_q: "+modify_q)
print("modify_q_mult: "+modify_q_mult)      

; write to namelist
WHICHSED="sed"
system("cp -v "+pthi+" "+pthi+".BAK")
;system(WHICHSED+" -i 's?.*rp=.*?rp="+rp+"?' "+pthi)
;system(WHICHSED+" -i 's?.*dp=.*?dp="+dp+"?' "+pthi)
;system(WHICHSED+" -i 's?.*zp=.*?zp="+zp+"?' "+pthi)
;system(WHICHSED+" -i 's?.*exppr=.*?exppr="+exppr+"?' "+pthi)
;system(WHICHSED+" -i 's?.*psminlat=.*?psminlat="+psminlat+"?' "+pthi)
;system(WHICHSED+" -i 's?.*psminlon=.*?psminlon="+psminlon+"?' "+pthi)

print("deleting vars from exiting namelist...")
; if we have added rp, dp, zp, or exppr, delete!
system(WHICHSED+" -i '/rp/d' "+pthi)
system(WHICHSED+" -i '/dp/d' "+pthi)
system(WHICHSED+" -i '/zp/d' "+pthi)
system(WHICHSED+" -i '/exppr/d' "+pthi)

print("forcing a newline at end of the file...")
; this line *forces* a newline at the end of the file if one does not exist so we can echo correctly
system(WHICHSED+" -i -e '$a\' "+pthi)

print("add new variables...")
; add the variables we have figured out!
system("echo 'rp = "+rp+"' >> "+pthi)
system("echo 'dp = "+dp+"' >> "+pthi)
system("echo 'zp = "+zp+"' >> "+pthi)
system("echo 'exppr = "+exppr+"' >> "+pthi)

print("overwrite ps cen lat and lon...")
; overwrite psminlat and psminlon since we needed those to begin with!
system(WHICHSED+" -i 's?.*psminlat=.*?psminlat="+psminlat+"?' "+pthi)
system(WHICHSED+" -i 's?.*psminlon=.*?psminlon="+psminlon+"?' "+pthi)
  
print("************* --------------------------")

; BEST CONFIG PROFILE CALCULATION

do ii=0,dimsizes(mylon)-1
  do jj=0,nlev-1
    tmp = tctestcase(cen_lon,cen_lat,dp,rp,zp,exppr,gamma_,mylon(ii),cen_lat,lev(jj)*100.,-999,0,psl_amb,0.,0.,Tenv(jj),0.0,invert_vortex,modify_q,modify_q_mult)
    anlT(jj,ii) = totype(tmp(3),typeof(anlT))
  end do
  anlPS(ii) = totype(tmp(4),typeof(anlPS))
end do
anlT_anom = anlT(:,0) - anlT(:,dimsizes(mylon)-1)


; ===================PLOTTING==============

wks   = gsn_open_wks ("png","best_fit_profiles")               ; send graphics to PNG file
plot = new(2,graphic)                          ; create a plot array

res                  = True                     ; plot mods desired
res@gsnDraw  = False                          ; don't draw
res@gsnFrame = False                          ; don't advance frame

res1=res
res2=res
res1@trYReverse = True 
res2@trYReverse = True 
res1@trXMinF = -10.
res1@trXMaxF = 10.
res1@trYMinF = 100.
res1@trYMaxF = 1010.
res2@trXMinF = res1@trXMinF
res2@trXMaxF = res1@trXMaxF
res2@trYMinF = res1@trYMinF
res2@trYMaxF = res1@trYMaxF
res1@xyLineColor = "red"
plot(0) = gsn_csm_x2y2(wks,anlT_anom,Tavg(:,0),lev,lev,res1,res2)

res3=res
res4=res
res3@trXMinF = 0.
res3@trXMaxF = 900.
res3@trYMinF = 95000.
res3@trYMaxF = 102000.
res4@trXMinF = res3@trXMinF
res4@trXMaxF = res3@trXMaxF
res4@trYMinF = res3@trYMinF
res4@trYMaxF = res3@trYMaxF
res3@xyLineColor           = "red"
plot(1) = gsn_csm_x2y2(wks,rad_psl&radius,rad_psl&radius,anlPS,rad_psl,res3,res4)

gsn_panel(wks,plot,(/2,1/),False)             ; now draw as one plot

print("PSL total RSMD: "+dim_rmsd_n(anlPS,rad_psl,0))

status_exit(9)

end 