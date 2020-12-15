; Inputs
; iniyear (1999)
; inimon (12)
; iniday (25)
; inisec (24000)
; filename

load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "/storage/home/cmz5202/sw/betacast/plotting_ncl/functions/rainratetodbz.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/esmf/ESMF_regridding.ncl"

begin

;ncl weatherplot.ncl inisec=43200 iniday=25 inimon=05 iniyear=2018 'filename="/glade/scratch/zarzycki/forecast_natlantic_30_x4_CAM5_L30/run/2018052512/_forecast_natlantic_30_x4_CAM5_L30.cam.h0.2018-05-27-00000.nc"' 'wgt_file="/glade/p/work/zarzycki/maps/forecast_plot_maps/map_natlantic_30_x4_to_0.25x0.25glob_bilinear.nc"'

; now read as command line option as of 11/17/17
;wgt_file = "/glade/p/work/zarzycki/maps/forecast_plot_maps/map_natlantic_30_x4_to_0.25x0.25glob_bilinear.nc"

iniyearstr = sprinti("%04d",iniyear)
inimonstr = sprinti("%02d",inimon)
inidaystr = sprinti("%02d",iniday)
inisecstr = sprinti("%05d",inisec)
inihourstr = sprinti("%02d",toint(inisec)/3600)
inihoursince = cd_inv_calendar(toint(iniyearstr),toint(inimonstr),toint(inidaystr),toint(inihourstr),0,0,"hours", 0)

out_type = "png"

;datafiledir = "/home/zarzycki/tcforecast_68_x4/run/test/proc/"
;fils = systemfunc("ls " + datafiledir + "*cam.h1*")
;numfiles = dimsizes(fils)

zz  = ispan(3,960,3)
zzs = sprinti("%0.3i",zz)   ; "000", "001", ... "100"

;do time = 0,numfiles-1

;filename = fils(time)

datafile = addfile(filename,"r")
;print("Using: "+fils(time))

strs = str_split(filename,".")
; strs(0) = path to grid
; strs(1) = model
; strs(2) = file (h1, etc)
; strs(3) = date
; strs(4) = extension

datestrs = str_split(strs(3),"-")
yearstr = datestrs(0)
monstr = datestrs(1)
daystr = datestrs(2)
secstr = datestrs(3)
hourstr = sprinti("%02d",toint(secstr)/3600)

hoursince = cd_inv_calendar(toint(yearstr),toint(monstr),toint(daystr),toint(hourstr),0,0,"hours", 0)
forecastHour = toint(hoursince - inihoursince)
forecastHourStr = sprinti("%03d",forecastHour)
print(forecastHourStr)

if (mod(forecastHour,6) .ne. 0) then
  print("For now, we are only printing data at 6 hr increments")
  exit
end if

print("Forecast for: "+hourstr+"Z "+monstr+" "+daystr+" "+yearstr)
print("Initialized at: "+inihourstr+"Z "+inimonstr+" "+inidaystr+" "+iniyearstr)

;---Get two variables to regrid and the source lat/lon grid
prect_se = datafile->PRECT            ; time (1) x nlev x ncol
;u250_se = datafile->U200
;v250_se = datafile->V200
u500_se = datafile->U500
v500_se = datafile->V500
z500_se = datafile->Z500
tmq_se = datafile->TMQ
flut_se = datafile->FLUT
u850_se = datafile->U850
v850_se = datafile->V850
mslp_se = datafile->PSL/100
preclav_se=datafile->PRECLav
preccav_se=datafile->PRECCav
precbsn_se=datafile->PRECBSN
ts_se=datafile->TREFHT
ptype_se=datafile->PTYPE

delete([/datafile/])

prectav_se=preclav_se
prectav_se=preclav_se+preccav_se
delete(preclav_se)
delete(preccav_se)

print(max(mslp_se))
print(min(mslp_se))

print("Begin interp")
Opt         = True
prect  = ESMF_regrid_with_weights(prect_se,wgt_file,Opt)
;u250  = ESMF_regrid_with_weights(u250_se,wgt_file,Opt)
;v250  = ESMF_regrid_with_weights(v250_se,wgt_file,Opt)
u500  = ESMF_regrid_with_weights(u500_se,wgt_file,Opt)
v500  = ESMF_regrid_with_weights(v500_se,wgt_file,Opt)
z500  = ESMF_regrid_with_weights(z500_se,wgt_file,Opt)
tmq  = ESMF_regrid_with_weights(tmq_se,wgt_file,Opt)
flut  = ESMF_regrid_with_weights(flut_se,wgt_file,Opt)
u850  = ESMF_regrid_with_weights(u850_se,wgt_file,Opt)
v850  = ESMF_regrid_with_weights(v850_se,wgt_file,Opt)
mslp  = ESMF_regrid_with_weights(mslp_se,wgt_file,Opt)
prectav= ESMF_regrid_with_weights(prectav_se,wgt_file,Opt)
precbsn= ESMF_regrid_with_weights(precbsn_se,wgt_file,Opt)
ts= ESMF_regrid_with_weights(ts_se,wgt_file,Opt)
ptype= ESMF_regrid_with_weights(ptype_se,wgt_file,Opt)

delete([/prect_se,u500_se,v500_se,z500_se,tmq_se,flut_se,u850_se,v850_se,mslp_se,prectav_se,precbsn_se,ts_se,ptype_se/])


guess     = 1                ; use zonal means
is_cyclic = True             ; cyclic [global]
nscan     = 500             ; usually much less than this
eps       = 1.e-2            ; variable dependent
relc      = 0.5              ; relaxation coefficient
opt       = 0                ; not used

;poisson_grid_fill( u250, is_cyclic, guess, nscan, eps, relc, opt)
;poisson_grid_fill( v250, is_cyclic, guess, nscan, eps, relc, opt)
poisson_grid_fill( u500, is_cyclic, guess, nscan, eps, relc, opt)
poisson_grid_fill( v500, is_cyclic, guess, nscan, eps, relc, opt)
poisson_grid_fill( prect, is_cyclic, guess, nscan, eps, relc, opt)
poisson_grid_fill( z500, is_cyclic, guess, nscan, eps, relc, opt)
poisson_grid_fill( tmq, is_cyclic, guess, nscan, eps, relc, opt)
poisson_grid_fill( flut, is_cyclic, guess, nscan, eps, relc, opt)
poisson_grid_fill( u850, is_cyclic, guess, nscan, eps, relc, opt)
poisson_grid_fill( v850, is_cyclic, guess, nscan, eps, relc, opt)
poisson_grid_fill( mslp, is_cyclic, guess, nscan, eps, relc, opt)
poisson_grid_fill( prectav, is_cyclic, guess, nscan, eps, relc, opt)
poisson_grid_fill( precbsn, is_cyclic, guess, nscan, eps, relc, opt)
poisson_grid_fill( ts, is_cyclic, guess, nscan, eps, relc, opt)
poisson_grid_fill( ptype, is_cyclic, guess, nscan, eps, relc, opt)

print(max(mslp))
print(min(mslp))

print("end interp")

dims=dimsizes(mslp)
ntime=dims(0)
nlat=dims(1)
nlon=dims(2)
lat = mslp&lat
lon = mslp&lon

;lat!0 = "lat"
;lon!0 = "lon"
;lat&lat = lat
;lon&lon = lon
;lev  = datafile->lev

;RAINDBZ = rainratetodbz(prect,1)

prect_mmday = prect
prect_mmday = prect_mmday*3600000
prect_mmday@units="mm/hr"

prectav_in = prectav
prectav_in = prectav_in*10800 ; seconds in 3 hours, now gives us meters per timestep
prectav_in = prectav_in*39.3701  ; now inches per timestep
;prectcum_in = dim_cumsum_n_Wrap(prectav_in,0)
;printVarSummary(prectcum_in)
prectav_in@units="inches"
prectav_in=where(prectav_in.lt.0.1,prectav_in@_FillValue,prectav_in)

precbsn_in = precbsn
precbsn_in = precbsn_in*10800 ; seconds in 3 hours, now gives us meters per timestep
precbsn_in = precbsn_in*39.3701  ; now inches per timestep
precbsn_in = precbsn_in*10  ; now snow water ratio
precbsn_in@units="inches"
precbsn_in=where(precbsn_in.lt.0.1,precbsn_in@_FillValue,precbsn_in)

;vort250 = uv2vrF_Wrap(u250,v250)
;vort250 = vort250/1e-5
;vort250@units="1E-5 1/s"

;vort500 = uv2vrF_Wrap(u500,v500)
;vort500 = vort500/1e-5
;vort500@units="1E-5 1/s"

vort500 = uv2vr_cfd (u500,v500,lat,lon, 2) 
copy_VarMeta(u850,vort500)
vort500 = vort500/1e-5
vort500@units="1E-5 1/s"
vort500@long_name="Vorticity"
;printVarSummary(vort500)

wind = u850 
wind = sqrt(u850^2 + v850^2)

tsf = ts
tsf = (ts-273.15)*9.0/5.0+32.0
tsf@units="degrees F"

ptype = toint(ptype)
copy_VarMeta(u850,ptype)
ptype@units="-"
ptype@long_name="Precipitation type"
ptype = where(prect_mmday .lt. 0.1,ptype@_FillValue,ptype)

;shear850250 = u850
;shear850250 = 1.943*sqrt((u850-u250)^2 + (v850-v250)^2)

;*****************************************************
; ANOMS
;*****************************************************

print("doing anoms")
f_anom=addfile("/storage/home/cmz5202/sw/betacast/plotting_ncl/climo-files/TS-1980-2009-NCEP.nc","r")
ANOMTIME=f_anom->time
ANOMDATE=cd_inv_calendar(2000,toint(monstr),toint(daystr),toint(hourstr),0,0,ANOMTIME@units,0)
TSANOM=f_anom->TS({ANOMDATE},:,:)
TSANOM=lonFlip(TSANOM)
TSANOM_CAM=linint2(TSANOM&lon,TSANOM&lat,TSANOM,True,tsf&lon,tsf&lat,0)
TSANOM_CAM_F=(TSANOM_CAM-273.15)*9.0/5.0+32.0
tsf_anom=tsf
tsf_anom(0,:,:)=tsf(0,:,:)-TSANOM_CAM_F(:,:)
delete([/f_anom,ANOMTIME,ANOMDATE,TSANOM,TSANOM_CAM,TSANOM_CAM_F/])


;*****************************************************
; PLOT
;*****************************************************

do plotindex = 0,2

; NORTH ATLANTIC
;if (plotindex .eq. 0)
;  plotname = "natl"
;  maxlat = 55.           ; NATL         
;  minlat = 5.
;  maxlon = -5. 
;  minlon = -100.
if (plotindex .eq. 0)
  plotname = "natl"
  maxlat = 53.           ; CONUS         
  minlat = 21.
  maxlon = 295.-360. 
  minlon = 231.-360.
;else if (plotindex .eq. 1)    ; ORIG EPAC
;  plotname = "epac"
;  maxlat = 40.             
;  minlat = 5.
;  maxlon = -90.
;  minlon = -165.
else if (plotindex .eq. 1)     ; Modif EPAC for NE USA
  plotname = "epac"
  maxlat = 49.              
  minlat = 37.
  maxlon = 295.-360. 
  minlon = 275.-360.
;else if (plotindex .eq. 1)     ; Modif EPAC for Caribbean
;  plotname = "epac"
;  maxlat = 30.              
;  minlat = 10.
;  maxlon = -55. 
;  minlon = -85.
;else if (plotindex .eq. 1)     ; Modif EPAC for Gulf
;  plotname = "epac"
;  maxlat = 35.              
;  minlat = 10.
;  maxlon = -70. 
;  minlon = -105.
;else if (plotindex .eq. 1)     ; Modif EPAC for SE USA
;  plotname = "epac"
;  maxlat = 37.              
;  minlat = 20.
;  maxlon = 289. 
;  minlon = 269.
;else if (plotindex .eq. 2)     ; for Gulf + Carrib
;  plotname = "float1"
;  maxlat = 33.              
;  minlat = 8.
;  maxlon = -58. 
;  minlon = -100.
;else if (plotindex .eq. 1)     ; Modif EPAC for CONUS
;  plotname = "epac"
;  maxlat = 53.              
;  minlat = 25.
;  maxlon = 300. 
;  minlon = 230.
;else if (plotindex .eq. 2)     ; Modif EPAC for Colorado
;  plotname = "float1"
;  maxlat = 41.5              
;  minlat = 36.5
;  maxlon = 258.5-360. 
;  minlon = 250.5-360.
;else if (plotindex .eq. 2)     ; Modif EPAC for WPAC
;  plotname = "float1"
;  maxlat = 40.0              
;  minlat = 5.0
;  maxlon = 170.0
;  minlon = 75.0
else if (plotindex .eq. 2)     ; Modif GLOB
  plotname = "float1"
  maxlat = 80.0
  minlat = -80.0
  maxlon = 180. 
  minlon = -180.
else
  plotname = "glob"
  maxlat = 90.           ; choose subregion           
  minlat = -90.
  maxlon = 180. 
  minlon = -180.
end if
end if
end if

print("Plotting: "+plotname)


;----------- Global resources -----------------------------------------

resbase                 = True
resbase@gsnDraw      = False                        ; don't draw yet
resbase@gsnFrame     = False                        ; don't advance frame yet
resbase@gsnMaximize = True
resbase@gsnLeftString = ""
resbase@gsnCenterString=""
resbase@gsnRightString = ""
resbase@tiMainString=""

;turn raster smoothing on
resbase@cnFillMode="RasterFill"
resbase@cnRasterSmoothingOn = True

res = resbase

res@gsnAddCyclic = False                    ; Don't add a cyclic point.

res@cnFillOn		  = True
res@cnLinesOn   = False 

res@mpFillOn                    = False        ; turn off gray fill
res@mpGeophysicalLineColor      = "Black"      ; color of cont. outlines
res@mpGeophysicalLineThicknessF = 3            ; thickness of outlines

res@mpMaxLatF                   = maxlat           ; choose subregion           
res@mpMinLatF                   = minlat
res@mpMinLonF = minlon 
res@mpMaxLonF = maxlon 
res@mpCenterLonF = (res@mpMinLonF+res@mpMaxLonF)/2.

; this adds county level data
if (plotindex .eq. 9)
  res@mpOutlineBoundarySets = "AllBoundaries"
  res@mpDataBaseVersion         = "mediumres"             ; select database
  res@mpDataSetName             = "Earth..4"
  res@mpUSStateLineThicknessF = 4.0
else
  res@mpOutlineBoundarySets     = "geophysicalandusstates"; turn on states
  res@mpDataBaseVersion         = "mediumres"             ; select database
  res@mpDataSetName             = "Earth..2"
end if

; OLD: removed 12/5/17
;res@tiMainString = "+"+forecastHourStr
;res@gsnCenterString = "VALID "+hourstr+"Z "+monstr+" "+daystr+" "+yearstr
; NEW: as of 12/15/17
res@gsnLeftString = "VALID "+hourstr+"Z "+monstr+" "+daystr+" "+yearstr+" [+"+forecastHourStr+"]"

res@tmXBMajorLengthF=0.0005
res@tmXBMajorOutwardLengthF=0.0005
res@tmXTMajorLengthF=res@tmXBMajorLengthF
res@tmXTMajorOutwardLengthF=res@tmXBMajorOutwardLengthF

res@tmYLMajorLengthF=0.0005
res@tmYLMajorOutwardLengthF=0.0005
res@tmYRMajorLengthF=res@tmYLMajorLengthF
res@tmYRMajorOutwardLengthF=res@tmYLMajorOutwardLengthF

res@tmXBLabelFontHeightF=0.008
res@tmXBLabelDeltaF=0.05

res@lbLabelFontHeightF = 0.01
res@pmLabelBarHeightF = 0.05
res@pmLabelBarWidthF = 0.6

out_type = "png"
; previously was 900x900 before 9/12/17
out_type@wkWidth = 825
out_type@wkHeight = 825

res_contour = True
res_contour@gsnDraw          = False           ; don't draw
res_contour@gsnFrame         = False           ; don't advance frame
res_contour@cnFillOn		  = False
res_contour@cnLinesOn   = True
res_contour@gsnRightString   = ""    
res_contour@gsnLeftString    = ""
res_contour@gsnCenterString  = ""
res_contour@cnInfoLabelOn = False

res_mslp = res_contour
if (plotindex .eq. 0 .or. plotindex .eq. 1)
  res_mslp@cnLowLabelsOn = True
else
  res_mslp@cnLowLabelsOn = False
end if
res_mslp@cnLowLabelBackgroundColor = -1
res_mslp@cnLabelMasking = False
res_mslp@cnLowLabelFontColor = "red"
res_mslp@cnLineLabelsOn = False
res_mslp@cnLevelSelectionMode = "ManualLevels"    ; set contour intervals
res_mslp@cnLevelSpacingF      =  4.0                ;
res_mslp@cnMinLevelValF       =  900.0           ;
res_mslp@cnMaxLevelValF       =  1060.0 ; Use 80.0 -> 60 for this

res_z500 = res_contour
res_z500@cnLevelSelectionMode = "ManualLevels"    ; set contour intervals
res_z500@cnLevelSpacingF      =  60.0                ;
res_z500@cnMinLevelValF       =  4800.0           ;
res_z500@cnMaxLevelValF       =  6200.0 ; Use 80.0 -> 60 for this
res_z500@gsnRightString   = ""    
res_z500@gsnLeftString    = ""
res_z500@gsnCenterString  = ""

; =============================================================================================================================================

print("Plotting instant precip...")

wks = gsn_open_wks(out_type,"prect_"+plotname+"_"+yearstr+monstr+daystr+hourstr)

;;;; RAIN, SIMULATED RADAR
; res_rain                      = res
; gsn_define_colormap(wks,"rainsnowradar")          ; choose colormap
; res_rain@gsnRightString = "Simulated Reflectivity (dBZ)"
; res_rain@gsnSpreadColors      = False              ; use full colormap
; res_rain@cnFillColors = (/-1,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32/)
; res_rain@cnLevelSelectionMode = "ExplicitLevels"    			; set explicit contour levels
; res_rain@cnLevels    = (/0,5,10,15,20,25,30,35,40,45,50,55,60,65,70/)	; set levels
; res_rain@lbLabelStride = 2
; res_rain@lbOrientation = "vertical"
; plotrain = gsn_csm_contour_map(wks,RAINDBZ(0,:,:),res_rain)
; draw(plotrain)
; frame(wks)
; delete(wks)

;;;; RAIN, MM/DAY
res_rain                      = res
cmap = read_colormap_file("prcp_1")
cmap_dims = dimsizes(cmap)
ncolors = cmap_dims(0)
npairs = cmap_dims(1)
nnewcolors = 200
cmap_interp = new((/nnewcolors,npairs/),typeof(cmap))
do ii = 0,npairs-1
  ;cmap_interp(3:nnewcolors-1,ii) = linint1 (fspan(0,1,ncolors-3), cmap(3:ncolors-1,ii), False, fspan(0,1,nnewcolors-3), 0)
  ;cmap_interp(2:nnewcolors-1,ii) = linint1 (fspan(0,1,ncolors-2), cmap(2:ncolors-1,ii), False, fspan(0,1,nnewcolors-2), 0)
  cmap_interp(:,ii) = linint1 (fspan(0,1,ncolors), cmap(:,ii), False, fspan(0,1,nnewcolors), 0)
end do
cmap_interp(0,:) = cmap(0,:)
cmap_interp(1,:) = cmap(1,:)
;cmap_interp((nnewcolors-1),:) = cmap((ncolors-1),:)
;cmap_interp((nnewcolors-2),:) = cmap((ncolors-2),:)

res_rain@cnFillPalette=cmap

res_rain@cnLevelSelectionMode = "ExplicitLevels"    			; set explicit contour levels
res_rain@cnLevels    = (/0.1,0.25,0.5,1.0,2.0,3.0,4.0,5.0,7.5,10.0,15.0,20.0,25.0,35.0/)	; set levels

res_rain@gsnRightString = "Instant. precip. rate (color, mm/hr), MSLP (black contour, hPa)"

plotrain = gsn_csm_contour_map(wks,prect_mmday(0,:,:),res_rain)

; overlay
plotmslp = gsn_csm_contour(wks,mslp(0,:,:),res_mslp)
overlay(plotrain,plotmslp)

draw(plotrain)
frame(wks)

delete(wks)
delete(plotmslp)
delete(cmap)
delete(cmap_interp)

; =============================================================================================================================================

print("Plotting ptype...")

wks = gsn_open_wks(out_type,"ptype_"+plotname+"_"+yearstr+monstr+daystr+hourstr)

;;;; ptype, MM/DAY
res_rainptype = res
cmap = read_colormap_file("MPL_Greys")
res_rainptype@cnFillPalette=cmap

res_rainptype@cnLevelSelectionMode = "ExplicitLevels"    			; set explicit contour levels
res_rainptype@cnLevels    = (/0.1,0.15,0.25,0.5,0.75,1.0,2.0,3.0,4.0,5.0,6.0,8.0,10.0,15.0/)	; set levels

res_rainptype@gsnRightString = "Instant. precip. rate (shading, mm/hr), precip. type (color)"

plotptype = gsn_csm_contour_map(wks,prect_mmday(0,:,:),res_rainptype)

res_ptype = resbase
res_ptype@cnFillOn		  = True
res_ptype@cnLinesOn   = False
res_ptype@lbLabelBarOn = False
res_ptype@cnInfoLabelOn = False

delete(cmap)
cmap = (/"white","dodgerblue","orange","green","red","white"/)
res_ptype@cnFillPalette=cmap
res_ptype@cnLevelSelectionMode = "ExplicitLevels"    			; set explicit contour levels
res_ptype@cnLevels    = (/0,1,2,3,4/)	; set levels
res_ptype@cnFillOpacityF=0.2

plotptype2 = gsn_csm_contour(wks,ptype(0,:,:),res_ptype)
overlay(plotptype,plotptype2)

; overlay
plotmslp = gsn_csm_contour(wks,mslp(0,:,:),res_mslp)
overlay(plotptype,plotmslp)

draw(plotptype)
frame(wks)

delete(wks)
delete(cmap)

delete(plotmslp)
delete(plotptype2)

; =============================================================================================================================================

print("Plotting accumulated precip...")
wks = gsn_open_wks(out_type,"sumprect_"+plotname+"_"+yearstr+monstr+daystr+hourstr)
res_sumprecip                    = res

cmap = read_colormap_file("precip2_17lev")
cmap_dims = dimsizes(cmap)
ncolors = cmap_dims(0)
npairs = cmap_dims(1)
nnewcolors = 200
cmap_interp = new((/nnewcolors,npairs/),typeof(cmap))
do ii = 0,npairs-1
  ;cmap_interp(3:nnewcolors-1,ii) = linint1 (fspan(0,1,ncolors-3), cmap(3:ncolors-1,ii), False, fspan(0,1,nnewcolors-3), 0)
  ;cmap_interp(2:nnewcolors-1,ii) = linint1 (fspan(0,1,ncolors-2), cmap(2:ncolors-1,ii), False, fspan(0,1,nnewcolors-2), 0)
  cmap_interp(:,ii) = linint1 (fspan(0,1,ncolors), cmap(:,ii), False, fspan(0,1,nnewcolors), 0)

end do
cmap_interp(0,:) = cmap(0,:)
cmap_interp(1,:) = cmap(1,:)
;cmap_interp((nnewcolors-1),:) = cmap((ncolors-1),:)
;cmap_interp((nnewcolors-2),:) = cmap((ncolors-2),:)

res_sumprecip@cnFillPalette=cmap_interp

res_sumprecip@cnLevelSelectionMode = "ExplicitLevels"    			; set explicit contour levels
res_sumprecip@cnLevels    = (/0.05,0.1,0.15,0.2,0.3,0.4,0.5,0.75,1,1.25,1.5,1.75,2,2.5,3,3.5,4,4.5,5,5.5,6,7,8,9,10,12,14,16,20,24/)	; set levels
res_sumprecip@lbLabelStride=1
res_sumprecip@lbLabelAngleF=50                 ; angle labels

res_sumprecip@gsnRightString = "Accum. precipitation (in.)"

plotprecipsum = gsn_csm_contour_map(wks,prectav_in(0,:,:),res_sumprecip)
draw(plotprecipsum)
frame(wks)

delete(wks)
delete(cmap)
delete(cmap_interp)

; =============================================================================================================================================

print("Plotting accumulated snow...")

wks = gsn_open_wks(out_type,"sumsnow_"+plotname+"_"+yearstr+monstr+daystr+hourstr)
res_sumsnow   = res

cmap = read_colormap_file("/storage/home/cmz5202/sw/betacast/plotting_ncl/colormaps/wxbell_snow.rgb")
res_sumsnow@cnFillPalette = cmap(::-1,:)       ; reverse color map
res_sumsnow@cnFillOpacityF = 0.75

res_sumsnow@cnLevelSelectionMode = "ExplicitLevels"    			; set explicit contour levels
res_sumsnow@cnLevels    = (/ 0.25, 0.5, 1., 1.5, 2., 3., 4., 5., 6., 7., 8., 9., 10., 11., 12., 13., 14., 15., 16., 17., 18., 19., 20., 21., 22., 23., 24., 26., 28., 30., 32., 34., 36., 38., 40., 42., 44., 46., 48., 52., 56., 60., 64., 68., 72./)

res_sumsnow@gsnRightString = "Accum. Bourgouin snow (10:1 in.)"

plotsnowsum = gsn_csm_contour_map(wks,precbsn_in(0,:,:),res_sumsnow)

txres               = True  
txres@txFontHeightF = 0.011             ; Set the font height
txres@txFont   = "helvetica-bold"
text = new((/nlat,nlon/),"graphic")

if (plotindex .eq. 0) then
  stride=18
else if (plotindex .eq. 1) then
  stride=6
else
  stride=60
end if
end if

do ii = 0,nlat-1,stride
do jj = 0,nlon-1,stride
  if (lat(ii) .ge. minlat .and. lat(ii) .le. maxlat .and. lon(jj) .ge. minlon .and. lon(jj) .le. maxlon) then
    if (.not.ismissing(precbsn_in(0,ii,jj)) .and. precbsn_in(0,ii,jj) .ge. 1.0) then
      text(ii,jj) = gsn_add_text(wks,plotsnowsum,tostring(toint(precbsn_in(0,ii,jj))),lon(jj),lat(ii),txres)
    end if
  end if
end do
end do

draw(plotsnowsum)
frame(wks)
delete(wks)
delete(cmap)

delete(text)
delete(stride)
delete(txres)


; =============================================================================================================================================

; >>>out_type = "png"
; >>>wks = gsn_open_wks(out_type,"250vort_"+plotname+"_"+yearstr+monstr+daystr+hourstr)
; >>>gsn_define_colormap(wks,"BlueWhiteOrangeRed")          ; choose colormap
; >>>res_250vort = res
; >>> res_250vort@cnLevelSelectionMode = "ManualLevels"    ; set contour intervals
; >>> res_250vort@cnLevelSpacingF      =  2.5              ;
; >>> res_250vort@cnMinLevelValF       =  -40.0           ;
; >>> res_250vort@cnMaxLevelValF       =  40.0 ; Use 80.0 -> 60 for this
; >>>
; >>>plot250vort = gsn_csm_contour_map(wks,vort250(0,:,:),res_250vort)
; >>>draw(plot250vort)
; >>>frame(wks)
; >>>delete(wks)

; =============================================================================================================================================

print("Plotting 500 mb vort + z500 + U/V 500...")


wks = gsn_open_wks(out_type,"500vort_"+plotname+"_"+yearstr+monstr+daystr+hourstr)

res_500vort = res

cmap = read_colormap_file("WhiteYellowOrangeRed")
res_500vort@cnFillPalette=cmap
;gsn_define_colormap(wks,"WhiteYellowOrangeRed")          ; choose colormap

res_500vort@cnLevelSelectionMode = "ManualLevels"    ; set contour intervals
res_500vort@cnLevelSpacingF      =  2.0              ;
res_500vort@cnMinLevelValF       =  0.0           ;
res_500vort@cnMaxLevelValF       =  40.0 ; Use 80.0 -> 60 for this

res_500vort@gsnRightString = "500mb vort/Z/wind"


plot500vort = gsn_csm_contour_map(wks,vort500(0,:,:),res_500vort)
plotz500 = gsn_csm_contour(wks,z500(0,:,:),res_z500)

vecres                  = True            ; vector only resources
vecres@gsnDraw          = False           ; don't draw
vecres@gsnFrame         = False           ; don't advance frame
vecres@vcGlyphStyle     = "CurlyVector"   ; curly vectors
vecres@vcRefMagnitudeF  = 20              ; define vector ref mag
vecres@vcRefLengthF     = 0.045           ; define length of vec ref
vecres@gsnRightString   = ""             ; turn off right string
vecres@gsnLeftString    = ""             ; turn off left string
vecres@tiXAxisString    = ""             ; turn off axis label 
;vecres@vcRefAnnoOrthogonalPosF = -.535    ; move ref vector into plot

; this is 8 for 0.25deg data, 16 for 0.125deg data

if (plotindex .eq. 0) then
  stride=16
else if (plotindex .eq. 1) then
  stride=16
else
  stride=40
end if
end if

plot500winds   = gsn_csm_vector(wks,u500(0,::stride,::stride),v500(0,::stride,::stride),vecres)
overlay(plot500vort,plotz500)
overlay(plot500vort,plot500winds)

draw(plot500vort)
frame(wks)
delete(wks)
delete(cmap)

; ====== Plot Z500 ============

; >>>out_type = "png"
; >>>wks = gsn_open_wks(out_type,"z500_"+plotname+"_"+yearstr+monstr+daystr+hourstr)
; >>>
; >>>res_z500 = res
; >>> res_z500@cnLevelSelectionMode = "ManualLevels"    ; set contour intervals
; >>>; res_z500@cnLevelSpacingF      =  2.5              ;
; >>>; res_z500@cnMinLevelValF       =  -40.0           ;
; >>>; res_z500@cnMaxLevelValF       =  40.0 ; Use 80.0 -> 60 for this
; >>>
; >>>plotz500 = gsn_csm_contour_map(wks,z500(0,:,:),res_z500)
; >>>draw(plotz500)
; >>>frame(wks)
; >>>delete(wks)

; ====== Plot TMQ ============

print("Plotting TMQ...")

wks = gsn_open_wks(out_type,"tmq_"+plotname+"_"+yearstr+monstr+daystr+hourstr)

res_tmq = res

cmap = read_colormap_file("nice_gfdl")
res_tmq@cnFillPalette=cmap

res_tmq@cnLevelSelectionMode = "ManualLevels"    ; set contour intervals
res_tmq@cnLevelSpacingF      =  2.0              ;
res_tmq@cnMinLevelValF       =  20.0           ;
res_tmq@cnMaxLevelValF       =  66.0 ; Use 80.0 -> 60 for this

res_tmq@gsnRightString = "Total precipitable water (mm)"

plottmq = gsn_csm_contour_map(wks,tmq(0,:,:),res_tmq)
draw(plottmq)
frame(wks)

delete(wks)
delete(cmap)

; ====== Plot TS ============

print("Plotting TS...")

wks = gsn_open_wks(out_type,"ts_"+plotname+"_"+yearstr+monstr+daystr+hourstr)

res_ts = res

cmap = read_colormap_file("amwg256")
res_ts@cnFillPalette=cmap

res_ts@cnLevelSelectionMode = "ManualLevels"    ; set contour intervals
res_ts@cnLevelSpacingF      =  2.5             ;
res_ts@cnMinLevelValF       =  -5.0           ;
res_ts@cnMaxLevelValF       =  110.0 ; Use 80.0 -> 60 for this

res_ts@gsnRightString = "Surface temperature (degF)"

plotts = gsn_csm_contour_map(wks,tsf(0,:,:),res_ts)

res_fzTS=res_contour
res_fzTS@cnLevelSelectionMode = "ExplicitLevels"     ; set manual contour levels
res_fzTS@cnLevels    = (/ 32. /)      ; set levels
res_fzTS@cnLineLabelsOn       = False
res_fzTS@cnInfoLabelOn = False
res_fzTS@cnLineColor="red3"
res_fzTS@cnLineThicknessF = 4.             ; doubles thickness of zero contour
plotzeroline = gsn_csm_contour(wks,tsf(0,::2,::2),res_fzTS)
overlay(plotts,plotzeroline)

txres               = True  
txres@txFontHeightF = 0.011             ; Set the font height
txres@txFont   = "helvetica-bold"
text = new((/nlat,nlon/),"graphic")

if (plotindex .eq. 0) then
  stride=25
else if (plotindex .eq. 1) then
  stride=8
else
  stride=120
end if
end if

;if (plotindex .eq. 0) then
do ii = 0,nlat-1,stride
do jj = 0,nlon-1,stride
  if (lat(ii) .ge. minlat .and. lat(ii) .le. maxlat .and. lon(jj) .ge. minlon .and. lon(jj) .le. maxlon) then
    if (.not.ismissing(tsf(0,ii,jj))) then
      text(ii,jj) = gsn_add_text(wks,plotts,tostring(toint(tsf(0,ii,jj))),lon(jj),lat(ii),txres)
    end if
  end if
end do
end do
;end if

draw(plotts)
frame(wks)

delete(wks)
delete(cmap)

; ====== Plot TS ANOM ============

print("Plotting TS ANOM...")

wks = gsn_open_wks(out_type,"tsanom_"+plotname+"_"+yearstr+monstr+daystr+hourstr)

res_tsanom = res

cmap = read_colormap_file("MPL_RdBu")
res_tsanom@cnFillPalette=cmap(::-1,:)

res_tsanom@cnLevelSelectionMode = "ManualLevels"    ; set contour intervals
res_tsanom@cnLevelSpacingF      =  2.0             ;
res_tsanom@cnMinLevelValF       =  -20.0           ;
res_tsanom@cnMaxLevelValF       =  20.0 ; Use 80.0 -> 60 for this

res_tsanom@gsnRightString = "Sfc. temp. anom. (1980-2009 NCEP base) (degF)"

;plotts = gsn_csm_contour_map(wks,tsf_anom(0,:,:),res_tsanom)
tsf_tmp=area_hi2lores_Wrap(tsf_anom&lon,tsf_anom&lat,tsf_anom,True,1,tsf_anom&lon(::4),tsf_anom&lat(::4), False)
plotts = gsn_csm_contour_map(wks,tsf_tmp(0,:,:),res_tsanom)

draw(plotts)
frame(wks)

delete(tsf_tmp)

delete(wks)
delete(cmap)

; ====== Plot FLUT ============

print("Plotting FLUT...")
wks = gsn_open_wks(out_type,"flut_"+plotname+"_"+yearstr+monstr+daystr+hourstr)

res_flut = res

cmap = read_colormap_file("MPL_gist_yarg")
res_flut@cnFillPalette=cmap

res_flut@cnLevelSelectionMode = "ManualLevels"
res_flut@cnLevelSpacingF      =  10.
res_flut@cnMinLevelValF       =  120.0
res_flut@cnMaxLevelValF       =  250.0

res_flut@mpFillOn = True
res_flut@mpOceanFillColor       = "black"
res_flut@mpLandFillColor        = "darkgreen"
res_flut@mpInlandWaterFillColor = "black"

res_flut@gsnRightString = "Outgoing longwave rad. (W/m2)"

res_flut@mpGeophysicalLineColor="green"

;flut = where(flut.ge.245,flut@_FillValue,flut)
plotflut = gsn_csm_contour_map(wks,flut(0,:,:),res_flut)
draw(plotflut)
frame(wks)

delete(wks)
delete(res_flut)
delete(cmap)

; ====== Plot 850 WIND ============

print("Plotting 850mb wind + MSLP...")

wks = gsn_open_wks(out_type,"wind_"+plotname+"_"+yearstr+monstr+daystr+hourstr)

res_wind = res

cmap = read_colormap_file("WhBlGrYeRe")
res_wind@cnFillPalette=cmap

res_wind@cnLevelSelectionMode = "ManualLevels"    ; set contour intervals
res_wind@cnLevelSpacingF      =  2.0                ;
res_wind@cnMinLevelValF       =  4.0           ;
res_wind@cnMaxLevelValF       =  46.0 ; Use 80.0 -> 60 for this

res_wind@gsnRightString = "MSLP/850mb wind"

plotwind = gsn_csm_contour_map(wks,wind(0,:,:),res_wind)
plotmslp = gsn_csm_contour(wks,mslp(0,:,:),res_mslp)
overlay(plotwind,plotmslp)

draw(plotwind)
frame(wks)

delete(wks)
delete(res_wind)
delete(cmap)

;system("rm wind_"+plotname+"_"+yearstr+monstr+daystr+hourstr+".000001.png")
;system("mv wind_"+plotname+"_"+yearstr+monstr+daystr+hourstr+".000002.png wind_"+plotname+"_"+yearstr+monstr+daystr+hourstr+".png")

; ========================= Plot 250-850 SHEAR ===============================

;print("Plotting 250-850 shear...")
;out_type = "png"
;wks = gsn_open_wks(out_type,"shear850250_"+plotname+"_"+yearstr+monstr+daystr+hourstr)
;
;
;res_shear850250 = res
;
;cmap = read_colormap_file("BlAqGrYeOrReVi200")
;res_shear850250@cnFillPalette=cmap
;;gsn_define_colormap(wks,"BlAqGrYeOrReVi200")
;
;res_shear850250@cnLevelSelectionMode = "ManualLevels"    ; set contour intervals
; res_shear850250@cnLevelSpacingF      =  5.0                ;
; res_shear850250@cnMinLevelValF       =  0.0           ;
; res_shear850250@cnMaxLevelValF       =  100.0 ; Use 80.0 -> 60 for this
;
;res_shear850250@gsnLeftString = "850-250mb shear mag."
;
;plotshear850250 = gsn_csm_contour_map(wks,shear850250(0,:,:),res_shear850250)
;draw(plotshear850250)
;frame(wks)
;delete(wks)
;delete(res_shear850250)
;delete(cmap)




; delete res at end
delete(res)

end do

end
