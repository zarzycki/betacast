begin

  ; Check for NLEV as commanline argument
  ;----------------------------------------
  if(isdefined("NLEV")) then
   print("Number of Vertical Levels: nlev:")
   nlev = NLEV
   print("nlev="+nlev)
   print("=====")
  else
   print("Using Default number of Vertical Levels: nlev = 58")
   nlev = 58
  end if
  printVarSummary(nlev)

  use_nclev=True
  lev_template = "../templates/L"+nlev+"template.nc"
  lt = addfile(lev_template,"r")
  nclev = lt->lev

  ; Read in Window paramters from namelist
  ;----------------------------------------
  Nudge_Hwin_lo=0.
  Nudge_Hwin_hi=0.
  Nudge_Vwin_lo=0.
  Nudge_Vwin_hi=0.
  Nudge_Hwin_lat0=0.
  Nudge_Hwin_lon0=0.
  Nudge_Hwin_latWidth=0.
  Nudge_Hwin_lonWidth=0.
  Nudge_Hwin_latDelta=0.
  Nudge_Hwin_lonDelta=0.
  Nudge_Vwin_Lindex=0.
  Nudge_Vwin_Hindex=0.
  Nudge_Vwin_Ldelta=0.
  Nudge_Vwin_Hdelta=0.
  Nudge_Hwin_Invert = "False"
  Nudge_Vwin_Invert = "False"

  pthi="user_nl_cam"

   Nudge_Hwin_Invert = keyword_values( pthi, "Nudge_Hwin_Invert", "string")
   Nudge_Vwin_Invert = keyword_values( pthi, "Nudge_Vwin_Invert", "string")

   Nudge_Hwin_lat0    = keyword_values( pthi, "Nudge_Hwin_lat0", "float")
   Nudge_Hwin_lon0    = keyword_values( pthi, "Nudge_Hwin_lon0", "float")

   Nudge_Vwin_lo      = keyword_values( pthi, "Nudge_Vwin_lo", "float")
   Nudge_Vwin_hi      = keyword_values( pthi, "Nudge_Vwin_hi", "float")

   Nudge_Hwin_lo      = keyword_values( pthi, "Nudge_Hwin_lo", "float")
   Nudge_Hwin_hi      = keyword_values( pthi, "Nudge_Hwin_hi", "float")

   Nudge_Hwin_latWidth  = keyword_values( pthi, "Nudge_Hwin_latWidth", "float")
   Nudge_Hwin_lonWidth  = keyword_values( pthi, "Nudge_Hwin_lonWidth", "float")

   Nudge_Hwin_latDelta  = keyword_values( pthi, "Nudge_Hwin_latDelta", "float")
   Nudge_Hwin_lonDelta  = keyword_values( pthi, "Nudge_Hwin_lonDelta", "float")

   Nudge_Vwin_Lindex  = keyword_values( pthi, "Nudge_Vwin_Lindex", "float")
   Nudge_Vwin_Hindex  = keyword_values( pthi, "Nudge_Vwin_Hindex", "float")

   Nudge_Vwin_Ldelta  = keyword_values( pthi, "Nudge_Vwin_Ldelta", "float")
   Nudge_Vwin_Hdelta  = keyword_values( pthi, "Nudge_Vwin_Hdelta", "float")

; Old FORTRAN wrapper from patc
;  Read_Namelist::read_namelist(Nudge_Hwin_Invert,                      \
;                               Nudge_Vwin_Invert,                      \
;                               Nudge_Hwin_lat0,Nudge_Hwin_lon0,        \
;                               Nudge_Hwin_latWidth,Nudge_Hwin_lonWidth,\
;                               Nudge_Hwin_latDelta,Nudge_Hwin_lonDelta,\
;                               Nudge_Vwin_Lindex,Nudge_Vwin_Hindex,    \
;                               Nudge_Vwin_Ldelta,Nudge_Vwin_Hdelta )


; Here, we check to see if Nudge_Hwin_Invert (a string) is returned as "missing"
; if it is, it means it wasn't in the namelist. If missing, ignore the var
; If any other string, it means NCL read something in user_nl. Convert that to a bool
; and then invert appropriately.
; NOTE: Nudge_Hwin_Invert/Nudge_Vwin_Invert is only used in CAM as of 9/2023
; so this block should never be evaluated for E3SM

if (Nudge_Hwin_Invert .ne. "missing")

  val=str_strip(Nudge_Hwin_Invert)
  if (val.eq.".true." .or. val.eq.".TRUE." .or.  val.eq."TRUE" .or. \
      val.eq.".T." .or. val.eq.".t."       .or.  val.eq."T"    .or. \
      val.eq."True".or. val.eq.".True.") then
    Nudge_Hwin_Invert := True
  else
    Nudge_Hwin_Invert := False
  end if
  delete(val)

  if(Nudge_Hwin_Invert) then
    Nudge_Hwin_lo= 1.0
    Nudge_Hwin_hi= 0.0
  else
    Nudge_Hwin_lo= 0.0
    Nudge_Hwin_hi= 1.0
  end if
end if

if (Nudge_Vwin_Invert .ne. "missing")

  val=str_strip(Nudge_Vwin_Invert)
  if (val.eq.".true." .or. val.eq.".TRUE." .or.  val.eq."TRUE" .or. \
      val.eq.".T." .or. val.eq.".t."       .or.  val.eq."T"    .or. \
      val.eq."True".or. val.eq.".True.") then
    Nudge_Vwin_Invert := True
  else
    Nudge_Vwin_Invert := False
  end if
  delete(val)

  if(Nudge_Vwin_Invert) then
    Nudge_Vwin_lo= 1.0
    Nudge_Vwin_hi= 0.0
  else
    Nudge_Vwin_lo= 0.0
    Nudge_Vwin_hi= 1.0
  end if
end if

  print("Nudge_Hwin_lo="+Nudge_Hwin_lo)
  print("Nudge_Hwin_hi="+Nudge_Hwin_hi)
  print("Nudge_Vwin_lo="+Nudge_Vwin_lo)
  print("Nudge_Vwin_hi="+Nudge_Vwin_hi)
  print("Nudge_Hwin_lat0="+Nudge_Hwin_lat0)
  print("Nudge_Hwin_lon0="+Nudge_Hwin_lon0)
  print("Nudge_Hwin_latWidth="+Nudge_Hwin_latWidth)
  print("Nudge_Hwin_lonWidth="+Nudge_Hwin_lonWidth)
  print("Nudge_Hwin_latDelta="+Nudge_Hwin_latDelta)
  print("Nudge_Hwin_lonDelta="+Nudge_Hwin_lonDelta)
  print("Nudge_Vwin_Lindex="+Nudge_Vwin_Lindex)
  print("Nudge_Vwin_Hindex="+Nudge_Vwin_Hindex)
  print("Nudge_Vwin_Ldelta="+Nudge_Vwin_Ldelta)
  print("Nudge_Vwin_Hdelta="+Nudge_Vwin_Hdelta)

  Hlo        = Nudge_Hwin_lo
  Hhi        = Nudge_Hwin_hi
  Vlo        = Nudge_Vwin_lo
  Vhi        = Nudge_Vwin_hi
  lat0       = Nudge_Hwin_lat0
  lat_width  = Nudge_Hwin_latWidth
  lat_delta  = Nudge_Hwin_latDelta
  lon0       = Nudge_Hwin_lon0
  lon_width  = Nudge_Hwin_lonWidth
  lon_delta  = Nudge_Hwin_lonDelta
  levH       = Nudge_Vwin_Hindex
  levH_delta = Nudge_Vwin_Hdelta
  levL       = Nudge_Vwin_Lindex
  levL_delta = Nudge_Vwin_Ldelta

;----------------------
; HORIZONTAL WINDOW
;----------------------
  ; Create a Horizontal test array
  ;-------------------------------
  nlon=360
  nlat=180
  lon       = fspan(0.,(nlon*360./(nlon+1)),nlon)
  lon@units = "degrees_east"
  lon!0     = "lon"
  lon&lon   = lon
  lat       = fspan(-90.,90.,nlat)
  lat@units = "degrees_north"
  lat!0     = "lat"
  lat&lat   = lat
  Hcoef     = new((/ nlat, nlon /), "float")
  Hcoef!0   = "lat"
  Hcoef!1   = "lon"
  Hcoef&lat = lat
  Hcoef&lon = lon

  ; Set lat/lon profiles for window function
  ;-----------------------------------------------------
  lonx     = (lon - lon0)
  lonx     = where((lonx.le.-180.),(lonx+360.),lonx)
  lonx     = where((lonx.gt. 180.),(lonx-360.),lonx)
  lon0_min = -(lon_width/2.)
  lon0_max =  (lon_width/2.)
  lon_lo   = (lonx - lon0_min)/lon_delta
  lon_hi   = (lon0_max - lonx)/lon_delta

  lat0_min = lat0 - (lat_width/2.)
  lat0_max = lat0 + (lat_width/2.)
  lat_lo   = (lat - lat0_min)/lat_delta
  lat_hi   = (lat0_max - lat)/lat_delta

  ; Calculate min/max of RAW window function
  ;------------------------------------------
  Val1_p=((1.+tanh(( 180. - lon0_min)/lon_delta))/2.)
  Val1_0=((1.+tanh((   0. - lon0_min)/lon_delta))/2.)
  Val1_n=((1.+tanh((-179. - lon0_min)/lon_delta))/2.)
  Val2_p=((1.+tanh((lon0_max - 180. )/lon_delta))/2.)
  Val2_0=((1.+tanh((lon0_max -   0. )/lon_delta))/2.)
  Val2_n=((1.+tanh((lon0_max + 179. )/lon_delta))/2.)

  Val3_p=((1.+tanh(( 90. - lat0_min)/lat_delta))/2.)
  Val3_0=((1.+tanh((lat0 - lat0_min)/lat_delta))/2.)
  Val3_n=((1.+tanh((-90. - lat0_min)/lat_delta))/2.)
  Val4_p=((1.+tanh((lat0_max -  90.)/lat_delta))/2.)
  Val4_0=((1.+tanh((lat0_max - lat0)/lat_delta))/2.)
  Val4_n=((1.+tanh((lat0_max +  90.)/lat_delta))/2.)

  Hmax = Val1_0*Val2_0*Val3_0*Val4_0
  Htest = (/ Val1_p*Val2_p*Val3_n*Val4_n, \
             Val1_p*Val2_p*Val3_p*Val4_p, \
             Val1_n*Val2_n*Val3_n*Val4_n, \
             Val1_n*Val2_n*Val3_p*Val4_p /)
  Hmin = min(Htest)

  ; Loop over lat/lon and load the RAW window function
  ;---------------------------------------------------
  do ilat=0,(nlat-1)
    do ilon=0,(nlon-1)
      Hcoef(ilat,ilon)=((1.+tanh(lon_lo(ilon)))/2.)  \
                      *((1.+tanh(lon_hi(ilon)))/2.)  \
                      *((1.+tanh(lat_lo(ilat)))/2.)  \
                      *((1.+tanh(lat_hi(ilat)))/2.)
    end do ; ilon=0,(nlon-1)
  end do ; ilat=0,(nlat-1)

  ; Scale the Window function to span the values between Hlo and Hhi:
  ;-----------------------------------------------------------------
  if(Hmax.le.Hmin) then
    Hcoef=1.0
  else
    Hcoef = (Hcoef-Hmin)/(Hmax-Hmin)
  end if
  Hcoef = Hlo + Hcoef*(Hhi-Hlo)

  print(" ")
  print("lon0="+lon0+" lon_width="+lon_width+" lon_delta="+lon_delta)
  print("lat0="+lat0+" lat_width="+lat_width+" lat_delta="+lat_delta)
  print("Hlo="+Hlo+" Hhi="+Hhi)
  print(" ")
  print(" RAW(Hmin)="+Hmin+" RAW(Hmax)="+Hmax)
  print(" min(Hcoef)="+min(Hcoef)+" max(Hcoef)="+max(Hcoef))
  print(" ")

;----------------------
; VERTICAL WINDOW
;----------------------
  ; Create a Vertical test array
  ;-------------------------------
  print(" nlev="+nlev)
;  nlev=32
;  lev       = fspan(0.,(nlev-1),nlev)
  lev       = fspan(1.,nlev,nlev)
  lev@units = "model_level"
  lev!0     = "lev"
  lev&lev   = lev
  Vcoef     = new((/ nlev /), "float")
  Vcoef!0   = "lev"
  Vcoef&lev = lev

  ; Set level profiles for window function
  ;-----------------------------------------------------
  lev0   = (levL + levH)/2.
  ilev   = minind(abs(lev(:)-lev0))
  lev_lo = (lev - levL)/levL_delta
  lev_hi = (levH - lev)/levH_delta

  ; Loop over lat/lon and load the RAW window function
  ;---------------------------------------------------
  do ilev=0,(nlev-1)
    Vcoef(ilev)=((1.+tanh(lev_lo(ilev)))/2.)  \
               *((1.+tanh(lev_hi(ilev)))/2.)
  end do ; ilev=0,(nlev-1)

  ; Scale the Window function to span the values between Vlo and Vhi:
  ;-----------------------------------------------------------------
  Vmax = max(Vcoef)
  Vmin = min(Vcoef)
  if(Vmax.le.Vmin) then
    Vcoef=1.0
  else
    Vcoef = (Vcoef-Vmin)/(Vmax-Vmin)
  end if
  Vcoef = Vlo + Vcoef*(Vhi-Vlo)

  print(" ")
  print("levH="+levH+" levH_delta="+levH_delta)
  print("levL="+levL+" levL_delta="+levL_delta)
  print("Vlo="+Vlo+" Vhi="+Vhi)
  print(" ")
  print(" RAW(Vmin)="+Vmin+" RAW(Vmax)="+Vmax)
  print(" min(Vcoef)="+min(Vcoef)+" max(Vcoef)="+max(Vcoef))
  print(" ")

;-----------------
; GRAPH results
;-----------------

; cmz, removed the paneling by default because it looks ugly, followed...
; https://www.ncl.ucar.edu/Applications/Scripts/panel_20.ncl

  do_panel=False

  wks = gsn_open_wks("x11","Wcoef")
;  wks = gsn_open_wks("X11","Wcoef")

  gsn_define_colormap(wks,"BlAqGrYeOrRe")

  if (do_panel) then
    Tplot = new(4,graphic)
  end if

  ilat = minind(abs(lat(:)-lat0))
  ilon = minind(abs(lon(:)-lon0))
  Hprof_lon = Hcoef(ilat,:)
  Hprof_lat = Hcoef(:,ilon)
  Xres = True
  Xres@gsnDraw      = False
  Xres@gsnFrame     = False

  ; common font heights/tick settings
  Xres@tiMainFontHeightF = 0.012
  Xres@tmEqualizeXYSizes = True
  Xres@tmXBLabelFontHeightF = 0.01
  Xres@tmXBMajorLengthF  = 0.008
  Xres@tmXBMinorLengthF = 0.005

  Xres@trYMinF      = 0.
  Xres@trYMaxF      = 1.1

  Xres@tiMainString = "Lat~B~0~N~ Zonal Window Profile"
  Xres@trXMinF      = min(lon)
  Xres@trXMaxF      = max(lon)
  if (.not. do_panel) then
    Xres@vpXF             = 0.15
    Xres@vpYF             = 0.96
    Xres@vpWidthF         = 0.4
    Xres@vpHeightF        = 0.2
    Tplot0 = gsn_csm_xy(wks,Hprof_lon&lon,Hprof_lon,Xres)
  else
    Tplot(0) = gsn_csm_xy(wks,Hprof_lon&lon,Hprof_lon,Xres)
  end if

  Xres@tiMainString = "Lon~B~0~N~ Meridional Window Profile"
  Xres@trXMinF      = min(lat)
  Xres@trXMaxF      = max(lat)
  if (.not. do_panel) then
    Xres@vpXF             = 0.65
    Xres@vpYF             = 0.96
    Xres@vpWidthF         = 0.2
    Xres@vpHeightF        = 0.2
    Tplot1 = gsn_csm_xy(wks,Hprof_lat&lat,Hprof_lat,Xres)
  else
    Tplot(1) = gsn_csm_xy(wks,Hprof_lat&lat,Hprof_lat,Xres)
  end if



  Cres = True
  Cres@gsnSpreadColors = True
  Cres@gsnDraw         = False
  Cres@gsnFrame        = False
  Cres@cnFillOn        = True
  Cres@cnLevelSelectionMode  = "ManualLevels"
  Cres@cnMinLevelValF  = 0.01
  Cres@cnMaxLevelValF  = 1.
  Cres@cnLevelSpacingF = 0.100
  Cres@mpCenterLonF    = 180.
  Cres@mpGridAndLimbOn=True
  Cres@mpGridLatSpacingF=10.0
  Cres@mpGridLonSpacingF=10.0

  if (.not. do_panel) then
    Cres@vpXF             = 0.25
    Cres@vpYF             = 0.71
    Cres@vpWidthF         = 0.5
    Cres@vpHeightF        = 0.25
    Tplot2 = gsn_csm_contour_map_ce(wks,Hcoef,Cres)
  else
    Tplot(2) = gsn_csm_contour_map_ce(wks,Hcoef,Cres)
  end if



  Vres = True
  Vres@gsnDraw      = False
  Vres@gsnFrame     = False
  Vres@tiMainFontHeightF = 0.012
  Vres@tiMainString = "Vertical Window Profile"
  Vres@tiXAxisString= "Nudging strength"
  Vres@trXMinF      = -0.05
  Vres@trXMaxF      = 1.05
  if (use_nclev) then
    Vres@trYMinF = nclev(0)
    Vres@trYMaxF = nclev(nlev-1)
    Vres@tmYLMode = "Manual"
    Vres@tmYLTickSpacingF  = 100.0
    Vres@tmYLTickStartF  = 0.0
    Vres@tiYAxisString= "Approx model plev"
  else
    Vres@trYMinF = lev(0)
    Vres@trYMaxF = lev(nlev-1)
    Vres@tiYAxisString= "Model Level Index"
  end if
  Vres@trYReverse   = True

  ; Marker settings
  Vres@xyMarkLineMode = "MarkLines"
  Vres@xyMarker      = 16
  Vres@xyMarkerColor = "red"

  ; Underlying grid
  Vres@tmXMajorGrid                = True          ; implement x grid
  Vres@tmXMajorGridThicknessF      = 1.0           ; 2.0 is default
  Vres@tmXMajorGridLineDashPattern = 2             ; select short dash lines
  Vres@tmYMajorGrid                = True          ; implement y grid
  Vres@tmYMajorGridThicknessF      = 1.0           ; 2.0 is default
  Vres@tmYMajorGridLineDashPattern = 2             ; select short dash lines
  Vres@tmXMinorGrid                = True          ; implement x grid
  Vres@tmXMinorGridThicknessF      = 1.0           ; 2.0 is default
  Vres@tmXMinorGridLineDashPattern = 2             ; select short dash lines
  Vres@tmYMinorGrid                = True          ; implement y grid
  Vres@tmYMinorGridThicknessF      = 1.0           ; 2.0 is default
  Vres@tmYMinorGridLineDashPattern = 2             ; select short dash lines

  if (.not. do_panel) then
    Vres@vpXF             = 0.36
    Vres@vpYF             = 0.34
    Vres@vpWidthF         = 0.28
    Vres@vpHeightF        = 0.28
    if (use_nclev) then
      Tplot3 = gsn_csm_xy(wks,Vcoef,nclev,Vres)
    else
      Tplot3 = gsn_csm_xy(wks,Vcoef,Vcoef&lev,Vres)
    end if
  else
    if (use_nclev) then
      Tplot(3) = gsn_csm_xy(wks,Vcoef,nclev,Vres)
    else
      Tplot(3) = gsn_csm_xy(wks,Vcoef,Vcoef&lev,Vres)
    end if
  end if



;drawNDCGrid(wks)
if (.not. do_panel) then
 psres = True
 maximize_output(wks,psres)  ; calls draw and frame for you
else
  Pres                 = True
  Pres@gsnPanelRowSpec = True
  Pres@gsnPanelCenter  = True
  Pres@gsnPanelDebug  = True
  Pres@gsnMaximize  = True
  gsn_panel(wks,Tplot,(/2,1,1/),Pres)
end if

end
