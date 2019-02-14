;**************************************
;  meteo_2.ncl
;
; Concepts illustrated:
;   - Drawing a meteogram
;   - Creating a color map using RGB triplets
;   - Reversing the Y axis
;   - Explicitly setting tickmarks and labels on the bottom X axis
;   - Increasing the thickness of contour lines
;   - Drawing wind barbs
;   - Drawing a bar chart
;   - Changing the width and height of a plot
;   - Overlaying wind barbs and line contours on filled contours
;   - gsnAttachPlotsXAxis = True
;
;***********************************************
; These files are loaded by default in NCL V6.2.0 and newer
; load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
; load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
;***********************************************
begin

all_files = systemfunc("ls ~/scratch/forecast_conus_30_x8_CAM5_L30/run/2018122600/*h0*nc")
fall = addfiles (all_files, "r")

TESTCIT=(/"STATECOLLEGE"/)
TESTSTR=(/"KUNV"/)
;TESTLON=(/-77.848611/)
;TESTLAT=(/ 40.849167/)
TESTLON=(/-96.7898/)
TESTLAT=(/ 46.8772/)

f1=addfile(all_files(0),"r")
lat=f1->lat
lon=f1->lon

do ii = 0,dimsizes(TESTLAT)-1
  gcdist = gc_latlon(TESTLAT(ii),TESTLON(ii),lat,lon,2,4)
  ix = minind(gcdist)
  print(ix+"")
  print(lat(ix)+"")
  print(lon(ix)+"")
end do

tempht = fall[:]->TREFHT(:,ix)
tempht = (tempht-273.15)*9.0/5.0+32.0
rain03 = fall[:]->PRECT(:,ix)
rain03 = rain03*141732*3
flut = fall[:]->FLUT(:,ix)
ptype = fall[:]->PTYPE(:,ix)
;z3tmp = fall[:]->Z3(:,:,ix)
;z3= z3tmp(lev|:,time|:)  
taus = fall[:]->time
levels = fall[0]->lev    
print(taus)
print(levels)

;***********************************************
;  open a workstation and define colormap using rgb triplets
;***********************************************
  wks = gsn_open_wks("x11","meteo") ; send graphics to PNG file
  colors = (/                                                        \
             (/255,255,255/),  \
             (/255,255,255/),   (/255,255,255/),   (/240,255,240/),  \
             (/220,255,220/),   (/190,255,190/),   (/120,255,120/),  \
             (/ 80,255, 80/),   (/ 50,200, 50/),   (/ 20,150, 20/) /) / 255.
;***********************************************
; Create a different resource list for each of the 5 plots
;***********************************************
  rh_res      = True
  temp_res    = True
  uv_res      = True
  rain_res    = True
  tempsfc_res = True
  flut_res    = True

;***********************************************
; Here are the resources that rh_res, temp_res, and uv_res share.
;***********************************************
  flut_res@gsnDraw      = False    ; Don't draw individual plot.
  flut_res@gsnFrame     = False    ; Don't advance frame.
  flut_res@vpXF         = 0.15     ; x location
  flut_res@vpYF         = 0.90     ; y location
  flut_res@vpWidthF     = 0.7      ; width
  flut_res@vpHeightF    = 0.40     ; height

  flut_res@trXMaxF         = max(taus)   ; max value on x-axis
  flut_res@trYMaxF         = max(flut)+10.   ; max value on y-axis
  flut_res@trYMinF         = min(flut)-10.   ; min value on y-axis
  
  flut_res@tiXAxisString   = ""            ; turn off x-axis string
  flut_res@tiYAxisString   = "FLUT"  ; set y-axis string
  flut_res@tmXBValues      = taus       ; location of explicit labels
  flut_res@tmXBLabelAngleF  = 90.            ; change label angle
  flut_res@tmXBMinorOn      = False          ; No minor tick marks.
  flut_res@tmXBLabelJust    = "CenterCenter" ; label justification
  flut_res@tmXBLabelFontHeightF = .016       ; Font size
 
  flut_res@tmXTOn           = False       ; turn off the top tickmarks
  flut_res@xyLineThicknesses= 2           ; increase line thickness
  flut_res@xyLineColor      =  "black"      ; set line color
                                     
;***********************************************
; rain histogram only resources
;***********************************************
  rain_res@vpXF            = 0.15   ; The left side of the box location
  rain_res@vpYF            = 0.34   ; The top side of the plot box location
  rain_res@vpWidthF        = 0.70   ; The Width of the plot box 
  rain_res@vpHeightF       = 0.10   ; The height of the plot box
  rain_res@trYMinF         = 0.0    ; min value on y-axis
  rain_res@trYMaxF         = 0.6    ; max value on y-axis
  rain_res@trXMinF         = 0.0    ; min value on x-axis
  rain_res@trXMaxF         = max(taus)   ; max value on x-axis

  rain_res@tmYLMode        = "Explicit"            ; Define own tick mark labels.
  rain_res@tmYLValues      = (/0.00, 0.10 , 0.20 , 0.30 , 0.40 , 0.50 , 0.60 /)  ; location of explicit labels
  rain_res@tmYLLabels      = (/ "" ,"0.10","0.20","0.30","0.40","0.50",  ""  /)
 ;rain_res@tmYLFormat      = "f"              ; no unnecessary 0s
  rain_res@tmYLFormat      = "0f"             ; equal number of significant digits

  rain_res@tiXAxisString  = ""                ; X axis label.
  rain_res@tiYAxisString  = "3hr rain total"  ; Y axis label.
  rain_res@tmXBMode       = "Explicit"        ; Define own tick mark labels.
  rain_res@tmXBValues     = taus              ; location of explicit labels
  rain_res@tmXTOn         = False             ; turn off the top tickmarks
  rain_res@tmXBMinorOn    = False             ; No minor tick marks.

  rain_res@gsnDraw         = False            ; Don't draw individual plot.
  rain_res@gsnFrame        = False            ; Don't advance frame.
  rain_res@gsnYRefLine     = 0.0              ; create a reference line
  rain_res@gsnXYBarChart   = True             ; turn on bar chart
  
  prect_bar_colors = new(dimsizes(taus),"string")
  prect_bar_cmap = (/"dodgerblue","","green","orange","red"/)
  print(ptype)
  do ii = 0,dimsizes(taus)-1
    prect_bar_colors(ii) = prect_bar_cmap(toint(ptype(ii)))
  end do
  rain_res@gsnAboveYRefLineBarColors = prect_bar_colors    ; above ref line fill green


;***********************************************
; temperature at surface only resources
;***********************************************
  tempsfc_res@vpXF            = 0.15   ; The left side of the box location
  tempsfc_res@vpYF            = 0.18   ; The top side of the plot box loc
  tempsfc_res@vpWidthF        = 0.70   ; The Width of the plot box
  tempsfc_res@vpHeightF       = 0.10   ; The height of the plot box
  tempsfc_res@trXMaxF         = max(taus)   ; max value on x-axis
  tempsfc_res@trYMaxF         = max(tempht)+2.   ; max value on y-axis
  tempsfc_res@trYMinF         = min(tempht)-2.   ; min value on y-axis

  tempsfc_res@tmYLValues      = ispan(59,64,1)*1.0   ; location of explicit labels
  tempsfc_res@tmYLLabels      = (/"59", "60", "61", "62", "63", "" /)
  tempsfc_res@tmYLFormat      = "f"          ; no unneccessary 0s
  tempsfc_res@tmYLFormat      = "0f"         ; equal number of significant digits
  
  tempsfc_res@tiXAxisString   = ""            ; turn off x-axis string
  tempsfc_res@tiYAxisString   = "Temp at 2m"  ; set y-axis string
  tempsfc_res@tmXBValues      = taus       ; location of explicit labels
  tempsfc_res@tmXBLabelAngleF  = 90.            ; change label angle
  tempsfc_res@tmXBMinorOn      = False          ; No minor tick marks.
  tempsfc_res@tmXBLabelJust    = "CenterCenter" ; label justification
  tempsfc_res@tmXBLabelFontHeightF = .016       ; Font size
 
  tempsfc_res@tmXTOn           = False       ; turn off the top tickmarks
  tempsfc_res@xyLineThicknesses= 2           ; increase line thickness
  tempsfc_res@xyLineColor      =  "red"      ; set line color

  tempsfc_res@gsnDraw          = False       ; Don't draw individual plot.
  tempsfc_res@gsnFrame         = False       ; Don't advance frame.

;********************************************************
; Create the five plots (they won't get drawn here, because gsnDraw
; was set to False for all of them.
;********************************************************
  rhfill    = gsn_csm_xy(wks,taus,flut,flut_res)
  rainhist  = gsn_csm_xy (wks,taus,rain03,rain_res)
  temptmsz  = gsn_csm_xy (wks,taus,tempht,tempsfc_res)

;********************************************************
; three of the above plots become one overlay plot
;********************************************************
  ;overlay(rhfill,templine)   ; result is templine over rhfill
  ;overlay(rhfill,windlayer)  ; result is windlayer over the other two

;********************************************************
; resources for attaching the plots                      
;********************************************************
; Set up resource lists for attaching the plot.  The res1 will apply to the base plot, 
; and the res2 to the plots being attached.  These resources lists are *not* for 
; changing things like line color, but for changing things like whether the plots
; are maximized, and which axis they are attached on.
;********************************************************
 
  res1                     = True
  res2                     = True
  res2@gsnAttachPlotsXAxis = True
  amid = gsn_attach_plots(rhfill,(/rainhist,temptmsz/),res1,res2)
  draw(rhfill)            ; 'base' plot   
  frame(wks)
end

