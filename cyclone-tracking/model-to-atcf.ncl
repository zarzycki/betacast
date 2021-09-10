;========================================================
; This script reads in a TempestExtremes trajectory file
; and produces an ATCF-formatted output given matching tcvitals
;========================================================
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"
load "./ncl-functions/getTrajectoriesGeneric.ncl"

;========================================================
begin

;thefile="trajectories.txt.NATLANTIC30X4-CYCLONES"
;vitalsfile="tcvitals"
;outfile="ATCF"

thefile=getenv("TRAJFILE")
vitalsfile=getenv("TCVITFILE")
outfile=getenv("ATCFFILE")
TECH=getenv("ATCFTECH")
print(thefile)
print(vitalsfile)
print(outfile)
print(TECH)

;--------------------------------------------------------

TIMETHRESH=12.  ; max number of hours init in tempest can be off from vitals
DISTTHRESH=5.0   ; max horiz distance (GC deg) init in tempest can be off from vitals

output_type="png"
cutwind=0.0

ncols = 9                     ; assign number of columns
maxLines = -1
headerDelimStr = "start"
delimStr="tab"
isHeader=False
traj_data = getTrajectories(thefile,ncols,maxLines,headerDelimStr,delimStr,isHeader)

lon = traj_data(:,1,:)
lat = traj_data(:,2,:)
wind = traj_data(:,4,:)
pres = traj_data(:,3,:)/100.
yyyy=  traj_data(:,5,:)
mm= traj_data(:,6,:)
dd= traj_data(:,7,:)
hh= traj_data(:,8,:)

ndims = dimsizes(lon)
nstorms = ndims(0)
ntimes = ndims(1)

tempesthit = new(nstorms,"string")
tempesthit = ""
;--------------------------------------------------------

data   = asciiread(vitalsfile,-1,"string")
v_id   = str_get_cols(data,5,7)
v_name  = str_get_cols(data,9,18)
v_yyyy = stringtoint(str_get_cols(data,19,22))
v_mm= stringtoint(str_get_cols(data,23,24))
v_dd= stringtoint(str_get_cols(data,25,26))
v_hh= stringtoint(str_get_cols(data,28,29))
v_lat=stringtofloat(str_get_cols(data,33,35))/10.
v_latdir= str_get_cols(data,36,36)
v_lon=stringtofloat(str_get_cols(data,38,41))/10.
v_londir= str_get_cols(data,42,42)

print(v_mm+" "+v_lon+" "+v_londir)

;NHC  12L JOSE      20170920 1200 388N 0702W 035 036 0976 1012 0556 28 093 0334 0334 0297 0297 D

num_vit = dimsizes(v_id)

do ii = 0,num_vit-1
  if (v_latdir(ii) .eq. "S") then
    v_lat(ii) = -v_lat(ii)
  end if
  if (v_londir(ii) .eq. "W") then
    v_lon(ii) = 360.-v_lon(ii)
  end if
end do

print(v_lat+" "+v_lon)

timeunits="hours since 1900-01-01 00:00:00"
opt=0
do ii = 0,num_vit-1
  print("    "+ii+" "+v_name(ii))
  v_single_time=cd_inv_calendar(v_yyyy(ii),v_mm(ii),v_dd(ii),v_hh(ii),0,0,timeunits,opt)
  do jj = 0,nstorms-1
    init_single_time=cd_inv_calendar(yyyy(jj,0),mm(jj,0),dd(jj,0),hh(jj,0),0,0,timeunits,opt)
    diff_single_time=init_single_time-v_single_time   ; hours b/w tempest init and tcvitals
    diff_single_dist=gc_latlon(v_lat(ii),v_lon(ii),lat(jj,0),lon(jj,0),2,2) ; dist b/w tempest init and tcvitals
    ;print(jj+"   difftime "+diff_single_time+" hrs --- distsep "+diff_single_dist+" deg")
    if (diff_single_time .le. TIMETHRESH .and. diff_single_dist .lt. DISTTHRESH) then
      print("Storm "+jj+"  in tempest file is "+v_name(ii))
      tempesthit(jj) = v_id(ii)
      do zz = 0,ntimes-1
      if (.not. ismissing(lat(jj,zz))) then
        if (zz .ne. 0)
          disttrav = gc_latlon(lat(jj,zz-1),lon(jj,zz-1),lat(jj,zz),lon(jj,zz),2,2)
        else
          disttrav = 0.
        end if
        print(v_id(ii)+" "+v_name(ii)+" "+lat(jj,zz)+" "+lon(jj,zz)+" "+yyyy(jj,zz)+" "+mm(jj,zz)+" "+dd(jj,zz)+" "+hh(jj,zz)+" dist: "+disttrav)
      end if
      end do
      break
    end if
  end do
end do

; WRITE ATCF FILE
foundstorms = ind(tempesthit .ne. "")
if (.not. all(ismissing(foundstorms)))
INITDATE=v_yyyy(0)*1000000+v_mm(0)*10000+v_dd(0)*100+v_hh(0)
inittime=cd_inv_calendar(v_yyyy(0),v_mm(0),v_dd(0),v_hh(0),0,0,timeunits,opt)
ms_to_kts = 1.94384
delim=", "
TECHNUM="03"
TY="XX"
FILLER="34, NEQ,    0,    0,    0,    0, "
quote = inttochar(34)
do jj = 0,dimsizes(foundstorms)-1
  do zz = 0,ntimes-1
      ix = foundstorms(jj)
      if (.not. ismissing(hh(ix,zz)) .and. mod(hh(ix,zz),6) .eq. 0)

        ; WORK HERE
        temptime=cd_inv_calendar(yyyy(ix,zz),mm(ix,zz),dd(ix,zz),hh(ix,zz),0,0,timeunits,opt)
        LEADTIME = toint(temptime - inittime)

        tmp = str_split_by_length(tempesthit(ix),(/2,1/))
        STORMID = toint(tmp(0))
        tmpbasin = tmp(1)
        if (tmpbasin .eq. "L") then
          STORMBASIN="AL"
        else if (tmpbasin .eq. "E") then
          STORMBASIN="EP"
        else if (tmpbasin .eq. "W") then
          STORMBASIN="WP"
        else if (tmpbasin .eq. "C") then
          STORMBASIN="CP"
        else if (tmpbasin .eq. "B" .or. tmpbasin .eq. "A") then
          STORMBASIN="IO"
        else if (tmpbasin .eq. "Q" .or. tmpbasin .eq. "P" .or. tmpbasin .eq. "S") then
          STORMBASIN="SH"
        else
          STORMBASIN="XX"
        end if       
        end if
        end if
        end if       
        end if
        end if

        templat=lat(ix,zz)
        if (templat .ge. 0.0) then
          LATDIR="N"
        else
          LATDIR="S"
        end if      
        LAT = toint(abs(templat)*10.)

        templon=lon(ix,zz)
        if (templon .ge. 180.0) then
          LONDIR="W"
        else
          LONDIR="E"
        end if
        if (templon .ge. 180.) then
          templon = templon-360.
        end if
        LON = toint(abs(templon)*10.)

        WIND = round((wind(ix,zz)*ms_to_kts),3)
        PRES = round(pres(ix,zz),3)

        ;AL, 01, 2012052106, 03, AP06,  60, 375N,  706W,  18, 1013, XX,  34, NEQ,    0,    0,    0,    0,
        ATCFSTRING=STORMBASIN+delim+sprinti("%0.2i",STORMID)+delim+sprinti("%0.10i",INITDATE)+delim+ \ 
          TECHNUM+delim+TECH+delim+sprinti("%0.3i",LEADTIME)+delim + \
          sprinti("%0.3i",LAT)+LATDIR+delim+sprinti("%0.4i",LON)+LONDIR+delim + \
          sprinti("%0.3i",WIND)+delim+sprinti("%0.4i",PRES)+delim+TY+delim+FILLER

        system("echo "+quote+ATCFSTRING+quote+" >> "+outfile)
      end if
  end do
end do
else
  print("Found no matches to tempest + vitals, cannot write ATCF")
end if

;--------------------------------------------------------
print("Beginning plot sequence")
wks  = gsn_open_wks(output_type,"trajs."+thefile)

numpanels=4
map = new(numpanels,graphic)
dum = new((/numpanels,nstorms,ntimes/),graphic)
text = new((/numpanels,nstorms/),graphic)

res  = True
res@gsnDraw    = False
res@gsnFrame   = False
res@gsnMaximize = True
res@mpFillOn             = True              ; turn off gray continents
res@mpLandFillColor = "Tan"
res@mpOceanFillColor = "LightBlue1"
res@mpOutlineOn = True
res@mpInlandWaterFillColor = res@mpOceanFillColor
res@mpGreatCircleLinesOn = True
 
do zz = 0,3

if (zz .eq. 0)  ;atlantic
  res@mpMinLatF             = 5
  res@mpMaxLatF             = 45
  res@mpMinLonF             = 260
  res@mpMaxLonF             = 355
else if (zz .eq. 1) ; epac
  res@mpMinLatF  = 5
  res@mpMaxLatF = 35
  res@mpMinLonF  = -165
  res@mpMaxLonF = -90
else if (zz .eq. 2) ; epac
  res@mpMinLatF  = 5
  res@mpMaxLatF = 45
  res@mpMinLonF  = 100
  res@mpMaxLonF = 195
else if (zz .eq. 3)
  res@mpMinLatF             = 0
  res@mpMaxLatF             = 75
  res@mpMinLonF             = 50
  res@mpMaxLonF             = 350
else
  print("invalid basin")
end if
end if
end if
end if

map(zz) = gsn_csm_map_ce(wks,res)

gsn_define_colormap(wks,"rainbow+white+gray")
  polyres               = True          ; poly marker mods desired 
  polyres@gsMarkerIndex = 1            ; choose circle as polymarker
  lineres = True

do i = 0,nstorms-1
if (tempesthit(i) .ne. "") then
do j = 0,ntimes-1

thisWind=wind(i,j)

;lineres@gsLineThicknessF = 8.0
; We need to use thinner lines for the global map
  linethicknesses = (/0.5,1.0,1.3,1.9,2.3,2.6/)
  linethicknesses = 4.0*linethicknesses

if (.not.ismissing(thisWind))

;  polyres@gsMarkerSizeF = 2.0           ; select size to avoid streaking
  if (thisWind .lt. 17)
    polyres@gsMarkerColor = (/"blue"/)
    lineres@gsLineColor = (/"blue"/)
    lineres@gsLineThicknessF = linethicknesses(0)
  else if (thisWind .ge. 17 .and. thisWind .lt. 33)
    polyres@gsMarkerColor = (/"green3"/)
    lineres@gsLineColor = (/"green3"/)
    lineres@gsLineThicknessF = linethicknesses(1)
  else if (thisWind .ge. 33 .and. thisWind .lt. 43)
    polyres@gsMarkerColor = (/"yellow2"/)
    lineres@gsLineColor = (/"yellow2"/)
    lineres@gsLineThicknessF = linethicknesses(2)
  else if (thisWind .ge. 43 .and. thisWind .lt. 49)
    polyres@gsMarkerColor = (/"orange"/)
    lineres@gsLineColor = (/"orange"/)
    lineres@gsLineThicknessF = linethicknesses(3)
  else if (thisWind .ge. 49 .and. thisWind .lt. 58)
    polyres@gsMarkerColor = (/"darkorange3"/)
    lineres@gsLineColor = (/"darkorange3"/)
    lineres@gsLineThicknessF  = linethicknesses(4)
  else
    polyres@gsMarkerColor = (/"red"/) 
    lineres@gsLineColor = (/"red"/)
    lineres@gsLineThicknessF  = linethicknesses(5)
  end if
  end if
  end if
  end if
  end if
  
 txres               = True                            
 txres@txFontHeightF = 0.01             ; Set the font height
  
  if (j .ne. 0) then
    if(thisWind .ge. cutwind) then
      dum(zz,i,j) = gsn_add_polyline(wks,map(zz),(/lon(i,j),lon(i,j-1)/),(/lat(i,j),lat(i,j-1)/),lineres)
    end if
  end if
  
end if

end do

txres = True         
txres@amJust = "TopLeft"                   
label = tempesthit(i)
text(zz,i) = gsn_add_text(wks,map(zz),label,lon(i,0)+0.5,lat(i,0)-0.5,txres)

end if
end do

; if we are doing a varres grid, let's push an outline of the atl30x4 grid
;if (modelgrid .eq. "vr") then
;  map = plot_outline_fcn(wks,map,"black",3.0)
;end if

end do

;******************
; Plot a legend   *
;******************

lgres                    = True

lgres@lgAutoManage       = False

lgres@vpWidthF           = 0.08       ; was 0.08        ; width of legend (NDC)
lgres@vpHeightF          = 0.12       ; was 0.08        ; height of legend (NDC)
; lgres@lgBottomMarginF    = 0.17     ; was 0.25

lgres@lgPerimFill        = 0                     ; Use solid fill (0) instead of the default hollow fill
lgres@lgPerimFillColor   = "Background"

;lgres@lgOrientation = "Horizontal"

; lgres@lgBoxMajorExtentF  = 0.4
lgres@lgBoxMinorExtentF  = 0.2	; controls how wide the box holding the legend items (lines and markers) can be in relation to legend

; lgres@lgBoxBackground    = "PaleTurquoise3"

lgres@lgMonoItemType        = False                 ; indicates that we wish to set the item types individually
lgres@lgMonoMarkerIndex     = False
lgres@lgMonoLineThickness   = False
lgres@lgMonoMarkerThickness = False
lgres@lgMonoMarkerSize      = False

;lgres@lgLabelFont             = 0
lgres@lgLabelFontHeightF      = 0.07
lgres@lgLabelFontAspectF      = 1.5
lgres@lgLabelConstantSpacingF = 0.0

legendType = "Markers"
markerSize = .008
markerIndex = 16

if (cutwind .eq. 0) then
lgres@lgItemCount        = 6
;lgres@lgLineDashSegLenF  = 0.8
lgres@lgItemTypes        = (/"Markers","Markers","Markers","Markers","Markers","Markers"/)
lgres@lgMarkerIndexes    = (/      16,      16,      16,      16,      16,      16/)
lgres@lgMarkerThicknesses  = (/    markerSize,    markerSize,    markerSize,    markerSize,    markerSize,    markerSize/)
lgres@lgMarkerColors       = (/ "red",     "darkorange3",   "orange",  "yellow",   "green",    "blue"    /) ; colors for legend lines
;lgres@lgDashIndexes      = (/      0,          0,          0,         0,               0,      0      /) ; dash indexes
legend_labels = (/"Cat. 4/5","Cat. 3","Cat. 2",\
               "Cat. 1","TS","TD"/)
else if (cutwind .eq. 17) then
  lgres@lgItemCount        = 5
  lgres@lgItemTypes        = (/legendType,legendType,legendType,legendType,legendType/)
  lgres@lgMarkerIndexes    = (/markerIndex,markerIndex,markerIndex,markerIndex,markerIndex/)
  lgres@lgMarkerThicknesses  = (/    markerSize,    markerSize,    markerSize,    markerSize,    markerSize/)
  lgres@lgMarkerColors       = (/ "red",     "darkorange3",   "orange",  "yellow2",   "green3" /) ; colors for legend lines
  legend_labels = (/"Category 4/5","Category 3","Category 2",\
                 "Category 1","Tropical Storm"/)
else
end if
end if


legend = gsn_create_legend(wks,lgres@lgItemCount,legend_labels,lgres)

amres = True

amres@amJust="TopLeft"
amres@amParallelPosF   = -.49		    ; move legend to the right
amres@amOrthogonalPosF = -.49             ; move the legend down


annoid1 = gsn_add_annotation(map(3),legend,amres)   ; add legend to plot
  
print("Drawing plot")

resP = True
resP@amJust           = "TopRight"
resP@gsnPanelFigureStringsFontHeightF = 0.027
gsn_panel(wks,map,(/4,1/),resP)             ; now draw as one plot

end
