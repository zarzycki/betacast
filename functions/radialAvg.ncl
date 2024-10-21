undef ( "radialAvg3D" )
function radialAvg3D ( \
  data[*][*][*]    : numeric, \
  lat[*]           : numeric, \
  lon[*]           : numeric, \
  lev[*]           : numeric, \
  psminlat         : numeric, \
  psminlon         : numeric, \
  outerRad         : numeric, \
  mergeInnerBins   : logical \
)
local None

; lat (deg)
; lon (deg)
; center_lat(deg)
; center_lon(deg)
; outerRad (km)

begin

nlat = dimsizes(lat)
nlon = dimsizes(lon)
nlev = dimsizes(lev)

deltaLat = abs(lat(0:nlat-2)-lat(1:nlat-1))
deltaLon = abs(lon(0:nlon-2)-lon(1:nlon-1))
deltaLatMax = max(deltaLat)
deltaLatMin = min(deltaLat)
deltaLonMax = max(deltaLon)
deltaLonMin = min(deltaLon)
kmInDeg = 111.32
kmGrid = kmInDeg * deltaLatMax
print("The max lat/lon km grid spacing at equator is "+kmGrid+" km")
print("The min lat/lon km grid spacing at equator is "+kmInDeg * deltaLatMin+" km")

pi = 3.14159
d2r = pi/180.
r2d = 180./pi
lonr = lon*d2r
latr = lat*d2r
copy_VarCoords(lon,lonr)
copy_VarCoords(lat,latr)

psminLocY = ind_nearest_coord(psminlat,lat,0)
psminLocX = ind_nearest_coord(psminlon,lon,0)

;;;;; DO RADIAL AVERAGES

loopWidthW = toint(ceil((outerRad*1.10) / (kmInDeg * deltaLonMin)))
loopWidthE = loopWidthW
loopWidthN = toint(ceil((outerRad*1.10) / (kmInDeg * deltaLatMin)))
loopWidthS = loopWidthN

; loopWidth checks
if (psminLocY-loopWidthS .lt. 0) then
  print("Correcting loopWidthS")
  loopWidthS = psminLocY
end if
if (psminLocY+loopWidthN .ge. nlat) then
  print("Correcting loopWidthN")
  loopWidthN = nlat - 1 - psminLocY
end if
if (psminLocX-loopWidthW .lt. 0) then
  print("Correcting loopWidthW")
  loopWidthW = 0
end if
if (psminLocX+loopWidthE .ge. nlon) then
  print("Correcting loopWidthE")
  loopWidthE = nlon - 1 - psminLocX
end if

print("nlat "+nlat+"  nlon "+nlon+"  nlev "+nlev)
; for standard lat/lon grid
;1.25 gives n ~= 5 in inner grid if numMerge = 2
;1.4 gives n ~= 9 in inner grid if numMerge = 2
timesGrid = 1.1 ; we want each radius bin to be timeGrid times kmGrid

nx = toint(outerRad/(timesGrid*kmGrid))
print("Number of bins is equal to "+nx)

if (mergeInnerBins) then
;  origRadiusArr = fspan(0,outerRad,nx)
;  radiusArr = new(dimsizes(origRadiusArr)-1,double)
;  radiusArr(0) = origRadiusArr(0)
;  radiusArr(1:dimsizes(radiusArr)-1) = origRadiusArr(2:dimsizes(origRadiusArr)-1)
  numMerge = 2 ; numMerge is the number of innermost radial bins we want to merge
  print("Merging innermost "+numMerge+" bins because radial average is so small.")
  numMergeMinusOne = numMerge-1
  origRadiusArr = fspan(0,outerRad,nx+numMergeMinusOne)
  radiusArr = new((dimsizes(origRadiusArr)-numMergeMinusOne),double)
  radiusArr(0) = origRadiusArr(0)
  printVarSummary(radiusArr)
  printVarSummary(origRadiusArr)
  radiusArr(1:dimsizes(radiusArr)-1) = origRadiusArr(numMerge:dimsizes(origRadiusArr)-1)
else
  print("Not merging any innermost bins -- be careful that your inner bins have > 1 pt.")
  radiusArr = fspan(0,outerRad,nx)
end if

numRadBins = dimsizes(radiusArr)

rad_thevar_hit  = new((/nlev,numRadBins/),integer)
rad_thevar_cum  = new((/nlev,numRadBins/),double)

; NEEDS to be 1 or 0
rad_thevar_hit!1 = "radius"
rad_thevar_cum!1 = "radius"
rad_thevar_hit&radius = radiusArr
rad_thevar_cum&radius = radiusArr

rad_thevar_hit!0 = "lev"
rad_thevar_cum!0 = "lev"
rad_thevar_hit&lev = lev
rad_thevar_cum&lev = lev

rad_thevar_hit = 0
rad_thevar_cum = 0

; NEED TO CORRECT NLEV
print("Starting loop")
do i = psminLocY-loopWidthS,psminLocY+loopWidthN
  print("Doing outer loop #: "+(i-psminLocY+loopWidthS)+" of "+(loopWidthS*2)+"   ("+((tofloat(i-psminLocY+loopWidthS)/tofloat(loopWidthS*2))*100.)+"%)")
  do j = psminLocX-loopWidthW,psminLocX+loopWidthE
    gcdist = tofloat(gc_latlon(psminlat,psminlon,lat(i),lon(j),2,4))
    if (gcdist .le. outerRad)
      ;do k = 0,dimsizes(lev)-1
      ;  rad_thevar_hit(k,{gcdist}) = rad_thevar_hit(k,{gcdist}) + 1
      ;  rad_thevar_cum(k,{gcdist}) = rad_thevar_cum(k,{gcdist}) + data(k,i,j)
      ;end do
      ;if (.not. any(ismissing(data(:,i,j)))) then
        rad_thevar_hit(:,{gcdist}) = rad_thevar_hit(:,{gcdist}) + 1
        rad_thevar_cum(:,{gcdist}) = rad_thevar_cum(:,{gcdist}) + data(:,i,j)
      ;end if
    end if
  end do
end do

; check for points where there were no hits
rad_thevar_hit = where(rad_thevar_hit.eq.0,rad_thevar_hit@_FillValue,rad_thevar_hit)

print("Minimum number of hits per gridbox: "+min(rad_thevar_hit))
print("Maximum number of hits per gridbox: "+max(rad_thevar_hit))
rad_thevar = rad_thevar_cum/rad_thevar_hit

copy_VarCoords(rad_thevar_hit,rad_thevar)

rad_thevar@units = data@units
rad_thevar@long_name=data@long_name

;rad_thevar_hit&radius=rad_thevar_hit&radius*1000.
return(rad_thevar)

end

; ====================================================================================================

undef ( "radialAvg2D" )
function radialAvg2D ( \
  data[*][*]       : numeric, \
  lat[*]           : numeric, \
  lon[*]           : numeric, \
  psminlat         : numeric, \
  psminlon         : numeric, \
  outerRad         : numeric, \
  mergeInnerBins   : logical \
)
local None

; lat (deg)
; lon (deg)
; center_lat(deg)
; center_lon(deg)
; outerRad (km)

begin

deltaLat = lat(1)-lat(0)
deltaLon = lon(1)-lon(0)
deltaMax = max((/deltaLon, deltaLat/))
deltaMin = min((/deltaLon, deltaLat/))
kmInDeg = 111.32
kmGrid = kmInDeg * deltaMax
print("The max lat/lon km grid spacing at equator is "+kmGrid+" km")

nlat = dimsizes(lat)
nlon = dimsizes(lon)

pi = 3.14159
d2r = pi/180.
r2d = 180./pi
lonr = lon*d2r
latr = lat*d2r
copy_VarCoords(lon,lonr)
copy_VarCoords(lat,latr)

psminLocY = ind_nearest_coord(psminlat,lat,0)
psminLocX = ind_nearest_coord(psminlon,lon,0)

;;;;; DO RADIAL AVERAGES

loopWidth = toint(ceil((outerRad*1.10) / (kmInDeg * deltaMin)))
; for standard lat/lon grid
;1.25 gives n ~= 5 in inner grid if numMerge = 2
;1.4 gives n ~= 9 in inner grid if numMerge = 2
timesGrid = 1.1 ; we want each radius bin to be timeGrid times kmGrid

nx = toint(outerRad/(timesGrid*kmGrid))
print("Number of bins is equal to "+nx)

if (mergeInnerBins) then
;  origRadiusArr = fspan(0,outerRad,nx)
;  radiusArr = new(dimsizes(origRadiusArr)-1,double)
;  radiusArr(0) = origRadiusArr(0)
;  radiusArr(1:dimsizes(radiusArr)-1) = origRadiusArr(2:dimsizes(origRadiusArr)-1)
  numMerge = 2 ; numMerge is the number of innermost radial bins we want to merge
  print("Merging innermost "+numMerge+" bins because radial average is so small.")
  numMergeMinusOne = numMerge-1
  origRadiusArr = fspan(0,outerRad,nx+numMergeMinusOne)
  radiusArr = new((dimsizes(origRadiusArr)-numMergeMinusOne),double)
  radiusArr(0) = origRadiusArr(0)
  printVarSummary(radiusArr)
  printVarSummary(origRadiusArr)
  radiusArr(1:dimsizes(radiusArr)-1) = origRadiusArr(numMerge:dimsizes(origRadiusArr)-1)
else
  print("Not merging any innermost bins -- be careful that your inner bins have > 1 pt.")
  radiusArr = fspan(0,outerRad,nx)
end if

numRadBins = dimsizes(radiusArr)

rad_thevar_hit  = new((/numRadBins/),integer)
rad_thevar_cum  = new((/numRadBins/),double)

; NEEDS to be 1 or 0
rad_thevar_hit!0 = "radius"
rad_thevar_cum!0 = "radius"
rad_thevar_hit&radius = radiusArr
rad_thevar_cum&radius = radiusArr

rad_thevar_hit = 0
rad_thevar_cum = 0

; NEED TO CORRECT NLEV
print("Starting loop")
do i = psminLocY-loopWidth,psminLocY+loopWidth
  print("Doing outer loop #: "+(i-psminLocY+loopWidth)+" of "+(loopWidth*2)+"   ("+((tofloat(i-psminLocY+loopWidth)/tofloat(loopWidth*2))*100.)+"%)")
  do j = psminLocX-loopWidth,psminLocX+loopWidth
      gcdist = tofloat(gc_latlon(psminlat,psminlon,lat(i),lon(j),2,4))
      if (gcdist .le. outerRad)
        rad_thevar_hit({gcdist}) = rad_thevar_hit({gcdist}) + 1
        rad_thevar_cum({gcdist}) = rad_thevar_cum({gcdist}) + data(i,j)
      end if
  end do
end do

print("Minimum number of hits per gridbox: "+min(rad_thevar_hit))
print("Maximum number of hits per gridbox: "+max(rad_thevar_hit))
rad_thevar = rad_thevar_cum/rad_thevar_hit

copy_VarCoords(rad_thevar_hit,rad_thevar)

rad_thevar@units = data@units
rad_thevar@long_name=data@long_name

return(rad_thevar)

end

; ====================================================================================================

undef ( "calcWindComponents" )
function calcWindComponents ( \
  U[*][*][*]       : numeric, \
  V[*][*][*]       : numeric, \
  lat[*]           : numeric, \
  lon[*]           : numeric, \
  psminlat         : numeric, \
  psminlon         : numeric, \
  movespd          : numeric, \
  movedir          : numeric \
)
local None
begin

; NOTE, MOVESPD is M/S
; NOTE, MOVEDIR is METEO DEG, 0 = N, 90 = E, 180 = S, 270 = W

pi = 3.14159
d2r = pi/180.
r2d = 180./pi
lonr = lon*d2r
latr = lat*d2r
copy_VarCoords(lon,lonr)
copy_VarCoords(lat,latr)

deltalong = lonr(:)-(psminlon*d2r)

arr1 = conform(U,sin(deltalong),2) * conform(U,cos(latr),1)
arr2 = cos(psminlat*d2r)*conform(U,sin(latr),1) - sin(psminlat*d2r)*conform(U,cos(latr),1)*conform(U,cos(deltalong),2)
dir_angle_r = atan2(arr1,arr2)
; need to conform here
dir_angle_d = dir_angle_r*r2d

copy_VarCoords(U,dir_angle_r)
copy_VarCoords(U,dir_angle_d)

eps = 10.e-6
phi_r = atan2(U,V)
phi_d = phi_r*r2d

copy_VarCoords(U,phi_r)
copy_VarCoords(U,phi_d)

; remove forward speed of TC
if (movespd .gt. 0)
  print("calcWindComponents: Adjusting U and V fields to remove TC forward motion")
  U = U - movespd*sin(movedir*d2r)
  V = V - movespd*cos(movedir*d2r)
end if

WIND = U
WIND = sqrt(U^2+V^2)

v_theta = WIND*sin(dir_angle_r - phi_r)
v_rad = -WIND*cos(dir_angle_r - phi_r)

copy_VarCoords(U,v_theta)
copy_VarCoords(U,v_rad)

v_theta@units = U@units
v_theta@long_name="Tangential wind"
v_rad@units = U@units
v_rad@long_name="Radial wind"

return( [/ v_rad, v_theta /] )

end


undef ( "radialInt3D" )
function radialInt3D ( \
  data[*][*][*]    : numeric, \
  lat[*]           : numeric, \
  lon[*]           : numeric, \
  lev[*]           : numeric, \
  psminlat         : numeric, \
  psminlon         : numeric, \
  outerRad         : numeric \
)
local None

; lat (deg)
; lon (deg)
; center_lat(deg)
; center_lon(deg)
; outerRad (km)

begin

deltaLat = lat(1)-lat(0)
deltaLon = lon(1)-lon(0)
deltaMax = max((/deltaLon, deltaLat/))
deltaMin = min((/deltaLon, deltaLat/))
kmInDeg = 111.32
kmGrid = kmInDeg * deltaMax
print("The max lat/lon km grid spacing at equator is "+kmGrid+" km")

nlat = dimsizes(lat)
nlon = dimsizes(lon)
nlev = dimsizes(lev)

pi = 3.14159
d2r = pi/180.
r2d = 180./pi
lonr = lon*d2r
latr = lat*d2r
copy_VarCoords(lon,lonr)
copy_VarCoords(lat,latr)

psminLocY = ind_nearest_coord(psminlat,lat,0)
psminLocX = ind_nearest_coord(psminlon,lon,0)

; Get cosine weights
rad    = 4.0*atan(1.0)/180.0
clat   = cos(lat*rad)

loopWidth = toint(ceil((outerRad*1.25) / (kmInDeg * deltaMin)))
; make sure we are more than 1 to do a few rows of buffer cells

rad_thevar_hit  = new((/nlev/),float)
rad_thevar_cum  = new((/nlev/),float)

rad_thevar_hit!0 = "lev"
rad_thevar_cum!0 = "lev"
rad_thevar_hit&lev = lev
rad_thevar_cum&lev = lev

rad_thevar_hit = 0
rad_thevar_cum = 0

; NEED TO CORRECT NLEV
print("Starting loop")
do i = psminLocY-loopWidth,psminLocY+loopWidth
  print("Doing outer loop #: "+(i-psminLocY+loopWidth)+" of "+(loopWidth*2)+"   ("+((tofloat(i-psminLocY+loopWidth)/tofloat(loopWidth*2))*100.)+"%)")
  do j = psminLocX-loopWidth,psminLocX+loopWidth
    do k = 0,dimsizes(lev)-1
      gcdist = tofloat(gc_latlon(psminlat,psminlon,lat(i),lon(j),2,4))
      if (gcdist .le. outerRad)
        rad_thevar_hit(k) = rad_thevar_hit(k) + data(k,i,j)*clat(i)
        rad_thevar_cum(k) = rad_thevar_cum(k) + clat(i)
      end if
    end do
  end do
end do

rad_thevar = rad_thevar_hit/rad_thevar_cum

copy_VarCoords(rad_thevar_hit,rad_thevar)

rad_thevar@units = data@units
rad_thevar@long_name=data@long_name

return(rad_thevar)

end



undef ( "calc_vertmax_from_rad" )
function calc_vertmax_from_rad ( \
  vrad[*][*]       : numeric \
)
local None

begin

lev = vrad&lev
nlev = dimsizes(lev)
rad = vrad&radius
nrad = dimsizes(rad)

rmw = new((/nlev/),"float")

do k = 0,dimsizes(lev)-1
  maxix = maxind(vrad(k,:))
  rmw(k) = tofloat(rad(maxix))
end do

rmw!0 = "lev"
rmw&lev = lev

return(rmw)

end














undef ( "radialAvg2D_unstruc" )
function radialAvg2D_unstruc ( \
  data[*]          : numeric, \
  lat[*]           : numeric, \
  lon[*]           : numeric, \
  deltaMax         : numeric, \
  psminlat         : numeric, \
  psminlon         : numeric, \
  outerRad         : numeric, \
  mergeInnerBins   : logical \
)
local None

; lat (deg)
; lon (deg)
; center_lat(deg)
; center_lon(deg)
; outerRad (km)

begin

kmInDeg = 111.32
kmGrid = kmInDeg * deltaMax
print("The max lat/lon km grid spacing at equator is "+kmGrid+" km")

ncol = dimsizes(lat)

pi = 3.14159
d2r = pi/180.
r2d = 180./pi
lonr = lon*d2r
latr = lat*d2r
copy_VarCoords(lon,lonr)
copy_VarCoords(lat,latr)

;;;;; DO RADIAL AVERAGES

; for standard lat/lon grid
;1.25 gives n ~= 5 in inner grid if numMerge = 2
;1.4 gives n ~= 9 in inner grid if numMerge = 2
timesGrid = 1.1 ; we want each radius bin to be timeGrid times kmGrid

nx = toint(outerRad/(timesGrid*kmGrid))
print("Number of bins is equal to "+nx)

if (mergeInnerBins) then
;  origRadiusArr = fspan(0,outerRad,nx)
;  radiusArr = new(dimsizes(origRadiusArr)-1,double)
;  radiusArr(0) = origRadiusArr(0)
;  radiusArr(1:dimsizes(radiusArr)-1) = origRadiusArr(2:dimsizes(origRadiusArr)-1)
  numMerge = 2 ; numMerge is the number of innermost radial bins we want to merge
  print("Merging innermost "+numMerge+" bins because radial average is so small.")
  numMergeMinusOne = numMerge-1
  origRadiusArr = fspan(0,outerRad,nx+numMergeMinusOne)
  radiusArr = new((dimsizes(origRadiusArr)-numMergeMinusOne),double)
  radiusArr(0) = origRadiusArr(0)
  ;printVarSummary(radiusArr)
  ;printVarSummary(origRadiusArr)
  radiusArr(1:dimsizes(radiusArr)-1) = origRadiusArr(numMerge:dimsizes(origRadiusArr)-1)
else
  print("Not merging any innermost bins -- be careful that your inner bins have > 1 pt.")
  radiusArr = fspan(0,outerRad,nx)
end if

numRadBins = dimsizes(radiusArr)

rad_thevar_hit  = new((/numRadBins/),integer)
rad_thevar_cum  = new((/numRadBins/),double)

; NEEDS to be 1 or 0
rad_thevar_hit!0 = "radius"
rad_thevar_cum!0 = "radius"
rad_thevar_hit&radius = radiusArr
rad_thevar_cum&radius = radiusArr

rad_thevar_hit = 0
rad_thevar_cum = 0

; NEED TO CORRECT NLEV
print("Starting loop")
do i = 0,ncol-1
  gcdist = tofloat(gc_latlon(psminlat,psminlon,lat(i),lon(i),2,4))
  if (gcdist .le. outerRad)
    rad_thevar_hit({gcdist}) = rad_thevar_hit({gcdist}) + 1
    rad_thevar_cum({gcdist}) = rad_thevar_cum({gcdist}) + data(i)
  end if
end do

print("Minimum number of hits per gridbox: "+min(rad_thevar_hit))
print("Maximum number of hits per gridbox: "+max(rad_thevar_hit))
rad_thevar = rad_thevar_cum/rad_thevar_hit

copy_VarCoords(rad_thevar_hit,rad_thevar)

rad_thevar@units = data@units
rad_thevar@long_name=data@long_name
rad_thevar@rad_thevar_hit=rad_thevar_hit

return(rad_thevar)

end




undef ( "radialAvg3D_unstruc" )
function radialAvg3D_unstruc ( \
  data[*][*]       : numeric, \
  lat[*]           : numeric, \
  lon[*]           : numeric, \
  lev[*]           : numeric, \
  deltaMax         : numeric, \
  psminlat         : numeric, \
  psminlon         : numeric, \
  outerRad         : numeric, \
  mergeInnerBins   : logical \
)
local None

; lat (deg)
; lon (deg)
; center_lat(deg)
; center_lon(deg)
; outerRad (km)

begin

kmInDeg = 111.32
kmGrid = kmInDeg * deltaMax
print("The max lat/lon km grid spacing at equator is "+kmGrid+" km")

ncol = dimsizes(lat)
nlev = dimsizes(lev)

pi = 3.14159
d2r = pi/180.
r2d = 180./pi
lonr = lon*d2r
latr = lat*d2r
copy_VarCoords(lon,lonr)
copy_VarCoords(lat,latr)

print("ncol "+ncol+"  nlev "+nlev)
; for standard lat/lon grid
;1.25 gives n ~= 5 in inner grid if numMerge = 2
;1.4 gives n ~= 9 in inner grid if numMerge = 2
timesGrid = 1.1 ; we want each radius bin to be timeGrid times kmGrid

nx = toint(outerRad/(timesGrid*kmGrid))
print("Number of bins is equal to "+nx)

if (mergeInnerBins) then
;  origRadiusArr = fspan(0,outerRad,nx)
;  radiusArr = new(dimsizes(origRadiusArr)-1,double)
;  radiusArr(0) = origRadiusArr(0)
;  radiusArr(1:dimsizes(radiusArr)-1) = origRadiusArr(2:dimsizes(origRadiusArr)-1)
  numMerge = 2 ; numMerge is the number of innermost radial bins we want to merge
  print("Merging innermost "+numMerge+" bins because radial average is so small.")
  numMergeMinusOne = numMerge-1
  origRadiusArr = fspan(0,outerRad,nx+numMergeMinusOne)
  radiusArr = new((dimsizes(origRadiusArr)-numMergeMinusOne),double)
  radiusArr(0) = origRadiusArr(0)
  ;printVarSummary(radiusArr)
  ;printVarSummary(origRadiusArr)
  radiusArr(1:dimsizes(radiusArr)-1) = origRadiusArr(numMerge:dimsizes(origRadiusArr)-1)
else
  print("Not merging any innermost bins -- be careful that your inner bins have > 1 pt.")
  radiusArr = fspan(0,outerRad,nx)
end if

numRadBins = dimsizes(radiusArr)

rad_thevar_hit  = new((/nlev,numRadBins/),integer)
rad_thevar_cum  = new((/nlev,numRadBins/),double)

; NEEDS to be 1 or 0
rad_thevar_hit!1 = "radius"
rad_thevar_cum!1 = "radius"
rad_thevar_hit&radius = radiusArr
rad_thevar_cum&radius = radiusArr

rad_thevar_hit!0 = "lev"
rad_thevar_cum!0 = "lev"
rad_thevar_hit&lev = lev
rad_thevar_cum&lev = lev

rad_thevar_hit = 0
rad_thevar_cum = 0

; NEED TO CORRECT NLEV
print("Starting loop")
do i = 0,ncol-1
  gcdist = tofloat(gc_latlon(psminlat,psminlon,lat(i),lon(i),2,4))
  if (gcdist .le. outerRad)
    ;do k = 0,dimsizes(lev)-1
    ;  rad_thevar_hit(k,{gcdist}) = rad_thevar_hit(k,{gcdist}) + 1
    ;  rad_thevar_cum(k,{gcdist}) = rad_thevar_cum(k,{gcdist}) + data(k,i,j)
    ;end do
    ;if (.not. any(ismissing(data(:,i,j)))) then
      rad_thevar_hit(:,{gcdist}) = rad_thevar_hit(:,{gcdist}) + 1
      rad_thevar_cum(:,{gcdist}) = rad_thevar_cum(:,{gcdist}) + data(:,i)
    ;end if
  end if
end do

; check for points where there were no hits
rad_thevar_hit = where(rad_thevar_hit.eq.0,rad_thevar_hit@_FillValue,rad_thevar_hit)

print("Minimum number of hits per gridbox: "+min(rad_thevar_hit))
print("Maximum number of hits per gridbox: "+max(rad_thevar_hit))
rad_thevar = rad_thevar_cum/rad_thevar_hit

copy_VarCoords(rad_thevar_hit,rad_thevar)

rad_thevar@units = data@units
rad_thevar@long_name=data@long_name

;rad_thevar_hit&radius=rad_thevar_hit&radius*1000.
return(rad_thevar)

end