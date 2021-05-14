undef("area_global_rectilinear_grid")
function area_global_rectilinear_grid (lat[*]:numeric, lon[*], opt[1]:logical)
; Area (approximate) of a global rectilinear grid cells.
; Actually, the rectilinear grid must be equally spaced: 1x1, 1x2, 3x5
; This could be modified for non-global grids by adding
;     the left and right and left dx.
local rearth, nlat, mlon, pi, rad, rr, dlon, dx, dy
    , area_lat, area_grid, area_total, area_sphere
begin
  rearth = 6371.220d0            ; default; kilometers
  if (opt .and. isatt(opt, "rearth")) then
      rearth = opt@rearth
  end if

  nlat   = dimsizes(lat)
  mlon   = dimsizes(lon)

  pi     = 4d0*atan(1.0d0)
  rad    = pi/180d0
  rr     = rearth*rad
  dlon   = rr*(lon(2)-lon(1))

  dx     = dlon*cos(lat*rad)   ; lon spacing at each latitude

  if (lat(0).lt.-89.9999d0) then ; force exact 0.0, rounding issue
      dx(0)  = 0.0d0
  end if
  if (lat(nlat-1).gt.89.9999d0) then
      dx(nlat-1)  = 0.0d0
  end if
                               ; lat can have variable spacing
  dy           = new (nlat, "double", "No_FillValue")
  dy(0)        = (lat(1)-lat(0))*rr
  dy(1:nlat-2) = (lat(2:)-lat(1:nlat-2))*rr
  dy(nlat-1)   = (lat(nlat-1)-lat(nlat-2))*rr

  area_lat      = dx*dy        ; [*]
  area_grid     = conform_dims( (/nlat,mlon/), area_lat, 0)
  area_total    = sum( area_grid )

  area_sphere   = 4d0*pi*rearth^2

  area_grid@long_name  = "area of each grid cell"
  area_grid@units      = "km^2"
  area_grid@area_total = area_total      ; scalar
  ;area_grid@area_lat   = area_lat        ; [*]
  area_grid@rearth     = rearth          ; scalar
  area_grid@area_sphere= area_sphere     ; scalar
  area_grid@area_ratio = area_total/area_sphere    ; scalar

  area_grid!0   = "lat"
  area_grid!1   = "lon"
  area_grid&lat =  lat
  area_grid&lon =  lon

  return(area_grid)
end
