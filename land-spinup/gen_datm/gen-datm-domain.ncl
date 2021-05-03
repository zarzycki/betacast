load "functions.ncl"

g=addfile("/glade/u/home/zarzycki/rda/ds633.0/e5.oper.an.sfc/198608/e5.oper.an.sfc.128_165_10u.ll025sc.1986080100_1986083123.nc","r")

;f=addfile("/glade/p/cgd/tss/CTSM_datm_forcing_data/atm_forcing.datm7.GSWP3.0.5d.v1.c170516/domain.lnd.360x720_gswp3.0v1.c170606.nc","r")
;print(sum(f->area))

;netcdf domain.CRUJRAV1.1.0.5x0.5.c190606 {
;dimensions:
;	ni = 720 ;
;	nj = 360 ;
;	nv = 4 ;
;variables:
;	float xc(nj, ni) ;
;		xc:units = "degrees_east" ;
;		xc:long_name = "longitude" ;
;		xc:mode = "time-invariant" ;
;	float yc(nj, ni) ;
;		yc:units = "degrees_north" ;
;		yc:long_name = "latitude" ;
;		yc:mode = "time-invariant" ;
;	float xv(nv, nj, ni) ;
;		xv:units = "degrees_east" ;
;		xv:long_name = "longitude of grid cell vertices" ;
;		xv:mode = "time-invariant" ;
;	float yv(nv, nj, ni) ;
;		yv:units = "degrees_north" ;
;		yv:long_name = "latitude of grid cell vertices" ;
;		yv:mode = "time-invariant" ;
;	float mask(nj, ni) ;
;		mask:units = "unitless" ;
;		mask:long_name = "domain mask" ;
;		mask:mode = "time-invariant" ;
;	float area(nj, ni) ;
;		area:units = "radians squared" ;
;		area:long_name = "area of grid cell in radians squared" ;
;		area:mode = "time-invariant" ;


lat = tofloat(g->latitude(::-1))
lon = tofloat(g->longitude)

nlat = dimsizes(lat)
nlon = dimsizes(lon)

xc = conform_dims((/nlat,nlon/),lon,1)
yc = conform_dims((/nlat,nlon/),lat,0)

opt=True
opt@rearth=1.0
area = tofloat(area_global_rectilinear_grid(lat,lon,opt))

xc!0="nj"
xc!1="ni"

yc!0="nj"
yc!1="ni"

area!0="nj"
area!1="ni"

thismask=area
thismask=1.0

fullfilename="era5-domain.nc"

system("/bin/rm -f "+fullfilename)   ; remove any pre-existing file
ncdf = addfile(fullfilename ,"c")  ; open output netCDF file

;===================================================================
; create global attributes of the file (optional)
;===================================================================
fAtt               = True            ; assign file attributes
fAtt@title         = "NCL Simple Approach to netCDF Creation"
fAtt@source_file   =  "original-file.nc"
fAtt@Conventions   = "None"
fAtt@creation_date = systemfunc ("date")
fileattdef( ncdf, fAtt )            ; copy file attributes

;===================================================================
; output variables directly; NCL will call appropriate functions
; to write the meta data associated with each variable
;===================================================================
ncdf->xc = xc
ncdf->yc = yc
ncdf->area = area
ncdf->mask = thismask
