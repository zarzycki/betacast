f=addfile("/glade/p/cgd/tss/CTSM_datm_forcing_data/atm_forcing.datm7.cruncep_qianFill.0.5d.v7.c160715/domain.lnd.360x720_cruncep.130305.nc","r")
g = addfile("/glade/u/home/zarzycki/rda/ds633.0/e5.oper.an.sfc/198608/e5.oper.an.sfc.128_165_10u.ll025sc.1986080100_1986083123.nc","r")

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

xc = f->xc
xv = f->xv
yc = f->yc
yv = f->yv

nj=14
ni=24

print(xc(nj,ni))
;print(xv(nj,ni,:))

print(yc(nj,ni))
;print(yv(:,nj,ni))

;print(xc(0,0))


print("xc: "+min(f->xc)+" "+max(f->xc))
print("yc: "+min(f->yc)+" "+max(f->yc))
print("mask: "+min(f->mask)+" "+max(f->mask))
print("frac: "+min(f->frac)+" "+max(f->frac))
print("area: "+min(f->area)+" "+max(f->area))

