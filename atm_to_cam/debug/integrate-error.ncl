f = addfile("diff.nc","r")

rho = f->rho(0,:,:,:)
lat = f->lat

rad    = 4.0*atan(1.0)/180.0
re     = 6371220.0
rr     = re*rad
clat   = cos(lat*rad)


rho_clat = wgt_areaave_Wrap(rho, clat, 1.0, 1)

print(""+rho_clat)