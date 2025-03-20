CAM_TO_CAM=false
dryrun=false

STYR=2021
ENYR=2021
HR_RES=6
STMON="jun"
STDAY=20
ENMON="jun"
ENDAY=28

####### PATHS ETC PER MACHINE

# BETACASTDIR=/glade/u/home/zarzycki/betacast/
# OUTDIR=/glade/p/univ/upsu0032/MPAS/ndg/
# RDADIR=/glade/campaign/collections/rda/data/ds633.0/

BETACASTDIR=/global/homes/c/czarzyck/betacast/
OUTDIR=/pscratch/sd/c/czarzyck/ndg/
RDADIR=/global/cfs/projectdirs/m3522/cmip6/ERA5/

####### MODEL CONFIGS

# DYCORE="se"
# GRIDSTR=ne30pg3
# BNDTOPO=/glade/p/cesmdata/cseg/inputdata/atm/cam/topo/se/ne30pg3_nc3000_Co060_Fi001_PF_nullRR_Nsw042_20171014.nc
# WGTNAME=/glade/u/home/zarzycki/betacast/remapping/map_era5_0.25x0.25_TO_ne30pg3_patc.nc
# NUMLEVS=58

# DYCORE="se"
# GRIDSTR=ne30pg3
# BNDTOPO=/glade/p/cesmdata/cseg/inputdata/atm/cam/topo/se/ne30pg3_nc3000_Co060_Fi001_PF_nullRR_Nsw042_20171014.nc
# WGTNAME=/glade/u/home/zarzycki/betacast/remapping/map_era5_0.25x0.25_TO_ne30pg3_patc.nc
# NUMLEVS=58

# DYCORE="fv"
# GRIDSTR=f09
# BNDTOPO=/glade/p/cesmdata/cseg/inputdata/atm/cam/topo/fv_0.9x1.25_nc3000_Nsw042_Nrs008_Co060_Fi001_ZR_sgh30_24km_GRNL_c170103.nc
# WGTNAME=/glade/u/home/zarzycki/work/maps/gfsmaps/map_era5_0.25_TO_fv0.9x1.25_patc.nc
# NUMLEVS=32

# DYCORE="mpas"
# GRIDSTR=mp120a
# BNDTOPO=/glade/p/cesmdata/cseg/inputdata/atm/cam/inic/mpas/mpasa120.CFSR.L32.nc
# WGTNAME=/glade/u/home/zarzycki/betacast/remapping/map_gfs_0.25x0.25_TO_mpasa120_patc.nc
# NUMLEVS=32

# DYCORE="mpas"
# GRIDSTR=mpasa3-60-florida
# BNDTOPO=/glade/p/univ/upsu0032/MPAS/3km_florida/x20.835586.florida.init.nc
# WGTNAME=/glade/u/home/zarzycki/betacast/remapping/map_era5_0.25x0.25_TO_mpasa3-60-florida_patc.nc
# NUMLEVS=32

# DYCORE="se"
# GRIDSTR=ne0np4natlanticref.ne30x4
# BNDTOPO=/glade/u/home/zarzycki/work/unigridFiles/ne0np4natlanticref.ne30x4/topo/topo_ne0np4natlanticref.ne30x4_smooth.nc
# WGTNAME=/glade/u/home/zarzycki/work/maps/gfsmaps/map_era5-0.25_TO_ne0np4natlanticref.ne30x4_patc.nc
# NUMLEVS=32

# DYCORE="se"
# GRIDSTR=conus_30_x8
# BNDTOPO=/global/homes/c/czarzyck/m2637/betacast/cesmfiles/topo/topo_conus_30_x8_smooth.nc
# WGTNAME=/global/homes/c/czarzyck/m2637/betacast/sewx/maps/map_era0.25_TO_conus_30_x8_patc.nc
# NUMLEVS=72

DYCORE="se"
GRIDSTR=ne30np4
BNDTOPO=/global/cfs/cdirs/e3sm/inputdata/atm/cam/topo/USGS-gtopo30_ne30np4_16xdel2-PFC-consistentSGH.nc
WGTNAME=/global/homes/c/czarzyck/m2637/betacast/sewx/maps/map_era0.25_TO_ne30np4_patc.nc
NUMLEVS=72
