CAM_TO_CAM=true
dryrun=false

STYR=2013
ENYR=2013
STMON="sep"
STDAY=9
ENMON="sep"
ENDAY=19
HR_RES=6
DESCSTR="CAM5"

BETACASTDIR=/glade/u/home/zarzycki/betacast/
OUTDIR=/glade/p/univ/upsu0032/MPAS/ndg/

DYCORE="mpas"
GRIDSTR=mpasa3-60-florida
BNDTOPO=/glade/u/home/zarzycki/scratch/MPAS_OLD/3km_florida/x20.835586.florida.init.nc
WGTNAME=/glade/u/home/zarzycki/betacast/remapping/map_era5_0.25x0.25_TO_mpasa3-60-florida_patc.nc
NUMLEVS=32

## Specific settings for CAM_TO_CAM
SUBNAME=CHEY.VR28.NATL.WAT.CAM5.4CLM5.0.dtime900.002
BINLIST="/glade/u/home/zarzycki/scratch/cam_to_cam/raw_files/storm_1354/CHEY.VR28.NATL.WAT.CAM5.4CLM5.0.dtime900.002.cam.h2.2013-08-26-00000.nc"
MODREMAPFILE=/glade/u/home/zarzycki/betacast/remapping/map_ne0np4natlantic\$GRIDLOWER.ne30x4_TO_era5_0.25x0.25_patc.nc
MODINTOPO=/glade/u/home/zarzycki/work/unigridFiles/ne0np4natlantic\$GRIDLOWER.ne30x4/topo/topo_ne0np4natlantic\$GRIDLOWER.ne30x4_smooth.nc
