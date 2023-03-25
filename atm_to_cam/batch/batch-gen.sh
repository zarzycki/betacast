#!/bin/bash

BASEDIR=/glade/u/home/zarzycki/scratch/cam_to_cam/raw_files/

dates=( 20080831 20000918 20050818 19911122 19890822 20110806 19950806 20141013 20130910 )
grids=( ext ref ref wat ref ext ref ref wat )
exps=( 001 003 001 001 002 003 002 002 002 )
stormids=( 1279 0755 1048 0310 0236 1307 0528 1521 1354 )
files=( storm_1279/CHEY.VR28.NATL.EXT.CAM5.4CLM5.0.dtime900.cam.h2.2008-08-22-00000.nc \
  storm_0755/CORI.VR28.NATL.REF.CAM5.4CLM5.0.dtime900.003.cam.h2.2000-09-03-00000.nc \
  storm_1048/CHEY.VR28.NATL.REF.CAM5.4CLM5.0.dtime900.cam.h2.2005-08-08-00000.nc \
  storm_0310/CHEY.VR28.NATL.WAT.CAM5.4CLM5.0.dtime900.cam.h2.1991-11-20-00000.nc \
  storm_0236/CHEY.VR28.NATL.REF.CAM5.4CLM5.0.dtime900.002.cam.h2.1989-08-02-00000.nc \
  storm_1307/CORI.VR28.NATL.EXT.CAM5.4CLM5.0.dtime900.003.cam.h2.2011-07-08-00000.nc \
  storm_0528/CHEY.VR28.NATL.REF.CAM5.4CLM5.0.dtime900.002.cam.h2.1995-08-01-00000.nc \
  storm_1521/CHEY.VR28.NATL.REF.CAM5.4CLM5.0.dtime900.002.cam.h2.2014-09-20-00000.nc \
  storm_1354/CHEY.VR28.NATL.WAT.CAM5.4CLM5.0.dtime900.002.cam.h2.2013-08-26-00000.nc \
)

iter=0   # user iter to reference array indices
for DATE in "${dates[@]}"
do
  # Set things we need
  YYYYMMDDHH="$DATE"00
  RAWFILE=${BASEDIR}/${files[$iter]}
  CONFIG=${grids[$iter]}
  EXP=${exps[$iter]}
  STORMID=${stormids[$iter]}
  echo $iter"    "$YYYYMMDDHH" "$RAWFILE" "$CONFIG" "$EXP" "$STORMID

  ncl -n atm_to_cam.ncl 'datasource="CAM"' \
    numlevels=32 \
    YYYYMMDDHH=${YYYYMMDDHH} \
    'dycore="mpas"' \
    'data_filename = "'${RAWFILE}'"' \
    'wgt_filename="/glade/p/univ/upsu0032/MPAS/maps/map_era5_0.25x0.25_TO_mpasa3-60-florida_patc.nc"' \
    mpas_as_cam=False \
    'adjust_config=""' \
    compress_file=False \
    write_floats=False \
    add_cloud_vars=False \
    'mod_remap_file ="/glade/p/univ/upsu0032/MPAS/maps/map_ne0np4natlantic'${CONFIG}'.ne30x4_TO_era5_0.25x0.25_patc.nc"' \
    'mod_in_topo ="/glade/u/home/zarzycki/work/unigridFiles/ne0np4natlantic'${CONFIG}'.ne30x4/topo/topo_ne0np4natlantic'${CONFIG}'.ne30x4_smooth.nc"' \
    'model_topo_file="/glade/p/univ/upsu0032/MPAS/3km_florida/x20.835586.florida.init.nc"' \
    'se_inic = "/glade/p/univ/upsu0032/MPAS/3km_florida/x20.835586.florida.init.'${YYYYMMDDHH}'.L32.nc"'

    ((iter++))
done
