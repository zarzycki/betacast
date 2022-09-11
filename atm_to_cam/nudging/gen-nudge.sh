#!/bin/bash
#>################################################################
#>#PBS -N gen_nudge_betacast 
#>#PBS -A P93300642 
#>#PBS -l walltime=11:59:00
#>#PBS -q regular
#>#PBS -j oe
#>#PBS -l select=1:ncpus=36:mem=109GB
#>################################################################

#PBS -N gen_nudge_betacast
#PBS -A P93300642
#PBS -l select=1:ncpus=36:mem=220GB
#PBS -l walltime=24:00:00
#PBS -q casper
#PBS -j oe

### CMZ NOTES:
# 4/2/21: Cannot find a benefit to running over multiple Casper nodes
# Cost penalty is 2x per file for 2 nodes, so double cost for same
# analysis throughput.
# Current recommendation is 36 CPUs, 36 GnuP tasks, scale mem as needed w/ gridres

CAM_TO_CAM=true   # Set to false for reanalysis -> CAM, otherwise true
dryrun=false       # Use for debugging -- will generate parallel file but will not execute

if ${CAM_TO_CAM} ; then
  STYR=1989
  ENYR=1989
  STMON="aug"
  STDAY=8
  ENMON="aug"
  ENDAY=22
  HR_RES=6
else
  # NOTE, DO NOT NEED TO TOUCH!
  STYR=2005
  ENYR=2005
  STDAY=1
  ENDAY=31
  STMON="jan"
  ENMON="dec"
  HR_RES=1
fi

#OUTDIR=/glade/u/home/zarzycki/scratch/nudge-E3SM/
OUTDIR=/glade/scratch/zarzycki/ndg/

BETACASTDIR=/glade/u/home/zarzycki/betacast/

#DYCORE="se"
#GRIDSTR=ne30pg3
#BNDTOPO=/glade/p/cesmdata/cseg/inputdata/atm/cam/topo/se/ne30pg3_nc3000_Co060_Fi001_PF_nullRR_Nsw042_20171014.nc
#WGTNAME=/glade/u/home/zarzycki/betacast/remapping/map_era5_0.25x0.25_TO_ne30pg3_patc.nc
#NUMLEVS=58

#DYCORE="fv"
#GRIDSTR=f09
#BNDTOPO=/glade/p/cesmdata/cseg/inputdata/atm/cam/topo/fv_0.9x1.25_nc3000_Nsw042_Nrs008_Co060_Fi001_ZR_sgh30_24km_GRNL_c170103.nc
#WGTNAME=/glade/u/home/zarzycki/work/maps/gfsmaps/map_era5_0.25_TO_fv0.9x1.25_patc.nc
#NUMLEVS=32

#DYCORE="mpas"
#GRIDSTR=mp120a
#BNDTOPO=/glade/p/cesmdata/cseg/inputdata/atm/cam/inic/mpas/mpasa120.CFSR.L32.nc
#WGTNAME=/glade/u/home/zarzycki/betacast/remapping/map_gfs_0.25x0.25_TO_mpasa120_patc.nc
#NUMLEVS=32

# Example of CAM->CAM
DYCORE="mpas"
GRIDSTR=mpasa3-60-florida
BNDTOPO=/glade/u/home/zarzycki/scratch/MPAS/3km_florida/x20.835586.florida.init.nc
####L58 BNDTOPO=/glade/u/home/zarzycki/work/cesmfiles/inic/x20.835586.florida.init.nc
NUMLEVS=32
WGTNAME=/glade/u/home/zarzycki/betacast/remapping/map_era5_0.25x0.25_TO_mpasa3-60-florida_patc.nc

#INFILE=/glade/scratch/zarzycki/cam_to_cam/CHEY.VR28.NATL.EXT.CAM5.4CLM5.0.dtime900.cam.h2.2008-08-22-00000.nc
#MODREMAPFILE=/glade/u/home/zarzycki/betacast/remapping/map_ne0np4natlanticext.ne30x4_TO_era5_0.25x0.25_patc.nc
#MODINTOPO=/glade/u/home/zarzycki/work/unigridFiles/ne0np4natlanticext.ne30x4/topo/topo_ne0np4natlanticext.ne30x4_smooth.nc

#INFILE=/glade/scratch/zarzycki/cam_to_cam/CHEY.VR28.NATL.EXT.CAM5.4CLM5.0.dtime900.cam.h2.2008-08-22-00000.nc
#MODREMAPFILE=/glade/u/home/zarzycki/betacast/remapping/map_ne0np4natlanticref.ne30x4_TO_era5_0.25x0.25_patc.nc
#MODINTOPO=/glade/u/home/zarzycki/work/unigridFiles/ne0np4natlanticref.ne30x4/topo/topo_ne0np4natlanticref.ne30x4_smooth.nc

INFILE=/glade/u/home/zarzycki/scratch/cam_to_cam/CORI.VR28.NATL.WAT.CAM5.4CLM5.0.dtime900.003.cam.h2.1989-08-02-00000.nc
MODREMAPFILE=/glade/u/home/zarzycki/betacast/remapping/map_ne0np4natlanticwat.ne30x4_TO_era5_0.25x0.25_patc.nc
MODINTOPO=/glade/u/home/zarzycki/work/unigridFiles/ne0np4natlanticwat.ne30x4/topo/topo_ne0np4natlanticwat.ne30x4_smooth.nc

THISDIR=${PWD}

## Append some subdir information for folders
OUTDIR=${OUTDIR}/${DYCORE}_${GRIDSTR}_L${NUMLEVS}/

# GNUPARALLEL SETTINGS
module load parallel
module load peak_memusage
NUMCORES=36
TIMESTAMP=`date +%s%N`
COMMANDFILE=commands.${TIMESTAMP}.txt

mkdir -p $OUTDIR

starttime=$(date -u +"%s")

#------------------------------------------------------

ENDHOUR=$((24-HR_RES))
TIME_INCREMENT=$((3600*HR_RES))
echo "ENDHOUR: "$ENDHOUR"    TIME_INCREMENT: "$TIME_INCREMENT

## now loop through the above array
for jj in $(seq $STYR $ENYR);
do

  YYYY=${jj}

  #### Get dates
  start=$(date -u --date ''${STDAY}' '${STMON}' '${YYYY}' 0:00' +%s)
  stop=$(date -u --date ''${ENDAY}' '${ENMON}' '${YYYY}' '${ENDHOUR}':00' +%s)

  shopt -s nullglob
  for t in $(seq ${start} ${TIME_INCREMENT} ${stop})
  do
    start=`date +%s`
    f=`date -u --date @${t} +'%Y%m%d%H'`
    echo $f
    YYYY=${f:0:4}
    MM=${f:4:2}
    DD=${f:6:2}
    HH=${f:8:2}
    HHBASE=10#${HH}
    SS=$(( 3600*HHBASE ))
    printf -v SSSSS "%05d" $SS
  
    ## Begin betacast block
    YYYYMMDDHH=${f}

    # OUTFILETMP is what betacast writes, but then quickly finalize.
    # This gives us a check if the run stops due to wallclock, etc. and the file is only partially written...
    OUTFILE=${OUTDIR}/ndg.ERA5.${GRIDSTR}.L${NUMLEVS}.cam2.i.$YYYY-$MM-$DD-$SSSSS.nc
    OUTFILETMP=${OUTFILE}.TMP.nc
    
    #4/4/22 After CESM2.2, nudging.F90 moved to PIO which doesn't support compression (I don't think...) ... supports floats, which are 25GB uncompressed vs 18GB compressed   
    if ${CAM_TO_CAM} ; then
      echo "Running CAM->CAM options"
      NCLCOMMAND="cd ${BETACASTDIR}/atm_to_cam/ ; 
          ncl -n atm_to_cam.ncl 
          'datasource=\"CAM\"' 
          write_floats=True 
          add_cloud_vars=False 
          compress_file=False 
          mpas_as_cam=True 
          numlevels=${NUMLEVS} 
          YYYYMMDDHH=${YYYYMMDDHH} 
          'dycore = \"'${DYCORE}'\"' 
          'data_filename = \"'${INFILE}'\"' 
          'wgt_filename=\"'${WGTNAME}'\"' 
          'model_topo_file=\"'${BNDTOPO}'\"' 
          'mod_remap_file=\"'${MODREMAPFILE}'\"' 
          'mod_in_topo=\"'${MODINTOPO}'\"' 
          'adjust_config=\"\"' 
          'se_inic = \"'${OUTFILETMP}'\"' ; 
          mv -v ${OUTFILETMP} ${OUTFILE} "
    else
      echo "Running ERA5 -> CAM options"
      RDADIR=/glade/collections/rda/data/ds633.0/
      INFILE=${RDADIR}/e5.oper.invariant/197901/e5.oper.invariant.128_129_z.ll025sc.1979010100_1979010100.nc
      NCLCOMMAND="cd ${BETACASTDIR}/atm_to_cam/ ; 
          ncl -n atm_to_cam.ncl 
          'datasource=\"ERA5RDA\"' 
          write_floats=True 
          add_cloud_vars=False 
          compress_file=False 
          mpas_as_cam=True 
          numlevels=${NUMLEVS} 
          YYYYMMDDHH=${YYYYMMDDHH} 
          'dycore = \"'${DYCORE}'\"' 
          'data_filename = \"'${INFILE}'\"' 
          'wgt_filename=\"'${WGTNAME}'\"' 
          'model_topo_file=\"'${BNDTOPO}'\"' 
          'adjust_config=\"-\"' 
          'se_inic = \"'${OUTFILETMP}'\"' ; 
          mv -v ${OUTFILETMP} ${OUTFILE} "
    fi
    # If file doesn't exist, we want to create, but if it does let's skip
    if [ ! -f ${OUTFILE} ]; then
      echo ${NCLCOMMAND} >> ${COMMANDFILE}
    fi
  done

done

if [ $dryrun = false ] ; then
  # Launch GNU parallel
  #parallel --jobs ${NUMCORES} -u --sshloginfile $PBS_NODEFILE --workdir $PWD < ${COMMANDFILE}
  peak_memusage.exe parallel --jobs ${NUMCORES} --workdir $PWD < ${COMMANDFILE}

  # Cleanup
  rm ${COMMANDFILE}
fi

endtime=$(date -u +"%s")
tottime=$(($endtime-$starttime))
printf "${tottime}\n" >> timing.txt
