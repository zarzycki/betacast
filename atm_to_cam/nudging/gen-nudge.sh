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
#PBS -A UNSB0017
#PBS -l select=1:ncpus=36:mem=220GB
#PBS -l walltime=24:00:00
#PBS -q casper
#PBS -j oe

### CMZ NOTES:
# 4/2/21: Cannot find a benefit to running over multiple Casper nodes
# Cost penalty is 2x per file for 2 nodes, so double cost for same
# analysis throughput.
# Current recommendation is 36 CPUs, 36 GnuP tasks, scale mem as needed w/ gridres

STYR=2014
ENYR=2014

#OUTDIR=/glade/u/home/zarzycki/scratch/nudge-E3SM/
OUTDIR=/glade/scratch/zarzycki/ndg/

BETACASTDIR=/glade/u/home/zarzycki/betacast/

DYCORE="se"
GRIDSTR=ne30
BNDTOPO=/glade/p/cesmdata/cseg/inputdata/atm/cam/topo/se/ne30np4_nc3000_Co060_Fi001_PF_nullRR_Nsw042_20171020.nc
WGTNAME=/glade/u/home/zarzycki/work/maps/gfsmaps/map_era5-0.25_TO_ne30np4_patc.nc
NUMLEVS=32

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

## now loop through the above array
for jj in $(seq $STYR $ENYR);
do

  YYYY=${jj}

  #### Get dates
  start=$(date -u --date '1 jan '${YYYY}' 0:00' +%s)
  stop=$(date -u --date '31 dec '${YYYY}' 23:00' +%s)

  shopt -s nullglob
  for t in $(seq ${start} 3600 ${stop})
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
  
    RDADIR=/glade/collections/rda/data/ds633.0/
    INFILE=${RDADIR}/e5.oper.invariant/197901/e5.oper.invariant.128_129_z.ll025sc.1979010100_1979010100.nc
    
    # OUTFILETMP is what betacast writes, but then quickly finalize.
    # This gives us a check if the run stops due to wallclock, etc. and the file is only partially written...
    OUTFILE=${OUTDIR}/ndg.ERA5.${GRIDSTR}.L${NUMLEVS}.cam2.i.$YYYY-$MM-$DD-$SSSSS.nc
    OUTFILETMP=${OUTFILE}.TMP.nc
    
    #4/4/22 After CESM2.2, nudging.F90 moved to PIO which doesn't support compression (I don't think...) ... supports floats, which are 25GB uncompressed vs 18GB compressed   
    NCLCOMMAND="cd ${BETACASTDIR}/atm_to_cam/ ; ncl -n atm_to_cam.ncl 'datasource=\"ERA5RDA\"' 'RDADIR = \"'${RDADIR}'\"' write_floats=True add_cloud_vars=False compress_file=False mpas_as_cam=True numlevels=${NUMLEVS} YYYYMMDDHH=${YYYYMMDDHH} 'dycore = \"'${DYCORE}'\"' 'data_filename = \"'${INFILE}'\"' 'wgt_filename=\"'${WGTNAME}'\"' 'model_topo_file=\"'${BNDTOPO}'\"' 'adjust_config=\"-\"' 'se_inic = \"'${OUTFILETMP}'\"' ; mv -v ${OUTFILETMP} ${OUTFILE} "
    
    # If file doesn't exist, we want to create, but if it does let's skip
    if [ ! -f ${OUTFILE} ]; then
      echo ${NCLCOMMAND} >> ${COMMANDFILE}
    fi
  done

done

# Launch GNU parallel
#parallel --jobs ${NUMCORES} -u --sshloginfile $PBS_NODEFILE --workdir $PWD < ${COMMANDFILE}
peak_memusage.exe parallel --jobs ${NUMCORES} --workdir $PWD < ${COMMANDFILE}

# Cleanup
rm ${COMMANDFILE}

endtime=$(date -u +"%s")
tottime=$(($endtime-$starttime))
printf "${tottime}\n" >> timing.txt
