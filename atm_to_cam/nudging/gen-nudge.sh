#!/bin/bash
##=======================================================================
#PBS -N gen_nudge_betacast 
#PBS -A UPSU0032 
#PBS -l walltime=11:59:00
#PBS -q regular
#PBS -j oe
#PBS -l select=1:ncpus=36:mem=109GB
################################################################

YYYY=2018
OUTDIR=/glade/u/home/zarzycki/scratch/nudge-E3SM/
BETACASTDIR=/glade/u/home/zarzycki/betacast/

DYCORE="se"
GRIDSTR=ne30
BNDTOPO=/glade/p/cesmdata/cseg/inputdata/atm/cam/topo/se/ne30np4_nc3000_Co060_Fi001_PF_nullRR_Nsw042_20171020.nc
WGTNAME=/glade/u/home/zarzycki/work/maps/gfsmaps/map_era5-0.25_TO_ne30np4_patc.nc
NUMLEVS=72

#DYCORE="fv"
#GRIDSTR=f09
#BNDTOPO=/glade/p/cesmdata/cseg/inputdata/atm/cam/topo/fv_0.9x1.25_nc3000_Nsw042_Nrs008_Co060_Fi001_ZR_sgh30_24km_GRNL_c170103.nc
#WGTNAME=/glade/u/home/zarzycki/work/maps/gfsmaps/map_era5_0.25_TO_fv0.9x1.25_patc.nc

#------------------------------------------------------
# DO NOT EDIT BELOW THIS LINE  
#------------------------------------------------------

THISDIR=${PWD}

# GNUPARALLEL SETTINGS
module load parallel
NUMCORES=12
TIMESTAMP=`date +%s%N`
COMMANDFILE=commands.${TIMESTAMP}.txt

#------------------------------------------------------

mkdir -p $OUTDIR

#### Get dates
start=$(date -u --date '1 jan '${YYYY}' 0:00' +%s)
stop=$(date -u --date '31 dec '${YYYY}' 18:00' +%s)

shopt -s nullglob
for t in $(seq ${start} 21600 ${stop})
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
  
  RDADIR=/glade/collections/rda/data/
  INFILE=${RDADIR}/ds633.0/e5.oper.invariant/197901/e5.oper.invariant.128_129_z.ll025sc.1979010100_1979010100.nc
  OUTFILE=${OUTDIR}/ndg.ERA5.${GRIDSTR}.L${NUMLEVS}.cam2.i.$YYYY-$MM-$DD-$SSSSS.nc
      
  NCLCOMMAND="cd ${BETACASTDIR}/atm_to_cam/ ; ncl -n atm_to_cam.ncl 'datasource=\"ERA5RDA\"' compress_file=True numlevels=${NUMLEVS} YYYYMMDDHH=${YYYYMMDDHH} 'dycore = \"'${DYCORE}'\"' 'data_filename = \"'${INFILE}'\"' 'wgt_filename=\"'${WGTNAME}'\"' 'model_topo_file=\"'${BNDTOPO}'\"' 'adjust_config=\"-\"' 'se_inic = \"'${OUTFILE}'\"'     "

  echo ${NCLCOMMAND} >> ${COMMANDFILE}
  
done

# Launch GNU parallel
parallel --jobs ${NUMCORES} -u --sshloginfile $PBS_NODEFILE --workdir $PWD < ${COMMANDFILE}

# Cleanup
rm ${COMMANDFILE}
