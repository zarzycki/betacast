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
#PBS -l walltime=4:00:00
#PBS -q casper
#PBS -j oe

### CMZ NOTES:
# 4/2/21: Cannot find a benefit to running over multiple Casper nodes
# Cost penalty is 2x per file for 2 nodes, so double cost for same
# analysis throughput.
# Current recommendation is 36 CPUs, 36 GnuP tasks, scale mem as needed w/ gridres

set -e

THISDIR=$PWD
source ${THISDIR}/../../utils.sh   # Source external bash functions

### Check if namelist was passed in as an argument, if yes, overwrite
if [[ $# -gt 0 ]] ; then
  NLFILE=${1}
fi

if [[ "$NLFILE" != /* ]] && [[ "$NLFILE" != ~* ]]; then NLFILE=${PWD}/${NLFILE}; fi
echo $NLFILE
read_bash_nl "${NLFILE}"

if ${CAM_TO_CAM} ; then
  echo "Doing CAM_TO_CAM, everything needs to be set in the namelist!"
else
  echo "Doing reanalysis, writing default input if not passed in via $NLFILE"
  if [ -z ${STDAY+x} ]; then STDAY=1; fi
  if [ -z ${ENDAY+x} ]; then ENDAY=31; fi
  if [ -z ${STMON+x} ]; then STMON="jan"; fi
  if [ -z ${ENMON+x} ]; then ENMON="dec"; fi
  if [ -z ${HR_RES+x} ]; then HR_RES=1; fi
  if [ -z ${DESCSTR+x} ]; then DESCSTR="ERA5"; fi
fi

# This block prints all timing inputs.
# The set -/+u causes the code to exit with unbound vars
set -u
echo "CHECKING FOR UNBOUND VARIABLES!"
echo "STYR: "$STYR"  ENYR: "$ENYR
echo "STDAY: "$STDAY"  ENDAY: "$ENDAY"  STMON: "$STMON"  ENMON: "$ENMON"  HR_RES: "$HR_RES
echo "DESCSTR: "$DESCSTR
set +u

# Special block of code to handle Hyperion data
if [[ -v SUBNAME ]] ; then
  GRID=`echo $SUBNAME | cut -d. -f 4`
  GRIDLOWER=$(echo $GRID | sed 's/./\L&/g' )
  MODREMAPFILE="${MODREMAPFILE//\$GRIDLOWER/$GRIDLOWER}"
  MODINTOPO="${MODINTOPO//\$GRIDLOWER/$GRIDLOWER}"
fi

## Append some subdir information for folders
OUTDIR=${OUTDIR}/${DYCORE}_${GRIDSTR}_L${NUMLEVS}/

if [[ -v SUBNAME ]] ; then
  OUTDIR=${OUTDIR}/${SUBNAME}/
else
  OUTDIR=${OUTDIR}/ERA5/
fi

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
    OUTFILE=${OUTDIR}/ndg.${DESCSTR}.${GRIDSTR}.L${NUMLEVS}.cam2.i.$YYYY-$MM-$DD-$SSSSS.nc
    OUTFILETMP=${OUTFILE}.TMP.nc
    
    #4/4/22 After CESM2.2, nudging.F90 moved to PIO which doesn't support compression (I don't think...) ... supports floats, which are 25GB uncompressed vs 18GB compressed   
    if ${CAM_TO_CAM} ; then
      echo "Running CAM->CAM options"
      for INFILE in $BINLIST; do
        echo $INFILE
        
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
        # If file doesn't exist, we want to create, but if it does let's skip
        if [ ! -f ${OUTFILE} ]; then
          echo ${NCLCOMMAND} >> ${COMMANDFILE}
        fi
      done
    else
      echo "Running ERA5 -> CAM options"
      INFILE=${RDADIR}/e5.oper.invariant/197901/e5.oper.invariant.128_129_z.ll025sc.1979010100_1979010100.nc
      NCLCOMMAND="cd ${BETACASTDIR}/atm_to_cam/ ; 
          ncl -n atm_to_cam.ncl 
          'datasource=\"ERA5RDA\"' 
          'RDADIR=\"'${RDADIR}'\"' 
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
      # If file doesn't exist, we want to create, but if it does let's skip
      if [ ! -f ${OUTFILE} ]; then
        echo ${NCLCOMMAND} >> ${COMMANDFILE}
      fi
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
