#!/bin/bash

################################################################
#### Cheyenne
################################################################
#>#PBS -N gen_nudge_betacast
#>#PBS -A P93300642
#>#PBS -l walltime=11:59:00
#>#PBS -q regular
#>#PBS -j oe
#>#PBS -l select=1:ncpus=36:mem=109GB

################################################################
#### Casper
################################################################
#PBS -N gen_nudge_betacast
#PBS -A P93300042
#PBS -l select=1:ncpus=12:mem=200GB
#PBS -l walltime=23:00:00
#PBS -q casper@casper-pbs
#PBS -j oe
################################################################

################################################################
#### PM-CPU interactive
################################################################
#>salloc --nodes 1 --qos interactive --time 01:00:00 --constraint cpu

################################################################
#### PMCPU
################################################################
#>#SBATCH --qos=regular
#>#SBATCH --time=5:00:00
#>#SBATCH --nodes=1
#>#SBATCH --ntasks-per-node=128
#>#SBATCH --constraint=cpu
################################################################

### CMZ NOTES:
# 4/2/21: Cannot find a benefit to running over multiple Casper nodes
# Cost penalty is 2x per file for 2 nodes, so double cost for same
# analysis throughput.
# Current recommendation is 36 CPUs, 36 GnuP tasks, scale mem as needed w/ gridres

set -e

SERVER_NAME=$(hostname -A)
echo $SERVER_NAME

case "$SERVER_NAME" in
  *"perlmutter"* | *"nid0"* | *"nid2"*)
    echo "Using pm-cpu"
    SERVER_CASE="pm-cpu"
    export BETACAST=/global/homes/c/czarzyck/betacast
    module load conda && conda activate betacast
    NUMCORES=12
    ;;
  *"casper"* | *"crhtc"*)
    echo "Using Casper"
    SERVER_CASE="casper"
    export BETACAST=/glade/u/home/zarzycki/betacast
    module load conda && conda activate npl
    NUMCORES=12
    ;;
  *)
    echo "Unrecognized server. Exiting."
    #exit 1
    ;;
esac

THISDIR=$PWD
echo "PWD: $PWD"
source ${THISDIR}/../../utils.sh   # Source external bash functions

### Check if namelist was passed in as an argument, if yes, overwrite
if [[ $# -gt 0 ]]; then
  NLFILE=${1}
fi

### Check if a second argument (input dates file) was provided
if [[ $# -gt 1 ]]; then
  input_dates_file=${2}
fi

echo "Read in: "
echo "NLFILE: $NLFILE"
echo "input_dates_file: $input_dates_file"

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
TIMESTAMP=`date +%s%N`
COMMANDFILE=commands.${TIMESTAMP}.txt

mkdir -p $OUTDIR

starttime=$(date -u +"%s")

#------------------------------------------------------

## Initialize an empty array to store YYYYMMDDHH values
DATES_ARRAY=()

# Number of seconds to skip for each nudging file
TIME_INCREMENT=$((3600*HR_RES))

#------------------------------------------------------

## If a dates file was passed in, use this to generate nudging for the specific forecast dates
## Will create nudging every HR_RES from init to NDAYS_PER_DATE

if [ -n "$input_dates_file" ] && [ -f "$input_dates_file" ]; then

  echo "Found input dates file: $input_dates_file"

  # Read the original dates from file
  ORIGINAL_DATES=($(cat "$input_dates_file"))
  echo "Read ${#ORIGINAL_DATES[@]} original dates from file"

  # Process each original date
  for original_date in "${ORIGINAL_DATES[@]}"; do
    # Add the original date to our final array
    DATES_ARRAY+=("$original_date")

    # Extract year, month, day, hour for date manipulation
    YYYY=${original_date:0:4}
    MM=${original_date:4:2}
    DD=${original_date:6:2}
    HH=${original_date:8:2}

    # Convert to seconds since epoch for easy date math
    base_time=$(date -u --date "${YYYY}-${MM}-${DD} ${HH}:00" +%s)

    # Calculate and add extended dates
    for ((i=1; i<=(NDAYS_PER_DATE*24/HR_RES); i++)); do
      # Calculate new time in seconds (increment by HR_RES hours each time)
      new_time=$((base_time + i*TIME_INCREMENT))

      # Convert back to YYYYMMDDHH format
      new_date=$(date -u --date "@$new_time" +"%Y%m%d%H")

      # Add to array
      DATES_ARRAY+=("$new_date")
    done
  done

  # After generating all dates in DATES_ARRAY, sort and remove duplicates
  # Convert array to sorted, unique array
  echo "Generated ${#DATES_ARRAY[@]} dates before deduplication"
  DATES_ARRAY=($(printf '%s\n' "${DATES_ARRAY[@]}" | sort -n | uniq))
  echo "Final array contains ${#DATES_ARRAY[@]} unique, sorted dates"

else

  ENDHOUR=$((24-HR_RES))

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
      YYYYMMDDHH=`date -u --date @${t} +'%Y%m%d%H'`
      echo $YYYYMMDDHH
      ## Append the current YYYYMMDDHH to the array
      DATES_ARRAY+=("$YYYYMMDDHH")
    done

  done
fi

#------------------------------------------------------

echo "DATES_ARRAY contains ${#DATES_ARRAY[@]} dates"
for date in "${DATES_ARRAY[@]}"; do
  echo "  $date"
done

#------------------------------------------------------

for f in "${DATES_ARRAY[@]}"; do

  echo $f
  YYYY=${f:0:4}
  MM=${f:4:2}
  DD=${f:6:2}
  HH=${f:8:2}
  HHBASE=10#${HH}
  SS=$(( 3600*HHBASE ))
  printf -v SSSSS "%05d" $SS

  # OUTFILETMP is what betacast writes, but then quickly finalize.
  # This gives us a check if the run stops due to wallclock, etc. and the file is only partially written...
  OUTFILE=${OUTDIR}/ndg.${DESCSTR}.${GRIDSTR}.L${NUMLEVS}.cam2.i.$YYYY-$MM-$DD-$SSSSS.nc
  OUTFILETMP=${OUTFILE}.TMP.nc

  #4/4/22 After CESM2.2, nudging.F90 moved to PIO which doesn't support compression (I don't think...) ... supports floats, which are 25GB uncompressed vs 18GB compressed
  if ${CAM_TO_CAM} ; then
    echo "Running CAM->CAM options"
    for INFILE in $BINLIST; do
      echo $INFILE

      TXTCOMMAND="python ${BETACAST}/py_atm_to_cam/atm_to_cam.py
          --datasource CAM \
          --numlevels ${NUMLEVS} \
          --YYYYMMDDHH ${f} \
          --data_filename ${INFILE} \
          --wgt_filename ${WGTNAME} \
          --dycore ${DYCORE} \
          --mpas_as_cam \
          --write_floats \
          --adjust_config \"\" \
          --model_topo_file ${BNDTOPO} \
          --mod_remap_file ${MODREMAPFILE} \
          --mod_in_topo ${MODINTOPO} \
          --se_inic ${OUTFILETMP} \
          --verbose ;
          mv -v ${OUTFILETMP} ${OUTFILE} "

      # If file doesn't exist, we want to create, but if it does let's skip
      if [ ! -f ${OUTFILE} ]; then
        echo ${TXTCOMMAND} >> ${COMMANDFILE}
      fi
    done
  else
    echo "Running ERA5 -> CAM options"
    INFILE=${RDADIR}/e5.oper.invariant/197901/e5.oper.invariant.128_129_z.ll025sc.1979010100_1979010100.nc

    TXTCOMMAND="python ${BETACAST}/py_atm_to_cam/atm_to_cam.py
        --datasource ERA5RDA \
        --numlevels ${NUMLEVS} \
        --YYYYMMDDHH ${f} \
        --data_filename ${INFILE} \
        --wgt_filename ${WGTNAME} \
        --dycore ${DYCORE} \
        --RDADIR ${RDADIR} \
        --mpas_as_cam \
        --write_floats \
        --adjust_config \"\" \
        --model_topo_file ${BNDTOPO} \
        --se_inic ${OUTFILETMP} \
        --verbose ;
        mv -v ${OUTFILETMP} ${OUTFILE} "

    # If file doesn't exist, we want to create, but if it does let's skip
    if [ ! -f ${OUTFILE} ]; then
      echo ${TXTCOMMAND} >> ${COMMANDFILE}
    fi
  fi

done

#------------------------------------------------------

if [ $dryrun = false ] ; then
  echo "Launching GNU parallel"
  # Launch GNU parallel
  echo "Running GNU parallel on $SERVER_CASE"
  case "$SERVER_CASE" in
    "pm-cpu")
      # Set up for NERSC pm-cpu
      parallel --jobs ${NUMCORES} < ${COMMANDFILE}
      ;;
    "casper")
      # Set up for NSF/NCAR Casper
      parallel --jobs ${NUMCORES} --workdir "$PWD" < ${COMMANDFILE}
      ;;
    *)
      echo "Unknown server. Exiting."
      #parallel --jobs ${NUMCORES} -u --sshloginfile $PBS_NODEFILE --workdir $PWD < ${COMMANDFILE}
      exit 1
      ;;
  esac
  # Cleanup
  rm ${COMMANDFILE}
fi

endtime=$(date -u +"%s")
tottime=$(($endtime-$starttime))
printf "${tottime}\n" >> timing.txt
