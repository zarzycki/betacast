#!/bin/bash
### CMZ NOTES:
# 4/2/21: Cannot find a benefit to running over multiple Casper nodes
# Cost penalty is 2x per file for 2 nodes, so double cost for same
# analysis throughput.
# Current recommendation is 36 CPUs, 36 GnuP tasks, scale mem as needed w/ gridres

echo "Command used: $0 \"$@\""
echo "PID       : $$"
echo "Host      : $(hostname)"
echo "User      : $(whoami)"
echo "Start time: $(date -u +"%Y-%m-%d %H:%M:%S UTC")"
echo "Shell     : $SHELL"

set -e

# Machine config: expect BETACAST, OUTDIR, RDADIR, NUMCORES from
# environment (set by submit wrapper). Falls back to namelist values
# if not set (backward compat for direct invocation).
if [ -z "${BETACAST+x}" ]; then
  echo "WARNING: BETACAST not set in environment, will use BETACASTDIR from namelist"
fi

if [ -n "${BETACAST:-}" ]; then
  THISDIR="${BETACAST}/py_atm_to_cam/nudging"
else
  THISDIR=$PWD
fi
echo "PWD: $PWD"
echo "THISDIR: $THISDIR"
source ${THISDIR}/../../utils.sh   # Source external bash functions

### Check if namelist was passed in as an argument, if yes, overwrite
if [[ $# -gt 0 ]]; then
  NLFILE=${1}
fi

### Check if a second argument (input dates file) was provided
if [[ $# -gt 1 ]]; then
  input_dates_file=${2}
fi

### Check if a third argument (INDEX) was provided
if [[ $# -gt 2 ]]; then
  INDEX=${3}
fi

echo "Raw args: $@"
echo "Parsed:"
echo "  NLFILE=${NLFILE:-<not set>}"
echo "  input_dates_file=${input_dates_file:-<not set>}"
echo "  INDEX=${INDEX:-<not set>}"

if [[ "$NLFILE" != /* ]] && [[ "$NLFILE" != ~* ]]; then NLFILE=${PWD}/${NLFILE}; fi
echo $NLFILE

# Save env vars before sourcing namelist so submit wrapper values take precedence
_ENV_BETACAST="${BETACAST:-}"  ; _ENV_BETACAST_SET="${BETACAST+x}"
_ENV_OUTDIR="${OUTDIR:-}"      ; _ENV_OUTDIR_SET="${OUTDIR+x}"
_ENV_RDADIR="${RDADIR:-}"      ; _ENV_RDADIR_SET="${RDADIR+x}"
_ENV_NUMCORES="${NUMCORES:-}"  ; _ENV_NUMCORES_SET="${NUMCORES+x}"

# Unset RDADIR so we can detect if the namelist explicitly sets it
unset RDADIR

read_bash_nl "${NLFILE}"

# Track whether RDADIR was explicitly set in the namelist
_NL_RDADIR_SET="${RDADIR+x}"

# Restore env vars: submit wrapper values override namelist values
# Exception: RDADIR — namelist wins if explicitly set, otherwise fall back to env
if [ -n "${_ENV_BETACAST_SET}" ]; then BETACAST="${_ENV_BETACAST}"; fi
if [ -n "${_ENV_OUTDIR_SET}" ];   then OUTDIR="${_ENV_OUTDIR}"; fi
if [ -z "${_NL_RDADIR_SET}" ] && [ -n "${_ENV_RDADIR_SET}" ]; then RDADIR="${_ENV_RDADIR}"; fi
if [ -n "${_ENV_NUMCORES_SET}" ]; then NUMCORES="${_ENV_NUMCORES}"; fi

# Set defaults
if [ -z "${add_scream+x}" ]; then add_scream="false"; fi
if [ -z "${SOURCE+x}" ]; then SOURCE="ERA5"; fi
if [ -z "${dryrun+x}" ]; then dryrun="false"; fi

# Backward compat: CAM_TO_CAM=true -> SOURCE=CAM
if [ -n "${CAM_TO_CAM+x}" ] && [ "$CAM_TO_CAM" = "true" ]; then
  echo "WARNING: CAM_TO_CAM is deprecated, use SOURCE=CAM instead"
  SOURCE="CAM"
fi

# Apply namelist fallback: BETACASTDIR -> BETACAST (backward compat)
if [ -z "${BETACAST+x}" ] && [ -n "${BETACASTDIR+x}" ]; then
  export BETACAST="${BETACASTDIR}"
fi
if [ -z "${BETACAST+x}" ]; then
  echo "ERROR: BETACAST must be set (via environment or BETACASTDIR in namelist)"; exit 1
fi
if [ -z "${OUTDIR+x}" ]; then
  echo "ERROR: OUTDIR must be set (via environment or namelist)"; exit 1
fi
if [ -z "${NUMCORES+x}" ]; then NUMCORES=12; fi

if [ "$SOURCE" = "CAM" ]; then
  # Check for required vars
  required_vars=(SUBNAME BINLIST MODREMAPFILE MODINTOPO)
  for v in "${required_vars[@]}"; do
    if [ -z "${!v+x}" ]; then
      echo "ERROR: Required SOURCE=CAM variable '$v' is not set."
      exit 1
    fi
  done
  echo "Doing CAM->CAM: STDAY, ENDAY, STMON, ENMON, HR_RES, DESCSTR need to be set in the namelist!"
elif [ "$SOURCE" = "ERA5" ] || [ "$SOURCE" = "CR20V3" ]; then
  if [ "$SOURCE" = "CR20V3" ] && [ -z "${_NL_RDADIR_SET}" ]; then
    echo "ERROR: SOURCE=CR20V3 requires RDADIR to be set in the namelist (default RDADIR points to ERA5 data)"; exit 1
  fi
  echo "Doing reanalysis (${SOURCE}), writing default input if not passed in via $NLFILE"
  if [ -z "${STDAY+x}" ]; then STDAY=1; fi
  if [ -z "${ENDAY+x}" ]; then ENDAY=31; fi
  if [ -z "${STMON+x}" ]; then STMON="jan"; fi
  if [ -z "${ENMON+x}" ]; then ENMON="dec"; fi
  if [ -z "${HR_RES+x}" ]; then HR_RES=1; fi
  if [ -z "${DESCSTR+x}" ]; then DESCSTR="${SOURCE}"; fi
else
  echo "ERROR: Unrecognized SOURCE='${SOURCE}'. Valid options: ERA5, CR20V3, CAM"; exit 1
fi

# This block prints all timing inputs.
# The set -/+u causes the code to exit with unbound vars
set -u
echo "CHECKING FOR UNBOUND VARIABLES!"
echo "STYR: "$STYR"  ENYR: "$ENYR
echo "STDAY: "$STDAY"  ENDAY: "$ENDAY"  STMON: "$STMON"  ENMON: "$ENMON"  HR_RES: "$HR_RES
echo "DESCSTR: "$DESCSTR"  DYCORE: "$DYCORE"  GRIDSTR: "$GRIDSTR"  NUMLEVS: "$NUMLEVS
echo "BNDTOPO: "$BNDTOPO
echo "WGTNAME: "$WGTNAME
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
  OUTDIR=${OUTDIR}/${SOURCE}/
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

  # These commands are needed to reorder scream nudging things
  if [ "$add_scream" = true ]; then
    # SCREAM needs dims reordered and variables renamed to match convention.
    SCREAM_CMDS="ncpdq -O -a ncol,lev ${OUTFILETMP} ${OUTFILETMP} ; ncrename -v T,T_mid -v Q,qv ${OUTFILETMP} ;"
  else
    SCREAM_CMDS=""
  fi

  #4/4/22 After CESM2.2, nudging.F90 moved to PIO which doesn't support compression (I don't think...) ... supports floats, which are 25GB uncompressed vs 18GB compressed
  if [ "$SOURCE" = "CAM" ]; then
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
          ${SCREAM_CMDS}
          mv -v ${OUTFILETMP} ${OUTFILE} "

      # If file doesn't exist, we want to create, but if it does let's skip
      if [ ! -f ${OUTFILE} ]; then
        echo ${TXTCOMMAND} >> ${COMMANDFILE}
      fi
    done
  elif [ "$SOURCE" = "ERA5" ]; then
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
        ${SCREAM_CMDS}
        mv -v ${OUTFILETMP} ${OUTFILE} "

    # If file doesn't exist, we want to create, but if it does let's skip
    if [ ! -f ${OUTFILE} ]; then
      echo ${TXTCOMMAND} >> ${COMMANDFILE}
    fi
  elif [ "$SOURCE" = "CR20V3" ]; then
    echo "Running CR20V3 -> CAM options"
    INFILE="NULL"  # not used by CR20V3 loader; RDADIR provides all paths

    TXTCOMMAND="python ${BETACAST}/py_atm_to_cam/atm_to_cam.py
        --datasource CR20V3 \
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
        ${SCREAM_CMDS}
        mv -v ${OUTFILETMP} ${OUTFILE} "

    # If file doesn't exist, we want to create, but if it does let's skip
    if [ ! -f ${OUTFILE} ]; then
      echo ${TXTCOMMAND} >> ${COMMANDFILE}
    fi
  fi

done

#------------------------------------------------------

if [ $dryrun = false ] ; then
  echo "Launching GNU parallel with ${NUMCORES} cores"
  # --halt now,fail=1: kill all jobs immediately if any single job fails
  parallel --jobs ${NUMCORES} --workdir "$PWD" --halt now,fail=1 < ${COMMANDFILE}
  rm ${COMMANDFILE}
fi

endtime=$(date -u +"%s")
tottime=$(($endtime-$starttime))
printf "${tottime}\n" >> timing.txt
