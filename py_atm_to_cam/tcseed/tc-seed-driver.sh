#!/bin/bash

conda activate betacast

SCRIPTPATH=$(dirname "$(realpath "$0")")
## Drop last two folders of path to get back to base betacast dir where utils.sh lives
SCRIPTPATH=`echo $SCRIPTPATH | rev | cut -d'/' -f3- | rev`
echo "Our script path is $SCRIPTPATH"
source ${SCRIPTPATH}/utils.sh   # Source external bash functions

if [[ $# -eq 0 ]] ; then echo "no namelist, exiting..." ; exit ; fi
if [ ! -s $1 ]; then echo "File is empty, exiting..." ; exit ; fi

CLEAN_MODE=false
if [[ $# -eq 2 && "${2,,}" == "clean" ]]; then
    CLEAN_MODE=true
    echo "Clean mode enabled - will perform cleanup operations"
fi

CONFIG_NAMELIST=$1

TIMEOUT=86400 # 24 hours

read_bash_nl "$CONFIG_NAMELIST"

############################### SET UP SOME THINGS ###############################

TEMPESTTMP=${TEMPESTTMP}.$(date +"%s%N")

#misc settings
WHICHSED="sed"

if [ -z ${CIMEsubstring+x} ]; then CIMEsubstring=""; fi
if [ -z ${CIMEbatchargs+x} ]; then CIMEbatchargs=""; fi
if [ -z ${NSEEDS+x} ]; then NSEEDS=0; fi
if [ -z ${MAXTIME+x} ]; then MAXTIME=99999999; fi

############################### RUN MODEL ###############################

if $CLEAN_MODE; then
  echo "Performing cleanup operations..."
  rm -fv ${path_to_rundir}/*.nc
  rm -fv ${path_to_rundir}/*.log.*
  rm -fv ${TEMPESTTMP}
  pushd "${path_to_case}" || { echo "Failed to change to case directory" >&2; exit 1; }
  ./xmlchange -v CONTINUE_RUN=FALSE
  ARCHIVE_PATH=$(./xmlquery DOUT_S_ROOT | cut -d: -f2 | tr -d ' ')
  rm -rfv $ARCHIVE_PATH
  popd
fi

# Turn on error checking...
set -e
set -u

if [ -f .${KILLSTR} ]; then
  echo "Exiting because of done..."
  echo "If this is a surprise to you, you need to rm .SEED_DONE in the driver folder"
  exit
else
  echo "If you want to kill this process, run..."
  echo "touch .$KILLSTR"
fi

echo "Setting up settings!"

cd ${path_to_case}

# turn off auto ST archiving
./xmlchange -v DOUT_S=FALSE

echo "Begin call to model run"

# Get number of log .gz files for sleeping
echo "Running again!" > ${path_to_rundir}/testrunning.gz
numlogfiles=`ls ${path_to_rundir}/*.gz | wc -l`
echo "numlogfiles: $numlogfiles"

echo "SUBMITTING FORECAST RUN"
set +e ; ./case.submit ${CIMEsubstring} --batch-args "${CIMEbatchargs}" ; set -e

## Set up NUKEing
if [ -f "${path_to_rundir}/NUKE" ] ; then rm -v ${path_to_rundir}/NUKE ; sleep 5 ; fi
echo "To NUKE, run \"touch ${path_to_rundir}/NUKE\" "

## Hold script while log files from filter run haven't been archived yet
START_TIME=$(date +%s)
check_timeout() {
  local current_time=$(date +%s)
  local elapsed_time=$((current_time - START_TIME))
  if [ $elapsed_time -gt $TIMEOUT ]; then
    echo "Error: Process timed out after $TIMEOUT seconds" >&2
    return 1
  fi
  return 0
}
count_log_files() {
  find "${path_to_rundir}" -maxdepth 1 -name "*.gz" -type f | wc -l
}

# Wait until job finishes
while true; do
  # Get number of files so we can see if they differ from log
  current_files=$(count_log_files)
  if [ "$current_files" != "$numlogfiles" ]; then
    echo "New log files detected. Continuing..."
    break
  fi
  if [ -f "${path_to_rundir}/NUKE" ]; then
    echo "Nuke sequence initiated, exiting betacast" >&2
    exit 1
  fi
  if ! check_timeout; then
    echo "Warning: Process took too long, exiting" >&2
    exit 1
  fi

  current_time=$(date +%s)
  elapsed=$((current_time - START_TIME))
  remaining=$((TIMEOUT - elapsed))
  printf "Sleeping... Elapsed: %3dm %2ds, Timeout in: %4dm %2ds --- %s\n" \
       $((elapsed/60)) $((elapsed%60)) $((remaining/60)) $((remaining%60)) "$(date '+%Y%m%d %H:%M:%S')"
  sleep 20
done

echo "Run over done sleeping ($(date '+%Y%m%d %H:%M:%S')) will hold for 20 more sec to make sure files moved"
sleep 20

############################### RESUB SCRIPT ###############################

cd $DRIVER_DIRECTORY

# Get most recent file for tracking
MOSTRECENT=`ls --color=never ${path_to_rundir}/*.cam.${HTRACKSTR}.*.nc -r | head -n 1`

# Get most recent YYYYMMDD
YYYYMMDD=$(get_YYYYMMDD_from_hfile "$MOSTRECENT" "$HTRACKSTR")

# Do TE to find candidates in final timeslice
STR_DETECT="--verbosity 0 --in_connect ${CONNECTDAT} --out ${TEMPESTFILE} --closedcontourcmd PSL,300.0,5.0,0;_DIFF(Z300,Z500),-6.0,5.0,1.0 --minlat -25.0 --maxlat 25.0 --mergedist 6.0 --searchbymin PSL --outputcmd PSL,min,0"
${TEMPESTEXTREMESDIR}/bin/DetectNodes --in_data ${MOSTRECENT} ${STR_DETECT}

# Get number of lines in Tempest file. NOTE, since TE includes a header line NLINES=1 means no TCs and >1 means candidates exist...
NLINES=`head ${TEMPESTFILE} | wc -l`
echo "The tempestExtremes candidate file has $NLINES lines"

# If we have more than one storm *or* no storms but do seed...
if (( $NLINES > 1 )) || ${do_seed} ; then

  # Get most recent restart file
  RESTARTFILE=`ls --color=never ${path_to_rundir}/*.cam.r.*.nc -r | head -n 1`
  cp -v $RESTARTFILE ${path_to_rundir}/ORIG.nc

  # If NLINES>1, we have some storms in the TMP file, so load them into bash arrays
  if (( $NLINES > 1 )) ; then
    tail -n +2 ${TEMPESTFILE} > ${TEMPESTTMP}
    declare stormLat=()
    declare stormLon=()
    while IFS=$'\t' read col1 col2 col3 col4
    do
      echo "$col2"
      echo "$col3"
      stormLon+=($col2)
      stormLat+=($col3)
    done < ${TEMPESTTMP}
    echo "${stormLat[*]}"
    echo "${stormLon[*]}"
    NSTORMS=`echo "${#stormLat[@]}"`
  else
    NSTORMS=0
  fi

  # If we are seeding, the loop is number of seeds to add.
  # If we are unseeding, loop is number of storms to remove from TMP file
  if ${do_seed} ; then
    NLOOP=$NSEEDS
  else
    NLOOP=$NSTORMS
  fi
  echo "NLOOP is: $NLOOP"

  # Loop over either how many seeds we need to add *or* how many TCs we need to unseed
  for ((ii=0; ii<NLOOP; ii++)); do
    echo "doing loop index: $ii"

    # First figure out if seeding, where this seed goes
    # ... or if not seeding, get rid of ii'th storm
    if ${do_seed} ; then
      echo "WE ARE SEEDING!"
      python random-seed.py --filename ${TEMPESTTMP} --pthi ${vortex_namelist}
    else
      echo "sedding lat/lon: ${stormLat[$ii]} ${stormLon[$ii]}"
      $WHICHSED -i "s?.*psminlat=.*?psminlat=${stormLat[$ii]}?" ${vortex_namelist}
      $WHICHSED -i "s?.*psminlon=.*?psminlon=${stormLon[$ii]}?" ${vortex_namelist}
    fi

    # Find TC parameters
    python find-tc-fill-params.py \
        --inic_file ${RESTARTFILE} \
        --vortex_namelist ${vortex_namelist}

    # Seed TC in data
    python py-seed-tc-in-ncdata.py \
        --se_inic ${RESTARTFILE} \
        --vortex_namelist ${vortex_namelist}

  done

  # Remove TMP traj file
  rm -v ${TEMPESTTMP}
else
  # If here it means we are *not* seeding and there were no detected storms (i.e., NLINES=1) that need filling
  echo "No seeding requested and/or no unseeding required... doing nothing!"
fi

############################### go back and do some things... ###############################

cd ${path_to_case}

## Force run to be "continued" now since we've run at least one cycle
./xmlchange -v CONTINUE_RUN=TRUE

## Run the ST archiver
if ${st_archive_ontheway} ; then
  set +e ; ./case.st_archive ; set -e
fi

## Check to see if we've exceeded our max time
if [[ ${YYYYMMDD} -ge ${MAXTIME} ]] ; then
  echo "Model time ${YYYYMMDD} has exceeded MAXTIME: ${MAXTIME}"
  echo "... EXITING!"
  exit
else
  echo "Model time ${YYYYMMDD} is before MAXTIME: ${MAXTIME}"
fi


############################### RESUB SCRIPT ###############################

cd $DRIVER_DIRECTORY ; pwd

# cleanup any BAK files
rm -fv ${vortex_namelist}.BAK

exec ./tc-seed-driver.sh "$CONFIG_NAMELIST"
