#!/bin/bash

# Script to automatically create a case, configure, build, and run I compset to spin up
# CLM or ELM for betacast runs

# Turn on error checking
set -e
source ../utils.sh

# Usage:
#./auto-script.sh MODELSYSTEM DOERA5 DATE_YYYYMMDD NMONTHS NCYCLES ANOMYEAR NORMYEAR NAMELIST
# CESM
#./auto-script.sh 0 0 20200103 36 1 -1 -1 NAMELIST.MACHINE
# E3SM
#./auto-script.sh 1 0 19960113 12 1 2018 1920 NAMELIST.MACHINE
# Model
#./auto-script.sh 1 2 19840101 0 1 -1 -1 nl.landspinup.pm-cpu

if [[ $# -ne 8 ]] ; then echo "Need 8 inputs, got $#, exiting..." ; exit ; fi
if [ ! -s $8 ]; then echo "Namelist is empty, exiting..." ; exit ; fi

### User settings
modelSystem=${1}         # 0 = CESM/E3SMv1, 1 = E3SMv2
doERA5=${2}              # Use ERA5 DATM? 0 = yes, 1 = no (internal CRUNCEP), 2 CESM/E3SM
FORECASTDATE=${3}        # What is the date you want to spin up for (00Z)
NMONTHSSPIN=${4}         # Duration of spinup (somewhere b/w 3-12 months seems reasonable)
NCYCLES=${5}
BETACAST_ANOMYEAR=${6}
BETACAST_REFYEAR=${7}
NAMELISTFILE=${8}

### Check if BETACAST_ANOMYEAR is positive integer -- if yes, add deltas, if no, nothing
if [ $BETACAST_ANOMYEAR -lt 1 ]; then
  addDeltas=1 ; echo "We are NOT adding deltas..."
else
  addDeltas=0 ; echo "We ARE adding deltas..."
fi

# Read the namelist
read_bash_nl "${NAMELISTFILE}"

# Set user things to empty strings or internal defaults if not provided
if [ -z "${USER_FSURDAT+x}" ]; then USER_FSURDAT=""; fi
if [ -z "${USER_FINIDAT+x}" ]; then USER_FINIDAT=""; fi
if [ -z "${USER_ICOMPSET+x}" ]; then USER_ICOMPSET=""; fi
if [ -z "${USER_JOB_PRIORITY+x}" ]; then USER_JOB_PRIORITY=""; fi
if [ -z "${BUILD_ONLY+x}" ]; then BUILD_ONLY=false; fi

# Purge settings
if [ -z "${FORCE_PURGE+x}" ]; then FORCE_PURGE=false; fi
if [ -z "${RUN_DIR_BASE+x}" ]; then RUN_DIR_BASE=""; fi

# Set directories to empty strings for use DATM
if [ -z "${BETACAST_DATM_FORCING_BASE+x}" ]; then BETACAST_DATM_FORCING_BASE=""; fi
if [ -z "${BETACAST_DATM_ANOMALY_BASE+x}" ]; then BETACAST_DATM_ANOMALY_BASE=""; fi

# Check bools
check_bool "BUILD_ONLY" $BUILD_ONLY
check_bool "FORCE_PURGE" $FORCE_PURGE

# Derived settings that should be same between all machines
BETACAST_DATMDOMAIN=${BETACAST}/land-spinup/gen_datm/gen-datm/
BETACAST_ANOMALIGN=1920

if [ $modelSystem -eq 0 ]; then
  echo "Using CESM"
  EXTRAFLAGS="--run-unsupported"
  COMPSET=I2000Clm50Sp
elif [ $modelSystem -eq 1 ]; then
  echo "Using E3SM"
  EXTRAFLAGS=""
  COMPSET=IELM
else
  echo "Unknown modeling system set for modelSystem: $modelSystem"
  exit 1
fi
if [[ -n "$USER_ICOMPSET" ]]; then
  echo "*** OVERRIDING ICOMPSET WITH USER_ICOMPSET!"
  COMPSET=$USER_ICOMPSET
fi

### Do not edit below this line!
### Date logic
FORECASTYEAR="${FORECASTDATE:0:4}"
FORECASTYEARM1=$((FORECASTYEAR-1))
FORECASTYEARP1=$((FORECASTYEAR+1))
echo "Trying to get CLM restart for "${FORECASTDATE}
echo "Year plus one equals: "${FORECASTYEARP1}
echo "Year minus one equals: "${FORECASTYEARM1}

if [ $NCYCLES -lt 2 ]; then
  ### NO CYCLE
  echo "CYCLE: no cycle since NCYCLES = ${NCYCLES}"
else
  ### CYCLE
  echo "CYCLE: since NCYCLES = ${NCYCLES}, update NMONTHSSPIN from: ${NMONTHSSPIN} to..."
  NMONTHSSPIN=$((NMONTHSSPIN*NCYCLES))
  echo "....... ${NMONTHSSPIN}"
fi

### Print diagnostics
STARTDATE=`date -d "${FORECASTDATE} - ${NMONTHSSPIN} months" "+%Y-%m-%d"`
STARTYEAR=`date -d "${FORECASTDATE} - ${NMONTHSSPIN} months" "+%Y"`
echo "Starting at: "${STARTDATE}
if [ $doERA5 -eq 0 ]; then
  echo "Using ERA5 DATM"
  ERA5STYR=1990
  ERA5ENYR=2022
elif [ $doERA5 -eq 2 ]; then
  echo "Using Hyperion DATM"
  ERA5STYR=1984
  ERA5ENYR=2014
else
  echo "Using CRUNCEP DATM"
fi

if [ ${#FORECASTDATE} -ne 8 ]; then
  echo "Incorrect string length for FORECASTDATE, needs to be YYYYMMDD"
  echo "You provided $FORECASTDATE"
  echo "STOP"
  exit 1
fi

### Configure, build, run land model w/ DATM
if [ $doERA5 -eq 1 ] && (( FORECASTYEAR > 2016 )); then
  echo "No default DATM files beyond 2016"
  echo "Need to either find a different DATM set or spin up when coupled"
  echo "STOP"
  exit 1
fi

### ERA5 and other external data-checking
if [ $doERA5 -eq 0 ] && (( STARTYEAR < ${ERA5STYR} )); then
  echo "No ERA5 DATM files earlier than ${ERA5STYR}"
  echo "You provided $STARTYEAR for a start year when accounting for spinup"
  echo "Need to either find a different DATM set or spin up when coupled"
  echo "STOP"
  exit 1
elif [ $doERA5 -eq 0 ] && (( FORECASTYEAR > ${ERA5ENYR} )); then
  echo "No ERA5 DATM files later than ${ERA5ENYR}"
  echo "You provided $FORECASTYEAR for an end year (forecast)"
  echo "Need to either find a different DATM set or spin up when coupled"
  echo "STOP"
  exit 1
fi

ICASENAME=${ICASENAME/RESSTRING/$RESOL}
ICASENAME=${ICASENAME}_${FORECASTDATE}_$(printf "%04d\n" $NMONTHSSPIN)
if [ $addDeltas -eq 0 ]; then
  ICASENAME=${ICASENAME}_${BETACAST_ANOMYEAR}
fi

# Check if the case already exists and decide what to do
if [ -d "${PATHTOCASE}/${ICASENAME}" ]; then
  echo "Directory ${PATHTOCASE}/${ICASENAME} already exists!"
  if [ "$FORCE_PURGE" = "true" ]; then
    echo "FORCE_PURGE is true. Deleting the directory..."
    rm -vrf "${PATHTOCASE}/${ICASENAME}"
    # CMZ note, this will also purge the scratch dir. Ideally we'd auto get this
    # path from CIME somehow, but for now ask user to input it and if they don't
    # rely on manual CIME prompts to purge.
    if [ -n "$RUN_DIR_BASE" ]; then
      echo "RUN_DIR_BASE is set to $RUN_DIR_BASE. Deleting the run directory..."
      rm -vrf "${RUN_DIR_BASE}/${ICASENAME}"
    else
      echo "RUN_DIR_BASE is not set. Doing nothing."
    fi
  else
    echo "ERROR: Directory ${PATHTOCASE}/${ICASENAME} exists and FORCE_PURGE is not true (it is $FORCE_PURGE)."
    exit 1
  fi
else
  echo "Directory ${PATHTOCASE}/${ICASENAME} does not exist. Proceeding."
fi

### Put a block to check everything here?
echo "--------------------------------------------"
echo "modelSystem: "${modelSystem}
echo "FORECASTDATE: "${FORECASTDATE}
echo "NMONTHSSPIN: "${NMONTHSSPIN}
echo "NCYCLES: "${NCYCLES}
echo "STARTDATE: "${STARTDATE}
echo "addDeltas: "${addDeltas}
echo "BETACAST_ANOMYEAR: "${BETACAST_ANOMYEAR}
echo "BETACAST_ANOMALIGN: "${BETACAST_ANOMALIGN}
echo "BETACAST_DATMDOMAIN: "${BETACAST_DATMDOMAIN}
echo "BETACAST_DATM_FORCING_BASE: "${BETACAST_DATM_FORCING_BASE}
echo "BETACAST_DATM_ANOMALY_BASE: "${BETACAST_DATM_ANOMALY_BASE}
echo "CIMEROOT: "${CIMEROOT}
echo "PATHTOCASE: "${PATHTOCASE}
echo "ICASENAME: "${ICASENAME}
echo "PROJECT: "${PROJECT}
echo "MACHINE: "${MACHINE}
echo "NNODES: "${NNODES}
echo "RESOL: "${RESOL}
echo "RUNQUEUE: "${RUNQUEUE}
echo "USER_JOB_PRIORITY: "${USER_JOB_PRIORITY}
echo "WALLCLOCK: "${WALLCLOCK}
echo "BUILD_ONLY: "${BUILD_ONLY}
echo "ICASENAME: "${ICASENAME}
echo "USER_FSURDAT: "${USER_FSURDAT}
echo "USER_FINIDAT: "${USER_FINIDAT}
echo "COMPSET: "${COMPSET}
echo "--------------------------------------------"
sleep 10  # sleep to hold this on the interactive window for 10 sec

cd ${CIMEROOT}/cime/scripts
./create_newcase --case ${PATHTOCASE}/${ICASENAME} --compset ${COMPSET} --res ${RESOL} --mach ${MACHINE} --project ${PROJECT} ${EXTRAFLAGS}
cd ${PATHTOCASE}/${ICASENAME}
./xmlchange NTASKS=-${NNODES}
./xmlchange NTASKS_ATM=-$((NNODES-1))   # NOTE: weird errors on Cheyenne w/ equal nodes for all components, but this works?
./xmlchange NTASKS_ESP=1
./xmlchange NTASKS_IAC=1
./xmlchange DATM_MODE=CLMCRUNCEPv7
./xmlchange STOP_N=10
./xmlchange STOP_OPTION='nyears'
./xmlchange DATM_CLMNCEP_YR_START=${FORECASTYEARM1}
./xmlchange DATM_CLMNCEP_YR_END=${FORECASTYEAR}
./xmlchange DATM_CLMNCEP_YR_ALIGN=${FORECASTYEARM1}
./xmlchange RUN_STARTDATE=${STARTDATE}
./xmlchange REST_OPTION='end'
./xmlchange DOUT_S=FALSE
# If NMONTHSSPIN is 0, doesn't make sense to stop model on same day
if [ $NMONTHSSPIN -gt 0 ]; then
  ./xmlchange STOP_DATE=${FORECASTDATE}
else
  echo "NMONTHSSPIN is zero (--> $NMONTHSSPIN), no update to STOP_DATE"
fi

### If using ERA5, add the stream files and reset DATM_CLMNCEP_YR_START, etc.
if [ $doERA5 -eq 0 ]; then
  echo "Injecting ERA5 DATM streams"
  cp ${BETACAST}/land-spinup/streams/user_datm.streams.txt.CLMCRUNCEPv7* .
  #REPLACEDIR
  sed -i "s?\${BETACAST_STREAMBASE}?${BETACAST_DATM_FORCING_BASE}?g" user_datm.streams.txt.CLMCRUNCEPv7.Solar
  sed -i "s?\${BETACAST_DATMDOMAIN}?${BETACAST_DATMDOMAIN}?g" user_datm.streams.txt.CLMCRUNCEPv7.Solar
  sed -i "s?\${BETACAST_STREAMBASE}?${BETACAST_DATM_FORCING_BASE}?g" user_datm.streams.txt.CLMCRUNCEPv7.Precip
  sed -i "s?\${BETACAST_DATMDOMAIN}?${BETACAST_DATMDOMAIN}?g" user_datm.streams.txt.CLMCRUNCEPv7.Precip
  sed -i "s?\${BETACAST_STREAMBASE}?${BETACAST_DATM_FORCING_BASE}?g" user_datm.streams.txt.CLMCRUNCEPv7.TPQW
  sed -i "s?\${BETACAST_DATMDOMAIN}?${BETACAST_DATMDOMAIN}?g" user_datm.streams.txt.CLMCRUNCEPv7.TPQW
  ./xmlchange DATM_CLMNCEP_YR_ALIGN=${ERA5STYR}
elif [ $doERA5 -eq 2 ]; then

  VARS=("Precip" "Solar" "TPQW")
  for VAR in "${VARS[@]}"; do

    NEW_STREAM=./user_datm.streams.txt.CLMCRUNCEPv7.${VAR}
    DATMVARFOLDER="${BETACAST_DATM_FORCING_BASE}/${VAR}/"

    cp -v ${BETACAST}/land-spinup/streams/GENERIC/user_datm.streams.txt.CLMCRUNCEPv7.${VAR} ${NEW_STREAM}

    TMPFILELIST=${VAR}_filelist.txt

    ls -1 "$DATMVARFOLDER" > ${TMPFILELIST}
    # Temporary files for before and after the placeholder
    BEFORE_TEMP=$(mktemp)
    AFTER_TEMP=$(mktemp)

    # Split the XML file at "FILES_HERE"
    awk "/FILES_HERE/{exit} 1" "$NEW_STREAM" > "$BEFORE_TEMP"
    awk "x==1{print} /FILES_HERE/{x=1}" "$NEW_STREAM" > "$AFTER_TEMP"

    # Concatenate: before the placeholder, file list, and after the placeholder
    cat "$BEFORE_TEMP" > "$NEW_STREAM"
    cat ${TMPFILELIST} >> "$NEW_STREAM"
    cat "$AFTER_TEMP" >> "$NEW_STREAM"

    # Remove any trailing newlines at the end of the file
    sed -i '/^$/d' "$NEW_STREAM"

    # Clean up temporary files
    rm -v "$BEFORE_TEMP" "$AFTER_TEMP" "$TMPFILELIST"

    #### Other stuff
    sed -i "s?\${BETACAST_STREAMBASE}?${DATMVARFOLDER}?g" "$NEW_STREAM"
    sed -i "s?\${BETACAST_DATMDOMAIN}?${BETACAST_DATMDOMAIN}?g" "$NEW_STREAM"
  done

  ./xmlchange DATM_CLMNCEP_YR_ALIGN=${ERA5STYR}

fi

if { [ $doERA5 -eq 0 ] || [ $doERA5 -eq 2 ]; } && [ $NCYCLES -lt 2 ]; then
  ### NO CYCLE
  ./xmlchange DATM_CLMNCEP_YR_START=${ERA5STYR}
  ./xmlchange DATM_CLMNCEP_YR_END=${ERA5ENYR}
  # Update general vars in case needed for anom stream overwrite
  FORECASTYEARM1=${ERA5STYR}
  FORECASTYEAR=${ERA5ENYR}
fi


if [ $addDeltas -eq 0 ]; then
  echo "Injecting anomaly DATM streams"
  cp ${BETACAST}/land-spinup/streams/user_datm.streams.txt.Anomaly.* .
  #REPLACEDIR
  sed -i "s?\${BETACAST_DATM_ANOMALY_BASE}?${BETACAST_DATM_ANOMALY_BASE}?g" user_datm.streams.txt.Anomaly.Forcing.Humidity
  sed -i "s?\${BETACAST_DATMDOMAIN}?${BETACAST_DATMDOMAIN}?g" user_datm.streams.txt.Anomaly.Forcing.Humidity
  sed -i "s?\${BETACAST_DATM_ANOMALY_BASE}?${BETACAST_DATM_ANOMALY_BASE}?g" user_datm.streams.txt.Anomaly.Forcing.Temperature
  sed -i "s?\${BETACAST_DATMDOMAIN}?${BETACAST_DATMDOMAIN}?g" user_datm.streams.txt.Anomaly.Forcing.Temperature
  sed -i "s?\${BETACAST_DATM_ANOMALY_BASE}?${BETACAST_DATM_ANOMALY_BASE}?g" user_datm.streams.txt.Anomaly.Forcing.Longwave
  sed -i "s?\${BETACAST_DATMDOMAIN}?${BETACAST_DATMDOMAIN}?g" user_datm.streams.txt.Anomaly.Forcing.Longwave
  sed -i "s?\${BETACAST_DATM_ANOMALY_BASE}?${BETACAST_DATM_ANOMALY_BASE}?g" user_datm.streams.txt.Anomaly.Forcing.Precip
  sed -i "s?\${BETACAST_DATMDOMAIN}?${BETACAST_DATMDOMAIN}?g" user_datm.streams.txt.Anomaly.Forcing.Precip

  if [ $BETACAST_REFYEAR -gt 0 ]; then
    # run ncl to normalize things
    echo "Running with normalized deltas"
    ncl ${BETACAST}/land-spinup/normalize-datm-deltas.ncl 'current_year='${BETACAST_REFYEAR}'' 'basedir="'${BETACAST_DATM_ANOMALY_BASE}'"'
    sed -i "s?ens_QBOT_anom.nc?ens_QBOT_${BETACAST_REFYEAR}ref_anom.nc?g"   user_datm.streams.txt.Anomaly.Forcing.Humidity
    sed -i "s?ens_TBOT_anom.nc?ens_TBOT_${BETACAST_REFYEAR}ref_anom.nc?g"   user_datm.streams.txt.Anomaly.Forcing.Temperature
    sed -i "s?ens_PRECT_anom.nc?ens_PRECT_${BETACAST_REFYEAR}ref_anom.nc?g" user_datm.streams.txt.Anomaly.Forcing.Longwave
    sed -i "s?ens_QBOT_anom.nc?ens_QBOT_${BETACAST_REFYEAR}ref_anom.nc?g"   user_datm.streams.txt.Anomaly.Forcing.Precip
  else
    sed -i "s?ens_QBOT_anom.nc?ens_QBOT_anom.nc?g"   user_datm.streams.txt.Anomaly.Forcing.Humidity
    sed -i "s?ens_TBOT_anom.nc?ens_TBOT_anom.nc?g"   user_datm.streams.txt.Anomaly.Forcing.Temperature
    sed -i "s?ens_PRECT_anom.nc?ens_PRECT_anom.nc?g" user_datm.streams.txt.Anomaly.Forcing.Longwave
    sed -i "s?ens_QBOT_anom.nc?ens_QBOT_anom.nc?g"   user_datm.streams.txt.Anomaly.Forcing.Precip
  fi

  # Need to replace pres aero stream in some cases where it is transient
  cp ${BETACAST}/land-spinup/streams/user_datm.streams.txt.presaero.clim_2000 .
  sed -i "s?\${BETACAST}?${BETACAST}?g" user_datm.streams.txt.presaero.clim_2000

  cp ${BETACAST}/land-spinup/streams/user_nl_datm .
  sed -i "s?\${FORECASTYEARM1}?${FORECASTYEARM1}?g" user_nl_datm
  sed -i "s?\${FORECASTYEAR}?${FORECASTYEAR}?g" user_nl_datm
  sed -i "s?\${BETACAST_ANOMALIGN}?${BETACAST_ANOMALIGN}?g" user_nl_datm
  sed -i "s?\${BETACAST_ANOMYEAR}?${BETACAST_ANOMYEAR}?g" user_nl_datm
fi

### USER! Edit this block if using ELM and need to inject any ELM specific mods (e.g., fsurdat, etc.)
cat > user_nl_elm <<EOF
do_transient_pfts = .false.
check_finidat_fsurdat_consistency = .false.
hist_avgflag_pertape='A','I'
hist_nhtfrq = 0,-1
hist_mfilt = 1,24
hist_fincl2 = 'TSA','TBOT','RAIN','SNOW','QBOT','WIND','SWdown','LWdown'
EOF

### USER! Edit this block if using CLM and need to inject any CLM specific mods (e.g., fsurdat, etc.)
cat > user_nl_clm <<EOF
!check_finidat_fsurdat_consistency = .false.
!use_init_interp = .true.
!do_transient_pfts = .false.
hist_avgflag_pertape='A','I'
hist_nhtfrq = 0,-1
hist_mfilt = 1,24
hist_fincl2 = 'TSA','TBOT','RAIN','SNOW','QBOT','WIND','SWdown','LWdown'
EOF

if [ $modelSystem -eq 0 ]; then
  rm -v user_nl_elm
  modelSystemString="clm"
elif [ $modelSystem -eq 1 ]; then
  rm -v user_nl_clm
  modelSystemString="elm"
else
  echo "Unknown modeling system set for modelSystem: $modelSystem"
  exit 1
fi

## Do any injection into the remaining user_nl* file
if [[ -n "$USER_FSURDAT" ]]; then
  sed -i '/.*fsurdat/d' user_nl_${modelSystemString}
  echo "fsurdat='${USER_FSURDAT}'" >> user_nl_${modelSystemString}
fi

if [[ -n "$USER_FINIDAT" ]]; then
  sed -i '/.*finidat/d' user_nl_${modelSystemString}
  echo "finidat='${USER_FINIDAT}'" >> user_nl_${modelSystemString}
else
  #echo "finidat=''" >> user_nl_${modelSystemString}
  ./xmlchange ${modelSystemString^^}_FORCE_COLDSTART="on"
fi

./case.setup

echo "Checking input data"
./preview_namelists
set +e ; ./check_input_data
RESULT=$?
if [ $RESULT -ne 0 ]; then
  echo "Something went wrong with the ERA5 input data!"
  exit 1
else
  echo "Data checks out!"
fi
set -e

./case.build
./xmlchange JOB_WALLCLOCK_TIME=${WALLCLOCK}
./xmlchange CHARGE_ACCOUNT=${PROJECT}
./xmlchange --force JOB_QUEUE=${RUNQUEUE}
if [[ -n "$USER_JOB_PRIORITY" ]]; then
  echo "Setting job priority based on user preference!"
  ./xmlchange --force JOB_PRIORITY=${USER_JOB_PRIORITY}
fi
if [ "$BUILD_ONLY" = false ]; then
  ./case.submit
fi

exit 0
