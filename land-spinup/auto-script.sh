#!/bin/bash

# Script to automatically create a case, configure, build, and run I compset to spin up
# CLM or ELM for betacast runs

# Turn on error checking
set -e
source ../utils.sh

# Usage:
#./auto-script.sh MODELSYSTEM DATAFORCING DATE_YYYYMMDD NMONTHS NCYCLES ANOMYEAR NORMYEAR NAMELIST
# CESM (historical)
#./auto-script.sh 0 0 20200103 36 1 -1 -1 nl.landspinup.derecho
# E3SM (historical)
#./auto-script.sh 1 0 20110827 12 1 -1 -1 nl.landspinup.philly.pm
# E3SM (perturbations)
#./auto-script.sh 1 0 19960113 12 1 2018 1920 nl.landspinup.pm-cpu
# Model
#./auto-script.sh 1 2 19840101 0 1 -1 -1 NAMELIST.MACHINE

if [[ $# -ne 8 ]] ; then echo "Need 8 inputs, got $#, exiting..." ; exit ; fi
if [ ! -s $8 ]; then echo "Namelist is empty, exiting..." ; exit ; fi

### User settings
modelSystem=${1}         # 0 = CESM/E3SMv1, 1 = E3SMv2
dataForcing=${2}         # DATM forcing? 0 = yes, 1 = no (internal CRUNCEP), 2 CESM/E3SM
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

# If user has specified FORCE_COLD=False, we won't force a cold start, otherwise we do!
if [ -z "${FORCE_COLD+x}" ]; then FORCE_COLD=true; fi

# Set directories to empty strings for use DATM
if [ -z "${BETACAST_DATM_FORCING_BASE+x}" ]; then BETACAST_DATM_FORCING_BASE=""; fi
if [ -z "${BETACAST_DATM_ANOMALY_BASE+x}" ]; then BETACAST_DATM_ANOMALY_BASE=""; fi

# Coupler default
if [ -z "${COUPLER+x}" ]; then COUPLER="mct"; fi

# Check bools
check_bool "BUILD_ONLY" $BUILD_ONLY
check_bool "FORCE_PURGE" $FORCE_PURGE

# Derived settings that should be same between all machines
BETACAST_DATMDOMAIN=${BETACAST}/land-spinup/gen_datm/gen-datm/
BETACAST_ANOMALIGN=1920

if [ $modelSystem -eq 0 ]; then
  echo "Using CESM"
  EXTRAFLAGS="--run-unsupported"
  #COMPSET=I2000Clm50Sp
  COMPSET=IHistClm60Sp
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

echo "CYCLE: NCYCLES = ${NCYCLES}"
NMONTHSSPIN_WITH_CYCLES=$((NMONTHSSPIN*NCYCLES))
echo "NMONTHSSPIN: ${NMONTHSSPIN} / NMONTHSSPIN_WITH_CYCLES: ${NMONTHSSPIN_WITH_CYCLES}"

### Print diagnostics
# This is the actual start date the model sees
MODEL_STARTDATE=`date -d "${FORECASTDATE} - ${NMONTHSSPIN_WITH_CYCLES} months" "+%Y-%m-%d"`
# This is the first year of the DATM stream required
DATM_STARTYEAR=`date -d "${FORECASTDATE} - ${NMONTHSSPIN} months" "+%Y"`
# Some other variables that may be helpful
FORECASTYEARM1=$((FORECASTYEAR-1))
FORECASTYEARP1=$((FORECASTYEAR+1))

echo "Trying to get CLM restart for "${FORECASTDATE}
echo "Starting model at: "${MODEL_STARTDATE}
echo "First year of DATM data (with $NCYCLES cycles) is: "${DATM_STARTYEAR}

if [ $dataForcing -eq 0 ]; then
  echo "Using ERA5 DATM"
  DATMMINYR=1980
  DATMMAXYR=2024
elif [ $dataForcing -eq 2 ]; then
  echo "Using Hyperion DATM"
  DATMMINYR=1984
  DATMMAXYR=2014
else
  echo "Using CRUNCEP DATM"
fi

if [ $NMONTHSSPIN -eq 0 ]; then
  echo "NMONTHSSPIN set to 0, we are using entire stream period $DATMMINYR to $DATMMAXYR"
  DATM_STARTYEAR=$DATMMINYR
  FORECASTYEAR=$DATMMAXYR
fi

### Error checking
if [ ${#FORECASTDATE} -ne 8 ]; then
  echo "Incorrect string length for FORECASTDATE, needs to be YYYYMMDD"
  echo "You provided $FORECASTDATE"
  echo "STOP"
  exit 1
fi
if [ $dataForcing -eq 1 ] && (( FORECASTYEAR > 2016 )); then
  echo "No default DATM files beyond 2016"
  echo "Need to either find a different DATM set or spin up when coupled"
  echo "STOP"
  exit 1
fi
if { [ $dataForcing -eq 0 ] || [ $dataForcing -eq 2 ]; } && (( DATM_STARTYEAR < ${DATMMINYR} )); then
  echo "No DATM files for dataset $dataForcing earlier than ${DATMMINYR}"
  echo "You provided $DATM_STARTYEAR for a start year when accounting for spinup"
  echo "Need to either find a different DATM set or spin up when coupled"
  echo "STOP"
  exit 1
elif { [ $dataForcing -eq 0 ] || [ $dataForcing -eq 2 ]; } && (( FORECASTYEAR > ${DATMMAXYR} )); then
  echo "No DATM files for dataset $dataForcing later than ${DATMMAXYR}"
  echo "You provided $FORECASTYEAR for an end year (forecast)"
  echo "Need to either find a different DATM set or spin up when coupled"
  echo "STOP"
  exit 1
fi

ICASENAME=${ICASENAME/RESSTRING/$RESOL}
ICASENAME=${ICASENAME}_${FORECASTDATE}_$(printf "%03d\n" $NMONTHSSPIN)_$(printf "%02d\n" $NCYCLES)
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
echo "dataForcing: "${dataForcing}
echo "FORECASTDATE: "${FORECASTDATE}
echo "NCYCLES: "${NCYCLES}
echo "NMONTHSSPIN: "${NMONTHSSPIN}
echo "NMONTHSSPIN_WITH_CYCLES: "${NMONTHSSPIN_WITH_CYCLES}
echo "MODEL_STARTDATE: "${MODEL_STARTDATE}
echo "DATM_STARTYEAR: "${DATM_STARTYEAR}
echo "FORECASTYEAR: "${FORECASTYEAR}
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
set +e ; ./xmlchange NTASKS_ESP=1 ; set -e
set +e ; ./xmlchange NTASKS_IAC=1 ; set -e
./xmlchange DATM_MODE=CLMCRUNCEPv7
./xmlchange STOP_N=20
./xmlchange STOP_OPTION='nyears'
# For now, let's try both with CLMNCEP and not in there...
# If CLMNCEP isn't there, try just DATM_ prefixes
./xmlchange DATM_CLMNCEP_YR_ALIGN=${DATM_STARTYEAR} || ./xmlchange DATM_YR_ALIGN=${DATM_STARTYEAR}
./xmlchange DATM_CLMNCEP_YR_START=${DATM_STARTYEAR} || ./xmlchange DATM_YR_START=${DATM_STARTYEAR}
./xmlchange DATM_CLMNCEP_YR_END=${FORECASTYEAR} || ./xmlchange DATM_YR_END=${FORECASTYEAR}
./xmlchange RUN_STARTDATE=${MODEL_STARTDATE}
./xmlchange REST_OPTION='end'
./xmlchange DOUT_S=FALSE
# If NMONTHSSPIN is 0, doesn't make sense to stop model on same day
if [ $NMONTHSSPIN -gt 0 ]; then
  ./xmlchange STOP_DATE=${FORECASTDATE}
else
  echo "NMONTHSSPIN is zero (--> $NMONTHSSPIN), no update to STOP_DATE"
fi

### If using ERA5, add the stream files and reset DATM_CLMNCEP_YR_START, etc.
if [ $dataForcing -eq 0 ]; then
  echo "Injecting ERA5 DATM streams"
  if [ "$COUPLER" == "mct" ]; then
    cp -v ${BETACAST}/land-spinup/streams/user_datm.streams.txt.CLMCRUNCEPv7* .
    sed -i "s?\${BETACAST_STREAMBASE}?${BETACAST_DATM_FORCING_BASE}?g" user_datm.streams.txt.CLMCRUNCEPv7.Solar
    sed -i "s?\${BETACAST_DATMDOMAIN}?${BETACAST_DATMDOMAIN}?g" user_datm.streams.txt.CLMCRUNCEPv7.Solar
    sed -i "s?\${BETACAST_STREAMBASE}?${BETACAST_DATM_FORCING_BASE}?g" user_datm.streams.txt.CLMCRUNCEPv7.Precip
    sed -i "s?\${BETACAST_DATMDOMAIN}?${BETACAST_DATMDOMAIN}?g" user_datm.streams.txt.CLMCRUNCEPv7.Precip
    sed -i "s?\${BETACAST_STREAMBASE}?${BETACAST_DATM_FORCING_BASE}?g" user_datm.streams.txt.CLMCRUNCEPv7.TPQW
    sed -i "s?\${BETACAST_DATMDOMAIN}?${BETACAST_DATMDOMAIN}?g" user_datm.streams.txt.CLMCRUNCEPv7.TPQW
  elif [ "$COUPLER" == "nuopc" ]; then
    echo "COUPLER is nuopc"
    cp -v ${BETACAST}/land-spinup/streams/nuopc/user_nl_datm_streams .
    sed -i "s?\${BETACAST_STREAMBASE}?${BETACAST_DATM_FORCING_BASE}?g" user_nl_datm_streams
    sed -i "s?\${BETACAST_DATMDOMAIN}?${BETACAST_DATMDOMAIN}?g" user_nl_datm_streams
  else
    echo "Error: COUPLER must be either 'mct' or 'nuopc'"
    exit 1
  fi
elif [ $dataForcing -eq 2 ]; then

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

  ./xmlchange DATM_CLMNCEP_YR_ALIGN=${DATMMINYR}

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
    # Run NCL to normalize things
    echo "Running with normalized deltas"
    set -x ; ncl ${BETACAST}/land-spinup/normalize-datm-deltas.ncl 'current_year='${BETACAST_REFYEAR}'' 'basedir="'${BETACAST_DATM_ANOMALY_BASE}'"' ; set +x
    # Replace default anomalies in the namelists with normalized ones
    sed -i "s?ens_QBOT_anom.nc?ens_QBOT_${BETACAST_REFYEAR}ref_anom.nc?g"   user_datm.streams.txt.Anomaly.Forcing.Humidity
    sed -i "s?ens_TBOT_anom.nc?ens_TBOT_${BETACAST_REFYEAR}ref_anom.nc?g"   user_datm.streams.txt.Anomaly.Forcing.Temperature
    sed -i "s?ens_FLDS_anom.nc?ens_FLDS_${BETACAST_REFYEAR}ref_anom.nc?g"   user_datm.streams.txt.Anomaly.Forcing.Longwave
    sed -i "s?ens_PRECT_anom.nc?ens_PRECT_${BETACAST_REFYEAR}ref_anom.nc?g" user_datm.streams.txt.Anomaly.Forcing.Precip
  fi

  # Need to replace pres aero stream in some cases where it is transient
  cp ${BETACAST}/land-spinup/streams/user_datm.streams.txt.presaero.clim_2000 .
  sed -i "s?\${BETACAST}?${BETACAST}?g" user_datm.streams.txt.presaero.clim_2000

  cp ${BETACAST}/land-spinup/streams/user_nl_datm .
  sed -i "s?\${FORECASTYEARM1}?${DATM_STARTYEAR}?g" user_nl_datm
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
flanduse_timeseries=''
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
  if [ "$FORCE_COLD" = "true" ]; then
    #echo "finidat=''" >> user_nl_${modelSystemString}
    ./xmlchange ${modelSystemString^^}_FORCE_COLDSTART="on"
  fi
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
