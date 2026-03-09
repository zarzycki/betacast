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
#./auto-script.sh 1 0 2011082612 12 1 -1 -1 nl.landspinup.philly.pm
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
  addDeltas=0 ; echo "We are NOT adding deltas..."
else
  addDeltas=1 ; echo "We ARE adding deltas..."
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

# Validate FORCE_PURGE requires RUN_DIR_BASE
if [ "$FORCE_PURGE" = "true" ] && [ -z "$RUN_DIR_BASE" ]; then
  echo "ERROR: FORCE_PURGE=true requires RUN_DIR_BASE to be set in the namelist."
  exit 1
fi

# Derived settings that should be same between all machines
BETACAST_DATMDOMAIN=${BETACAST}/land-spinup/gen_datm/gen-datm/
BETACAST_ANOMALIGN=1920

if [ $modelSystem -eq 0 ]; then
  echo "Using CESM"
  modelSystemString="clm"
  EXTRAFLAGS="--run-unsupported"
  #COMPSET=I2000Clm50Sp
  #COMPSET=IHistClm60Sp
  COMPSET="2000_DATM%GSWP3v1_CLM50%SP_SICE_SOCN_MOSART_CISM2%NOEVOLVE_SWAV"
elif [ $modelSystem -eq 1 ]; then
  echo "Using E3SM"
  modelSystemString="elm"
  EXTRAFLAGS=""
  #COMPSET=IELM
  COMPSET="2000_DATM%QIA_ELM%SP_SICE_SOCN_SROF_SGLC_SWAV_SIAC_SESP"
else
  echo "Unknown modeling system set for modelSystem: $modelSystem"
  exit 1
fi
if [[ -n "$USER_ICOMPSET" ]]; then
  echo "*** OVERRIDING ICOMPSET WITH USER_ICOMPSET!"
  COMPSET=$USER_ICOMPSET
fi



### Do not edit below this line!

### Check date passed in
date_length=${#FORECASTDATE}
if [[ $date_length -ne 8 && $date_length -ne 10 ]]; then
  echo "ERROR: FORECASTDATE must be 8 (YYYYMMDD) or 10 (YYYYMMDDHH) digits long"
  echo "STOP"
  exit 1
fi
### Date logic
FORECASTYEAR="${FORECASTDATE:0:4}"
FORECASTMON="${FORECASTDATE:4:2}"
FORECASTDAY="${FORECASTDATE:6:2}"
if [[ $date_length -eq 10 ]]; then
  FORECASTCYCLE="${FORECASTDATE:8:2}"
  FORECASTSSSSS=$((10#$FORECASTCYCLE * 3600))
else
  FORECASTCYCLE="00"
  FORECASTSSSSS="00000"
fi

echo "CYCLE: NCYCLES = ${NCYCLES}"
NMONTHSSPIN_WITH_CYCLES=$((NMONTHSSPIN*NCYCLES))
echo "NMONTHSSPIN: ${NMONTHSSPIN} / NMONTHSSPIN_WITH_CYCLES: ${NMONTHSSPIN_WITH_CYCLES}"

### Print diagnostics
# This is the actual start date the model sees
MODEL_STARTDATE=$(date -d "${FORECASTYEAR}${FORECASTMON}${FORECASTDAY} - ${NMONTHSSPIN_WITH_CYCLES} months" "+%Y-%m-%d")
# This is the first year of the DATM stream required
DATM_STARTYEAR=$(date -d "${FORECASTYEAR}${FORECASTMON}${FORECASTDAY} - ${NMONTHSSPIN} months" "+%Y")
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
  BETACAST_DATMDOMAIN_FILE="era5-domain.nc"
elif [ $dataForcing -eq 1 ]; then
  echo "Using CRUNCEP DATM"
  # Note, these are just dummy years
  DATMMINYR=1901
  DATMMAXYR=2016
  BETACAST_DATMDOMAIN_FILE=""
elif [ $dataForcing -eq 2 ]; then
  echo "Using Hyperion DATM"
  DATMMINYR=1984
  DATMMAXYR=2014
  BETACAST_DATMDOMAIN_FILE="era5-domain.nc"
elif [ $dataForcing -eq 3 ]; then
  echo "Using CR20V3"
  DATMMINYR=1850
  DATMMAXYR=2015
  BETACAST_DATMDOMAIN_FILE="cr20v3-domain.nc"
else
  echo "ERROR: Unknown dataForcing value: $dataForcing"
  exit 1
fi

if [ $NMONTHSSPIN -eq 0 ]; then
  if [ -z "${DATMMINYR}" ] || [ -z "${DATMMAXYR}" ]; then
    echo "ERROR: NMONTHSSPIN=0 requires DATMMINYR and DATMMAXYR to be set, but dataForcing=$dataForcing did not define them."
    exit 1
  fi
  echo "NMONTHSSPIN set to 0, we are using entire stream period $DATMMINYR to $DATMMAXYR"
  DATM_STARTYEAR=$DATMMINYR
  FORECASTYEAR=$DATMMAXYR
fi

### Error checking
if [ $dataForcing -eq 1 ] && (( FORECASTYEAR > 2016 )); then
  echo "No default DATM files beyond 2016"
  echo "Need to either find a different DATM set or spin up when coupled"
  echo "STOP"
  exit 1
fi
### DATM main forcing bounds checking
if { [ $dataForcing -eq 0 ] || [ $dataForcing -eq 2 ] || [ $dataForcing -eq 3 ]; } && (( DATM_STARTYEAR < ${DATMMINYR} )); then
  echo "No DATM files for dataset $dataForcing earlier than ${DATMMINYR}"
  echo "You provided $DATM_STARTYEAR for a start year when accounting for spinup"
  echo "Need to either find a different DATM set or spin up when coupled"
  echo "STOP"
  exit 1
elif { [ $dataForcing -eq 0 ] || [ $dataForcing -eq 2 ] || [ $dataForcing -eq 3 ]; } && (( FORECASTYEAR > ${DATMMAXYR} )); then
  echo "No DATM files for dataset $dataForcing later than ${DATMMAXYR}"
  echo "You provided $FORECASTYEAR for an end year (forecast)"
  echo "Need to either find a different DATM set or spin up when coupled"
  echo "STOP"
  exit 1
fi
if (( addDeltas == 1 && DATM_STARTYEAR < 1920 )); then
  echo "Anomaly in $DATM_STARTYEAR requested but no anomaly forcing (currently) before 1920."
  echo "Either edit the script or manually add your file:"
  echo "(set BETACAST_ANOMALIGN=1920 to clear this message)"
  echo "(set BUILD_ONLY=True to not submit so you can manually edit)"
  echo "STOP"
  exit 1
fi

ICASENAME=${ICASENAME/RESSTRING/$RESOL}
ICASENAME=${ICASENAME}_${FORECASTDATE}_$(printf "%03d\n" $NMONTHSSPIN)_$(printf "%02d\n" $NCYCLES)
if [ $addDeltas -eq 1 ]; then
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
echo "CIMEROOT: "${CIMEROOT}
echo "modelSystem: "${modelSystem}
echo "COUPLER: "${COUPLER}
echo "dataForcing: "${dataForcing}
echo "FORECASTDATE: "${FORECASTDATE}
echo "NCYCLES: "${NCYCLES}
echo "NMONTHSSPIN: "${NMONTHSSPIN}
echo "NMONTHSSPIN_WITH_CYCLES: "${NMONTHSSPIN_WITH_CYCLES}
echo "MODEL_STARTDATE: "${MODEL_STARTDATE}
echo "DATM_STARTYEAR: "${DATM_STARTYEAR}
echo "FORECASTYEAR: "$FORECASTYEAR
echo "FORECASTMON: "$FORECASTMON
echo "FORECASTDAY: "$FORECASTDAY
echo "FORECASTCYCLE: "$FORECASTCYCLE
echo "FORECASTSSSSS: "$FORECASTSSSSS
echo "addDeltas: "${addDeltas}
echo "BETACAST_ANOMYEAR: "${BETACAST_ANOMYEAR}
echo "BETACAST_ANOMALIGN: "${BETACAST_ANOMALIGN}
echo "BETACAST_DATMDOMAIN: "${BETACAST_DATMDOMAIN}
echo "BETACAST_DATM_FORCING_BASE: "${BETACAST_DATM_FORCING_BASE}
echo "BETACAST_DATM_ANOMALY_BASE: "${BETACAST_DATM_ANOMALY_BASE}
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
echo "USER_FSURDAT: "${USER_FSURDAT}
echo "USER_FINIDAT: "${USER_FINIDAT}
echo "COMPSET: "${COMPSET}
echo "--------------------------------------------"
sleep 10  # sleep to hold this on the interactive window for 10 sec

cd ${CIMEROOT}/cime/scripts
./create_newcase --case ${PATHTOCASE}/${ICASENAME} --compset ${COMPSET} --res ${RESOL} --mach ${MACHINE} --project ${PROJECT} ${EXTRAFLAGS}
cd ${PATHTOCASE}/${ICASENAME}
xmlchange_verbose "NTASKS" "-${NNODES}"
xmlchange_verbose "NTASKS_ATM" "-$((NNODES-1))" # NOTE: weird errors on Cheyenne w/ equal nodes for all components, but this works?
set +e ; xmlchange_verbose "NTASKS_ESP" "1" ; set -e
set +e ; xmlchange_verbose "NTASKS_IAC" "1" ; set -e
xmlchange_verbose "DATM_MODE" "CLMCRUNCEPv7"
xmlchange_verbose "STOP_N" "$NMONTHSSPIN_WITH_CYCLES"
xmlchange_verbose "STOP_OPTION" "nmonths"
# For now, let's try both with CLMNCEP and not in there...
# If CLMNCEP isn't there, try just DATM_ prefixes
xmlchange_verbose "DATM_CLMNCEP_YR_ALIGN" "$DATM_STARTYEAR" || xmlchange_verbose "DATM_YR_ALIGN" "$DATM_STARTYEAR"
xmlchange_verbose "DATM_CLMNCEP_YR_START" "$DATM_STARTYEAR" || xmlchange_verbose "DATM_YR_START" "$DATM_STARTYEAR"
xmlchange_verbose "DATM_CLMNCEP_YR_END" "$FORECASTYEAR" || xmlchange_verbose "DATM_YR_END" "$FORECASTYEAR"
xmlchange_verbose "RUN_STARTDATE" "$MODEL_STARTDATE"
xmlchange_verbose "START_TOD" "$FORECASTSSSSS"
xmlchange_verbose "REST_OPTION" "end"
xmlchange_verbose "DOUT_S" "FALSE"

# # If NMONTHSSPIN is 0, doesn't make sense to stop model on same day
# if [ $NMONTHSSPIN -gt 0 ]; then
#   ./xmlchange STOP_DATE=${FORECASTDATE}
# else
#   echo "NMONTHSSPIN is zero (--> $NMONTHSSPIN), no update to STOP_DATE"
# fi

# Only do user_datm files if "not supported" by default
if [ $dataForcing -ne 1 ] ; then
  # Now we need to do different things for different couplers
  if [ "$COUPLER" == "mct" ]; then
    # Variable streams
    VARS=("Precip" "Solar" "TPQW")

    for VAR in "${VARS[@]}"; do

      # Define output paths
      NEW_STREAM=./user_datm.streams.txt.CLMCRUNCEPv7.${VAR}
      DATMVARFOLDER="${BETACAST_DATM_FORCING_BASE}/${VAR}/"

      # Copy the generic stream over
      cp -v ${BETACAST}/land-spinup/streams/GENERIC/user_datm.streams.txt.CLMCRUNCEPv7.${VAR} ${NEW_STREAM}

      # Define a temporary filelist and store available forcing netCDF files
      TMPFILELIST=${VAR}_filelist.txt
      shopt -s nullglob

      full_file_list=false
      if [ "$full_file_list" = "true" ]; then
        # Include all NetCDF files (old method)
        files=("$DATMVARFOLDER"/*.nc)
      else
        # Include only files containing YYYY between DATM_STARTYEAR and FORECASTYEAR
        files=()
        for ((yr=DATM_STARTYEAR; yr<=FORECASTYEAR; yr++)); do
          files+=("$DATMVARFOLDER"/*.${yr}-*.nc)
        done
      fi

      # Exit if no NetCDF files are found
      if [ ${#files[@]} -eq 0 ]; then
        echo "ERROR: No NetCDF forcing files found in $DATMVARFOLDER" >&2
        exit 1
      fi

      # Write filenames (not full paths) to the file list
      printf "%s\n" "${files[@]##*/}" > "$TMPFILELIST"

      # Verify expected number of files (12 per year)
      expected_files=$(( (FORECASTYEAR - DATM_STARTYEAR + 1) * 12 ))
      actual_files=${#files[@]}
      if [ "$actual_files" -ne "$expected_files" ]; then
        echo "ERROR: Forcing file count mismatch in $DATMVARFOLDER" >&2
        echo "Expected files : $expected_files  (years ${DATM_STARTYEAR}-${FORECASTYEAR}, 12/month)" >&2
        echo "Actual files   : $actual_files" >&2
        exit 1
      else
        echo "OK: Forcing file count verified in $DATMVARFOLDER"
        echo "Files found    : $actual_files"
        echo "Expected files : $expected_files  (years ${DATM_STARTYEAR}-${FORECASTYEAR}, 12/month)"
      fi

      # Temporary files for splitting the user_datm stream before and after the placeholder
      BEFORE_TEMP=$(mktemp)
      AFTER_TEMP=$(mktemp)

      # Split the user_datm stream file at "FILES_HERE" (now we have two split files)
      awk "/FILES_HERE/{exit} 1" "$NEW_STREAM" > "$BEFORE_TEMP"
      awk "x==1{print} /FILES_HERE/{x=1}" "$NEW_STREAM" > "$AFTER_TEMP"

      # Concatenate! Glue three files together in order
      # BEFORE_TEMP, TMPFILELIST, AFTER_TEMP
      cat "$BEFORE_TEMP" > "$NEW_STREAM"    # Overwrite
      cat ${TMPFILELIST} >> "$NEW_STREAM"   # Append
      cat "$AFTER_TEMP" >> "$NEW_STREAM"    # Append

      # Remove any trailing newlines at the end of the file
      vsed -i '/^$/d' "$NEW_STREAM"

      # Clean up temporary files
      rm -v "$BEFORE_TEMP" "$AFTER_TEMP" "$TMPFILELIST"

      # Replace path placeholders
      vsed -i "s?\${BETACAST_STREAMBASE}?${DATMVARFOLDER}?g" "$NEW_STREAM"
      vsed -i "s?\${BETACAST_DATMDOMAIN}?${BETACAST_DATMDOMAIN}?g" "$NEW_STREAM"
      vsed -i "s?\${BETACAST_DATMDOMAIN_FILE}?${BETACAST_DATMDOMAIN_FILE}?g" "$NEW_STREAM"
    done

  elif [ "$COUPLER" == "nuopc" ]; then
    echo "COUPLER is nuopc"
    cp -v ${BETACAST}/land-spinup/streams/nuopc/user_nl_datm_streams .
    vsed -i "s?\${BETACAST_STREAMBASE}?${BETACAST_DATM_FORCING_BASE}?g" user_nl_datm_streams
    vsed -i "s?\${BETACAST_DATMDOMAIN}?${BETACAST_DATMDOMAIN}?g" user_nl_datm_streams
  else
    echo "Error: COUPLER must be either 'mct' or 'nuopc'"
    exit 1
  fi
fi

if [ $addDeltas -eq 1 ]; then
  echo "Injecting anomaly DATM streams"
  cp ${BETACAST}/land-spinup/streams/user_datm.streams.txt.Anomaly.* .
  #REPLACEDIR
  vsed -i "s?\${BETACAST_DATM_ANOMALY_BASE}?${BETACAST_DATM_ANOMALY_BASE}?g" user_datm.streams.txt.Anomaly.Forcing.Humidity
  vsed -i "s?\${BETACAST_DATMDOMAIN}?${BETACAST_DATMDOMAIN}?g" user_datm.streams.txt.Anomaly.Forcing.Humidity
  vsed -i "s?\${BETACAST_DATM_ANOMALY_BASE}?${BETACAST_DATM_ANOMALY_BASE}?g" user_datm.streams.txt.Anomaly.Forcing.Temperature
  vsed -i "s?\${BETACAST_DATMDOMAIN}?${BETACAST_DATMDOMAIN}?g" user_datm.streams.txt.Anomaly.Forcing.Temperature
  vsed -i "s?\${BETACAST_DATM_ANOMALY_BASE}?${BETACAST_DATM_ANOMALY_BASE}?g" user_datm.streams.txt.Anomaly.Forcing.Longwave
  vsed -i "s?\${BETACAST_DATMDOMAIN}?${BETACAST_DATMDOMAIN}?g" user_datm.streams.txt.Anomaly.Forcing.Longwave
  vsed -i "s?\${BETACAST_DATM_ANOMALY_BASE}?${BETACAST_DATM_ANOMALY_BASE}?g" user_datm.streams.txt.Anomaly.Forcing.Precip
  vsed -i "s?\${BETACAST_DATMDOMAIN}?${BETACAST_DATMDOMAIN}?g" user_datm.streams.txt.Anomaly.Forcing.Precip

  if [ $BETACAST_REFYEAR -gt 0 ]; then
    # Run NCL to normalize things
    echo "Running with normalized deltas"
    set -x ; ncl ${BETACAST}/land-spinup/normalize-datm-deltas.ncl 'current_year='${BETACAST_REFYEAR}'' 'bavsedir="'${BETACAST_DATM_ANOMALY_BASE}'"' ; set +x
    # Replace default anomalies in the namelists with normalized ones
    vsed -i "s?ens_QBOT_anom.nc?ens_QBOT_${BETACAST_REFYEAR}ref_anom.nc?g"   user_datm.streams.txt.Anomaly.Forcing.Humidity
    vsed -i "s?ens_TBOT_anom.nc?ens_TBOT_${BETACAST_REFYEAR}ref_anom.nc?g"   user_datm.streams.txt.Anomaly.Forcing.Temperature
    vsed -i "s?ens_FLDS_anom.nc?ens_FLDS_${BETACAST_REFYEAR}ref_anom.nc?g"   user_datm.streams.txt.Anomaly.Forcing.Longwave
    vsed -i "s?ens_PRECT_anom.nc?ens_PRECT_${BETACAST_REFYEAR}ref_anom.nc?g" user_datm.streams.txt.Anomaly.Forcing.Precip
  fi

  # Need to replace pres aero stream in some cases where it is transient
  cp ${BETACAST}/land-spinup/streams/user_datm.streams.txt.presaero.clim_2000 .
  vsed -i "s?\${BETACAST}?${BETACAST}?g" user_datm.streams.txt.presaero.clim_2000

  cp ${BETACAST}/land-spinup/streams/user_nl_datm .
  vsed -i "s?\${FORECASTYEARM1}?${DATM_STARTYEAR}?g" user_nl_datm
  vsed -i "s?\${FORECASTYEAR}?${FORECASTYEAR}?g" user_nl_datm
  vsed -i "s?\${BETACAST_ANOMALIGN}?${BETACAST_ANOMALIGN}?g" user_nl_datm
  vsed -i "s?\${BETACAST_ANOMYEAR}?${BETACAST_ANOMYEAR}?g" user_nl_datm
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
elif [ $modelSystem -eq 1 ]; then
  rm -v user_nl_clm
fi

## Do any injection into the remaining user_nl* file
if [[ -n "$USER_FSURDAT" ]]; then
  vsed -i '/.*fsurdat/d' user_nl_${modelSystemString}
  echo "fsurdat='${USER_FSURDAT}'" >> user_nl_${modelSystemString}
fi

if [[ -n "$USER_FINIDAT" ]]; then
  vsed -i '/.*finidat/d' user_nl_${modelSystemString}
  echo "finidat='${USER_FINIDAT}'" >> user_nl_${modelSystemString}
else
  if [ "$FORCE_COLD" = "true" ]; then
    #echo "finidat=''" >> user_nl_${modelSystemString}
    xmlchange_verbose "${modelSystemString^^}_FORCE_COLDSTART" "on"
  fi
fi

./case.setup

echo "Checking input data"
./preview_namelists
set +e ; ./check_input_data
RESULT=$?
if [ $RESULT -ne 0 ]; then
  echo "UH OH... Something went wrong with the input data!"
  exit 1
else
  echo "Data checks out!"
fi
set -e

./case.build
xmlchange_verbose "JOB_WALLCLOCK_TIME" "$WALLCLOCK"
xmlchange_verbose "CHARGE_ACCOUNT" "$PROJECT"
xmlchange_verbose "JOB_QUEUE" "$RUNQUEUE" "--force"
if [[ -n "$USER_JOB_PRIORITY" ]]; then
  echo "Setting job priority based on user preference!"
  xmlchange_verbose "JOB_PRIORITY" "$USER_JOB_PRIORITY" "--force"
fi
if [ "$BUILD_ONLY" = false ]; then
  ./case.submit
fi

exit 0
