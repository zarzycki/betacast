#!/bin/bash

# Script to automatically create a case, configure, build, and run I compset to spin up
# CLM or ELM for betacast runs

# Turn on error checking
set -e

# Usage:
#./auto-script.sh MODELSYSTEM DOERA5 DATE_YYYYMMDD NMONTHS ANOMYEAR
# CESM
#./auto-script.sh 0 0 20200821 12 2081
# E3SM
#./auto-script.sh 1 0 19960113 12 2081

### User settings
modelSystem=${1}         # 0 = CESM/E3SMv1, 1 = E3SMv2
doERA5=${2}               # Use ERA5 DATM? 0 = yes, 1 = no (internal CRUNCEP)
FORECASTDATE=${3}  #What is the date you want to spin up for (00Z)
NMONTHSSPIN=${4}      # Duration of spinup (somewhere b/w 3-12 months seems reasonable)
BETACAST_ANOMYEAR=${5}

### Check if BETACAST_ANOMYEAR is positive integer -- if yes, add deltas, if no, nothing
if [ $BETACAST_ANOMYEAR -lt 1 ]; then
  addDeltas=1
else
  addDeltas=0
fi

if [[ $HOSTNAME = cheyenne* ]]; then
  ### CLM on Cheyenne
  ### REQUIRED
  CIMEROOT=~/work/cesm2_2_0
  PATHTOCASE=~/I-compsets
  ICASENAME=TEST6-ICLM45-f09
  PROJECT=UPSU0032
  MACHINE=cheyenne
  NNODES=12
  RESOL=f09_g16
  RUNQUEUE=premium
  WALLCLOCK="00:29:00"
  
  ### Only required if doERA5 = 0
  BETACAST=/glade/u/home/zarzycki/betacast
  BETACAST_DATM_FORCING_BASE=/glade/scratch/zarzycki/ERA5-DATM/DATM_FORCING/
  
elif [[ $HOSTNAME = cori* || $HOSTNAME = nid* ]]; then
  ### ELM on Cori
  CIMEROOT=~/clean/E3SM-20210824/
  PATHTOCASE=~/I-compsets
  ICASENAME=RoS-ICLM45-ne0conus30x8-002
  PROJECT=m2637
  MACHINE=cori-knl
  NNODES=12
  RESOL=ne0conus30x8_ne0conus30x8_t12
  RUNQUEUE=regular
  WALLCLOCK="02:57:00"
  
  ### Only required if doERA5 = 0
  BETACAST=/global/homes/c/czarzyck/betacast/
  BETACAST_DATM_FORCING_BASE=/global/homes/c/czarzyck/scratch/DATM_FORCING/
  
else
  echo "Can't figure out hostname, need to add new host and config. Exiting for now..."
  exit
fi

if [ $addDeltas -eq 0 ]; then
  ICASENAME=${ICASENAME}_${BETACAST_ANOMYEAR}
fi

# Derived settings that should be same between all machines
BETACAST_DATMDOMAIN=${BETACAST}/land-spinup/gen_datm/gen-datm/
BETACAST_ANOMALIGN=1920
BETACAST_STREAMBASE=${BETACAST_DATM_FORCING_BASE}/ERA5/
BETACAST_ANOMBASE=${BETACAST_DATM_FORCING_BASE}/ANOM_LENS/

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

### Do not edit below this line!
### Date logic
FORECASTYEAR="${FORECASTDATE:0:4}"
FORECASTYEARM1=$((FORECASTYEAR-1))
FORECASTYEARP1=$((FORECASTYEAR+1))
echo "Trying to get CLM restart for "${FORECASTDATE}
echo "Year plus one equals: "${FORECASTYEARP1}
echo "Year minus one equals: "${FORECASTYEARM1}

### Print diagnostics
STARTDATE=`date -d "${FORECASTDATE} - ${NMONTHSSPIN} months" "+%Y-%m-%d"`
echo "Starting at: "${STARTDATE}
if [ $doERA5 -eq 0 ]; then
  echo "Using ERA5 DATM"
  ERA5STYR=1995
  ERA5ENYR=2020
else
  echo "Using CRUNCEP DATM"
fi

### Configure, build, run land model w/ DATM
if [ $doERA5 -ne 0 ] && (( FORECASTYEAR > 2016 )); then
  echo "No default DATM files beyond 2016"
  echo "Need to either find a different DATM set or spin up when coupled"
  echo "STOP"
  exit 1
elif [ $doERA5 -eq 0 ] && (( FORECASTYEARM1 < ${ERA5STYR} )); then
  echo "No ERA5 DATM files earlier than ${ERA5STYR}"
  echo "Need to either find a different DATM set or spin up when coupled"
  echo "STOP"
  exit 1
elif [ $doERA5 -eq 0 ] && (( FORECASTYEAR > ${ERA5ENYR} )); then
  echo "No ERA5 DATM files later than ${ERA5ENYR}"
  echo "Need to either find a different DATM set or spin up when coupled"
  echo "STOP"
  exit 1
elif [ ${#FORECASTDATE} -ne 8 ]; then
  echo "Incorrect string length for FORECASTDATE, needs to be YYYYMMDD"
  echo "STOP"
  exit 1
fi

### Put a block to check everything here?
echo "--------------------------------------------"
echo "modelSystem: "${modelSystem}
echo "FORECASTDATE: "${FORECASTDATE}
echo "NMONTHSSPIN: "${NMONTHSSPIN}
echo "addDeltas: "${addDeltas}
echo "BETACAST_ANOMYEAR: "${BETACAST_ANOMYEAR}
echo "BETACAST_ANOMALIGN: "${BETACAST_ANOMALIGN}
echo "BETACAST_DATMDOMAIN: "${BETACAST_DATMDOMAIN}
echo "BETACAST_STREAMBASE: "${BETACAST_STREAMBASE}
echo "BETACAST_ANOMBASE: "${BETACAST_ANOMBASE}
echo "CIMEROOT: "${CIMEROOT}
echo "PATHTOCASE: "${PATHTOCASE}
echo "ICASENAME: "${ICASENAME}
echo "PROJECT: "${PROJECT}
echo "MACHINE: "${MACHINE}
echo "NNODES: "${NNODES}
echo "RESOL: "${RESOL}
echo "RUNQUEUE: "${RUNQUEUE}
echo "WALLCLOCK: "${WALLCLOCK}
echo "ICASENAME: "${ICASENAME}
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
./xmlchange STOP_N=2
./xmlchange STOP_OPTION='nyears'
./xmlchange DATM_CLMNCEP_YR_START=${FORECASTYEARM1}
./xmlchange DATM_CLMNCEP_YR_END=${FORECASTYEAR}
./xmlchange DATM_CLMNCEP_YR_ALIGN=${FORECASTYEARM1}
./xmlchange RUN_STARTDATE=${STARTDATE}
./xmlchange STOP_DATE=${FORECASTDATE}
./xmlchange REST_OPTION='end'
./xmlchange DOUT_S=FALSE

### If using ERA5, add the stream files and reset DATM_CLMNCEP_YR_START, etc.
if [ $doERA5 -eq 0 ]; then
  echo "Injecting ERA5 DATM streams"
  cp ${BETACAST}/land-spinup/streams/user_datm.streams.txt.CLMCRUNCEPv7* .
  #REPLACEDIR
  sed -i "s?\${BETACAST_STREAMBASE}?${BETACAST_STREAMBASE}?g" user_datm.streams.txt.CLMCRUNCEPv7.Solar 
  sed -i "s?\${BETACAST_DATMDOMAIN}?${BETACAST_DATMDOMAIN}?g" user_datm.streams.txt.CLMCRUNCEPv7.Solar 
  sed -i "s?\${BETACAST_STREAMBASE}?${BETACAST_STREAMBASE}?g" user_datm.streams.txt.CLMCRUNCEPv7.Precip 
  sed -i "s?\${BETACAST_DATMDOMAIN}?${BETACAST_DATMDOMAIN}?g" user_datm.streams.txt.CLMCRUNCEPv7.Precip 
  sed -i "s?\${BETACAST_STREAMBASE}?${BETACAST_STREAMBASE}?g" user_datm.streams.txt.CLMCRUNCEPv7.TPQW 
  sed -i "s?\${BETACAST_DATMDOMAIN}?${BETACAST_DATMDOMAIN}?g" user_datm.streams.txt.CLMCRUNCEPv7.TPQW  
  ./xmlchange DATM_CLMNCEP_YR_ALIGN=${ERA5STYR}
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
  sed -i "s?\${BETACAST_ANOMBASE}?${BETACAST_ANOMBASE}?g" user_datm.streams.txt.Anomaly.Forcing.Humidity
  sed -i "s?\${BETACAST_DATMDOMAIN}?${BETACAST_DATMDOMAIN}?g" user_datm.streams.txt.Anomaly.Forcing.Humidity 
  sed -i "s?\${BETACAST_ANOMBASE}?${BETACAST_ANOMBASE}?g" user_datm.streams.txt.Anomaly.Forcing.Temperature
  sed -i "s?\${BETACAST_DATMDOMAIN}?${BETACAST_DATMDOMAIN}?g" user_datm.streams.txt.Anomaly.Forcing.Temperature
  sed -i "s?\${BETACAST_ANOMBASE}?${BETACAST_ANOMBASE}?g" user_datm.streams.txt.Anomaly.Forcing.Longwave
  sed -i "s?\${BETACAST_DATMDOMAIN}?${BETACAST_DATMDOMAIN}?g" user_datm.streams.txt.Anomaly.Forcing.Longwave
  sed -i "s?\${BETACAST_ANOMBASE}?${BETACAST_ANOMBASE}?g" user_datm.streams.txt.Anomaly.Forcing.Precip
  sed -i "s?\${BETACAST_DATMDOMAIN}?${BETACAST_DATMDOMAIN}?g" user_datm.streams.txt.Anomaly.Forcing.Precip

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
!fsurdat="/global/homes/c/czarzyck/m2637/betacast/cesmfiles/clm_surfdata_4_5/surfdata_conus_30_x8_simyr2000_c201027.nc"
!finidat=''
!do_transient_pfts = .false.
!check_finidat_fsurdat_consistency = .false.
!fsurdat="/global/cfs/cdirs/e3sm/inputdata/lnd/clm2/surfdata_map/surfdata_ne30np4_simyr2000_c190730.nc"
fsurdat="/global/homes/c/czarzyck/m2637/betacast/cesmfiles/clm_surfdata_4_5/surfdata_conus_30_x8_simyr2000_c201027.nc"
finidat="/global/homes/c/czarzyck/scratch/e3sm_scratch/cori-knl/RoS-F2010C5-ne0conus30x8-001-control/run//landstart//RoS-F2010C5-ne0conus30x8-001-control.elm.r.1996-01-15-00000.nc"
!finidat=''
do_transient_pfts = .false.
check_finidat_fsurdat_consistency = .false.
EOF

### USER! Edit this block if using CLM and need to inject any CLM specific mods (e.g., fsurdat, etc.)
cat > user_nl_clm <<EOF
!finidat=''
EOF

# remove any "incorrect" user nl files depending on model system
if [ $modelSystem -eq 0 ]; then
  rm -v user_nl_elm
elif [ $modelSystem -eq 1 ]; then
  rm -v user_nl_clm
else
  echo "Unknown modeling system set for modelSystem: $modelSystem"
  exit 1
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
./case.submit

exit 0
