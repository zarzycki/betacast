#!/bin/bash

# Turn on error checking
set -e

### User settings
FORECASTDATE=20160304
CIMEROOT=~/scream
PATHTOCASE=~/
CASENAME=YYY_ICLM45
PROJECT=m2637
MACHINE=cori-knl
NNODES=8
COMPSET=ICLM45
RESOL=ne30_ne30
RUNQUEUE=premium

### Date logic
FORECASTYEAR="${FORECASTDATE:0:4}"
FORECASTYEARM1=$((FORECASTYEAR-1))
FORECASTYEARP1=$((FORECASTYEAR+1))
echo "Trying to get CLM restart for "${FORECASTDATE}
echo "Year plus one equals: "${FORECASTYEARP1}
echo "Year minus one equals: "${FORECASTYEARM1}

### Configure, build, run land model w/ DATM
if (( FORECASTYEAR > 2016 )); then
  echo "No default DATM files beyond 2016"
  echo "Need to either find a different DATM set or spin up when coupled"
  echo "STOP"
elif [ ${#FORECASTDATE} -ne 8 ]; then
  echo "Incorrect string length for FORECASTDATE, needs to be YYYYMMDD"
  echo "STOP"
else
  cd ${CIMEROOT}/cime/scripts
  ./create_newcase --case ${PATHTOCASE}/${CASENAME} --compset ${COMPSET} --res ${RESOL} --mach ${MACHINE} --project ${PROJECT}
  cd ${PATHTOCASE}/${CASENAME}
  ./xmlchange CHARGE_ACCOUNT=${PROJECT}
  ./xmlchange NTASKS=-${NNODES}
  ./xmlchange NTASKS_ESP=1
  ./xmlchange --force JOB_QUEUE=${RUNQUEUE}
  ./case.setup
  ./xmlchange DATM_MODE=CLMCRUNCEPv7
  ./xmlchange STOP_N=3
  ./xmlchange STOP_OPTION='nyears'
  ./xmlchange DATM_CLMNCEP_YR_START=${FORECASTYEARM1}
  ./xmlchange DATM_CLMNCEP_YR_END=${FORECASTYEARP1}
  ./xmlchange DATM_CLMNCEP_YR_ALIGN=${FORECASTYEARM1}
  ./xmlchange RUN_STARTDATE=${FORECASTYEARM1}-01-01
  ./xmlchange STOP_DATE=${FORECASTDATE}
  ./xmlchange REST_OPTION='end'
  ./case.build
  ./case.submit
fi

exit
