#!/bin/bash

# Script to automatically create a case, configure, build, and run I compset to spin up
# CLM or ELM for betacast runs

# Turn on error checking
set -e

### User settings
FORECASTDATE=19960113
NMONTHSSPIN=12
CIMEROOT=~/E3SM-dev
PATHTOCASE=~/I-compsets
CASENAME=RoS-ICLM45-ne30
PROJECT=m2637
MACHINE=cori-knl
NNODES=12
COMPSET=IELM
RESOL=ne30_ne30
RUNQUEUE=premium

### Do not edit below this line!

### Date logic
FORECASTYEAR="${FORECASTDATE:0:4}"
FORECASTYEARM1=$((FORECASTYEAR-1))
FORECASTYEARP1=$((FORECASTYEAR+1))
echo "Trying to get CLM restart for "${FORECASTDATE}
echo "Year plus one equals: "${FORECASTYEARP1}
echo "Year minus one equals: "${FORECASTYEARM1}

STARTDATE=`date -d "${FORECASTDATE} - ${NMONTHSSPIN} months" "+%Y-%m-%d"`

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
  ./xmlchange NTASKS=-${NNODES}
  ./xmlchange NTASKS_ESP=1
  ./xmlchange DATM_MODE=CLMCRUNCEPv7
  ./xmlchange STOP_N=2
  ./xmlchange STOP_OPTION='nyears'
  ./xmlchange DATM_CLMNCEP_YR_START=${FORECASTYEARM1}
  ./xmlchange DATM_CLMNCEP_YR_END=${FORECASTYEAR}
  ./xmlchange DATM_CLMNCEP_YR_ALIGN=${FORECASTYEARM1}
  ./xmlchange RUN_STARTDATE=${STARTDATE}
  ./xmlchange STOP_DATE=${FORECASTDATE}
  ./xmlchange REST_OPTION='end'
  
  cat > user_nl_elm <<EOF
  !fsurdat="/global/homes/c/czarzyck/m2637/betacast/cesmfiles/clm_surfdata_4_5/surfdata_conus_30_x8_simyr2000_c201027.nc"
  !finidat=''
  !do_transient_pfts = .false.
  !check_finidat_fsurdat_consistency = .false.
  fsurdat="/global/cfs/cdirs/e3sm/inputdata/lnd/clm2/surfdata_map/surfdata_ne30np4_simyr2000_c190730.nc"
  finidat=''
  do_transient_pfts = .false.
  check_finidat_fsurdat_consistency = .false.
EOF

  ./case.setup
  ./case.build
  ./xmlchange JOB_WALLCLOCK_TIME="06:15:00"
  ./xmlchange CHARGE_ACCOUNT=${PROJECT}
  ./xmlchange --force JOB_QUEUE=${RUNQUEUE}
  ./case.submit
fi

exit
