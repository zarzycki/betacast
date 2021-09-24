#!/bin/bash

# Colin Zarzycki

# Input args
# $1 -- YYYYMMDDHH of cycle (ex: 2018082112)
# $2 -- TCVITFOLDER (ex: ./fin-tcvitals/)
# $3 -- betacast script dir (ex: ~/betacast/)

############ USER OPTIONS #####################

TCVITFOLDER=${2}
YYYYMMDDHH=${1}

echo "Getting TC vitals for ${YYYYMMDDHH} and putting them in ${TCVITFOLDER}"

# Split YYYYMMDDHH into components
yearstr=${YYYYMMDDHH:0:4}
monthstr=${YYYYMMDDHH:4:2}
daystr=${YYYYMMDDHH:6:2}
cyclestr=${YYYYMMDDHH:8:2}

TCVITFILE=${TCVITFOLDER}/tcvitals.${yearstr}${monthstr}${daystr}${cyclestr}
mkdir -p ${TCVITFOLDER}
mkdir -p ./fin-figs/
mkdir -p ./fin-atcf/
if [ ! -f ${TCVITFILE} ]; then   #if TCVITFILE doesn't exist, download
  COMBINEDVITFILE=combined_tcvitals.${yearstr}.dat
  echo "${TCVITFILE} doesn't exist, attempting to create from ${TCVITFOLDER}/combined/${COMBINEDVITFILE}"
  mkdir -v -p ${TCVITFOLDER}/combined/
  if [ ! -f ${TCVITFOLDER}/combined/${COMBINEDVITFILE} ]; then
    echo "${TCVITFOLDER}/combined/${COMBINEDVITFILE} doesn't exist, attempting to download"
    cd ${TCVITFOLDER}/combined/
    if [[ $yearstr -gt 2018 ]]; then
      echo "Year $yearstr is 2019 or later, try to get from NCAR RAL"
      wget http://hurricanes.ral.ucar.edu/repository/data/tcvitals_open/${COMBINEDVITFILE}
      CLEARCOMBINEDVIT=true
    else
      echo "Year $yearstr is 2018 or earlier, try to get from NOAA NHC archives"
      wget -q -nd -r --no-parent -P ./ -A dat "https://ftp.nhc.noaa.gov/atcf/archive/${yearstr}/tcvitals/"
      rm *tcvitals.dat
      cat *tcvitals-arch.dat > ${COMBINEDVITFILE}
      rm *-tcvitals-arch.dat
      CLEARCOMBINEDVIT=false
    fi
    sed -i $'s/\t/  /g' combined_tcvitals.${yearstr}.dat
    cd ../..
  fi
  grep "${yearstr}${monthstr}${daystr} ${cyclestr}00" ${TCVITFOLDER}/combined/combined_tcvitals.${yearstr}.dat > ${TCVITFILE}
  if [ "$CLEARCOMBINEDVIT" = true ]; then rm -v ${TCVITFOLDER}/combined/combined_tcvitals.${yearstr}.dat ; fi
fi

echo "******* TCVIT FILE INFO"
head -20 ${TCVITFILE}