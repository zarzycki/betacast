#!/bin/bash

# Colin Zarzycki

# Input args
# $1 -- YYYYMMDDHH of cycle (ex: 2018082112)
# $2 -- TCVITFOLDER (ex: ./fin-tcvitals/)
# $3 -- betacast script dir (ex: ~/betacast/)
# ./get-vitals.sh 2024070612 ./fin-tcvitals/ ~/work/betacast/

############ USER OPTIONS #####################

TCVITFOLDER=${2}
YYYYMMDDHH=${1}

echo "get-vitals: Getting TC vitals for ${YYYYMMDDHH} and putting them in ${TCVITFOLDER}"

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
  echo "get-vitals: ${TCVITFILE} doesn't exist, attempting to create from ${TCVITFOLDER}/combined/${COMBINEDVITFILE}"
  mkdir -v -p ${TCVITFOLDER}/combined/

  if [ ! -f ${TCVITFOLDER}/combined/${COMBINEDVITFILE} ]; then
    echo "get-vitals: ${TCVITFOLDER}/combined/${COMBINEDVITFILE} doesn't exist, attempting to download"
    cd ${TCVITFOLDER}/combined/

    # if get_from_nhc, get from NHC servers for this current year by gluing storm archived vitals
    get_from_nhc=true

    if [ "$get_from_nhc" = true ]; then
      BASE_URL="https://ftp.nhc.noaa.gov/atcf/com/"
      TEMP_DIR="./vit_tmp/"
      ORIGINAL_DIR=$(pwd)
      echo "Temporary Directory: $TEMP_DIR"
      echo "Original Directory: $ORIGINAL_DIR"
      mkdir -p $TEMP_DIR ; cd $TEMP_DIR

      # CURL STUFF
      curl -s $BASE_URL > index.html
      # Extract any tcvitals files we can find
      FILES=$(grep -oP 'href="\K[^"]*tcvitals[^"]*' index.html)
      # Loop through each file and download it
      for FILE in $FILES; do
        echo "get-vitals: downloading: $BASE_URL$FILE"
        curl -O "$BASE_URL$FILE"
      done

      # Concatenate all the downloaded files into a single combined file
      head *tcvitals*
      cat *tcvitals* > "$COMBINEDVITFILE"

      # Sort the combined file by the 4th and 5th columns (YYYYMMDD then HHNN)
      sort -k 4,4 -k 5,5 $COMBINEDVITFILE -o $COMBINEDVITFILE

      # Move the combined and sorted file to the original directory
      mv $COMBINEDVITFILE $ORIGINAL_DIR/$COMBINEDVITFILE

      # Clean up the temporary directory
      cd .. ; rm -rf $TEMP_DIR

      # Tell the code to clear this file when done
      CLEARCOMBINEDVIT=true
    else  # get from RAL server
      if [[ $yearstr -gt 2018 ]]; then
        echo "get-vitals: Year $yearstr is 2019 or later, try to get from NCAR RAL"
        wget http://hurricanes.ral.ucar.edu/repository/data/tcvitals_open/${COMBINEDVITFILE}
        CLEARCOMBINEDVIT=true
      else
        echo "get-vitals: Year $yearstr is 2018 or earlier, try to get from NOAA NHC archives"
        wget -q -nd -r --no-parent -P ./ -A dat "https://ftp.nhc.noaa.gov/atcf/archive/${yearstr}/tcvitals/"
        rm *tcvitals.dat
        cat *tcvitals-arch.dat > ${COMBINEDVITFILE}
        rm *-tcvitals-arch.dat
        CLEARCOMBINEDVIT=false
      fi
    fi
    sed -i $'s/\t/  /g' combined_tcvitals.${yearstr}.dat
    cd ../..
  fi
  echo "Last 10 lines of ${TCVITFOLDER}/combined/combined_tcvitals.${yearstr}.dat"
  tail -10 ${TCVITFOLDER}/combined/combined_tcvitals.${yearstr}.dat
  echo "-------"
  grep "${yearstr}${monthstr}${daystr} ${cyclestr}00" ${TCVITFOLDER}/combined/combined_tcvitals.${yearstr}.dat > ${TCVITFILE}
  if [ "$CLEARCOMBINEDVIT" = true ]; then rm -v ${TCVITFOLDER}/combined/combined_tcvitals.${yearstr}.dat ; fi
fi

# Remove duplicate lines if necessary
sort -u ${TCVITFILE} -o ${TCVITFILE}

echo "get-vitals: end of script"
echo "******* TCVIT FILE INFO"
head -20 ${TCVITFILE}
