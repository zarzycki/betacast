#!/bin/bash

# Driver script to search for cyclones in unstructured CAM-SE data
# After finding trajectories, uses model-to-atcf.ncl to collocate
# found storms with tcvitals, writes ATCF format which can be used
# by met-tools for track stats, etc.
#
# Colin Zarzycki

# Input args
# $1 -- YYYYMMDDHH of cycle (ex: 2018082112)
# $2 -- Casename (ex: conus_30_x8_forecast
# $3 -- Path to TCVitals file (input)
# $4 -- Path to ATCF file (output)
# $5 -- connectivity file (full path)
# $6 -- output base (everything before date)
# $7 -- send HTML (usually false)
# $8 -- which h stream are we tracking on? (ex: h1)
# $9 -- timestride, should be 2 for 3-hourly data, 1 for 6-hourly data
# $10 -- ATCFTECH
# $11 -- where is the TE noMPI binary on this system?
# $12 -- OPTIONAL: ARCHIVEDIR
############ USER OPTIONS #####################

if [ $# -eq 11 ] || [ $# -eq 12 ]; then
  archive_files=false
  YYYYMMDDHH=${1}
  UQSTR=${2}
  TCVITFILE=${3}
  ATCFFILE=${4}
  CONNECTFILE=${5}
  OUTPUTBASE=${6}
  SENDHTML=${7}
  HSTREAMTRACK=${8}
  TIMESTRIDE=${9}
  ATCFTECH=${10}
  TEMPESTEXTREMESDIR=${11}
  if [ $# -eq 12 ]; then
    ARCHIVEDIR=${12}
    archive_files=true
    echo "cyclone-tracking-drive: archive requested"
  fi
else
  echo "cyclone-tracking-driver: Oops... this script requires 11 or 12 arguments, not $#"
  exit 1
fi

# Currently, code purges this ATCFTECH's record in YYYYMMDDHH ATCF file.
# It keeps other models/ensembles, to allow for serial catting.
# Setting this to true nukes whole thing each time tracker is invoked
OVERWRITE_ATCF=false

# The TE trajectory file for this particular betacast case
TRAJFILE=trajectories.txt.${UQSTR}

# Call TE tracking
bash run-te-tracking.sh \
    "${CONNECTFILE}" \
    "${OUTPUTBASE}" \
    "${YYYYMMDDHH}" \
    "${UQSTR}" \
    "${HSTREAMTRACK}" \
    "${TEMPESTEXTREMESDIR}" \
    "${TIMESTRIDE}" \
    "${TRAJFILE}"

# Archive files if requested
if [ "$archive_files" = true ]; then
  echo "cyclone-tracking-driver: archiving files"
  mkdir -p "$ARCHIVEDIR/cyclone-tracking"
  cp -v "${TRAJFILE}" "$ARCHIVEDIR/cyclone-tracking"
  cp -v "cyc.${UQSTR}" "$ARCHIVEDIR/cyclone-tracking"
fi

# If a vitals file exists, we want to "match" up the corresponding trajectory from
# tempest with the +0hr vitals location and write an ATCF record
if [ -f ${TCVITFILE} ]; then
  export TRAJFILE="${TRAJFILE}"
  export TCVITFILE="${TCVITFILE}"
  export ATCFFILE="${ATCFFILE}"
  export ATCFTECH="${ATCFTECH}"
  ncl model-to-atcf.ncl
  export YYYYMMDDHH="${ATCFFILE:(-10)}"
  export ATCFFILEMERGE="./fin-atcf/atcf.tempest.${YYYYMMDDHH}"
  if [ "$OVERWRITE_ATCF" = true ]; then
    cp ${ATCFFILE} ${ATCFFILEMERGE}
  else
    # First check to see if the master file exists
    # Delete all matching lines from this model the master ATCF if they exist!
    if [[ -e "${ATCFFILEMERGE}" ]]; then
      sed -i "/${ATCFTECH}/d" ${ATCFFILEMERGE}
    fi
    # Now concat new lines!
    cat ${ATCFFILE} >> ${ATCFFILEMERGE}
  fi
  echo "cyclone-tracking-driver: plotting ATCF"
  ncl plot_ATCF.ncl
  if [ "$SENDHTML" = true ]; then
    echo "cyclone-tracking-driver: SSHings..."
    ssh colinzar@colinzarzycki.com "mkdir -p /home/colinzar/www/www/current2/${YYYYMMDDHH} "
    scp enstraj.${YYYYMMDDHH}.png colinzar@colinzarzycki.com:/home/colinzar/www/www/current2/${YYYYMMDDHH}/ens_traj.png
  fi
  if [ "$archive_files" = true ]; then
    mkdir -p "$ARCHIVEDIR/cyclone-tracking"
    cp -v "${ATCFFILE}" "$ARCHIVEDIR/cyclone-tracking"
  fi
  rm -v ${ATCFFILE}
fi

# Cleanup
rm cyc.${UQSTR}   #Delete candidate cyclone file
