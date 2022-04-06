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
############ USER OPTIONS #####################

if [ $# -eq 11 ] 
then
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
else
  echo "Oops... this script requires 11 arguments, not "$#
  exit
fi

TRAJFILE=trajectories.txt.${UQSTR}
OVERWRITE_ATCF=false #CMZ need to build logic that checks if ATCF file has TECH already on it, then overwrite

### Flag needed if using unstructured data (if not using unstruc data, empty string)
CONNECTFLAG="--in_connect ${CONNECTFILE}"

### Path + filelist of data to process
PATHTOFILES=${OUTPUTBASE}/${YYYYMMDDHH}/
FILES=`ls ${PATHTOFILES}/*.?am.${HSTREAMTRACK}.*.nc`

############ TRACKER MECHANICS #####################

# Loop over files to find candidate cyclones
# Find ALL local min in PSL, merge within 3deg
rm cyc.${UQSTR} trajectories.txt.${UQSTR}
touch cyc.${UQSTR}
for f in ${FILES[@]};
do
  echo "Processing $f..."
  ${TEMPESTEXTREMESDIR}/bin/DetectNodes --verbosity 0 --timestride ${TIMESTRIDE} --in_data "${f}" ${CONNECTFLAG} --out cyc_tempest.${UQSTR} --mergedist 5.0 --searchbymin PSL --outputcmd "PSL,min,0;_VECMAG(UBOT,VBOT),max,2"
  cat cyc_tempest.${UQSTR} >> cyc.${UQSTR}
  rm cyc_tempest.${UQSTR}
done

# Stitch candidate cyclones together
# Using range = 3.0 for 3hr, use 5.0 for 6hr.
# Using minlength = 5 = 12 hrs for 3hrly; minlength = 3 = 12 hrs for 6hrly
${TEMPESTEXTREMESDIR}/bin/StitchNodes --format "ncol,lon,lat,slp,wind" --range 5.0 --minlength 5 --maxgap 0 --in cyc.${UQSTR} --out ${TRAJFILE} --threshold "lat,<=,55.0,4;lat,>=,-55.0,4"

sleep 5 #brief delay to nip any small I/O issues in the bud

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
    cat ${ATCFFILE} >> ${ATCFFILEMERGE}  
  fi
  ncl plot_ATCF.ncl
  if [ "$SENDHTML" = true ]; then
    echo "SSHings..."
    ssh colinzar@colinzarzycki.com "mkdir -p /home/colinzar/www/www/current2/${YYYYMMDDHH} "
    scp enstraj.${YYYYMMDDHH}.png colinzar@colinzarzycki.com:/home/colinzar/www/www/current2/${YYYYMMDDHH}/ens_traj.png
  fi
  rm ${ATCFFILE}
fi

# Cleanup
rm cyc.${UQSTR}   #Delete candidate cyclone file
