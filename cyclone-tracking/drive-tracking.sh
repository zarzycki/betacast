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
############ USER OPTIONS #####################

## Path to serial TempestExtremes binaries
TEMPESTEXTREMESDIR=/storage/home/cmz5202/sw/tempestextremes_noMPI/

TCVITFILE=${3}
UQSTR=${2}
ATCFFILE=${4}
TRAJFILE=trajectories.txt.${UQSTR}
if [ ${UQSTR:(-3)} == "001" ]; then
  ATCFTECH="C501"
elif [ ${UQSTR:(-3)} == "002" ]; then
  ATCFTECH="C502"
elif [ ${UQSTR:(-3)} == "003" ]; then
  ATCFTECH="C503"
else
  ATCFTECH="CAM5"
fi

### Flag needed if using unstructured data (if not using unstruc data, empty string)
CONNECTFLAG="--in_connect ./nhemitc30x4.connect.dat"

### Path + filelist of data to process
PATHTOFILES=/storage/home/cmz5202/scratch/output/${2}/run/${1}/
FILES=`ls ${PATHTOFILES}/*.cam.h0.*.nc`

############ TRACKER MECHANICS #####################

# Loop over files to find candidate cyclones
# Find ALL local min in PSL, merge within 3deg
rm cyc.${UQSTR} trajectories.txt.${UQSTR}
touch cyc.${UQSTR}
for f in ${FILES[@]};
do
  echo "Processing $f..."
  ${TEMPESTEXTREMESDIR}/bin/DetectCyclonesUnstructured --verbosity 0 --timestride 2 --in_data "${f}" ${CONNECTFLAG} --out cyc_tempest.${UQSTR} --mergedist 3.0 --searchbymin PSL --outputcmd "PSL,min,0;_VECMAG(UBOT,VBOT),max,2"
  cat cyc_tempest.${UQSTR} >> cyc.${UQSTR}
  rm cyc_tempest.${UQSTR}
done

# Stitch candidate cyclones together
# Using range = 3.0 for 3hr, use 5.0 for 6hr.
# Using minlength = 5 = 12 hrs for 3hrly; minlength = 3 = 12 hrs for 6hrly
${TEMPESTEXTREMESDIR}/bin/StitchNodes --format "ncol,lon,lat,slp,wind" --range 3.0 --minlength 5 --maxgap 0 --in cyc.${UQSTR} --out ${TRAJFILE} --threshold "lat,<=,55.0,4;lat,>=,-55.0,4"

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
  cat ${ATCFFILE} >> ${ATCFFILEMERGE}
  ncl plot_ATCF.ncl
  echo "SSHings..."
  ssh colinzar@colinzarzycki.com "mkdir -p /home/colinzar/www/www/current2/${YYYYMMDDHH} "
  scp enstraj.${YYYYMMDDHH}.png colinzar@colinzarzycki.com:/home/colinzar/www/www/current2/${YYYYMMDDHH}/ens_traj.png
  rm ${ATCFFILE}
fi

# Cleanup
rm cyc.${UQSTR}   #Delete candidate cyclone file
