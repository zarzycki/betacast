#!/bin/bash

# Usage:
#   run_te_tracking.sh CONNECTFILE OUTPUTBASE YYYYMMDDHH UQSTR HSTREAMTRACK TEMPESTEXTREMESDIR TIMESTRIDE TRAJFILE

set -e

if [ "$#" -ne 8 ]; then
  echo "Usage: $0 CONNECTFILE OUTPUTBASE YYYYMMDDHH UQSTR HSTREAMTRACK TEMPESTEXTREMESDIR TIMESTRIDE TRAJFILE"
  exit 1
fi

CONNECTFILE=$1
OUTPUTBASE=$2
YYYYMMDDHH=$3
UQSTR=$4
HSTREAMTRACK=$5
TEMPESTEXTREMESDIR=$6
TIMESTRIDE=$7

# The TE trajectory file for this particular betacast case
TRAJFILE=$8

### Flag needed if using unstructured data (if not using unstruc data, empty string)
CONNECTFLAG="--in_connect ${CONNECTFILE}"

### Path + filelist of data to process
PATHTOFILES=${OUTPUTBASE}/${YYYYMMDDHH}/
FILES=$(ls ${PATHTOFILES}/*.?am.${HSTREAMTRACK}.*.nc)

DN_MERGEDIST=5.0
SN_RANGE=5.0
SN_MINLENGTH=5
SN_MAXGAP=0
SN_LATMAX=55.0
SN_LATMAX_N=4

[ -f cyc.${UQSTR} ] && rm -v cyc.${UQSTR}
[ -f "${TRAJFILE}" ] && rm -v "${TRAJFILE}"
touch cyc.${UQSTR}

# Loop over files to find candidate cyclones
# Find ALL local min in PSL, merge within DN_MERGEDIST
for f in ${FILES[@]}; do
  echo "run_te_tracking.sh: Processing $f..."
  ${TEMPESTEXTREMESDIR}/bin/DetectNodes \
    --verbosity 0 \
    --timestride ${TIMESTRIDE} \
    --in_data "${f}" \
    ${CONNECTFLAG} \
    --out cyc_tempest.${UQSTR} \
    --mergedist ${DN_MERGEDIST} \
    --searchbymin PSL \
    --outputcmd "PSL,min,0;_VECMAG(UBOT,VBOT),max,2"

  cat cyc_tempest.${UQSTR} >> cyc.${UQSTR}
  rm -v cyc_tempest.${UQSTR}
done

# Stitch candidate cyclones together
# Using range = 3.0 for 3hr, use 5.0 for 6hr.
# Using minlength = 5 = 12 hrs for 3hrly; minlength = 3 = 12 hrs for 6hrly
${TEMPESTEXTREMESDIR}/bin/StitchNodes \
  --format "ncol,lon,lat,slp,wind" \
  --range ${SN_RANGE} \
  --minlength ${SN_MINLENGTH} \
  --maxgap ${SN_MAXGAP} \
  --in cyc.${UQSTR} \
  --out "${TRAJFILE}" \
  --threshold "lat,<=,${SN_LATMAX},${SN_LATMAX_N};lat,>=,-${SN_LATMAX},${SN_LATMAX_N}"

sleep 5 #brief delay to nip any small I/O issues in the bud