#!/bin/bash

CASESTR="ERA5-ELMinic-ne128pg2_v2_20020101_000_01"
DIRECTORY="/pscratch/sd/c/czarzyck/e3sm_scratch/pm-cpu/${CASESTR}/run/"
IGNORELAST=1
FINALDIR="/global/cfs/cdirs/m2637/E3SM_SCREAM_files/finidat/ne128pg2_v2/ERA5/"

#####

source /global/homes/c/czarzyck/betacast/utils.sh
module load conda ; conda activate betacast
START_TIME=$(date +%s)
START_TIME_FORMATTED=$(date '+%Y-%m-%d %H:%M:%S')

cdv $DIRECTORY

# Get all elm.r.*nc files, sort them
files=($(ls ${CASESTR}*elm.r.*nc 2>/dev/null | sort))
count=${#files[@]}

echo "Got $count files"

if (( count <= IGNORELAST )); then
  echo "Not enough restart files (found $count), exiting."
  exit 0
fi

# Loop through all but last $IGNORELAST
for ((i=0; i<count-$IGNORELAST; i++)); do
  f=${files[$i]}
  timestring="${f##*.elm.r.}"  # gets 2002-09-23-43200.nc
  timestring="${timestring%.nc}"  # strips off .nc
  echo $timestring
  # Check if file exists first before compressing or removing
  [[ -f "${CASESTR}.elm.r.${timestring}.nc" ]] && compress_file "${CASESTR}.elm.r.${timestring}.nc" zstd
  [[ -f "${CASESTR}.mosart.r.${timestring}.nc" ]] && compress_file "${CASESTR}.mosart.r.${timestring}.nc" zstd
  [[ -f "${CASESTR}.mosart.rh0.${timestring}.nc" ]] && rm -fv "${CASESTR}.mosart.rh0.${timestring}.nc"
  [[ -f "${CASESTR}.elm.rh0.${timestring}.nc" ]] && rm -fv "${CASESTR}.elm.rh0.${timestring}.nc"
done

# Create archival dir if it doesn't exist
mkdir -p $FINALDIR

# compgen -G checks to see if any files exist first, if they do move.
if compgen -G "${CASESTR}.elm.r.*.zst" > /dev/null; then
  mv -v ${CASESTR}.elm.r.*.zst "$FINALDIR"
fi
if compgen -G "${CASESTR}.mosart.r.*.zst" > /dev/null; then
  mv -v ${CASESTR}.mosart.r.*.zst "$FINALDIR"
fi

# Calculate and display completion time and duration
END_TIME=$(date +%s)
END_TIME_FORMATTED=$(date '+%Y-%m-%d %H:%M:%S')
ELAPSED=$((END_TIME - START_TIME))
HOURS=$((ELAPSED / 3600))
MINUTES=$(((ELAPSED % 3600) / 60))
SECONDS=$((ELAPSED % 60))
echo "Job completed at $END_TIME_FORMATTED, time to complete: $HOURS hours, $MINUTES minutes, $SECONDS seconds"
