#!/bin/bash

set -e

NAMELISTFILE=$1

source ../../utils.sh

# Read the namelist
read_bash_nl "${NAMELISTFILE}"

# Build relevant directories
CASEDIR=$CASESRC/$CASENAME
PATH_TO_RUNDIR=$BASERUN/$CASENAME/run/
SCRIPTDIR=$BETACAST/land-spinup/drive-dates/

echo $CASEDIR
echo $PATH_TO_RUNDIR
echo $SCRIPTDIR

## Set defaults
## Set these in your namelistfile if need to
if [ -z "${CIMEsubstring+x}" ]; then CIMEsubstring=""; fi
if [ -z "${CIMEbatchargs+x}" ]; then CIMEbatchargs=""; fi

## Get to this directory
cd $SCRIPTDIR

# Parse dates file
echo "Using dates in: "${datesfile}
longdate=$(get_top_line_from_dates "${datesfile}")
echo "Getting parsed time from $longdate"
parse_YYYYMMDDHH $longdate
echo "From datesfile, read in: "$yearstr' '$monthstr' '$daystr' '$cyclestr'Z'

### Go to the case dir and do things
cd $CASEDIR
./xmlchange JOB_WALLCLOCK_TIME=${WALLCLOCK}
./xmlchange --force JOB_QUEUE=${RUNQUEUE}
./xmlchange REST_OPTION="end"
./xmlchange STOP_N="86400"
./xmlchange STOP_OPTION="date"
./xmlchange STOP_DATE="$yearstr$monthstr$daystr"
run_CIME2 "$PATH_TO_RUNDIR" "$CIMEsubstring" "$CIMEbatchargs" true

# Copy restart files to some directory for stashing purposes
cd $PATH_TO_RUNDIR

# Assuming run was succesful, delete all other restarts except the last one
delete_except_last "*.clm2.r.*,*.elm.r.*,*.elm.rh0.*,*.elm.rh1.*"
delete_except_last "*.mosart.r.*,*.mosart.rh0.*"
delete_except_last "*.cpl.r.*"

echo "Done with delete_except_last"

## Make stash directories if they don't exist
if [ ! -d "$DIRSTASH" ]; then mkdir -pv "$DIRSTASH"; fi
if [ ! -d "$DIRSTASH/logs/" ]; then mkdir -pv "$DIRSTASH/logs/"; fi

## Copy relevant files over to dirstash, compress any uncompressed, and cleanup
safe_cp2 "*.clm2.r.*,*.elm.r.*" $DIRSTASH
safe_cp2 "*.mosart.r.*" $DIRSTASH
for file in $DIRSTASH/*.nc; do
  [ -e "$file" ] || continue # Skip if no files match
  compress_file "$file" zstd
done
mv -fv *.log.*.gz $DIRSTASH/logs/
rm -fv *.bin

cd $CASEDIR
# Assuming run_CIME2 was successful, let's set CONTINUE_RUN to TRUE
# this is because for a cold start we want FALSE, but every other run is TRUE
# So there is no harm in doing this as long as the first run succeeds.
./xmlchange CONTINUE_RUN=TRUE

# Return to the script dir to do things and resubmit
cd $SCRIPTDIR
remove_top_line_from_dates ${datesfile}

echo "*-*-*-* Automatically resubbing next date!"
exec ./drive-dates.sh "$NAMELISTFILE"
