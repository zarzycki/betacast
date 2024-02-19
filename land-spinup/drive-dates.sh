#!/bin/bash

set -e
source ../utils.sh

### Set vars here for now
datesfile="test-dates.txt"
CASEDIR="/glade/u/home/zarzycki/I-compsets/INC-spinup_19900101_0000"
path_to_rundir="/glade/derecho/scratch/zarzycki/INC-spinup_19900101_0000/run"
CIMEsubstring=""
CIMEbatchargs=""
betacast="/glade/u/home/zarzycki/betacast"

cd $betacast/land-spinup/
# Parse dates file
echo "Using dates in: "${datesfile}
longdate=$(get_top_line_from_dates "${datesfile}")
echo "Getting parsed time from $longdate"
parse_YYYYMMDDHH $longdate
echo "From datesfile, read in: "$yearstr' '$monthstr' '$daystr' '$cyclestr'Z'

### Go to the case dir and do things
cd $CASEDIR
./xmlchange REST_OPTION="end"
./xmlchange STOP_DATE="$yearstr$monthstr$daystr"
run_CIME2 "$path_to_rundir" "$CIMEsubstring" "$CIMEbatchargs" true

# Assuming run_CIME2 was successful, let's set CONTINUE_RUN to TRUE
# this is because for a cold start we want FALSE, but every other run is TRUE
# So there is no harm in doing this as long as the first run succeeds.
./xmlchange CONTINUE_RUN=TRUE

# Copy restart files to some directory for stashing purposes
cd $path_to_rundir

# Assuming run was succesful, delete all other restarts except the last one
delete_except_last "*.clm2.r.*"
delete_except_last "*.mosart.r.*"

# Copy this last restart to a stash folder
DIRSTASH=$path_to_rundir/landstart/

if [ ! -d "$DIRSTASH" ]; then mkdir -v "$DIRSTASH"; fi
if [ ! -d "$DIRSTASH/logs/" ]; then mkdir -v "$DIRSTASH/logs/"; fi

cp -v *.clm2.r.* $DIRSTASH
cp -v *.mosart.r.* $DIRSTASH
mv -v *log*.gz $DIRSTASH/logs/

# Return to the script dir to do things and resubmit
cd $betacast/land-spinup/
remove_top_line_from_dates ${datesfile}

echo "*-*-*-* Automatically resubbing next date!"
exec ./drive-dates.sh