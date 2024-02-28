#!/bin/bash

set -e
source ../utils.sh

### Set vars here for now
datesfile="dates.CHEY.VR28.NATL.REF.CAM5.4CLM5.0.dtime900.txt"
CASEDIR="/global/homes/c/czarzyck/I-runs/ne128pg2-spinup_19840101_0000"
path_to_rundir="/pscratch/sd/c/czarzyck/e3sm_scratch/pm-cpu/ne128pg2-spinup_19840101_0000/run"
CIMEsubstring=""
CIMEbatchargs=""
betacast="/global/homes/c/czarzyck/betacast"
WALLCLOCK="02:00:00"
RUNQUEUE="regular"
# Copy this last restart to a stash folder
DIRSTASH=/global/cfs/cdirs/m2637/E3SM_SCREAM_files/finidat/ne128pg2/CHEY.VR28.NATL.REF.CAM5.4CLM5.0.dtime900/

#--------

cd $betacast/land-spinup/
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
./xmlchange STOP_DATE="$yearstr$monthstr$daystr"
run_CIME2 "$path_to_rundir" "$CIMEsubstring" "$CIMEbatchargs" true

# Assuming run_CIME2 was successful, let's set CONTINUE_RUN to TRUE
# this is because for a cold start we want FALSE, but every other run is TRUE
# So there is no harm in doing this as long as the first run succeeds.
./xmlchange CONTINUE_RUN=TRUE

# Copy restart files to some directory for stashing purposes
cd $path_to_rundir

# Assuming run was succesful, delete all other restarts except the last one
delete_except_last "*.clm2.r.*,*.elm.r.*,*.elm.rh0.*,*.elm.rh1.*"
delete_except_last "*.mosart.r.*,*.mosart.rh0.*"
delete_except_last "*.cpl.r.*"

echo "Done with delete_except_last"

if [ ! -d "$DIRSTASH" ]; then mkdir -v "$DIRSTASH"; fi
if [ ! -d "$DIRSTASH/logs/" ]; then mkdir -v "$DIRSTASH/logs/"; fi

function safe_cp2() {

  local dest="$2"
  # split on multiple input patterns if necessary
  IFS=',' read -r -a patterns <<< "$1"

  # Loop over all matching patterns
  for pattern in "${patterns[@]}"; do

    # Trim leading and trailing spaces from the pattern
    pattern=$(echo "$pattern" | xargs)

    # Collect all files matching the pattern, sorted
    files=($(ls $pattern 2>/dev/null | sort))

    if [ ${#files[@]} -eq 0 ]; then
      echo "Error: ZERO files found matching pattern $pattern."
    else
      cp -v -- "${files[@]}" "$dest"
    fi
  done
}

safe_cp2 "*.clm2.r.*,*.elm.r.*" $DIRSTASH
safe_cp2 "*.mosart.r.*" $DIRSTASH
for file in $DIRSTASH/*.nc; do
  [ -e "$file" ] || continue # Skip if no files match
  compress_file "$file" zstd
done
mv -fv *tar.gz $DIRSTASH/logs/
rm -fv *.bin

# Return to the script dir to do things and resubmit
cd $betacast/land-spinup/
remove_top_line_from_dates ${datesfile}

echo "*-*-*-* Automatically resubbing next date!"
exec ./drive-dates.sh