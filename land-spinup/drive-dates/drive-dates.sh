#!/bin/bash

echo "Command used: $0 \"$@\""
echo "PID       : $$"
echo "Host      : $(hostname)"
echo "User      : $(whoami)"
echo "Start time: $(date -u +"%Y-%m-%d %H:%M:%S UTC")"
echo "Shell     : $SHELL"

set -e
# Print the failing command and line number before exiting on error
trap 'echo "Error on line ${LINENO}: ${BASH_COMMAND}"; exit 1' ERR

NAMELISTFILE=$1

source ../../utils.sh

# Read the namelist
read_bash_nl "${NAMELISTFILE}"

# Build relevant directories
CASEDIR=$CASESRC/$CASENAME
PATH_TO_RUNDIR=$BASERUN/$CASENAME/run/
SCRIPTDIR=$BETACAST/land-spinup/drive-dates/

echo "CASEDIR: $CASEDIR"
echo "PATH_TO_RUNDIR: $PATH_TO_RUNDIR"
echo "SCRIPTDIR: $SCRIPTDIR"

## Set defaults
## Set these in your namelistfile if need to
if [ -z "${CIMEsubstring+x}" ]; then CIMEsubstring=""; fi
if [ -z "${CIMEbatchargs+x}" ]; then CIMEbatchargs=""; fi
if [ -z "${CIMEMAXTRIES+x}" ]; then CIMEMAXTRIES=3; fi
if [ -z "${SHORTCLOCK+x}" ]; then SHORTCLOCK=$WALLCLOCK; fi
if [ -z "${SHORTQUEUE+x}" ]; then SHORTQUEUE=$RUNQUEUE; fi
if [ -z "${SHORTCUTOFFHRS+x}" ]; then SHORTCUTOFFHRS=120; fi

## Get to this directory
cd $SCRIPTDIR

## ---------------------------------------------------------------------------------------

# Parse dates file
echo "Using dates in: "${datesfile}
longdate=$(get_top_line_from_dates "${datesfile}")
echo "Getting parsed time from $longdate"
parse_YYYYMMDDHH $longdate
echo "From datesfile, read in: "$yearstr' '$monthstr' '$daystr' '$cyclestr'Z'

### New logic for figure out how many **hours** we have to march forward before restart

echo "Target date: ${yearstr}-${monthstr}-${daystr} ${cyclestr}:00"

# Need to go to casedir to see what's up
cd $CASEDIR

# Get DRV_RESTART_POINTER and trim
output=$(./xmlquery DRV_RESTART_POINTER | xargs)
echo "Raw xmlquery output: $output"

restart_pointer="${output##*.}"
echo "Initial restart_pointer (from output): $restart_pointer"

# Check if restart_pointer is YYYY-MM-DD-SSSSS, if not try rpointer
if [[ ! "$restart_pointer" =~ ^[0-9]{4}-[0-9]{2}-[0-9]{2}-[0-9]{5}$ ]]; then

  # Check for rpointer.lnd file
  if [[ -f "$PATH_TO_RUNDIR/rpointer.lnd" ]]; then
    restart_pointer=$(sed -E 's/.*([0-9]{4}-[0-9]{2}-[0-9]{2}-[0-9]{5}).*/\1/' "$PATH_TO_RUNDIR/rpointer.lnd")
    echo "Second restart_pointer (from rpointer): $restart_pointer"
  else
    echo "No rpointer file found at: $PATH_TO_RUNDIR/rpointer.lnd, continuing..."
  fi

  # Check if we still don't have YYYY-MM-DD-SSSSS, if not use STARTDATE
  if [[ ! "$restart_pointer" =~ ^[0-9]{4}-[0-9]{2}-[0-9]{2}-[0-9]{5}$ ]]; then
    echo "Restart pointer format invalid. Falling back to RUN_STARTDATE + START_TOD..."

    run_startdate=$(./xmlquery RUN_STARTDATE | awk -F': ' '{print $2}' | xargs)
    echo "RUN_STARTDATE: $run_startdate"

    start_tod=$(./xmlquery START_TOD | awk -F': ' '{print $2}' | xargs)
    echo "START_TOD: $start_tod"

    start_tod_padded=$(printf "%05d" "$start_tod")
    echo "START_TOD padded: $start_tod_padded"

    restart_pointer="${run_startdate}-${start_tod_padded}"
  fi
fi

echo "Final restart pointer: $restart_pointer"

# Function: convert YYYY-MM-DD to days since 0000-01-01 (no leap years)
days_since_epoch() {
  local y=$1 m=$2 d=$3
  local -a month_days=(0 31 59 90 120 151 181 212 243 273 304 334)
  local total_days=$(( y * 365 + month_days[m - 1] + (d - 1) ))
  echo "$total_days"
}

# Parse restart_pointer
restart_date=${restart_pointer%-*}
restart_secs=${restart_pointer##*-}
restart_year=${restart_date:0:4}
restart_month=${restart_date:5:2}
restart_day=${restart_date:8:2}

echo "Parsed restart date: $restart_year-$restart_month-$restart_day"
echo "Restart seconds since midnight: $restart_secs"

# Strip leading zeros
restart_month=$((10#$restart_month))
restart_day=$((10#$restart_day))
monthstr=$((10#$monthstr))
daystr=$((10#$daystr))

echo "Parsed restart date (int): $restart_year-$restart_month-$restart_day"
echo "Parsed future date (int): $yearstr-$monthstr-$daystr"

# Compute total hours since epoch
restart_days=$(days_since_epoch "$restart_year" "$restart_month" "$restart_day")
echo "Restart days since epoch: $restart_days"

restart_hours=$(( restart_days * 24 + restart_secs / 3600 ))
echo "Restart total hours: $restart_hours"

future_days=$(days_since_epoch "$yearstr" "$monthstr" "$daystr")
echo "Future days since epoch: $future_days"

future_hours=$(( future_days * 24 + cyclestr ))
echo "Future total hours: $future_hours"

# Final difference how many hours to go from restart -> target
diff_hours=$(( future_hours - restart_hours ))
echo "Hour difference: $diff_hours"

## ---------------------------------------------------------------------------------------

if (( diff_hours != 0 )); then

  cd $CASEDIR
  if (( diff_hours <= SHORTCUTOFFHRS )); then
    # If we only have a few days to go, allow user to specify debug or other tiny queue
    echo "Short run: $diff_hours hours (<= $SHORTCUTOFFHRS)"
    xmlchange_verbose "JOB_WALLCLOCK_TIME" "$SHORTCLOCK"
    xmlchange_verbose "JOB_QUEUE" "$SHORTQUEUE" "--force"
  else
    # If longer than a few days, go to a regular queue
    echo "Long run: $diff_hours hours (> $SHORTCUTOFFHRS)"
    xmlchange_verbose "JOB_WALLCLOCK_TIME" "$WALLCLOCK"
    xmlchange_verbose "JOB_QUEUE" "$RUNQUEUE" "--force"
  fi

  if [[ -n "${RUNPRIORITY:-}" ]]; then
    xmlchange_verbose "JOB_PRIORITY" "$RUNPRIORITY" "--force"
  fi
  xmlchange_verbose "REST_OPTION" "end"
  #./xmlchange STOP_N="86400"
  #./xmlchange STOP_OPTION="date"
  #./xmlchange STOP_DATE="$yearstr$monthstr$daystr"

  xmlchange_verbose "STOP_N" "$diff_hours"
  xmlchange_verbose "STOP_OPTION" "nhours"
  xmlchange_verbose "STOP_DATE" "-99999"

  ## ---------------------------------------------------------------------------------------
  # Run the model

  # New formulation
  CIMESTATUS=1; CIMEITER=0
  while [ $CIMESTATUS != 0 ] ; do
    if [ "$CIMEITER" -ge "$CIMEMAXTRIES" ]; then
      echo "Exceeded the max tries: $CIMEMAXTRIES ... exiting"
      exit 1
    fi
    CIMEITER=$((CIMEITER+1))
    echo "CIME try $CIMEITER of $CIMEMAXTRIES"
    set +e
    run_CIME2 "$PATH_TO_RUNDIR" "$CIMEsubstring" "$CIMEbatchargs" false
    CIMESTATUS=$?
    set -e
  done
  echo "Returned status $CIMESTATUS"

else
  echo "WARNING: Hour difference is zero. Restart time equals target time."
  echo "  Restart: $restart_year-$restart_month-$restart_day ${restart_secs} seconds"
  echo "  Target: $yearstr-$monthstr-$daystr ${cyclestr}:00"
fi

## ---------------------------------------------------------------------------------------

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

## Move log files to DIRSTASH/logs
shopt -s nullglob
logfiles=(*.log.*.gz)
if [ ${#logfiles[@]} -gt 0 ]; then
  mv -fv "${logfiles[@]}" "$DIRSTASH/logs/"
else
  echo "WARNING: No log files matching *.log.*.gz found to move." >&2
fi

## Cleanup
rm -fv *.bin
rm -fv mpibind.*.log

cd $CASEDIR
# Assuming run_CIME2 was successful, let's set CONTINUE_RUN to TRUE
# this is because for a cold start we want FALSE, but every other run is TRUE
# So there is no harm in doing this as long as the first run succeeds.
xmlchange_verbose "CONTINUE_RUN" "TRUE"

# Return to the script dir to do things and resubmit
cd $SCRIPTDIR
remove_top_line_from_dates ${datesfile}

echo "*-*-*-* Automatically resubbing next date!"
exec ./drive-dates.sh "$NAMELISTFILE"
