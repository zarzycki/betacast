#!/bin/bash

### -----------------------------------------------------------------------------------

# This function allows you to create "checkpoints" in the driver script to see
# what vars are being defined globally in each block of code.
# Usage: dump_vars 1
# Creates vars001.txt
function dump_vars {
  local N=$(printf "%03d" $1)
  set | awk -F= '/^\w/ {print $1}' > "vars${N}.txt"
}

### -----------------------------------------------------------------------------------

check_python_dependency() {
  local PKG_NAME=$1
  if python -c "import pkgutil; exit(not pkgutil.find_loader('$PKG_NAME'))"; then
    echo "CHECK_PYTHON_DEPS: $PKG_NAME found"
  else
    echo "CHECK_PYTHON_DEPS: $PKG_NAME NOT found"
    echo "$> conda install $PKG_NAME"
    echo "or otherwise... exiting..."
    exit 24
  fi
}

### -----------------------------------------------------------------------------------

replace_betacast_string() {
  local original_var=$1
  local replacement=$2
  local placeholder=$3
  local NEWVAR

  # Check if any of the inputs are empty
  if [ -z "$original_var" ] || [ -z "$replacement" ] || [ -z "$placeholder" ]; then
    echo "Error: One or more input strings are empty." >&2
    exit 1
  fi

  # Swap placeholder with replacement in original_var
  NEWVAR=${original_var/${placeholder}/${replacement}}

  echo "NAMELIST UPDATED: $original_var --> $NEWVAR" >&2

  # Return the new var, which is captured by the code
  echo "$NEWVAR"
}

### -----------------------------------------------------------------------------------

check_required_vars() {
  local missing_vars=()
  for var_name in "$@"; do
    if [[ -z ${!var_name+x} ]]; then
      missing_vars+=("$var_name")
    fi
  done

  if [[ ${#missing_vars[@]} -ne 0 ]]; then
    echo "The following required variables are not set:" >&2
    printf ' - %s\n' "${missing_vars[@]}" >&2
    exit 1
  else
    echo "All required variables are set!"
    return 0
  fi
}

### -----------------------------------------------------------------------------------

# Strip surrounding quotes from string [$1: variable name]
function strip_quotes() {
  local -n var="$1"
  [[ "${var}" == \"*\" || "${var}" == \'*\' ]] && var="${var:1:-1}"
}

### -----------------------------------------------------------------------------------

# -----------------------------------------------------------------------------
# check_shell_flags
#
# Displays the current status of key Bash execution flags
# Useful for debugging script, just drop "check_shell_flags" in various locales
#
# Function combines the output of two helpers:
#   - check_direct_flags: Inspects shell flags stored in the special $- variable
#       - errexit (-e): Exit immediately on command failure
#       - nounset (-u): Treat unset variables as an error
#       - xtrace  (-x): Print each command before execution
#       - noexec  (-n): Parse commands but do not execute them
#   - check_pipefail: Uses `set -o` to determine the state of the `pipefail` option
# -----------------------------------------------------------------------------

# Helper function #1
check_direct_flags() {
  echo "Current shell flags in \$-: '$-'"
  [[ $- == *e* ]] && echo "  errexit (-e)  : ON" || echo "  errexit (-e)  : OFF"
  [[ $- == *u* ]] && echo "  nounset (-u)  : ON" || echo "  nounset (-u)  : OFF"
  [[ $- == *x* ]] && echo "  xtrace  (-x)  : ON" || echo "  xtrace  (-x)  : OFF"
  [[ $- == *n* ]] && echo "  noexec  (-n)  : ON" || echo "  noexec  (-n)  : OFF"
}


# Helper function #2
check_pipefail() {
  pipefail_status=$(set -o | grep pipefail | awk '{print $2}')
  echo "  pipefail      : $pipefail_status"
}


# Use this function as the one invoked
check_shell_flags() {
  echo "Bash flag status:"
  check_direct_flags
  check_pipefail
}

### -----------------------------------------------------------------------------------

# Check if a boolean is set via 0/1 and correct.
# Usage: check_bool "do_tracking" $do_tracking
check_bool() {

  local var_string=$1
  local var_totest=$2
  local cborig=$2
  local cbnew

  # Convert to lowercase
  var_totest=$(echo ${var_totest,,})

  if [ "${var_totest}" == "0" ] ; then
    #echo "$var_totest"
    echo "* WARNING: Setting $var_string from 0 (as set in namelist) to false."
    echo "WARNING: This will be deprecated in the future!"
    echo "WARNING: To fix, update your namelist for $var_string from 0 to false."
    cbnew=false
  elif [ "${var_totest}" == "1" ] ; then
    echo "* WARNING: Setting $var_string from 1 (as set in namelist) to true."
    echo "WARNING: This will be deprecated in the future!"
    echo "WARNING: To fix, update your namelist for $var_string from 1 to true."
    cbnew=true
  elif [ "${var_totest}" == "f" ] ; then
    echo "* WARNING: Setting $var_string from f to false."
    echo "WARNING: To fix, update your namelist for $var_string from f to false."
    cbnew=false
  elif [ "${var_totest}" == "t" ] ; then
    echo "* WARNING: Setting $var_string from t to true."
    echo "WARNING: To fix, update your namelist for $var_string from t to true."
    cbnew=true
  else
    if [ "$var_totest" != "false" ] && [ "$var_totest" != "true" ] ; then
      echo "ERROR: $var_string is set to $var_totest and not a valid true/false or 0/1, exiting."
      exit
    fi
    # Set cbnew to our lowercased var from earlier since we didn't have to fix.
    cbnew=$var_totest
  fi
  export ${var_string}=${cbnew}
  echo "CHECK_BOOL: ${var_string}     in: $cborig   out: $cbnew"
}
### -----------------------------------------------------------------------------------

function read_bash_nl() {

  local FILETOREAD=$1
  # Sanitize namelist files (add carriage return to end)
  sanitize_file "${FILETOREAD}"

  # Note, ___ will be converted to a space. Namelists cannot have whitespace due to
  # parsing on whitespaces...
  echo "Reading namelist ${FILETOREAD}..."
  local inputstream
  inputstream=$(grep -v "^#" "${FILETOREAD}")
  inputstream="${inputstream//=/ = }"
  #echo $inputstream
  set -- $inputstream
  while [ $# -ge 3 ]; do
    local var="$1"
    local eq="$2"
    local val="${3//___/}"

    if [ "$eq" != "=" ]; then
      echo "NAMELIST ERR: Uh oh, in $FILETOREAD bad format: $1 $2 $3"
      exit 1
    fi

    if [[ "$var" =~ ^[a-zA-Z_][a-zA-Z0-9_]*$ ]]; then
      # Expand references (e.g., $INDEX) in value
      # Note, $INDEX must be defined higher in the namelist
      val=$(eval echo "\"$val\"")
      echo "NAMELIST: setting ${var} to ${val}"
      declare -g "$var=$val"
    else
      echo "NAMELIST ERR: In $FILETOREAD, invalid variable name: $var"
      exit 1
    fi

    shift 3
  done
}

### -----------------------------------------------------------------------------------

# Check if a file exists, if not exit
# Usage: exit_file_no_exist $OUTPUTSTREAMS
function exit_file_no_exist() {
  if [ ! -f "$1" ]; then
    echo "CHECK_EXIST: File $1 does not exist, exiting..."
    exit 1
  fi
}


# Usage: exit_files_no_exist $path_to_case/SourceMods/src.${lndName}/lnd_comp_mct.F90 $path_to_case/SourceMods/src.${lndName}/lnd_comp_nuopc.F90
function exit_files_no_exist() {
  local checkthisfile
  for checkthisfile in "$@"; do
    if [ -f "$checkthisfile" ]; then
      return 0  # Exit the function, not the script, as soon as we find an existing file
    fi
  done
  # If we've got here, none of the files exist
  echo "CHECK_EXIST: None of the specified files exist, exiting..."
  for checkthisfile in "$@"; do
    echo "  - $checkthisfile"
  done
  exit 1
}


function check_ncl_exit() {
  local ncl_script_name=$1
  local ncl_exit_status=$2
  if [[ $ncl_exit_status -ne 9 ]]; then
    echo "${ncl_script_name} exited with a non-9 error code: ${ncl_exit_status}"
    exit 240
  else
    echo "${ncl_script_name} exited successfully"
  fi
}


### -----------------------------------------------------------------------------------

function get_YYYYMMDD_from_hfile() {
  local FILE=$1
  local DELIM=$2
  local var2=${FILE#*${DELIM}.}
  local THEDATE=$(echo $var2 | cut -c1-4)$(echo $var2 | cut -c6-7)$(echo $var2 | cut -c9-10)
  echo $THEDATE
}

### -----------------------------------------------------------------------------------

sanitize_file () {
  echo "Sanitizing $1"
  sed -i -e '$a\' $1
}

### -----------------------------------------------------------------------------------

# cdv: Verbose cd with error checking.
# Usage: cdv /path/to/dir
# Description:
#   Changes directory to the specified path, printing the destination before switching.
#   Exits the script if the directory change fails.
#
# Example:
#   cdv "$outputdir"
cdv() {
  echo "cd to: $1"
  cd "$1" || { echo "Failed to cd to $1"; exit 1; }
}

### -----------------------------------------------------------------------------------

# Usage--> remove_top_line_from_dates "dates.txt"
remove_top_line_from_dates() {
  local datesfile="$1"

  if [ ! -f "${datesfile}" ]; then
    echo "remove_top_line_from_dates: File does not exist: ${datesfile}"
    exit 1
  fi

  # Remove the top line from the file
  tail -n +2 "${datesfile}" > "${datesfile}.tmp" && mv -v "${datesfile}.tmp" "${datesfile}"
}


# Usage--> longdate=$(get_top_line_from_dates "dates.txt")
get_top_line_from_dates() {
  local datesfile="$1"  # The first argument to the function is the path to the file

  if [ ! -f "${datesfile}" ]; then
    # The >&2 redirects to stderr
    echo "ERROR: get_top_line_from_dates: File does not exist: ${datesfile}" >&2
    exit 1
  fi

  local top_line=$(head -n 1 "${datesfile}")

  if [ -z "$top_line" ]; then
    echo "ERROR: get_top_line_from_dates: ${datesfile} is empty (or the top line is empty)." >&2
    exit 1
  fi

  # If not empty, echo back to the main program to be stored in variable
  echo "$top_line"
}

### -----------------------------------------------------------------------------------

## Input is YYYYMMDDHH
## Sets ${yearstr}, ${monthstr}, ${daystr}, ${hourstr}, ${cyclestrsec}
## Ex: parse_time 2022092812
parse_YYYYMMDDHH () {
  echo "Getting time!"
  local thisDate=$1

  # Make sure 1 (and only 1) args are passed in
  if [ $# -ne 1 ]; then { echo "Error: Exactly one argument is required in function parse_YYYYMMDDHH." ; exit 90; } ; fi

  # Do some simple error trapping on date string to ensure validity
  if [ -z "$thisDate" ]; then { echo "Date string passed in is empty, exiting..." ; exit 91; } ; fi
  if [ ${#thisDate} -ne 10 ]; then { echo "Malformed date string, $thisDate is ${#thisDate} characters, needs 10 (YYYYMMDDHH). Exiting..." ; exit 92; } ; fi
  if [[ -n $(echo $thisDate | tr -d '[0-9]') ]]; then { echo "Malformed date string, $thisDate contains non-numeric values. Exiting..." ; exit 93; } ; fi

  yearstr=${thisDate:0:4}
  monthstr=${thisDate:4:2}
  daystr=${thisDate:6:2}
  cyclestr=${thisDate:8:2}

  # Do some error trapping on returned time values
  if (( yearstr > 3000 || yearstr < 1 )); then { echo "Year set to $yearstr, this sounds wrong, exiting..." ; exit 94; } ; fi
  if (( cyclestr > 23 )); then { echo "Cycle string set to $cyclestr Z, this sounds wrong, exiting..." ; exit 95; } ; fi
}


## Input is HH int
## Sets $cyclestrsec which is zero padded 3600*hr
get_cyclestrsec () {
  local hrstr=$1
  ## Figure out the seconds which correspond to the cycle and zero pad if neces
  cyclestrsec=$(($hrstr*3600))
  while [ ${#cyclestrsec} -lt 5 ];
  do
    cyclestrsec="0"$cyclestrsec
  done
}


function getSSTtime() {
  local islive=$1
  local currtime_UTC_HHMM=$2
  local monthstr=$3
  local daystr=$4
  local yearstr=$5
  local cyclestr=$6

  if [ $islive = true ] ; then
    ## Use currtime_UTC_HHMM to figure out what SST we can download
    ## Current GDAS SST appears 555 after cycle time
    ## First guess is today's dates!
    sstmonthstr=$monthstr
    sstdaystr=$daystr
    sstyearstr=$yearstr
    if [ $currtime_UTC_HHMM -lt 0555 ] ; then
      echo "SSTs are previous days 18Z"
      sstmonthstr=$(date --date="yesterday" -u +%m)
      sstdaystr=$(date --date="yesterday" -u +%d)
      sstyearstr=$(date --date="yesterday" -u +%Y)
      sstcyclestr=18
    elif [ $currtime_UTC_HHMM -lt 1155 ] ; then
      echo "SSTs are current days 00Z"
      sstcyclestr=00
    elif [ $currtime_UTC_HHMM -lt 1755 ] ; then
      echo "SSTs are current days 06Z"
      sstcyclestr=06
    elif [ $currtime_UTC_HHMM -lt 2355 ] ; then
      echo "SSTs are current days 12Z"
      sstcyclestr=12
    elif [ $currtime_UTC_HHMM -ge 2355 ] ; then
      echo "SSTs are current days 18Z"
      sstcyclestr=18
    else
      echo "Can't figure out start time"
      exit 1
    fi
  else
    sstmonthstr=$monthstr
    sstdaystr=$daystr
    sstyearstr=$yearstr
    sstcyclestr=$cyclestr
  fi
  # Global variables set:
  # sstmonthstr
  # sstdaystr
  # sstyearstr
  # sstcyclestr
}


### This is the OG version of this that figures out what the status is by
### counting the number of gzipped files in the run directory
### when the gzip files increments, it means the job is terminated.
run_CIME () {
  # 1 path_to_run_dir
  # 2 ${CIMEsubstring}
  # 3 ${CIMEbatchargs}

  echo "run_CIME: path_to_rundir: "$1
  echo "run_CIME: CIMEsubstring: "$2
  echo "run_CIME: CIMEbatchargs: "$3

  # Get number of log .gz files for sleeping
  echo "Running again!" > $1/testrunning.gz
  local numlogfiles=$(ls ${1}/*.gz | wc -l)
  echo "run_CIME: numlogfiles: $numlogfiles"

  echo "SUBMITTING FORECAST RUN"
  set +e ; ./case.submit $2 --batch-args "${3}" ; set -e

  ## Set up NUKEing
  if [ -f "${1}/NUKE" ] ; then rm -v $1/NUKE ; sleep 5 ; fi
  echo "To NUKE, run \"touch ${1}/NUKE\" "

  ## Hold script while log files from filter run haven't been archived yet
  while [ $(ls ${1}/*.gz | wc -l) == $numlogfiles ]
  do
    if [ -f "${1}/NUKE" ] ; then echo "Nuke sequence initiated, exiting betacast" ; exit ; fi
    sleep 10 ; echo "Sleeping... $(date '+%Y%m%d %H:%M:%S')"
  done
  echo "Run over done sleeping ($(date '+%Y%m%d %H:%M:%S')) will hold for 10 more sec to make sure files moved"
  sleep 10
}


### This is the v2 version of this that figures out what the status is by reading
### the last valid entry in CaseStatus
run_CIME2 () {
  # 1: path_to_run_dir
  # 2: CIMEsubstring
  # 3: CIMEbatchargs
  # 4: exit_on_error (true/false)

  echo "run_CIME2: path_to_rundir: $1"
  echo "run_CIME2: CIMEsubstring: $2"
  echo "run_CIME2: CIMEbatchargs: $3"
  echo "run_CIME2: exit_on_error: $4"

  local run_dir="$1"
  local substr="$2"
  local batchargs="$3"
  local exit_on_error="$4"
  local case_status_str=""
  local case_status=1

  echo "RUN_CIME2: SUBMITTING FORECAST RUN"
  ./case.submit "$substr" --batch-args "$batchargs"

  # Reset NUKE file if it exists
  if [ -f "$run_dir/NUKE" ]; then rm -v "$run_dir/NUKE"; sleep 5; fi
  echo "RUN_CIME2: To NUKE, run \"touch $run_dir/NUKE\" "

  # Wait for model run to complete or fail
  while [ $case_status -eq 1 ]; do
    if [ -f "$run_dir/NUKE" ]; then
      echo "RUN_CIME2: Nuke sequence initiated, exiting betacast"
      exit 1
    fi

    case_status_str=$(grep -a "^20" CaseStatus | tail -2)

    if echo "$case_status_str" | grep -q "case.run success"; then
      case_status=0
    elif echo "$case_status_str" | grep -q "case.run error"; then
      case_status=99
    elif echo "$case_status_str" | grep -q "ERROR: Model did not complete"; then
      case_status=99
    else
      case_status=1
    fi

    sleep 10
    echo "RUN_CIME2: Sleeping... status=$case_status -- $(date '+%Y%m%d %H:%M:%S')"
  done

  if [ $case_status -eq 0 ]; then
    echo "RUN_CIME2: Run completed successfully, sleeping 30s to allow file flush..."
    sleep 30
    return 0
  else
    echo "RUN_CIME2: Run failed"
    if [ "$exit_on_error" = true ]; then
      echo "RUN_CIME2: Exiting due to error in model run; check logs in $run_dir"
      exit 1
    fi
    return 1
  fi
}


# Example usage of the function (assuming COMPRESS_SW and a file path are set)
# COMPRESS_SW is determined by the find_compression_software function
# file_to_decompress="/path/to/your/compressed_file"

# Uncomment the line below to use the function after setting COMPRESS_SW and file_to_decompress appropriately
# decompress_file "$file_to_decompress" "$COMPRESS_SW"

# Function to determine the available compression software
find_compression_software() {
  # Define an array of desired compression algorithms
  local compression_tools=("xz" "zstd" "pigz" "lz4" "gzip")

  # Loop through the array to find the first available tool
  for tool in "${compression_tools[@]}"; do
    if command -v "$tool" > /dev/null 2>&1; then
      # Return the found tool
      echo "$tool"
      return
    fi
  done

  # If we reach this point, no tool was found
  echo "null"
}


# Function to compress a file using the specified compression software
compress_file() {
  local file="$1"
  local compress_sw="$2"
  local start_time end_time original_size compressed_size time_taken speed compression_percentage

  start_time=$(date +%s.%N)
  original_size=$(stat --format=%s "$file")

  case "$compress_sw" in
    zstd)
      zstd -3 -T4 -q --rm "$file"
      compressed_size=$(stat --format=%s "${file}.zst")
      ;;
    pigz)
      pigz "$file"
      compressed_size=$(stat --format=%s "${file}.gz")
      ;;
    gzip)
      gzip "$file"
      compressed_size=$(stat --format=%s "${file}.gz")
      ;;
    xz)
      xz -0 -T8 -f -q "$file"
      compressed_size=$(stat --format=%s "${file}.xz")
      ;;
    lz4)
      lz4 -1 --rm -q "$file"
      compressed_size=$(stat --format=%s "${file}.lz4")
      ;;
    *)
      echo "Unsupported compression software: $compress_sw"
      return 1
      ;;
  esac

  end_time=$(date +%s.%N)
  time_taken=$(echo "$end_time - $start_time" | bc)
  compression_percentage=$(echo "scale=2; (100 * $compressed_size / $original_size)" | bc)
  speed=$(echo "scale=2; $original_size / 1048576 / $time_taken" | bc)

  echo "COMPRESS: Compressed using $compress_sw: ${file}"
  echo "COMPRESS: -- Original size: $original_size bytes, Compressed size: $compressed_size bytes, Compression percentage: ${compression_percentage}%, Time taken: ${time_taken}s, Speed: ${speed}MB/s"
}


uncompress_file() {
  local file="$1"
  local compress_sw="$2"
  local start_time end_time time_taken

  start_time=$(date +%s.%N)

  case "$compress_sw" in
    zstd)
      zstd -d -q "$file" --rm
      ;;
    pigz)
      pigz -d "$file"
      ;;
    gzip)
      gzip -d "$file"
      ;;
    xz)
      xz -d -q "$file"
      ;;
    lz4)
      lz4 -d -q --rm "$file"
      ;;
    *)
      echo "Unsupported compression software for decompression: $compress_sw"
      return 1
      ;;
  esac

  end_time=$(date +%s.%N)
  time_taken=$(echo "$end_time - $start_time" | bc)

  echo "UNCOMPRESS: file uncompressed using $compress_sw: ${file%.*}"
  echo "UNCOMPRESS: -- Time taken: ${time_taken}s"
}


try_uncompress() {
  local file="$1"
  local extension="${file##*.}"

  case "$extension" in
    gz)  uncompress_file "$file" "gzip" ;;
    zst) uncompress_file "$file" "zstd" ;;
    xz)  uncompress_file "$file" "xz" ;;
    *)
      echo "Unsupported extension for decompression: $extension. Exiting."
      exit 1
      ;;
  esac
}


#shopt -s nullglob
#Usage: compress_history $DIR
#Attempt to compress all *.nc in $DIR
#Will try using ncks deflation first
#If that fails, will use pigz
#If that fails, we'll just move on but should deal with it...

#xz -0 -T4 -f -v LANDFALL_20050828_0024.elm.r.2005-08-28-00000.nc
compress_history () {
  echo "Compressing model history files..."
  cd "$1" || exit
  local start_time end_time duration f original_size compressed_size compression_percentage
  for f in *.nc; do
    if [ ! -f "$f" ]; then
      echo "compress_history: No .nc files found. Moving on..."
      continue
    fi
    echo "compress_history: compressing $f with ncks"
    start_time=$(date +%s)
    original_size=$(stat --format=%s "$f")
    if ! ncks -4 -L 1 --rad --no_abc -O "$f" "$f"; then
      rm -v *.ncks.tmp
      echo "compress_history: ncks failed for $f, attempting xz..."
      if ! xz -2 -T8 -q -f "$f"; then
        echo "compress_history: xz failed for $f, attempting zstd..."
        if ! zstd --adapt -T8 -q --rm "$f"; then
          echo "compress_history: zstd failed for $f, attempting pigz..."
          if ! pigz "$f"; then
            echo "compress_history: error: Failed to compress $f with ncks, zstd, and pigz. Moving on..."
          fi #pigz
        fi   #zstd
      fi   #xz
    fi  #ncks
    compressed_size=$(stat --format=%s "$f")
    compression_percentage=$(echo "scale=2; (100 * $compressed_size / $original_size)" | bc)
    end_time=$(date +%s)
    duration=$((end_time - start_time))
    echo "compress_history: ... final percentage: $compression_percentage%, took $duration seconds."
  done
}


#$path_to_nc_files $2
#$tmparchivecdir $1
#$compress_history $3
archive_nudging () {
  (
  echo "Archiving nudging files..."
  cd $2
  mkdir -p $1/nudging
  mv -v *.i.*.nc $1/nudging
  if [ $3 = true ]; then
    # Compress files using lossless compression
    compress_history $1/nudging
  fi
  )
}


#tmparchivecdir = 1
#path_to_case =2
#compress_history_nc = 3
#atmName = 4
#lndName = 5
#rofName = 6
#sstFileIC = 7
archive_inic () {
local compress_sw

  echo "Archiving initial condition files..."
  mkdir -p $1/inic

  # Copy LND initial conditions
  ARCFILE=$(grep ^finidat $2/user_nl_$5 | cut -d "=" -f2)
  if [[ -n "$ARCFILE" ]]; then
    strip_quotes ARCFILE
    echo "Found land initial file: $ARCFILE"
    cp -v "$ARCFILE" "$1/inic"
  else
    echo "Land ARCFILE not found in $2/user_nl_$5"
  fi

  # Copy ATM initial conditions
  ARCFILE=$(grep ^ncdata $2/user_nl_$4 | cut -d "=" -f2)
  if [[ -n "$ARCFILE" ]]; then
    strip_quotes ARCFILE
    echo "Found atm initial file: $ARCFILE"
    cp -v "$ARCFILE" "$1/inic"
  else
    echo "Atm ARCFILE not found in $2/user_nl_$4"
  fi

  # Copy ROF initial conditions
  if [ $do_runoff = true ]; then
    ARCFILE=$(grep ^finidat_rtm $2/user_nl_$6 | cut -d "=" -f2)
    if [[ -n "$ARCFILE" ]]; then
      strip_quotes ARCFILE
      echo "Found rof initial file: $ARCFILE"
      cp -v "$ARCFILE" "$1/inic"
    else
      echo "rof ARCFILE not found in $2/user_nl_$6"
    fi
  fi

  # Copy SST conditions
  cp -v "$7" "$1/inic" || { echo "WARNING: No SST file found in archive_inic"; true; }

  if [ $3 = true ]; then
    # ncks compression
    #compress_history "$1/inic"

    # Since these are just INIC, we can hammer them with a compression algo
    compress_sw=$(find_compression_software)
    if [[ $compress_sw != "null" ]]; then
      echo "Using compression software: $compress_sw"
      for file in "$1"/inic/*; do
        compress_file "$file" "$compress_sw"
      done
    else
      echo "No suitable compression software found."
    fi
  fi
}


## Usage: main_archive $archive_dir $atmName $lndName $rofName
main_archive () {
  echo "Creating archive folder data structure"
  mkdir -p $1
  mkdir -p $1/images
  mkdir -p $1/text
  mkdir -p $1/nl_files
  mkdir -p $1/logs
  mkdir -p $1/betacast

  echo "Moving relevant files to archive folder (${1})"
  mv -v *.$2.h*.nc $1 || true
  mv -v *.$3*.h*.nc $1 || true
  mv -v *.$4*.h*.nc $1 || true
  cp -v *_in seq_maps.rc *_modelio.nml docn.streams.txt.prescribed $1/nl_files || true
  mv -v *.txt $1/text || true
  mv -v *.log.* $1/logs || true
  mv -v timing.*.gz $1/logs || true
  mv -v atm_chunk_costs*.gz $1/logs || true

  mv -v timing/ $1/ || true
}


#Usage: delete_leftovers $rundir $atmName $lndName $rofName
delete_leftovers () {
  (
  echo "Deleting restart/misc. files produced by CESM that aren't needed"
  cd $1
  rm -f -v *.$3*.rh0.*.nc
  rm -f -v *.docn.rs1.*.bin
  rm -f -v *.$2.r.*.nc
  rm -f -v *.$2.rs.*.nc
  rm -f -v *.cpl.r.*.nc
  rm -f -v *.$2.rh3.*.nc
  rm -f -v *.$4.rh0.*.nc
  rm -f -v *.cice.r.*.nc
  rm -f -v rpointer.*
  rm -f -v *.bin
  rm -f -v *.h.*.nc
  rm -f -v *initial_hist*.nc
  )
}


#Usage example:
#delete_except_last "*.clm2.r.*"
#delete_except_last "*.clm2.r.*,*.elm.r.*"
function delete_except_last() {

  # split on multiple input patterns if necessary
  IFS=',' read -r -a patterns <<< "$1"

  # Loop over all matching patterns
  for pattern in "${patterns[@]}"; do

    # Trim leading and trailing spaces from the pattern
    pattern=$(echo "$pattern" | xargs)

    # Collect all files matching the pattern, sorted
    files=($(ls $pattern 2>/dev/null | sort))

    # Count the number of matched files
    file_count=${#files[@]}

    # If there are more than one file, delete all except the last one
    if [ "$file_count" -gt 1 ]; then
      for ((i=0; i<file_count-1; i++)); do
        echo "Deleting: ${files[i]}"
        rm -v "${files[i]}"
      done
    fi
  done
}


# Copy files given a set of patterns (comma-separated)
# Usage: safe_cp2 "*.clm2.r.*,*.elm.r.*" /path/to/destination
safe_cp2() {
  local dest="$2"

  # Split on commas
  IFS=',' read -r -a patterns <<< "$1"

  for pattern in "${patterns[@]}"; do
    # Trim whitespace
    pattern=$(echo "$pattern" | xargs)

    # Enable nullglob to safely handle unmatched globs
    shopt -s nullglob
    files=( $pattern )
    shopt -u nullglob

    if [ ${#files[@]} -eq 0 ]; then
      echo "Warning: No files found matching pattern: $pattern"
    else
      # Sort the files
      IFS=$'\n' sorted_files=($(printf "%s\n" "${files[@]}" | sort))
      cp -v -- "${sorted_files[@]}" "$dest"
    fi
  done
}


# Copy files only if non-empty and exist
# Usage: safe_cp_files file1 file2 ... destdir
safe_cp_files() {
  local dest="${@: -1}"   # Last argument is destination
  local files=("${@:1:$#-1}")  # All but last

  for f in "${files[@]}"; do
    if [ -n "$f" ] && [ -f "$f" ]; then
      echo "Copying: $f → $dest"
      cp -v -- "$f" "$dest"
    else
      echo "Skipping: '$f' (empty or not a file)"
    fi
  done
}


# Usage:
# xmlchange_verbose "JOB_QUEUE" "$RUNQUEUE" "--force"
# xmlchange_verbose "ATM_NCPL" "$ATM_NCPL"
xmlchange_verbose() {
  local variable_name="$1"
  local variable_value="$2"
  shift 2
  local additional_args="$@"

  echo "XMLCHANGE: Setting ${variable_name} to ${variable_value} ${additional_args}"

  if [ -z "$additional_args" ]; then
    ./xmlchange "${variable_name}=${variable_value}"
  else
    ./xmlchange "${additional_args}" "${variable_name}=${variable_value}"
  fi
}


# Usage:
# script_start=$(date +%s)
# script_end=$(date +%s)
# print_elapsed_time "$script_start" "$script_end"
function print_elapsed_time() {
  local time_start=$1
  local time_end=$2

  local time_elapsed=$((time_end - time_start))
  local elapsed_days=$((time_elapsed/86400))
  local elapsed_hours=$((time_elapsed/3600%24))
  local elapsed_minutes=$((time_elapsed/60%60))
  local elapsed_seconds=$((time_elapsed%60))

  echo "TIME ELAPSED: $elapsed_days days, $elapsed_hours hours, $elapsed_minutes minutes, $elapsed_seconds seconds"
}


function compress_single_file() {
  local filename="$1"
  if ! type ncks &> /dev/null ; then
    echo "ncks was NOT FOUND, no compression performed"
  else
    echo "ncks compressing: $filename"
    ncks -O -4 -L 1 "${filename}" "${filename}"
  fi
}


# Timer function to measure execution time of another function
# Usage: timer compress_history /path/to/directory
timer() {
  local start_time=$(date +%s)

  # Call the function passed as argument with all its arguments
  "$@"

  local end_time=$(date +%s)
  local elapsed=$((end_time - start_time))
  echo "Time elapsed for $1: $elapsed seconds"
}

#### NCO functions!
# ncdmnsz $dmn_nm $fl_nm : What is dimension size?
function ncdmnsz { ncks --trd -m -M ${2} | grep -E -i ": ${1}, size =" | cut -f 7 -d ' ' | uniq ; }


# -----------------------------------------------------------------------------
# process_model_times
#
# Determines appropriate forecast cycle times and sets all time-related variables
# needed by the main driver script. Handles both live forecasting and historical
# cases from dates files.
#
# Inputs (passed as parameters):
#   $1 - islive (true/false)
#   $2 - casename (required for non-live runs)
#   $3 - datestemplate (optional, for non-live runs)
#   $4 - numHoursSEStart (for SE filtering calculations)
#   $5 - doFilter (true/false, affects output messaging)
#
# Sets the following global variables:
#   - yearstr, monthstr, daystr, cyclestr, cyclestrsec
#   - se_yearstr, se_monthstr, se_daystr, se_cyclestr, se_cyclestrsec
#   - sstyearstr, sstmonthstr, sstdaystr, sstcyclestr (via getSSTtime)
#   - yestmonthstr, yestdaystr, yestyearstr
#   - datesfile (for non-live runs)
#   - currtime_UTC_HHMM (current UTC time in HHMM format)
# -----------------------------------------------------------------------------
process_model_times() {
  local islive="$1"
  local casename="$2"
  local datestemplate="$3"
  local numHoursSEStart="$4"
  local doFilter="$5"
  local datesbase
  local longdate

  # Get the current time
  currtime_UTC_HHMM=$(date -u +%H%M)

  if [ "$islive" = true ]; then
    echo "Processing live forecast times..."
    # Get current UTC date components
    monthstr=$(date -u +%m)
    daystr=$(date -u +%d)
    yearstr=$(date -u +%Y)

    # Determine latest available GFS cycle based on current time
    # (GFS output lags by ~3.5 hours)
    if [ "$currtime_UTC_HHMM" -lt 0328 ]; then
      echo "12Z cycle (previous day)"
      monthstr=$(date --date="yesterday" -u +%m)
      daystr=$(date --date="yesterday" -u +%d)
      yearstr=$(date --date="yesterday" -u +%Y)
      cyclestr=12
    elif [ "$currtime_UTC_HHMM" -lt 0928 ]; then
      echo "00Z cycle"
      cyclestr=00
    elif [ "$currtime_UTC_HHMM" -lt 1528 ]; then
      echo "00Z cycle"
      cyclestr=00
    elif [ "$currtime_UTC_HHMM" -lt 2128 ]; then
      echo "12Z cycle"
      cyclestr=12
    elif [ "$currtime_UTC_HHMM" -ge 2128 ]; then
      echo "12Z cycle"
      cyclestr=12
    else
      echo "ERROR: Can't determine appropriate start time for current time: $currtime_UTC_HHMM"
      exit 1
    fi
  else
    echo "Processing historical times from dates file..."
    datesbase="dates.${casename}.txt"
    # Locate the dates file
    if [ -f "./dates/${datesbase}" ]; then
      echo "$datesbase exists in dates subdirectory"
      datesfile="./dates/${datesbase}"
    elif [ -f "./${datesbase}" ]; then
      echo "WARNING! $datesbase isn't in dates subdir but exists in home dir!"
      echo "This is allowed for backwards compatibility but is less organized!"
      echo "You should create a subdir called 'dates' and put the dates.CASE.txt file there!"
      datesfile="$datesbase"
    else
      # Try to use template if provided
      if [[ -n "${datestemplate}" && -f "./dates/${datestemplate}" ]]; then
        echo "Didn't find case dates, but you specified a datestemplate..."
        echo "So I'm copying $datestemplate to this case and using that..."
        cp -v "./dates/${datestemplate}" "./dates/${datesbase}"
        datesfile="./dates/${datesbase}"
      else
        echo "ERROR: Can't find a dates file OR template AND run isn't live. Exiting..."
        exit 1
      fi
    fi

    echo "Using dates file: $datesfile"

    # Parse the top line from dates file
    longdate=$(get_top_line_from_dates "$datesfile")
    echo "Getting parsed time from $longdate"
    parse_YYYYMMDDHH "$longdate"
    echo "From datesfile, read in: $yearstr $monthstr $daystr ${cyclestr}Z"

  fi

  # Calculate cycle seconds
  get_cyclestrsec "$cyclestr"

  # Calculate SE (Spectral Element) start times for filtering
  if [ "$numHoursSEStart" -lt 6 ]; then
    # Calculate SE start time (3 hours after main cycle)
    se_cyclestr=$((10#$cyclestr + 3))

    # Zero pad if necessary
    while [ ${#se_cyclestr} -lt 2 ]; do
      se_cyclestr="0$se_cyclestr"
    done

    # For now, assume same day (this logic might need enhancement for edge cases)
    se_monthstr="$monthstr"
    se_daystr="$daystr"
    se_yearstr="$yearstr"

    # Calculate seconds and zero pad
    se_cyclestrsec=$((10#$se_cyclestr * 3600))
    while [ ${#se_cyclestrsec} -lt 5 ]; do
      se_cyclestrsec="0$se_cyclestrsec"
    done
  else
    echo "ERROR: SE forecast lead time too long, 18Z cycle causes trouble"
    echo "Not supported."
    exit 1
  fi

  # Get SST times using existing function
  getSSTtime "$islive" "$currtime_UTC_HHMM" "$monthstr" "$daystr" "$yearstr" "$cyclestr"

  # Set yesterday's date strings
  yestmonthstr=$(date --date="yesterday" -u +%m)
  yestdaystr=$(date --date="yesterday" -u +%d)
  yestyearstr=$(date --date="yesterday" -u +%Y)

  # Print summary
  echo "=== TIME PROCESSING SUMMARY ==="
  echo "ATM init data: $yearstr $monthstr $daystr $cyclestr Z ($cyclestrsec seconds)"
  echo "SST init data: $sstyearstr $sstmonthstr $sstdaystr $sstcyclestr Z"

  if [ "$doFilter" = true ]; then
    echo "Filter: True model init will occur at $se_yearstr $se_monthstr $se_daystr $se_cyclestr Z ($se_cyclestrsec seconds)"
  else
    echo "No filter: True model init will occur at $se_yearstr $se_monthstr $se_daystr $cyclestr Z ($cyclestrsec seconds)"
  fi
  echo "==============================="
}