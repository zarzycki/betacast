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

# Strip surrounding quotes from string [$1: variable name]
function strip_quotes() {
  local -n var="$1"
  [[ "${var}" == \"*\" || "${var}" == \'*\' ]] && var="${var:1:-1}"
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
  sanitize_file ${FILETOREAD}

  # Note, ___ will be converted to a space. Namelists cannot have whitespace due to
  # parsing on whitespaces...
  echo "Reading namelist ${FILETOREAD}..."
  local inputstream=$(cat ${FILETOREAD} | grep -v "^#")
  inputstream="${inputstream//=/ = }"
  #echo $inputstream
  set -- $inputstream
  while [ $1 ]
   do
     if [ "${2}" != "=" ] ; then echo "Uh oh, $1, $2, $3!" ; exit ; fi
     echo "NAMELIST: setting ${1} to ${3//___/ }"
     #eval $1=$3
     eval $1="${3//___/ }"
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
  local currtime=$2
  local monthstr=$3
  local daystr=$4
  local yearstr=$5
  local cyclestr=$6

  if [ $islive = true ] ; then
    ## Use currtime to figure out what SST we can download
    ## Current GDAS SST appears 555 after cycle time
    ## First guess is today's dates!
    sstmonthstr=$monthstr
    sstdaystr=$daystr
    sstyearstr=$yearstr
    if [ $currtime -lt 0555 ] ; then
      echo "SSTs are previous days 18Z"
      sstmonthstr=$(date --date="yesterday" -u +%m)
      sstdaystr=$(date --date="yesterday" -u +%d)
      sstyearstr=$(date --date="yesterday" -u +%Y)
      sstcyclestr=18
    elif [ $currtime -lt 1155 ] ; then
      echo "SSTs are current days 00Z"
      sstcyclestr=00
    elif [ $currtime -lt 1755 ] ; then
      echo "SSTs are current days 06Z"
      sstcyclestr=06
    elif [ $currtime -lt 2355 ] ; then
      echo "SSTs are current days 12Z"
      sstcyclestr=12
    elif [ $currtime -ge 2355 ] ; then
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
  # 1 path_to_run_dir
  # 2 ${CIMEsubstring}
  # 3 ${CIMEbatchargs}
  # 4 exit_on_error

  echo "run_CIME2: path_to_rundir: "$1
  echo "run_CIME2: CIMEsubstring: "$2
  echo "run_CIME2: CIMEbatchargs: "$3
  echo "run_CIME2: exit_on_error: "$4

  # Get number of log .gz files for sleeping

  local CASESTR=""

  echo "RUN_CIME2: SUBMITTING FORECAST RUN"
  set +e ; ./case.submit $2 --batch-args "${3}" ; set -e

  ## Set up NUKEing
  if [ -f "${1}/NUKE" ] ; then rm -v $1/NUKE ; sleep 5 ; fi
  echo "RUN_CIME2: To NUKE, run \"touch ${1}/NUKE\" "

  ## Hold script while log files from filter run haven't been archived yet
  CIMESTATUS=1
  while [ $CIMESTATUS == 1 ]
  do
    if [ -f "${1}/NUKE" ] ; then echo "RUN_CIME2: Nuke sequence initiated, exiting betacast" ; exit ; fi

    # Get the last valid line from the CaseStatus file...
    CASESTR=$(grep "^20" CaseStatus | tail -1)

    if [[ "$CASESTR" == *"case.run success"* ]]; then
      CIMESTATUS=0
    elif [[ "$CASESTR" == *"case.run error"* ]]; then
      CIMESTATUS=99
    else
      CIMESTATUS=1
    fi

    sleep 10 ; echo "RUN_CIME2: Sleeping... CIMESTATUS: $CIMESTATUS -- $(date '+%Y%m%d %H:%M:%S')"
  done

  if [ $CIMESTATUS -eq 0 ]; then
    echo "RUN_CIME2: Run over done sleeping ($(date '+%Y%m%d %H:%M:%S')) will hold for 30 more sec to make sure files moved"
    sleep 30
  else
    echo "RUN_CIME2: Uh oh, something wrong!"
    if [ $4 = true ]; then
      echo "RUN_CIME2: Exiting because of error with model run, check log files in $path_to_rundir"
      exit 1
    fi
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
    local compression_tools=("xz" "zstd" "pigz" "gzip")

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
            zstd -3 -T8 --rm "$file"
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
            xz -0 -T8 -f "$file"
            compressed_size=$(stat --format=%s "${file}.xz")
            ;;
        *)
            echo "Unsupported compression software: $compress_sw"
            return 1
            ;;
    esac

    end_time=$(date +%s.%N)
    time_taken=$(echo "$end_time - $start_time" | bc)
    compression_percentage=$(echo "scale=2; ($compressed_size / $original_size) * 100" | bc)
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
            zstd -d "$file" --rm
            ;;
        pigz)
            pigz -d "$file"
            ;;
        gzip)
            gzip -d "$file"
            ;;
        xz)
            xz -d "$file"
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
  for f in *.nc; do
    if [ ! -f "$f" ]; then
      echo "compress_history: No .nc files found. Moving on..."
      continue
    fi
    echo "compress_history: compressing $f with ncks"
    if ! ncks -4 -L 1 --rad --no_abc -O "$f" "$f"; then
      rm -v *.ncks.tmp
      echo "compress_history: ncks failed for $f, attempting xz..."
      if ! xz -2 -T8 -q -f "$f"; then
        echo "compress_history: xz failed for $f, attempting zstd..."
        if ! zstd --adapt -T8 -q --rm "$f"; then
          echo "compress_history: zstd failed for $f, attempting pigz..."
          if ! pigz "$f"; then
            echo "compress_history: error: Failed to compress $f with ncks, zstd, and pigz. Moving on..."
          fi
        fi
      fi
    fi
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
  cp -v $7 $1/inic

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

  echo "Moving relevant files to archive folder"
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
  /usr/bin/time -f "Time elapsed for $1: %E ---> Maximum memory used: %M kilobytes" "$@" 2>&1
}

#### NCO functions!
# ncdmnsz $dmn_nm $fl_nm : What is dimension size?
function ncdmnsz { ncks --trd -m -M ${2} | grep -E -i ": ${1}, size =" | cut -f 7 -d ' ' | uniq ; }

