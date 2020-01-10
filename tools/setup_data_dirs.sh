#!/bin/bash
# Usage ./setup_data_dirs.sh ../machine_files/machine.cheyenne
# This script just builds data structure from machine file so scripts can write appropriately

set -e

# Set files, in reality order doesn't matter
MACHINEFILE=${1}
# If relative path, convert to absolute path
if [[ "$MACHINEFILE" != /* ]] && [[ "$MACHINEFILE" != ~* ]]; then MACHINEFILE=${PWD}/${MACHINEFILE}; fi
echo $MACHINEFILE

# Sanitize namelist files (add carriage return to end)
sed -i -e '$a\' ${MACHINEFILE}
#'

# Note, ___ will be converted to a space. Namelists cannot have whitespace due to
# parsing on whitespaces...
echo "Reading namelist ${NAMELISTFILE}..."
inputstream=`cat ${NAMELISTFILE} ${MACHINEFILE} | grep -v "^#"`
set -- $inputstream
while [ $1 ]
do
  echo "NAMELIST: setting ${1} to ${3/___/ }"
  eval $1="${3/___/ }"  
  shift 3
done

set -u  # turn on crashes for unbound variables in bash

### MAKE DIRS
mkdir -vp $path_to_inputdata
mkdir -vp $pathToINICfiles
mkdir -vp $pathToSSTfiles
mkdir -vp $gfs_files_path
mkdir -vp $era_files_path
mkdir -vp $sst_files_path

exit 0

