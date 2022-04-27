#!/bin/bash

### -----------------------------------------------------------------------------------

# Strip surrounding quotes from string [$1: variable name]
function strip_quotes() {
  local -n var="$1"
  [[ "${var}" == \"*\" || "${var}" == \'*\' ]] && var="${var:1:-1}"
}
export -f strip_quotes

### -----------------------------------------------------------------------------------

# Check if a boolean is set via 0/1 and correct.
# Usage: check_bool "do_tracking" $do_tracking
check_bool() {

  local cborig
  local cbnew
  var_string=$1
  var_totest=$2
  cborig=$2

  # Convert to lowercase
  var_totest=`echo ${var_totest,,}`
  
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
export -f check_bool

### -----------------------------------------------------------------------------------
