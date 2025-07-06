#!/bin/bash

# Additional configuration for filter test
export filterHourLength=6
export filtTcut=6
export filtfile_name="${TEST_FILES_DIR}/F-betacast-F2010-CICE.eam.h0.1996-12-29-00000.nc"
export sePostFilterIC="${TEST_FILES_DIR}/F-betacast-F2010-CICE_INIC.nc"

# Function to run Python test
run_python_test() {
  # Copy the initial condition file
  cp -v ${sePostFilterIC} ${DEBUG_FILE_DIR}/python_filtered.nc

  # Run the Python filter
  python ${BETACAST}/py_atm_to_cam/filter/lowmemfilter.py \
    --endhour ${filterHourLength} \
    --tcut ${filtTcut} \
    --filtfile_name "${filtfile_name}" \
    --writefile_name "${DEBUG_FILE_DIR}/python_filtered.nc"

  return 0
}

# Function to run NCL test
run_ncl_test() {
  # Copy the initial condition file
  cp -v ${sePostFilterIC} ${DEBUG_FILE_DIR}/ncl_filtered.nc

  # Run the NCL filter
  ncl ${BETACAST}/atm_to_cam/filter/lowmemfilter.ncl \
    endhour=${filterHourLength} \
    tcut=${filtTcut} \
    'filtfile_name = "'${filtfile_name}'"' \
    'writefile_name = "'${DEBUG_FILE_DIR}'/ncl_filtered.nc"'

  return 9
}

# Function to validate results
run_validation() {
  echo "Validating Python filtered file with NCL filtered file..."
  python ${BETACAST}/py_testing/check_same.py ${DEBUG_FILE_DIR}/python_filtered.nc ${DEBUG_FILE_DIR}/ncl_filtered.nc || return 1

  return 0
}

# Optional cleanup function
cleanup() {
  rm -f ${DEBUG_FILE_DIR}/python_filtered.nc ${DEBUG_FILE_DIR}/ncl_filtered.nc
}
