#!/bin/bash

# Environment setup
export TEST_FILES_DIR="/glade/derecho/scratch/zarzycki/test_files/"
export DEBUG_FILE_DIR="/glade/derecho/scratch/zarzycki/tmp_betacast/"
export BETACAST="/glade/derecho/scratch/zarzycki/betacast/"

DATE=2024081400
NLEV=10
DATASOURCE="GFS"

# Function to run Python test
run_python_test() {
    python ${BETACAST}/py_atm_to_cam/atm_to_cam.py \
        --datasource $DATASOURCE \
        --numlevels $NLEV \
        --YYYYMMDDHH $DATE \
        --data_filename "${TEST_FILES_DIR}/gfs_atm_2024081400.grib2" \
        --wgt_filename "${TEST_FILES_DIR}/map_gfs_0.25x0.25_TO_ne30_patc.nc" \
        --dycore "se" \
        --write_floats \
        --add_cloud_vars \
        --adjust_config "a" \
        --model_topo_file "${TEST_FILES_DIR}/ne30np4_nc3000_Co060_Fi001_PF_nullRR_Nsw042_20171020.nc" \
        --se_inic "${DEBUG_FILE_DIR}/py_final.nc" \
        --verbose \
        --write_debug_files \
        --write_debug_dir "${DEBUG_FILE_DIR}"
}

# Function to run NCL test
run_ncl_test() {
    ncl -n ${BETACAST}/atm_to_cam/atm_to_cam.ncl \
        'datasource="'${DATASOURCE}'"' \
        numlevels=$NLEV \
        YYYYMMDDHH=$DATE \
        'dycore="se"' \
        'data_filename="'${TEST_FILES_DIR}'/gfs_atm_2024081400.grib2"' \
        'wgt_filename="'${TEST_FILES_DIR}'/map_gfs_0.25x0.25_TO_ne30_patc.nc"' \
        'model_topo_file="'${TEST_FILES_DIR}'/ne30np4_nc3000_Co060_Fi001_PF_nullRR_Nsw042_20171020.nc"' \
        'adjust_config="a"' \
        write_floats=True \
        add_cloud_vars=True \
        write_debug_files=True \
        'write_debug_dir="'${DEBUG_FILE_DIR}'"' \
        'se_inic = "'${DEBUG_FILE_DIR}'/ncl_final.nc"'
}

# Function to validate results
run_validation() {
    python ${BETACAST}/py_testing/check_same.py ${DEBUG_FILE_DIR}/py_final.nc ${DEBUG_FILE_DIR}/ncl_final.nc
}

# Optional cleanup function
cleanup() {
    rm -f ${DEBUG_FILE_DIR}/py_final.nc ${DEBUG_FILE_DIR}/ncl_final.nc
}
