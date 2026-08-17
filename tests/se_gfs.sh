#!/bin/bash

DATE=2024081400
NLEV=10
DATASOURCE="GFS"

# Function to run Python test
run_python_test() {
    python ${BETACAST}/atm_to_cam/atm_to_cam.py \
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

# Optional cleanup function
cleanup() {
    rm -f ${DEBUG_FILE_DIR}/py_final.nc
}
