#!/bin/bash

# Function to run Python test
run_python_test() {
    python ${BETACAST}/py_atm_to_cam/atm_to_cam.py \
        --datasource "ERA5RDA" \
        --numlevels 32 \
        --YYYYMMDDHH 2005082800 \
        --data_filename "${TEST_FILES_DIR}/ds633.0/e5.oper.invariant/197901/e5.oper.invariant.128_129_z.ll025sc.1979010100_1979010100.nc" \
        --wgt_filename "${TEST_FILES_DIR}/map_gfs_0.25x0.25_TO_mpasa120_patc.nc" \
        --dycore "mpas" \
        --RDADIR "${TEST_FILES_DIR}/ds633.0/" \
        --mpas_as_cam \
        --compress_file \
        --write_floats \
        --add_cloud_vars \
        --adjust_config "" \
        --model_topo_file "${TEST_FILES_DIR}/mpasa120.CFSR.L32.nc" \
        --se_inic "${DEBUG_FILE_DIR}/py_final.nc" \
        --verbose \
        --write_debug_files \
        --write_debug_dir "${DEBUG_FILE_DIR}"
}

# Function to run NCL test
run_ncl_test() {
    ncl -n ${BETACAST}/atm_to_cam/atm_to_cam.ncl \
        'datasource="ERA5RDA"' \
        numlevels=32 \
        YYYYMMDDHH=2005082800 \
        'dycore="mpas"' \
        'data_filename="'${TEST_FILES_DIR}'/ds633.0/e5.oper.invariant/197901/e5.oper.invariant.128_129_z.ll025sc.1979010100_1979010100.nc"' \
        'wgt_filename="'${TEST_FILES_DIR}'/map_gfs_0.25x0.25_TO_mpasa120_patc.nc"' \
        'RDADIR="'${TEST_FILES_DIR}'/ds633.0/"' \
        mpas_as_cam=True \
        write_floats=True \
        'model_topo_file="'${TEST_FILES_DIR}'/mpasa120.CFSR.L32.nc"' \
        'adjust_config=""' \
        write_debug_files=True \
        'write_debug_dir="'${DEBUG_FILE_DIR}'"' \
        'se_inic="'${DEBUG_FILE_DIR}'/ncl_final.nc"'
}

# Function to validate results
run_validation() {
    python ${BETACAST}/py_testing/check_same.py ${DEBUG_FILE_DIR}/py_final.nc ${DEBUG_FILE_DIR}/ncl_final.nc
}

# Optional cleanup function
cleanup() {
    rm -f ${DEBUG_FILE_DIR}/py_final.nc ${DEBUG_FILE_DIR}/ncl_final.nc
}
