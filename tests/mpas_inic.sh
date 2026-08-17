#!/bin/bash

# Function to run Python test
run_python_test() {
    python ${BETACAST}/atm_to_cam/atm_to_cam.py \
        --datasource "ERA5RDA" \
        --numlevels 32 \
        --YYYYMMDDHH 2005082800 \
        --data_filename "${TEST_FILES_DIR}/ds633.0/e5.oper.invariant/197901/e5.oper.invariant.128_129_z.ll025sc.1979010100_1979010100.nc" \
        --wgt_filename "${TEST_FILES_DIR}/map_gfs_0.25x0.25_TO_mpasa120_patc.nc" \
        --dycore "mpas" \
        --RDADIR "${TEST_FILES_DIR}/ds633.0/" \
        --compress_file \
        --write_floats \
        --add_cloud_vars \
        --adjust_config "" \
        --model_topo_file "${TEST_FILES_DIR}/mpasa120.CFSR.L32_no-met.nc" \
        --se_inic "${DEBUG_FILE_DIR}/py_final.nc" \
        --verbose \
        --write_debug_files \
        --write_debug_dir "${DEBUG_FILE_DIR}"
}

# Optional cleanup function
cleanup() {
    rm -f ${DEBUG_FILE_DIR}/py_final.nc
}
