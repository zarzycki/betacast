#!/bin/bash

# Function to run Python test
run_python_test() {
    python ${BETACAST}/py_atm_to_cam/atm_to_cam.py \
        --datasource "ERA5RDA" \
        --numlevels 10 \
        --YYYYMMDDHH 2005082800 \
        --data_filename "${TEST_FILES_DIR}/ds633.0/e5.oper.invariant/197901/e5.oper.invariant.128_129_z.ll025sc.1979010100_1979010100.nc" \
        --wgt_filename "${TEST_FILES_DIR}/map_era5_0.25_TO_fv0.9x1.25_patc.nc" \
        --dycore "fv" \
        --RDADIR "${TEST_FILES_DIR}/ds633.0/" \
        --compress_file \
        --write_floats \
        --add_cloud_vars \
        --adjust_config "a" \
        --model_topo_file "${TEST_FILES_DIR}/fv_0.9x1.25_nc3000_Nsw042_Nrs008_Co060_Fi001_ZR_sgh30_24km_GRNL_c170103.nc" \
        --se_inic "${DEBUG_FILE_DIR}/py_final.nc" \
        --verbose \
        --write_debug_files \
        --write_debug_dir "${DEBUG_FILE_DIR}"
}

# Function to run NCL test
run_ncl_test() {
    ncl -n ${BETACAST}/atm_to_cam/atm_to_cam.ncl \
        'datasource="ERA5RDA"' \
        numlevels=10 \
        YYYYMMDDHH=2005082800 \
        'dycore="fv"' \
        'data_filename="'${TEST_FILES_DIR}'/ds633.0/e5.oper.invariant/197901/e5.oper.invariant.128_129_z.ll025sc.1979010100_1979010100.nc"' \
        'wgt_filename="'${TEST_FILES_DIR}'/map_era5_0.25_TO_fv0.9x1.25_patc.nc"' \
        'RDADIR="'${TEST_FILES_DIR}'/ds633.0/"' \
        'model_topo_file="'${TEST_FILES_DIR}'/fv_0.9x1.25_nc3000_Nsw042_Nrs008_Co060_Fi001_ZR_sgh30_24km_GRNL_c170103.nc"' \
        'adjust_config="a"' \
        compress_file=True \
        write_floats=True \
        add_cloud_vars=True \
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
