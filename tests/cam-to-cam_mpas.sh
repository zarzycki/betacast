#!/bin/bash

# Function to run Python test
run_python_test() {
    python ${BETACAST}/py_atm_to_cam/atm_to_cam.py \
        --datasource "CAM" \
        --numlevels 10 \
        --YYYYMMDDHH 2005081800 \
        --data_filename "${TEST_FILES_DIR}/cam_to_cam/CHEY.VR28.NATL.REF.CAM5.4CLM5.0.dtime900.cam.h2.2005-08-08-00000.nc" \
        --wgt_filename "${TEST_FILES_DIR}/map_gfs_0.25x0.25_TO_mpasa120_patc.nc" \
        --dycore "mpas" \
        --compress_file \
        --write_floats \
        --add_cloud_vars \
        --adjust_config "a" \
        --mod_remap_file "${TEST_FILES_DIR}/map_ne0np4natlanticref.ne30x4_TO_era5_0.25x0.25_patc.nc" \
        --mod_in_topo "${TEST_FILES_DIR}/topo_ne0np4natlanticref.ne30x4_smooth.nc" \
        --model_topo_file "${TEST_FILES_DIR}/mpasa120.CFSR.L32.nc" \
        --se_inic "${DEBUG_FILE_DIR}/py_final.nc" \
        --verbose \
        --write_debug_files \
        --write_debug_dir "${DEBUG_FILE_DIR}"
}

# Function to run NCL test
run_ncl_test() {
    ncl -n ${BETACAST}/atm_to_cam/atm_to_cam.ncl \
        'datasource="CAM"' \
        numlevels=10 \
        YYYYMMDDHH=2005081800 \
        'dycore="mpas"' \
        'data_filename="'${TEST_FILES_DIR}'/cam_to_cam/CHEY.VR28.NATL.REF.CAM5.4CLM5.0.dtime900.cam.h2.2005-08-08-00000.nc"' \
        'wgt_filename="'${TEST_FILES_DIR}'/map_gfs_0.25x0.25_TO_mpasa120_patc.nc"' \
        'adjust_config="a"' \
        compress_file=False \
        write_floats=True \
        'mod_remap_file="'${TEST_FILES_DIR}'/map_ne0np4natlanticref.ne30x4_TO_era5_0.25x0.25_patc.nc"' \
        'mod_in_topo="'${TEST_FILES_DIR}'/topo_ne0np4natlanticref.ne30x4_smooth.nc"' \
        'model_topo_file="'${TEST_FILES_DIR}'/mpasa120.CFSR.L32.nc"' \
        write_debug_files=True \
        'write_debug_dir="'${DEBUG_FILE_DIR}'"' \
        'se_inic="'${DEBUG_FILE_DIR}'/ncl_final.nc"'
}

# Function to validate results
run_validation() {
    python ${BETACAST}/py_testing/check_same.py ${DEBUG_FILE_DIR}/py_final.nc ${DEBUG_FILE_DIR}/ncl_final.nc u,qv,rho,theta
}

# Optional cleanup function
cleanup() {
    rm -f ${DEBUG_FILE_DIR}/py_final.nc ${DEBUG_FILE_DIR}/ncl_final.nc
}
