#!/bin/bash
# Environment setup
export TEST_FILES_DIR="/glade/u/home/zarzycki/scratch/test_files/"
export DEBUG_FILE_DIR="/glade/u/home/zarzycki/scratch/tmp_betacast/"
export BETACAST="/glade/u/home/zarzycki/scratch/betacast/"

# Additional variables specific to this test
DATE="2023052300"
RES="180x360"
DATASRC="NOAAOI"

# Function to run Python test
run_python_test() {
    python ${BETACAST}/py_sst_to_cam/sst_to_cam.py \
        --initdate ${DATE} \
        --predict_docn 0 \
        --inputres ${RES} \
        --datasource ${DATASRC} \
        --sstDataFile "${TEST_FILES_DIR}/sst.day.mean.2023.nc" \
        --iceDataFile "${TEST_FILES_DIR}/icec.day.mean.2023.nc" \
        --SST_write_file "${DEBUG_FILE_DIR}/sst_python_final.nc" \
        --smooth_ice \
        --smooth_iter 3 \
        --verbose
}

# Function to run NCL test
run_ncl_test() {
    ncl ${BETACAST}/sst_to_cam/sst_interp.ncl \
        'initdate="'${DATE}'"' \
        predict_docn=0 \
        'inputres="'${RES}'"' \
        'datasource="'${DATASRC}'"' \
        'sstDataFile="'${TEST_FILES_DIR}'/sst.day.mean.2023.nc"' \
        'iceDataFile="'${TEST_FILES_DIR}'/icec.day.mean.2023.nc"' \
        'SST_write_file="'${DEBUG_FILE_DIR}'/sst_ncl_final.nc"'
}

# Function to validate results
run_validation() {
    python ${BETACAST}/py_testing/check_same.py ${DEBUG_FILE_DIR}/sst_python_final.nc ${DEBUG_FILE_DIR}/sst_ncl_final.nc
}

# Optional cleanup function
cleanup() {
    rm -f ${DEBUG_FILE_DIR}/sst_python_final.nc ${DEBUG_FILE_DIR}/sst_ncl_final.nc
}
