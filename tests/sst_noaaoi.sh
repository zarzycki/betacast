#!/bin/bash

# Additional variables specific to this test
DATE="2023052300"
RES="180x360"
DATASRC="NOAAOI"

# Function to run Python test
run_python_test() {
    python ${BETACAST}/sst_to_cam/sst_to_cam.py \
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

# Optional cleanup function
cleanup() {
    rm -f ${DEBUG_FILE_DIR}/sst_python_final.nc
}
