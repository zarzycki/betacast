#!/bin/bash

# Function to run Python test
run_python_test() {

    mkdir -pv $DEBUG_FILE_DIR/ERA5_PY/

    python ${BETACAST}/land-spinup/gen_datm/gen-datm/gen-forcing.py \
        --era5_file="${TEST_FILES_DIR}/out.2003.01.nc" \
        --year=2003 \
        --month=01 \
        --outdirbase="${DEBUG_FILE_DIR}/ERA5_PY/" \
        --do_q \
        --do_flds \
        --greg_to_noleap \
        --convert_nc3

    return 0
}

# Optional cleanup function
cleanup() {
    rm -frv $DEBUG_FILE_DIR/ERA5_PY/
}
