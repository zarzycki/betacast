#!/bin/bash

# Test: gen-sst-domain Python vs NCL
# Generates domain + SCRIP files at 180x360 resolution and compares outputs.
# Fully self-contained — no external data required.

RES="180x360"

# Function to run Python test
run_python_test() {
    python ${BETACAST}/py_sst_to_cam/gen-sst-domain.py --inputres ${RES}

    # Copy outputs to DEBUG_FILE_DIR for comparison
    cp ${BETACAST}/grids/domains/domain.ocn.${RES}.nc ${DEBUG_FILE_DIR}/domain.ocn.py.nc
    cp ${BETACAST}/grids/domains/scrip.ocn.${RES}.nc ${DEBUG_FILE_DIR}/scrip.ocn.py.nc
}

# Function to run NCL test
run_ncl_test() {
    # NCL writes to ./domains/ relative to CWD
    pushd ${DEBUG_FILE_DIR} > /dev/null
    ncl ${BETACAST}/sst_to_cam/gen-sst-domain.ncl 'inputres="'${RES}'"'
    local ncl_status=$?
    popd > /dev/null

    # Copy NCL outputs for comparison
    cp ${DEBUG_FILE_DIR}/domains/domain.ocn.${RES}.nc ${DEBUG_FILE_DIR}/domain.ocn.ncl.nc
    cp ${DEBUG_FILE_DIR}/domains/scrip.ocn.${RES}.nc ${DEBUG_FILE_DIR}/scrip.ocn.ncl.nc

    return ${ncl_status}
}

# Function to validate results
run_validation() {
    echo "Comparing domain files..."
    python ${BETACAST}/py_testing/check_same.py \
        ${DEBUG_FILE_DIR}/domain.ocn.py.nc \
        ${DEBUG_FILE_DIR}/domain.ocn.ncl.nc || return 1

    echo "Comparing SCRIP files..."
    python ${BETACAST}/py_testing/check_same.py \
        ${DEBUG_FILE_DIR}/scrip.ocn.py.nc \
        ${DEBUG_FILE_DIR}/scrip.ocn.ncl.nc || return 1

    return 0
}

# Optional cleanup function
cleanup() {
    rm -f ${DEBUG_FILE_DIR}/domain.ocn.py.nc ${DEBUG_FILE_DIR}/domain.ocn.ncl.nc
    rm -f ${DEBUG_FILE_DIR}/scrip.ocn.py.nc ${DEBUG_FILE_DIR}/scrip.ocn.ncl.nc
    rm -rf ${DEBUG_FILE_DIR}/domains
}
