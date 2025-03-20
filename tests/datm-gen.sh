#!/bin/bash

# Environment setup
export TEST_FILES_DIR="/glade/campaign/univ/upsu0032/betacast_test_files/"
export DEBUG_FILE_DIR="/glade/derecho/scratch/zarzycki/tmp_betacast/"
export BETACAST="/glade/u/home/zarzycki/betacast/"

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

# Function to run NCL test
run_ncl_test() {

    mkdir -pv $DEBUG_FILE_DIR/ERA5_NCL/

    ncl ${BETACAST}/land-spinup/gen_datm/gen-datm/gen-forcing.ncl YYYY=2003 \
        'MM="01"' \
        'RAWERA5FILE="'${TEST_FILES_DIR}'/out.2003.01.nc"' \
        'outdirbase="'${DEBUG_FILE_DIR}'/ERA5_NCL/"'

    return 9

}

# Function to validate results
run_validation() {
    python ${BETACAST}/py_testing/check_same.py $DEBUG_FILE_DIR/ERA5_PY/TPQW/CMZERA5.v0.c2021.0.5d.TPQWL.2003-01.nc $DEBUG_FILE_DIR/ERA5_NCL/TPQW/CMZERA5.v0.c2021.0.5d.TPQWL.2003-01.nc
    python ${BETACAST}/py_testing/check_same.py $DEBUG_FILE_DIR/ERA5_PY/Precip/CMZERA5.v0.c2021.0.5d.Prec.2003-01.nc $DEBUG_FILE_DIR/ERA5_NCL/Precip/CMZERA5.v0.c2021.0.5d.Prec.2003-01.nc
    python ${BETACAST}/py_testing/check_same.py $DEBUG_FILE_DIR/ERA5_PY/Solar/CMZERA5.v0.c2021.0.5d.Solar.2003-01.nc $DEBUG_FILE_DIR/ERA5_NCL/Solar/CMZERA5.v0.c2021.0.5d.Solar.2003-01.nc
}

# Optional cleanup function
cleanup() {
    rm -frv $DEBUG_FILE_DIR/ERA5_PY/
    rm -frv $DEBUG_FILE_DIR/ERA5_NCL/
}
