#!/bin/bash
# Environment setup
export TEST_FILES_DIR="/glade/u/home/zarzycki/scratch/test_files/"
export DEBUG_FILE_DIR="/glade/u/home/zarzycki/scratch/tmp_betacast/"
export BETACAST="/glade/u/home/zarzycki/scratch/betacast/"

# Additional configuration
export NCL_BEFORE_ATM="${DEBUG_FILE_DIR}/ncl_final.nc"
export PY_BEFORE_ATM="${DEBUG_FILE_DIR}/py_final.nc"
export PERTURB_NAMELIST="/glade/u/home/zarzycki/scratch/betacast/namelists/perturb.101.plus4K.nl"
export MODEL_GRID="/glade/p/cesmdata/inputdata/share/scripgrids/ne16np4_scrip_171002.nc"

# Function to run NCL test
run_ncl_test() {
    # First generate initial conditions
    ncl -n ${BETACAST}/atm_to_cam/atm_to_cam.ncl \
        'datasource="ERA5RDA"' \
        numlevels=10 \
        YYYYMMDDHH=2005082800 \
        'dycore="se"' \
        'data_filename="'${TEST_FILES_DIR}'/ds633.0/e5.oper.invariant/197901/e5.oper.invariant.128_129_z.ll025sc.1979010100_1979010100.nc"' \
        'wgt_filename="'${TEST_FILES_DIR}'/map_era5_0.25x0.25_TO_ne16_patc.nc"' \
        'RDADIR="'${TEST_FILES_DIR}'/ds633.0/"' \
        mpas_as_cam=True \
        'model_topo_file="'${TEST_FILES_DIR}'/ne16np4_nc3000_Co120_Fi001_PF_nullRR_Nsw084_20171012.nc"' \
        'adjust_config="a"' \
        compress_file=True \
        write_floats=True \
        add_cloud_vars=True \
        write_debug_files=True \
        'write_debug_dir="'${DEBUG_FILE_DIR}'"' \
        'se_inic="'${NCL_BEFORE_ATM}'"'

    # Then add perturbations
    ncl -n ${BETACAST}/atm_to_cam/perturb/add_perturbations_to_cam.ncl \
        'BEFOREPERTFILE="'${NCL_BEFORE_ATM}'"' \
        'AFTERPERTFILE="'${NCL_BEFORE_ATM}'_PERT.nc"' \
        'gridfile="'${MODEL_GRID}'"' \
        'MAPFILEPATH="'${DEBUG_FILE_DIR}'"' \
        'pthi="'${PERTURB_NAMELIST}'"'

    return 9
}

# Function to run Python test
run_python_test() {
    # Generate initial conditions
    python ${BETACAST}/py_atm_to_cam/atm_to_cam.py \
        --datasource "ERA5RDA" \
        --numlevels 10 \
        --YYYYMMDDHH 2005082800 \
        --data_filename "${TEST_FILES_DIR}/ds633.0/e5.oper.invariant/197901/e5.oper.invariant.128_129_z.ll025sc.1979010100_1979010100.nc" \
        --wgt_filename "${TEST_FILES_DIR}/map_era5_0.25x0.25_TO_ne16_patc.nc" \
        --dycore "se" \
        --RDADIR "${TEST_FILES_DIR}/ds633.0/" \
        --mpas_as_cam \
        --compress_file \
        --write_floats \
        --add_cloud_vars \
        --adjust_config "a" \
        --model_topo_file "${TEST_FILES_DIR}/ne16np4_nc3000_Co120_Fi001_PF_nullRR_Nsw084_20171012.nc" \
        --se_inic "${PY_BEFORE_ATM}" \
        --verbose \
        --write_debug_files \
        --write_debug_dir "${DEBUG_FILE_DIR}"

    # Add perturbations
    python ${BETACAST}/py_atm_to_cam/perturb/add_perturbations_to_cam.py \
        --BEFOREPERTFILE ${PY_BEFORE_ATM} \
        --AFTERPERTFILE "${PY_BEFORE_ATM}_PERT.nc" \
        --gridfile ${MODEL_GRID} \
        --MAPFILEPATH ${DEBUG_FILE_DIR} \
        --pthi ${PERTURB_NAMELIST}
}

# Function to validate results
run_validation() {
    # First compare the initial conditions
    echo "Comparing initial conditions..."
    python ${BETACAST}/py_testing/check_same.py ${PY_BEFORE_ATM} ${NCL_BEFORE_ATM} || return 1

    # Then compare the perturbed files
    echo "Comparing perturbed files..."
    python ${BETACAST}/py_testing/check_same.py "${PY_BEFORE_ATM}_PERT.nc" "${NCL_BEFORE_ATM}_PERT.nc" || return 1

    return 0
}

# Optional cleanup function
cleanup() {
    rm -f ${PY_BEFORE_ATM} ${NCL_BEFORE_ATM}
    rm -f ${PY_BEFORE_ATM}_PERT.nc ${NCL_BEFORE_ATM}_PERT.nc
}