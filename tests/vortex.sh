#!/bin/bash
# Environment setup
export TEST_FILES_DIR="/glade/u/home/zarzycki/scratch/test_files/"
export DEBUG_FILE_DIR="/glade/u/home/zarzycki/scratch/tmp_betacast/"
export BETACAST="/glade/u/home/zarzycki/scratch/betacast/"

# Set vortex namelist path
export VORTEX_NAMELIST="/glade/u/home/zarzycki/scratch/betacast/namelists/unseed.default.nl"

# Function to run Python test
run_python_test() {
    # Initial ERA5 conversion
    python ${BETACAST}/py_atm_to_cam/atm_to_cam.py \
        --datasource "ERA5RDA" \
        --numlevels 10 \
        --YYYYMMDDHH 2005082800 \
        --data_filename "${TEST_FILES_DIR}/ds633.0/e5.oper.invariant/197901/e5.oper.invariant.128_129_z.ll025sc.1979010100_1979010100.nc" \
        --wgt_filename "${TEST_FILES_DIR}/map_era5_0.25x0.25_TO_ne120np4_SCRIP_patc.nc" \
        --dycore "se" \
        --RDADIR "${TEST_FILES_DIR}/ds633.0/" \
        --mpas_as_cam \
        --compress_file \
        --write_floats \
        --add_cloud_vars \
        --adjust_config "a" \
        --model_topo_file "${TEST_FILES_DIR}/ne120np4_nc3000_Co015_Fi001_PF_nullRR_Nsw010_20171011.nc" \
        --se_inic "${DEBUG_FILE_DIR}/py_final.nc" \
        --verbose \
        --write_debug_files \
        --write_debug_dir "${DEBUG_FILE_DIR}"

    # Find TC parameters
    python ${BETACAST}/py_atm_to_cam/tcseed/find-tc-fill-params.py \
        --inic_file "${DEBUG_FILE_DIR}/py_final.nc" \
        --vortex_namelist ${VORTEX_NAMELIST}

    # Seed TC in data
    python ${BETACAST}/py_atm_to_cam/tcseed/py-seed-tc-in-ncdata.py \
        --se_inic "${DEBUG_FILE_DIR}/py_final.nc" \
        --vortex_namelist ${VORTEX_NAMELIST}
}

# Function to run NCL test
run_ncl_test() {
    # Initial ERA5 conversion
    ncl -n ${BETACAST}/atm_to_cam/atm_to_cam.ncl \
        'datasource="ERA5RDA"' \
        numlevels=10 \
        YYYYMMDDHH=2005082800 \
        'dycore="se"' \
        'data_filename="'${TEST_FILES_DIR}'/ds633.0/e5.oper.invariant/197901/e5.oper.invariant.128_129_z.ll025sc.1979010100_1979010100.nc"' \
        'wgt_filename="'${TEST_FILES_DIR}'/map_era5_0.25x0.25_TO_ne120np4_SCRIP_patc.nc"' \
        'RDADIR="'${TEST_FILES_DIR}'/ds633.0/"' \
        mpas_as_cam=True \
        'model_topo_file="'${TEST_FILES_DIR}'/ne120np4_nc3000_Co015_Fi001_PF_nullRR_Nsw010_20171011.nc"' \
        'adjust_config="a"' \
        compress_file=True \
        write_floats=True \
        add_cloud_vars=True \
        write_debug_files=True \
        'write_debug_dir="'${DEBUG_FILE_DIR}'"' \
        'se_inic="'${DEBUG_FILE_DIR}'/ncl_final.nc"'

    # Find TC parameters
    ncl -n ${BETACAST}/atm_to_cam/tcseed/find-tc-fill-params.ncl \
        'inic_file="'${DEBUG_FILE_DIR}'/ncl_final.nc"' \
        'pthi="'${VORTEX_NAMELIST}'"'

    # Seed TC in data
    ncl -n ${BETACAST}/atm_to_cam/tcseed/seed-tc-in-ncdata.ncl \
        'seedfile="'${DEBUG_FILE_DIR}'/ncl_final.nc"' \
        'pthi="'${VORTEX_NAMELIST}'"'

    return 9
}

# Function to validate results
run_validation() {
    python ${BETACAST}/py_testing/check_same.py ${DEBUG_FILE_DIR}/py_final.nc ${DEBUG_FILE_DIR}/ncl_final.nc
}

# Optional cleanup function
cleanup() {
    rm -f ${DEBUG_FILE_DIR}/py_final.nc ${DEBUG_FILE_DIR}/ncl_final.nc
}
