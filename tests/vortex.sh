#!/bin/bash

# Set vortex namelist path
export VORTEX_NAMELIST="/glade/u/home/zarzycki/scratch/betacast/namelists/unseed.default.nl"

# Function to run Python test
run_python_test() {
    # Initial ERA5 conversion
    python ${BETACAST}/atm_to_cam/atm_to_cam.py \
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
    python ${BETACAST}/atm_to_cam/tcseed/find-tc-fill-params.py \
        --inic_file "${DEBUG_FILE_DIR}/py_final.nc" \
        --vortex_namelist ${VORTEX_NAMELIST}

    # Seed TC in data
    python ${BETACAST}/atm_to_cam/tcseed/py-seed-tc-in-ncdata.py \
        --se_inic "${DEBUG_FILE_DIR}/py_final.nc" \
        --vortex_namelist ${VORTEX_NAMELIST}
}

# Optional cleanup function
cleanup() {
    rm -f ${DEBUG_FILE_DIR}/py_final.nc
}
