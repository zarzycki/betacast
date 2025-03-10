#!/bin/bash
# Environment setup
export TEST_FILES_DIR="/glade/u/home/zarzycki/scratch/test_files/"
export DEBUG_FILE_DIR="/glade/u/home/zarzycki/scratch/tmp_betacast/"
export BETACAST="/glade/u/home/zarzycki/scratch/betacast/"

# Additional variables
DATE="2023052300"
NLEV="30"

# Function to run Python test
run_python_test() {
    # First ERA5 background
    python ${BETACAST}/py_atm_to_cam/atm_to_cam.py \
        --datasource "ERA5RDA" \
        --numlevels ${NLEV} \
        --YYYYMMDDHH ${DATE} \
        --data_filename "${TEST_FILES_DIR}/ds633.0/e5.oper.invariant/197901/e5.oper.invariant.128_129_z.ll025sc.1979010100_1979010100.nc" \
        --wgt_filename "${TEST_FILES_DIR}/map_era5_0.25x0.25_TO_ne120np4_SCRIP_patc.nc" \
        --dycore "se" \
        --add_cloud_vars \
        --RDADIR "${TEST_FILES_DIR}/ds633.0/" \
        --adjust_config "a" \
        --model_topo_file "${TEST_FILES_DIR}/ne120np4_nc3000_Co015_Fi001_PF_nullRR_Nsw010_20171011.nc" \
        --mod_remap_file "" \
        --mod_in_topo "" \
        --se_inic "${DEBUG_FILE_DIR}/frank_py_final.nc" \
        --verbose \
        --write_debug_files \
        --write_debug_dir "${DEBUG_FILE_DIR}"

    # Generate SCRIP grid
    python ${BETACAST}/py_remapping/gen_reglatlon_SCRIP.py \
        --dstGridName "hwrf_storm_scrip.nc" \
        --dstDir "./" \
        --srcfilename "${TEST_FILES_DIR}/mawar02w.2023052300.hwrfprs.storm.0p015.f000.grb2"

    # Generate weight file
    python ${BETACAST}/py_remapping/gen_analysis_to_model_wgt_file.py \
        --ANLGRID "hwrf_storm" \
        --DSTGRIDNAME "modelgrid" \
        --DSTGRIDFILE "/glade/u/home/zarzycki/work/grids/scrip/ne120np4_SCRIP.nc" \
        --ANLGRIDPATH "./" \
        --WGTFILEDIR "./"

    # HWRF overlay
    python ${BETACAST}/py_atm_to_cam/atm_to_cam.py \
        --datasource "HWRF" \
        --numlevels ${NLEV} \
        --YYYYMMDDHH ${DATE} \
        --data_filename "${TEST_FILES_DIR}/mawar02w.2023052300.hwrfprs.storm.0p015.f000.grb2" \
        --wgt_filename "./map_hwrf_storm_TO_modelgrid_patc.nc" \
        --dycore "se" \
        --add_cloud_vars \
        --adjust_config "" \
        --se_inic "${DEBUG_FILE_DIR}/frank_py_reg.nc" \
        --verbose \
        --write_debug_files \
        --write_debug_dir "${DEBUG_FILE_DIR}"

    # Perform overlay
    python ${BETACAST}/py_atm_to_cam/overlay.py \
        "${DEBUG_FILE_DIR}/frank_py_final.nc" \
        "${DEBUG_FILE_DIR}/frank_py_reg.nc" \
        --maxLev 80.

    # Cleanup intermediate files
    rm -f hwrf_storm_scrip.nc map_hwrf_storm_TO_modelgrid_patc.nc PET0.RegridWeightGen.Log

    return 0
}

# Function to run NCL test
run_ncl_test() {
    # First ERA5 background
    ncl -n ${BETACAST}/atm_to_cam/atm_to_cam.ncl \
        'datasource="ERA5RDA"' \
        numlevels=${NLEV} \
        YYYYMMDDHH=${DATE} \
        'dycore="se"' \
        'data_filename="'${TEST_FILES_DIR}'/ds633.0/e5.oper.invariant/197901/e5.oper.invariant.128_129_z.ll025sc.1979010100_1979010100.nc"' \
        'wgt_filename="'${TEST_FILES_DIR}'/map_era5_0.25x0.25_TO_ne120np4_SCRIP_patc.nc"' \
        'RDADIR="'${TEST_FILES_DIR}'/ds633.0/"' \
        'model_topo_file="'${TEST_FILES_DIR}'/ne120np4_nc3000_Co015_Fi001_PF_nullRR_Nsw010_20171011.nc"' \
        'mod_remap_file=""' \
        'mod_in_topo=""' \
        'adjust_config="a"' \
        write_debug_files=True \
        'write_debug_dir="'${DEBUG_FILE_DIR}'"' \
        'se_inic="'${DEBUG_FILE_DIR}'/frank_ncl_final.nc"'

    # Generate SCRIP grid
    ncl ${BETACAST}/remapping/gen_reglatlon_SCRIP.ncl \
        'DSTGRIDNAME="hwrf_storm_scrip.nc"' \
        'DSTDIR="./"' \
        'SRCFILENAME="'${TEST_FILES_DIR}'/mawar02w.2023052300.hwrfprs.storm.0p015.f000.grb2"'

    # Generate weight file
    ncl ${BETACAST}/remapping/gen_analysis_to_model_wgt_file.ncl \
        'ANLGRID="hwrf_storm"' \
        'DSTGRIDNAME="modelgrid"' \
        'ANLGRIDPATH="./"' \
        'WGTFILEDIR="./"' \
        'DSTGRIDFILE="'${TEST_FILES_DIR}'/ne120np4_SCRIP.nc"'

    # HWRF overlay
    ncl -n ${BETACAST}/atm_to_cam/atm_to_cam.ncl \
        'datasource="HWRF"' \
        'dycore="se"' \
        numlevels=${NLEV} \
        YYYYMMDDHH=${DATE} \
        'data_filename="'${TEST_FILES_DIR}'/mawar02w.2023052300.hwrfprs.storm.0p015.f000.grb2"' \
        'wgt_filename="./map_hwrf_storm_TO_modelgrid_patc.nc"' \
        'adjust_config=""' \
        write_debug_files=True \
        'write_debug_dir="'${DEBUG_FILE_DIR}'"' \
        'se_inic="'${DEBUG_FILE_DIR}'/frank_ncl_reg.nc"'

    # Perform overlay
    python ${BETACAST}/py_atm_to_cam/overlay.py \
        "${DEBUG_FILE_DIR}/frank_ncl_final.nc" \
        "${DEBUG_FILE_DIR}/frank_ncl_reg.nc" \
        --maxLev 80.

    # Cleanup intermediate files
    rm -f hwrf_storm_scrip.nc map_hwrf_storm_TO_modelgrid_patc.nc PET0.RegridWeightGen.Log

    return 9
}

# Function to validate results
run_validation() {
    python ${BETACAST}/py_testing/check_same.py ${DEBUG_FILE_DIR}/frank_py_final.nc ${DEBUG_FILE_DIR}/frank_ncl_final.nc
}

# Optional cleanup function
cleanup() {
    rm -f ${DEBUG_FILE_DIR}/frank_py_final.nc ${DEBUG_FILE_DIR}/frank_ncl_final.nc
    rm -f ${DEBUG_FILE_DIR}/frank_py_reg.nc ${DEBUG_FILE_DIR}/frank_ncl_reg.nc
}
