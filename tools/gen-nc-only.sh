#!/bin/bash

# Environment variables
export BETACAST=/glade/u/home/zarzycki/scratch/betacast/
TEST_FILES_DIR=/glade/campaign/univ/upsu0032/betacast_test_files/
DEBUG_FILE_DIR=/glade/u/home/zarzycki/scratch/tmp_betacast/
mapping_files_path=$PWD

# For E3SM, np is the gridfile, pg is the topofile
modelgridfile="/glade/u/home/zarzycki/work/CESM_files/grids/scrip/SnowEater_WUS_30x8_np4_scrip.nc"
MODEL_TOPO_FILE="/glade/u/home/zarzycki/work/CESM_files/topo/USGS-gtopo30_SnowEater_WUS_30x8np4pg2_x6t.nc"
SE_INIC="/glade/derecho/scratch/zarzycki/SnowEater_WUS_30x8np4_L32_inic.nc"

# Model configuration
DYCORE="se"
ADJUST_CONFIG="a"
NUM_LEVELS=32
YYYYMMDDHH=2005082900

# Source data configuration
DATASOURCE="SAMPLE"
RLLSOURCEGRID="era5_2deg"
DATA_FILENAME="${BETACAST}/grids/samples/ERA5_2005082900_strided.nc"
RDA_DIR=""

# Python paths and directories
PY_REMAPPING_PATH=${BETACAST}/py_remapping
GRIDS_PATH=${BETACAST}/grids/anl_scrip
ATM_TO_CAM_PATH=${BETACAST}/py_atm_to_cam

# Change to working directory
cd $ATM_TO_CAM_PATH
echo "cd'ing to interpolation directory: $ATM_TO_CAM_PATH"

# Use existing weights file if valid, otherwise generate from modelgridfile

if [[ -n "${anl2mdlWeights}" && -e "${anl2mdlWeights}" ]]; then
  echo "User has provided anl2mdlWeights --> ${anl2mdlWeights}, using those!"
  write_weights="false"
else
  echo "User did not explicitly specify anl2mdlWeights, trying to generate from SCRIP grid"

  # Ensure modelgridfile is defined and exists
  if [[ -z "${modelgridfile}" || ! -f "${modelgridfile}" ]]; then
    echo "modelgridfile is not defined or file --> ${modelgridfile} does not exist, exiting"
    echo "specify this as a SCRIP file in the namelist or anl2mdlWeights"
    exit 19
  fi

  # Prepare for weight file generation
  modelgridshortname=$(basename "${modelgridfile%.*}")
  anl2mdlWeights=${mapping_files_path}/map_${RLLSOURCEGRID}_TO_${modelgridshortname}_patc.nc

  if [ ! -f ${anl2mdlWeights} ]; then
    echo "Writing anl2mdlWeights --> ${anl2mdlWeights}"
    write_weights="true"
    python ${PY_REMAPPING_PATH}/gen_analysis_to_model_wgt_file.py \
        --ANLGRID "${RLLSOURCEGRID}" \
        --DSTGRIDNAME "${modelgridshortname}" \
        --DSTGRIDFILE "${modelgridfile}" \
        --ANLGRIDPATH "${GRIDS_PATH}" \
        --WGTFILEDIR "${mapping_files_path}"
  else
    echo "Betacast-generated anl2mdlWeights --> ${anl2mdlWeights} already exists, using those!"
    write_weights="false"
  fi
fi

# Run main processing with Python
time python ${ATM_TO_CAM_PATH}/atm_to_cam.py \
    --datasource "${DATASOURCE}" \
    --numlevels ${NUM_LEVELS} \
    --YYYYMMDDHH ${YYYYMMDDHH} \
    --data_filename "${DATA_FILENAME}" \
    --wgt_filename "${anl2mdlWeights}" \
    --dycore "${DYCORE}" \
    --RDADIR "${RDA_DIR}" \
    --mpas_as_cam \
    --compress_file \
    --write_floats \
    --add_cloud_vars \
    --adjust_config "${ADJUST_CONFIG}" \
    --model_topo_file "${MODEL_TOPO_FILE}" \
    --se_inic "${SE_INIC}" \
    --verbose \
    --write_debug_files \
    --write_debug_dir "${DEBUG_FILE_DIR}"

if [ "$write_weights" = "true" ]; then
  echo "Removing temporary weights..."
  rm -fv "$anl2mdlWeights"
fi
