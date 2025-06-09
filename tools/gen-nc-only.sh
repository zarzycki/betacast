#!/bin/bash

# Environment variables
export BETACAST=/glade/u/home/zarzycki/betacast/
mapping_files_path=$PWD

# For CESM, np is the gridfile and topofile
# For E3SM, np is the gridfile, pg is the topofile
modelgridfile="/glade/u/home/zarzycki/work/grids/scrip/ne0np4natlanticref.ne30x4.g_scrip.nc"
MODEL_TOPO_FILE="/glade/work/zarzycki/CESM_files/topo/ne0np4natlanticref.ne30x4_np4_gmted2010_modis_bedmachine_nc3000_Laplace0100_noleak_20250515.nc"
SE_INIC="/glade/work/zarzycki/FG_ne0np4natlanticref.ne30x4_np4_L32_inic.nc"
NUM_LEVELS=32

# Model configuration
DYCORE="se"
ADJUST_CONFIG="a"
YYYYMMDDHH=2017010700

# Source data configuration
# DATASOURCE="SAMPLE"
# RLLSOURCEGRID="era5_2deg"
# RDA_DIR=""
# DATA_FILENAME="${BETACAST}/grids/samples/ERA5_2deg_L21_sample.nc"
DATASOURCE="ERA5RDA"
RLLSOURCEGRID="era5_0.25x0.25"
RDA_DIR="/glade/u/home/zarzycki/rda/ds633.0/"
#RDA_DIR="/global/cfs/projectdirs/m3522/cmip6/ERA5/"
DATA_FILENAME="${RDA_DIR}/e5.oper.invariant/197901/e5.oper.invariant.128_129_z.ll025sc.1979010100_1979010100.nc"

do_frankengrid=true
REGDATASRC="HRRR"
regionalName="hrrr_3km"
regional_src="/glade/u/home/zarzycki/work/regional/HRRR_wrfprs_f00_2017010700.grib2"
DEBUG_FILE_DIR="/glade/work/zarzycki/"

#----------------------------------------------------------------------------------------

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
    --verbose

if [ "$write_weights" = "true" ]; then
  echo "Removing temporary weights..."
  rm -fv "$anl2mdlWeights"
fi

echo "Done with primary generation of: ${SE_INIC}"

if [ "$do_frankengrid" = true ] ; then

  if [ "$REGDATASRC" == "HRRR" ]; then
    regionalName="hrrr_3km"
  elif [ "$REGDATASRC" == "HWRF" ]; then
    regionalName="hwrf"
  elif [ "$REGDATASRC" == "RAP" ]; then
    regionalName="rap_13km"
  else
    echo "Error: Unknown REGDATASRC value '$REGDATASRC'" >&2
    exit 1
  fi

  echo "Doing Frankengrid $REGDATASRC/$regionalName with $regional_src"

  if [ "$REGDATASRC" == "HWRF" ]; then
    (set -x ; time python ${PY_REMAPPING_PATH}/py_remapping/gen_reglatlon_SCRIP.py \
        --dstGridName "${regionalName}_storm_scrip.nc" \
        --dstDir "${mapping_files_path}" \
        --srcfilename "${regional_src}"
    )
    TMPWGTFILE="${mapping_files_path}/map_hwrf_storm_TO_modelgrid_patc.nc"
  else
    TMPWGTFILE="${mapping_files_path}/map_${regionalName}_TO_modelgrid_patc.nc"
  fi

  (set -x ; time python ${PY_REMAPPING_PATH}/gen_analysis_to_model_wgt_file.py \
      --ANLGRID "${regionalName}" \
      --DSTGRIDNAME "modelgrid" \
      --DSTGRIDFILE "${modelgridfile}" \
      --ANLGRIDPATH "${GRIDS_PATH}" \
      --WGTFILEDIR "${mapping_files_path}"
  )

  (set -x ; time python ${ATM_TO_CAM_PATH}/atm_to_cam.py \
      --datasource "${REGDATASRC}" \
      --numlevels ${NUM_LEVELS} \
      --YYYYMMDDHH ${YYYYMMDDHH} \
      --data_filename "${regional_src}" \
      --wgt_filename "${TMPWGTFILE}" \
      --dycore "${DYCORE}" \
      --compress_file \
      --write_floats \
      --add_cloud_vars \
      --adjust_config "${ADJUST_CONFIG}" \
      --model_topo_file "${MODEL_TOPO_FILE}" \
      --se_inic "${SE_INIC}_reg.nc" \
      --verbose \
      --write_debug_files \
      --write_debug_dir "${DEBUG_FILE_DIR}"
  )

  echo "Overlay regional file on top of basefile"
  # Make a copy to archive
  cp -v ${SE_INIC} ${SE_INIC}_base.nc
  # Overlay the reg file to the OG base file
  (set -x ; time python ${ATM_TO_CAM_PATH}/overlay.py \
      "${SE_INIC}" \
      "${SE_INIC}_reg.nc" \
      --maxLev 80.
  )

  echo "Cleaning up temporary ESMF files"
  rm -fv "$TMPWGTFILE"
  rm -fv hwrf_storm_scrip.nc
  rm -fv ${SE_INIC}_reg.nc

fi
