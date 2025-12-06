#!/bin/bash

# export BETACAST atm_to_cam_path mapping_files_path \
#        modelgridfile m2m_gridfile model_scrip m2m_parent_source \
#        anl2mdlWeights m2m_remap_file adjust_topo m2m_topo_in adjust_flags \
#        yearstr monthstr daystr cyclestr uniqtime \
#        DO_PYTHON DYCORE atmDataType numLevels ERA5RDA \
#        do_frankengrid standalone_vortex add_noise add_perturbs \
#        modelSystem sstDataType \
#        sePreFilterIC sstFileIC perturb_namelist vortex_namelist regional_src \
#        RDADIR gfs_files_path era_files_path \
#        AUGMENT_STR VORTEX_STR

set -euo pipefail

if [ -z "${BETACAST}" ]; then
  echo "BETACAST must be set in environment" >&2
  exit 1
fi

source "${BETACAST}/utils.sh"
source "${BETACAST}/datahelpers.sh"

echo "ATM runner host: $(hostname)"
echo "ATM runner start time: $(date -u +"%Y-%m-%d %H:%M:%S UTC")"

############################### ARRAYS ###############################

declare -A atm_data_sources=(
  ["1"]="GFS"
  ["2"]="ERAI"
  ["3"]="CFSR"
  ["4"]="ERA5"
  ["9"]="CAM"
)
declare -A atm_data_glob_anl=(
  ["1"]="gfs_0.25x0.25"
  ["2"]=""
  ["3"]="gfs_0.50x0.50"
  ["4"]="era5_0.25x0.25"
  ["9"]=""
)
declare -A atm_file_paths=(
  ["1"]="${gfs_files_path}/gfs_atm_${yearstr}${monthstr}${daystr}${cyclestr}.grib2"
  ["2"]="${era_files_path}/ERA-Int_${yearstr}${monthstr}${daystr}${cyclestr}.nc"
  ["3"]="${gfs_files_path}/cfsr_atm_${yearstr}${monthstr}${daystr}${cyclestr}.grib2"
  ["4"]="${era_files_path}/ERA5_${yearstr}${monthstr}${daystr}${cyclestr}.nc"
  ["9"]=""
)

# If ERA5RDA flag toggled, set value w/ key to RDA data
if [ "$ERA5RDA" -eq 1 ] ; then
  atm_data_sources["4"]="ERA5RDA"
  atm_file_paths["4"]="${RDADIR}/e5.oper.invariant/197901/e5.oper.invariant.128_129_z.ll025sc.1979010100_1979010100.nc"
fi

############################### FLAGS ###############################

if [ "${modelSystem}" -eq 1 ]; then
  ADDCHEM_STR="--add_chemistry"
fi

############################### ATM NCL / PYTHON ###############################

echo "cd'ing to interpolation directory: $atm_to_cam_path"
cdv "$atm_to_cam_path"

# Figure out which anl2mdlWeights we want to use.
if [[ -z "${anl2mdlWeights}" || ! -e "${anl2mdlWeights}" ]]; then
  echo "User did not explicitly specify anl2mdlWeights, trying to generate from SCRIP grid"
  if [ ! -f "${modelgridfile}" ]; then
    echo "modelgridfile --> ${modelgridfile} does not exist, exiting"
    echo "specify this as a SCRIP file in the namelist or anl2mdlWeights"
    exit 19
  fi
  modelgridshortname=$(basename "${modelgridfile%.*}")
  if [[ "$atmDataType" -eq 9 ]]; then
    RLLSOURCEGRID="era5_0.25x0.25"
  else
    RLLSOURCEGRID="${atm_data_glob_anl[$atmDataType]}"
  fi
  anl2mdlWeights=${mapping_files_path}/map_${RLLSOURCEGRID}_TO_${modelgridshortname}_patc.nc

  if [ ! -f "${anl2mdlWeights}" ]; then
    echo "Writing anl2mdlWeights --> ${anl2mdlWeights}"
    if [ "$DO_PYTHON" = true ]; then
      (set -x; python ../py_remapping/gen_analysis_to_model_wgt_file.py \
        --ANLGRID "${RLLSOURCEGRID}" \
        --DSTGRIDNAME "${modelgridshortname}" \
        --DSTGRIDFILE "${modelgridfile}" \
        --ANLGRIDPATH "../grids/anl_scrip/" \
        --WGTFILEDIR "${mapping_files_path}"
      )
    else
      set +e
      (set -x; ncl ../remapping/gen_analysis_to_model_wgt_file.ncl \
        'ANLGRID="'${RLLSOURCEGRID}'"' \
        'ANLGRIDPATH="../grids/anl_scrip/"' \
        'DSTGRIDNAME="'${modelgridshortname}'"' \
        'DSTGRIDFILE="'${modelgridfile}'"' \
        'WGTFILEDIR="'${mapping_files_path}'"' \
      ) ; exit_status=$?
      check_ncl_exit "gen_analysis_to_model_wgt_file.ncl" $exit_status
      set -e
    fi
  else
    echo "Betacast-generated anl2mdlWeights --> ${anl2mdlWeights} already exists, using those!"
  fi
else
  echo "User has provided anl2mdlWeights --> ${anl2mdlWeights}, using those!"
fi

if [[ "$atmDataType" -eq 9 ]]; then
  if [[ -z "${m2m_remap_file}" || ! -e "${m2m_remap_file}" ]]; then
    echo "User did not explicitly specify m2m_remap_file, trying to generate from SCRIP grid"
    if [ ! -f "${m2m_gridfile}" ]; then
      echo "m2m_gridfile --> ${m2m_gridfile} does not exist, exiting"
      echo "specify this as a SCRIP file in the namelist or m2m_remap_file"
      exit 19
    fi
    m2mgridshortname=$(basename "${m2m_gridfile%.*}")
    m2m_remap_file=${mapping_files_path}/map_${m2mgridshortname}_TO_era5_0.25x0.25_patc.nc

    if [ ! -f "${m2m_remap_file}" ]; then
      echo "Writing m2m_remap_file --> ${m2m_remap_file}"
      set +e
      (set -x; ncl ../remapping/gen_analysis_to_model_wgt_file.ncl \
        'ANLGRID="era5_0.25x0.25"' \
        'ANLGRIDPATH="../grids/anl_scrip/"' \
        'DSTGRIDNAME="'${m2mgridshortname}'"' \
        'DSTGRIDFILE="'${m2m_gridfile}'"' \
        'WGTFILEDIR="'${mapping_files_path}'"' \
        FLIP_MODEL_AND_ANALYSIS=True \
      ) ; exit_status=$?
      check_ncl_exit "gen_analysis_to_model_wgt_file.ncl" $exit_status
      set -e
    else
      echo "Betacast-generated m2m_remap_file --> ${m2m_remap_file} already exists, using those!"
    fi
  else
    echo "User has provided m2m_remap_file --> ${m2m_remap_file}, using those!"
  fi
fi

if [[ "$atmDataType" -eq 9 ]]; then
  if [[ -d "$m2m_parent_source" ]]; then
    echo "m2m_parent_source ($m2m_parent_source) is provided as a dir of nc files."
    if [ "$DO_PYTHON" = true ]; then
      (set -x; python find-time-file.py \
        --DIR "${m2m_parent_source}" \
        --YYYYMMDDHH ${yearstr}${monthstr}${daystr}${cyclestr} \
        --UQSTR "${uniqtime}"
      )
    else
      set +e
      (set -x; ncl find-time-file.ncl \
        'DIR="'${m2m_parent_source}'"' \
        YYYYMMDDHH=${yearstr}${monthstr}${daystr}${cyclestr} \
        'UQSTR="'${uniqtime}'"' \
      ) ; exit_status=$?
      check_ncl_exit "find-time-file.ncl" $exit_status
      set -e
    fi
    while IFS= read -r line || [[ -n "$line" ]]; do
      atm_file_paths["9"]="$line"
      break
    done < "m2mfile.$uniqtime"
    [[ ! -f "${atm_file_paths["9"]}" ]] && { echo "File does not exist."; exit 1; }
  elif [[ -f "$m2m_parent_source" ]]; then
    echo "m2m_parent_source ($m2m_parent_source) is provided as a file."
    atm_file_paths["9"]="$m2m_parent_source"
  else
    echo "m2m_parent_source ($m2m_parent_source) is not a file or directory. Exiting."
    exit 1
  fi
fi

echo "Doing atm_to_cam"

if [ "$DO_PYTHON" = true ]; then
  (set -x; python atm_to_cam.py \
    --datasource "${atm_data_sources[$atmDataType]}" \
    --numlevels ${numLevels} \
    --YYYYMMDDHH ${yearstr}${monthstr}${daystr}${cyclestr} \
    --data_filename "${atm_file_paths[$atmDataType]}" \
    --wgt_filename "${anl2mdlWeights}" \
    --dycore "${DYCORE}" \
    --add_cloud_vars \
    --add_chemistry \
    --RDADIR "${RDADIR}" \
    --adjust_config "${adjust_flags-}" \
    --model_topo_file "${adjust_topo-}" \
    --mod_remap_file "${m2m_remap_file-}" \
    --mod_in_topo "${m2m_topo_in-}" \
    --se_inic "${sePreFilterIC}" \
    ${ADDCHEM_STR:+$ADDCHEM_STR} ${AUGMENT_STR:+$AUGMENT_STR} ${VORTEX_STR:+$VORTEX_STR}
  )
else
  set +e
  (set -x; ncl -n atm_to_cam.ncl \
      'datasource="'${atm_data_sources[$atmDataType]}'"' \
      numlevels=${numLevels} \
      YYYYMMDDHH=${yearstr}${monthstr}${daystr}${cyclestr} \
      'dycore="'${DYCORE}'"' \
      'data_filename="'${atm_file_paths[$atmDataType]}'"' \
      'RDADIR="'${RDADIR}'"' \
      'wgt_filename="'${anl2mdlWeights}'"' \
      'model_topo_file="'${adjust_topo-}'"' \
      'mod_remap_file="'${m2m_remap_file-}'"' \
      'mod_in_topo="'${m2m_topo_in-}'"' \
      'adjust_config="'${adjust_flags-}'"' \
      'se_inic = "'${sePreFilterIC}'"'
  ) ; exit_status=$?
  check_ncl_exit "atm_to_cam.ncl" $exit_status
  set -e
fi

############################### FRANKENGRID ###############################

if [ "${do_frankengrid}" = true ] ; then
  regional_src=${regional_src/YYYY/$yearstr}
  regional_src=${regional_src/MM/$monthstr}
  regional_src=${regional_src/DD/$daystr}
  regional_src=${regional_src/HH/$cyclestr}

  echo "Doing Frankengrid with $regional_src"

  TMPWGTFILE="./map_hwrf_storm_TO_modelgrid_patc.nc"

  if [ "$DO_PYTHON" = true ]; then
    (set -x; python ../py_remapping/gen_reglatlon_SCRIP.py \
      --dstGridName "hwrf_storm_scrip.nc" \
      --dstDir "./" \
      --srcfilename "${regional_src}"
    )
  else
    set +e
    echo "Generating a temporary SCRIP file for HWRF"
    (set -x; ncl ../remapping/gen_reglatlon_SCRIP.ncl \
      'DSTGRIDNAME="hwrf_storm_scrip.nc"' \
      'DSTDIR="./"' \
      'SRCFILENAME="'${regional_src}'"'
    ) ; exit_status=$?
    check_ncl_exit "gen_reglatlon_SCRIP.ncl" $exit_status
  fi

  if [ "$DO_PYTHON" = true ]; then
    (set -x; python ../py_remapping/gen_analysis_to_model_wgt_file.py \
      --ANLGRID "hwrf_storm" \
      --DSTGRIDNAME "modelgrid" \
      --DSTGRIDFILE "${model_scrip}" \
      --ANLGRIDPATH "./" \
      --WGTFILEDIR "./"
    )
  else
    echo "Generating a temporary map file for HWRF"
    (set -x; ncl ../remapping/gen_analysis_to_model_wgt_file.ncl \
      'ANLGRID="hwrf_storm"' \
      'DSTGRIDNAME="modelgrid"' \
      'ANLGRIDPATH="./"' \
      'WGTFILEDIR="./"' \
      'DSTGRIDFILE="'${model_scrip}'"'
    ) ; exit_status=$?
    check_ncl_exit "gen_analysis_to_model_wgt_file.ncl" $exit_status
  fi

  if [ "$DO_PYTHON" = true ]; then
    (set -x; python atm_to_cam.py \
      --datasource "HWRF" \
      --numlevels ${numLevels} \
      --YYYYMMDDHH ${yearstr}${monthstr}${daystr}${cyclestr} \
      --data_filename "${regional_src}" \
      --wgt_filename "${TMPWGTFILE}" \
      --dycore "${DYCORE}" \
      --add_cloud_vars \
      --adjust_config "" \
      --se_inic "${sePreFilterIC}_reg.nc"
    )
  else
    echo "Generating regional Frankengrid for HWRF"
    (set -x; ncl -n atm_to_cam.ncl 'datasource="HWRF"' \
      'dycore="'${DYCORE}'"' \
      numlevels=${numLevels} \
      YYYYMMDDHH=${yearstr}${monthstr}${daystr}${cyclestr} \
      'data_filename = "'${regional_src}'"' \
      'wgt_filename="'${TMPWGTFILE}'"' \
      'adjust_config=""' \
      'se_inic = "'${sePreFilterIC}_reg.nc'"'
    ) ; exit_status=$?
    check_ncl_exit "atm_to_cam.ncl" $exit_status
    set -e
  fi

  echo "Overlay regional file on top of basefile"
  cp -v ${sePreFilterIC} ${sePreFilterIC}_base.nc
  (set -x; python overlay.py "${sePreFilterIC}" "${sePreFilterIC}_reg.nc" --maxLev 80. )

  echo "Cleaning up temporary ESMF files"
  rm -v "$TMPWGTFILE"
  rm -v hwrf_storm_scrip.nc
  rm -v ${sePreFilterIC}_reg.nc
fi

############################### VORTEX / NOISE / PERTURBS ###############################

if [ "${standalone_vortex}" = true ] ; then
  cdv "$atm_to_cam_path/tcseed"
  set +e
  echo "Adding or removing a TC from initial condition based on ${vortex_namelist}"

  echo "... finding fill parameters"
  if [ "$DO_PYTHON" = true ]; then
    (set -x; python find-tc-fill-params.py \
        --inic_file "${sePreFilterIC}" \
        --vortex_namelist ${vortex_namelist}
    )
  else
    (set -x; ncl -n find-tc-fill-params.ncl 'inic_file= "'${sePreFilterIC}'"' 'pthi = "'${vortex_namelist}'"' ) ; exit_status=$?
    check_ncl_exit "find-tc-fill-params.ncl" $exit_status
  fi

  echo "... seeding or unseeding TC"
  if [ "$DO_PYTHON" = true ]; then
    (set -x; python py-seed-tc-in-ncdata.py \
        --inic_file "${sePreFilterIC}" \
        --vortex_namelist ${vortex_namelist}
    )
  else
    (set -x; ncl -n seed-tc-in-ncdata.ncl 'seedfile = "'${sePreFilterIC}'"' 'pthi = "'${vortex_namelist}'"' ) ; exit_status=$?
    check_ncl_exit "seed-tc-in-ncdata.ncl" $exit_status
  fi
  set -e
fi

if [ "${add_noise}" = true ] ; then
  set +e
  echo "Adding white noise to initial condition"
  cdv "$atm_to_cam_path"
  (set -x; ncl -n perturb_white_noise.ncl 'basFileName = "'${sePreFilterIC}'"' ) ; exit_status=$?
  check_ncl_exit "perturb_white_noise.ncl" $exit_status
  set -e
fi

if [ "${add_perturbs}" = true ] ; then
  echo "Adding perturbations"
  cdv "$atm_to_cam_path/perturb"
  set +e

  sstFileIC_WPERT=${sstFileIC}_PERT.nc
  (set -x; ncl -n add_perturbations_to_sst.ncl 'BEFOREPERTFILE="'${sstFileIC}'"' \
     'AFTERPERTFILE = "'${sstFileIC_WPERT}'"' \
     'pthi="'${perturb_namelist}'"'
  ) ; exit_status=$?
  check_ncl_exit "add_perturbations_to_sst.ncl" $exit_status
  echo "SST perturbations added successfully"

  sePreFilterIC_WPERT=${sePreFilterIC}_PERT.nc
  (set -x; ncl -n add_perturbations_to_cam.ncl 'BEFOREPERTFILE="'${sePreFilterIC}'"'  \
     'AFTERPERTFILE = "'${sePreFilterIC_WPERT}'"' \
     'gridfile = "'${modelgridfile}'"' \
     'MAPFILEPATH = "'${mapping_files_path}'"' \
     'pthi="'${perturb_namelist}'"'
  ) ; exit_status=$?
  check_ncl_exit "add_perturbations_to_cam.ncl" $exit_status
  echo "ATM NCL completed successfully"

  set -e
  mv ${sstFileIC_WPERT} ${sstFileIC}
  mv ${sePreFilterIC_WPERT} ${sePreFilterIC}
fi

################################# E3SM LINOZ CHEM ###################################

if [ "${modelSystem}" -eq 1 ]; then
  echo "Adding chemistry variables to E3SMv3"
  check_bash_dependency ncap2 "E3SM chemistry addition"
  # Set Linoz chem variables for E3SMv3
  # H2OLNZ is just mass mixing ratio of water, ~Q
  # CH4LNZ assumes 1.8 ppm in free trop converted to MMR
  # N2OLNZ assumes 300 ppbv in free trop converted to MMR
  # NOYLNZ was "derived" from an existing input file with different values below and above 50mb
  # Note, other vars use Q*0 to create var shape/dims
  ncap2 -O \
    -s 'H2OLNZ=Q' \
    -s 'CH4LNZ=(Q*0)+1e-6' \
    -s 'N2OLNZ=(Q*0)+4.5e-7' \
    -s 'NOYLNZ=(Q*0)+6.0e-11' \
    -s 'where(lev<50) NOYLNZ=5.0e-09' \
    "${sePreFilterIC}" "${sePreFilterIC}_LINOZ.nc"
  mv -v "${sePreFilterIC}_LINOZ.nc" "${sePreFilterIC}"
  echo "... done adding chemistry variables to E3SMv3"
fi

############################### SCREAM CDF5 CONVERSION ###############################

if [ "${modelSystem}" -eq 2 ]; then
  echo "SCREAM, converting to CDF5"
  check_bash_dependency nccopy "SCREAM CDF5 conversion"
  if [[ "${sstDataType}" -ne 9 ]]; then
    timer nccopy_convert 5 "${sstFileIC}"
  fi
  timer nccopy_convert 5 "${sePreFilterIC}"
  echo "... done with converting to CDF5"
fi

echo "ATM runner complete."