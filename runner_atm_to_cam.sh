#!/bin/bash

# export BETACAST atm_to_cam_path mapping_files_path \
#        modelgridfile m2m_gridfile model_scrip m2m_parent_source \
#        anl2mdlWeights m2m_remap_file adjust_topo m2m_topo_in adjust_flags \
#        yearstr monthstr daystr cyclestr uniqtime \
#        DYCORE atmDataType numLevels ERA5RDA \
#        cr20v3_mean_dir cr20v3_member cr20v3_member_dir cr20v3_orog_file \
#        cr20v3_blend cr20v3_blend_lev cr20v3_blend_taper \
#        do_frankengrid standalone_vortex add_noise add_perturbs \
#        modelSystem sstDataType \
#        sePreFilterIC sstFileIC perturb_namelist vortex_namelist regional_name regional_src regional_fail_if_missing \
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
  ["5"]="CR20V3"
  ["9"]="CAM"
)
declare -A atm_data_glob_anl=(
  ["1"]="gfs_0.25x0.25"
  ["2"]=""
  ["3"]="gfs_0.50x0.50"
  ["4"]="era5_0.25x0.25"
  ["5"]="cr20v3_0.70x0.70"
  ["9"]=""
)
declare -A atm_file_paths=(
  ["1"]="${gfs_files_path}/gfs_atm_${yearstr}${monthstr}${daystr}${cyclestr}.grib2"
  ["2"]="${era_files_path}/ERA-Int_${yearstr}${monthstr}${daystr}${cyclestr}.nc"
  ["3"]="${gfs_files_path}/cfsr_atm_${yearstr}${monthstr}${daystr}${cyclestr}.grib2"
  ["4"]="${era_files_path}/ERA5_${yearstr}${monthstr}${daystr}${cyclestr}.nc"
  ["5"]="${cr20v3_orog_file}"
  ["9"]=""
)

# If ERA5RDA flag toggled, set value w/ key to RDA data
if [ "$ERA5RDA" -eq 1 ] ; then
  atm_data_sources["4"]="ERA5RDA"
  atm_file_paths["4"]="${RDADIR}/e5.oper.invariant/197901/e5.oper.invariant.128_129_z.ll025sc.1979010100_1979010100.nc"
fi

# Most datasources read their RDA-style data out of RDADIR as-is. CR20V3 is the
# exception: the mean and the per-member data live in separate archives, so point
# RDADIR at whichever one we are initializing from. The blend below goes back to
# cr20v3_mean_dir directly when it needs the mean.
if [[ "$atmDataType" -eq 5 ]]; then
  if [ -n "${cr20v3_member}" ] ; then
    # We are requesting member code (mean to come later if blended)
    RDADIR="${cr20v3_member_dir}"
    atm_data_sources["5"]="CR20V3-${cr20v3_member}"
    if [[ -z "${cr20v3_orog_file}" || ! -f "${cr20v3_orog_file}" ]]; then
      # Orography ships with the mean data, but fall back to the member archive
      # in case the user keeps everything under one root
      if [[ -f "${cr20v3_mean_dir}/invariants/surface_height.nc" ]]; then
        atm_file_paths["5"]="${cr20v3_mean_dir}/invariants/surface_height.nc"
      else
        atm_file_paths["5"]="${cr20v3_member_dir}/invariants/surface_height.nc"
      fi
    fi
  else
    # We are requesting mean only
    RDADIR="${cr20v3_mean_dir}"
  fi
  if [[ ! -f "${atm_file_paths[5]}" ]]; then
    echo "Orography file for CR20V3 (${atm_file_paths[5]}) doesn't exist, exiting!"
    exit 1
  fi
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
    (set -x; python ../remapping/gen_analysis_to_model_wgt_file.py \
      --ANLGRID "${RLLSOURCEGRID}" \
      --DSTGRIDNAME "${modelgridshortname}" \
      --DSTGRIDFILE "${modelgridfile}" \
      --ANLGRIDPATH "../grids/anl_scrip/" \
      --WGTFILEDIR "${mapping_files_path}"
    )
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
      (set -x; python ../remapping/gen_analysis_to_model_wgt_file.py \
        --ANLGRID "era5_0.25x0.25" \
        --ANLGRIDPATH "../grids/anl_scrip/" \
        --DSTGRIDNAME "${m2mgridshortname}" \
        --DSTGRIDFILE "${m2m_gridfile}" \
        --WGTFILEDIR "${mapping_files_path}" \
        --FLIP_MODEL_AND_ANALYSIS
      )
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
    (set -x; python find-time-file.py \
      --DIR "${m2m_parent_source}" \
      --YYYYMMDDHH ${yearstr}${monthstr}${daystr}${cyclestr} \
      --UQSTR "${uniqtime}"
    )
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

##################################### VERTICAL BLEND ###############################

# CR20V3 per-member files stop at 200hPa, so extrapolating above that is unreliable.
# Instead, build a second IC from the ensemble mean (which has full upper-level data)
# and splice its upper levels onto the member file. Conceptually the same idea as
# Frankengrid, but blending in the vertical rather than overlaying in the horizontal.
if [[ "$atmDataType" -eq 5 && -n "${cr20v3_member}" && "${cr20v3_blend}" = true ]] ; then

  echo "Blending CR20V3 mean above ${cr20v3_blend_lev} hPa onto member ${cr20v3_member}"

  CR20MEANFILE="${sePreFilterIC}_${yearstr}${monthstr}${daystr}${cyclestr}_cr20mean.nc"

  # Mean is on the same 0.70x0.70 grid as the members, so anl2mdlWeights is reused
  (set -x; time python atm_to_cam.py \
    --datasource "CR20V3" \
    --numlevels ${numLevels} \
    --YYYYMMDDHH ${yearstr}${monthstr}${daystr}${cyclestr} \
    --data_filename "${atm_file_paths[5]}" \
    --wgt_filename "${anl2mdlWeights}" \
    --dycore "${DYCORE}" \
    --add_cloud_vars \
    --RDADIR "${cr20v3_mean_dir}" \
    --adjust_config "${adjust_flags-}" \
    --model_topo_file "${adjust_topo-}" \
    --se_inic "${CR20MEANFILE}"
  )

  # Updates sePreFilterIC in place with mean conditions aloft, CR20MEANFILE is read-only
  (set -x; time python vertical_blend.py \
      "${sePreFilterIC}" \
      "${CR20MEANFILE}" \
      --blendLev ${cr20v3_blend_lev} \
      --taperScale ${cr20v3_blend_taper}
  )

  echo "Cleaning up temporary CR20V3 mean file"
  rm -fv "${CR20MEANFILE}"
fi

############################### FRANKENGRID ###############################

if [ "${do_frankengrid}" = true ] ; then
  regional_src=${regional_src/YYYY/$yearstr}
  regional_src=${regional_src/MM/$monthstr}
  regional_src=${regional_src/DD/$daystr}
  regional_src=${regional_src/HH/$cyclestr}

  echo "Doing Frankengrid $regional_name with $regional_src"

  # Check to see if regional file exists and decide what to do if it doesn't.
  if [ ! -f "${regional_src}" ]; then
    if [ "${regional_fail_if_missing}" = true ]; then
      echo "!!! ERROR: regional_src file '${regional_src}' does not exist and regional_fail_if_missing=true, exiting !!!" >&2
      exit 1
    else
      echo "!!! WARNING: regional_src file '${regional_src}' does not exist and regional_fail_if_missing=false, skipping Frankengrid !!!" >&2
    fi
  else
    case "$regional_name" in
      hwrf)        regFrankenDataSource="HWRF"   ;;
      hrrr_3km)    regFrankenDataSource="HRRR"   ;;
      hrrr_3km_ml) regFrankenDataSource="HRRRml" ;;
      rap_13km)    regFrankenDataSource="RAP"    ;;
      *) echo "Unknown regional_name '${regional_name}' for do_frankengrid, exiting"; exit 1 ;;
    esac

    # hrrr_3km_ml (native/hybrid levels) shares its horizontal grid with hrrr_3km
    # (pressure levels), so it reuses the same SCRIP file and weights.
    case "$regional_name" in
      hrrr_3km_ml) REG_ANLGRID="hrrr_3km" ;;
      *)           REG_ANLGRID="${regional_name}" ;;
    esac

    if [ "$regional_name" == "hwrf" ]; then
      (set -x; time python ../remapping/gen_reglatlon_SCRIP.py \
          --dstGridName "hwrf_storm_scrip.nc" \
          --dstDir "${mapping_files_path}" \
          --srcfilename "${regional_src}"
      )
      TMPWGTFILE="${mapping_files_path}/map_hwrf_storm_TO_modelgrid_patc.nc"
      REG_ANLGRID="hwrf_storm"
      REG_ANLGRIDPATH="${mapping_files_path}"
    else
      TMPWGTFILE="${mapping_files_path}/map_${REG_ANLGRID}_TO_modelgrid_patc.nc"
      REG_ANLGRIDPATH="../grids/anl_scrip/"
    fi

    (set -x; time python ../remapping/gen_analysis_to_model_wgt_file.py \
        --ANLGRID "${REG_ANLGRID}" \
        --DSTGRIDNAME "modelgrid" \
        --DSTGRIDFILE "${modelgridfile}" \
        --ANLGRIDPATH "${REG_ANLGRIDPATH}" \
        --WGTFILEDIR "${mapping_files_path}"
    )

    (set -x; time python atm_to_cam.py \
        --datasource "${regFrankenDataSource}" \
        --numlevels ${numLevels} \
        --YYYYMMDDHH ${yearstr}${monthstr}${daystr}${cyclestr} \
        --data_filename "${regional_src}" \
        --wgt_filename "${TMPWGTFILE}" \
        --dycore "${DYCORE}" \
        --compress_file \
        --write_floats \
        --add_cloud_vars \
        --adjust_config "${adjust_flags-}" \
        --model_topo_file "${adjust_topo-}" \
        --se_inic "${sePreFilterIC}_reg.nc"
    )

    echo "Overlay regional file on top of basefile"
    # Make a copy to archive
    cp -v ${sePreFilterIC} ${sePreFilterIC}_base.nc
    # Overlay the reg file to the OG base file
    (set -x; time python overlay.py \
        "${sePreFilterIC}" \
        "${sePreFilterIC}_reg.nc" \
        --maxLev 80.
    )

    echo "Cleaning up temporary ESMF files"
    rm -fv "$TMPWGTFILE"
    rm -fv "${mapping_files_path}/hwrf_storm_scrip.nc"
    rm -fv ${sePreFilterIC}_reg.nc
  fi
fi

############################### VORTEX / NOISE / PERTURBS ###############################

if [ "${standalone_vortex}" = true ] ; then
  cdv "$atm_to_cam_path/tcseed"
  set +e
  echo "Adding or removing a TC from initial condition based on ${vortex_namelist}"

  echo "... finding fill parameters"

  (set -x; python find-tc-fill-params.py \
      --inic_file "${sePreFilterIC}" \
      --vortex_namelist ${vortex_namelist}
  ) ; exit_status=$?
  check_python_exit "find-tc-fill-params.py" $exit_status

  echo "... seeding or unseeding TC"

  (set -x; python py-seed-tc-in-ncdata.py \
      --se_inic "${sePreFilterIC}" \
      --vortex_namelist ${vortex_namelist}
  ) ; exit_status=$?
  check_python_exit "py-seed-tc-in-ncdata.py" $exit_status

  set -e
fi

if [ "${add_noise}" = true ] ; then
  set +e
  echo "Adding white noise to initial condition"
  cdv "$atm_to_cam_path"
  (set -x; python perturb_white_noise.py "${sePreFilterIC}") ; exit_status=$?
  check_python_exit "perturb_white_noise.py" $exit_status
  set -e
fi

if [ "${add_perturbs}" = true ] ; then
  echo "Adding perturbations"
  cdv "$atm_to_cam_path/perturb"
  set +e

  sstFileIC_WPERT=${sstFileIC}_PERT.nc

  (set -x; python add_perturbations_to_sst.py \
     --BEFOREPERTFILE "${sstFileIC}" \
     --AFTERPERTFILE "${sstFileIC_WPERT}" \
     --pthi "${perturb_namelist}") ; exit_status=$?
  check_python_exit "add_perturbations_to_sst.py" $exit_status

  echo "SST perturbations added successfully"

  sePreFilterIC_WPERT=${sePreFilterIC}_PERT.nc

  (set -x; python add_perturbations_to_cam.py \
     --BEFOREPERTFILE "${sePreFilterIC}" \
     --AFTERPERTFILE "${sePreFilterIC_WPERT}" \
     --gridfile "${modelgridfile}" \
     --MAPFILEPATH "${mapping_files_path}" \
     --pthi "${perturb_namelist}") ; exit_status=$?
  check_python_exit "add_perturbations_to_cam.py" $exit_status

  echo "ATM perturbations added successfully"

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