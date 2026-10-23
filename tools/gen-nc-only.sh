#!/bin/bash

# Environment variables
#export BETACAST=/glade/u/home/zarzycki/betacast/
export BETACAST=/global/homes/c/czarzyck/betacast/
mapping_files_path=$PWD
cleanup_new_weights=false  # delete any map_*.nc this script generates in mapping_files_path

# Model configuration
DYCORE="se"
modelSystem=1   # 0 = CESM, 1 = E3SM (adds Linoz chem), 2 = SCREAM (CDF5 conversion)
adjust_flags="-"
YYYYMMDDHH=1979010100

# For CESM, np is the gridfile and topofile
# For E3SM, np is the gridfile, pg is the topofile
modelgridfile="/global/homes/c/czarzyck/m2637/E3SM_SCREAM_files/grids/scrip/ne0np4alcc.ne30x4_np4_scrip.nc"
adjust_topo="/global/homes/c/czarzyck/m2637/E3SM_SCREAM_files/topo/USGS-gtopo30_ne0np4alcc.ne30x4np4pg2_x6t.nc"
sePreFilterIC="/global/homes/c/czarzyck/m2637/E3SM_SCREAM_files/ncdata/IC_ERA5_ne0np4alcc.ne30x4np4pg2_eam_${YYYYMMDDHH}.nc"
numLevels=80
#anl2mdlWeights=""   # optional, generated from modelgridfile if empty

# Source data configuration
# atmDataType: 1 = GFS, 2 = ERAI, 3 = CFSR, 4 = ERA5, 5 = CR20V3, 9 = CAM (m2m)
# Local files are read from gfs_files_path/era_files_path, ERA5 with ERA5RDA=1 from RDADIR
atmDataType=4
ERA5RDA=1
#RDADIR="/glade/u/home/zarzycki/rda/ds633.0/"
RDADIR="/global/cfs/projectdirs/m3522/cmip6/ERA5/"
#gfs_files_path=""
#era_files_path=""

# Frankengrid (regional_name: hrrr_3km, hrrr_3km_ml, rap_13km, hwrf)
do_frankengrid=false
regional_name="hrrr_3km"
regional_src="/global/homes/c/czarzyck/m2637/betacast/sewx/HRRR/HRRR_wrfprs_f00_YYYYMMDDHH.grib2"
regional_fail_if_missing=true
keep_regional_file=true

# Flags appended to every atm_to_cam.py call
ATM_TO_CAM_EXTRA_ARGS="--compress_file --write_floats --verbose"
#ATM_TO_CAM_EXTRA_ARGS+=" --mpas_as_cam"
#ATM_TO_CAM_EXTRA_ARGS+=" --write_debug_files --write_debug_dir /pscratch/sd/c/czarzyck/"

# Anything else runner_atm_to_cam.sh understands (cr20v3_*, m2m_*, standalone_vortex,
# vortex_namelist, add_noise, add_perturbs, AUGMENT_STR, VORTEX_STR, ...) can be set here too

# EDIT ABOVE HERE!
#----------------------------------------------------------------------------------------

set -e

yearstr=${YYYYMMDDHH:0:4}
monthstr=${YYYYMMDDHH:4:2}
daystr=${YYYYMMDDHH:6:2}
cyclestr=${YYYYMMDDHH:8:2}

# Export everything set above (and any extra runner options) to the runner
export BETACAST mapping_files_path DYCORE modelSystem adjust_flags adjust_topo \
       modelgridfile sePreFilterIC numLevels anl2mdlWeights \
       atmDataType ERA5RDA RDADIR gfs_files_path era_files_path \
       do_frankengrid regional_name regional_src regional_fail_if_missing keep_regional_file \
       ATM_TO_CAM_EXTRA_ARGS yearstr monthstr daystr cyclestr \
       cr20v3_mean_dir cr20v3_member cr20v3_member_dir cr20v3_orog_file \
       cr20v3_blend cr20v3_blend_lev cr20v3_blend_taper \
       m2m_gridfile m2m_remap_file m2m_parent_source m2m_topo_in \
       standalone_vortex vortex_namelist add_noise add_perturbs perturb_namelist \
       sstFileIC sstDataType AUGMENT_STR VORTEX_STR

shopt -s nullglob
weights_before=("${mapping_files_path}"/map_*.nc)

bash "${BETACAST}/runner_atm_to_cam.sh"

if [ "$cleanup_new_weights" = true ]; then
  for f in "${mapping_files_path}"/map_*.nc; do
    if [[ ! " ${weights_before[*]} " == *" ${f} "* ]]; then
      rm -fv "$f"
    fi
  done
fi

echo "Done with generation of: ${sePreFilterIC}"
