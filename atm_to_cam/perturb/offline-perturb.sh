#!/bin/bash

SERVER_NAME=$(hostname -A) ; echo $SERVER_NAME
if [[ $SERVER_NAME == *"perlmutter"* ]] || [[ $SERVER_NAME == *"nid0"* ]]; then
  echo "Using pm-cpu"
  export NCARG_ROOT=/global/homes/c/czarzyck/.conda/pkgs/ncl-6.6.2-h3fdc804_41/
  PATHTONCL=/global/homes/c/czarzyck/.conda/envs/e3sm_unified_1.8.1_nompi/bin/
  module load parallel
  NUMCORES=32
elif [[ $SERVER_NAME == *"casper"* ]]; then
  echo "Using Casper"
else
  echo "Unrecognized server. exiting"
  exit
fi

THISDIR=$PWD
source ${THISDIR}/../../utils.sh   # Source external bash functions

PERTNAME="Plus2K"
path_to_inputdata=/global/cfs/cdirs/m2637/betacast/sewx/
perturb_namelist=../../namelists/PNW_2021/perturb.101.plus2K.nl
modelgridfile=/global/homes/c/czarzyck/m2637/betacast/cesmfiles/grids/conus_30_x8.g_scrip.nc
NUDGINGBASEDIR=~/scratch/ndg/se_conus_30_x8_L72/ERA5/

# Automated stuff done once.
mapping_files_path=${path_to_inputdata}/mapping/ ; mkdir -p ${mapping_files_path}

# sePreFilterIC passed in
#modelgridfile --> SCRIP

for sePreFilterIC in ${NUDGINGBASEDIR}/*nc; do
  echo "Doing: $sePreFilterIC"
  # Create a new subfolder with PERTNAME where the base files are located
  sePreFilterIC_WPERT="${sePreFilterIC%/*}/${PERTNAME}/${sePreFilterIC##*/}"
  # Make the subdirectory
  mkdir -pv "${sePreFilterIC_WPERT%/*}"
  echo "Writing to: $sePreFilterIC_WPERT"
  (set -x; ncl -n add_perturbations_to_cam.ncl 'BEFOREPERTFILE="'${sePreFilterIC}'"'  \
     'AFTERPERTFILE = "'${sePreFilterIC_WPERT}'"' \
     'gridfile = "'${modelgridfile}'"' \
     'MAPFILEPATH = "'${mapping_files_path}'"' \
     'pthi="'${perturb_namelist}'"' ) ; exit_status=$?
  check_ncl_exit "add_perturbations_to_cam.ncl" $exit_status
done