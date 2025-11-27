#!/bin/bash

set -euo pipefail

# BETACAST must be set so we can find utils, etc.
if [ -z "${BETACAST}" ]; then
  echo "BETACAST must be set in environment" >&2
  exit 1
fi

source "${BETACAST}/utils.sh"
source "${BETACAST}/datahelpers.sh"

echo "SST runner host: $(hostname)"
echo "SST runner start time: $(date -u +"%Y-%m-%d %H:%M:%S UTC")"

############################### SST TO CAM ###############################

# Switch bash bool to int for NCL / Python input
if [ "${predict_docn}" = true ]; then
  INT_PREDICT_DOCN=1
else
  INT_PREDICT_DOCN=0
fi

# Check if domain or SCRIP exist, if one is missing create both domain and scrip
if [ ! -f "${sst_domain_file}" ] || [ ! -f "${sst_scrip_file}" ]; then
  echo "Creating SST domain file for: ${docnres}"
  if [ "${DO_PYTHON}" = true ]; then
    # TODO: Python equivalent if desired
    echo "Python SST domain generation not implemented (DO_PYTHON=true)"; exit 1
  else
    set +e
    (set -x; ncl gen-sst-domain.ncl 'inputres="'${docnres}'"' ) ; exit_status=$?
    check_ncl_exit "gen-sst-domain.ncl" $exit_status
    set -e
  fi
  compress_single_file "${sst_scrip_file}"
fi

# If the CIME coupler is nuopc, we need to generate an ESMF file
if [ "${cime_coupler}" = "nuopc" ] && [ ! -f "${sst_ESMF_file}" ]; then
  ESMF_Scrip2Unstruct "${sst_scrip_file}" "${sst_ESMF_file}" 0
  rm -fv PET0.ESMF_LogFile || true
  compress_single_file "${sst_ESMF_file}"
fi

# Now generate the SST/ice datastream
if [ "${DO_PYTHON}" = true ]; then
  (set -x; python sst_to_cam.py \
      --initdate "${yearstr}${monthstr}${daystr}${cyclestr}" \
      --predict_docn "${INT_PREDICT_DOCN}" \
      --inputres "${docnres}" \
      --datasource "${SSTTYPE}" \
      --sstDataFile "${sst_files_path}/${sstFile}" \
      --iceDataFile "${sst_files_path}/${iceFile}" \
      --SST_write_file "${sstFileIC}" \
      --smooth_ice \
      --smooth_iter 3 \
      --verbose
  )
else
  set +e
  (set -x; ncl sst_interp.ncl \
      'initdate="'${yearstr}${monthstr}${daystr}${cyclestr}'"' \
      predict_docn=${INT_PREDICT_DOCN} \
      'inputres="'${docnres}'"' \
      'datasource="'${SSTTYPE}'"' \
      'sstDataFile = "'${sst_files_path}/${sstFile}'"' \
      'iceDataFile = "'${sst_files_path}/${iceFile}'"' \
      'SST_write_file = "'${sstFileIC}'"' \
  ) ; exit_status=$?
  check_ncl_exit "sst_interp.ncl" $exit_status
  set -e
fi

if [ "$do_slab" = true ] ; then
  ## Add slab ocean parameters
  (set -x; python ${BETACAST}/slab-ocn/h-qdp.py \
      --sst_file "${sstFileIC}" \
      --output_file "${sstFileIC}"
  )
fi

echo "SST runner complete."
