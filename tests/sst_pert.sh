#!/bin/bash

# Test: add_perturbations_to_sst Python vs NCL
# Self-contained test that generates synthetic delta SST data,
# creates an SST file via NOAAOI, then applies PGW SST perturbations.

# SST generation settings
DATE="2023052300"
RES="180x360"
DATASRC="NOAAOI"

# Test-specific paths
SYNTH_BASEDIR="${DEBUG_FILE_DIR}/sst_pert_testdata"
PERTURB_NAMELIST="${SYNTH_BASEDIR}/test_perturb.nl"

# Generate synthetic delta SST data and namelist for testing
generate_test_data() {
    echo "Generating synthetic delta SST data for perturbation test..."
    mkdir -p "${SYNTH_BASEDIR}/CESMLENS_plev"

    python -c "
import numpy as np
from netCDF4 import Dataset

# Create a synthetic ens_SST_anom.nc with a smooth warming pattern
# Use current_year=1921, comp_year=1922 so we need minimal time steps
# Index for 1921: (1921-1920)*12 + 0 = 12
# Index for 1922: (1922-1920)*12 + 0 = 24
# So we need at least 25 time steps

ntime = 25
nlat = 73   # 2.5 degree
nlon = 144  # 2.5 degree

lat = np.linspace(-90, 90, nlat)
lon = np.linspace(0, 357.5, nlon)

# Create SST anomaly: cos(lat) pattern scaled to ~2K at equator
# Vary slightly over time so current and comp years differ
lat_rad = np.deg2rad(lat)
base_pattern = 2.0 * np.cos(lat_rad)[:, None] * np.ones((1, nlon))

ts = np.zeros((ntime, nlat, nlon))
for t in range(ntime):
    # Scale grows with time index to create meaningful delta
    ts[t, :, :] = base_pattern * (0.5 + 0.1 * t)

ds = Dataset('${SYNTH_BASEDIR}/CESMLENS_plev/ens_SST_anom.nc', 'w')
ds.createDimension('time', ntime)
ds.createDimension('lat', nlat)
ds.createDimension('lon', nlon)

v = ds.createVariable('TS', 'f8', ('time', 'lat', 'lon'))
v[:] = ts

v = ds.createVariable('lat', 'f8', ('lat',))
v.units = 'degrees_north'
v[:] = lat

v = ds.createVariable('lon', 'f8', ('lon',))
v.units = 'degrees_east'
v[:] = lon

ds.close()
print('Synthetic delta file written.')
"

    # Write test namelist
    cat > "${PERTURB_NAMELIST}" << 'NLEOF'
! Test perturbation namelist (synthetic data)
case         = "CESMLENS"
NLEOF
    echo "basedir      = \"${SYNTH_BASEDIR}\"" >> "${PERTURB_NAMELIST}"
    cat >> "${PERTURB_NAMELIST}" << 'NLEOF'
start_month  = 1
end_month    = 1
current_year = 1921
comp_year    = 1922
adjust_ice        = True
output_sst_diag   = False
NLEOF

    echo "Test namelist written to ${PERTURB_NAMELIST}"
    cat "${PERTURB_NAMELIST}"
}

# Generate SST IC file (shared input for both Python and NCL perturbation runs)
generate_sst_input() {
    echo "Generating SST input file for perturbation test..."
    python ${BETACAST}/py_sst_to_cam/sst_to_cam.py \
        --initdate ${DATE} \
        --predict_docn 0 \
        --inputres ${RES} \
        --datasource ${DATASRC} \
        --sstDataFile "${TEST_FILES_DIR}/sst.day.mean.2023.nc" \
        --iceDataFile "${TEST_FILES_DIR}/icec.day.mean.2023.nc" \
        --SST_write_file "${DEBUG_FILE_DIR}/sst_pert_input.nc" \
        --smooth_ice \
        --smooth_iter 3 \
        --verbose
}

# Function to run Python test
run_python_test() {
    generate_test_data
    generate_sst_input

    python ${BETACAST}/py_atm_to_cam/perturb/add_perturbations_to_sst.py \
        --BEFOREPERTFILE "${DEBUG_FILE_DIR}/sst_pert_input.nc" \
        --AFTERPERTFILE "${DEBUG_FILE_DIR}/sst_pert_python.nc" \
        --pthi "${PERTURB_NAMELIST}"
}

# Function to run NCL test
run_ncl_test() {
    # NCL script uses relative load paths — must run from atm_to_cam/perturb/
    pushd ${BETACAST}/atm_to_cam/perturb > /dev/null
    ncl -n add_perturbations_to_sst.ncl \
        'BEFOREPERTFILE="'${DEBUG_FILE_DIR}'/sst_pert_input.nc"' \
        'AFTERPERTFILE="'${DEBUG_FILE_DIR}'/sst_pert_ncl.nc"' \
        'pthi="'${PERTURB_NAMELIST}'"'
    local ncl_status=$?
    popd > /dev/null
    return ${ncl_status}
}

# Function to validate results
run_validation() {
    python ${BETACAST}/py_testing/check_same.py \
        ${DEBUG_FILE_DIR}/sst_pert_python.nc \
        ${DEBUG_FILE_DIR}/sst_pert_ncl.nc
}

# Optional cleanup function
cleanup() {
    rm -f ${DEBUG_FILE_DIR}/sst_pert_input.nc
    rm -f ${DEBUG_FILE_DIR}/sst_pert_python.nc
    rm -f ${DEBUG_FILE_DIR}/sst_pert_ncl.nc
    rm -rf ${SYNTH_BASEDIR}
    rm -f deltas_sst.nc
}
