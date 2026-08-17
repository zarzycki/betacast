#!/bin/bash

# Test: gen-sst-domain
# Generates domain + SCRIP files at 180x360 resolution.
# Fully self-contained — no external data required.

RES="180x360"

# Function to run Python test
run_python_test() {
    python ${BETACAST}/sst_to_cam/gen-sst-domain.py --inputres ${RES}
}
