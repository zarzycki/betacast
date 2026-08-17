#!/bin/bash

STYR=1980
ENYR=2001

echo "Processing years from ${STYR} to ${ENYR}"

for year in $(seq ${STYR} ${ENYR}); do
  echo "----------------------------------------"
  echo "Running for year: ${year}"
  echo "----------------------------------------"
  python convert-ebt-to-tcv.py ${year}
done

echo "All years processed!"
