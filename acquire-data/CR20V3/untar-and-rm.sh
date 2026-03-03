#!/bin/bash
#PBS -N casper-untar
#PBS -A P93300042
#PBS -l select=1:ncpus=8:mem=200GB
#PBS -l walltime=23:59:00
#PBS -q casper@casper-pbs
#PBS -j oe

module load parallel
NCPUS=8
TARDIR="/glade/derecho/scratch/zarzycki/globus/"

set -euo pipefail
shopt -s nullglob

cd ${TARDIR}

# Define the function
untar_and_rm() {
  local tarfile="$1"
  echo "======================================="
  echo "Processing: $tarfile"
  echo "Start time: $(date)"

  if tar -xf "$tarfile"; then
    echo "Extraction successful for $tarfile"
    echo "Deleting $tarfile"
    rm -f "$tarfile"
  else
    echo "ERROR extracting $tarfile"
    echo "Skipping delete."
  fi

  echo "End time: $(date)"
  echo
}
export -f untar_and_rm

# Now call the function
parallel -j ${NCPUS} untar_and_rm ::: *.tar

echo "All done!"
