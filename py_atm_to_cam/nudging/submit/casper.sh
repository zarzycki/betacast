#!/bin/bash
#PBS -N gen_nudge_betacast
#PBS -A P93300042
#PBS -l select=1:ncpus=12:mem=200GB
#PBS -l walltime=01:59:00
#PBS -q casper@casper-pbs
#PBS -j oe

module load conda && conda activate npl

export BETACAST=/glade/u/home/zarzycki/betacast
export OUTDIR=/glade/derecho/scratch/zarzycki/ndg
export RDADIR=/glade/campaign/collections/rda/data/d633000
export NUMCORES=12

SCRIPTDIR="${BETACAST}/py_atm_to_cam/nudging"

# When run interactively, use positional args; via qsub, use -v env vars
# e.g., qsub -v NLFILE="ndg.era5.nl",input_dates_file="dates.txt",INDEX="001" submit/casper.sh
if [ $# -gt 0 ]; then
  exec "${SCRIPTDIR}/gen-nudge.sh" "$@"
else
  ARGS=()
  [ -n "${NLFILE:-}" ] && ARGS+=("$NLFILE")
  [ -n "${input_dates_file:-}" ] && ARGS+=("$input_dates_file")
  [ -n "${INDEX:-}" ] && ARGS+=("$INDEX")
  exec "${SCRIPTDIR}/gen-nudge.sh" "${ARGS[@]}"
fi
