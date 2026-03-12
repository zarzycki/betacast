#!/bin/bash
#SBATCH --qos=regular
#SBATCH --time=12:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=128
#SBATCH --constraint=cpu

module load conda && conda activate betacast

export BETACAST=/global/homes/c/czarzyck/betacast
export OUTDIR=/pscratch/sd/c/czarzyck/ndg
export RDADIR=/global/cfs/projectdirs/m3522/cmip6/ERA5
export NUMCORES=24

SCRIPTDIR="${BETACAST}/py_atm_to_cam/nudging"

# sbatch passes trailing args directly; also supports interactive use
# e.g., sbatch submit/pm-cpu.sh ndg.era5.nl dates.txt 001
exec "${SCRIPTDIR}/gen-nudge.sh" "$@"
