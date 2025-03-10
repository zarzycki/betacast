#!/bin/bash
################################################################
#### PMCPU
################################################################
#SBATCH --qos=premium
#SBATCH --time=11:30:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=128
#SBATCH --constraint=cpu
################################################################

dryrun=false
MODEL_DATASTREAM_DIR=/pscratch/sd/c/czarzyck/hyperion/CORI.VR28.NATL.WAT.CAM5.4CLM5.0.dtime900.003/h5/
MAPFILE=~/m2637/betacast/regrid_maps/map_ne0np4natlanticwat.ne30x4_TO_era5_0.25x0.25_patc.nc
OUTDIRBASE="/global/homes/c/czarzyck/m2637/DATM_FORCING/"

SCRIPTDIR=$PWD

SERVER_NAME=$(hostname -A)
echo $SERVER_NAME

if [[ $SERVER_NAME == *"perlmutter"* ]] || [[ $SERVER_NAME == *"nid0"* ]]; then
  echo "Using pm-cpu"
  export NCARG_ROOT=/global/homes/c/czarzyck/.conda/pkgs/ncl-6.6.2-h3fdc804_41/
  PATHTONCL=/global/homes/c/czarzyck/.conda/envs/e3sm_unified_1.8.1_nompi/bin/
  source /global/common/software/e3sm/anaconda_envs/load_latest_e3sm_unified_pm-cpu.sh
  #source /global/common/software/e3sm/anaconda_envs/load_e3sm_unified_1.9.3_pm-cpu.sh
  module load parallel
  NUMCORES=32
elif [[ $SERVER_NAME == *"casper"* ]] || [[ $SERVER_NAME == *"crhtc"* ]]; then
  echo "Using Casper"
  module load parallel
  module load peak_memusage
  NUMCORES=36
else
  echo "Unrecognized server. exiting"
  exit
fi

if ! type ncks &> /dev/null ; then
  echo "ERROR: ncks does not exist. Make sure ncks is in your path when betacast is invoked"
  exit 1
fi

# GNUPARALLEL SETTINGS
TIMESTAMP=`date +%s%N`
COMMANDFILE=commands.${TIMESTAMP}.txt

FILES=`find ${MODEL_DATASTREAM_DIR} -type f -name "*.nc"`

for f in $FILES; do
  COMMAND="bash $SCRIPTDIR/single-file-to-datm.sh $f $OUTDIRBASE $MAPFILE"
  echo ${COMMAND} >> ${COMMANDFILE}
done

if [ $dryrun = false ] ; then
  echo "Launching GNU parallel"
  # Launch GNU parallel
  if [[ $SERVER_NAME == *"perlmutter"* ]] || [[ $SERVER_NAME == *"nid0"* ]]; then
    parallel --jobs ${NUMCORES} < ${COMMANDFILE}
  elif [[ $SERVER_NAME == *"casper"* ]] || [[ $SERVER_NAME == *"crhtc"* ]]; then
    peak_memusage.exe parallel --jobs ${NUMCORES} --workdir $PWD < ${COMMANDFILE}
  else
    echo "unknown" ; exit
  fi

  # Cleanup
  rm ${COMMANDFILE}
fi

