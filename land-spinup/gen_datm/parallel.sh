#!/bin/bash -l

#SBATCH --job-name=MERRA_mpi
#SBATCH --account=P05010048
#SBATCH --ntasks=4
#SBATCH --ntasks-per-node=4
#SBATCH --time=18:00:00
#SBATCH --mem=60G
#SBATCH --partition=dav
#SBATCH --output=MERRA_mpi.out.%j

module load parallel
module load ncl
module load nco

# use numcores = 32 for 1deg runs, 16 for high-res runs
NUMCORES=4
TIMESTAMP=`date +%s%N`
COMMANDFILE=commands.${TIMESTAMP}.txt

declare -a months=("01" "02" "03" "04" "05" "06" "07" "08" "09" "10" "11" "12")

for DATA_YEAR in {2020..2020}
do
  for DATA_MONTH in "${months[@]}"
  do
    LINECOMMAND="ncl gen-forcing.ncl YYYY=$DATA_YEAR 'MM=\"$DATA_MONTH\"'"
    echo ${LINECOMMAND} >> ${COMMANDFILE}
  done
done

parallel --jobs ${NUMCORES} --workdir $PWD < ${COMMANDFILE}

rm ${COMMANDFILE}

