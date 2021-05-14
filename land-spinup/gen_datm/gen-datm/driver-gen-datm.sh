#!/bin/bash -l

#PBS -N forcing-gen 
#PBS -A UNSB0017
#PBS -l select=1:ncpus=4:mpiprocs=4:mem=80GB
#PBS -l walltime=24:00:00
#PBS -q casper

module load parallel
module load ncl

STYR=1995
ENYR=1996
NUMCORES=4
TIMESTAMP=`date +%s%N`
COMMANDFILE=commands.${TIMESTAMP}.txt

declare -a months=("01" "02" "03" "04" "05" "06" "07" "08" "09" "10" "11" "12")

for DATA_YEAR in `seq $STYR $ENYR`
do
  for DATA_MONTH in "${months[@]}"
  do
    LINECOMMAND="ncl gen-forcing.ncl YYYY=$DATA_YEAR 'MM=\"$DATA_MONTH\"'"
    echo ${LINECOMMAND} >> ${COMMANDFILE}
  done
done

parallel --jobs ${NUMCORES} --workdir $PWD < ${COMMANDFILE}

rm ${COMMANDFILE}

