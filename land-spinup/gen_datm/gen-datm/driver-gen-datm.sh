#!/bin/bash -l

################################################################
#### Casper
################################################################
#PBS -N forcing-gen
#PBS -A P93300642
#PBS -l select=1:ncpus=4:mpiprocs=4:mem=80GB
#PBS -l walltime=24:00:00
#PBS -q casper@casper-pbs
#PBS -j oe
################################################################

#module load parallel
#module load ncl

NUMCORES=4

STYR=1992
ENYR=1992
PATHTORAWERA5=/glade/derecho/scratch/${LOGNAME}/ERA5-DATM/

TIMESTAMP=`date +%s%N`
COMMANDFILE=commands.${TIMESTAMP}.txt

declare -a months=("01" "02" "03" "04" "05" "06" "07" "08" "09" "10" "11" "12")

ERR="0"

for DATA_YEAR in `seq $STYR $ENYR`
do
  for DATA_MONTH in "${months[@]}"
  do
    THISFILE=${PATHTORAWERA5}/out.${DATA_YEAR}.${DATA_MONTH}.nc
    if [ -f "$THISFILE" ]; then
      LINECOMMAND="ncl gen-forcing.ncl YYYY=$DATA_YEAR 'MM=\"$DATA_MONTH\"' 'RAWERA5FILE=\"$THISFILE\"'"
      echo ${LINECOMMAND} >> ${COMMANDFILE}
    else
      echo "$THISFILE does not exist!"
      ERR="yes"
    fi
  done
done

if [[ ${ERR} != "0" ]] ; then
  echo "Found errors, not running gen script..."
else
  echo "Good to go!"
  parallel --jobs ${NUMCORES} --workdir $PWD < ${COMMANDFILE}
fi

rm -v ${COMMANDFILE}