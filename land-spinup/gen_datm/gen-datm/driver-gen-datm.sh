#!/bin/bash -l

################################################################
#### Casper
################################################################
#PBS -N forcing-gen
#PBS -A UPSU0032
#PBS -l select=1:ncpus=4:mpiprocs=4:mem=80GB
#PBS -l walltime=24:00:00
#PBS -q casper@casper-pbs
#PBS -j oe
################################################################

module load parallel
#module load ncl
module load conda
conda activate npl

echo "parallel location: $(which parallel)"
#echo "ncl location: $(which ncl)"

NUMCORES=4

STYR=2023
ENYR=2024
PATHTORAWERA5=/glade/campaign/cgd/amp/zarzycki/DATM_FORCING/raw-ERA5/
OUTDIRBASE=/glade/campaign/cgd/amp/zarzycki/DATM_FORCING/ERA5/

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
      #LINECOMMAND="ncl gen-forcing.ncl YYYY=$DATA_YEAR 'MM=\"$DATA_MONTH\"' 'RAWERA5FILE=\"$THISFILE\"' 'outdirbase=\"$OUTDIRBASE\"' "
      LINECOMMAND="python gen-forcing.py --year=$DATA_YEAR --month=$DATA_MONTH --era5_file=\"$THISFILE\" --outdirbase=\"$OUTDIRBASE\" --do_q --do_flds --greg_to_noleap --convert_nc3"
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
