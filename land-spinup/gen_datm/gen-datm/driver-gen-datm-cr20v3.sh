#!/bin/bash -l

################################################################
#### Casper
################################################################
#PBS -N cr20v3-datm-gen
#PBS -A P93300042
#PBS -l select=1:ncpus=4:mpiprocs=4:mem=80GB
#PBS -l walltime=24:00:00
#PBS -q casper@casper-pbs
#PBS -j oe
################################################################

module load parallel
module load conda
conda activate npl

echo "parallel location: $(which parallel)"

NUMCORES=4

STYR=1850
ENYR=1875
RAWDATADIR=/glade/u/home/zarzycki/rda/d131003/anl
#OUTDIRBASE=/glade/derecho/scratch/zarzycki/CR20V3-DATM/
OUTDIRBASE=/glade/campaign/cgd/amp/zarzycki/DATM_FORCING/CR20V3/

TIMESTAMP=$(date +%s%N)
COMMANDFILE=commands.${TIMESTAMP}.txt

ERR="0"

for DATA_YEAR in $(seq $STYR $ENYR)
do
    THISFILE=${RAWDATADIR}/anl_mean_${DATA_YEAR}_TMP_2m.nc
    if [ -f "$THISFILE" ]; then
        LINECOMMAND="python gen-forcing-cr20v3.py --year=${DATA_YEAR} --rawdatadir=\"${RAWDATADIR}\" --outdirbase=\"${OUTDIRBASE}\" --greg_to_noleap --time_stride 6 --convert_nc3"
        echo ${LINECOMMAND} >> ${COMMANDFILE}
    else
        echo "$THISFILE does not exist!"
        ERR="yes"
    fi
done

if [[ ${ERR} != "0" ]]; then
    echo "Found missing files, not running gen script..."
else
    echo "Good to go! Submitting ${NUMCORES} parallel jobs."
    parallel --jobs ${NUMCORES} --workdir $PWD < ${COMMANDFILE}
fi

rm -v ${COMMANDFILE}
