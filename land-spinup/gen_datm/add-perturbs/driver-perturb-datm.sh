#!/bin/bash -l

#PBS -N forcing-gen 
#PBS -A UNSB0017
#PBS -l select=1:ncpus=8:mpiprocs=8:mem=60GB
#PBS -l walltime=24:00:00
#PBS -q casper

ORIGPATH=/glade/u/home/zarzycki/scratch/ERA5-DATM/DATM-perturb3/
PERTURBPATH=/glade/u/home/zarzycki/scratch/ERA5-DATM/DATM-perturb4/
DIRTOSCRIPTS=./

NUMCORES=4
TIMESTAMP=`date +%s%N`
COMMANDFILE=commands.${TIMESTAMP}.txt

module load parallel
module load ncl

echo "Cleaning up old perturb files..."
rm -rfv ${PERTURBPATH}
echo "Copying reference DATM files..."
cp -v -rf ${ORIGPATH} ${PERTURBPATH}

for f in `find ${PERTURBPATH} -name "*TPQWL*nc"`
do
  echo $f
  TPQWLf=$f
  tmp1=${f/TPHWL6Hrly/Precip6Hrly}
  Precf=${tmp1/TPQWL/Prec}
  tmp1=${f/TPHWL6Hrly/Solr}
  Solarf=${tmp1/TPQWL/Solar}
  
  #ls ${TPQWLf}
  #ls ${Precf}
  #ls ${Solarf}
  
  LINECOMMAND="ncl ${DIRTOSCRIPTS}/add_perturbations_to_datm.ncl 'datm_file_name=\"'${TPQWLf}'\"' 'datm2_file_name=\"'${Precf}'\"' 'datm3_file_name=\"'${Solarf}'\"'    "
  echo ${LINECOMMAND} >> ${COMMANDFILE}
done

exit

parallel --jobs ${NUMCORES} --workdir $PWD < ${COMMANDFILE}

rm ${COMMANDFILE}
