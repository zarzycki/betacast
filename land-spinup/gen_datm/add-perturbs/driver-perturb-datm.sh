#!/bin/bash -l

#PBS -N forcing-gen 
#PBS -A UNSB0017
#PBS -l select=1:ncpus=8:mpiprocs=8:mem=60GB
#PBS -l walltime=24:00:00
#PBS -q casper

ORIGPATH=/glade/u/home/zarzycki/scratch/ERA5-DATM/DATM-1996/
PERTURBPATH=/glade/u/home/zarzycki/scratch/ERA5-DATM/DATM-1996-Plus5K/
DIRTOSCRIPTS=./

NUMCORES=4
# true: after perturbing (in place, full precision), re-quantize + deflate the modified TPQWL and Prec files
# (netCDF4, 1 time step per chunk). Solar is not perturbed (no --do_sw_fluxes) so it is left as is.
LOSSY_COMPRESS=true
CMPSTR="ncks -O -4 -L 1 --cnk_plc=g2d --cnk_map=rd1"
TIMESTAMP=`date +%s%N`
COMMANDFILE=commands.${TIMESTAMP}.txt

module load parallel
module load conda
conda activate npl

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
  
  LINECOMMAND="python ${DIRTOSCRIPTS}/add_perturbations_to_datm.py --datm_file_name=\"${TPQWLf}\" --datm2_file_name=\"${Precf}\" --datm3_file_name=\"${Solarf}\""
  if [ "$LOSSY_COMPRESS" = true ]; then
    LINECOMMAND="${LINECOMMAND} && ${CMPSTR} --qnt 'TBOT=.2#PSRF=.0#QBOT=4#WIND=.2#FLDS=.1' \"${TPQWLf}\" \"${TPQWLf}\" && ${CMPSTR} --qnt PRECTmms=3 \"${Precf}\" \"${Precf}\""
  fi
  echo ${LINECOMMAND} >> ${COMMANDFILE}
done

parallel --jobs ${NUMCORES} --workdir $PWD < ${COMMANDFILE}

rm ${COMMANDFILE}
