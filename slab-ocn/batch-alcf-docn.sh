#!/bin/bash

SLABDIR=/glade/u/home/zarzycki/betacast/slab-ocn
OUTDIR=/glade/u/home/zarzycki/scratch/CONUSClimateSims

for f in ${OUTDIR}/sst_cice*nc ; do

  ORIGFILE=$(basename $f)

  cd $SLABDIR
  ln -s ${OUTDIR}/${ORIGFILE} .

  cd $SLABDIR
  python external-to-docn.py --infile $ORIGFILE

  cd /glade/u/home/zarzycki/work/cesm2_2_0/components/cam/tools/icesst/bcgen
  ./bcgen -i $SLABDIR/${ORIGFILE}_CIMEified.nc -c sstice_clim.nc -t $SLABDIR/${ORIGFILE}_diddled.nc < namelist
  rm -v sstice_clim.nc

  cd $SLABDIR
  python h-qdp.py --sst_file ${ORIGFILE}_diddled.nc --output_file ${OUTDIR}/final_${ORIGFILE}
  rm -v $SLABDIR/${ORIGFILE}_CIMEified.nc
  rm -v $SLABDIR/${ORIGFILE}_diddled.nc
  rm -v $SLABDIR/${ORIGFILE}

done