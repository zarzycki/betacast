#!/bin/bash

NCLDIR=/glade/u/home/zarzycki/betacast/slab-ocn
OUTDIR=/glade/u/home/zarzycki/scratch/CONUSClimateSims

for f in ${OUTDIR}/sst_cice*nc ; do

  ORIGFILE=$(basename $f)

  cd $NCLDIR
  ln -s ${OUTDIR}/${ORIGFILE} .

  cd $NCLDIR
  ncl external-to-docn.ncl 'infile="'$ORIGFILE'"'

  cd /glade/u/home/zarzycki/work/cesm2_2_0/components/cam/tools/icesst/bcgen
  ./bcgen -i $NCLDIR/${ORIGFILE}_CIMEified.nc -c sstice_clim.nc -t $NCLDIR/${ORIGFILE}_diddled.nc < namelist
  rm -v sstice_clim.nc

  cd $NCLDIR
  ncl add-h-qdp-to-sst.ncl 'sfile="'$ORIGFILE'_diddled.nc"'

  mv simple.nc ${OUTDIR}/final_${ORIGFILE}
  rm -v $NCLDIR/${ORIGFILE}_CIMEified.nc
  rm -v $NCLDIR/${ORIGFILE}_diddled.nc
  rm -v $NCLDIR/${ORIGFILE}

done