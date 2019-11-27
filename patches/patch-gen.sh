#!/bin/bash

  CESMROOT=/global/homes/c/czarzyck/E3SM
  filename="lnd_comp_mct"
  SUBPATH="components/clm/src/cpl/"
  CESMCOMPONENT="clm"
  SOURCEMODDIR=~/betacast/patches/
  mkdir -p ${SOURCEMODDIR}
  diff -b -U 3 ${CESMROOT}/${SUBPATH}/${filename}.F90 /global/homes/c/czarzyck/F-betacast-FC5AV1C/SourceMods/src.${CESMCOMPONENT}/${filename}.F90 > ${SOURCEMODDIR}/${filename}.patch
