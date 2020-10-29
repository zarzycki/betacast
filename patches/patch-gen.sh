#!/bin/bash

CESMROOT=/global/homes/c/czarzyck/E3SM-dev/
filename="rof_comp_mct"
SUBPATH="components/mosart/src/cpl/"
CESMCOMPONENT="mosart"
SOURCEMODDIR=~/betacast/patches/
mkdir -p ${SOURCEMODDIR}
diff -b -U 3 ${CESMROOT}/${SUBPATH}/${filename}.F90 /global/homes/c/czarzyck/v1.2-CONUS-F2010C5-v3/SourceMods/src.${CESMCOMPONENT}/${filename}.F90 > ${SOURCEMODDIR}/${filename}.patch
