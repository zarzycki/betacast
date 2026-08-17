#!/bin/bash

CONNECTFILE=/global/homes/c/czarzyck/scratch/philly128connect.dat
UQSTR=Philly128x8
PATHTOFILES=/global/homes/c/czarzyck/scratch/e3sm_scratch/pm-cpu/Philly-F2010SCREAMHRDYAMOND2SP-Philly128x8-101-control/run/2075090400_uv3/
HSTREAMTRACK="h2"
TIMESTRIDE=1
TEMPESTEXTREMESDIR=~/sw/tempestextremes_noMPI/

# The TE trajectory file for this particular betacast case
TRAJFILE=trajectories.txt.${UQSTR}

### Flag needed if using unstructured data (if not using unstruc data, empty string)
CONNECTFLAG="--in_connect ${CONNECTFILE}"

### Path + filelist of data to process
FILES=`ls ${PATHTOFILES}/*.?am.${HSTREAMTRACK}.*.nc`

############ TRACKER MECHANICS #####################

### THIS ONE!
DCU_PSLFOMAG=100.0
DCU_PSLFODIST=5.5
DCU_WCFOMAG=-6.0
DCU_WCFODIST=6.5
DCU_WCMAXOFFSET=1.0
DCU_MERGEDIST=3.0
SN_TRAJRANGE=3.0
SN_TRAJMINLENGTH=5
SN_TRAJMAXGAP=1
SN_MAXLAT=60.0
SN_MINLEN=5
SN_MINENDPOINTDIST=2.0

rm -v cyc.${UQSTR} trajectories.txt.${UQSTR}
touch cyc.${UQSTR}
for f in ${FILES[@]}; do
  echo "Processing $f..."
  ${TEMPESTEXTREMESDIR}/bin/DetectNodes --verbosity 0 --timestride ${TIMESTRIDE} --in_data "${f}" ${CONNECTFLAG} --out cyc_tempest.${UQSTR} --closedcontourcmd "PSL,${DCU_PSLFOMAG},${DCU_PSLFODIST},0;_DIFF(Z300,Z500),${DCU_WCFOMAG},${DCU_WCFODIST},${DCU_WCMAXOFFSET}" --mergedist ${DCU_MERGEDIST} --searchbymin PSL --outputcmd "PSL,min,0;_VECMAG(UBOT,VBOT),max,2"
  cat cyc_tempest.${UQSTR} >> cyc.${UQSTR}
  rm -v cyc_tempest.${UQSTR}
done

${TEMPESTEXTREMESDIR}/bin/StitchNodes --format "ncol,lon,lat,slp,wind" --min_endpoint_dist ${SN_MINENDPOINTDIST} --range ${SN_TRAJRANGE} --minlength ${SN_MINLEN} --maxgap ${SN_TRAJMAXGAP} --in cyc.${UQSTR} --out ${TRAJFILE} --threshold "lat,<=,${SN_MAXLAT},4;lat,>=,-${SN_MAXLAT},4"

