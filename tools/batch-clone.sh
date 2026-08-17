#!/bin/bash

FARRAY=( "RoS-F2010CICEROF-ne0wus30x32-301-control" )

IARRAY=( "RoS-ICLM45-ne0wus30x32-301_19961225_0060" )

# FARRAY=( "RoS-F2010CICEROF-ne0wus30x32-101-plus1K" \
#   "RoS-F2010CICEROF-ne0wus30x32-101-plus2K" \
#   "RoS-F2010CICEROF-ne0wus30x32-101-plus3K" \
#   "RoS-F2010CICEROF-ne0wus30x32-101-plus4K" \
#   "RoS-F2010CICEROF-ne0wus30x32-101-PI" )
#
# IARRAY=( "RoS-ICLM45-ne0wus30x32-201_19961225_0060_2019" \
#   "RoS-ICLM45-ne0wus30x32-201_19961225_0060_2044" \
#   "RoS-ICLM45-ne0wus30x32-201_19961225_0060_2064" \
#   "RoS-ICLM45-ne0wus30x32-201_19961225_0060_2081" \
#   "RoS-ICLM45-ne0wus30x32-201_19961225_0060_1921" )

# FARRAY=( "RoS-F2010C5-ne0conus30x8-101-plus1K" \
#   "RoS-F2010C5-ne0conus30x8-101-plus2K" \
#   "RoS-F2010C5-ne0conus30x8-101-plus3K" \
#   "RoS-F2010C5-ne0conus30x8-101-plus4K" \
#   "RoS-F2010C5-ne0conus30x8-101-PI" )
#
# IARRAY=( "RoS-ICLM45-ne0conus30x8-101_19961225_0036_2019" \
#   "RoS-ICLM45-ne0conus30x8-101_19961225_0036_2044" \
#   "RoS-ICLM45-ne0conus30x8-101_19961225_0036_2064" \
#   "RoS-ICLM45-ne0conus30x8-101_19961225_0036_2081" \
#   "RoS-ICLM45-ne0conus30x8-101_19961225_0036_1921" )

# FARRAY=( "RoS-F2010C5-ne0conus30x8-006-control" \
#   "RoS-F2010C5-ne0conus30x8-006-plus1K" \
#   "RoS-F2010C5-ne0conus30x8-006-plus2K" \
#   "RoS-F2010C5-ne0conus30x8-006-plus3K" \
#   "RoS-F2010C5-ne0conus30x8-006-plus4K" \
#   "RoS-F2010C5-ne0conus30x8-006-PI" )
#
# IARRAY=( "RoS-ICLM45-ne0conus30x8-006_19960113_0012" \
#   "RoS-ICLM45-ne0conus30x8-006_19960113_0012_2019" \
#   "RoS-ICLM45-ne0conus30x8-006_19960113_0012_2044" \
#   "RoS-ICLM45-ne0conus30x8-006_19960113_0012_2064" \
#   "RoS-ICLM45-ne0conus30x8-006_19960113_0012_2081" \
#   "RoS-ICLM45-ne0conus30x8-006_19960113_0012_1921" )

#MASTERCASENAME=RoS-F2010C5-ne0conus30x8-101-control
MASTERCASENAME=RoS-F2010CICEROF-ne0wus30x32-101-control
E3SMPATH=/global/homes/c/czarzyck/E3SM-20220616/

BETACAST=/global/homes/c/czarzyck/betacast/
PATHTOCASE=~/F-compsets

for i in $(seq 0 $(( ${#FARRAY[@]} - 1 )))
do
  echo "----"
  echo ${FARRAY[$i]}" <----> "${IARRAY[$i]}

  CASENAME=${FARRAY[$i]}
  cd ${E3SMPATH}/cime/scripts
  ./create_clone --case ${PATHTOCASE}/${CASENAME} --clone ${PATHTOCASE}/${MASTERCASENAME} --keepexe
  cd ${PATHTOCASE}/${CASENAME}
  ./case.setup --silent
  ./case.build --silent

  echo "now doing land stuff"

  LANDSPINUPDIR=${IARRAY[$i]}
  cd ~/scratch/e3sm_scratch/cori-knl/${CASENAME}/run
  mkdir landstart
  cd landstart
  cp -v ~/scratch/e3sm_scratch/cori-knl/${LANDSPINUPDIR}/run/*.elm.r.*.nc /global/homes/c/czarzyck/m2637/betacast/elm-inic/1997wus
  cp -v ~/scratch/e3sm_scratch/cori-knl/${LANDSPINUPDIR}/run/*.mosart.r.*.nc /global/homes/c/czarzyck/m2637/betacast/elm-inic/1997wus
  cp -v ~/scratch/e3sm_scratch/cori-knl/${LANDSPINUPDIR}/run/*.elm.r.*.nc .
  cp -v ~/scratch/e3sm_scratch/cori-knl/${LANDSPINUPDIR}/run/*.mosart.r.*.nc .
  rename -v ${LANDSPINUPDIR} ${CASENAME} *nc

  cd ${BETACAST}

  echo "****** DONE!"
done

