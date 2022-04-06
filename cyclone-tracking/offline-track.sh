#!/bin/bash

## User settings
casename=t37_v1.2-NATL-F2010C5-v2_L71
YYYYMMDDHH=2011082512
connectfile="/global/homes/c/czarzyck/betacast/cyclone-tracking/ne0np4natlanticref.ne30x4.connect_v2.dat"
OUTPUTBASE="/global/cscratch1/sd/jjbenedi/e3sm_scratch/cori-knl/${casename}.${YYYYMMDDHH}.ens001/run/"
sendHTML=false
hstream="h1"
TEMPESTEXTREMESDIR=~/sw/tempestextremes_noMPI/
TIMESTRIDE=2         # should be 2 for 3-hourly data, 1 for 6-hourly data
ATCFTECH="EAM001"
 
# Define some betacast specific stuff
ATCFFILE=atcf.${casename}.${YYYYMMDDHH}
TCVITFOLDER=./fin-tcvitals/

## First try and get observed vitals files for reference
./get-vitals.sh ${YYYYMMDDHH} ${TCVITFOLDER}

## Now track using tempest and match up tracked TCs with those in observations + put in ATCF format (for use in met tools, etc.)
./drive-tracking.sh ${YYYYMMDDHH} \
    ${casename} \
    ${TCVITFOLDER}/tcvitals.${YYYYMMDDHH} \
    ${ATCFFILE} \
    ${connectfile} \
    ${OUTPUTBASE} \
    ${sendHTML} \
    ${hstream} \
    ${TIMESTRIDE} \
    ${ATCFTECH} \
    ${TEMPESTEXTREMESDIR}