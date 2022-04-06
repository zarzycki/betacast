#!/bin/bash

## User settings
casename=t37_v1.2-NATL-F2010C5-v2_L71  # Betacast case
YYYYMMDDHH=2011082512                  # Init cycle
ATCFTECH="EAM001"                      # This is the name for this model for ATCF (if ensemble, include member number)
hstream="h1"                           # Which h stream contains the tracker data PSL, UBOT, VBOT
TIMESTRIDE=2                           # 2 for 3-hourly data, 1 for 6-hourly data (in h stream)
# Where are the output files archived after Betacast is done running?
OUTPUTBASE="/global/cscratch1/sd/jjbenedi/e3sm_scratch/cori-knl/${casename}.${YYYYMMDDHH}.ens001/run/"

sendHTML=false
connectfile="/global/homes/c/czarzyck/betacast/cyclone-tracking/ne0np4natlanticref.ne30x4.connect_v2.dat"
TEMPESTEXTREMESDIR=/global/homes/c/czarzyck/sw/tempestextremes_noMPI/
 
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
