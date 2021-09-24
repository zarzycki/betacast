#!/bin/bash

## User settings
casename=t35_v1.2-NATL-F2010C5-v2_L71
YYYYMMDDHH=2011082512
betacastdir=~/betacast/
connectfile="/global/homes/c/czarzyck/betacast/cyclone-tracking/ne0np4natlanticref.ne30x4.connect_v2.dat"
OUTPUTBASE="/global/cscratch1/sd/jjbenedi/e3sm_scratch/cori-knl/"
sendHTML=false

# Define some betacast specific stuff
ATCFFILE=atcf.${casename}.${YYYYMMDDHH}
TCVITFOLDER=./fin-tcvitals/

## First try and get observed vitals files for reference
./get-vitals.sh ${YYYYMMDDHH} ${TCVITFOLDER} ${betacastdir}

## Now track
./drive-tracking.sh ${YYYYMMDDHH} ${casename} ${TCVITFOLDER}/tcvitals.${YYYYMMDDHH} ${ATCFFILE} ${connectfile} ${OUTPUTBASE} ${sendHTML}