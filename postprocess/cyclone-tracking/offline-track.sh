#!/bin/bash

## User settings
casename="ERA5-tclf002-ne30x4-L32"  # Betacast case
YYYYMMDDHH=2002100200                  # Init cycle
ATCFTECH="SE28L32"                     # This is the name for this model for ATCF (if ensemble, include member number)
hstream="h3i"                           # Which h stream contains the tracker data PSL, UBOT, VBOT
TIMESTRIDE=2                           # 2 for 3-hourly data, 1 for 6-hourly data (in h stream)
# Where are the output files archived after Betacast is done running?
OUTPUTBASE="/glade/derecho/scratch/zarzycki/${casename}/run/"

sendHTML=false
#connectfile="/glade/u/home/zarzycki/work/grids/connect/mpasa3-60-tclf002_connect.dat"
connectfile="/glade/u/home/zarzycki/tempest-scripts/hyperion/ne0np4natlanticref.ne30x4.connect_v2.dat"
TEMPESTEXTREMESDIR=/glade/work/zarzycki/tempestextremes_noMPI/

# Define some betacast specific stuff
ATCFFILE=atcf.${casename}.${YYYYMMDDHH}
TCVITFOLDER=./fin-tcvitals/

## First try and get observed vitals files for reference
./get-vitals.sh ${YYYYMMDDHH} ${TCVITFOLDER}

## Now track using tempest and match up tracked TCs with those in observations + put in ATCF format (for use in met tools, etc.)
./cyclone-tracking-driver.sh ${YYYYMMDDHH} \
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
