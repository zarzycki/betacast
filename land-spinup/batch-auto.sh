#!/bin/bash

#YYYYMMDD=19960113
YYYYMMDD=19961225
NMONTHS=60

./auto-script.sh 1 0 ${YYYYMMDD} ${NMONTHS} 1 1921
#./auto-script.sh 1 0 ${YYYYMMDD} ${NMONTHS} 1 -1
./auto-script.sh 1 0 ${YYYYMMDD} ${NMONTHS} 1 2019
./auto-script.sh 1 0 ${YYYYMMDD} ${NMONTHS} 1 2044
./auto-script.sh 1 0 ${YYYYMMDD} ${NMONTHS} 1 2064
./auto-script.sh 1 0 ${YYYYMMDD} ${NMONTHS} 1 2081

# MPAS stuff for 3km TCs
#./auto-script.sh 0 1 20050818 12 1 -1 nl.landspinup.mpas
#./auto-script.sh 0 1 19911122 12 1 -1 nl.landspinup.mpas
#./auto-script.sh 0 1 19890822 12 1 -1 nl.landspinup.mpas
#./auto-script.sh 0 1 20080831 12 1 -1 nl.landspinup.mpas
#./auto-script.sh 0 1 20000918 12 1 -1 nl.landspinup.mpas
#./auto-script.sh 0 1 20110806 12 1 -1 nl.landspinup.mpas
#./auto-script.sh 0 1 19950806 12 1 -1 nl.landspinup.mpas
#./auto-script.sh 0 1 20141013 12 1 -1 nl.landspinup.mpas
#./auto-script.sh 0 1 20130910 12 1 -1 nl.landspinup.mpas
#./auto-script.sh 0 1 20120825 12 1 -1 -1 nl.landspinup.mpas
#./auto-script.sh 0 1 20080818 12 1 -1 -1 nl.landspinup.mpas

# Testing Derecho
#./auto-script.sh 0 1 20000113 48 3 2080 1996 nl.landspinup.derecho
#./auto-script.sh 0 1 20000113 48 3 -1 -1 nl.landspinup.derecho
#./auto-script.sh 0 0 20050902 36 1 2080 1996 nl.landspinup.derecho
#./auto-script.sh 0 0 20050902 36 1 -1 -1 nl.landspinup.derecho
