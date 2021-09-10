#!/bin/bash

casename=RHEL7_forecast_nhemitc_30_x4_CAM5_L30
YYYYMMDDHH=2021082100
ATCFFILE=atcf.${casename}.${YYYYMMDDHH}

./drive-tracking.sh ${YYYYMMDDHH} ${casename} ./fin-tcvitals/tcvitals.${YYYYMMDDHH} ${ATCFFILE}

