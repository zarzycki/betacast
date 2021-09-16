#!/bin/bash

casename=RHEL7_forecast_nhemitc_30_x4_CAM5_L30.002
YYYYMMDDHH=2021082200
ATCFFILE=atcf.${casename}.${YYYYMMDDHH}

./drive-tracking.sh ${YYYYMMDDHH} ${casename} ./fin-tcvitals/tcvitals.${YYYYMMDDHH} ${ATCFFILE}

