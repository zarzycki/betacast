#!/bin/bash
# python3 must be loaded with cdsapi package
# usage: ./prestage-ERA5.sh 1996011800
# where you want the data file for Jan 18, 1996 at 00Z

DATE=$1

yearstr=${DATE:0:4}
monthstr=${DATE:4:2}
daystr=${DATE:6:2}
cyclestr=${DATE:8:2}

set -e
python getERA5.py ${yearstr}${monthstr}${daystr} ${cyclestr}
ncrename -v z,phis ERA5_sfc_${yearstr}${monthstr}${daystr}${cyclestr}.nc
ncks -A ERA5_ml_${yearstr}${monthstr}${daystr}${cyclestr}.nc ERA5_sfc_${yearstr}${monthstr}${daystr}${cyclestr}.nc
mv -v ERA5_sfc_${yearstr}${monthstr}${daystr}${cyclestr}.nc ERA5_${yearstr}${monthstr}${daystr}${cyclestr}.nc
rm -f ERA5_ml_${yearstr}${monthstr}${daystr}${cyclestr}.nc
