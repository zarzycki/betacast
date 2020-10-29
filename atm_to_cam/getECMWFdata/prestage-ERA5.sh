#!/bin/bash
#python3 must be loaded with cdsapi package

yearstr=1996
monthstr=01
daystr=15
cyclestr=12

set -e
python getERA5.py ${yearstr}${monthstr}${daystr} ${cyclestr}
ncrename -v z,phis ERA5_sfc_${yearstr}${monthstr}${daystr}${cyclestr}.nc
ncks -A ERA5_ml_${yearstr}${monthstr}${daystr}${cyclestr}.nc ERA5_sfc_${yearstr}${monthstr}${daystr}${cyclestr}.nc
mv -v ERA5_sfc_${yearstr}${monthstr}${daystr}${cyclestr}.nc ERA5_${yearstr}${monthstr}${daystr}${cyclestr}.nc
rm -f ERA5_ml_${yearstr}${monthstr}${daystr}${cyclestr}.nc
