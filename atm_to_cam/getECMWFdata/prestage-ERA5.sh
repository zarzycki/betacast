#!/bin/bash -l
# python3 must be loaded with cdsapi package
# on cheyenne...
# $> module load python
# $> ncar_pylib
#
# usage: ./prestage-ERA5.sh /glade/work/${LOGNAME}/sewx/ECMWF/ 1996011800
# where you want the data file for Jan 18, 1996 at 00Z

conda activate meteo

DATE=$2
DIR=$1
THISDIR=$PWD

yearstr=${DATE:0:4}
monthstr=${DATE:4:2}
daystr=${DATE:6:2}
cyclestr=${DATE:8:2}

set -e
mkdir -p $DIR
cd $DIR
python $THISDIR/getERA5.py ${yearstr}${monthstr}${daystr} ${cyclestr}
ncrename -v z,phis ERA5_sfc_${yearstr}${monthstr}${daystr}${cyclestr}.nc
ncks -A ERA5_ml_${yearstr}${monthstr}${daystr}${cyclestr}.nc ERA5_sfc_${yearstr}${monthstr}${daystr}${cyclestr}.nc
mv -v ERA5_sfc_${yearstr}${monthstr}${daystr}${cyclestr}.nc ERA5_${yearstr}${monthstr}${daystr}${cyclestr}.nc
rm -f ERA5_ml_${yearstr}${monthstr}${daystr}${cyclestr}.nc
