#!/bin/bash

# Name of the job - You'll probably want to customize this.
#SBATCH -J sewx_driver

# Standard out and Standard Error output files with the job number in the name.
#SBATCH -o CAM-%j.output
#SBATCH -e CAM-%j.error

#SBATCH -n 1

################################################################

## 1 = Yellowstone, 2 = Flux, 3 = Agri

machineid=3

gridname=tcforecast_68_x4
numHoursSEStart=3

## Here we get two digit strings for UTC time for month, day, year
## We also get current time in hoursminutes (because the GFS output lags by 3.5 hours)
monthstr=`date -u +%m`
daystr=`date -u +%d`
yearstr=`date -u +%Y`
currtime=`date -u +%H%M`

## Use currtime to figure out what is the latest cycle we have access to
if [ $currtime -lt 0330 ]
then
  echo "18Z cycle"
  monthstr=`date --date="yesterday" -u +%m`
  daystr=`date --date="yesterday" -u +%d`
  yearstr=`date --date="yesterday" -u +%Y`
  cyclestr=18
elif [ $currtime -lt 0930 ]
then
  echo "00Z cycle"
  cyclestr=00
elif [ $currtime -lt 1530 ]
then
  echo "06Z cycle"
  cyclestr=06
elif [ $currtime -lt 2130 ]
then
  echo "12Z cycle"
  cyclestr=12
elif [ $currtime -ge 2130 ]
then
  echo "18Z cycle"
  cyclestr=18
else
  echo "Can't figure out start time"
#  exit
fi

## Figure out the seconds which correspond to the cycle and zero pad if neces
let cyclestrsec=$cyclestr*3600
while [ ${#cyclestrsec} -lt 5 ];
do
  cyclestrsec="0"$cyclestrsec
done

## Figure out what the SE start time will be after filter
if [ $numHoursSEStart -lt 6 ]
then
  let se_cyclestr=$cyclestr+03
  while [ ${#se_cyclestr} -lt 2 ];
  do
    se_cyclestr="0"$se_cyclestr
  done
  echo $se_cyclestr

  se_monthstr=$monthstr
  se_daystr=$daystr
  se_yearstr=$yearstr
  let se_cyclestrsec=$se_cyclestr*3600
  while [ ${#se_cyclestrsec} -lt 5 ];
  do
    se_cyclestrsec="0"$se_cyclestrsec
  done
else
  echo "SE forecast lead time too long, 18Z cycle causes trouble"
  echo "Not yet supported."
  exit 1
fi

## Use currtime to figure out what SST we can download
## Current GDAS SST appears 555 after cycle time
## First guess is today's dates!
sstmonthstr=$monthstr
sstdaystr=$daystr
sstyearstr=$yearstr
if [ $currtime -lt 0555 ]
then
  echo "SSTs are previous days 18Z"
  sstmonthstr=`date --date="yesterday" -u +%m`
  sstdaystr=`date --date="yesterday" -u +%d`
  sstyearstr=`date --date="yesterday" -u +%Y`
  sstcyclestr=18
elif [ $currtime -lt 1155 ]
then
  echo "SSTs are current days 00Z"
  sstcyclestr=00
elif [ $currtime -lt 1755 ]
then
  echo "SSTs are current days 06Z"
  sstcyclestr=06
elif [ $currtime -lt 2355 ]
then
  echo "SSTs are current days 12Z"
  sstcyclestr=12
elif [ $currtime -ge 2355 ]
then
  echo "SSTs are current days 18Z"
  sstcyclestr=18
else
  echo "Can't figure out start time"
  exit 1
fi

echo "The current time is $currtime"
echo "We are using $yearstr $monthstr $daystr $cyclestr Z ($cyclestrsec seconds) for GFS data"
echo "We are using $sstyearstr $sstmonthstr $sstdaystr $sstcyclestr Z for SST data"
echo "SE initialization will occur at $se_yearstr $se_monthstr $se_daystr $se_cyclestr Z ($se_cyclestrsec seconds)"

# Set input data type
# 1 = GFS forecast ICs
# 2 = ECMWF ERA40 reanalysis
inputdatatype=1

if [ $machineid -eq 1 ]
then
  echo "Using Yellowstone"
  backupdir=/glade/u/home/$LOGNAME/scriptbackups
  gfs_files_path=/glade/u/home/$LOGNAME/sewx/GFS
  gfs_to_cam_path=/glade/u/home/$LOGNAME/sewx/gfs_to_cam
  era_to_cam_path=/glade/home/zarzycki/obs_data_init/reanal_to_cam
  filter_path=/glade/u/home/$LOGNAME/sewx/filter
  path_to_se_build=/glade/u/home/$LOGNAME
  path_to_gz_touch=/glade/scratch/$LOGNAME/archive/$gridname/atm/logs
  path_to_nc_files=/glade/scratch/$LOGNAME/archive/$gridname/atm/hist
elif [ $machineid -eq 2 ]
then
  echo "Using UMich Flux"
  exit 1
elif [ $machineid -eq 3 ]
then
  echo "Using UCDavis Agri"
  backupdir=/home/$LOGNAME/scriptbackups
  gfs_files_path=/home/$LOGNAME/sewx/GFS
  gfs_to_cam_path=/home/$LOGNAME/sewx/gfs_to_cam
  era_to_cam_path=/home/$LOGNAME/sewx/reanal_to_cam
  filter_path=/home/$LOGNAME/sewx/filter
  path_to_se_build=/home/$LOGNAME
  path_to_gz_touch=$path_to_se_build/$gridname/run
  path_to_nc_files=$path_to_se_build/$gridname/run
else
  echo "Machine"
fi

#### Now scrape for NCL files

COUNTER=0
MAXCOUNTER=900
DTIME=90
outputdir=/home/zarzycki/tcforecast_68_x4/run
procdir=$outputdir/proc

mkdir -p $procdir
cd $outputdir

while [  $COUNTER -lt 900 ]; do
  echo The starting counter is $COUNTER
  filenames=`ls tc*h1*.nc`
  numfiles=`ls tc*h1*.nc | wc -l`
  echo $numfiles
  if [ $numfiles -eq 0 ]
  then
    echo "Nothing found"
    let COUNTER=COUNTER+DTIME
  else
    sleep 10
    echo "Found at least one file"
    rename 's/^/_/' tc*h1*.nc
    newfiles=`ls _tc*h1*.nc`
    for f in $newfiles
    do
	    echo "Processing"
	    ncl /home/zarzycki/ncl/sewx/weatherplot.ncl inisec=$cyclestrsec iniday=$daystr inimon=$monthstr iniyear=$yearstr 'filename="'$f'"'
    done
    mv $newfiles $procdir
    let COUNTER=0
  fi
  echo Now sleeping $DTIME seconds
  sleep $DTIME
done

exit 0
