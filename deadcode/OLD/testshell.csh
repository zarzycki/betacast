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

gridname=atlsquad_30_x8
monthstr=10
daystr=24
yearstr=2012
cyclestr=00

se_monthstr=10
se_daystr=24
se_yearstr=2012
se_cyclestr=03
se_cyclestrsec=10800

# Set input data type
# 1 = GFS forecast ICs
# 2 = ECMWF ERA40 reanalysis
inputdatatype=1

if [ $machineid -eq 1 ]
then
  echo "Using Yellowstone"
  backupdir=/glade/u/home/$LOGNAME/scriptbackups
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
  gfs_to_cam_path=/home/$LOGNAME/sewx/gfs_to_cam
  era_to_cam_path=/home/$LOGNAME/sewx/reanal_to_cam
  filter_path=/home/$LOGNAME/sewx/filter
  path_to_se_build=/home/$LOGNAME
  path_to_gz_touch=$path_to_se_build/$gridname/run
  path_to_nc_files=$path_to_se_build/$gridname/run
else
  echo "Unsupported start time"
fi

if [ $cyclestr -eq 00 ]
then
  cyclestrsec=00000
  echo $cyclestrsec
elif [ $cyclestr -eq 06 ]
then
  cyclestrsec=21600
  echo $cyclestrsec
elif [ $cyclestr -eq 12 ]
then
  cyclestrsec=43200
  echo $cyclestrsec
elif [ $cyclestr -eq 18 ]
then
  cyclestrsec=64800
  echo $cyclestrsec
else
  echo "Unsupported start time"
fi

# 30 -> CAM5 physics
# 26 -> CAM4 physics
numlevels=30

##set string=$month$day
##echo $string

# Set timestamp for backing up files, etc.
timestamp=`date +%Y%m%d.%H%M`
mkdir -p $backupdir

if [ $inputdatatype -eq 1 ]
then  
    echo "Using GFS forecast ICs"
    echo "Cding to GFS interpolation directory"
    cd $gfs_to_cam_path 
elif [$inputdatatype -eq 2 ]
then  
    echo "Using ECMWF forecast ICs"
    echo "Cding to ECMWF interpolation directory"
    cd $era_to_cam_path
else
    echo "Incorrect model IC entered"
    exit 1
fi

echo "Sedding commands in NCL interpolation scripts"
cp se_interp.ncl $backupdir/se_interp.BAK.$timestamp
sed 's?^numlevels.*?numlevels='$numlevels'?' se_interp.ncl > test.ncl
mv test.ncl se_interp.ncl
sed 's?^monthstr.*?monthstr="'$monthstr'"?' se_interp.ncl > test.ncl
mv test.ncl se_interp.ncl
sed 's?^yearstr.*?yearstr="'$yearstr'"?' se_interp.ncl > test.ncl
mv test.ncl se_interp.ncl
sed 's?^daystr.*?daystr="'$daystr'"?' se_interp.ncl > test.ncl
mv test.ncl se_interp.ncl
sed 's?^cyclestr.*?cyclestr="'$cyclestr'"?' se_interp.ncl > test.ncl
mv test.ncl se_interp.ncl

##ncl se_interp.ncl machineid=$machineid 'gridname = "'$gridname'"'
## ncl sst_interp.ncl

cd $path_to_se_build/$gridname

echo "Change date in FV namelist"
./xmlchange -file env_run.xml -id RUN_STARTDATE -val $yearstr-$monthstr-$daystr
./xmlchange -file env_run.xml -id START_TOD -val $cyclestrsec
./xmlchange -file env_run.xml -id STOP_OPTION -val nhours
./xmlchange -file env_run.xml -id STOP_N -val 6
cp user_nl_cam_filter user_nl_cam

## Get number of log files archived
numlogfiles=`ls $path_to_gz_touch/*.gz | wc -l`

echo $numlogfiles

echo "Begin call to filter-run"
if [ $machineid -eq 1 ]
then
  bsub < $gridname.run
elif [ $machineid -eq 2 ]
then
  echo "Using UMich Flux"
  exit 1
elif [ $machineid -eq 3 ]
then
  sbatch $gridname.run
  echo "Skip me"
else
  echo "Unsupported start time"
fi

## Hold script while log files from filter run haven't been archived yet
while [ `ls $path_to_gz_touch/*.gz | wc -l` == $numlogfiles ]
do
  sleep 20
  echo "Sleeping"
done
echo "Run over done sleeping will hold for 60 more sec to make sure files moved"

sleep 60

## Run NCL filter
cd $filter_path
filtfile_name=$gridname.cam.h1.$yearstr-$monthstr-$daystr-$cyclestrsec.nc
ncl digifilter.ncl machineid=$machineid 'filtfile_name = "'$filtfile_name'"' 'gridname = "'$gridname'"'

echo "Removing filter files"
mkdir -p $path_to_nc_files/filtered
mv $path_to_nc_files/*.nc $path_to_nc_files/filtered

echo "Make changes in CESM-SE namelist"
cd $path_to_se_build/$gridname
./xmlchange -file env_run.xml -id RUN_STARTDATE -val $se_yearstr-$se_monthstr-$se_daystr
./xmlchange -file env_run.xml -id START_TOD -val $se_cyclestrsec
./xmlchange -file env_run.xml -id STOP_OPTION -val ndays
./xmlchange -file env_run.xml -id STOP_N -val 2
cp user_nl_cam_run user_nl_cam

echo "Begin call to forecast run"
if [ $machineid -eq 1 ]
then
  bsub < $gridname.run
elif [ $machineid -eq 2 ]
then
  echo "Using UMich Flux"
  exit 1
elif [ $machineid -eq 3 ]
then
  sbatch $gridname.run
else
  echo "Unsupported start time"
fi

exit 0