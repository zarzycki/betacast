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
numHoursSEStart=3

monthstr=`date -u +%d`
daystr=`date -u +%m`
yearstr=`date -u +%Y`
currtime=`date -u +%H%M`

if [ $currtime -lt 0330 ]
then
  echo "18 Z cycle"
  cyclestr=18
elif [ $currtime -lt 0930 ]
then
  echo "00 Z cycle"
  cyclestr=00
elif [ $currtime -lt 1530 ]
then
  echo "06 Z cycle"
  cyclestr=06
elif [ $currtime -lt 2130 ]
then
  echo "12 Z cycle"
  cyclestr=12
elif [ $currtime -ge 2130 ]
then
  echo "18 Z cycle"
  cyclestr=18
else
  echo "Can't figure out start time"
#  exit
fi

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
else
  echo "SE forecast lead time too long, 18Z cycle causes trouble"
  echo "Not yet supported."
#  exit
fi

monthstr=06
daystr=25
yearstr=2013
cyclestr=12
dayminusonestr=24

se_monthstr=monthstr
se_daystr=daystr
se_yearstr=yearstr
se_cyclestr=15
se_cyclestrsec=54000

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
    echo "Getting GFS conditions"
    cd $gfs_files_path
    
    ## Pull atmospheric conditions
    ## Need to write an if call depending on how far back the GFS files are located
    ## The NCEP server only holds them for a few days -- otherwise go to nomads at NCDC
    # gfsFtpPath=http://nomads.ncdc.noaa.gov/data/gfsanl/$yearstr$monstr/$yearstr$monstr$daystr/
    # gfs.t12z.pgrb2f00.grb2
    echo "Getting Atmo file"
    rm gfs.t*pgrb2f00*
    gfsFTPPath=ftp://ftp.ncep.noaa.gov/pub/data/nccf/com/gfs/prod/gfs.$yearstr$monthstr$daystr$cyclestr/
    gfsFTPFile='gfs.t'$cyclestr'z.pgrb2f00'
    wget $gfsFTPPath$gfsFTPFile
    mv $gfsFTPFile 'gfs_atm_'$yearstr$monthstr$daystr$cyclestr'.grib2'
    
    ## Pull sea surface temps, need to rename and delete (if necess)
    echo "Getting SST data"
    rm gdas1*sstgrb*
    sstFTPPath=ftp://ftp.ncep.noaa.gov/pub/data/nccf/com/gfs/prod/gdas.$yearstr$monthstr$dayminusonestr/
    sstFTPFile='gdas1.t'$cyclestr'z.sstgrb.grib2'
    wget $sstFTPPath$sstFTPFile
    mv $sstFTPFile 'gfs_sst_'$yearstr$monthstr$daystr$cyclestr'.grib2'
    
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

ncl se_interp.ncl machineid=$machineid 'gridname = "'$gridname'"' 'gfs_filename = "'$gfs_files_path'/gfs_atm_'$yearstr$monthstr$daystr$cyclestr'.grib2"'
ncl sst_interp.ncl 'sst_file_full = "'$gfs_files_path'/gfs_sst_'$yearstr$monthstr$daystr$cyclestr'.grib2"'

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
./xmlchange -file env_run.xml -id STOP_N -val 10
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
