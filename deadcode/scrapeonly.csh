#!/bin/bash

# Name of the job - You'll probably want to customize this.
#SBATCH -J scrape_only

# Standard out and Standard Error output files with the job number in the name.
#SBATCH -o CAM-%j.output
#SBATCH -e CAM-%j.error

#SBATCH -n 1

################################################################

## 1 = Yellowstone, 2 = Flux, 3 = Agri

machineid=3

gridname=tcforecast_60_x4
numHoursSEStart=3

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
outputdir=/home/zarzycki/tcforecast_60_x4/run
procdir=$outputdir/proc

mkdir -p $procdir
cd $outputdir

while [  $COUNTER -lt 1500 ]; do
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
	    ncl /home/zarzycki/ncl/sewx/weatherplot.ncl inisec=43200 iniday=26 inimon=07 iniyear=2013 'filename="'$f'"'
    done
    
    ## Trim whitespace around pngs
    genpngs=`ls *png`
    for g in $genpngs
    do
	    convert -trim +repage $g $g
    done
    
    
    ## Get file list in txt file for FLANIS viewer
    outflds=(wind tmq flut prect)

    # Loop over items in outflds
    for item in ${outflds[*]}
    do
        #printf "   %s\n" $item
        
        if [ ! -f ${item}files.txt ]; then
          echo "File not found!"
          ls -1 ${item}*png > ${item}files_nopath.txt
        else
          mv ${item}files_nopath.txt orig_${item}files.txt
          ls -1 ${item}*png > tocat_${item}files.txt
          cat orig_${item}files.txt tocat_${item}files.txt > ${item}files_nopath.txt
          rm tocat_${item}files.txt orig_${item}files.txt
        fi
        
        cp ${item}files_nopath.txt ${item}files.txt
        #sed -i.bak 's/^/##/' ${item}files.txt        
      
    done


    
    ## Move files to server
    ## Google create remote directory if not existant
    ## use rysnc?
    echo "Moving files to remote server"
    scp *.png *.txt balloflight@s483.sureserver.com:/home/balloflight/www/weather/current/
    
    mv *.png $procdir

    mv $newfiles $procdir

    let COUNTER=0
  fi
  echo Now sleeping $DTIME seconds
  sleep $DTIME
done

exit 0
