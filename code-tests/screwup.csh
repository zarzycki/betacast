#!/bin/bash

# Name of the job - You'll probably want to customize this.
#SBATCH -J sewx_driver

# Standard out and Standard Error output files with the job number in the name.
#SBATCH -o CAM-%j.output
#SBATCH -e CAM-%j.error

#SBATCH -n 1
#SBATCH --exclusive

################################################################

## 1 = Yellowstone, 2 = Flux, 3 = Agri

debug=0

machineid=3

gridname=newtcforecast_60_x4
numHoursSEStart=3

date

## Here we get two digit strings for UTC time for month, day, year
## We also get current time in hoursminutes (because the GFS output lags by 3.5 hours)
monthstr=`date -u +%m`
daystr=`date -u +%d`
yearstr=`date -u +%Y`
currtime=`date -u +%H%M`
machzone=`date +%z`
twodaysago=`date --date='2 days ago' -u +"%Y%m%d"`

## Use currtime to figure out what is the latest cycle we have access to
if [ $currtime -lt 0328 ]
then
  echo "18Z cycle"
  monthstr=`date --date="yesterday" -u +%m`
  daystr=`date --date="yesterday" -u +%d`
  yearstr=`date --date="yesterday" -u +%Y`
  twodaysago=`date --date='3 days ago' -u +"%Y%m%d"`
  cyclestr=18
elif [ $currtime -lt 0928 ]
then
  echo "00Z cycle"
  cyclestr=00
elif [ $currtime -lt 1528 ]
then
  echo "06Z cycle"
  cyclestr=06
elif [ $currtime -lt 2128 ]
then
  echo "12Z cycle"
  cyclestr=12
elif [ $currtime -ge 2128 ]
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
  sewxscriptsdir=/glade/u/home/$LOGNAME/sewx/
  backupdir=/glade/u/home/$LOGNAME/scriptbackups
  gfs_files_path=/glade/u/home/$LOGNAME/sewx/GFS
  gfs_to_cam_path=/glade/u/home/$LOGNAME/sewx/gfs_to_cam
  era_to_cam_path=/glade/home/zarzycki/obs_data_init/reanal_to_cam
  filter_path=/glade/u/home/$LOGNAME/sewx/filter
  path_to_se_build=/glade/u/home/$LOGNAME
  path_to_gz_touch=/glade/scratch/$LOGNAME/archive/$gridname/atm/logs
  path_to_nc_files=/glade/scratch/$LOGNAME/archive/$gridname/atm/hist
  outputdir=/glade/scratch/$LOGNAME/${gridname}/run
elif [ $machineid -eq 2 ]
then
  echo "Using UMich Flux"
  exit 1
elif [ $machineid -eq 3 ]
then
  echo "Using UCDavis Agri"
  sewxscriptsdir=/home/$LOGNAME/sewx/
  backupdir=/home/$LOGNAME/scriptbackups
  gfs_files_path=/home/$LOGNAME/sewx/GFS
  gfs_to_cam_path=/home/$LOGNAME/sewx/gfs_to_cam
  era_to_cam_path=/home/$LOGNAME/sewx/reanal_to_cam
  filter_path=/home/$LOGNAME/sewx/filter
  path_to_se_build=/home/$LOGNAME
  path_to_gz_touch=$path_to_se_build/$gridname/run
  path_to_nc_files=$path_to_se_build/$gridname/run
  outputdir=/home/$LOGNAME/${gridname}/run
else
  echo "Machine"
fi

# 30 -> CAM5 physics
# 26 -> CAM4 physics
numlevels=30

##set string=$month$day
##echo $string

  echo "Done sleeping, time to start scraping for files to plot"
  
  #### Now scrape for NCL files
      
  ## rm -rf `ls -tr | head -n 1`
  #ssh balloflight@s483.sureserver.com 'rm /home/balloflight/www/weather/current/*.png ; rm /home/balloflight/www/weather/current/*.txt'
  echo "SSHings..."
  ssh balloflight@s483.sureserver.com "mkdir -p /home/balloflight/www/weather/current/${yearstr}${monthstr}${daystr}${cyclestr} ; \ 
       rm -rf /home/balloflight/www/weather/current/${twodaysago}${cyclestr} ; \
       cd /home/balloflight/www/weather/current/ ; \
       cp *cfg *html ${yearstr}${monthstr}${daystr}${cyclestr} ; \
       rm ${yearstr}${monthstr}${daystr}${cyclestr}/index.html "
  
  COUNTER=0
  MAXCOUNTER=630
  DTIME=270
  procdir=$outputdir/proc
  
  mkdir -p $procdir
  mkdir -p $procdir/images
  mkdir -p $procdir/text
  mkdir -p $procdir/nl_files
  mkdir -p $procdir/logs
  cd $outputdir
  
  ## UPDATE html page
  if [ ! -f index.html ]; then
    cp index.HOLD index.html
    sed '/<!--FORECASTHEAD-->/ r htmltemplate.html' index.html > _index1.html
    sed -e "/$twodaysago${cyclestr}/d" _index1.html > _index2.html
    sed -e "s/YYYYMMDDHH/${yearstr}${monthstr}${daystr}${cyclestr}/" _index2.html > _index3.html
    mv _index3.html index.html
    rm _index*.html
    scp index.html balloflight@s483.sureserver.com:/home/balloflight/www/weather/current
  fi
  
  while [  $COUNTER -lt 630 ]; do
    echo The starting counter is $COUNTER
    filenames=`ls newtc*h1*.nc`
    numfiles=`ls newtc*h1*.nc | wc -l`
    echo $numfiles
    if [ $numfiles -eq 0 ]
    then
      echo "Nothing found"
      let COUNTER=COUNTER+DTIME
    else
      sleep 8
      echo "Found at least one file"
      rename 's/^/_/' newtc*h1*.nc
      newfiles=`ls _newtc*h1*.nc`
      for f in $newfiles
      do
        echo "Processing"
        ncl /home/zarzycki/ncl/sewx/weatherplot.ncl inisec=$cyclestrsec iniday=$daystr inimon=$monthstr iniyear=$yearstr 'filename="'$f'"' > ncl.output 2>&1
        if [ grep FileReadVar ncl.output ]; then
          sleep 8
          echo "Found an error"
          ncl /home/zarzycki/ncl/sewx/weatherplot.ncl inisec=$cyclestrsec iniday=$daystr inimon=$monthstr iniyear=$yearstr 'filename="'$f'"'
        fi
        rm ncl.output
      done
      
      ## Trim whitespace around pngs
      genpngs=`ls *png`
      for g in $genpngs
      do
        convert -trim +repage $g $g
      done
      
      ## Get file list in txt file for FLANIS viewer
      basins=(natl epac)
      outflds=(wind tmq flut prect 500vort shear850250)
  
      # Loop over items in outflds
      for basin in ${basins[*]}
      do
        for item in ${outflds[*]}
        do
            #printf "   %s\n" $item
            
            if [ ! -f ${item}_${basin}_files.txt ]; then
              echo "File not found!"
              ls -1 ${item}*${basin}*png > ${item}_${basin}_files_nopath.txt
            else
              mv ${item}_${basin}_files_nopath.txt orig_${item}_${basin}_files.txt
              ls -1 ${item}*${basin}*png > tocat_${item}_${basin}_files.txt
              cat orig_${item}_${basin}_files.txt tocat_${item}_${basin}_files.txt > ${item}_${basin}_files_nopath.txt
              rm tocat_${item}_${basin}_files.txt orig_${item}_${basin}_files.txt
            fi
            
            cp ${item}_${basin}_files_nopath.txt ${item}_${basin}_files.txt
            #sed -i.bak 's/^/##/' ${item}files.txt        
          
        done
      done
      
      ## Move files to server
      ## Google create remote directory if not existant
      ## use rysnc?
      echo "Moving files to remote server"
      scp *.png *.txt balloflight@s483.sureserver.com:/home/balloflight/www/weather/current/${yearstr}${monthstr}${daystr}${cyclestr}
      
      mv *.png $procdir/images
  
      mv $newfiles $procdir
  
      let COUNTER=0
    fi
    echo Now sleeping $DTIME seconds
    sleep $DTIME
  done

fi

### Finish uploading index.html
printtime=`date -u`
sed -e 's/\"red/\"green/' index.html > _index4.html
sed -e "s/CURRENTLY UPDATING/COMPLETED AT $printtime/" _index4.html > _index5.html
mv _index5.html index.html
rm _index*html
scp index.html balloflight@s483.sureserver.com:/home/balloflight/www/weather/current
mv index.html index.HOLD

## Move other cam files to dir
mv *.cam.h*.nc $procdir
cp *_in seq_maps.rc *_modelio.nml docn.streams.txt.prescribed $procdir/nl_files
mv *.txt $procdir/text
mv *.log.* $procdir/logs
mv timing/ $procdir/

## Move land files to new restart location
echo "Moving land restart files"
cd $path_to_nc_files
landdir=$path_to_nc_files/clmstart
mkdir $landdir 
mv *.clm2.r.*nc $landdir

echo "Deleting boring restart files"
# Delete boring restart files
cd $path_to_nc_files
rm *.clm2.rh0.*.nc
rm *.docn.rs1.*.bin
rm *.cam.r.*.nc
rm *.cam.rs.*.nc
rm *.cpl.r.*.nc
rm *.clm2.rh0.*.nc
rm *.cam.rh3.*.nc
rm *.cice.r.*.nc
rm rpointer.*

cp /home/zarzycki/sewx/INIC/${gridname}_INIC_filter.nc $procdir
cp /home/zarzycki/sewx/SST/sst_1x1.nc $procdir

mv $procdir $outputdir/${yearstr}${monthstr}${daystr}${cyclestr}

date

cd $sewxscriptsdir

### Still need to figure out how to convert machzone properly
if [ $machzone -eq -0000 ]
then
  machzone=0
elif [ $machzone -eq -0300 ]
then
  machzone=-300
elif [ $machzone -eq -0400 ]
then
  machzone=-400
elif [ $machzone -eq -0500 ]
then
  machzone=-500
elif [ $machzone -eq -0600 ]
then
  machzone=-600
elif [ $machzone -eq -0700 ]
then
  machzone=-700
elif [ $machzone -eq -0800 ]
then
  machzone=-800
elif [ $machzone -eq -0900 ]
then
  machzone=-900
else
  echo "Can't figure out machzone, you probably need to manually add it"
  echo "at least until Colin figures out bash arithmetic"
fi

## 
if [ $cyclestr -eq 00 ]
then
  let "resubhour = 1500 + $machzone"
elif [ $cyclestr -eq 06 ]
then
  let "resubhour = 2100 + $machzone"
elif [ $cyclestr -eq 12 ]
then
  let "resubhour = 300 + $machzone"
elif [ $cyclestr -eq 18 ]
then
  let "resubhour = 900 + $machzone"
else
  echo "Can't figure out start time"
fi

if [ $resubhour -lt 0 ]
then
  let "resubhour = 2400 + $resubhour"
fi
let "resubhour = $resubhour/100"

echo "The resub hour is: $resubhour"
echo "Therefore restart time is: ${resubhour}:30"

#Sed into testshell.csh
sed -i 's?^#SBATCH --begin.*?#SBATCH --begin='$resubhour':30?' testshell.csh

sbatch testshell.csh

exit 0

