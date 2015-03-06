#!/bin/bash

##=======================================================================
#BSUB -a poe                     # use LSF openmp elim
#BSUB -N
#BSUB -n 2                      # yellowstone setting
#BSUB -o out.%J                  # output filename
#BSUB -e out.%J                  # error filename
#BSUB -q geyser                 # queue
#BSUB -J sewx_driver
#BSUB -W 23:59                    # wall clock limit
#BSUB -P P54048000               # account number
##=======================================================================

set -e
#set -v

debug=0     # 0 = no debug, 1 = debug
echo "debug set to $debug"

# Set whether we are running forecast in realtime
# 0 no, we are using historical data --- 1 yes, we are running live
islive=0
echo "islive set to $islive"

# Set whether we want to send plots to external server
sendplots=0         # 0 no, 1 yes
echo "sendplots set to $sendplots"

## 1 = Yellowstone, 2 = Flux, 3 = Agri
machineid=1

numHoursSEStart=3
filterHourLength=6

#gridname=uniform_240
#gfs2seWeights=/glade/u/home/zarzycki/scratch/unigridFiles/uniform_240/maps/map_gfs0.50_TO_uniform240_patc.141127.nc
#sePreFilterIC=/glade/u/home/zarzycki/scratch/unigridFiles/uniform_240/inic/inic_uniform_240_INIC.nc
#sePostFilterIC=/glade/u/home/zarzycki/scratch/unigridFiles/uniform_240/inic/inic_uniform_240_INIC.nc

gridname=uniform_60
gfs2seWeights=/glade/p/work/zarzycki/maps/gfsmaps/map_gfs0.50_TO_uniform_60_patc.nc
sePreFilterIC=/glade/p/work/zarzycki/sewx/INIC/uniform_60_INIC.nc
sePostFilterIC=/glade/p/work/zarzycki/sewx/INIC/uniform_60_INIC_filter.nc

#gridname=newgulf_30_x4
#gfs2seWeights=/glade/p/work/zarzycki/maps/gfsmaps/map_gfs0.50_TO_newgulf_30_x4_patc.nc
#sePreFilterIC=/glade/p/work/zarzycki/sewx/INIC/newgulf_30_x4_INIC.nc
#sePostFilterIC=/glade/p/work/zarzycki/sewx/INIC/newgulf_30_x4_INIC_filter.nc

#gridname=haiyan_48_x8
#gfs2seWeights=/glade/p/work/zarzycki/maps/gfsmaps/map_gfs0.50_TO_haiyan_48_x8_patc.nc
#sePreFilterIC=/glade/p/work/zarzycki/sewx/INIC/haiyan_48_x8_INIC.nc
#sePostFilterIC=/glade/p/work/zarzycki/sewx/INIC/haiyan_48_x8_INIC_filter.nc

echo "We are using $gridname for the grid"
echo "The formal SE run will start at +$numHoursSEStart hours from actual init time"

# Set input data type
# 1 = GFS forecast ICs
# 2 = ERA-interim reanalysis
# 3 = CFSR
inputdatatype=2

# 30 -> CAM5 physics, 26 -> CAM4 physics
numlevels=30

#forecast length (in days)
numdays=3

if [ $islive -ne 0 ]    # Find most recent GFS forecast
then
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
    exit 1
  fi
else     # if not live, draw from dates.txt file
  longdate=$(head -n 1 dates.txt)
  
  yearstr=${longdate:0:4}
  monthstr=${longdate:4:2}
  daystr=${longdate:6:2}
  cyclestr=${longdate:8:2}
  
  echo $yearstr
  echo $monthstr
  echo $daystr
  echo $cyclestr
  
  #Remove top line from dates file
  tail -n +2 dates.txt > dates2.txt
  mv -v dates2.txt dates.txt
fi

## Figure out the seconds which correspond to the cycle and zero pad if neces
cyclestrsec=$(($cyclestr*3600))
echo $cyclestrsec
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

if [ $islive -ne 0 ]
then
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
else
  sstmonthstr=$monthstr
  sstdaystr=$daystr
  sstyearstr=$yearstr
  sstcyclestr=$cyclestr
fi

echo "The current time is $currtime"
echo "We are using $yearstr $monthstr $daystr $cyclestr Z ($cyclestrsec seconds) for GFS data"
echo "We are using $sstyearstr $sstmonthstr $sstdaystr $sstcyclestr Z for SST data"
echo "SE initialization will occur at $se_yearstr $se_monthstr $se_daystr $se_cyclestr Z ($se_cyclestrsec seconds)"

if [ $machineid -eq 1 ]
then
  echo "Using Yellowstone"
  sewxscriptsdir=/glade/u/home/$LOGNAME/sewx/
  backupdir=/glade/u/home/$LOGNAME/sewx/scriptbackups
  gfs_files_path=/glade/p/work/$LOGNAME/sewx/GFS
  era_files_path=/glade/p/work/$LOGNAME/getECMWFdata/
  sst_files_path=/glade/p/work/$LOGNAME/sewx/SST
  gfs_to_cam_path=/glade/u/home/$LOGNAME/sewx/gfs_to_cam
  era_to_cam_path=/glade/u/home/$LOGNAME/sewx/interim_to_cam
  atm_to_cam_path=/glade/u/home/$LOGNAME/sewx/atm_to_cam
  filter_path=/glade/u/home/$LOGNAME/sewx/filter
  path_to_se_build=/glade/u/home/$LOGNAME
  path_to_rundir=/glade/scratch/$LOGNAME/$gridname/run
  path_to_nc_files=/glade/scratch/$LOGNAME/$gridname/run
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
  era_files_path=/home/$LOGNAME/getECMWFdata
  gfs_to_cam_path=/home/$LOGNAME/sewx/gfs_to_cam
  era_to_cam_path=/home/$LOGNAME/sewx/interim_to_cam
  filter_path=/home/$LOGNAME/sewx/filter
  path_to_se_build=/home/$LOGNAME
  path_to_gz_touch=$path_to_se_build/$gridname/run
  path_to_nc_files=$path_to_se_build/$gridname/run
  outputdir=/home/$LOGNAME/${gridname}/run
else
  echo "Machine not supported"
  exit 1
fi

# Set timestamp for backing up files, etc.
timestamp=`date +%Y%m%d.%H%M`
mkdir -p $backupdir

if [ $debug -ne 1 ]
then
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
      if [ $islive -ne 0 ]
      then
        rm -f gfs.t*pgrb2f00*
        gfsFTPPath=ftp://ftp.ncep.noaa.gov/pub/data/nccf/com/gfs/prod/gfs.$yearstr$monthstr$daystr$cyclestr/
        gfsFTPFile='gfs.t'$cyclestr'z.pgrb2f00'
        echo "Attempting to download ${gfsFTPPath}${gfsFTPFile}"
      else
        rm -f gfs.t*pgrb2f00*
        gfsFTPPath=ftp://nomads.ncdc.noaa.gov/GFS/analysis_only/${yearstr}${monthstr}/${yearstr}${monthstr}${daystr}/
        gfsFTPFile='gfsanl_4_'${yearstr}${monthstr}${daystr}'_'${cyclestr}'00_000.grb2'
        ##gfsFTPPath=ftp://nomads.ncdc.noaa.gov/GFS/Grid4/${yearstr}${monthstr}/${yearstr}${monthstr}${daystr}/
        ##gfsFTPFile='gfs_4_'${yearstr}${monthstr}${daystr}'_'${cyclestr}'00_000.grb2'
        echo "Attempting to download ${gfsFTPPath}${gfsFTPFile}"
      fi
      ## Scrape for files
      error=1
      while [ $error != 0 ]
      do
        wget -nv $gfsFTPPath$gfsFTPFile
        error=`echo $?`
        if [ $error -ne 0 ]
        then
          echo "Cannot get file, will wait 2 min and scrape again"
          sleep 120
        fi
      done
      mv $gfsFTPFile 'gfs_atm_'$yearstr$monthstr$daystr$cyclestr'.grib2'
      
      ## Pull sea surface temps, need to rename and delete (if necess)
      echo "Getting SST data"
      if [ $islive -ne 0 ]
      then
        # Here is where we get the "live" GDAS SST file
        rm -f gdas1*sstgrb*
        sstFTPPath=ftp://ftp.ncep.noaa.gov/pub/data/nccf/com/gfs/prod/gdas.$sstyearstr$sstmonthstr$sstdaystr/
        sstFTPFile='gdas1.t'$sstcyclestr'z.sstgrb.grib2'
        echo "Attempting to download ${sstFTPPath}${sstFTPFile}"
      else
        # Here is where we get the archived historical SST file
        rm -f gdas1*sstgrb*
        sstFTPPath=http://nomads.ncdc.noaa.gov/data/gdas/${yearstr}${monthstr}/${yearstr}${monthstr}${daystr}/
        #sstFTPFile='gdas1.t'$sstcyclestr'z.sstgrb.grib2'
        sstFTPFile='gdas-sstgrb_3_'${yearstr}${monthstr}${daystr}'_'$sstcyclestr'00_000.grb2'
        #sstFTPFile='gdas-sstgrb_3_'${yearstr}${monthstr}${daystr}'_0600_000.grb2'
        echo "Attempting to download ${sstFTPPath}${sstFTPFile}"
      fi
      ## Scrape for files
      error=1
      while [ $error != 0 ]
      do
        wget -nv $sstFTPPath$sstFTPFile
        error=`echo $?`
        if [ $error -ne 0 ]
        then
          echo "Cannot get file, will wait 2 min and scrape again"
          sleep 120
        fi
      done
      mv $sstFTPFile 'gfs_sst_'$yearstr$monthstr$daystr$cyclestr'.grib2'
  elif [ $inputdatatype -eq 2 ]
  then  
      echo "Using ERA-Interim forecast ICs"
      echo "Cding to ERA-Interim interpolation directory"
      cd $era_files_path
      python getInterim.py ${yearstr}${monthstr}${daystr} ${cyclestr}
      ncks -A ERA-Int_ml_${yearstr}${monthstr}${daystr}${cyclestr}.nc ERA-Int_sfc_${yearstr}${monthstr}${daystr}${cyclestr}.nc
      mv -v ERA-Int_sfc_${yearstr}${monthstr}${daystr}${cyclestr}.nc ERA-Int_${yearstr}${monthstr}${daystr}${cyclestr}.nc
      rm -f ERA-Int_ml_${yearstr}${monthstr}${daystr}${cyclestr}.nc
  elif [ $inputdatatype -eq 3 ]
  then  
      echo "Using CFSR ICs"
      echo "Cding to GFS interpolation directory since they are practically the same thing"
      cd $gfs_files_path
      STCUTARR=(26 21 16 11 06 01)
      ENCUTARR=(99 25 20 15 10 05)
      zero=0
      index=$zero
      for FILEDAY in "${STCUTARR[@]}"
      do
        if [ "$daystr" -ge "$FILEDAY" ]
        then
      #    echo $FILEDAY
          break
        fi
        index=$((index+1))
      done
      #echo $index
      if [[ "$index" -eq "$zero" ]]; then
        ENCUTARR[${zero}]=`date -d "$monthstr/1 + 1 month - 1 day" "+%d"`
        echo "Last day of month ($monstr) is $ENCUTARR[${zero}]"
      fi
      ## NEED TO IMPLEMENT LEAP YEAR FIX
      CFSRFILENAME=pgbhnl.gdas.${yearstr}${monthstr}${FILEDAY}-${yearstr}${monthstr}${ENCUTARR[$index]}.tar
      echo "Getting file: ${CFSRFILENAME}"
      wget -q --load-cookies ~/.thecookies http://rda.ucar.edu/data/ds093.0/${yearstr}/${CFSRFILENAME}
      
      tar -xvf $CFSRFILENAME
      mv pgbhnl.gdas.${yearstr}${monthstr}${daystr}${cyclestr}.grb2 'cfsr_atm_'$yearstr$monthstr$daystr$cyclestr'.grib2'
      rm pgbhnl.gdas.*
  else
      echo "Incorrect model IC entered"
      exit 1
  fi

  set +e #Need to turn off error checking b/c NCL returns 0 even if fatal
  ### We can probably clean this up by merging the above sed commands into command line arguments
  ### then put this if/else statement up inside the whole get data structure above
  if [ $inputdatatype -eq 1 ] # GFS
  then
    echo "Cding to GFS interpolation directory"
    cd $atm_to_cam_path 
    echo "Doing NCL"    
    ncl -n atm_to_cam.ncl 'datasource="GFS"'     \
        numlevels=$numlevels \
        YYYYMMDDHH=${yearstr}${monthstr}${daystr}${cyclestr} \
       'gridname = "'$gridname'"' \
       'data_filename = "'$gfs_files_path'/gfs_atm_'$yearstr$monthstr$daystr$cyclestr'.grib2"'  \
       'wgt_filename="'${gfs2seWeights}'"' \
       'se_inic = "'${sePreFilterIC}'"'
       
    ncl sst_interp.ncl initdate=${yearstr}${monthstr}${daystr}${cyclestr} 'sst_file_full = "'$gfs_files_path'/gfs_sst_'$yearstr$monthstr$daystr$cyclestr'.grib2"'
  elif [ $inputdatatype -eq 2 ] # ERA
  then
    echo "CD ing to ERA-interim interpolation directory"
    cd $atm_to_cam_path

    ncl -n atm_to_cam.ncl 'datasource="ERAI"'     \
        numlevels=$numlevels \
        YYYYMMDDHH=${yearstr}${monthstr}${daystr}${cyclestr} \
       'gridname = "'$gridname'"' \
       'data_filename = "/glade/p/work/zarzycki/getECMWFdata/ERA-Int_'$yearstr$monthstr$daystr$cyclestr'.nc"'  \
       'wgt_filename="/glade/p/work/zarzycki/getECMWFdata/ERA_to_uniform_60_patch.nc"' \
       'se_inic = "'${sePreFilterIC}'"'

#    ncl -n se_interp.ncl machineid=$machineid \
#      numlevels=$numlevels \
#      YYYYMMDDHH=${yearstr}${monthstr}${daystr}${cyclestr} \
#      'gridname = "'$gridname'"' \
#      'ERA_dir="'${era_files_path}'"'

#    ncl sst_interp.ncl 'YYYYMMDDHH="'${yearstr}${monthstr}${daystr}${cyclestr}'"' \
#      'dataSource="ERAI"' \
#      'sstDataFile = "'${era_files_path}'/ERA-Int_sfc_'$yearstr$monthstr$daystr$cyclestr'.nc"' \
#      'iceDataFile = "'${era_files_path}'/ERA-Int_sfc_'$yearstr$monthstr$daystr$cyclestr'.nc"'  

  elif [ $inputdatatype -eq 3 ] # CFSR
  then
      echo "CD ing to interpolation directory"
      cd $atm_to_cam_path 
    
     ncl -n atm_to_cam.ncl 'datasource="CFSR"'     \
        numlevels=$numlevels \
        YYYYMMDDHH=${yearstr}${monthstr}${daystr}${cyclestr} \
       'gridname = "'$gridname'"' \
       'data_filename = "'$gfs_files_path'/cfsr_atm_'$yearstr$monthstr$daystr$cyclestr'.grib2"'  \
       'wgt_filename="'${gfs2seWeights}'"' \
       'se_inic = "'${sePreFilterIC}'"'
  else
      echo "Incorrect model IC entered"
      exit 1
  fi
  # Since NCL doesn't return non-zero codes, I have NCL returning a non-zero code
  # if successful! However, this means we have to check if code is successful with
  # something other than zero. Generally, if NCL fails expect a 0 return, but lets
  # be safe and call everything non-9.
  if [[ $? -ne 9 ]]
  then
    echo "NCL exited with non-9 error code"
    exit 240
  fi
  echo "NCL completed successfully"
  set -e # Turn error checking back on
  
fi #End debug if statement

#cd /home/zarzycki/sewx/INIC
#cp tcforecast_60_x4_INIC_filter.nc tcforecast_60_x4_INIC.nc
cd $path_to_se_build/$gridname

echo "Turning off archiving and restart file output in env_run.xml"
./xmlchange -v -file env_run.xml -id DOUT_S -val FALSE
./xmlchange -v -file env_run.xml -id REST_OPTION -val nyears
./xmlchange -v -file env_run.xml -id REST_N -val 9999
echo "Update env_run.xml with runtime parameters"
./xmlchange -v -file env_run.xml -id RUN_STARTDATE -val $yearstr-$monthstr-$daystr
./xmlchange -v -file env_run.xml -id START_TOD -val $cyclestrsec
./xmlchange -v -file env_run.xml -id STOP_OPTION -val nhours
./xmlchange -v -file env_run.xml -id STOP_N -val ${filterHourLength}
cp -v user_nl_cam_filter user_nl_cam

echo "Setting input land dataset"
# Copy dummy lnd namelist over with commented "!finidat" line
# If input land DOES exist, we'll sed in the file and remove ! comment
# If it does not exist, we'll do nothing and let CESM use arbitrary ICs

#cp -v user_nl_clm_dummy user_nl_clm
if [ -f /home/zarzycki/${gridname}/run/clmstart/${gridname}.clm2.r.${se_yearstr}-${se_monthstr}-${se_daystr}-${se_cyclestrsec}.nc ]
then
    echo "Land file exists, sedding that into CLM namelist"
    sed -i 's?!finidat.*?finidat='"'"/home/zarzycki/"${gridname}"/run/clmstart/"${gridname}".clm2.r."${se_yearstr}"-"${se_monthstr}"-"${se_daystr}"-"${se_cyclestrsec}".nc"'"'?' user_nl_clm
else
    echo "Land file DOES NOT EXIST, will use arbitrary CESM spinup"
    #echo "OK, Colin is cheating and using a different land file"
    #echo "He really should specify 3-4 files by month as dummies instead of CESM cold starts"
    #sed -i 's?!finidat.*?finidat='"'"/home/zarzycki/"${gridname}"/run/clmstart/"${gridname}".clm2.r.2012-08-24-10800.nc"'"'?' user_nl_clm
fi

echo "Running again!" > ${path_to_rundir}/testrunning.gz

## Get number of log files archived
numlogfiles=`ls ${path_to_rundir}/*.gz | wc -l`

echo $numlogfiles

#cp -v ${sst_files_path}/sst_1x1.nc ${path_to_rundir}/INIC
## IF GFS
#ncap2 -s "SST_cpl=SST_cpl-273.15" ${path_to_rundir}/INIC/sst_1x1.nc ${path_to_rundir}/INIC/sst_1x1.nc.new
#mv -v ${path_to_rundir}/INIC/sst_1x1.nc.new ${path_to_rundir}/INIC/sst_1x1.nc
#./xmlchange -v -file env_run.xml -id SSTICE_DATA_FILENAME -val "$DIN_LOC_ROOT/atm/cam/sst/sst_HadOIBl_bc_1x1_1850_2009_c101029.nc"
#./xmlchange -v -file env_run.xml -id SSTICE_DATA_FILENAME -val "${path_to_rundir}/INIC/sst_1x1.nc"

if [ $debug -ne 1 ]
then
  echo "Begin call to filter-run"
  if [ $machineid -eq 1 ]
  then
    bsub < $gridname.run
    echo "Skip me"
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

  while [ `ls ${path_to_rundir}/*.gz | wc -l` == $numlogfiles ]
  do
    sleep 20
    echo "Sleeping"
  done
  echo "Run over done sleeping will hold for 60 more sec to make sure files moved" 
  sleep 40


#  echo "Exiting because Colin doesn't want to run the full forecast sim"
#  exit 0

  ## Run NCL filter
  cd $filter_path
  echo "Running filter"
  filtfile_name=$gridname.cam.h0.$yearstr-$monthstr-$daystr-$cyclestrsec.nc
  ncl lowmemfilter.ncl \
   endhour=${filterHourLength} tcut=6 \
  'filtfile_name = "'${path_to_rundir}'/'${filtfile_name}'"' \
  'writefile_name = "'${sePostFilterIC}'"'
fi

exit

echo "done with filter, removing filter files"
mkdir -p $path_to_nc_files/filtered
mv -v $path_to_nc_files/*.nc $path_to_nc_files/filtered
## Delete filter files that aren't h0 since those are the only ones we care about.
find $path_to_nc_files/filtered/ -type f -not -name '*.cam.h0.*.nc' | xargs rm
## For now, I'm going to go back and delete all the h0 files in filtered
rm -v $path_to_nc_files/filtered/*.cam.h0.*.nc

echo "Make changes in CESM-SE namelist"
cd $path_to_se_build/$gridname
./xmlchange -v -file env_run.xml -id RUN_STARTDATE -val $se_yearstr-$se_monthstr-$se_daystr
./xmlchange -v -file env_run.xml -id START_TOD -val $se_cyclestrsec
./xmlchange -v -file env_run.xml -id STOP_OPTION -val ndays
./xmlchange -v -file env_run.xml -id STOP_N -val $numdays
cp -v user_nl_cam_run user_nl_cam

if [ $debug -ne 1 ]
then
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
  
  procdir=$outputdir/proc
  
  mkdir -p $procdir
  mkdir -p $procdir/images
  mkdir -p $procdir/text
  mkdir -p $procdir/nl_files
  mkdir -p $procdir/logs
  cd $outputdir
  
  echo "Running again!" > ${path_to_rundir}/testrunning.gz

  ## Get number of log files archived
  numlogfiles=`ls ${path_to_rundir}/*.gz | wc -l`

  echo $numlogfiles

  ## Hold script while log files from filter run haven't been archived yet
  while [ `ls ${path_to_rundir}/*.gz | wc -l` == $numlogfiles ]
  do
    sleep 20
    echo "Sleeping"
  done
  echo "Run over done sleeping will hold for 30 more sec to make sure files moved"
  
  sleep 30
fi

set +e

## Move other cam files to dir
mv *.cam.h*.nc $procdir
cp *_in seq_maps.rc *_modelio.nml docn.streams.txt.prescribed $procdir/nl_files
mv *.txt $procdir/text
mv *.log.* $procdir/logs
mv timing/ $procdir/

## Move land files to new restart location
cd $path_to_nc_files
landdir=$path_to_nc_files/clmstart
mkdir $landdir
echo "Removing 06Z and 18Z land restart files"
rm -v *.clm2.r.*32400.nc
rm -v *.clm2.r.*75600.nc
echo "Moving land restart files"
mv -v *.clm2.r.*nc $landdir
rm -v *.clm2.r.*.nc

echo "Deleting boring restart files"
# Delete boring restart files
cd $path_to_nc_files
rm -v *.clm2.rh0.*.nc
rm -v *.docn.rs1.*.bin
rm -v *.cam.r.*.nc
rm -v *.cam.rs.*.nc
rm -v *.cpl.r.*.nc
rm -v *.clm2.rh0.*.nc
rm -v *.cam.rh3.*.nc
rm -v *.cice.r.*.nc
rm -v rpointer.*

#cp -v /home/zarzycki/sewx/INIC/${gridname}_INIC_filter.nc $procdir
#cp -v /home/zarzycki/sewx/SST/sst_1x1.nc $procdir

mv -v $procdir $outputdir/${yearstr}${monthstr}${daystr}${cyclestr}

date

cd $sewxscriptsdir

if [ $islive -ne 0 ]
then
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
  #sed -i 's?^#SBATCH --begin.*?#SBATCH --begin='$resubhour':30?' testshell.csh
fi

bsub < notliveshell.csh

exit 0

