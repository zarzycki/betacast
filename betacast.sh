#!/bin/bash

##=======================================================================
#PBS -N wx-driver
#PBS -A P54048000 
#PBS -l walltime=6:00:00
#PBS -q share
#PBS -k oe
#PBS -m a 
#PBS -M zarzycki@ucar.edu
#PBS -l select=1:ncpus=8
##=======================================================================

###################################################################################
# Colin Zarzycki (czarzycki@psu.edu)
#
# Driver script for running CESM/CAM in "forecast" or "hindcast" mode.
# This code will either read UTC unix clock or read in a specified list of dates,
# download dataset, map to CAM grid, and run a forecast.
#
# Generally can be executed on login nodes or in the background, ex:
# $> nohup ./main-realtime.sh nl.neconus30x8 &
# but see above for example of PBS options to submit to batch nodes
#
# Details can be found in:
# C. M. Zarzycki and C. Jablonowski (2015), Experimental tropical cyclone forecasts 
# using a variable-resolution global model. Mon. Weat. Rev., 143(10), 4012–4037.
# doi:10.1175/MWR-D-15-0159.1.
###################################################################################

set -e
#set -v

SCRIPTPATH=$(dirname "$(realpath "$0")")
echo "Our script path is $SCRIPTPATH"
source ${SCRIPTPATH}/utils.sh   # Source external bash functions

# Set files
MACHINEFILE=${1}
NAMELISTFILE=${2}
OUTPUTSTREAMS=${3}
# If relative path, convert to absolute path
if [[ "$MACHINEFILE" != /* ]] && [[ "$MACHINEFILE" != ~* ]]; then MACHINEFILE=${PWD}/${MACHINEFILE}; fi
if [[ "$NAMELISTFILE" != /* ]] && [[ "$NAMELISTFILE" != ~* ]]; then NAMELISTFILE=${PWD}/${NAMELISTFILE}; fi
if [[ "$OUTPUTSTREAMS" != /* ]] && [[ "$OUTPUTSTREAMS" != ~* ]]; then OUTPUTSTREAMS=${PWD}/${OUTPUTSTREAMS}; fi
echo $MACHINEFILE; echo $NAMELISTFILE; echo $OUTPUTSTREAMS

# Read namelists
read_bash_nl "${NAMELISTFILE}"
read_bash_nl "${MACHINEFILE}"

set -u  # turn on crashes for unbound variables in bash

###################################################################################
############### OPTIONAL TO BE SET BY USER ########################################
path_to_nc_files=${path_to_rundir}              # Path where .nc files are
outputdir=${path_to_rundir}                     # Path where .nc files are being written
tmparchivecdir=${path_to_rundir}/proc/              # Path to temporarily stage final data
landdir=${path_to_rundir}/landstart/            # Path to store land restart files
###################################################################################
### THESE COME WITH THE REPO, DO NOT CHANGE #######################################
gfs_to_cam_path=${sewxscriptsdir}/gfs_to_cam
era_to_cam_path=${sewxscriptsdir}/interim_to_cam
atm_to_cam_path=${sewxscriptsdir}/atm_to_cam
sst_to_cam_path=${sewxscriptsdir}/sst_to_cam
filter_path=${sewxscriptsdir}/filter
###################################################################################

### setting variables not included in namelist for backwards compat
if [ -z ${CIMEbatchargs+x} ]; then CIMEbatchargs=""; fi
if [ -z ${do_runoff+x} ]; then do_runoff=false; fi
if [ -z ${keep_land_restarts+x} ]; then keep_land_restarts=true; fi
if [ -z ${perturb_namelist+x} ]; then perturb_namelist=""; fi
if [ -z ${predict_docn+x} ]; then predict_docn=false; fi
if [ -z ${archive_inic+x} ]; then archive_inic=false; fi
if [ -z ${compress_history_nc+x} ]; then compress_history_nc=true; fi
if [ -z ${add_vortex+x} ]; then add_vortex=false; fi
if [ -z ${vortex_namelist+x} ]; then vortex_namelist=""; fi
if [ -z ${save_nudging_files+x} ]; then save_nudging_files=false; fi
### Some defaults infrequently set
if [ -z ${doFilter+x} ]; then doFilter=false; fi
if [ -z ${filterOnly+x} ]; then filterOnly=false; fi
if [ -z ${numHoursSEStart+x} ]; then numHoursSEStart=3; fi
if [ -z ${filterHourLength+x} ]; then filterHourLength=6; fi
if [ -z ${filtTcut+x} ]; then filtTcut=6; fi
if [ -z ${FILTERWALLCLOCK+x} ]; then FILTERWALLCLOCK="00:29:00"; fi
if [ -z ${FILTERQUEUE+x} ]; then FILTERQUEUE="batch"; fi
if [ -z ${use_nsplit+x} ]; then use_nsplit="true"; fi

### Set correct E3SM/CESM split
if [ -z ${modelSystem+x} ]; then modelSystem=0; fi
if [ $modelSystem -eq 0 ]; then
  echo "Using CESM"
  atmName="cam"
  lndName="clm"
  lndSpecialName="clm2"
  rofName="mosart"
  rofSpecialName="mosart"
elif [ $modelSystem -eq 1 ]; then
  echo "Using E3SM"
  atmName="eam"
  lndName="elm"
  lndSpecialName="elm"
  rofName="mosart"
  rofSpecialName="mosart"
else
  echo "Unknown modeling system set for modelSystem: $modelSystem"
  exit 1
fi

# Figure out where to archive
if [ -z ${ARCHIVEDIR+x} ] || [[ -z "${ARCHIVEDIR// }" ]] ; then
  ARCHIVEDIR=${outputdir}/
else
  ARCHIVEDIR=${ARCHIVEDIR}/${casename}/
  mkdir -v -p ${ARCHIVEDIR}
fi
echo "Files will be archived in ${ARCHIVEDIR}/YYYYMMDDHH/"

### ERROR CHECKING BLOCK! #########################################################

# Exit if add_perturbs is turned on but no namelist is passed in with perturbation config settings
if ${add_perturbs} && { [ -z "$perturb_namelist" ] || [ ! -f "$perturb_namelist" ]; } ; then
  echo "add_perturbs is true but can't find namelist: "$perturb_namelist
  exit 1
fi

# Exit if add_perturbs is turned on but no namelist is passed in with perturbation config settings
if ${add_vortex} && { [ -z "$vortex_namelist" ] || [ ! -f "$vortex_namelist" ]; } ; then
  echo "add_vortex is true but can't find namelist: "$vortex_namelist
  exit 1
fi

# Check for deprecated namelist options
if [ -z ${anl2mdlWeights+x} ] && [ ${gfs2seWeights+x} ] ; then
  echo "WARNING: Setting anl2mdlWeights to ${gfs2seWeights}"
  echo "WARNING: This is deprecated and will be removed in the future! To fix, change 'gfs2seWeights' to 'anl2mdlWeights' in ${NAMELISTFILE}"
  anl2mdlWeights=$gfs2seWeights
fi

# Check if ncks exists
if ! type ncl &> /dev/null ; then
  echo "ERROR: ncl does not exist. Make sure ncl is in your path when betacast is invoked"
  exit 1
fi

# Check if ncks exists for compression
if ! type ncks &> /dev/null ; then
  #echo "ERROR: ncks does not exist. Make sure ncks is in your path when betacast is invoked"
  #exit 1
  echo "WARNING: ncks does not exist, cannot compress. Setting compress_history_nc to 0 (false)"
  echo "WARNING: if you'd like remedy, make sure ncks is in your path when betacast.sh is invoked"
  compress_history_nc=0
fi

# Adjust bools (for backwards compatibility, 0 = false and 1 = true)
check_bool "islive" $islive
check_bool "debug" $debug
check_bool "do_runoff" $do_runoff
check_bool "keep_land_restarts" $keep_land_restarts
check_bool "predict_docn" $predict_docn
check_bool "archive_inic" $archive_inic
check_bool "compress_history_nc" $compress_history_nc

###################################################################################

# do some stability calcs
# if USERSTAB is 0, use internal calcs.
# if USERSTAB is <0, use se_nsplit=-1
# If USERSTAB >0, use the value in seconds for dt_dyn and calculate nsplit accordingly from DTIME
USERSTABTF=`python -c "print('TRUE' if ${USERSTAB} > 0 else 'FALSE')"`
if [ ${USERSTABTF} == 'FALSE' ] ; then
  if [ `python -c "print('TRUE' if ${USERSTAB} < -0.001 else 'FALSE')"` == 'FALSE' ]; then
    STABILITY=`python -c "print(30./${FINERES}*450.)"`
    VALIDSTABVAL=true
    echo "Dynamic stability for ne${FINERES} to be ${STABILITY} seconds"
  else
    STABILITY=-1
    VALIDSTABVAL=false   # Setting to false means we won't calculate nsplit for SE/HOMME later
    echo "Dynamic stability for ne${FINERES} to be ${STABILITY} and internal STABILITY setting"
  fi
else
  STABILITY=${USERSTAB}
  VALIDSTABVAL=true
  echo "User defined stability set to ${STABILITY}"
fi
echo "VALIDSTABVAL is set to "${VALIDSTABVAL}

## Create paths to generate initial files if they don't exist...
mkdir -p ${pathToINICfiles}
mkdir -p ${pathToSSTfiles}

# Set timestamp for backing up files, etc.
timestamp=`date +%Y%m%d.%H%M`
uniqtime=`date +"%s%N"`

echo "We are using ${casename} for the case"

if [ $islive = true ] ; then    # Find most recent GFS forecast
  ## Here we get two digit strings for UTC time for month, day, year
  ## We also get current time in hoursminutes (because the GFS output lags by 3.5 hours)
  monthstr=`date -u +%m`
  daystr=`date -u +%d`
  yearstr=`date -u +%Y`
  currtime=`date -u +%H%M`
  machzone=`date +%z`
  twodaysago=`date --date='2 days ago' -u +"%Y%m%d"`
  
  ## Use currtime to figure out what is the latest cycle we have access to
  if [ $currtime -lt 0328 ] ; then
    echo "12Z cycle"
    monthstr=`date --date="yesterday" -u +%m`
    daystr=`date --date="yesterday" -u +%d`
    yearstr=`date --date="yesterday" -u +%Y`
    twodaysago=`date --date='3 days ago' -u +"%Y%m%d"`
    cyclestr=12
  elif [ $currtime -lt 0928 ] ; then
    echo "00Z cycle"
    cyclestr=00
  elif [ $currtime -lt 1528 ] ; then
    echo "00Z cycle"
    cyclestr=00
  elif [ $currtime -lt 2128 ] ; then
    echo "12Z cycle"
    cyclestr=12
  elif [ $currtime -ge 2128 ] ; then
    echo "12Z cycle"
    cyclestr=12
  else
    echo "Can't figure out start time"
    exit 1
  fi

else     # if not live, draw from head of dates.txt file

  # Figure out if dates.CASE.txt exists and where it is located
  datesbase=dates.${casename}.txt
  if [ -f "./dates/${datesbase}" ]; then
    # Preferred location is dates subdir as of 4/13/2022
    echo "$datesbase exists in dates subdirectory"
    datesfile=./dates/${datesbase}
  else
    if [ -f "./${datesbase}" ]; then
      echo "WARNING! $datesbase isn't in dates subdir but exists in home dir!"
      echo "This is allowed for backwards compatibility but is less organized!"
      echo "You should create a subdir called 'dates' and put the dates.CASE.txt file there!"
      datesfile=${datesbase}
    else
      if [[ ${datestemplate+x} && -f ./dates/${datestemplate} ]]; then
        echo "Didn't find case dates, but you specified a datestemplate in the namelist..."
        echo "So I'm copying $datestemplate to this case and using that..."
        cp -v ./dates/${datestemplate} ./dates/${datesbase}
        datesfile=./dates/${datesbase}
      else
        echo "Uh, oh. Can't find a dates file OR template AND run isn't live. Exiting..."
        exit 1
      fi
    fi
  fi

  echo "Using dates in: "${datesfile}
  longdate=$(head -n 1 ${datesfile})
  
  # Do some simple error trapping on date string to ensure validity
  if [ -z "$longdate" ]; then { echo "Date string passed in is empty, exiting..." ; exit 91; } ; fi
  if [ ${#longdate} -ne 10 ]; then { echo "Malformed date string, $longdate is ${#longdate} characters, needs 10 (YYYYMMDDHH). Exiting..." ; exit 92; } ; fi
  if [[ -n $(echo $longdate | tr -d '[0-9]') ]]; then { echo "Malformed date string, $longdate contains non-numeric values. Exiting..." ; exit 93; } ; fi
  
  echo "Getting parsed time from $longdate"
  parse_YYYYMMDDHH $longdate    
  echo "From datesfile, read in: "$yearstr' '$monthstr' '$daystr' '$cyclestr'Z'
  # Do some error trapping on returned time values
  if (( yearstr > 3000 || yearstr < 1 )); then { echo "Year set to $yearstr, this sounds wrong, exiting..." ; exit 94; } ; fi
  if (( cyclestr > 23 )); then { echo "Cycle string set to $cyclestr Z, this sounds wrong, exiting..." ; exit 95; } ; fi
fi

## Figure out the seconds which correspond to the cycle and zero pad if neces
get_cyclestrsec "$cyclestr"

## Figure out what the SE start time will be after filter
if [ $numHoursSEStart -lt 6 ] ; then
  let se_cyclestr=$cyclestr+03
  while [ ${#se_cyclestr} -lt 2 ];
  do
    se_cyclestr="0"$se_cyclestr
  done
  se_monthstr=$monthstr
  se_daystr=$daystr
  se_yearstr=$yearstr
  let se_cyclestrsec=$((10#$se_cyclestr))*3600
  while [ ${#se_cyclestrsec} -lt 5 ];
  do
    se_cyclestrsec="0"$se_cyclestrsec
  done
else
  echo "SE forecast lead time too long, 18Z cycle causes trouble"
  echo "Not supported."
  exit 1
fi

if [ $islive = true ] ; then
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

yestmonthstr=`date --date="yesterday" -u +%m`
yestdaystr=`date --date="yesterday" -u +%d`
yestyearstr=`date --date="yesterday" -u +%Y`

echo "We are using $yearstr $monthstr $daystr $cyclestr Z ($cyclestrsec seconds) for ATM init. data"
echo "We are using $sstyearstr $sstmonthstr $sstdaystr $sstcyclestr Z for SST init. data"
if $doFilter ; then
  echo "Filter: True model init. will occur at $se_yearstr $se_monthstr $se_daystr $se_cyclestr Z ($se_cyclestrsec seconds)"
else
  echo "No filter: True model init. will occur at $se_yearstr $se_monthstr $se_daystr $cyclestr Z ($cyclestrsec seconds)"
fi

if $runmodel ; then

############################### GET DYCORE INFO ############################### 

cd $path_to_case
DYCORE=`./xmlquery CAM_DYCORE | sed 's/^[^\:]\+\://' | xargs`
echo "DYCORE: "$DYCORE

if [ $debug = false ] ; then
############################### GET ATM DATA ############################### 

  if [ $atmDataType -eq 1 ] ; then
    echo "Using GFS forecast ICs"
    echo "Getting GFS conditions"
    mkdir -p $gfs_files_path
    cd $gfs_files_path
    LOCALGFSFILE='gfs_atm_'$yearstr$monthstr$daystr$cyclestr'.grib2'
    ## Pull atmospheric conditions
    ## Need to write an if call depending on how far back the GFS files are located
    ## The NCEP server only holds them for a few days -- otherwise go to nomads at NCDC
    # gfsFtpPath=http://nomads.ncdc.noaa.gov/data/gfsanl/$yearstr$monstr/$yearstr$monstr$daystr/
    # gfs.t12z.pgrb2f00.grb2
    if [ ! -f ${LOCALGFSFILE} ]; then
      echo "Getting Atmo file"
      if [ $islive = true ] ; then
        rm -f gfs.t*pgrb2f00*
        gfsFTPPath=ftp://ftp.ncep.noaa.gov/pub/data/nccf/com/gfs/prod/gfs.$yearstr$monthstr$daystr/$cyclestr/atmos/
        #gfsFTPFile='gfs.t'$cyclestr'z.pgrb2f00'
        gfsFTPFile='gfs.t'$cyclestr'z.pgrb2.0p25.anl'
        #gfsFTPFile='gfs.t'$cyclestr'z.pgrb2.0p50.anl'
        echo "Attempting to download ${gfsFTPPath}${gfsFTPFile}"
        ## Scrape for files
        error=1
        while [ $error != 0 ] ; do
          wget -nv --read-timeout=30 $gfsFTPPath$gfsFTPFile
          error=`echo $?`
          if [ $error -ne 0 ] ; then
            echo "Cannot get file, will wait 2 min and scrape again"
            sleep 120
          fi
        done
      else                  # Copy GFS data from RDA archive at NCAR
        rm -f gfs.t*pgrb2f00*
        gfsFTPPath=/glade/collections/rda/data/ds084.1/${yearstr}/${yearstr}${monthstr}${daystr}/
        gfsFTPFile=gfs.0p25.${yearstr}${monthstr}${daystr}${cyclestr}.f000.grib2
        cp ${gfsFTPPath}/${gfsFTPFile} .
        echo "Attempting to copy ${gfsFTPPath}${gfsFTPFile}"
      fi
      mv $gfsFTPFile ${LOCALGFSFILE}
    fi
  elif [ $atmDataType -eq 2 ] ; then  
    echo "Using ERA-Interim forecast ICs"
    echo "Cding to ERA-Interim interpolation directory"
    mkdir -p $era_files_path
    cd $era_files_path
    LOCALGFSFILE=ERA-Int_${yearstr}${monthstr}${daystr}${cyclestr}.nc
    if [ ! -f ${LOCALGFSFILE} ]; then
      echo "support broken for auto download ERA, please prestage"
      exit
      python getInterim.py ${yearstr}${monthstr}${daystr} ${cyclestr}
      ncks -A ERA-Int_ml_${yearstr}${monthstr}${daystr}${cyclestr}.nc ERA-Int_sfc_${yearstr}${monthstr}${daystr}${cyclestr}.nc
      mv -v ERA-Int_sfc_${yearstr}${monthstr}${daystr}${cyclestr}.nc ERA-Int_${yearstr}${monthstr}${daystr}${cyclestr}.nc
      rm -f ERA-Int_ml_${yearstr}${monthstr}${daystr}${cyclestr}.nc
    fi
  elif [ $atmDataType -eq 3 ] ; then  
    echo "Using CFSR ICs"
    echo "Cding to GFS interpolation directory since they are practically the same thing"
    mkdir -p $gfs_files_path
    cd $gfs_files_path

    LOCALCFSRFILE='cfsr_atm_'$yearstr$monthstr$daystr$cyclestr'.grib2'
    if [ ! -f ${LOCALCFSRFILE} ]; then
      STCUTARR=(26 21 16 11 06 01)
      ENCUTARR=(99 25 20 15 10 05)
      zero=0
      index=$zero
      for FILEDAY in "${STCUTARR[@]}" ; do
        if [ "$daystr" -ge "$FILEDAY" ] ; then
      #    echo $FILEDAY
          break
        fi
        index=$((index+1))
      done
      #echo $index
      if [[ "$index" -eq "$zero" ]] ; then
        ENCUTARR[${zero}]=`date -d "$monthstr/1 + 1 month - 1 day" "+%d"`
        echo "Last day of month ($monthstr) is $ENCUTARR[${zero}]"
      fi
      ## NEED TO IMPLEMENT LEAP YEAR FIX

      echo "Getting file: ${CFSRFILENAME}"
      #Register with RDA, then do following command to get wget cookies
      #wget --save-cookies ~/.thecookies --post-data="email=your_email_address&passwd=your_password&action=login" https://rda.ucar.edu/cgi-bin/login
      CFSRFILENAME=pgbhnl.gdas.${yearstr}${monthstr}${FILEDAY}-${yearstr}${monthstr}${ENCUTARR[$index]}.tar
      if [[ $(hostname -s) = cheyenne* ]]; then
        cp /glade/collections/rda/data/ds093.0/${yearstr}/${CFSRFILENAME} .
      else
        wget -nv --load-cookies ~/.thecookies http://rda.ucar.edu/data/ds093.0/${yearstr}/${CFSRFILENAME}
      fi
      tar -xvf $CFSRFILENAME
      mv pgbhnl.gdas.${yearstr}${monthstr}${daystr}${cyclestr}.grb2 ${LOCALCFSRFILE}
      rm pgbhnl.gdas.*
    fi
  elif [ $atmDataType -eq 4 ] ; then  
    echo "Using ERA5 forecast ICs"
    echo "Cding to ERA5 interpolation directory"
    mkdir -p $era_files_path
    cd $era_files_path
    LOCALGFSFILE=ERA5_${yearstr}${monthstr}${daystr}${cyclestr}.nc
    ERA5RDA=0  # Set whether or not ERA5 is local (0 = local, 1 = RDA)
    if [ ! -f ${LOCALGFSFILE} ]; then
      echo "cannot find: ${era_files_path}/${LOCALGFSFILE}"
      if [[ "$MACHINEFILE" == *cheyenne* ]]; then
        echo "We are on Cheyenne, so even though we lack a local file, we can use RDA"
        ERA5RDA=1
        RDADIR=/glade/collections/rda/data/ds633.0/
      elif [[ "$MACHINEFILE" == *cori* ]]; then
        echo "We are on Cori, so even though we lack a local file, we can use RDA"
        ERA5RDA=1
        RDADIR=/global/cfs/projectdirs/m3522/cmip6/ERA5/    
      else
        echo "support broken for auto download ERA, please prestage!"
        exit
      fi
    fi
  else
    echo "Incorrect model IC entered"
    exit 1
  fi

############################### GET SST / NCL ############################### 

  if [ ${sstDataType} -eq 1 ] ; then
    SSTTYPE=GDAS
    ## Pull sea surface temps, need to rename and delete (if necess)
    #echo "Live GFS SST hasn't been updated in a while..." ; exit
    echo "Getting SST data"
    mkdir -p ${sst_files_path}
    cd ${sst_files_path}
    if [ $islive = true ] ; then
      # Here is where we get the "live" GDAS SST file
      rm -f gdas1*sstgrb*
      #sstFTPPath=ftp://ftp.ncep.noaa.gov/pub/data/nccf/com/gfs/prod/sst.${yestyearstr}${yestmonthstr}${yestdaystr}/
      sstFTPPath=ftp://ftp.ncep.noaa.gov//pub/data/nccf/com/nsst/v1.2/nsst.${yestyearstr}${yestmonthstr}${yestdaystr}/
      sstFTPFile='rtgssthr_grb_0.5.grib2'
      echo "Attempting to download ${sstFTPPath}${sstFTPFile}"
    else
      echo "NCEP broke support for historical GDAS, use NOAAOI instead."
      exit
    fi
    ## Scrape for files
    error=1
    while [ $error != 0 ]
    do
      wget -nv --read-timeout=30 -nv $sstFTPPath$sstFTPFile
      error=`echo $?`
      if [ $error -ne 0 ]
      then
        echo "Cannot get file, will wait 2 min and scrape again"
        sleep 120
      fi
    done
    sstFile='gfs_sst_'$yearstr$monthstr$daystr$cyclestr'.grib2'
    mv ${sstFTPFile} ${sstFile}
    iceFile=''   # do not need icefile since ice stored on sstfile
  elif [ ${sstDataType} -eq 2 ] ; then  
    echo "ERA SST not quite supported yet..." ; exit 1
  elif [ ${sstDataType} -eq 3 ] ; then
    SSTTYPE=NOAAOI
    echo "Using NOAAOI SSTs"
    mkdir -p ${sst_files_path}
    cd ${sst_files_path}
    sstFile=sst.day.mean.${yearstr}.v2.nc
    if [ ! -f ${sst_files_path}/${sstFile} ] ; then
      echo "NOAAOI file doesn't exist, need to download"
      sstFTPPath=ftp://ftp.cdc.noaa.gov/Datasets/noaa.oisst.v2.highres/
      error=1
      while [ $error != 0 ]
      do
        wget -nv ${sstFTPPath}/${sstFile}
        error=`echo $?`
        if [ $error -ne 0 ]
        then
          echo "Cannot get file, will wait 2 min and scrape again"
          sleep 120
        fi
      done        
    fi
    iceFile=icec.day.mean.${yearstr}.v2.nc
    if [ ! -f ${sst_files_path}/${iceFile} ] ; then
      echo "NOAAOI file doesn't exist, need to download"
      sstFTPPath=ftp://ftp.cdc.noaa.gov/Datasets/noaa.oisst.v2.highres/
      error=1
      while [ $error != 0 ] ; do
        wget -nv ${sstFTPPath}/${iceFile}
        error=`echo $?`
        if [ $error -ne 0 ]
        then
          echo "Cannot get file, will wait 2 min and scrape again"
          sleep 120
        fi
      done        
    fi
  else
      echo "Incorrect SST data type entered" ; exit 1
  fi

  # Switch bash bool to int for NCL input
  if [ $predict_docn = true ]; then INT_PREDICT_DOCN=1; else INT_PREDICT_DOCN=0; fi

  set +e
  cd ${sst_to_cam_path}
  (set -x; ncl sst_interp.ncl 'initdate="'${yearstr}${monthstr}${daystr}${cyclestr}'"' \
    predict_docn=${INT_PREDICT_DOCN} \
    'datasource="'${SSTTYPE}'"' \
    'sstDataFile = "'${sst_files_path}/${sstFile}'"' \
    'iceDataFile = "'${sst_files_path}/${iceFile}'"' \
    'SST_write_file = "'${sstFileIC}'"' )
  if [[ $? -ne 9 ]] ; then
    echo "NCL exited with non-9 error code"
    exit 240
  fi
  echo "SST NCL completed successfully"
  set -e # Turn error checking back on
  # Removed 4/18/22 since this can be done inside sst_interp.ncl
  #ncatted -O -a units,time,o,c,"days since 0001-01-01 00:00:00" ${sstFileIC} ${sstFileIC}

  ############################### ATM NCL ############################### 

  set +e #Need to turn off error checking b/c NCL returns 0 even if fatal
  ### We can probably clean this up by merging the above sed commands into command line arguments
  ### then put this if/else statement up inside the whole get data structure above
  if [ $atmDataType -eq 1 ] ; then
    echo "Cding to GFS interpolation directory"
    cd $atm_to_cam_path 
    echo "Doing NCL"
    (set -x; ncl -n atm_to_cam.ncl 'datasource="GFS"'     \
        numlevels=${numLevels} \
        YYYYMMDDHH=${yearstr}${monthstr}${daystr}${cyclestr} \
         'dycore="'${DYCORE}'"' \
       'data_filename = "'$gfs_files_path'/gfs_atm_'$yearstr$monthstr$daystr$cyclestr'.grib2"'  \
       'wgt_filename="'${anl2mdlWeights}'"' \
       'model_topo_file="'${adjust_topo-}'"' \
       'adjust_config="'${adjust_flags-}'"' \
       'se_inic = "'${sePreFilterIC}'"' )
  elif [ $atmDataType -eq 2 ] ; then
    echo "CD ing to ERA-interim interpolation directory"
    cd $atm_to_cam_path
    (set -x; ncl -n atm_to_cam.ncl 'datasource="ERAI"'     \
        numlevels=${numLevels} \
        YYYYMMDDHH=${yearstr}${monthstr}${daystr}${cyclestr} \
         'dycore="'${DYCORE}'"' \
       'data_filename = "/glade/p/work/zarzycki/getECMWFdata/ERA-Int_'$yearstr$monthstr$daystr$cyclestr'.nc"'  \
       'wgt_filename="/glade/p/work/zarzycki/getECMWFdata/ERA_to_uniform_60_patch.nc"' \
       'se_inic = "'${sePreFilterIC}'"' )
  elif [ $atmDataType -eq 3 ] ; then
    echo "CD ing to interpolation directory"
    cd $atm_to_cam_path 
    (set -x; ncl -n atm_to_cam.ncl 'datasource="CFSR"'     \
      numlevels=${numLevels} \
      YYYYMMDDHH=${yearstr}${monthstr}${daystr}${cyclestr} \
     'dycore="'${DYCORE}'"' \
     'data_filename = "'$gfs_files_path'/cfsr_atm_'$yearstr$monthstr$daystr$cyclestr'.grib2"'  \
     'wgt_filename="'${anl2mdlWeights}'"' \
     'model_topo_file="'${adjust_topo-}'"' \
     'adjust_config="'${adjust_flags-}'"' \
     'se_inic = "'${sePreFilterIC}'"' )
  elif [ $atmDataType -eq 4 ] ; then
    echo "CD ing to ERA5 interpolation directory"
    cd $atm_to_cam_path
    if [ $ERA5RDA -eq 1 ] ; then
      (set -x; ncl -n atm_to_cam.ncl 'datasource="ERA5RDA"'     \
          numlevels=${numLevels} \
          YYYYMMDDHH=${yearstr}${monthstr}${daystr}${cyclestr} \
         'dycore="'${DYCORE}'"' \
         'data_filename = "'${RDADIR}'/e5.oper.invariant/197901/e5.oper.invariant.128_129_z.ll025sc.1979010100_1979010100.nc"'  \
         'RDADIR="'${RDADIR}'"' \
         'wgt_filename="'${anl2mdlWeights}'"' \
         'model_topo_file="'${adjust_topo-}'"' \
         'adjust_config="'${adjust_flags-}'"' \
         'se_inic = "'${sePreFilterIC}'"' )
    else
      (set -x; ncl -n atm_to_cam.ncl 'datasource="ERA5"'     \
          numlevels=${numLevels} \
          YYYYMMDDHH=${yearstr}${monthstr}${daystr}${cyclestr} \
         'dycore="'${DYCORE}'"' \
         'data_filename = "'$era_files_path'/ERA5_'$yearstr$monthstr$daystr$cyclestr'.nc"'  \
         'wgt_filename="'${anl2mdlWeights}'"' \
         'model_topo_file="'${adjust_topo-}'"' \
         'adjust_config="'${adjust_flags-}'"' \
         'se_inic = "'${sePreFilterIC}'"' )
    fi
  else
    echo "Incorrect model IC entered"
    exit 1
  fi
  # Since NCL doesn't return non-zero codes, I have NCL returning a non-zero code
  # if successful! However, this means we have to check if code is successful with
  # something other than zero. Generally, if NCL fails expect a 0 return, but lets
  # be safe and call everything non-9.
  if [[ $? -ne 9 ]] ; then echo "NCL exited with non-9 error code" ; exit 240 ; fi
  echo "ATM NCL completed successfully"
  set -e # Turn error checking back on
  
fi #End debug if statement

############################### #### ############################### 
##### ADD OR REMOVE VORTEX

if ${add_vortex} ; then
  cd $atm_to_cam_path/tcseed 
  set +e
  echo "Adding or removing a TC from initial condition based on ${vortex_namelist}"
  
  (set -x; ncl -n find-tc-fill-params.ncl 'inic_file= "'${sePreFilterIC}'"' 'pthi = "'${vortex_namelist}'"' )
  if [[ $? -ne 9 ]] ; then echo "NCL exited with non-9 error code" ; exit 240 ; fi
  
  (set -x; ncl -n seed-tc-in-ncdata.ncl   'seedfile = "'${sePreFilterIC}'"' 'pthi = "'${vortex_namelist}'"' )
  if [[ $? -ne 9 ]] ; then echo "NCL exited with non-9 error code" ; exit 240 ; fi
  
  set -e
fi

############################### #### ############################### 
##### ADD WHITE NOISE PERTURBATIONS

if ${add_noise} ; then
  set +e
  echo "Adding white noise to initial condition"
  cd $atm_to_cam_path 
  (set -x; ncl -n perturb_white_noise.ncl 'basFileName = "'${sePreFilterIC}'"' )
  set -e
fi

############################### #### ############################### 
##### ADD PERTURBATIONS

if ${add_perturbs} ; then
  echo "Adding perturbations"

#   cp /glade/scratch/zarzycki/apply-haiyan-perturb/sst_1x1_Nat-Hist-CMIP5-est1-v1-0.nc ${sstFileIC}
# 
#   cd $atm_to_cam_path
# 
#   sePreFilterIC_WPERT=${sePreFilterIC}_PERT.nc
# 
#   set +e
#   ncl -n add_perturbations_to_cam.ncl 'BEFOREPERTFILE="'${sePreFilterIC}'"'  \
#     'AFTERPERTFILE = "'${sePreFilterIC_WPERT}'"'
#   if [[ $? -ne 9 ]]
#   then
#     echo "NCL exited with non-9 error code"
#     exit 240
#   fi
#   echo "ATM NCL completed successfully"
#   set -e # Turn error checking back on
# 
#   mv ${sePreFilterIC_WPERT} ${sePreFilterIC}
  
  cd $atm_to_cam_path/perturb
  set +e
  
  ## Add perturbations to SST file
  sstFileIC_WPERT=${sstFileIC}_PERT.nc
  (set -x; ncl -n add_perturbations_to_sst.ncl 'BEFOREPERTFILE="'${sstFileIC}'"' \
     'AFTERPERTFILE = "'${sstFileIC_WPERT}'"' \
     'pthi="'${perturb_namelist}'"' )

  if [[ $? -ne 9 ]]
  then
    echo "NCL exited with non-9 error code"
    exit 240
  fi
  echo "SST perturbations added successfully"

  ## Add perturbations to ATM file
  sePreFilterIC_WPERT=${sePreFilterIC}_PERT.nc
  (set -x; ncl -n add_perturbations_to_cam.ncl 'BEFOREPERTFILE="'${sePreFilterIC}'"'  \
     'AFTERPERTFILE = "'${sePreFilterIC_WPERT}'"' \
     'pthi="'${perturb_namelist}'"' )
     
  if [[ $? -ne 9 ]]
  then
   echo "NCL exited with non-9 error code"
   exit 240
  fi
  echo "ATM NCL completed successfully"

  set -e # Turn error checking back on
 
  # move temp perturb files to overwrite unperturbed files 
  mv ${sstFileIC_WPERT} ${sstFileIC}
  mv ${sePreFilterIC_WPERT} ${sePreFilterIC}
fi

cd $path_to_case
############################### SETUP AND QUERY ############################### 

############################### CISM SETUP ############################### 

if [ -f user_nl_cism ]; then
  sed -i '/dt_count/d' user_nl_cism
  echo "dt_count = 8" >> user_nl_cism
fi

############################### CLM SETUP ############################### 

## TEMPORARY CLM -> LAND FIX
## Check if clmstart exists but landstart doesn't, move if that is the case
if [[ ! -d ${landdir} ]] && [[ -d ${path_to_rundir}/clmstart/ ]] ; then
  echo "Moving ${path_to_rundir}/clmstart/ to ${landdir}"
  mv -v ${path_to_rundir}/clmstart/ ${landdir}
else
  echo "No need to modify land directory, either appropriately exists or will be created later"
fi
## END TEMP FIX

echo "Setting input land dataset"
# Clean up file to delete any special interp lines that may be needed later (but aren't needed for native init)
#sed -i '/finidat/d' user_nl_${lndName}
sed -i '/init_interp_fill_missing_with_natveg/d' user_nl_${lndName}
sed -i '/use_init_interp/d' user_nl_${lndName}
#echo "finidat=''" >> user_nl_${lndName}

# We want to check ${landdir} for land restart files. If so, use those.
landrestartfile=${landdir}/${casename}.${lndSpecialName}.r.${yearstr}-${monthstr}-${daystr}-${cyclestrsec}.nc
echo $landrestartfile

# Check to see if file exists on native SE land grid
if [ -f ${landrestartfile} ]; then
   echo "File exists at exact time"
else
   echo "File does not exist at exact time"
   landrestartfile=${landdir}/${casename}.${lndSpecialName}.r.${yearstr}-${monthstr}-${daystr}-00000.nc
   if [ -f ${landrestartfile} ]; then
     echo "File exists at 00Z"
   else
     echo "No restart file exists, setting to empty string."
     landrestartfile=
   fi
fi
echo "landrestartfile: ${landrestartfile}"

## Now modify user_nl_${lndName}
if [ ${landrestartfile} ] ; then
  sed -i '/.*finidat/d' user_nl_${lndName}
  echo "finidat='${landrestartfile}'" >> user_nl_${lndName}
  #sed -i 's?.*finidat.*?finidat='"'${landrestartfile}'"'?' user_nl_${lndName}
else
  rawlandrestartfile=`ls ${landrawdir}/*.${lndSpecialName}.r.${yearstr}-${monthstr}-${daystr}-${cyclestrsec}.nc || true`   # check for file, suppress failed ls error if true
  echo "rawlandrestartfile: ${rawlandrestartfile}"
  if [ ! -z ${rawlandrestartfile} ]; then   # if rawlandrestartfile string is NOT empty, add it.
    sed -i '/.*finidat/d' user_nl_${lndName}
    echo "finidat='${rawlandrestartfile}'" >> user_nl_${lndName}
    echo "use_init_interp = .true." >> user_nl_${lndName}
    echo "init_interp_fill_missing_with_natveg = .true." >> user_nl_${lndName}
  else
    if [ -f user_nl_${lndName}_presave ] ; then
      echo "Using pre-written user_nl_${lndName} file"
      cp -v user_nl_${lndName}_presave user_nl_${lndName}
    else
      echo "WARNING: Land file DOES NOT EXIST, will use arbitrary user_nl_${lndName} already in folder"
      echo "!!!!!!!!!!!!!"
      sed -i '/.*finidat/d' user_nl_${lndName}
      ./xmlchange CLM_FORCE_COLDSTART=on
      #exit
      #sed -i 's?.*finidat.*?!finidat='"''"'?' user_nl_${lndName}
    fi
  fi
fi

## Append a check error to ignore inconsistencies in the dataset
if [ $modelSystem -eq 0 ]; then   # CLM/CTSM
  sed -i '/check_finidat_pct_consistency/d' user_nl_${lndName}
  sed -i '/check_finidat_year_consistency/d' user_nl_${lndName}
  echo "check_finidat_pct_consistency = .false." >> user_nl_${lndName}
  echo "check_finidat_year_consistency = .false." >> user_nl_${lndName}
elif [ $modelSystem -eq 1 ]; then   # ELM
  sed -i '/check_finidat_fsurdat_consistency/d' user_nl_${lndName}
  echo "check_finidat_fsurdat_consistency = .false." >> user_nl_${lndName}
else
  echo "Unknown modeling system set for modelSystem: $modelSystem"
  exit 1
fi

############################### ROF SETUP ############################### 

if [ $do_runoff = true ]; then

  echo "Setting input rof dataset"

  # Delete any existing input data
  sed -i '/.*finidat_rtm/d' user_nl_${rofName}
  
  # We want to check ${landdir} for land restart files. If so, use those.
  rofrestartfile=${landdir}/${casename}.${rofSpecialName}.r.${yearstr}-${monthstr}-${daystr}-${cyclestrsec}.nc
  echo $rofrestartfile

  # Check to see if file exists on native SE land grid
  if [ -f ${rofrestartfile} ]; then
     echo "File exists at exact time"
  else
     echo "File does not exist at exact time"
     rofrestartfile=${landdir}/${casename}.${rofSpecialName}.r.${yearstr}-${monthstr}-${daystr}-00000.nc
     if [ -f ${rofrestartfile} ]; then
       echo "File exists at 00Z"
     else
       echo "No restart file exists, setting to empty string."
       rofrestartfile=
     fi
  fi
  echo "rofrestartfile: ${rofrestartfile}"

  ## Now modify user_nl_${lndName}
  if [ ${rofrestartfile} ] ; then
    echo "finidat_rtm='${rofrestartfile}'" >> user_nl_${rofName}
  fi

fi

############################### GENERIC CAM SETUP ############################### 

./xmlchange PROJECT=${PROJECTID}

echo "Turning off archiving and restart file output in env_run.xml"
./xmlchange DOUT_S=FALSE,REST_OPTION=nyears,REST_N=9999
echo "Setting SST from default to our SST"
./xmlchange SSTICE_DATA_FILENAME="${sstFileIC}"
echo "Setting GLC coupling to handle forecasts across calendar years"
./xmlchange GLC_AVG_PERIOD="glc_coupling_period"
echo "Standardizing streams for SST"
./xmlchange SSTICE_YEAR_START=1,SSTICE_YEAR_END=1
echo "Setting projectID and queue"
./xmlchange --force JOB_QUEUE=${RUNQUEUE},PROJECT=${PROJECTID}
####### 
if [ "$land_spinup" = true ] ; then
  ./xmlchange REST_OPTION=ndays,REST_N=1
fi
#######

cp user_nl_${atmName} user_nl_${atmName}.BAK

echo "Update env_run.xml with runtime parameters"
./xmlchange RUN_STARTDATE=$yearstr-$monthstr-$daystr,START_TOD=$cyclestrsec,STOP_OPTION=ndays,STOP_N=$numdays

SEINIC=${sePreFilterIC}

############################### (IF) FILTER SETUP ############################### 

if $doFilter ; then
  # If filtering, need to change these options

  ./xmlchange STOP_OPTION=nhours,STOP_N=${filterHourLength}

  sed -i '/.*ncdata/d' user_nl_${atmName}
  echo "ncdata='${SEINIC}'" >> user_nl_${atmName}

  # Do filter timestepping stability
  ATM_NCPL=192
  echo "ATM_NCPL: $ATM_NCPL  DTIME: $DTIME"
  if [[ "$DYCORE" == "se" ]]; then
    if [ "$use_nsplit" = true ]; then # if trad. SE nsplit timestepping
      sed -i '/.*se_nsplit/d' user_nl_${atmName}
      if [ "$VALIDSTABVAL" = false ]; then
        # Use se_nsplit = -1, which is internal
        echo "SE_NSPLIT: -1"
        echo "se_nsplit=-1" >> user_nl_${atmName}
      else
        # Add se_nsplit based off of dtime=450 and stability
        SE_NSPLIT=`python -c "from math import ceil; print(int(ceil(450/${STABILITY})))"`
        echo "SE_NSPLIT: $SE_NSPLIT   STABILITY: $STABILITY   "
        echo "se_nsplit=${SE_NSPLIT}" >> user_nl_${atmName}
      fi
    else
      sed -i '/.*se_tstep/d' user_nl_${atmName}
      if [ "$VALIDSTABVAL" = false ]; then
        echo "Betacast cannot handle when $VALIDSTABVAL is false and use_nsplit is false"
        echo "You need to set USERSTAB -> desired se_tstep in the namelist. Exiting..."
        exit
      else
        # Add se_tstep directly from stability
        echo "SE_TSTEP --> STABILITY: $STABILITY   "
        echo "se_tstep=${STABILITY}" >> user_nl_${atmName}
      fi
    fi
  else
    echo "non-SE core, make sure timestepping is happy!"
  fi
  ./xmlchange ATM_NCPL=${ATM_NCPL}

  # Set NHTFRQ, MFILT, and FINCL fields
  sed -i '/.*nhtfrq/d' user_nl_${atmName}
  sed -i '/.*mfilt/d' user_nl_${atmName}
  sed -i '/.*fincl/d' user_nl_${atmName}
  sed -i '/.*empty_htapes/d' user_nl_${atmName}
  sed -i '/.*inithist/d' user_nl_${atmName}
  echo "nhtfrq=1" >> user_nl_${atmName}
  echo "mfilt=200" >> user_nl_${atmName}
  echo "fincl1='U:I','V:I','T:I','PS:I','Q:I','CLDICE:I','CLDLIQ:I','NUMICE:I','NUMLIQ:I','ICEFRAC:I','SNOWHICE:I'" >> user_nl_${atmName}
  echo "empty_htapes=.TRUE." >> user_nl_${atmName}
  echo "inithist='6-HOURLY'" >> user_nl_${atmName}

  ./xmlchange JOB_WALLCLOCK_TIME=${FILTERWALLCLOCK}
  ./xmlchange --force JOB_QUEUE=${FILTERQUEUE}

  if [ $debug = false ] ; then
  
    echo "Begin call to filter-run"
    run_CIME2 "$path_to_rundir" "$CIMEsubstring" "$CIMEbatchargs"

    ## Run NCL filter
    cd $filter_path
    echo "Running filter"
    cp ${sePreFilterIC} ${sePostFilterIC}
    filtfile_name=${casename}.${atmName}.h0.$yearstr-$monthstr-$daystr-$cyclestrsec.nc
    set +e 
    (set -x; ncl lowmemfilter.ncl \
     endhour=${filterHourLength} tcut=${filtTcut} \
    'filtfile_name = "'${path_to_rundir}'/'${filtfile_name}'"' \
    'writefile_name = "'${sePostFilterIC}'"' )
    if [[ $? -ne 9 ]] ; then
      echo "NCL exited with non-9 error code"
      exit 240
    fi
    set -e
  fi  # debug

  echo "done with filter, removing filter files"
  mkdir -p $path_to_nc_files/filtered
  mv -v $path_to_nc_files/*.nc $path_to_nc_files/filtered
  ## Delete filter files that aren't h0 since those are the only ones we care about.
  find $path_to_nc_files/filtered/ -type f -not -name '*.${atmName}.h0.*.nc' | xargs rm

  ## For now, I'm going to go back and delete all the h0 files in filtered
  rm -v $path_to_nc_files/filtered/*.${atmName}.h0.*.nc

  ### Output filter files only, can be used for ensemble or other initialization after the fact
  if ${filterOnly} ; then
    FILTONLYDIR=${path_to_nc_files}/FILT_INIC/${se_yearstr}-${se_monthstr}-${se_daystr}-${se_cyclestrsec}
    mkdir -p ${FILTONLYDIR}
    cp ${sePostFilterIC} ${FILTONLYDIR}/${casename}_FILTERED_${se_yearstr}-${se_monthstr}-${se_daystr}-${se_cyclestrsec}.nc
    cp ${sstFileIC} ${FILTONLYDIR}/${casename}_SST_1x1_${se_yearstr}-${se_monthstr}-${se_daystr}-${se_cyclestrsec}.nc
    mkdir ${FILTONLYDIR}/config_files
    mv ${path_to_nc_files}/*.gz ${FILTONLYDIR}/config_files
    mv ${path_to_nc_files}/*.nml ${FILTONLYDIR}/config_files
    mv ${path_to_nc_files}/*_in ${FILTONLYDIR}/config_files
    exit
  fi

  # Special things we have to do to "reset" CESM after filtering
  echo "Pushing CESM start back from $cyclestrsec to $se_cyclestrsec seconds..."
  cd $path_to_case
  ./xmlchange RUN_STARTDATE=$se_yearstr-$se_monthstr-$se_daystr,START_TOD=$se_cyclestrsec,STOP_OPTION=ndays,STOP_N=$numdays

  # Set new inic to point to post filter file
  SEINIC=${sePostFilterIC}
fi

############################### "ACTUAL" FORECAST RUN ############################### 

## Initial modification of user_nl_atm
sed -i '/.*nhtfrq/d' user_nl_${atmName}
sed -i '/.*mfilt/d' user_nl_${atmName}
sed -i '/.*fincl/d' user_nl_${atmName}
sed -i '/.*empty_htapes/d' user_nl_${atmName}
sed -i '/.*collect_column_output/d' user_nl_${atmName}
sed -i '/.*avgflag_pertape/d' user_nl_${atmName}  # Note, we delete this and user either specifies as :A, :I for each var or as a sep var
echo "empty_htapes=.TRUE." >> user_nl_${atmName}
sed -i '/.*inithist/d' user_nl_${atmName}
if ${save_nudging_files} ; then
  echo "inithist='6-HOURLY'" >> user_nl_${atmName}
else
  echo "inithist='NONE'" >> user_nl_${atmName}
fi

# Concatenate output streams to end of user_nl_${atmName}
cat ${OUTPUTSTREAMS} >> user_nl_${atmName}

# Calculate timestep stability criteria
ATM_NCPL=`python -c "print(int(86400/${DTIME}))"`
echo "ATM_NCPL: $ATM_NCPL  DTIME: $DTIME"
if [[ "$DYCORE" == "se" ]]; then
  if [ "$use_nsplit" = true ]; then # if trad. SE nsplit timestepping
    # Delete existing se_nsplit
    sed -i '/.*se_nsplit/d' user_nl_${atmName}
    if [ "$VALIDSTABVAL" = false ]; then
      echo "SE_NSPLIT: -1 "
      echo "se_nsplit=-1" >> user_nl_${atmName}
    else
      # Add se_nsplit based off of dtime and stability
      SE_NSPLIT=`python -c "from math import ceil; print(int(ceil(${DTIME}/${STABILITY})))"`
      echo "SE_NSPLIT: $SE_NSPLIT   STABILITY: $STABILITY   "
      echo "se_nsplit=${SE_NSPLIT}" >> user_nl_${atmName}
    fi
  else # otherwise, use E3SM se_tstep
    # Delete existing tstep
    sed -i '/.*se_tstep/d' user_nl_${atmName}
    if [ "$VALIDSTABVAL" = false ]; then
      echo "Betacast cannot handle when $VALIDSTABVAL is false and use_nsplit is false"
      echo "You need to set USERSTAB -> desired se_tstep in the namelist. Exiting..."
      exit
    else
      # Add se_tstep directly from stability
      echo "SE_TSTEP --> STABILITY: $STABILITY   "
      echo "se_tstep=${STABILITY}" >> user_nl_${atmName}
    fi
  fi
else
  echo "non-SE core, make sure timestepping is happy!"
fi

./xmlchange ATM_NCPL=${ATM_NCPL}

./xmlchange JOB_WALLCLOCK_TIME=${RUNWALLCLOCK}
./xmlchange --force JOB_QUEUE=${RUNQUEUE}

sed -i '/.*ncdata/d' user_nl_${atmName}
echo "ncdata='${SEINIC}'" >> user_nl_${atmName}

if [ $debug = false ] ; then
  echo "Begin call to forecast run"
  run_CIME2 "$path_to_rundir" "$CIMEsubstring" "$CIMEbatchargs"
fi

fi # end run model

cd $outputdir

# Generate folder structure and move NetCDF files
main_archive "$tmparchivecdir" "$atmName" "$lndName" "$rofName"

# Copy betacast configs to archive directory for posterity
cp -v $MACHINEFILE $NAMELISTFILE $OUTPUTSTREAMS $perturb_namelist $tmparchivecdir/betacast

# Copy user_nl* files to archive dir as well...
cp -v $path_to_case/user* $tmparchivecdir/nl_files

# Archive initial conditions?
if [ $archive_inic = true ]; then
  archive_inic "$tmparchivecdir" "$path_to_case" "$compress_history_nc" "$atmName" "$lndName" "$rofName" "$sstFileIC"
fi

# Compress model output streams
if [ $compress_history_nc = true ]; then
  compress_history "$tmparchivecdir"
fi

# Archive nudging files generated by hindcasts
if [ $save_nudging_files = true ] ; then
  archive_nudging "$tmparchivecdir" "$path_to_nc_files" "$compress_history_nc"
fi

## Move land files to new restart location
cd $path_to_nc_files
if [ $keep_land_restarts = true ]; then
  echo "Archiving land restart files"
  mkdir -p $landdir
  echo "Removing 06Z and 18Z land restart files if they exist"
  rm -v *.${lndName}*.r.*32400.nc || true
  rm -v *.${lndName}*.r.*75600.nc || true
  echo "Moving land restart files for future runs"
  mv -v *.${lndName}*.r.*.nc $landdir || true
  ## Move runoff files to land dir if doing runoff
  if [ $do_runoff = true ]; then
    echo "Removing 06Z and 18Z runoff restart files if they exist"
    rm -v *.${rofName}*.r.*32400.nc || true
    rm -v *.${rofName}*.r.*75600.nc || true
    echo "Moving runoff restart files for future runs"
    mv -v *.${rofName}*.r.*.nc $landdir || true
  fi
else
  echo "Removing all land restart files!"
  rm -v *.${lndName}*.r.*.nc || true
  if [ $do_runoff = true ]; then
    rm -v *.${rofName}*.r.*.nc || true
  fi
fi

## Delete any leftover files in the run dir that we don't need/want anymore
delete_leftovers "$path_to_nc_files" "$atmName" "$lndName" "$rofName"

cd $path_to_nc_files
echo "Moving tmp archive directory to ARCHIVEDIR/YYYYMMDDHH"
if [ -d "${ARCHIVEDIR}/${yearstr}${monthstr}${daystr}${cyclestr}" ]
then
  echo "${ARCHIVEDIR}/${yearstr}${monthstr}${daystr}${cyclestr} already exists! Appending date string!"
  mv -v $tmparchivecdir ${ARCHIVEDIR}/${yearstr}${monthstr}${daystr}${cyclestr}_$(date +\%Y\%m\%d\%H\%M)
else
  mv -v $tmparchivecdir ${ARCHIVEDIR}/${yearstr}${monthstr}${daystr}${cyclestr}
fi

if $dotracking ; then
  
  # Go to cyclone tracking folder...
  cd ${sewxscriptsdir}/cyclone-tracking/
  
  # Set a few things that are hardcoded
  TCVITFOLDER=./fin-tcvitals/
  TCVITFILE=${TCVITFOLDER}/tcvitals.${yearstr}${monthstr}${daystr}${cyclestr}
  ATCFFILE=atcf.${casename}.${yearstr}${monthstr}${daystr}${cyclestr}

  ### Get TC vitals for YYYYMMDDHH and store in $TCVITFOLDER
  (set -x
  /bin/bash ./get-vitals.sh ${yearstr}${monthstr}${daystr}${cyclestr} ${TCVITFOLDER} \
  )
  
  ### If we have vitals, run tracking
  if [ -f ${TCVITFILE} ]; then  # if file does exist (i.e., it was downloaded), run tracker
    (set -x; /bin/bash ./drive-tracking.sh ${yearstr}${monthstr}${daystr}${cyclestr} \
      ${casename} \
      ${TCVITFILE} \
      ${ATCFFILE} \
      ${track_connectfile} \
      ${path_to_rundir} \
      ${track_sendhtml} \
      ${track_hstream} \
      ${track_stride} \
      ${track_ATCFTECH} \
      ${TE_SERIAL_DIR} )
    cp -v trajs.trajectories.txt.${casename}.png ./fin-figs/trajs.trajectories.txt.${casename}.${yearstr}${monthstr}${daystr}${cyclestr}.png
  fi
  
  # Return to where we were...
  cd $path_to_nc_files
fi

if $sendplots ; then
  ### Begin output calls
  echo "Sending plots!"
  cp ${upload_ncl_script} ${upload_ncl_script}.${uniqtime}.ncl
  sed -i 's?.*yearstr=.*?yearstr='${yearstr}'?' ${upload_ncl_script}.${uniqtime}.ncl
  sed -i 's?.*monthstr=.*?monthstr='${monthstr}'?' ${upload_ncl_script}.${uniqtime}.ncl
  sed -i 's?.*daystr=.*?daystr='${daystr}'?' ${upload_ncl_script}.${uniqtime}.ncl
  sed -i 's?.*cyclestrsec=.*?cyclestrsec='${cyclestrsec}'?' ${upload_ncl_script}.${uniqtime}.ncl
  sed -i 's?.*cyclestr=.*?cyclestr='${cyclestr}'?' ${upload_ncl_script}.${uniqtime}.ncl
  (set -x; /bin/bash ${upload_ncl_script}.${uniqtime}.ncl ${nclPlotWeights} ${outputdir}/${yearstr}${monthstr}${daystr}${cyclestr} )
  # Cleanup
  rm ${upload_ncl_script}.${uniqtime}.ncl
fi

### If not live and the run has made it here successively, delete top line of datesfile
if [ $islive = false ] ; then
  cd ${sewxscriptsdir}
  #Remove top line from dates file
  tail -n +2 ${datesfile} > ${datesfile}.2
  mv -v ${datesfile}.2 ${datesfile}
  
  AUTORESUB="yes"
  if [ $AUTORESUB == "yes" ]; then
    echo "*-*-*-* Automatically resubbing next date!" 
    exec ./betacast.sh ${MACHINEFILE} ${NAMELISTFILE} ${OUTPUTSTREAMS}
  fi
fi

echo "DONE"

exit 0
