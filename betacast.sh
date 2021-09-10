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

date

set -e
#set -v

# Set files, in reality order doesn't matter
MACHINEFILE=${1}
NAMELISTFILE=${2}
OUTPUTSTREAMS=${3}
# If relative path, convert to absolute path
if [[ "$MACHINEFILE" != /* ]] && [[ "$MACHINEFILE" != ~* ]]; then MACHINEFILE=${PWD}/${MACHINEFILE}; fi
if [[ "$NAMELISTFILE" != /* ]] && [[ "$NAMELISTFILE" != ~* ]]; then NAMELISTFILE=${PWD}/${NAMELISTFILE}; fi
if [[ "$OUTPUTSTREAMS" != /* ]] && [[ "$OUTPUTSTREAMS" != ~* ]]; then OUTPUTSTREAMS=${PWD}/${OUTPUTSTREAMS}; fi
echo $MACHINEFILE; echo $NAMELISTFILE; echo $OUTPUTSTREAMS

# Sanitize namelist files (add carriage return to end)
sed -i -e '$a\' ${MACHINEFILE}
sed -i -e '$a\' ${NAMELISTFILE}
sed -i -e '$a\' ${OUTPUTSTREAMS}
#'

# Note, ___ will be converted to a space. Namelists cannot have whitespace due to
# parsing on whitespaces...
echo "Reading namelist ${NAMELISTFILE}..."
inputstream=`cat ${NAMELISTFILE} ${MACHINEFILE} | grep -v "^#"`
#echo $inputstream
set -- $inputstream
while [ $1 ]
 do
  echo "NAMELIST: setting ${1} to ${3//___/ }"
  #eval $1=$3
  eval $1="${3//___/ }"  
  shift 3
 done

set -u  # turn on crashes for unbound variables in bash

###################################################################################
############### OPTIONAL TO BE SET BY USER ########################################
path_to_nc_files=${path_to_rundir}              # Path where .nc files are
outputdir=${path_to_rundir}                     # Path where .nc files are being written
archivedir=${path_to_rundir}/proc/              # Path to temporarily stage final data
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
if [ -z ${do_runoff+x} ]; then do_runoff=0; fi
if [ -z ${keep_land_restarts+x} ]; then keep_land_restarts=1; fi
if [ -z ${perturb_namelist+x} ]; then perturb_namelist=""; fi

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

### ERROR CHECKING BLOCK! #########################################################

if [ "${add_perturbs}" = true ] && { [ -z "$perturb_namelist" ] || [ ! -f "$perturb_namelist" ]; } ; then
  echo "add_perturbs is true but can't find namelist: "$perturb_namelist
  exit 1
fi

###################################################################################

# do some stability calcs
# if USERSTAB is negative, use internal calcs.
# If positive, use the value in seconds for dt_dyn
USERSTABTF=`python -c "print('TRUE' if ${USERSTAB} > 0 else 'FALSE')"`
if [ ${USERSTABTF} == 'FALSE' ] ; then
  STABILITY=`python -c "print(30./${FINERES}*450.)"`
else
  STABILITY=${USERSTAB}
fi
echo "Dynamic stability for ne${FINERES} to be ${STABILITY} seconds"

## Create paths to generate initial files if they don't exist...
mkdir -p ${pathToINICfiles}
mkdir -p ${pathToSSTfiles}

# Set timestamp for backing up files, etc.
timestamp=`date +%Y%m%d.%H%M`
uniqtime=`date +"%s%N"`

echo "We are using ${casename} for the case"

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
    echo "12Z cycle"
    monthstr=`date --date="yesterday" -u +%m`
    daystr=`date --date="yesterday" -u +%d`
    yearstr=`date --date="yesterday" -u +%Y`
    twodaysago=`date --date='3 days ago' -u +"%Y%m%d"`
    cyclestr=12
  elif [ $currtime -lt 0928 ]
  then
    echo "00Z cycle"
    cyclestr=00
  elif [ $currtime -lt 1528 ]
  then
    echo "00Z cycle"
    cyclestr=00
  elif [ $currtime -lt 2128 ]
  then
    echo "12Z cycle"
    cyclestr=12
  elif [ $currtime -ge 2128 ]
  then
    echo "12Z cycle"
    cyclestr=12
  else
    echo "Can't figure out start time"
    exit 1
  fi

##CMZ 
#cyclestr=00
else     # if not live, draw from head of dates.txt file
  datesfile=dates.${casename}.txt
  longdate=$(head -n 1 ${datesfile})
  yearstr=${longdate:0:4}
  monthstr=${longdate:4:2}
  daystr=${longdate:6:2}
  cyclestr=${longdate:8:2}
  echo "From datesfile, read in: "$yearstr' '$monthstr' '$daystr' '$cyclestr'Z'

  # Do some simple error trapping on date string to ensure validity
  if [ -z "$longdate" ]; then
    echo "Date string passed in is empty, exiting..." ; exit 91
  fi
  if [ ${#longdate} -ne 10 ]; then 
    echo "Malformed date string, $longdate is ${#longdate} characters, needs 10 (YYYYMMDDHH). Exiting..." ; exit 92
  fi
  if [[ -n $(echo $longdate | tr -d '[0-9]') ]] ; then
    echo "Malformed date string, $longdate contains non-numeric values. Exiting..." ; exit 93
  fi
  if (( yearstr > 3000 || yearstr < 1 )); then
    echo "Year set to $yearstr, this sounds wrong, exiting..." ; exit 94
  fi
  if (( cyclestr > 23 )); then
    echo "Cycle string set to $cyclestr Z, this sounds wrong, exiting..." ; exit 95
  fi
fi

## Figure out the seconds which correspond to the cycle and zero pad if neces
cyclestrsec=$(($cyclestr*3600))
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

if [ $debug -ne 1 ] ; then

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
      if [ $islive -ne 0 ] ; then
        rm -f gfs.t*pgrb2f00*
        gfsFTPPath=ftp://ftp.ncep.noaa.gov/pub/data/nccf/com/gfs/prod/gfs.$yearstr$monthstr$daystr/$cyclestr/atmos/
        #gfsFTPFile='gfs.t'$cyclestr'z.pgrb2f00'
        gfsFTPFile='gfs.t'$cyclestr'z.pgrb2.0p25.anl'
        #gfsFTPFile='gfs.t'$cyclestr'z.pgrb2.0p50.anl'
        echo "Attempting to download ${gfsFTPPath}${gfsFTPFile}"
        ## Scrape for files
        error=1
        while [ $error != 0 ] ; do
          wget --read-timeout=30 $gfsFTPPath$gfsFTPFile
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
        wget --load-cookies ~/.thecookies http://rda.ucar.edu/data/ds093.0/${yearstr}/${CFSRFILENAME}
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
      if [[ "$MACHINEFILE" != *cheyenne* ]]; then
        echo "support broken for auto download ERA, please prestage!"
        exit
      else
        echo "We are on Cheyenne, so even though we lack a local file, we can use RDA"
        ERA5RDA=1
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
    if [ $islive -ne 0 ] ; then
      # Here is where we get the "live" GDAS SST file
      rm -f gdas1*sstgrb*
      sstFTPPath=ftp://ftp.ncep.noaa.gov/pub/data/nccf/com/gfs/prod/sst.${yestyearstr}${yestmonthstr}${yestdaystr}/
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
      wget --read-timeout=30 -nv $sstFTPPath$sstFTPFile
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

  set +e
  cd ${sst_to_cam_path}
  (set -x; ncl sst_interp.ncl 'initdate="'${yearstr}${monthstr}${daystr}${cyclestr}'"' \
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
  ncatted -O -a units,time,o,c,"days since 0001-01-01 00:00:00" ${sstFileIC} ${sstFileIC}

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
       'data_filename = "'$gfs_files_path'/gfs_atm_'$yearstr$monthstr$daystr$cyclestr'.grib2"'  \
       'wgt_filename="'${gfs2seWeights}'"' \
       'model_topo_file="'${adjust_topo-}'"' \
       'adjust_config="'${adjust_flags-}'"' \
       'se_inic = "'${sePreFilterIC}'"' )
  elif [ $atmDataType -eq 2 ] ; then
    echo "CD ing to ERA-interim interpolation directory"
    cd $atm_to_cam_path
    (set -x; ncl -n atm_to_cam.ncl 'datasource="ERAI"'     \
        numlevels=${numLevels} \
        YYYYMMDDHH=${yearstr}${monthstr}${daystr}${cyclestr} \
       'data_filename = "/glade/p/work/zarzycki/getECMWFdata/ERA-Int_'$yearstr$monthstr$daystr$cyclestr'.nc"'  \
       'wgt_filename="/glade/p/work/zarzycki/getECMWFdata/ERA_to_uniform_60_patch.nc"' \
       'se_inic = "'${sePreFilterIC}'"' )
  elif [ $atmDataType -eq 3 ] ; then
    echo "CD ing to interpolation directory"
    cd $atm_to_cam_path 
    (set -x; ncl -n atm_to_cam.ncl 'datasource="CFSR"'     \
      numlevels=${numLevels} \
      YYYYMMDDHH=${yearstr}${monthstr}${daystr}${cyclestr} \
     'data_filename = "'$gfs_files_path'/cfsr_atm_'$yearstr$monthstr$daystr$cyclestr'.grib2"'  \
     'wgt_filename="'${gfs2seWeights}'"' \
     'model_topo_file="'${adjust_topo-}'"' \
     'adjust_config="'${adjust_flags-}'"' \
     'se_inic = "'${sePreFilterIC}'"' )
  elif [ $atmDataType -eq 4 ] ; then
    echo "CD ing to ERA5 interpolation directory"
    cd $atm_to_cam_path
    if [ $ERA5RDA -eq 1 ] ; then
      RDADIR=/glade/collections/rda/data/
      (set -x; ncl -n atm_to_cam.ncl 'datasource="ERA5RDA"'     \
          numlevels=${numLevels} \
          YYYYMMDDHH=${yearstr}${monthstr}${daystr}${cyclestr} \
         'data_filename = "'${RDADIR}'/ds633.0/e5.oper.invariant/197901/e5.oper.invariant.128_129_z.ll025sc.1979010100_1979010100.nc"'  \
         'wgt_filename="'${gfs2seWeights}'"' \
         'model_topo_file="'${adjust_topo-}'"' \
         'adjust_config="'${adjust_flags-}'"' \
         'se_inic = "'${sePreFilterIC}'"' )
    else
      (set -x; ncl -n atm_to_cam.ncl 'datasource="ERA5"'     \
          numlevels=${numLevels} \
          YYYYMMDDHH=${yearstr}${monthstr}${daystr}${cyclestr} \
         'data_filename = "'$era_files_path'/ERA5_'$yearstr$monthstr$daystr$cyclestr'.nc"'  \
         'wgt_filename="'${gfs2seWeights}'"' \
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
  if [[ $? -ne 9 ]]
  then
    echo "NCL exited with non-9 error code"
    exit 240
  fi
  echo "ATM NCL completed successfully"
  set -e # Turn error checking back on
  
fi #End debug if statement

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

if [ "${add_perturbs}" = true ] ; then
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

############################### SETUP AND QUERY ############################### 

cd $path_to_case
DYCORE=`./xmlquery CAM_DYCORE | sed 's/^[^\:]\+\://' | xargs`
echo "DYCORE: "$DYCORE

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
      cp user_nl_${lndName}_presave user_nl_${lndName}
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

if [ $do_runoff -ne 0 ]; then

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
  SE_NSPLIT=`python -c "from math import ceil; print(int(ceil(450/${STABILITY})))"`
  echo "ATM_NCPL: $ATM_NCPL  DTIME: $DTIME"
  if [[ "$DYCORE" == "se" ]]; then
    echo "SE_NSPLIT: $SE_NSPLIT   STABILITY: $STABILITY   "
    sed -i '/.*se_nsplit/d' user_nl_${atmName}
    echo "se_nsplit=${SE_NSPLIT}" >> user_nl_${atmName}
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

  if [ $debug -ne 1 ] ; then
    echo "Begin call to filter-run"

    # Get number of log .gz files for sleeping
    echo "Running again!" > ${path_to_rundir}/testrunning.gz
    numlogfiles=`ls ${path_to_rundir}/*.gz | wc -l`
    echo "numlogfiles: $numlogfiles"

    echo "SUBMITTING FILTER RUN"
    if $usingCIME ; then
      set +e ; ./case.submit ${CIMEsubstring} --batch-args "${CIMEbatchargs}" ; set -e
    else
      # needs to be modified for your machine if not using CIME
      bsub < ${casename}.run
    fi

    ## Hold script while log files from filter run haven't been archived yet
    while [ `ls ${path_to_rundir}/*.gz | wc -l` == $numlogfiles ]
    do
      sleep 20 ; echo "Sleeping... $(date '+%Y%m%d %H:%M:%S')"
    done
    echo "Run over done sleeping ($(date '+%Y%m%d %H:%M:%S')) will hold for 30 more sec to make sure files moved" ; sleep 40

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

sed -i '/.*nhtfrq/d' user_nl_${atmName}
sed -i '/.*mfilt/d' user_nl_${atmName}
sed -i '/.*fincl/d' user_nl_${atmName}
sed -i '/.*empty_htapes/d' user_nl_${atmName}
sed -i '/.*inithist/d' user_nl_${atmName}
echo "empty_htapes=.TRUE." >> user_nl_${atmName}
echo "inithist='NONE'" >> user_nl_${atmName}

# Concatenate output streams to end of user_nl_${atmName}
cat ${OUTPUTSTREAMS} >> user_nl_${atmName}

# Calculate timestep criteria
ATM_NCPL=`python -c "print(int(86400/${DTIME}))"`
SE_NSPLIT=`python -c "from math import ceil; print(int(ceil(${DTIME}/${STABILITY})))"`
echo "ATM_NCPL: $ATM_NCPL  DTIME: $DTIME"
if [[ "$DYCORE" == "se" ]]; then
  echo "SE_NSPLIT: $SE_NSPLIT   STABILITY: $STABILITY   "
  sed -i '/.*se_nsplit/d' user_nl_${atmName}
  echo "se_nsplit=${SE_NSPLIT}" >> user_nl_${atmName}
else
  echo "non-SE core, make sure timestepping is happy!"
fi
./xmlchange ATM_NCPL=${ATM_NCPL}

./xmlchange JOB_WALLCLOCK_TIME=${RUNWALLCLOCK}
./xmlchange --force JOB_QUEUE=${RUNQUEUE}

sed -i '/.*ncdata/d' user_nl_${atmName}
echo "ncdata='${SEINIC}'" >> user_nl_${atmName}

if [ $debug -ne 1 ]
then
  echo "Begin call to forecast run"

  # Get number of log .gz files for sleeping
  echo "Running again!" > ${path_to_rundir}/testrunning.gz
  numlogfiles=`ls ${path_to_rundir}/*.gz | wc -l`
  echo "numlogfiles: $numlogfiles"

  echo "SUBMITTING FORECAST RUN"
  if $usingCIME ; then
    set +e ; ./case.submit ${CIMEsubstring} --batch-args "${CIMEbatchargs}" ; set -e
  else
    # needs to be modified for your machine if not using CIME
    bsub < ${casename}.run
  fi
  
  ## Hold script while log files from filter run haven't been archived yet
  while [ `ls ${path_to_rundir}/*.gz | wc -l` == $numlogfiles ]
  do
    sleep 20 ; echo "Sleeping... $(date '+%Y%m%d %H:%M:%S')"
  done
  echo "Run over done sleeping ($(date '+%Y%m%d %H:%M:%S')) will hold for 30 more sec to make sure files moved"
  sleep 30
fi

fi # end run model

cd $outputdir

echo "Creating archive folder data structure"
mkdir -p $archivedir
mkdir -p $archivedir/images
mkdir -p $archivedir/text
mkdir -p $archivedir/nl_files
mkdir -p $archivedir/logs

set +e

echo "Moving relevant files to archive folder"
mv -v *.${atmName}.h*.nc $archivedir
mv -v *.${lndName}*.h*.nc $archivedir
mv -v *.${rofName}*.h*.nc $archivedir
cp -v *_in seq_maps.rc *_modelio.nml docn.streams.txt.prescribed $archivedir/nl_files
mv -v *.txt $archivedir/text
mv -v *.log.* $archivedir/logs
mv -v timing.*.gz $archivedir/logs
mv -v atm_chunk_costs*.gz $archivedir/logs

mv -v timing/ $archivedir/

## Move land files to new restart location
cd $path_to_nc_files
if [ $keep_land_restarts -eq 1 ]; then
  echo "Archiving land restart files"
  mkdir $landdir
  echo "Removing 06Z and 18Z land restart files if they exist"
  rm -v *.${lndName}*.r.*32400.nc
  rm -v *.${lndName}*.r.*75600.nc
  echo "Moving land restart files for future runs"
  mv -v *.${lndName}*.r.*nc $landdir
  rm -v *.${lndName}*.r.*.nc
  ## Move runoff files to land dir if doing runoff
  if [ $do_runoff -ne 0 ]; then
    echo "Removing 06Z and 18Z runoff restart files if they exist"
    rm -v *.${rofName}*.r.*32400.nc
    rm -v *.${rofName}*.r.*75600.nc
    echo "Moving runoff restart files for future runs"
    mv -v *.${rofName}*.r.*nc $landdir
    rm -v *.${rofName}*.r.*.nc
  fi
else
  echo "Removing all land restart files!"
  rm -v *.${lndName}*.r.*.nc
  rm -v *.${lndName}*.r.*.nc
  if [ $do_runoff -ne 0 ]; then
    rm -v *.${rofName}*.r.*.nc
    rm -v *.${rofName}*.r.*.nc
  fi
fi

echo "Deleting restart/misc. files produced by CESM that aren't needed"
cd $path_to_nc_files
rm -v *.${lndName}*.rh0.*.nc
rm -v *.docn.rs1.*.bin
rm -v *.${atmName}.r.*.nc
rm -v *.${atmName}.rs.*.nc
rm -v *.cpl.r.*.nc
rm -v *.${atmName}.rh3.*.nc
rm -v *.${rofName}.rh0.*.nc
rm -v *.cice.r.*.nc
rm -v rpointer.*
rm -v *.bin
rm -v *.h.*.nc
rm -v *initial_hist*.nc

echo "Renaming tmp archive directory to YYYYMMDDHH"
mv -v $archivedir ${outputdir}/${yearstr}${monthstr}${daystr}${cyclestr}

if $dotracking ; then
  #yearstr="2020"
  #monthstr="07"
  #daystr="31"
  #cyclestr="12"
  #casename="forecast_nhemitc_30_x4_CAM5_L30"
  #sewxscriptsdir=~/sw/betacast/
  cd ${sewxscriptsdir}/cyclone-tracking/
  TCVITFOLDER=./fin-tcvitals/
  TCVITFILE=${TCVITFOLDER}/tcvitals.${yearstr}${monthstr}${daystr}${cyclestr}
  mkdir -p ${TCVITFOLDER}
  mkdir -p ./fin-figs/
  mkdir -p ./fin-atcf/
  ATCFFILE=atcf.${casename}.${yearstr}${monthstr}${daystr}${cyclestr}
  if [ ! -f ${TCVITFILE} ]; then   #if TCVITFILE doesn't exist, download
    wget http://hurricanes.ral.ucar.edu/repository/data/tcvitals_open/${yearstr}
    grep "${yearstr}${monthstr}${daystr} ${cyclestr}00" ${yearstr} > ${TCVITFILE}
    rm ${yearstr}
  fi
  if [ -f ${TCVITFILE} ]; then  # if file does exist (i.e., it was downloaded), run tracker
    /bin/bash ./drive-tracking.sh ${yearstr}${monthstr}${daystr}${cyclestr} ${casename} ${TCVITFILE} ${ATCFFILE}
    cp trajs.trajectories.txt.${casename}.png ./fin-figs/trajs.trajectories.txt.${casename}.${yearstr}${monthstr}${daystr}${cyclestr}.png
  fi
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
  /bin/bash ${upload_ncl_script}.${uniqtime}.ncl ${nclPlotWeights} ${outputdir}/${yearstr}${monthstr}${daystr}${cyclestr}
  # Cleanup
  rm ${upload_ncl_script}.${uniqtime}.ncl
fi

### If not live and the run has made it here successively, delete top line of datesfile
if [ $islive -eq 0 ] ; then
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

exit 0
