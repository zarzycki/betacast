## Input is YYYYMMDDHH
## Sets ${yearstr}, ${monthstr}, ${daystr}, ${hourstr}, ${cyclestrsec}
## Ex: parse_time 2022092812
parse_YYYYMMDDHH () {
  echo "Getting time!"
  local thisDate=$1
  yearstr=${thisDate:0:4}
  monthstr=${thisDate:4:2}
  daystr=${thisDate:6:2}
  cyclestr=${thisDate:8:2}
}


## Input is HH int
## Sets $cyclestrsec which is zero padded 3600*hr
get_cyclestrsec () {
  local hrstr=$1
  ## Figure out the seconds which correspond to the cycle and zero pad if neces
  cyclestrsec=$(($hrstr*3600))
  while [ ${#cyclestrsec} -lt 5 ];
  do
    cyclestrsec="0"$cyclestrsec
  done
}


### This is the OG version of this that figures out what the status is by
### counting the number of gzipped files in the run directory
### when the gzip files increments, it means the job is terminated.
run_CIME () {
  # 1 path_to_run_dir
  # 2 ${CIMEsubstring}
  # 3 ${CIMEbatchargs}
  
  echo "run_CIME: path_to_rundir: "$1
  echo "run_CIME: CIMEsubstring: "$2
  echo "run_CIME: CIMEbatchargs: "$3

  # Get number of log .gz files for sleeping
  echo "Running again!" > $1/testrunning.gz
  local numlogfiles=`ls ${1}/*.gz | wc -l`
  echo "run_CIME: numlogfiles: $numlogfiles"

  echo "SUBMITTING FORECAST RUN"
  set +e ; ./case.submit $2 --batch-args "${3}" ; set -e

  ## Set up NUKEing
  if [ -f "${1}/NUKE" ] ; then rm -v $1/NUKE ; sleep 5 ; fi
  echo "To NUKE, run \"touch ${1}/NUKE\" "

  ## Hold script while log files from filter run haven't been archived yet
  while [ `ls ${1}/*.gz | wc -l` == $numlogfiles ]
  do
    if [ -f "${1}/NUKE" ] ; then echo "Nuke sequence initiated, exiting betacast" ; exit ; fi
    sleep 10 ; echo "Sleeping... $(date '+%Y%m%d %H:%M:%S')"
  done
  echo "Run over done sleeping ($(date '+%Y%m%d %H:%M:%S')) will hold for 10 more sec to make sure files moved"
  sleep 10
}


### This is the v2 version of this that figures out what the status is by reading
### the last valid entry in CaseStatus
run_CIME2 () {
  # 1 path_to_run_dir
  # 2 ${CIMEsubstring}
  # 3 ${CIMEbatchargs}
  
  echo "run_CIME: path_to_rundir: "$1
  echo "run_CIME: CIMEsubstring: "$2
  echo "run_CIME: CIMEbatchargs: "$3

  # Get number of log .gz files for sleeping

  local CASESTR=""
  
  echo "SUBMITTING FORECAST RUN"
  set +e ; ./case.submit $2 --batch-args "${3}" ; set -e

  ## Set up NUKEing
  if [ -f "${1}/NUKE" ] ; then rm -v $1/NUKE ; sleep 5 ; fi
  echo "To NUKE, run \"touch ${1}/NUKE\" "

  ## Hold script while log files from filter run haven't been archived yet
  STATUS=1
  while [ $STATUS == 1 ]
  do
    if [ -f "${1}/NUKE" ] ; then echo "Nuke sequence initiated, exiting betacast" ; exit ; fi
    
    # Get the last valid line from the CaseStatus file...
    CASESTR=`grep "^20" CaseStatus | tail -1`
    
    if [[ "$CASESTR" == *"case.run success"* ]]; then
      STATUS=0
    elif [[ "$CASESTR" == *"case.run error"* ]]; then
      STATUS=99
    else
      STATUS=1
    fi
    
    sleep 10 ; echo "Sleeping... STATUS: $STATUS -- $(date '+%Y%m%d %H:%M:%S')"
  done
  
  if [ $STATUS -eq 0 ]; then
    echo "Run over done sleeping ($(date '+%Y%m%d %H:%M:%S')) will hold for 30 more sec to make sure files moved"
    sleep 30
  else
    echo "Uh oh, something wrong!"
  fi
}


compress_history () {
  (
  # Compress files using lossless compression
  echo "Compressing model history files..."
  cd $1
  for f in *.nc ; do echo "Compressing $f" ; ncks -4 -L 1 -O $f $f ; done
  )
}


#$path_to_nc_files $2
#$tmparchivecdir $1
#$compress_history $3
archive_nudging () {
  (
  echo "Archiving nudging files..."
  cd $2
  mkdir -p $1/nudging
  mv -v *.i.*.nc $1/nudging
  if [ $3 = true ]; then
    # Compress files using lossless compression
    compress_history $1/nudging
  fi
  )
}

#tmparchivecdir = 1
#path_to_case =2 
#compress_history_nc = 3
#atmName = 4
#lndName = 5
#rofName = 6
#sstFileIC = 7

archive_inic () {

  echo "Archiving initial condition files..."
  mkdir -p $1/inic

  # Copy LND initial conditions
  ARCFILE=`grep ^finidat $2/user_nl_$5 | cut -d "=" -f2`
  strip_quotes ARCFILE
  echo "Found initial file: "$ARCFILE
  cp -v $ARCFILE $1/inic

  # Copy ATM initial conditions
  ARCFILE=`grep ^ncdata $2/user_nl_$4 | cut -d "=" -f2`
  strip_quotes ARCFILE
  echo "Found initial file: "$ARCFILE
  cp -v $ARCFILE $1/inic

  # Copy ROF initial conditions
  if [ $do_runoff = true ]; then
    ARCFILE=`grep ^finidat_rtm $2/user_nl_$6 | cut -d "=" -f2`
    strip_quotes ARCFILE
    echo "Found initial file: "$ARCFILE
    cp -v $ARCFILE $1/inic
  fi

  # Copy SST conditions
  cp -v $7 $1/inic

  if [ $3 = true ]; then
    compress_history "$1/inic"
  fi
}


#tmparchivecdir = 1
#atmName = 2
#lndName = 3
#rofName = 4
main_archive () {
  echo "Creating archive folder data structure"
  mkdir -p $1
  mkdir -p $1/images
  mkdir -p $1/text
  mkdir -p $1/nl_files
  mkdir -p $1/logs
  mkdir -p $1/betacast

  set +e

  echo "Moving relevant files to archive folder"
  mv -v *.$2.h*.nc $1
  mv -v *.$3*.h*.nc $1
  mv -v *.$4*.h*.nc $1
  cp -v *_in seq_maps.rc *_modelio.nml docn.streams.txt.prescribed $1/nl_files
  mv -v *.txt $1/text
  mv -v *.log.* $1/logs
  mv -v timing.*.gz $1/logs
  mv -v atm_chunk_costs*.gz $1/logs

  mv -v timing/ $1/
  
  set -e
}