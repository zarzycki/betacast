#!/bin/bash

# Case directory set
DRIVER_DIRECTORY=/glade/u/home/zarzycki/betacast/atm_to_cam/tcseed/
path_to_case=/glade/u/home/zarzycki/aqua/APE3-test/
path_to_rundir=/glade/u/home/zarzycki/scratch/APE3-test/run

# Tempeststuff
TEMPESTEXTREMESDIR=~/work/tempestextremes_noMPI/
CONNECTDAT="/glade/u/home/zarzycki/tempest-scripts/hyperion/ne30.connect_v2.dat"
TEMPESTFILE=${path_to_rundir}/cyclones_tempest
TEMPESTTMP=${DRIVER_DIRECTORY}/cyclones_tempest.tmp

#Vortex settings
vortex_namelist=/glade/u/home/zarzycki/betacast/atm_to_cam/tcseed/seed.ape.nl
do_seed=true
NSEEDS=3
#vortex_namelist=/glade/u/home/zarzycki/betacast/atm_to_cam/tcseed/unseed.ape.nl
#do_seed=false
#NSEEDS=0

#misc settings
WHICHSED="sed"

############################### RUN MODEL ############################### 

if [ -f .SEED_DONE ]; then
  echo "Exiting because of done..."
  exit
fi

############################### RUN MODEL ############################### 
  
echo "Setting up settings!"

cd ${path_to_case}

# turn off restarts
# turn off archiving
./xmlchange DOUT_S=FALSE

echo "Begin call to model run"

# Get number of log .gz files for sleeping
echo "Running again!" > ${path_to_rundir}/testrunning.gz
numlogfiles=`ls ${path_to_rundir}/*.gz | wc -l`
echo "numlogfiles: $numlogfiles"

echo "SUBMITTING FORECAST RUN"
set +e ; ./case.submit ${CIMEsubstring} --batch-args "${CIMEbatchargs}" ; set -e

## Set up NUKEing
if [ -f "${path_to_rundir}/NUKE" ] ; then rm -v ${path_to_rundir}/NUKE ; sleep 5 ; fi
echo "To NUKE, run \"touch ${path_to_rundir}/NUKE\" "

## Hold script while log files from filter run haven't been archived yet
while [ `ls ${path_to_rundir}/*.gz | wc -l` == $numlogfiles ]
do
  if [ -f "${path_to_rundir}/NUKE" ] ; then echo "Nuke sequence initiated, exiting betacast" ; exit ; fi
  sleep 10 ; echo "Sleeping... $(date '+%Y%m%d %H:%M:%S')"
done
echo "Run over done sleeping ($(date '+%Y%m%d %H:%M:%S')) will hold for 10 more sec to make sure files moved"
sleep 10

############################### RESUB SCRIPT ############################### 

cd $DRIVER_DIRECTORY

# Get most recent file for tracking
MOSTRECENT=`ls --color=never ${path_to_rundir}/*.cam.h0.*.nc -r | head -n 1`
# TempestExtremes
STR_DETECT="--verbosity 0 --in_connect ${CONNECTDAT} --out ${TEMPESTFILE} --closedcontourcmd PSL,300.0,5.0,0;_DIFF(Z300,Z500),-6.0,5.0,1.0 --minlat -25.0 --maxlat 25.0 --mergedist 6.0 --searchbymin PSL --outputcmd PSL,min,0"
${TEMPESTEXTREMESDIR}/bin/DetectNodes --in_data ${MOSTRECENT} ${STR_DETECT} 

# Get number of lines in Tempest file
NLINES=`head ${TEMPESTFILE} | wc -l`
echo $NLINES

# If we have more than one storm *or* no storms but do seed...
if (( $NLINES > 1 )) || ${do_seed} ; then

  # Get most recent restart file
  RESTARTFILE=`ls --color=never ${path_to_rundir}/*.cam.r.*.nc -r | head -n 1`
  cp -v $RESTARTFILE ${path_to_rundir}/ORIG.nc

  # If NLINES>1, we have some storms in the TMP file, so load them into bash arrays
  if (( $NLINES > 1 )) ; then
    tail -n +2 ${TEMPESTFILE} > ${TEMPESTTMP}
    declare stormLat=()
    declare stormLon=()
    while IFS=$'\t' read col1 col2 col3 col4
    do 
      echo "$col2"
      echo "$col3"
      stormLon+=($col2)
      stormLat+=($col3)
    done < ${TEMPESTTMP}
    echo "${stormLat[*]}"
    echo "${stormLon[*]}"
    NSTORMS=`echo "${#stormLat[@]}"`
  else
    NSTORMS=0
  fi
  
  # If we are seeding, the loop is number of seeds to add.
  # If we are unseeding, loop is number of storms to remove from TMP file
  if ${do_seed} ; then
    NLOOP=$NSEEDS
  else
    NLOOP=$NSTORMS
  fi
  echo "NLOOP is: $NLOOP"

  for ((ii=0; ii<NLOOP; ii++)); do
    echo "doing loop index: $ii"
    
    # First figure out if seeding, where this seed goes
    # ... or if not seeding, get rid of ii'th storm
    if ${do_seed} ; then
      echo "WE ARE SEEDING!"
      set +e
      (set -x; ncl -n random-seed.ncl 'filename= "'${TEMPESTTMP}'"' 'pthi = "'${vortex_namelist}'"' )
      if [[ $? -ne 9 ]] ; then echo "NCL exited with non-9 error code" ; exit 240 ; fi
      set -e
    else
      echo "sedding lat/lon: ${stormLat[$ii]} ${stormLon[$ii]}"
      $WHICHSED -i "s?.*psminlat=.*?psminlat=${stormLat[$ii]}?" ${vortex_namelist}
      $WHICHSED -i "s?.*psminlon=.*?psminlon=${stormLon[$ii]}?" ${vortex_namelist}
    fi

    # Run usual code
    set +e
    (set -x; ncl -n find-tc-fill-params.ncl 'inic_file= "'${RESTARTFILE}'"' 'pthi = "'${vortex_namelist}'"' )
    if [[ $? -ne 9 ]] ; then echo "NCL exited with non-9 error code" ; exit 240 ; fi
    (set -x; ncl -n seed-tc-in-ncdata.ncl   'seedfile = "'${RESTARTFILE}'"' 'pthi = "'${vortex_namelist}'"' )
    if [[ $? -ne 9 ]] ; then echo "NCL exited with non-9 error code" ; exit 240 ; fi
    set -e
    
  done
  
  # Remove TMP traj file
  rm -v ${TEMPESTTMP}
else
  # If here it means we are *not* seeding and there were no detected storms (i.e., NLINES=1) that need filling
  echo "Doing nothing!"
fi

############################### go back and do some things... ############################### 

cd ${path_to_case}

./xmlchange CONTINUE_RUN=TRUE
#set +e ; ./case.st_archive ; set -e

############################### RESUB SCRIPT ############################### 

cd $DRIVER_DIRECTORY ; pwd

exec ./master.sh
