#!/bin/bash

COUNTER=0
DTIME=90

outputdir=/home/zarzycki/tcforecast_68_x4/run/test

cd $outputdir

while [  $COUNTER -lt 900 ]; do
  echo The starting counter is $COUNTER
  filenames=`ls tc*h1*.nc`
  numfiles=`ls tc*h1*.nc | wc -l`
  if [ $numfiles -eq 0 ]
  then
    echo "Nothing found"
    let COUNTER=COUNTER+DTIME
  else
    echo "Found at least one file"
    rename 's/^/_/' tc*h1*.nc
    newfiles=`ls __tc*h1*.nc`
    for f in $newfiles
    do
	    echo "Processing $f"
	    #ncl $f
    done
    let COUNTER=0
  fi 
done