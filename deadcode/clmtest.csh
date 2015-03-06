#!/bin/bash

machineid=3

gridname=tcforecast_60_x4
numHoursSEStart=3

date

## Here we get two digit strings for UTC time for month, day, year
## We also get current time in hoursminutes (because the GFS output lags by 3.5 hours)
monthstr=`date -u +%m`
daystr=`date -u +%d`
yearstr=`date -u +%Y`
currtime=`date -u +%H%M`

## Use currtime to figure out what is the latest cycle we have access to
if [ $currtime -lt 0330 ]
then
  echo "18Z cycle"
  monthstr=`date --date="yesterday" -u +%m`
  daystr=`date --date="yesterday" -u +%d`
  yearstr=`date --date="yesterday" -u +%Y`
  cyclestr=18
elif [ $currtime -lt 0930 ]
then
  echo "00Z cycle"
  cyclestr=00
elif [ $currtime -lt 1530 ]
then
  echo "06Z cycle"
  cyclestr=06
elif [ $currtime -lt 2130 ]
then
  echo "12Z cycle"
  cyclestr=12
elif [ $currtime -ge 2130 ]
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


cd ~/tcforecast_60_x4

if [ -f /home/zarzycki/${gridname}/run/clmstart/${gridname}.clm2.r.${se_yearstr}-${se_monthstr}-${se_daystr}-${se_cyclestrsec}.nc ]
then
    echo "Land file exists, sedding that into CLM namelist"
    sed -i 's?!finidat.*?finidat='"'"/home/zarzycki/"${gridname}"/run/clmstart/"${gridname}".clm2.r."${se_yearstr}"-"${se_monthstr}"-"${se_daystr}"-"${se_cyclestrsec}".nc"'"'?' user_nl_clm
else
    echo "Land file DOES NOT EXIST, will use arbitrary CESM spinup"
fi


