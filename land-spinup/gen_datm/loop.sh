#!/bin/bash

## declare an array variable
declare -a months=("01" "02" "03" "04" "05" "06" "07" "08" "09" "10" "11" "12")

STYR=2020
ENYR=2020

## now loop through the above array
for ii in "${months[@]}"
do
  for jj in $(seq $STYR $ENYR);
  do
    echo $jj' '$ii
    python get-data.py "${jj}" "${ii}"
  done
done
