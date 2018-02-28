#!/bin/bash

array=(FLUT PRECT WIND U850 PRECL)

echo "Array size: ${#array[*]}"

echo "Array items:"
for item in ${array[*]}
do
    printf "   %s\n" $item
done

echo "Array indexes:"
for index in ${!array[*]}
do
    printf "   %d\n" $index
done

echo "Array items and indexes:"
for index in ${!array[*]}
do
    printf "%4d: %s\n" $index ${array[$index]}
done

  for item in ${array[*]}
    do
        printf "   %s\n" $item
        
        echo ${item}files.txt
              
      
    done
