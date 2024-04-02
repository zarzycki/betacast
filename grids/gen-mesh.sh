#!/bin/bash -l

BASERES=48
REFINE_LEVEL=2
#BASERES=128
#REFINE_LEVEL=3
SQUADGEN=~/Software/squadgen/

# Define from outer to inner
#LATWIDTH=(32 18 16)
#LONWIDTH=(36 30 24)
LATWIDTH=(36 24 18)
LONWIDTH=(40 28 22)

REFINE_RECT_STR=""
REFINE_FACTOR=$((2**REFINE_LEVEL))

# Assuming your CSV file is named data.csv
csv_file="tc_mesh_centers.csv"

# Loop over each line in the CSV file
while IFS=, read -r MESH_ID LONCENTER LATCENTER
do
  echo "ID: $MESH_ID, Longitude: $LONCENTER, Latitude: $LATCENTER"

  for (( i=0; i<REFINE_LEVEL; i++ )); do

      # Use bc for calculations involving decimals
      half_lat_width=$(echo "${LATWIDTH[i]} / 2" | bc -l)
      half_lon_width=$(echo "${LONWIDTH[i]} / 2" | bc -l)

      LATMIN=$(echo "$LATCENTER - $half_lat_width" | bc -l)
      LATMAX=$(echo "$LATCENTER + $half_lat_width" | bc -l)
      LONMIN=$(echo "$LONCENTER - $half_lon_width" | bc -l)
      LONMAX=$(echo "$LONCENTER + $half_lon_width" | bc -l)

      # Format output to limit decimal places, if desired
      LATMIN=$(printf "%.8f" $LATMIN)
      LATMAX=$(printf "%.8f" $LATMAX)
      LONMIN=$(printf "%.8f" $LONMIN)
      LONMAX=$(printf "%.8f" $LONMAX)

      REFINE_RECT_STR=${REFINE_RECT_STR}"${LONMIN},${LATMIN},${LONMAX},${LATMAX},$((i+1));"

  done

  # Remove final trailing semicolon
  REFINE_RECT_STR="${REFINE_RECT_STR%;}"
  echo "$REFINE_RECT_STR"

  GRIDNAME=TClandfall-${MESH_ID}_ne${BASERES}x${REFINE_FACTOR}
  ${SQUADGEN}/SQuadGen --output ${GRIDNAME}.g --refine_rect ${REFINE_RECT_STR} --refine_level ${REFINE_LEVEL} --lat_ref ${LATCENTER} --lon_ref ${LONCENTER} --resolution ${BASERES} --orient_ref -3 --smooth_type SPRING

  GenerateVolumetricMesh --in ${GRIDNAME}.g --out ${GRIDNAME}_pg2.g --np 2 --uniform
  ConvertMeshToSCRIP --in ${GRIDNAME}_pg2.g --out ${GRIDNAME}_pg2_scrip.nc

  ncl gridplot.ncl 'gridfile="'${GRIDNAME}'.g"'

  # Clear REFINE_RECT_STR
  unset REFINE_RECT_STR

  mkdir -p exodus
  mkdir -p scrip

  mv *_scrip.nc scrip
  mv *.g exodus

done < "$csv_file"
