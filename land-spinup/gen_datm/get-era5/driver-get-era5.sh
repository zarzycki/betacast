#!/bin/bash

# NOTE, need cdsapi Python
# On Cheyenne; run ncar_pylib
# conda activate meteo

## declare an array variable
declare -a months=("01" "02" "03" "04" "05" "06" "07" "08" "09" "10" "11" "12")

STYR=1988
ENYR=2023
OUTDIR=/glade/derecho/scratch/zarzycki/ERA5-DATM/
#OUTDIR="/global/homes/c/czarzyck/scratch/ERA5-DATM/"

mkdir -p -v $OUTDIR

## now loop through the above array
for jj in $(seq $STYR $ENYR);
do
  for ii in "${months[@]}"
  do
    echo $jj' '$ii
    OUTFILE=${OUTDIR}/out.${jj}.${ii}.nc
    if [ -f "$OUTFILE" ]; then
      echo "$OUTFILE exists, not downloading"
    else
      echo "$OUTFILE doesn't exist..."
      python get-data.py "${jj}" "${ii}" "${OUTDIR}"
      if command -v ncks >/dev/null 2>&1; then
        echo "Compressing file using ncks..."
        ncks -O -4 -L 1 "$OUTFILE" "$OUTFILE"
      else
        echo "ncks is not in PATH; skipping compression."
      fi
    fi
  done
done
