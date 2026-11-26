#!/bin/bash

#MAPFILE=~/m2637/betacast/regrid_maps/map_ne0np4natlanticref.ne30x4_TO_era5_0.25x0.25_patc.nc
MAPFILE=$3

#SOURCE=/pscratch/sd/c/czarzyck/hyperion/CHEY.VR28.NATL.REF.CAM5.4CLM5.0.dtime900/h5/CHEY.VR28.NATL.REF.CAM5.4CLM5.0.dtime900.cam.h5.1989-11-30-00000.nc
SOURCE=$1

#OUTDIRBASE="/pscratch/sd/c/czarzyck/hyperion/DATM/"
OUTDIRBASE=$2

# true: quantize + deflate to netCDF4 (1 time step per chunk); false: convert to nc3 (64-bit offset)
LOSSY_COMPRESS=${4:-true}

SOURCEDIR="${SOURCE%/*}/"
SOURCEFILE="${SOURCE##*/}"
SOURCEDATE="${SOURCEFILE##*.h?.}"
SOURCEDATE="${SOURCEDATE%.nc}"
SOURCECASE="${SOURCEFILE%%.?am.h*}"

echo "SOURCEDIR: $SOURCEDIR"
echo "SOURCEFILE: $SOURCEFILE"
echo "SOURCEDATE: $SOURCEDATE"
echo "SOURCECASE: $SOURCECASE"

OUTDIR=${OUTDIRBASE}/${SOURCECASE}

mkdir -pv $OUTDIR
mkdir -pv $OUTDIR/Solar/
mkdir -pv $OUTDIR/Precip/
mkdir -pv $OUTDIR/TPQW/

TMPsolar=$OUTDIR/Solar/${SOURCEDATE}_solar.nc
TMPprecip=$OUTDIR/Precip/${SOURCEDATE}_precip.nc
TMPstate=$OUTDIR/TPQW/${SOURCEDATE}_state.nc

FINALsolar=$OUTDIR/Solar/${SOURCECASE}_${SOURCEDATE}_Solar.nc
FINALprecip=$OUTDIR/Precip/${SOURCECASE}_${SOURCEDATE}_Precip.nc
FINALstate=$OUTDIR/TPQW/${SOURCECASE}_${SOURCEDATE}_TPQW.nc

if [[ -f "$FINALsolar" && -f "$FINALprecip" && -f "$FINALstate" ]]; then
  echo "All files for $SOURCEDATE already exist. Exiting script."
  exit 1
fi

echo "Remapping $SOURCEDATE"
echo "Writing ${TMPsolar}  ${TMPprecip}  ${TMPstate}"
ncremap -v FSDS -i ${SOURCE} -o ${TMPsolar} -m ${MAPFILE}
ncremap -v PRECT -i ${SOURCE} -o ${TMPprecip} -m ${MAPFILE}
ncremap -v FLDS,UBOT,VBOT,QBOT,TBOT,PS -i ${SOURCE} -o ${TMPstate} -m ${MAPFILE}

echo "Inverting latitude dimension because of annoying ERA5 N->S"
ncpdq -O -a '-lat' ${TMPsolar} ${TMPsolar}
ncpdq -O -a '-lat' ${TMPprecip} ${TMPprecip}
ncpdq -O -a '-lat' ${TMPstate} ${TMPstate}

echo "Creating ZBOT and WIND"
ncap2 -O -s 'ZBOT=UBOT*float(0.0)+float(50.0)' -s 'WIND=sqrt(pow(UBOT,2)+pow(VBOT,2))' ${TMPstate} ${TMPstate}

echo "Updating attributes"
ncatted -O -a long_name,ZBOT,o,c,"Lowest model level height" ${TMPstate}
ncatted -O -a units,ZBOT,o,c,"m" ${TMPstate}
ncatted -O -a long_name,WIND,o,c,"Lowest model level wind speed" ${TMPstate}

echo "Updating units"
ncap2 -O -s 'PRECT=PRECT*float(1000.0)' ${TMPprecip} ${TMPprecip}
ncatted -O -a units,PRECT,o,c,"mm/s" ${TMPprecip}

echo "Renaming variables to match ERA5 streams"
ncrename -v PS,PSRF ${TMPstate}
ncrename -v PRECT,PRECTmms ${TMPprecip}

echo "Deleting variables"
ncks -O -C -x -v area,gw ${TMPsolar} ${TMPsolar}
ncks -O -C -x -v UBOT,VBOT,area,gw ${TMPstate} ${TMPstate}
ncks -O -C -x -v area,gw ${TMPprecip} ${TMPprecip}

if [ "$LOSSY_COMPRESS" = true ]; then
  # Quantize (VAR=.n -> nearest 10^-n, VAR=n -> n sig. digits), then DEFLATE level 1.
  # Chunk as (1 time, full lat, full lon) since DATM reads one time slice at a time.
  echo "Lossy compressing to nc4"
  ncks -O -4 -L 1 --cnk_plc=g2d --cnk_map=rd1 --qnt FSDS=.1 ${TMPsolar} ${FINALsolar}
  ncks -O -4 -L 1 --cnk_plc=g2d --cnk_map=rd1 --qnt PRECTmms=3 ${TMPprecip} ${FINALprecip}
  ncks -O -4 -L 1 --cnk_plc=g2d --cnk_map=rd1 --qnt 'TBOT=.2#PSRF=.0#QBOT=4#WIND=.2#FLDS=.1' ${TMPstate} ${FINALstate}
else
  echo "Converting to nc3"
  nccopy -k nc6 ${TMPsolar} ${FINALsolar}
  nccopy -k nc6 ${TMPprecip} ${FINALprecip}
  nccopy -k nc6 ${TMPstate} ${FINALstate}
fi

echo "Deleting bad files"
rm -v ${TMPstate} ${TMPsolar} ${TMPprecip}

echo "all done"
