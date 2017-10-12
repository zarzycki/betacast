#!/bin/bash

FILES=/glade/scratch/zarzycki/forecast_natlantic_30_x4_CAM5/run/2017090412/*h0*.nc
for f in $FILES
do
  echo "Processing $f file..."
  ncl weatherplot.NEW.ncl inisec=43200 iniday=04 inimon=09 iniyear=2017 \
  'filename="'${f}'"'
done
