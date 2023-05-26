# Betacast and regional analyses/models

There may be situations where it is beneficial to overlay a regional analysis product over a global analysis. This appears to be particularly useful when running very high-resolution models such as MPAS or SCREAM at O(3km) where mesoscale features are poorly initialized from ERA5 or other traditional reanalysis products. Examples of regional analysis products are RAP, HRRR, or NARR.

Informally known as Frankengrid (patent pending!)

### RAP support

```
YYYYMMDD=20120630
HH=06
YYYYMMDDHH=${YYYYMMDD}${HH}
NLEV=128
RAPFILE=/glade/work/zarzycki/sewx//INIC/RAP_${YYYYMMDDHH}.nc
ERA5FILE=/glade/work/zarzycki/sewx//INIC/ERA5_${YYYYMMDDHH}.nc
RDADIR=/glade/collections/rda/data/ds633.0/

ncl -n atm_to_cam.ncl 'datasource="RAP"' \
  'dycore="se"' \
  numlevels=${NLEV} \
  YYYYMMDDHH=${YYYYMMDDHH} \
  'data_filename = "/glade/u/home/zarzycki/work/2012_derecho/rap_130_'${YYYYMMDD}'_'${HH}'00_000.grib2"' \
  'wgt_filename="/glade/u/home/zarzycki/betacast/remapping/map_rap_13km_TO_rrm.x4_ESMF.nc_patc.nc"' \
  'adjust_config=""' \
  'se_inic = "'${RAPFILE}'"' \
  'model_topo_file="/glade/u/home/zarzycki/work/2012_derecho/USGS-gtopo30_rrm_x4_np4pg2_16xdel2.nc"'

ncl -n atm_to_cam.ncl 'datasource="ERA5RDA"' \
  'dycore="se"' \
  numlevels=${NLEV} \
  YYYYMMDDHH=${YYYYMMDDHH} \
  'RDADIR="'${RDADIR}'"' \
  'data_filename = "/glade/collections/rda/data//ds633.0/e5.oper.invariant/197901/e5.oper.invariant.128_129_z.ll025sc.1979010100_1979010100.nc"' \
  'wgt_filename="/glade/u/home/zarzycki/betacast/remapping/map_era5_0.25x0.25_TO_rrm.x4_ESMF.nc_patc.nc"' \
  'adjust_config=""' \
  'se_inic = "'${ERA5FILE}'"' \
  'model_topo_file="/glade/u/home/zarzycki/work/2012_derecho/USGS-gtopo30_rrm_x4_np4pg2_16xdel2.nc"'

ncl overlay-file.ncl 'base_file = "'${ERA5FILE}'"' 'top_file = "'${RAPFILE}'"'
```

### HWRF support

```
### Settings
YYYYMMDD=20230523
HH=00
YYYYMMDDHH=${YYYYMMDD}${HH}
NLEV=32
HWRF_SRC=/glade/u/home/zarzycki/work/2023_mawar/mawar02w.${YYYYMMDDHH}.hwrfprs.storm.0p015.f000.grb2
CAM_GFSFILE=/glade/work/zarzycki/sewx//INIC/GFS_${YYYYMMDDHH}.nc
CAM_HWRFFILE=/glade/work/zarzycki/sewx//INIC/HWRF_${YYYYMMDDHH}.nc
WGTFILE=./map_hwrf_storm_TO_modelgrid_patc.nc
MODEL_SCRIP=/glade/work/zarzycki/grids/scrip/ne120np4_SCRIP.nc
MODEL_TOPO=/glade/p/cesmdata/cseg/inputdata/atm/cam/topo/se/ne120np4_nc3000_Co015_Fi001_PF_nullRR_Nsw010_20171011.nc

echo "Generating global initial conditions"
ncl -n atm_to_cam.ncl 'datasource="GFS"' \
  'dycore="se"' \
  numlevels=${NLEV} \
  YYYYMMDDHH=${YYYYMMDDHH} \
  'data_filename = "/glade/u/home/zarzycki/work/2023_mawar/gfs_atm_'${YYYYMMDDHH}'.grib2"' \
  'wgt_filename="/glade/u/home/zarzycki/betacast/remapping/map_gfs_0.25x0.25_TO_ne120np4_patc.nc"' \
  'adjust_config=""' \
  'se_inic = "'${CAM_GFSFILE}'"' \
  'model_topo_file="'${MODEL_TOPO}'"'

echo "Generating a temporary SCRIP file for HWRF"
ncl ../remapping/gen_reglatlon_SCRIP.ncl \
  'DSTGRIDNAME="hwrf_storm_scrip.nc"' \
  'DSTDIR="./"' \
  'SRCFILENAME="'${HWRF_SRC}'"'

echo "Generating a temporary map file for HWRF"
ncl ../remapping/gen_analysis_to_model_wgt_file.ncl \
  'ANLGRID="hwrf_storm"' \
  'DSTGRIDNAME="modelgrid"' \
  'ANLGRIDPATH="./"' \
  'WGTFILEDIR="./"' \
  'DSTGRIDFILE="'${MODEL_SCRIP}'"'

echo "Generating regional Frankengrid for HWRF"
ncl -n atm_to_cam.ncl 'datasource="HWRF"' \
  'dycore="se"' \
  numlevels=${NLEV} \
  YYYYMMDDHH=${YYYYMMDDHH} \
  'data_filename = "'${HWRF_SRC}'"' \
  'wgt_filename="'${WGTFILE}'"' \
  'adjust_config=""' \
  'se_inic = "'${CAM_HWRFFILE}'"'

echo "Overlay regional file on top of basefile"
cp -v ${CAM_GFSFILE} ${CAM_GFSFILE}_bak.nc
ncl overlay-file.ncl 'base_file = "'${CAM_GFSFILE}'"' 'top_file = "'${CAM_HWRFFILE}'"'

echo "Cleaning up temporary ESMF files"
rm -v $WGTFILE
rm -v hwrf_storm_scrip.nc
```

