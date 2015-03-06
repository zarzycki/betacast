



if (datasource .eq. "GFS") then
end if
if (datasource .eq. "ERAI") then
end if


ncl -n gfs.ncl 'datasource="GFS"'     \
       numlevels=30           \
       YYYYMMDDHH=2014080712  \
       'gridname = "uniform_60"' \
       'data_filename = "/glade/p/work/zarzycki/sewx/GFS/gfs_atm_2014080718.grib2"'  \
       'wgt_filename="/glade/p/work/zarzycki/maps/gfsmaps/map_gfs0.50_TO_uniform_60_patc.nc"' \
       'se_inic = "/glade/p/work/zarzycki/sewx/INIC/uniform_60_INIC.nc"'
       

    ncl -n gfs.ncl 'datasource="ERAI"'     \
       numlevels=30           \
       YYYYMMDDHH=2013102300  \
       'gridname = "tcforecast_60_x4"' \
       'data_filename = "/glade/scratch/zarzycki/getECMWFdata/ERA-Int_2013102300.nc"'  \
       'wgt_filename="/glade/u/home/zarzycki/scratch/getECMWFdata/map_ERA_TO_tcforecast_60_x4_patc.140530.nc"' \
       'se_inic = "/glade/p/work/zarzycki/sewx/INIC/tcforecast_60_x4_INIC.nc"'

    ncl -n gfs.ncl 'datasource="CFSR"'     \
       numlevels=30           \
       YYYYMMDDHH=2005082812  \
       'gridname = "uniform_60"' \
       'data_filename = "/glade/scratch/zarzycki/pgbh06.gdas.2005082812.grb2"'  \
       'wgt_filename="/glade/p/work/zarzycki/maps/gfsmaps/map_gfs0.50_TO_uniform_60_patc.nc"' \
       'se_inic = "/glade/p/work/zarzycki/sewx/INIC/uniform_60_INIC.nc"'