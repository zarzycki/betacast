#!/bin/bash

#PBS -N CESMLENS-delta
#PBS -A UNSB0017
#PBS -l select=1:ncpus=4:mpiprocs=4:mem=60GB
#PBS -l walltime=24:00:00
#PBS -q casper

module load nco
module load cdo
module load parallel

which ncremap

VAR="SST"
PINT=false #if a variable is 3D and converstion to pressure levels is needed we need PS first
TEMP_DIR=/glade/scratch/zarzycki/CESM_LENS_temp/${VAR}
if [ $VAR = "SST" ] ; then
  # Go to ocn dir
  LENS_DIR=/glade/collections/cdg/data/cesmLE/CESM-CAM5-BGC-LE/ocn/proc/tseries/monthly/${VAR}
else
  LENS_DIR=/glade/collections/cdg/data/cesmLE/CESM-CAM5-BGC-LE/atm/proc/tseries/monthly/${VAR}
fi

#HYPERSLABSTRING="-d lev,29,29"
HYPERSLABSTRING=""

# use numcores = 32 for 1deg runs, 16 for high-res runs
NUMCORES=4
TIMESTAMP=`date +%s%N`
COMMANDFILE=commands.${TIMESTAMP}.txt

#=================================================
#Calculate Mean Annual Cycle Over Control Run
#=================================================
mkdir -p ${TEMP_DIR}

rm ${TEMP_DIR}/1850_all_files.nc

if [ $VAR = 'SST' ] ; then
  ncrcat -v ${VAR} ${HYPERSLABSTRING} ${LENS_DIR}/b.e11.B20TRC5CNBDRD.f09_g16.005.pop.h.*.nc ${TEMP_DIR}/1850_all_files.nc
  #/glade/u/apps/ch/opt/nco/4.9.5/gnu/9.1.0/bin/ncremap -v ${VAR} --alg_typ=neareststod --esmf_typ=neareststod -i ${TEMP_DIR}/1850_all_files.nc -d /glade/scratch/zarzycki/CESM_LENS_temp/T/ens_T_anom.nc -o ${TEMP_DIR}/1850_all_files_remap.nc
  #mv ${TEMP_DIR}/1850_all_files_remap.nc ${TEMP_DIR}/1850_all_files.nc
else
  ncrcat -v ${VAR} ${HYPERSLABSTRING} ${LENS_DIR}/b.e11.B1850C5CN.f09_g16.005.cam.h0.*.nc ${TEMP_DIR}/1850_all_files.nc
fi

if ${PINT} ; then
  echo "Converting 3D ${VAR} to pressure levels for 1850"
  ncks -A -v PS /glade/scratch/zarzycki/CESM_LENS_temp/PS/1850_all_files.nc ${TEMP_DIR}/1850_all_files.nc
  ncatted -a bounds,lev,c,c,"ilev"  ${TEMP_DIR}/1850_all_files.nc
  cdo ml2pl,100000,92500,85000,70000,60000,50000,40000,30000,25000,20000,15000,10000,7000,5000,3000,2000,1000 ${TEMP_DIR}/1850_all_files.nc ${TEMP_DIR}/1850_all_files_pres.nc 
  rm ${TEMP_DIR}/1850_all_files.nc
  mv ${TEMP_DIR}/1850_all_files_pres.nc ${TEMP_DIR}/1850_all_files.nc
fi

if [ $VAR != 'SST' ] ; then
  ncrename -a .missing_value,_FillValue ${TEMP_DIR}/1850_all_files.nc
fi

rm ${TEMP_DIR}/1850_climo.nc
for month in {01..12}
  do
  echo "ncra -F -d time,$month,,12 ${TEMP_DIR}/1850_all_files.nc ${TEMP_DIR}/temp$month.nc" >> ${COMMANDFILE}
done
parallel --jobs ${NUMCORES} --workdir $PWD < ${COMMANDFILE}
rm ${COMMANDFILE}

ncrcat -h ${TEMP_DIR}/temp*nc ${TEMP_DIR}/1850_climo.nc
rm ${TEMP_DIR}/temp*nc

#=================================================
#Calculate Ensemble Mean for 1920 thru 2100
#=================================================

for n in {001..105} #SET number of ensembles
   do
     echo "echo $n; \
     rm ${TEMP_DIR}/ens_${n}.nc; \
     ls ${LENS_DIR}/*D.f09_g16.${n}.*.nc; \
     ncrcat -v ${VAR} ${HYPERSLABSTRING} -d time,'1920-01-01 0:00:0.1','2101-01-01 00:00:0.0' ${LENS_DIR}/*D.f09_g16.${n}.*.nc ${TEMP_DIR}/ens_${n}.nc; \
     if ${PINT}; \
     then echo 'Converting 3D ${VAR} to pressure levels for ensemble # ${n}'; \
       ncks -A -v PS /glade/scratch/zarzycki/CESM_LENS_temp/PS/ens_${n}.nc ${TEMP_DIR}/ens_${n}.nc; \
       ncatted -a bounds,lev,c,c,'ilev'  ${TEMP_DIR}/ens_${n}.nc; \
       cdo ml2pl,100000,92500,85000,70000,60000,50000,40000,30000,25000,20000,15000,10000,7000,5000,3000,2000,1000 ${TEMP_DIR}/ens_${n}.nc ${TEMP_DIR}/ens_${n}_pres.nc; \
       rm ${TEMP_DIR}/ens_${n}.nc; \
       mv ${TEMP_DIR}/ens_${n}_pres.nc ${TEMP_DIR}/ens_${n}.nc; \
     fi; \
     if [ $VAR != 'SST' ] ; \
       then ncrename -a .missing_value,_FillValue ${TEMP_DIR}/ens_${n}.nc;\
     fi" >> ${COMMANDFILE}
done

#       if [ $VAR = 'SST' ] ; \
#         then /glade/u/apps/ch/opt/nco/4.9.5/gnu/9.1.0/bin/ncremap -v ${VAR} --alg_typ=neareststod --esmf_typ=neareststod -i ${TEMP_DIR}/ens_${n}.nc -d /glade/scratch/zarzycki/CESM_LENS_temp/T/ens_T_anom.nc -o ${TEMP_DIR}/ens_${n}_remap.nc; \
#         mv ${TEMP_DIR}/ens_${n}_remap.nc ${TEMP_DIR}/ens_${n}.nc; \
#       fi; \

parallel --jobs ${NUMCORES} --workdir $PWD < ${COMMANDFILE}
rm ${COMMANDFILE}

rm ${TEMP_DIR}/ens_mean.nc
nces ${TEMP_DIR}/ens_*.nc ${TEMP_DIR}/ens_mean.nc

#=================================================
#Remove 1850 Annual Cycle from Ensemble Mean
#=================================================
rm ${TEMP_DIR}/ens_anom.nc
cdo sub ${TEMP_DIR}/ens_mean.nc ${TEMP_DIR}/1850_climo.nc ${TEMP_DIR}/ens_anom.nc
if ${PINT} ; then
  ncpdq -a '-plev' ${TEMP_DIR}/ens_anom.nc ${TEMP_DIR}/ens_${VAR}_anom.nc
else
  cp ${TEMP_DIR}/ens_anom.nc ${TEMP_DIR}/ens_${VAR}_anom.nc
fi
ncks -O -C -x -v lev_bnds ${TEMP_DIR}/ens_${VAR}_anom.nc ${TEMP_DIR}/ens_${VAR}_anom.nc
rm ${TEMP_DIR}/ens_anom.nc

#ncwa -a z_t ens_SST_anom.nc tmp.nc
#ncremap -v SST --alg_typ=neareststod --esmf_typ=neareststod -d /glade/scratch/zarzycki/CESM_LENS_temp/T/ens_T_anom.nc -i tmp.nc -o tmp2.nc

