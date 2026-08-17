################################################################
#### Casper
################################################################
#PBS -N betacast
#PBS -A P93300042
#PBS -l select=1:ncpus=16:mem=240GB
#PBS -l walltime=23:59:00
#PBS -q casper@casper-pbs
#PBS -j oe
################################################################

module load peak-memusage
module load conda
conda activate betacast

date

#python write.py

TEST_FILES_DIR=/glade/u/home/zarzycki/scratch/test_files/
DEBUG_FILE_DIR=/glade/u/home/zarzycki/scratch/tmp_betacast/
export BETACAST=/glade/u/home/zarzycki/scratch/betacast/

# peak_memusage python atm_to_cam.py \
#   --datasource "ERA5RDA" \
#   --numlevels 32 \
#   --YYYYMMDDHH 2005082800 \
#   --data_filename "${TEST_FILES_DIR}/ds633.0/e5.oper.invariant/197901/e5.oper.invariant.128_129_z.ll025sc.1979010100_1979010100.nc" \
#   --wgt_filename "${TEST_FILES_DIR}/map_gfs_0.25x0.25_TO_mpasa120_patc.nc" \
#   --dycore "mpas" \
#   --RDADIR "${TEST_FILES_DIR}/ds633.0/" \
#   --compress_file \
#   --write_floats \
#   --add_cloud_vars \
#   --adjust_config '""' \
#   --model_topo_file "${TEST_FILES_DIR}/mpasa120.CFSR.L32.nc" \
#   --se_inic "${DEBUG_FILE_DIR}/py_final.nc" \
#   --verbose \
#   --write_debug_files \
#   --write_debug_dir "${DEBUG_FILE_DIR}" \
#   --verbose

peak_memusage python atm_to_cam.py \
  --datasource "ERA5mlRDA" \
  --numlevels 93 \
  --YYYYMMDDHH 2000010100 \
  --data_filename "/glade/campaign/collections/rda/data//ds633.6/e5.oper.invariant/e5.oper.invariant.128_129_z.regn320sc.2016010100_2016010100.nc" \
  --wgt_filename "/glade/campaign/cgd/amp/aherring/betacast-pivot/py_atm_to_cam/ERA5ml_TO_mpasa3.75-bilin.nc" \
  --dycore "mpas" \
  --RDADIR "/glade/campaign/collections/rda/data/ds633.6/" \
  --compress_file \
  --write_floats \
  --adjust_config '""' \
  --model_topo_file "/lustre/desc1/scratch/bdobbins/ic_io_err_3.75km/fmthist_3.75km/fmthist_3.75km.init.nc" \
  --se_inic "/glade/derecho/scratch/zarzycki/betacast_oob_fhist_mpasa3.75km_c250110_v2.nc" \
  --verbose

date
