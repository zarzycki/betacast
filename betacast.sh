#!/bin/bash -l

#SBATCH -C cpu
#SBATCH -A m2637
#SBATCH --qos=premium
#SBATCH --time=02:20:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=0

module load conda
conda activate betacast

###################################################################################
# Colin Zarzycki (czarzycki@psu.edu)
#
# Driver script for running CESM/CAM in "forecast" or "hindcast" mode.
# This code will either read UTC unix clock or read in a specified list of dates,
# download dataset, map to CAM grid, and run a forecast.
#
# Generally can be executed on login nodes or in the background, ex:
# $> nohup ./main-realtime.sh nl.neconus30x8 &
# but see above for example of PBS options to submit to batch nodes
# $> sbatch betacast.sh machine_files/machine.pm-cpu namelists/nl.philly.pm-gpu output_streams/SCREAM.packed.yaml
#
# Details can be found in:
# C. M. Zarzycki and C. Jablonowski (2015), Experimental tropical cyclone forecasts
# using a variable-resolution global model. Mon. Weat. Rev., 143(10), 4012–4037.
# doi:10.1175/MWR-D-15-0159.1.
###################################################################################

echo "Command used: $0 \"$@\""
echo "PID       : $$"
echo "Host      : $(hostname)"
echo "User      : $(whoami)"
echo "Start time: $(date -u +"%Y-%m-%d %H:%M:%S UTC")"
echo "Shell     : $SHELL"

script_start=$(date +%s)
set -e
#set -v

if [[ -n "$SLURM_SUBMIT_DIR" ]]; then
  # Running under Slurm batch
  SCRIPTPATH="$SLURM_SUBMIT_DIR"
  USING_BATCH=true
  BATCH_PREFIX="sbatch"
  AVAIL_CPUS=${SLURM_CPUS_ON_NODE:-${SLURM_CPUS_PER_TASK:-1}}
  AVAIL_MEM_MB=${SLURM_MEM_PER_NODE:-0}  # In MB
  SUBMIT_EPOCH=$(scontrol show job $SLURM_JOB_ID | grep SubmitTime | sed 's/.*SubmitTime=\([^ ]*\).*/\1/' | xargs -I {} date -d {} +%s)
  QUEUE_MINUTES=$(( ($(date +%s) - SUBMIT_EPOCH) / 60 ))
  echo "Job $SLURM_JOB_ID currently running, was queued for $QUEUE_MINUTES minutes"
elif [[ -n "$PBS_O_WORKDIR" ]]; then
  # Running under PBS/Torque
  SCRIPTPATH="$PBS_O_WORKDIR"
  USING_BATCH=true
  BATCH_PREFIX="qsub"
  AVAIL_CPUS=${PBS_NUM_PPN:-$(nproc)}
  AVAIL_MEM_MB=$(free -m | awk '/^Mem:/ {print $2}')
  SUBMIT_EPOCH=$(qstat -f $PBS_JOBID | grep 'qtime = ' | sed 's/.*= //')
  QUEUE_MINUTES=$(( ($(date +%s) - SUBMIT_EPOCH) / 60 ))
  echo "Job $PBS_JOBID currently running, was queued for $QUEUE_MINUTES minutes"
else
  # Running interactively, nohup, by hand
  SCRIPTPATH=$(dirname "$(realpath "$0")")
  USING_BATCH=false
  BATCH_PREFIX=""
  AVAIL_CPUS=$(nproc)
  AVAIL_MEM_MB=$(free -m | awk '/^Mem:/ {print $2}')
fi
echo "Available CPUs: $AVAIL_CPUS"
echo "Available Memory: $AVAIL_MEM_MB MB"
if [[ "$USING_BATCH" == true ]] && command -v parallel &>/dev/null; then
  echo "GNU parallel is available"
  GNUPAR_AVAIL=true
else
  GNUPAR_AVAIL=false
fi

echo "Our script path is $SCRIPTPATH"
export BETACAST=${SCRIPTPATH}
echo "BETACAST exported as $BETACAST"
source ${BETACAST}/utils.sh
source ${BETACAST}/datahelpers.sh

# Set files
MACHINEFILE=${1}
NAMELISTFILE=${2}
OUTPUTSTREAMS=${3}
# If relative path, convert to absolute path
if [[ "$MACHINEFILE" != /* ]] && [[ "$MACHINEFILE" != ~* ]]; then MACHINEFILE=${PWD}/${MACHINEFILE}; fi
if [[ "$NAMELISTFILE" != /* ]] && [[ "$NAMELISTFILE" != ~* ]]; then NAMELISTFILE=${PWD}/${NAMELISTFILE}; fi
if [[ "$OUTPUTSTREAMS" != /* ]] && [[ "$OUTPUTSTREAMS" != ~* ]]; then OUTPUTSTREAMS=${PWD}/${OUTPUTSTREAMS}; fi
# If files don't exist, exit now
exit_file_no_exist "$MACHINEFILE"
exit_file_no_exist "$NAMELISTFILE"
exit_file_no_exist "$OUTPUTSTREAMS"
echo "$MACHINEFILE"; echo "$NAMELISTFILE"; echo "$OUTPUTSTREAMS"

# Read namelists
read_bash_nl "${MACHINEFILE}"
read_bash_nl "${NAMELISTFILE}"

# Replace Betacast-specific placeholders. This allows a user to override these in NAMELISTFILE
path_to_case=$(replace_betacast_string "$path_to_case" "$casename" "!CASENAME!")
path_to_rundir=$(replace_betacast_string "$path_to_rundir" "$casename" "!CASENAME!")
sePreFilterIC=$(replace_betacast_string "$sePreFilterIC" "$casename" "!CASENAME!")
sePostFilterIC=$(replace_betacast_string "$sePostFilterIC" "$casename" "!CASENAME!")
sstFileIC=$(replace_betacast_string "$sstFileIC" "$casename" "!CASENAME!")

set -u  # turn on crashes for unbound variables in bash

###################################################################################
############### OPTIONAL TO BE SET BY USER ########################################
path_to_nc_files=${path_to_rundir}              # Path where .nc files are
outputdir=${path_to_rundir}                     # Path where .nc files are being written
tmparchivecdir=${path_to_rundir}/proc/              # Path to temporarily stage final data
landdir=${path_to_rundir}/landstart/            # Path to store land restart files
###################################################################################
### THESE COME WITH THE REPO, DO NOT CHANGE #######################################
atm_to_cam_path=${BETACAST}/atm_to_cam
sst_to_cam_path=${BETACAST}/sst_to_cam
filter_path=${atm_to_cam_path}/filter
###################################################################################

### setting variables not included in namelist for backwards compat
if [ -z "${CIMEbatchargs+x}" ]; then CIMEbatchargs=""; fi
if [ -z "${do_runoff+x}" ]; then do_runoff=false; fi
if [ -z "${keep_land_restarts+x}" ]; then keep_land_restarts=true; fi
if [ -z "${perturb_namelist+x}" ]; then perturb_namelist=""; fi
if [ -z "${predict_docn+x}" ]; then predict_docn=false; fi
if [ -z "${archive_inic+x}" ]; then archive_inic=false; fi
if [ -z "${compress_history_nc+x}" ]; then compress_history_nc=true; fi
if [ -z "${standalone_vortex+x}" ]; then standalone_vortex=false; fi
if [ -z "${vortex_namelist+x}" ]; then vortex_namelist=""; fi
if [ -z "${augment_tcs+x}" ]; then augment_tcs=false; fi
if [ -z "${save_nudging_files+x}" ]; then save_nudging_files=false; fi
if [ -z "${override_rest_check+x}" ]; then override_rest_check=false; fi
if [ -z "${tararchivedir+x}" ]; then tararchivedir=true; fi
if [ -z "${docnres+x}" ]; then docnres="180x360"; fi
if [ -z "${modelgridfile+x}" ]; then modelgridfile=""; fi
if [ -z "${anl2mdlWeights+x}" ]; then anl2mdlWeights=""; fi
if [ -z "${CIMEMAXTRIES+x}" ]; then CIMEMAXTRIES=1; fi
if [ -z "${add_noise+x}" ]; then add_noise=false; fi
if [ -z "${runmodel+x}" ]; then runmodel=true; fi
if [ -z "${DO_PYTHON+x}" ]; then DO_PYTHON=true; fi
### Some defaults infrequently set
if [ -z "${debug+x}" ]; then debug=false; fi
if [ -z "${islive+x}" ]; then islive=false; fi
if [ -z "${datestemplate+x}" ]; then datestemplate=""; fi
if [ -z "${doFilter+x}" ]; then doFilter=false; fi
if [ -z "${filterOnly+x}" ]; then filterOnly=false; fi
if [ -z "${numHoursSEStart+x}" ]; then numHoursSEStart=3; fi
if [ -z "${filterHourLength+x}" ]; then filterHourLength=6; fi
if [ -z "${filtTcut+x}" ]; then filtTcut=6; fi
if [ -z "${add_perturbs+x}" ]; then add_perturbs=false; fi
if [ -z "${land_spinup+x}" ]; then land_spinup=false; fi
if [ -z "${FILTERWALLCLOCK+x}" ]; then FILTERWALLCLOCK="00:29:00"; fi
if [ -z "${FILTERQUEUE+x}" ]; then FILTERQUEUE="batch"; fi
if [ -z "${RUNWALLCLOCK+x}" ]; then RUNWALLCLOCK="12:00:00"; fi
if [ -z "${RUNQUEUE+x}" ]; then RUNQUEUE="regular"; fi
if [ -z "${RUNPRIORITY+x}" ]; then RUNPRIORITY=""; fi
if [ -z "${BATCH_SCHEDULER+x}" ]; then BATCH_SCHEDULER=""; fi
if [ -z "${use_nsplit+x}" ]; then use_nsplit="true"; fi
if [ -z "${cime_coupler+x}" ]; then cime_coupler="mct"; fi
if [ -z "${nclPlotWeights+x}" ]; then nclPlotWeights="NULL"; fi
if [ -z "${landrawdir+x}" ]; then landrawdir="NULL"; fi
if [ -z "${sendplots+x}" ]; then sendplots=false; fi
if [ -z "${dotracking+x}" ]; then dotracking=false; fi
if [ -z "${m2m_gridfile+x}" ]; then m2m_gridfile=""; fi
if [ -z "${m2m_remap_file+x}" ]; then m2m_remap_file=""; fi
### Runner code
if [ -z "${batch_sst_gen+x}" ]; then batch_sst_gen=false; fi
if [ -z "${batch_atm_gen+x}" ]; then batch_atm_gen=false; fi

echo "DO_PYTHON set to $DO_PYTHON"
if [ "$DO_PYTHON" = true ]; then
  # Update default paths to point to python folders
  atm_to_cam_path=${BETACAST}/py_atm_to_cam
  sst_to_cam_path=${BETACAST}/py_sst_to_cam
  # Update paths that relied on atm_to_cam_path or sst_to_cam_path
  filter_path=${atm_to_cam_path}/filter
fi

### Set correct model split
if [ -z "${modelSystem+x}" ]; then modelSystem=0; fi
if [ "$modelSystem" -eq 0 ]; then
  echo "Using CESM"
  atmName="cam"
  lndName="clm"
  lndSpecialName="clm2"
  rofName="mosart"
  rofSpecialName="mosart"
  # As of late 2023, _rtm was dropped in MOSART inic string
  # See: https://github.com/ESCOMP/MOSART/commit/cdd878de4cff9dc562093978e44aeb7122237983
  # <CESM3 dev tags need to set this to finidat_rtm
  rof_finidat="finidat"
elif [ "$modelSystem" -eq 1 ]; then
  echo "Using E3SM"
  atmName="eam"
  lndName="elm"
  lndSpecialName="elm"
  rofName="mosart"
  rofSpecialName="mosart"
  rof_finidat="finidat_rtm"
elif [ "$modelSystem" -eq 2 ]; then
  echo "Using SCREAMv1"
  atmName="scream"
  lndName="elm"
  lndSpecialName="elm"
  rofName="mosart"
  rofSpecialName="mosart"
  rof_finidat="finidat_rtm"
else
  echo "Unknown modeling system set for modelSystem: $modelSystem"
  exit 1
fi

# Are we running with frankengrid?
do_frankengrid=false
if [ -n "${regional_src+x}" ]; then do_frankengrid=true ; fi
echo "do_frankengrid set to: $do_frankengrid"

# Check if ARCHIVEDIR is unset or empty or only whitespace
if [ -z "${ARCHIVEDIR+x}" ] || [ -z "${ARCHIVEDIR//[[:space:]]/}" ]; then
  ARCHIVEDIR=${outputdir}/
else
  ARCHIVEDIR=${ARCHIVEDIR}/${casename}/
  mkdir -v -p ${ARCHIVEDIR}
fi
echo "Files will be archived in ${ARCHIVEDIR}/YYYYMMDDHH/"

# Update SST filename based on docnres
sstFileIC=${sstFileIC/DOCNRES/$docnres}
echo "Actual sstFileIC: ${sstFileIC}"

# Get times!
process_model_times "$islive" "$casename" "$datestemplate" "$numHoursSEStart" "$doFilter"
# Now we have globally set...
#   - yearstr, monthstr, daystr, cyclestr, cyclestrsec
#   - se_yearstr, se_monthstr, se_daystr, se_cyclestr, se_cyclestrsec
#   - sstyearstr, sstmonthstr, sstdaystr, sstcyclestr (via getSSTtime)
#   - yestmonthstr, yestdaystr, yestyearstr
#   - datesfile (for non-live runs)
#   - currtime_UTC_HHMM (current UTC time in HHMM format)
#   ${yearstr}${monthstr}${daystr}${cyclestr}

# Replace any time strings
if [[ -n "$vortex_namelist" ]]; then
  vortex_namelist=$(replace_betacast_string "$vortex_namelist" "${yearstr}${monthstr}${daystr}${cyclestr}" "!YYYYMMDDHH!")
fi

### ERROR CHECKING BLOCK! #########################################################

# Adjust bools (for backwards compatibility, 0 = false and 1 = true)
bools_to_check=("islive" "debug" "doFilter" "filterOnly"
    "do_runoff" "keep_land_restarts" "predict_docn" "archive_inic"
    "compress_history_nc" "override_rest_check" "batch_atm_gen" "batch_sst_gen"
    "standalone_vortex" "augment_tcs" "tararchivedir" "save_nudging_files"
    "add_noise" "runmodel" "DO_PYTHON" "add_perturbs" "land_spinup"
    "sendplots" "dotracking")
for bool_to_check in ${bools_to_check[@]}; do
  check_bool "$bool_to_check" ${!bool_to_check}
done

# Exit if add_perturbs is turned on but no namelist is passed in with perturbation config settings
if [ "$add_perturbs" = true ] && { [ -z "$perturb_namelist" ] || [ ! -f "$perturb_namelist" ]; } ; then
  echo "add_perturbs is true but can't find namelist: $perturb_namelist"
  exit 1
fi

# Exit vortex_namelist is defined, is not NULL, and doesn't exist
if [ -n "$vortex_namelist" ] && [ "$vortex_namelist" != "NULL" ] && [ ! -f "$vortex_namelist" ] ; then
  echo "Erroneous vortex namelist: $vortex_namelist"
  exit 1
fi

# Check for deprecated namelist options
if [ -z "${anl2mdlWeights+x}" ] && [ -n "${gfs2seWeights+x}" ] ; then
  echo "WARNING: Setting anl2mdlWeights to ${gfs2seWeights}"
  echo "WARNING: This is deprecated and will be removed in the future! To fix, change 'gfs2seWeights' to 'anl2mdlWeights' in ${NAMELISTFILE}"
  anl2mdlWeights="$gfs2seWeights"
fi

if [ "$save_nudging_files" = true ] && [ "$modelSystem" -eq 2 ] ; then
  echo "ERROR: save_nudging_files cannot be true with SCREAM model (2) set"
  exit 1
fi

# Check if ncl exists
if ! type ncl &> /dev/null ; then
  echo "ERROR: ncl does not exist. Make sure ncl is in your path when betacast is invoked"
  exit 1
fi

# Check if ncks exists for compression
ncks_exists=true
if ! type ncks &> /dev/null ; then
  # SCREAM requires NCO to convert filetypes, so this is a hard stop
  if [ "$modelSystem" -eq 2 ]; then
    echo "ERROR: ncks does not exist. Make sure ncks is in your path when betacast is invoked"
    exit 1
  fi
  echo "WARNING: ncks does not exist, cannot compress. Setting compress_history_nc to 0 (false)"
  echo "WARNING: if you'd like remedy, make sure ncks is in your path when betacast.sh is invoked"
  compress_history_nc=false
  ncks_exists=false
fi

# If we are using nuopc, we need to generate an ESMF mesh
if [ "${cime_coupler}" == "nuopc" ] && ! type ESMF_Scrip2Unstruct &> /dev/null; then
  echo "ERROR: ESMF_Scrip2Unstruct does not exist, which is needed when cime_coupler is: ${cime_coupler}"
  echo "ERROR: Install via conda/mamba or from source and add to PATH"
  exit 1
fi

# Check to make sure required ATM variables are provided by the user if we are using model-to-model
if [[ "$atmDataType" -eq 9 ]]; then
  if [[ -z "$m2m_topo_in" ]] || [[ -z "$m2m_parent_source" ]] || { [[ -z "$m2m_remap_file" ]] && [[ -z "$m2m_gridfile" ]]; }; then
    echo "ERROR: For m2m ATM, one or more required variables (m2m_topo_in, m2m_parent_source, m2m_remap_file/m2m_gridfile pair) are not defined."
    echo "ERROR: Ensure these are defined in the Betacast namelist"
    exit 1
  fi
fi

# Check to make sure required SST variables are provided by the user if we are using model-to-model
if [[ "$sstDataType" -eq 9 ]]; then
  if [[ -z "$m2m_sst_grid_filename" ]] || [[ -z "$m2m_sstice_data_filename" ]] || [[ -z "$m2m_sstice_year_align" ]] || [[ -z "$m2m_sstice_year_start" ]] || [[ -z "$m2m_sstice_year_end" ]]; then
    echo "ERROR: For m2m SST, one or more required variables (m2m_sst_grid_filename, m2m_sstice_data_filename, m2m_sstice_year_align, m2m_sstice_year_start, m2m_sstice_year_end) are not defined."
    echo "ERROR: Ensure these are defined in the Betacast namelist"
    exit 1
  fi
fi

# If batch mode is requested, scheduler must be defined
if { [ "${batch_atm_gen}" = true ] || [ "${batch_sst_gen}" = true ]; } && [ -z "${BATCH_SCHEDULER:-}" ]; then
  echo "ERROR: BATCH_SCHEDULER is not set, but batch_atm_gen or batch_sst_gen is true."
  echo "Please export BATCH_SCHEDULER in the machine file or turn off batch data generation."
  exit 1
fi

# Check Python version for 3+
python -c 'import sys; exit(sys.version_info.major != 3)' && echo "CHECK_PYTHON: Python 3 found" || { echo "CHECK_PYTHON: Please install Python 3+"; exit 25; }

# Check Python package dependencies
if [ "$do_frankengrid" = true ] ; then
  check_python_dependency numpy
  check_python_dependency netCDF4
  check_python_dependency sklearn
fi

if [ "$override_rest_check" = false ]; then
  echo "Checking for SourceMods permiting additional restart writes for land model"
  echo "This check can be ignored with override_rest_check = true in the namelist."
  exit_files_no_exist "$path_to_case/SourceMods/src.${lndName}/lnd_comp_mct.F90" "$path_to_case/SourceMods/src.${lndName}/lnd_comp_nuopc.F90"
  if [ "$do_runoff" = true ]; then
    exit_files_no_exist "$path_to_case/SourceMods/src.${rofName}/rof_comp_mct.F90" "$path_to_case/SourceMods/src.${rofName}/rof_comp_nuopc.F90"
  fi
fi

### Check required variables
betacast_required_vars=("casename" "atmDataType" "sstDataType" "numLevels" "numdays"
        "anl2mdlWeights" "PROJECTID" "DTIME" "FINERES" "USERSTAB")
check_required_vars "${betacast_required_vars[@]}"

# Boolean flags for Python that need to be appended to any command line calls
# TC augmentation in atm_to_cam
if [ "$augment_tcs" = true ]; then
  AUGMENT_STR="--augment_tcs"
else
  AUGMENT_STR=""
fi
# Vortex modification in atm_to_cam (if standalone_vortex = true, use old code)
if [ -f "$vortex_namelist" ] && [ "$standalone_vortex" = false ]; then
  VORTEX_STR="--vortex_namelist $vortex_namelist"
else
  VORTEX_STR=""
fi

###################################################################################

# do some stability calcs
# if USERSTAB is 0, use internal calcs.
# if USERSTAB is <0, use se_nsplit=-1
# If USERSTAB >0, use the value in seconds for dt_dyn and calculate nsplit accordingly from DTIME
USERSTABTF=$(python -c "print('TRUE' if ${USERSTAB} > 0 else 'FALSE')")
if [ ${USERSTABTF} == 'FALSE' ] ; then
  if [ $(python -c "print('TRUE' if ${USERSTAB} < -0.001 else 'FALSE')") == 'FALSE' ]; then
    if [ -z "$FINERES" ] || [ "$FINERES" -eq 0 ] 2>/dev/null; then
      echo "ERROR: FINERES is zero, empty, or non-numeric: '$FINERES'"
      exit 1
    fi
    STABILITY=$(python -c "print(30./${FINERES}*450.)")
    VALIDSTABVAL=true
    echo "Dynamic stability for ne${FINERES} to be ${STABILITY} seconds"
  else
    STABILITY=-1
    VALIDSTABVAL=false   # Setting to false means we won't calculate nsplit for SE/HOMME later
    echo "Dynamic stability for ne${FINERES} to be ${STABILITY} and internal STABILITY setting"
  fi
else
  STABILITY=${USERSTAB}
  VALIDSTABVAL=true
  echo "User defined stability set to ${STABILITY}"
fi
echo "VALIDSTABVAL is set to "${VALIDSTABVAL}

## Create paths to generate initial files if they don't exist...
mkdir -p ${pathToINICfiles}
mkdir -p ${pathToSSTfiles}

# Set mapping file directory and make if needed
mapping_files_path=${path_to_inputdata}/mapping/ ; mkdir -p ${mapping_files_path}

# Set timestamp for backing up files, etc.
timestamp=$(date +%Y%m%d.%H%M)
uniqtime=$(date +"%s%N")

echo "We are using ${casename} for the case"

if [ "$runmodel" = true ] ; then

############################### GET DYCORE INFO ###############################

cdv "$path_to_case"
DYCORE=$(./xmlquery CAM_DYCORE | sed 's/^[^\:]\+\://' | xargs)
if [ -z "$DYCORE" ] && [ "$modelSystem" -ne 2 ]; then
  echo "Couldn't automagically figure out dycore, assuming SE/HOMME"
  DYCORE="se"
elif [ "$modelSystem" -eq 2 ]; then
  echo "Using SCREAM, setting dycore to SCREAM/EAMXX"
  DYCORE="scream"
fi

echo "DYCORE: $DYCORE"

if [ "$debug" = false ] ; then

############################### GET ATM DATA ###############################

  # Initialize any "global" variables needed when getting atmospheric data
  RDADIR="" # Init to empty, but fill in if RDA available later
  ERA5RDA=0 # Set whether or not ERA5 is local (0 = local, 1 = RDA)

  echo "FLAGS at get_atm_data"
  check_shell_flags

  case "$atmDataType" in
    1)
      get_gfs_atm
      ;;
    2)
      get_era_interim_atm
      ;;
    3)
      get_cfsr_atm
      ;;
    4)
      get_era5_atm
      ;;
    9)
      echo "User has specified a CESM/E3SM data file"
      ;;
    *)
      echo "Incorrect model IC entered"
      exit 1
      ;;
  esac

############################### GET SST / NCL ###############################

  case "$sstDataType" in
    1)
      get_gdas_sst
      ;;
    2)
      get_erai_sst
      ;;
    3)
      get_noaaoi_sst
      ;;
    9)
      echo "User has specified a CESM/E3SM data stream"
      ;;
    *)
      echo "Incorrect SST data type entered"
      exit 1
      ;;
  esac

  # If not using data streams, we have to generate the SST forcing
  if [ "$sstDataType" -ne 9 ] ; then

    echo "Doing sst_to_cam"
    echo "cd'ing to interpolation directory: ${sst_to_cam_path}"
    cdv "${sst_to_cam_path}"

    sst_domain_file=${BETACAST}/grids/domains/domain.ocn.${docnres}.nc
    sst_scrip_file=${BETACAST}/grids/domains/scrip.ocn.${docnres}.nc
    sst_ESMF_file=${BETACAST}/grids/domains/ESMF.ocn.${docnres}.nc

    export BETACAST \
           sst_domain_file sst_scrip_file sst_ESMF_file \
           sstDataType predict_docn sst_to_cam_path \
           docnres DO_PYTHON cime_coupler \
           yearstr monthstr daystr cyclestr \
           SSTTYPE sst_files_path sstFile iceFile sstFileIC

    if [ "$batch_sst_gen" = true ]; then
      echo "Submitting sst_to_cam as batch job..."
      if [ "${BATCH_SCHEDULER}" = "sbatch" ]; then
        sbatch --wait \
          -C cpu \
          -A m2637 \
          --qos=debug \
          --time=00:30:00 \
          --nodes=1 \
          --ntasks-per-node=1 \
          --mem=0 \
          --output=${BETACAST}/sst_to_cam_%j.out \
          --error=${BETACAST}/sst_to_cam_%j.out \
          -J sst_to_cam \
          --export=ALL \
          ${BETACAST}/runner_sst_to_cam.sh
      elif [ "${BATCH_SCHEDULER}" = "qsub" ]; then
        qsub -V -W block=true \
          -A P93300042 \
          -q develop \
          -l walltime=00:30:00,select=1:ncpus=1:mpiprocs=1:ompthreads=1:mem=80GB \
          -N sst_to_cam \
          -o ${BETACAST}/sst_to_cam_"$(date +%Y%m%d%H%M%S)".out \
          -j oe \
          -- ${BETACAST}/runner_sst_to_cam.sh
      else
        echo "Error: Unknown BATCH_SCHEDULER='${BATCH_SCHEDULER}'. Must be 'sbatch' or 'qsub'." >&2
        exit 1
      fi
      echo "sst_to_cam batch job finished."
    else
      # Just run serially
      bash ${BETACAST}/runner_sst_to_cam.sh
    fi

  fi

  ############################### ATM NCL ###############################

  echo "Doing atm_to_cam"
  echo "cd'ing to interpolation directory: $atm_to_cam_path"
  cdv "$atm_to_cam_path"

  export BETACAST atm_to_cam_path mapping_files_path modelgridfile m2m_gridfile \
         m2m_parent_source uniqtime sePreFilterIC sstFileIC \
         atmDataType numLevels yearstr monthstr daystr cyclestr \
         RDADIR anl2mdlWeights DO_PYTHON DYCORE \
         do_frankengrid model_scrip standalone_vortex add_noise add_perturbs \
         perturb_namelist vortex_namelist modelSystem sstDataType \
         gfs_files_path era_files_path ERA5RDA

  if [ "$batch_atm_gen" = true ]; then
    echo "Submitting atm_to_cam as batch job..."
    if [ "${BATCH_SCHEDULER}" = "sbatch" ]; then
      sbatch --wait \
        -C cpu \
        -A m2637 \
        --qos=debug \
        --time=00:30:00 \
        --nodes=1 \
        --ntasks-per-node=1 \
        --mem=0 \
        --output=${BETACAST}/atm_to_cam_%j.out \
        --error=${BETACAST}/atm_to_cam_%j.out \
        -J atm_to_cam \
        --export=ALL \
        ${BETACAST}/runner_atm_to_cam.sh
    elif [ "${BATCH_SCHEDULER}" = "qsub" ]; then
      qsub -V -W block=true \
        -A P93300042 \
        -q develop \
        -l walltime=00:30:00,select=1:ncpus=1:mpiprocs=1:ompthreads=1:mem=80GB \
        -N atm_to_cam \
        -o ${BETACAST}/atm_to_cam_"$(date +%Y%m%d%H%M%S)".out \
        -j oe \
        -- ${BETACAST}/runner_atm_to_cam.sh
    else
      echo "Error: Unknown BATCH_SCHEDULER='${BATCH_SCHEDULER}'. Must be 'sbatch' or 'qsub'." >&2
      exit 1
    fi
    echo "atm_to_cam batch job finished."
  else
    # Just run serially
    bash ${BETACAST}/runner_atm_to_cam.sh
  fi

fi # debug

cdv "$path_to_case"

############################### CISM SETUP ###############################

if [ -f user_nl_cism ]; then
  sed -i '/dt_count/d' user_nl_cism
  echo "dt_count = 8" >> user_nl_cism
fi

############################### CLM SETUP ###############################

echo "................. configuring land and/or runoff"

## TEMPORARY CLM -> LAND FIX
## Check if clmstart exists but landstart doesn't, move if that is the case
if [[ ! -d "$landdir" && -d "${path_to_rundir}/clmstart" ]]; then
  echo "Moving ${path_to_rundir}/clmstart to ${landdir}"
  mv -v "${path_to_rundir}/clmstart" "$landdir"
else
  echo "No need to modify land directory; either it exists already or will be created later."
fi
## END TEMP FIX

echo "-------\nUSER_NL: Setting input ${lndName} dataset"
# Clean up file to delete any special interp lines that may be needed later (but aren't needed for native init)
sed -i '/init_interp_fill_missing_with_natveg/d' user_nl_${lndName}
sed -i '/use_init_interp/d' user_nl_${lndName}

# Create a temp directory for now
RUNTMPDIR="${path_to_nc_files}/tmp"
[ -z "$RUNTMPDIR" ] && { echo "RUNTMPDIR is not set. Exiting."; exit 1; }
# Create directory and make sure it was actually generated (to prevent user from not having write perms)
mkdir -vp "$RUNTMPDIR"
[ -d "$RUNTMPDIR" ] || { echo "Failed to create $RUNTMPDIR. Exiting."; exit 1; }
# Delete files in RUNTMPDIR
find "$RUNTMPDIR" -type f -exec rm -v {} +

# We want to check ${landdir} for land restart file with exact date match. If so, use that.
landrestartfile=$(find "${landdir}" -maxdepth 1 \( -type f -o -type l \) -name "${casename}.${lndSpecialName}.r.${yearstr}-${monthstr}-${daystr}-${cyclestrsec}.*" | head -n 1)

# Check to see if file exists on native SE land grid
if [ -n "${landrestartfile}" ] && [ -f "${landrestartfile}" ]; then
  echo "File ${landrestartfile} exists at exact time -- in landstart folder"
else
  # If we don't have an exact match from above, trying 00Z first, then check landrawdir
  echo "File ${landrestartfile} does not exist at exact time -- in landstart folder"
  landrestartfile=$(find "${landdir}" -maxdepth 1 -type f -name "${casename}.${lndSpecialName}.r.${yearstr}-${monthstr}-${daystr}-00000.*" | head -n 1)
  if [ -n "${landrestartfile}" ] && [ -f "${landrestartfile}" ]; then
    echo "File ${landrestartfile} exists at 00Z -- in landstart folder"
  else
    rawlandrestartfile=$(find "${landrawdir}" -maxdepth 1 -type f -name "*.${lndSpecialName}.r.${yearstr}-${monthstr}-${daystr}-${cyclestrsec}.*" | head -n 1)
    echo "rawlandrestartfile: ${rawlandrestartfile} -- in raw folder"
    if [ -n "${rawlandrestartfile}" ] && [ -f "${rawlandrestartfile}" ]; then
      landrestartfile=${rawlandrestartfile}
    else
      echo "No LND restart file exists, setting landrestartfile to empty string."
      landrestartfile=
    fi
  fi
fi
echo "--> FINAL landrestartfile: ${landrestartfile}"

if [ -f "${landrestartfile}" ] && [[ "${landrestartfile}" != *.nc ]]; then
  echo "${landrestartfile} was found, but does not have an *.nc extension. Assuming compressed. Copying..."
  cp -v ${landrestartfile} "$RUNTMPDIR"
  try_uncompress "$RUNTMPDIR/$(basename "${landrestartfile}")"
  landrestartfile=$(find "${RUNTMPDIR}" -maxdepth 1 -type f -name "*.${lndSpecialName}.r.${yearstr}-${monthstr}-${daystr}-*.nc" | head -n 1)
  echo "--> FINAL (uncompressed) landrestartfile: ${landrestartfile}"
fi

## Now modify user_nl_${lndName}
xmlchange_verbose "${lndName^^}_FORCE_COLDSTART" "off"

if [ "${landrestartfile}" ] ; then
  sed -i '/.*finidat/d' user_nl_${lndName}
  echo "finidat='${landrestartfile}'" >> user_nl_${lndName}
  if [ -n "${rawlandrestartfile-}" ] && [ -f "${rawlandrestartfile-}" ]; then
    echo "Adding interpolation lines to user_nl_${lndName}"
    echo "use_init_interp = .true." >> user_nl_${lndName}
    echo "init_interp_fill_missing_with_natveg = .true." >> user_nl_${lndName}
  fi
else
  if [ -f user_nl_${lndName}_presave ] ; then
    echo "Using pre-written user_nl_${lndName} file"
    cp -v user_nl_${lndName}_presave user_nl_${lndName}
  else
    echo "WARNING: Land file DOES NOT EXIST, will use arbitrary user_nl_${lndName} already in folder"
    echo "!!!!!!!!!!!!!"
    sed -i '/.*finidat/d' user_nl_${lndName}
    xmlchange_verbose "${lndName^^}_FORCE_COLDSTART" "on"
  fi
fi

## Append a check error to ignore inconsistencies in the dataset
if [ "$modelSystem" -eq 0 ]; then   # CLM/CTSM
  sed -i '/check_finidat_pct_consistency/d' user_nl_${lndName}
  sed -i '/check_finidat_year_consistency/d' user_nl_${lndName}
  echo "check_finidat_pct_consistency = .false." >> user_nl_${lndName}
  echo "check_finidat_year_consistency = .false." >> user_nl_${lndName}
elif [ "$modelSystem" -eq 1 ] || [ "$modelSystem" -eq 2 ]; then   # ELM
  sed -i '/check_finidat_fsurdat_consistency/d' user_nl_${lndName}
  echo "check_finidat_fsurdat_consistency = .false." >> user_nl_${lndName}
  # 2/25/24 CMZ added since ELM doesn't have use_init_interp support for rawlandrestartfile
  if [ -n "${rawlandrestartfile-}" ]; then # if rawlandrestartfile is SET *and* not empty...
    echo "WARNING USER_NL: We used a rawlandrestartfile, but E3SM/ELM doesn't support interpolation, removing relevant lines"
    echo "WARNING USER_NL: If your model fails, check to make sure the rawlandrestartfile supports your target land grid"
    sed -i '/use_init_interp/d' user_nl_${lndName}
    sed -i '/init_interp_fill_missing_with_natveg/d' user_nl_${lndName}
  fi
else
  echo "Unknown modeling system set for modelSystem: $modelSystem"
  exit 1
fi

############################### ROF SETUP ###############################

if [ "$do_runoff" = true ]; then

  echo "-------\nUSER_NL: Setting input ${rofName} dataset"

  # Delete any existing input data
  sed -i '/.*finidat/d' user_nl_${rofName}

  # We want to check ${landdir} for land restart files. If so, use those.
  rofrestartfile=$(find "${landdir}" -maxdepth 1 \( -type f -o -type l \) -name "${casename}.${rofSpecialName}.r.${yearstr}-${monthstr}-${daystr}-${cyclestrsec}.*" | head -n 1)
  echo "$rofrestartfile"

  # Check to see if file exists on native SE land grid
  if [ -n "${rofrestartfile}" ] && [ -f "${rofrestartfile}" ]; then
     echo "USER_NL: ${rofName} file exists at exact time -- in landstart folder"
  else
     echo "USER_NL: ${rofName} file does not exist at exact time -- in landstart folder"
     rofrestartfile=$(find "${landdir}" -maxdepth 1 -type f -name "${casename}.${rofSpecialName}.r.${yearstr}-${monthstr}-${daystr}-00000.*" | head -n 1)
     if [ -n "${rofrestartfile}" ] && [ -f "${rofrestartfile}" ]; then
       echo "USER_NL: ${rofName} file exists at 00Z -- in landstart folder"
    else
      rawrofrestartfile=$(find "${landrawdir}" -maxdepth 1 -type f -name "*.${rofSpecialName}.r.${yearstr}-${monthstr}-${daystr}-${cyclestrsec}.*" | head -n 1)
      echo "rawrofrestartfile: ${rawrofrestartfile} -- in raw folder"
      if [ -n "${rawrofrestartfile}" ] && [ -f "${rawrofrestartfile}" ]; then
        rofrestartfile=${rawrofrestartfile}
      else
        echo "No ROF restart file exists, setting rofrestartfile to empty string."
        rofrestartfile=
      fi
    fi
  fi
  echo "USER_NL: rofrestartfile: ${rofrestartfile}"

  if [ -f "${rofrestartfile}" ] && [[ "${rofrestartfile}" != *.nc ]]; then
    echo "${rofrestartfile} was found, but does not have an *.nc extension. Assuming compressed. Copying..."
    cp -v ${rofrestartfile} "$RUNTMPDIR"
    try_uncompress "$RUNTMPDIR/$(basename "$rofrestartfile")"
    rofrestartfile=$(find "${RUNTMPDIR}" -maxdepth 1 -type f -name "*.${rofSpecialName}.r.${yearstr}-${monthstr}-${daystr}-*.nc" | head -n 1)
    echo "Updated rofrestartfile: ${rofrestartfile}"
  fi

  ## Now modify user_nl_${rofName}
  if [ -n "${rofrestartfile-}" ]; then
    echo "USER_NL: Adding ${rofrestartfile} to user_nl_${rofName}"
    echo "${rof_finidat}='${rofrestartfile}'" >> user_nl_${rofName}
  else
    # Check to see if there is a raw ROF file in the landrawdir given by the user
    rawrofrestartfile=$(ls ${landrawdir}/*.${rofSpecialName}.r.${yearstr}-${monthstr}-${daystr}-${cyclestrsec}.nc || true)   # check for file, suppress failed ls error if true
    echo "rawrofrestartfile: ${rawrofrestartfile}"
    if [ ! -z "${rawrofrestartfile}" ]; then   # if rawrofrestartfile string is NOT empty, add it.
      echo "WARNING USER_NL: Adding rawrofrestartfile for runoff, although no check to see if grid is consistent!"
      echo "USER_NL: Adding ${rawrofrestartfile} to user_nl_${rofName}"
      echo "${rof_finidat}='${rawrofrestartfile}'" >> user_nl_${rofName}
    else
      echo "USER_NL: No ${rofName} restart file added, cold start"
    fi
  fi

fi  # if do_runoff

echo "................. done configuring land and/or runoff"

############################### GENERIC CAM SETUP ###############################

# Setting projectID and queue
xmlchange_verbose "PROJECT" "$PROJECTID"
xmlchange_verbose "JOB_QUEUE" "$RUNQUEUE" "--force"

# Turning off archiving and restart file output in env_run.xml
xmlchange_verbose "DOUT_S" "FALSE"
xmlchange_verbose "REST_OPTION" "nyears"
xmlchange_verbose "REST_N" "999"

if [ ${sstDataType} -ne 9 ] ; then
  # We are using some sort of analysis SST
  # Setting SST domain file
  if [ "${cime_coupler}" == "nuopc" ] ; then
    xmlchange_verbose "SSTICE_MESH_FILENAME" "${sst_ESMF_file}"
  else
    xmlchange_verbose "SSTICE_GRID_FILENAME" "${sst_domain_file}"
  fi
  # Setting SST from default to our SST
  xmlchange_verbose "SSTICE_DATA_FILENAME" "${sstFileIC}"
  # Standardizing streams for SST
  xmlchange_verbose "SSTICE_YEAR_START" "1"
  xmlchange_verbose "SSTICE_YEAR_END" "1"
else
  # We already have a DOCN stream available to use (sstDataType 9)
  # Reproducing previous DOCN configuration as specified by m2m_sst* in namelist
  if [ "${cime_coupler}" == "nuopc" ] ; then
    xmlchange_verbose "SSTICE_MESH_FILENAME" "${m2m_sst_grid_filename}"
  else
    xmlchange_verbose "SSTICE_GRID_FILENAME" "${m2m_sst_grid_filename}"
  fi
  xmlchange_verbose "SSTICE_DATA_FILENAME" "${m2m_sstice_data_filename}"
  xmlchange_verbose "SSTICE_YEAR_ALIGN" "${m2m_sstice_year_align}"
  xmlchange_verbose "SSTICE_YEAR_START" "${m2m_sstice_year_start}"
  xmlchange_verbose "SSTICE_YEAR_END" "${m2m_sstice_year_end}"

  # This is a cheat so that the m2m_sstice_data_filename is archived correctly.
  #sstFileIC=$m2m_sstice_data_filename
fi

# Getting GLC coupling to handle forecasts across calendar years
xmlchange_verbose "GLC_AVG_PERIOD" "glc_coupling_period"
#######
if [ "$land_spinup" = true ] ; then
  xmlchange_verbose "REST_OPTION" "ndays"
  xmlchange_verbose "REST_N" "1"
fi
#######

# Update name of atmospheric model if E3SM/SCREAM
# CMZ, may need to just do this for all modeling systems
#if [ $modelSystem -eq 1 ]; then
#  atmName=$(./xmlquery COMP_ATM | sed 's/^[^\:]\+\://' | xargs)
#  echo "Updating atmName to $atmName"
#fi

# Non-SCREAM models backup atm config file
if [ "$modelSystem" -ne 2 ]; then
  cp user_nl_${atmName} user_nl_${atmName}.BAK
fi

# Update env_run.xml with runtime parameters
xmlchange_verbose "RUN_STARTDATE" "$yearstr-$monthstr-$daystr"
xmlchange_verbose "START_TOD" "$cyclestrsec"
xmlchange_verbose "STOP_OPTION" "ndays"
xmlchange_verbose "STOP_N" "$numdays"

SEINIC=${sePreFilterIC}

############################### (IF) FILTER SETUP ###############################

if [ "$doFilter" = true ] ; then

  # Can we run filtering?
  if [ "$modelSystem" -ne 2 ]; then
    echo "SCREAM and digital filter not supported"
    exit
  fi

  # If filtering, need to change these options
  xmlchange_verbose "STOP_OPTION" "nhours"
  xmlchange_verbose "STOP_N" "$filterHourLength"

  sed -i '/.*ncdata/d' user_nl_${atmName}
  echo "ncdata='${SEINIC}'" >> user_nl_${atmName}

  # Do filter timestepping stability
  ATM_NCPL=192
  echo "ATM_NCPL: $ATM_NCPL  DTIME: $DTIME"
  if [[ "$DYCORE" == "se" ]]; then
    if [ "$use_nsplit" = true ]; then # if trad. SE nsplit timestepping
      echo "Using nsplit --> trad. SE nsplit timestepping"
      sed -i '/.*se_nsplit/d' user_nl_${atmName}
      if [ "$VALIDSTABVAL" = false ]; then
        # Use se_nsplit = -1, which is internal
        echo "SE_NSPLIT: -1"
        echo "se_nsplit=-1" >> user_nl_${atmName}
      else
        # Add se_nsplit based off of dtime=450 and stability
        SE_NSPLIT=$(python -c "from math import ceil; print(int(ceil(450/${STABILITY})))")
        echo "SE_NSPLIT: $SE_NSPLIT   STABILITY: $STABILITY   "
        echo "se_nsplit=${SE_NSPLIT}" >> user_nl_${atmName}
      fi
    else
      echo "Not using nsplit"
      sed -i '/.*se_tstep/d' user_nl_${atmName}
      if [ "$VALIDSTABVAL" = false ]; then
        echo "Betacast cannot handle when $VALIDSTABVAL is false and use_nsplit is false"
        echo "You need to set USERSTAB -> desired se_tstep in the namelist. Exiting..."
        exit
      else
        # Add se_tstep directly from stability
        echo "SE_TSTEP --> STABILITY: $STABILITY   "
        echo "se_tstep=${STABILITY}" >> user_nl_${atmName}
      fi
    fi
  else
    echo "non-SE core, make sure timestepping is happy!"
  fi
  xmlchange_verbose "ATM_NCPL" "$ATM_NCPL"

  # Set NHTFRQ, MFILT, and FINCL fields
  sed -i '/.*nhtfrq/d' user_nl_${atmName}
  sed -i '/.*mfilt/d' user_nl_${atmName}
  sed -i '/.*fincl/d' user_nl_${atmName}
  sed -i '/.*empty_htapes/d' user_nl_${atmName}
  sed -i '/.*inithist/d' user_nl_${atmName}
  echo "nhtfrq=1" >> user_nl_${atmName}
  echo "mfilt=200" >> user_nl_${atmName}
  echo "fincl1='U:I','V:I','T:I','PS:I','Q:I','CLDICE:I','CLDLIQ:I','NUMICE:I','NUMLIQ:I','ICEFRAC:I','SNOWHICE:I'" >> user_nl_${atmName}
  echo "empty_htapes=.TRUE." >> user_nl_${atmName}
  echo "inithist='6-HOURLY'" >> user_nl_${atmName}

  xmlchange_verbose "JOB_WALLCLOCK_TIME" "$FILTERWALLCLOCK"
  xmlchange_verbose "JOB_QUEUE" "$FILTERQUEUE" "--force"

  if [ "$debug" = false ] ; then

    echo "Begin call to filter-run"
    set +e
    run_CIME2 "$path_to_rundir" "$CIMEsubstring" "$CIMEbatchargs" true
    CIMESTATUS=$?
    set -e

    ## Run NCL filter
    echo "Running filter, currently in $PWD"
    cp -v ${sePreFilterIC} ${sePostFilterIC}
    filtfile_name=${casename}.${atmName}.h0.$yearstr-$monthstr-$daystr-$cyclestrsec.nc

    if [ "$DO_PYTHON" = true ]; then
      (set -x; python $filter_path/lowmemfilter.py \
          --endhour ${filterHourLength} \
          --tcut ${filtTcut} \
          --filtfile_name "${path_to_rundir}/${filtfile_name}" \
          --writefile_name "${sePostFilterIC}"
      )
    else
      set +e # Turn off error check
      (set -x; ncl $filter_path/lowmemfilter.ncl \
       endhour=${filterHourLength} tcut=${filtTcut} \
      'filtfile_name = "'${path_to_rundir}'/'${filtfile_name}'"' \
      'writefile_name = "'${sePostFilterIC}'"' ) ; exit_status=$?
      check_ncl_exit "lowmemfilter.ncl" $exit_status
      set -e # Turn on error check
    fi

  fi  # debug

  echo "done with filter, removing filter files"
  mkdir -p $path_to_nc_files/filtered
  mv -v $path_to_nc_files/*.nc $path_to_nc_files/filtered
  ## Delete filter files that aren't h0 since those are the only ones we care about.
  find $path_to_nc_files/filtered/ -type f -not -name '*.${atmName}.h0.*.nc' | xargs rm

  ## For now, I'm going to go back and delete all the h0 files in filtered
  rm -v $path_to_nc_files/filtered/*.${atmName}.h0.*.nc

  ### Output filter files only, can be used for ensemble or other initialization after the fact
  if [ "$filterOnly" = true ] ; then
    FILTONLYDIR=${path_to_nc_files}/FILT_INIC/${se_yearstr}-${se_monthstr}-${se_daystr}-${se_cyclestrsec}
    mkdir -p ${FILTONLYDIR}
    cp -v ${sePostFilterIC} ${FILTONLYDIR}/${casename}_FILTERED_${se_yearstr}-${se_monthstr}-${se_daystr}-${se_cyclestrsec}.nc
    cp -v ${sstFileIC} ${FILTONLYDIR}/${casename}_SST_1x1_${se_yearstr}-${se_monthstr}-${se_daystr}-${se_cyclestrsec}.nc
    mkdir ${FILTONLYDIR}/config_files
    mv -v ${path_to_nc_files}/*.gz ${FILTONLYDIR}/config_files
    mv -v ${path_to_nc_files}/*.nml ${FILTONLYDIR}/config_files
    mv -v ${path_to_nc_files}/*_in ${FILTONLYDIR}/config_files
    exit
  fi

  # Special things we have to do to "reset" model after filtering
  echo "Pushing model start back from $cyclestrsec to $se_cyclestrsec seconds..."
  cdv "$path_to_case"
  xmlchange_verbose "RUN_STARTDATE" "$se_yearstr-$se_monthstr-$se_daystr"
  xmlchange_verbose "START_TOD" "$se_cyclestrsec"
  xmlchange_verbose "STOP_OPTION" "ndays"
  xmlchange_verbose "STOP_N" "$numdays"

  # Set new inic to point to post filter file
  SEINIC=${sePostFilterIC}
fi

############################### "ACTUAL" FORECAST RUN ###############################

# Set NHTFRQ, MFILT, and FINCL fields for models that user user_nl_atm
if [ "$modelSystem" -eq 0 ] || [ "$modelSystem" -eq 1 ]; then
  ## Initial modification of user_nl_atm
  sed -i '/.*nhtfrq/d' user_nl_${atmName}
  sed -i '/.*mfilt/d' user_nl_${atmName}
  sed -i '/.*fincl/d' user_nl_${atmName}
  sed -i '/.*empty_htapes/d' user_nl_${atmName}
  sed -i '/.*collect_column_output/d' user_nl_${atmName}
  sed -i '/.*avgflag_pertape/d' user_nl_${atmName}  # Note, we delete this and user either specifies as :A, :I for each var or as a sep var
  echo "empty_htapes=.TRUE." >> user_nl_${atmName}
  sed -i '/.*inithist/d' user_nl_${atmName}
  if [ "$save_nudging_files" = true ] ; then
    echo "inithist='6-HOURLY'" >> user_nl_${atmName}
  else
    echo "inithist='NONE'" >> user_nl_${atmName}
  fi
  # Concatenate output streams to end of user_nl_${atmName}
  cat ${OUTPUTSTREAMS} >> user_nl_${atmName}
elif [ "$modelSystem" -eq 2 ]; then
  # Clear any existing output datastream YAML links
  xmlchange_verbose --atmchange "output_yaml_files" ""
  # Here, we unpack the packed output streams as separate yaml files
  unpack_files ${OUTPUTSTREAMS} "./"
  # Then append the streams
  scream_atmchange_from_packed ${OUTPUTSTREAMS}
fi

# Calculate timestep stability criteria
ATM_NCPL=$(python -c "print(int(86400/${DTIME}))")
echo "ATM_NCPL: $ATM_NCPL  DTIME: $DTIME"
if [[ "$DYCORE" == "se" ]]; then
  if [ "$use_nsplit" = true ]; then # if trad. SE nsplit timestepping
    # Delete existing se_nsplit
    sed -i '/.*se_nsplit/d' user_nl_${atmName}
    if [ "$VALIDSTABVAL" = false ]; then
      echo "SE_NSPLIT: -1 "
      echo "se_nsplit=-1" >> user_nl_${atmName}
    else
      # Add se_nsplit based off of dtime and stability
      SE_NSPLIT=$(python -c "from math import ceil; print(int(ceil(${DTIME}/${STABILITY})))")
      echo "SE_NSPLIT: $SE_NSPLIT   STABILITY: $STABILITY   "
      echo "se_nsplit=${SE_NSPLIT}" >> user_nl_${atmName}
    fi
  else # otherwise, use E3SM se_tstep
    # Delete existing tstep
    sed -i '/.*se_tstep/d' user_nl_${atmName}
    if [ "$VALIDSTABVAL" = false ]; then
      echo "Betacast cannot handle when $VALIDSTABVAL is false and use_nsplit is false"
      echo "You need to set USERSTAB -> desired se_tstep in the namelist. Exiting..."
      exit
    else
      # Add se_tstep directly from stability
      echo "SE_TSTEP --> STABILITY: $STABILITY   "
      echo "se_tstep=${STABILITY}" >> user_nl_${atmName}
    fi
  fi
elif [[ "$DYCORE" == "scream" ]]; then
  # SCREAM/EAMXX must have use_nsplit = false, so error if true
  if [ "$use_nsplit" = true ]; then
    echo "SCREAM/EAMXX cannot have use_nsplit be true, exiting"
    exit 1
  else
    if [ "$VALIDSTABVAL" = false ]; then
      echo "Betacast cannot handle when $VALIDSTABVAL is false and use_nsplit is false"
      echo "You need to set USERSTAB -> desired se_tstep in the namelist. Exiting..."
      exit
    else
      # Add se_tstep directly from stability
      echo "SE_TSTEP --> STABILITY: $STABILITY   "
      xmlchange_verbose --atmchange "se_tstep" "${STABILITY}"
    fi
  fi
else
  echo "non-SE/HOMME core, make sure timestepping is happy!"
fi

# Set dtime
xmlchange_verbose "ATM_NCPL" "$ATM_NCPL"

# Set queue settings
xmlchange_verbose "JOB_WALLCLOCK_TIME" "$RUNWALLCLOCK"
xmlchange_verbose "JOB_QUEUE" "$RUNQUEUE" "--force"
if [ -n "${RUNPRIORITY+x}" ] && [ -n "$RUNPRIORITY" ]; then
  xmlchange_verbose "JOB_PRIORITY" "$RUNPRIORITY" "--force"
fi

if [ "$modelSystem" -eq 0 ] || [ "$modelSystem" -eq 1 ]; then
  # Delete existing ncdata and inject new one into user_nl_atm file
  sed -i '/.*ncdata/d' user_nl_${atmName}
  echo "ncdata='${SEINIC}'" >> user_nl_${atmName}
elif [ "$modelSystem" -eq 2 ]; then
  xmlchange_verbose --atmchange "filename" "${SEINIC}"
else
  echo "Cannot inject initial conditions file!"
  exit
fi

if [ "$debug" = false ] ; then
  echo "Begin call to forecast run"
  #run_CIME2 "$path_to_rundir" "$CIMEsubstring" "$CIMEbatchargs" true
  CIMESTATUS=1; CIMEITER=0
  while [ "$CIMESTATUS" != 0 ] ; do
    if [ "$CIMEITER" -ge "$CIMEMAXTRIES" ]; then
      echo "Exceeded the max tries: $CIMEMAXTRIES ... exiting"
      exit 1
    fi
    CIMEITER=$((CIMEITER+1))
    echo "Running CIME on CIMEITER $CIMEITER"
    set +e
    run_CIME2 "$path_to_rundir" "$CIMEsubstring" "$CIMEbatchargs" false
    CIMESTATUS=$?
    set -e
  done
  echo "Returned status $CIMESTATUS"
fi

fi # end run model

cdv "$outputdir"
check_shell_flags

# Generate folder structure and move NetCDF files
# This also moves nl_files and logs
timer main_archive "$tmparchivecdir" "$atmName" "$lndName" "$rofName"

# Copy betacast configs to archive directory for posterity
safe_cp_files "$MACHINEFILE" "$NAMELISTFILE" "$OUTPUTSTREAMS" "$perturb_namelist" "$tmparchivecdir/betacast"
# Archive shell configuration for debugging
declare -p | sort > "$tmparchivecdir/betacast/shellvars_$(date +%Y%m%d_%H%M%S).txt"

# Copy user_nl* files to archive dir as well...
safe_cp2 "${path_to_case}/user_nl_*" "${tmparchivecdir}/nl_files/"
safe_cp2 "${path_to_case}/*.xml" "${tmparchivecdir}/nl_files/"

# Archive initial conditions?
if [ "$archive_inic" = true ]; then
  echo "BETACAST_USER: requesting initial condition files be archived"
  timer archive_inic "$tmparchivecdir" "$path_to_case" "$compress_history_nc" "$atmName" "$lndName" "$rofName" "$sstFileIC"
fi

# Archive nudging files generated by hindcasts
if [ "$save_nudging_files" = true ] ; then
  echo "BETACAST_USER: requesting nudging files be archived"
  timer archive_nudging "$tmparchivecdir" "$path_to_nc_files" "$compress_history_nc"
fi

## Move land files to new restart location
cdv "$path_to_nc_files" || { echo "Failed to change directory to $path_to_nc_files"; exit 1; }
if [ "$keep_land_restarts" = true ]; then
  echo "BETACAST_USER: Archiving land restart files"
  mkdir -p "$landdir" #make landdir if doesn't exist
  echo "Removing 06Z and 18Z land restart files if they exist"
  rm -v *."${lndName}"*.r.*32400.nc || true
  rm -v *."${lndName}"*.r.*75600.nc || true
  echo "Moving land restart files for future runs"
  for file in *."${lndName}"*.r.*.nc; do
    [ -e "$file" ] || continue # Skip if no files match
    compress_file "$file" zstd
    mv -v "${file}"* "$landdir" || true
  done
  ## Move runoff files to land dir if doing runoff
  if [ "$do_runoff" = true ]; then
    echo "Removing 06Z and 18Z runoff restart files if they exist"
    rm -v *."${rofName}"*.r.*32400.nc || true
    rm -v *."${rofName}"*.r.*75600.nc || true
    echo "Moving runoff restart files for future runs"
    for file in *."${rofName}"*.r.*.nc; do
      [ -e "$file" ] || continue # Skip if no files match
      compress_file "$file" zstd
      mv -v "${file}"* "$landdir" || true
    done
  fi
else
  echo "BETACAST_USER: Removing all land restart files!"
  rm -v *."${lndName}"*.r.*.nc || true
  if [ "$do_runoff" = true ]; then
    rm -v *."${rofName}"*.r.*.nc || true
  fi
fi

## Delete any leftover files in the run dir that we don't need/want anymore
timer delete_leftovers "$path_to_nc_files" "$atmName" "$lndName" "$rofName"

## Delete run temp directory
rm -rfv "$RUNTMPDIR"

## Move directory to YYYYMMDDHH, handing duplicates
cdv "$path_to_nc_files"
ARCHIVESUBDIR="${yearstr}${monthstr}${daystr}${cyclestr}"
echo "Moving tmp archive directory to ARCHIVEDIR/YYYYMMDDHH --> ${ARCHIVESUBDIR}"
if [ -d "${ARCHIVEDIR}/${ARCHIVESUBDIR}" ]; then
  echo "${ARCHIVEDIR}/${ARCHIVESUBDIR} already exists! Appending date string and moving!"
  mv -v "${ARCHIVEDIR}/${ARCHIVESUBDIR}" "${ARCHIVEDIR}/${ARCHIVESUBDIR}_$(date +%Y%m%d%H%M)"
fi
mv -v "$tmparchivecdir" "${ARCHIVEDIR}/${ARCHIVESUBDIR}"

if [ "$dotracking" = true ] ; then
  echo "BETACAST_USER: requesting cyclones to be tracked."

  # Go to cyclone tracking folder...
  pushd "${BETACAST}/cyclone-tracking/" > /dev/null

  # Set a few things that are hardcoded
  TCVITFOLDER="./fin-tcvitals/"
  TCVITFILE="${TCVITFOLDER}/tcvitals.${yearstr}${monthstr}${daystr}${cyclestr}"
  ATCFFILE="atcf.${casename}.${yearstr}${monthstr}${daystr}${cyclestr}"

  # Get TC vitals for YYYYMMDDHH and store in $TCVITFOLDER
  (
    set -x
    /bin/bash ./get-vitals.sh "${yearstr}${monthstr}${daystr}${cyclestr}" "${TCVITFOLDER}"
  )

  # If we have vitals, run tracking
  if [ -s "${TCVITFILE}" ] ; then  # if file does exist (i.e., it was downloaded) and has size >0, run tracker
    echo "Found a TCvitals with at least one storm, running tracking code"
    echo "${TCVITFILE}" ; head -100 "${TCVITFILE}"
    (
      set -x
      bash ./cyclone-tracking-driver.sh "${yearstr}${monthstr}${daystr}${cyclestr}" \
        "${casename}" \
        "${TCVITFILE}" \
        "${ATCFFILE}" \
        "${track_connectfile}" \
        "${path_to_rundir}" \
        "${track_sendhtml}" \
        "${track_hstream}" \
        "${track_stride}" \
        "${track_ATCFTECH}" \
        "${TE_SERIAL_DIR}" \
        "${ARCHIVEDIR}/${ARCHIVESUBDIR}"
    ) || { echo "Tracking code failed"; popd > /dev/null; }
    safe_cp_files "trajs.trajectories.txt.${casename}.png" "./fin-figs/trajs.trajectories.txt.${casename}.${yearstr}${monthstr}${daystr}${cyclestr}.png"
  else
    echo "No TCvitals file exists and/or no storms on said TCvitals file, no reason to run the tracking code"
  fi

  # Return to where we were...
  popd > /dev/null
fi

if [ "$sendplots" = true ]; then
  echo "BETACAST_USER: Sending plots to remote server!"
  upload_ncl_script="${BETACAST}/upload_ncl.sh"
  temp_upload_script="${upload_ncl_script}.${uniqtime}.ncl"
  cp "${upload_ncl_script}" "${temp_upload_script}"
  sed -i -e "s?.*yearstr=.*?yearstr='${yearstr}'?" \
         -e "s?.*monthstr=.*?monthstr='${monthstr}'?" \
         -e "s?.*daystr=.*?daystr='${daystr}'?" \
         -e "s?.*cyclestrsec=.*?cyclestrsec='${cyclestrsec}'?" \
         -e "s?.*cyclestr=.*?cyclestr='${cyclestr}'?" \
         "${temp_upload_script}" || { echo "Failed to update script"; exit 1; }
  # Call script to create plots and push to remove server
  (
    set -x
    /bin/bash "${temp_upload_script}" "${nclPlotWeights}" "${outputdir}/${yearstr}${monthstr}${daystr}${cyclestr}" "${BETACAST}"
  )
  # Cleanup
  rm -f "${temp_upload_script}"
fi

# Compress model output streams
# Let's do this last so all the above scripts can operate on uncompressed files
if [ "$compress_history_nc" = true ]; then
  echo "BETACAST_USER: Requesting history be compressed"
  timer compress_history "${ARCHIVEDIR}/${ARCHIVESUBDIR}" $GNUPAR_AVAIL
fi

# Tar archive files so they are contained in a single directory
if [ "$tararchivedir" = true ] ; then
  echo "BETACAST_USER: Requesting archive dir be tarred"
  tar -cvf ${ARCHIVEDIR}/${ARCHIVESUBDIR}.tar -C ${ARCHIVEDIR} ${ARCHIVESUBDIR}
  rm -rfv ${ARCHIVEDIR}/${ARCHIVESUBDIR} || true
fi

# Timing statistics
script_end=$(date +%s)
print_elapsed_time "$script_start" "$script_end"

### If not live and the run has made it here successively, delete top line of datesfile
if [ "$islive" = false ] ; then
  cdv ${BETACAST}
  #Remove top line from dates file
  remove_top_line_from_dates ${datesfile}

  AUTORESUB="yes"
  if [ "$AUTORESUB" == "yes" ]; then
    echo "*-*-*-* Automatically resubbing next date!"

    case "$BATCH_PREFIX" in
      sbatch)
        echo "(Re)submitting: sbatch $0 $MACHINEFILE $NAMELISTFILE $OUTPUTSTREAMS"
        exec sbatch "$0" "$MACHINEFILE" "$NAMELISTFILE" "$OUTPUTSTREAMS"
        ;;
      qsub)
        echo "(Re)submitting: qsub -v MACHINEFILE=$MACHINEFILE,NAMELISTFILE=$NAMELISTFILE,OUTPUTSTREAMS=$OUTPUTSTREAMS $0"
        exec qsub -v MACHINEFILE="$MACHINEFILE",NAMELISTFILE="$NAMELISTFILE",OUTPUTSTREAMS="$OUTPUTSTREAMS" "$0"
        ;;
      *)
        echo "(Re)submitting directly: $0 $MACHINEFILE $NAMELISTFILE $OUTPUTSTREAMS"
        exec "$0" "$MACHINEFILE" "$NAMELISTFILE" "$OUTPUTSTREAMS"
        ;;
    esac

  fi
fi

echo "DONE at $(date)"

exit 0
