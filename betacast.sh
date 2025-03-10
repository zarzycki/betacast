#!/bin/bash

##=======================================================================
#PBS -N wx-driver
#PBS -A P54048000
#PBS -l walltime=6:00:00
#PBS -q share
#PBS -k oe
#PBS -m a
#PBS -M zarzycki@ucar.edu
#PBS -l select=1:ncpus=8
##=======================================================================

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
#
# Details can be found in:
# C. M. Zarzycki and C. Jablonowski (2015), Experimental tropical cyclone forecasts
# using a variable-resolution global model. Mon. Weat. Rev., 143(10), 4012–4037.
# doi:10.1175/MWR-D-15-0159.1.
###################################################################################

script_start=$(date +%s)
set -e
#set -v

SCRIPTPATH=$(dirname "$(realpath "$0")")
echo "Our script path is $SCRIPTPATH"
source ${SCRIPTPATH}/utils.sh         # Source external bash functions
source ${SCRIPTPATH}/datahelpers.sh   # Source external bash functions

# Set files
MACHINEFILE=${1}
NAMELISTFILE=${2}
OUTPUTSTREAMS=${3}
# If relative path, convert to absolute path
if [[ "$MACHINEFILE" != /* ]] && [[ "$MACHINEFILE" != ~* ]]; then MACHINEFILE=${PWD}/${MACHINEFILE}; fi
if [[ "$NAMELISTFILE" != /* ]] && [[ "$NAMELISTFILE" != ~* ]]; then NAMELISTFILE=${PWD}/${NAMELISTFILE}; fi
if [[ "$OUTPUTSTREAMS" != /* ]] && [[ "$OUTPUTSTREAMS" != ~* ]]; then OUTPUTSTREAMS=${PWD}/${OUTPUTSTREAMS}; fi
exit_file_no_exist $MACHINEFILE
exit_file_no_exist $NAMELISTFILE
exit_file_no_exist $OUTPUTSTREAMS
echo $MACHINEFILE; echo $NAMELISTFILE; echo $OUTPUTSTREAMS

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
gfs_to_cam_path=${SCRIPTPATH}/gfs_to_cam
era_to_cam_path=${SCRIPTPATH}/interim_to_cam
atm_to_cam_path=${SCRIPTPATH}/atm_to_cam
sst_to_cam_path=${SCRIPTPATH}/sst_to_cam
filter_path=${SCRIPTPATH}/filter
###################################################################################

### setting variables not included in namelist for backwards compat
if [ -z "${CIMEbatchargs+x}" ]; then CIMEbatchargs=""; fi
if [ -z "${do_runoff+x}" ]; then do_runoff=false; fi
if [ -z "${keep_land_restarts+x}" ]; then keep_land_restarts=true; fi
if [ -z "${perturb_namelist+x}" ]; then perturb_namelist=""; fi
if [ -z "${predict_docn+x}" ]; then predict_docn=false; fi
if [ -z "${archive_inic+x}" ]; then archive_inic=false; fi
if [ -z "${compress_history_nc+x}" ]; then compress_history_nc=true; fi
if [ -z "${add_vortex+x}" ]; then add_vortex=false; fi
if [ -z "${vortex_namelist+x}" ]; then vortex_namelist=""; fi
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
if [ -z "${perturb_namelist+x}" ]; then perturb_namelist=""; fi
if [ -z "${land_spinup+x}" ]; then land_spinup=false; fi
if [ -z "${FILTERWALLCLOCK+x}" ]; then FILTERWALLCLOCK="00:29:00"; fi
if [ -z "${FILTERQUEUE+x}" ]; then FILTERQUEUE="batch"; fi
if [ -z "${RUNWALLCLOCK+x}" ]; then RUNWALLCLOCK="12:00:00"; fi
if [ -z "${RUNQUEUE+x}" ]; then RUNQUEUE="regular"; fi
if [ -z "${RUNPRIORITY+x}" ]; then RUNPRIORITY=""; fi
if [ -z "${use_nsplit+x}" ]; then use_nsplit="true"; fi
if [ -z "${cime_coupler+x}" ]; then cime_coupler="mct"; fi
if [ -z "${nclPlotWeights+x}" ]; then nclPlotWeights="NULL"; fi
if [ -z "${landrawdir+x}" ]; then landrawdir="NULL"; fi
if [ -z "${sendplots+x}" ]; then sendplots=false; fi
if [ -z "${dotracking+x}" ]; then dotracking=false; fi
if [ -z "${m2m_gridfile+x}" ]; then m2m_gridfile=""; fi
if [ -z "${m2m_remap_file+x}" ]; then m2m_remap_file=""; fi

echo "DO_PYTHON set to $DO_PYTHON"
if [ "$DO_PYTHON" = true ]; then
  atm_to_cam_path=${SCRIPTPATH}/py_atm_to_cam
  sst_to_cam_path=${SCRIPTPATH}/py_sst_to_cam
fi

### Set correct E3SM/CESM split
if [ -z "${modelSystem+x}" ]; then modelSystem=0; fi
if [ $modelSystem -eq 0 ]; then
  echo "Using CESM"
  atmName="cam"
  lndName="clm"
  lndSpecialName="clm2"
  rofName="mosart"
  rofSpecialName="mosart"
elif [ $modelSystem -eq 1 ]; then
  echo "Using E3SM"
  atmName="eam"
  lndName="elm"
  lndSpecialName="elm"
  rofName="mosart"
  rofSpecialName="mosart"
else
  echo "Unknown modeling system set for modelSystem: $modelSystem"
  exit 1
fi

# Are we running with frankengrid?
do_frankengrid=false
if [ -n "${regional_src+x}" ]; then do_frankengrid=true ; fi
echo "do_frankengrid set to: $do_frankengrid"

# Figure out where to archive
if [ -z ${ARCHIVEDIR+x} ] || [[ -z "${ARCHIVEDIR// }" ]] ; then
  ARCHIVEDIR=${outputdir}/
else
  ARCHIVEDIR=${ARCHIVEDIR}/${casename}/
  mkdir -v -p ${ARCHIVEDIR}
fi
echo "Files will be archived in ${ARCHIVEDIR}/YYYYMMDDHH/"

# Update SST filename based on docnres
sstFileIC=${sstFileIC/DOCNRES/$docnres}
echo "Actual sstFileIC: ${sstFileIC}"

### ERROR CHECKING BLOCK! #########################################################

# Exit if add_perturbs is turned on but no namelist is passed in with perturbation config settings
if [ "$add_perturbs" = true ] && { [ -z "$perturb_namelist" ] || [ ! -f "$perturb_namelist" ]; } ; then
  echo "add_perturbs is true but can't find namelist: "$perturb_namelist
  exit 1
fi

# Exit if add_perturbs is turned on but no namelist is passed in with perturbation config settings
if [ "$add_vortex" = true ] && { [ -z "$vortex_namelist" ] || [ ! -f "$vortex_namelist" ]; } ; then
  echo "add_vortex is true but can't find namelist: "$vortex_namelist
  exit 1
fi

# Check for deprecated namelist options
if [ -z ${anl2mdlWeights+x} ] && [ ${gfs2seWeights+x} ] ; then
  echo "WARNING: Setting anl2mdlWeights to ${gfs2seWeights}"
  echo "WARNING: This is deprecated and will be removed in the future! To fix, change 'gfs2seWeights' to 'anl2mdlWeights' in ${NAMELISTFILE}"
  anl2mdlWeights=$gfs2seWeights
fi

# Check if ncl exists
if ! type ncl &> /dev/null ; then
  echo "ERROR: ncl does not exist. Make sure ncl is in your path when betacast is invoked"
  exit 1
fi

# Check if ncks exists for compression
ncks_exists=true
if ! type ncks &> /dev/null ; then
  #echo "ERROR: ncks does not exist. Make sure ncks is in your path when betacast is invoked"
  #exit 1
  echo "WARNING: ncks does not exist, cannot compress. Setting compress_history_nc to 0 (false)"
  echo "WARNING: if you'd like remedy, make sure ncks is in your path when betacast.sh is invoked"
  compress_history_nc=0
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

# Check Python version for 3+
python -c 'import sys; exit(sys.version_info.major != 3)' && echo "CHECK_PYTHON: Python 3 found" || { echo "CHECK_PYTHON: Please install Python 3+"; exit 25; }

# Check Python package dependencies
if [ "$do_frankengrid" = true ] ; then
  check_python_dependency numpy
  check_python_dependency netCDF4
  check_python_dependency sklearn
fi

# Adjust bools (for backwards compatibility, 0 = false and 1 = true)
bools_to_check=("islive" "debug" "doFilter" "filterOnly" "do_runoff" "keep_land_restarts"
       "predict_docn" "archive_inic" "compress_history_nc" "override_rest_check"
       "tararchivedir" "save_nudging_files" "add_vortex")
for bool_to_check in ${bools_to_check[@]}; do
  check_bool $bool_to_check ${!bool_to_check}
done

if [ $override_rest_check = false ]; then
  echo "Checking for SourceMods permiting additional restart writes for land model"
  echo "This check can be ignored with override_rest_check = true in the namelist."
  exit_files_no_exist $path_to_case/SourceMods/src.${lndName}/lnd_comp_mct.F90 $path_to_case/SourceMods/src.${lndName}/lnd_comp_nuopc.F90
  if [ $do_runoff = true ]; then
    exit_files_no_exist $path_to_case/SourceMods/src.${rofName}/rof_comp_mct.F90 $path_to_case/SourceMods/src.${rofName}/rof_comp_nuopc.F90
  fi
fi

### Check required variables
betacast_required_vars=("casename" "atmDataType" "sstDataType" "numLevels" "numdays"
        "anl2mdlWeights" "PROJECTID" "DTIME" "FINERES" "USERSTAB")
check_required_vars "${betacast_required_vars[@]}"

###################################################################################

# do some stability calcs
# if USERSTAB is 0, use internal calcs.
# if USERSTAB is <0, use se_nsplit=-1
# If USERSTAB >0, use the value in seconds for dt_dyn and calculate nsplit accordingly from DTIME
USERSTABTF=$(python -c "print('TRUE' if ${USERSTAB} > 0 else 'FALSE')")
if [ ${USERSTABTF} == 'FALSE' ] ; then
  if [ $(python -c "print('TRUE' if ${USERSTAB} < -0.001 else 'FALSE')") == 'FALSE' ]; then
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

# Get the current time
currtime=$(date -u +%H%M)

if [ $islive = true ] ; then    # Find most recent GFS forecast
  ## Here we get two digit strings for UTC time for month, day, year
  ## We also get current time in hoursminutes (because the GFS output lags by 3.5 hours)
  monthstr=$(date -u +%m)
  daystr=$(date -u +%d)
  yearstr=$(date -u +%Y)

  ## Use currtime to figure out what is the latest cycle we have access to
  if [ $currtime -lt 0328 ] ; then
    echo "12Z cycle"
    monthstr=$(date --date="yesterday" -u +%m)
    daystr=$(date --date="yesterday" -u +%d)
    yearstr=$(date --date="yesterday" -u +%Y)
    cyclestr=12
  elif [ $currtime -lt 0928 ] ; then
    echo "00Z cycle"
    cyclestr=00
  elif [ $currtime -lt 1528 ] ; then
    echo "00Z cycle"
    cyclestr=00
  elif [ $currtime -lt 2128 ] ; then
    echo "12Z cycle"
    cyclestr=12
  elif [ $currtime -ge 2128 ] ; then
    echo "12Z cycle"
    cyclestr=12
  else
    echo "Can't figure out start time"
    exit 1
  fi

else     # if not live, draw from head of dates.txt file

  # Figure out if dates.CASE.txt exists and where it is located
  datesbase=dates.${casename}.txt
  if [ -f "./dates/${datesbase}" ]; then
    # Preferred location is dates subdir as of 4/13/2022
    echo "$datesbase exists in dates subdirectory"
    datesfile=./dates/${datesbase}
  else
    if [ -f "./${datesbase}" ]; then
      echo "WARNING! $datesbase isn't in dates subdir but exists in home dir!"
      echo "This is allowed for backwards compatibility but is less organized!"
      echo "You should create a subdir called 'dates' and put the dates.CASE.txt file there!"
      datesfile=${datesbase}
    else
      if [[ ${datestemplate+x} && -f ./dates/${datestemplate} ]]; then
        echo "Didn't find case dates, but you specified a datestemplate in the namelist..."
        echo "So I'm copying $datestemplate to this case and using that..."
        cp -v ./dates/${datestemplate} ./dates/${datesbase}
        datesfile=./dates/${datesbase}
      else
        echo "Uh, oh. Can't find a dates file OR template AND run isn't live. Exiting..."
        exit 1
      fi
    fi
  fi

  echo "Using dates in: "${datesfile}
  longdate=$(get_top_line_from_dates "${datesfile}")
  echo "Getting parsed time from $longdate"
  parse_YYYYMMDDHH $longdate
  echo "From datesfile, read in: "$yearstr' '$monthstr' '$daystr' '$cyclestr'Z'
fi

## Figure out the seconds which correspond to the cycle and zero pad if neces
get_cyclestrsec "$cyclestr"

## Figure out what the SE start time will be after filter
## These values are *only* used if do_filter is true
if [ $numHoursSEStart -lt 6 ] ; then
  let se_cyclestr=$cyclestr+03
  while [ ${#se_cyclestr} -lt 2 ];
  do
    se_cyclestr="0"$se_cyclestr
  done
  se_monthstr=$monthstr
  se_daystr=$daystr
  se_yearstr=$yearstr
  let se_cyclestrsec=$((10#$se_cyclestr))*3600
  while [ ${#se_cyclestrsec} -lt 5 ];
  do
    se_cyclestrsec="0"$se_cyclestrsec
  done
else
  echo "SE forecast lead time too long, 18Z cycle causes trouble"
  echo "Not supported."
  exit 1
fi

# Get time for pulling SST
getSSTtime $islive $currtime $monthstr $daystr $yearstr $cyclestr

yestmonthstr=$(date --date="yesterday" -u +%m)
yestdaystr=$(date --date="yesterday" -u +%d)
yestyearstr=$(date --date="yesterday" -u +%Y)

echo "We are using $yearstr $monthstr $daystr $cyclestr Z ($cyclestrsec seconds) for ATM init. data"
echo "We are using $sstyearstr $sstmonthstr $sstdaystr $sstcyclestr Z for SST init. data"
if [ "$doFilter" = true ] ; then
  echo "Filter: True model init. will occur at $se_yearstr $se_monthstr $se_daystr $se_cyclestr Z ($se_cyclestrsec seconds)"
else
  echo "No filter: True model init. will occur at $se_yearstr $se_monthstr $se_daystr $cyclestr Z ($cyclestrsec seconds)"
fi

if [ "$runmodel" = true ] ; then

############################### GET DYCORE INFO ###############################

cd $path_to_case
DYCORE=$(./xmlquery CAM_DYCORE | sed 's/^[^\:]\+\://' | xargs)
if [ -z "$DYCORE" ]; then
  echo "Couldn't automagically figure out dycore, assuming SE/HOMME"
  DYCORE="se"
fi
echo "DYCORE: "$DYCORE

if [ $debug = false ] ; then
############################### GET ATM DATA ###############################

  # The keys are $atmDataType
  declare -A atm_data_sources=(
    ["1"]="GFS"
    ["2"]="ERAI"
    ["3"]="CFSR"
    ["4"]="ERA5"
    ["9"]="CAM"
  )
  declare -A atm_data_glob_anl=(
    ["1"]="gfs_0.25x0.25"
    ["2"]=""
    ["3"]="gfs_0.50x0.50"
    ["4"]="era5_0.25x0.25"
    ["9"]=""
  )
  declare -A atm_file_paths=(
    ["1"]="${gfs_files_path}/gfs_atm_${yearstr}${monthstr}${daystr}${cyclestr}.grib2"
    ["2"]="${era_files_path}/ERA-Int_${yearstr}${monthstr}${daystr}${cyclestr}.nc"
    ["3"]="${gfs_files_path}/cfsr_atm_${yearstr}${monthstr}${daystr}${cyclestr}.grib2"
    ["4"]="${era_files_path}/ERA5_${yearstr}${monthstr}${daystr}${cyclestr}.nc"
    ["9"]=""
  )

  # Initialize any "global" variables needed when getting atmospheric data
  RDADIR="" # Init to empty, but fill in if RDA available later
  ERA5RDA=0 # Set whether or not ERA5 is local (0 = local, 1 = RDA)

  case $atmDataType in
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

  # If ERA5RDA flag toggled, set value w/ key to RDA data
  if [ $ERA5RDA -eq 1 ] ; then
    atm_data_sources["4"]="ERA5RDA"
    atm_file_paths["4"]="${RDADIR}/e5.oper.invariant/197901/e5.oper.invariant.128_129_z.ll025sc.1979010100_1979010100.nc"
  fi

############################### GET SST / NCL ###############################

  case $sstDataType in
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
  if [ ${sstDataType} -ne 9 ] ; then
    # Switch bash bool to int for NCL input
    if [ $predict_docn = true ]; then INT_PREDICT_DOCN=1; else INT_PREDICT_DOCN=0; fi

    cd $sst_to_cam_path ; echo "cd'ing to interpolation directory: $sst_to_cam_path"

    sst_domain_file=${SCRIPTPATH}/grids/domains/domain.ocn.${docnres}.nc
    sst_scrip_file=${SCRIPTPATH}/grids/domains/scrip.ocn.${docnres}.nc
    sst_ESMF_file=${SCRIPTPATH}/grids/domains/ESMF.ocn.${docnres}.nc

    # check if domain or SCRIP exist, if one is missing create both domain and scrip
    if [ ! -f "$sst_domain_file" ] || [ ! -f "$sst_scrip_file" ]; then
      echo "Creating SST domain file for: ${docnres}"
      if [ "$DO_PYTHON" = true ]; then
        echo "Oops" ; exit
      else
        set +e
        (set -x; ncl gen-sst-domain.ncl 'inputres="'${docnres}'"' ) ; exit_status=$?
        check_ncl_exit "gen-sst-domain.ncl" $exit_status
        set -e
      fi
      compress_single_file "${sst_scrip_file}"
    fi
    # if the CIME coupler is nuopc, we need to generate a scrip file
    if [ "${cime_coupler}" == "nuopc" ] && [ ! -f "$sst_ESMF_file" ]; then
      ESMF_Scrip2Unstruct ${sst_scrip_file} ${sst_ESMF_file} 0 ; rm -fv PET0.ESMF_LogFile
      compress_single_file "${sst_ESMF_file}"
    fi

    # Now generate the SST/ice datastream
    if [ "$DO_PYTHON" = true ]; then
      (set -x; python sst_to_cam.py \
          --initdate "${yearstr}${monthstr}${daystr}${cyclestr}" \
          --predict_docn ${INT_PREDICT_DOCN} \
          --inputres "${docnres}" \
          --datasource "${SSTTYPE}" \
          --sstDataFile "${sst_files_path}/${sstFile}" \
          --iceDataFile "${sst_files_path}/${iceFile}" \
          --SST_write_file "${sstFileIC}" \
          --smooth_ice \
          --smooth_iter 3 \
          --verbose
      )
    else
      set +e
      (set -x; ncl sst_interp.ncl \
          'initdate="'${yearstr}${monthstr}${daystr}${cyclestr}'"' \
          predict_docn=${INT_PREDICT_DOCN} \
          'inputres="'${docnres}'"' \
          'datasource="'${SSTTYPE}'"' \
          'sstDataFile = "'${sst_files_path}/${sstFile}'"' \
          'iceDataFile = "'${sst_files_path}/${iceFile}'"' \
          'SST_write_file = "'${sstFileIC}'"' \
      ) ; exit_status=$?
      check_ncl_exit "sst_interp.ncl" $exit_status
      set -e # Turn error checking back on
    fi
  fi

  ############################### ATM NCL ###############################

  cd $atm_to_cam_path ; echo "cd'ing to interpolation directory: $atm_to_cam_path"

  # Figure out which anl2mdlWeights we want to use. If the user gave us one
  # we will just use that, otherwise we'll hope they gave us modelgridfile (SCRIP)
  # and we will use that if not generated. Once it's been generated, keep reusing until
  # purged from the mapping directory
  if [[ -z "${anl2mdlWeights}" || ! -e "${anl2mdlWeights}" ]]; then
    echo "User did not explicitly specify anl2mdlWeights, trying to generate from SCRIP grid"
    if [ ! -f ${modelgridfile} ]; then
      echo "modelgridfile --> ${modelgridfile} does not exist, exiting"
      echo "specify this as a SCRIP file in the namelist or anl2mdlWeights"
      exit 19
    fi
    # Get name without suffix
    modelgridshortname=$(basename "${modelgridfile%.*}")
    # If we are doing m2m, specify the target as ERA5, otherwise let user decide
    if [[ "$atmDataType" -eq 9 ]]; then
      RLLSOURCEGRID="era5_0.25x0.25"
    else
      RLLSOURCEGRID=${atm_data_glob_anl[$atmDataType]}
    fi
    # Define new anl2mdlWeights
    anl2mdlWeights=${mapping_files_path}/map_${RLLSOURCEGRID}_TO_${modelgridshortname}_patc.nc
    # Check if anl2mdlWeights exist or not, if not try to generate them
    if [ ! -f ${anl2mdlWeights} ]; then
      echo "Writing anl2mdlWeights --> ${anl2mdlWeights}"
      if [ "$DO_PYTHON" = true ]; then
        (set -x; python ../py_remapping/gen_analysis_to_model_wgt_file.py \
          --ANLGRID "${RLLSOURCEGRID}" \
          --DSTGRIDNAME "${modelgridshortname}" \
          --DSTGRIDFILE "${modelgridfile}" \
          --ANLGRIDPATH "../grids/anl_scrip/" \
          --WGTFILEDIR "${mapping_files_path}"
        )
      else
        set +e
        (set -x; ncl ../remapping/gen_analysis_to_model_wgt_file.ncl \
          'ANLGRID="'${RLLSOURCEGRID}'"' \
          'ANLGRIDPATH="../grids/anl_scrip/"' \
          'DSTGRIDNAME="'${modelgridshortname}'"' \
          'DSTGRIDFILE="'${modelgridfile}'"' \
          'WGTFILEDIR="'${mapping_files_path}'"' \
        ) ; exit_status=$?
        check_ncl_exit "gen_analysis_to_model_wgt_file.ncl" $exit_status
        set -e
      fi
    else
      echo "Betacast-generated anl2mdlWeights --> ${anl2mdlWeights} already exists, using those!"
    fi
  else
    echo "User has provided anl2mdlWeights --> ${anl2mdlWeights}, using those!"
  fi # check if anl2mdlWeights is passed in

  if [[ "$atmDataType" -eq 9 ]]; then
    if [[ -z "${m2m_remap_file}" || ! -e "${m2m_remap_file}" ]]; then
      echo "User did not explicitly specify m2m_remap_file, trying to generate from SCRIP grid"
      if [ ! -f ${m2m_gridfile} ]; then
        echo "m2m_gridfile --> ${m2m_gridfile} does not exist, exiting"
        echo "specify this as a SCRIP file in the namelist or m2m_remap_file"
        exit 19
      fi
      # Get name without suffix
      m2mgridshortname=$(basename "${m2m_gridfile%.*}")

      m2m_remap_file=${mapping_files_path}/map_${m2mgridshortname}_TO_era5_0.25x0.25_patc.nc

      # Check if m2m_remap_file exist or not, if not try to generate them
      if [ ! -f ${m2m_remap_file} ]; then
        echo "Writing m2m_remap_file --> ${m2m_remap_file}"
        set +e
        (set -x; ncl ../remapping/gen_analysis_to_model_wgt_file.ncl \
          'ANLGRID="era5_0.25x0.25"' \
          'ANLGRIDPATH="../grids/anl_scrip/"' \
          'DSTGRIDNAME="'${m2mgridshortname}'"' \
          'DSTGRIDFILE="'${m2m_gridfile}'"' \
          'WGTFILEDIR="'${mapping_files_path}'"' \
          FLIP_MODEL_AND_ANALYSIS=True \
        ) ; exit_status=$?
        check_ncl_exit "gen_analysis_to_model_wgt_file.ncl" $exit_status
        set -e
      else
        echo "Betacast-generated m2m_remap_file --> ${m2m_remap_file} already exists, using those!"
      fi
    else
      echo "User has provided m2m_remap_file --> ${m2m_remap_file}, using those!"
    fi
  fi

  if [[ "$atmDataType" -eq 9 ]]; then
    # If do_model2model is true, then we need to find the file
    if [[ -d "$m2m_parent_source" ]]; then
      echo "m2m_parent_source ($m2m_parent_source) is provided as a dir of nc files."
      if [ "$DO_PYTHON" = true ]; then
        (set -x; python find-time-file.py \
          --DIR "${m2m_parent_source}" \
          --YYYYMMDDHH ${yearstr}${monthstr}${daystr}${cyclestr} \
          --UQSTR "${uniqtime}"
        )
      else
        set +e
        (set -x; ncl find-time-file.ncl \
          'DIR="'${m2m_parent_source}'"' \
          YYYYMMDDHH=${yearstr}${monthstr}${daystr}${cyclestr} \
          'UQSTR="'${uniqtime}'"' \
        ) ; exit_status=$?
        check_ncl_exit "find-time-file.ncl" $exit_status
        set -e
      fi
      while IFS= read -r line || [[ -n "$line" ]]; do
        atm_file_paths["9"]=$line
        break # Exit after reading the first line
      done < m2mfile.$uniqtime
      [[ ! -f "${atm_file_paths["9"]}" ]] && { echo "File does not exist."; exit 1; }
    elif [[ -f "$m2m_parent_source" ]]; then
      echo "m2m_parent_source ($m2m_parent_source) is provided as a file."
      atm_file_paths["9"]=$m2m_parent_source
    else
      echo "m2m_parent_source ($m2m_parent_source) is not a file or directory. Exiting."
      exit 1
    fi
  fi

  echo "Doing atm_to_cam"

  if [ "$DO_PYTHON" = true ]; then
    (set -x ; python atm_to_cam.py \
      --datasource "${atm_data_sources[$atmDataType]}" \
      --numlevels ${numLevels} \
      --YYYYMMDDHH ${yearstr}${monthstr}${daystr}${cyclestr} \
      --data_filename "${atm_file_paths[$atmDataType]}" \
      --wgt_filename "${anl2mdlWeights}" \
      --dycore "${DYCORE}" \
      --add_cloud_vars \
      --RDADIR "${RDADIR}" \
      --adjust_config "${adjust_flags-}" \
      --model_topo_file "${adjust_topo-}" \
      --mod_remap_file "${m2m_remap_file-}" \
      --mod_in_topo "${m2m_topo_in-}" \
      --model_topo_file "${adjust_topo-}" \
      --se_inic "${sePreFilterIC}" )
  else
    set +e #Need to turn off error checking b/c NCL returns 0 even if fatal
    (set -x; ncl -n atm_to_cam.ncl \
        'datasource="'${atm_data_sources[$atmDataType]}'"' \
        numlevels=${numLevels} \
        YYYYMMDDHH=${yearstr}${monthstr}${daystr}${cyclestr} \
        'dycore="'${DYCORE}'"' \
        'data_filename="'${atm_file_paths[$atmDataType]}'"' \
        'RDADIR="'${RDADIR}'"' \
        'wgt_filename="'${anl2mdlWeights}'"' \
        'model_topo_file="'${adjust_topo-}'"' \
        'mod_remap_file="'${m2m_remap_file-}'"' \
        'mod_in_topo="'${m2m_topo_in-}'"' \
        'adjust_config="'${adjust_flags-}'"' \
        'se_inic = "'${sePreFilterIC}'"'
    ) ; exit_status=$?
    check_ncl_exit "atm_to_cam.ncl" $exit_status
    set -e
  fi

  if [ "$do_frankengrid" = true ] ; then

    # Fill in "templated" information from regional_src
    regional_src=${regional_src/YYYY/$yearstr}
    regional_src=${regional_src/MM/$monthstr}
    regional_src=${regional_src/DD/$daystr}
    regional_src=${regional_src/HH/$cyclestr}

    echo "Doing Frankengrid with $regional_src"

    TMPWGTFILE="./map_hwrf_storm_TO_modelgrid_patc.nc"

    if [ "$DO_PYTHON" = true ]; then
      (set -x; python ../py_remapping/gen_reglatlon_SCRIP.py \
        --dstGridName "hwrf_storm_scrip.nc" \
        --dstDir "./" \
        --srcfilename "${regional_src}"
      )
    else
      set +e # Turn error checking off for NCL
      echo "Generating a temporary SCRIP file for HWRF"
      (set -x; ncl ../remapping/gen_reglatlon_SCRIP.ncl \
        'DSTGRIDNAME="hwrf_storm_scrip.nc"' \
        'DSTDIR="./"' \
        'SRCFILENAME="'${regional_src}'"'
      ) ; exit_status=$?
      check_ncl_exit "gen_reglatlon_SCRIP.ncl" $exit_status
    fi

    if [ "$DO_PYTHON" = true ]; then
      (set -x; python ../py_remapping/gen_analysis_to_model_wgt_file.py \
        --ANLGRID "hwrf_storm" \
        --DSTGRIDNAME "modelgrid" \
        --DSTGRIDFILE "${model_scrip}" \
        --ANLGRIDPATH "./" \
        --WGTFILEDIR "./"
      )
    else
      echo "Generating a temporary map file for HWRF"
      (set -x; ncl ../remapping/gen_analysis_to_model_wgt_file.ncl \
        'ANLGRID="hwrf_storm"' \
        'DSTGRIDNAME="modelgrid"' \
        'ANLGRIDPATH="./"' \
        'WGTFILEDIR="./"' \
        'DSTGRIDFILE="'${model_scrip}'"'
      ) ; exit_status=$?
      check_ncl_exit "gen_analysis_to_model_wgt_file.ncl" $exit_status
    fi

    if [ "$DO_PYTHON" = true ]; then
      (set -x ; python atm_to_cam.py \
        --datasource "HWRF" \
        --numlevels ${numLevels} \
        --YYYYMMDDHH ${yearstr}${monthstr}${daystr}${cyclestr} \
        --data_filename "${regional_src}" \
        --wgt_filename "${TMPWGTFILE}" \
        --dycore "${DYCORE}" \
        --add_cloud_vars \
        --adjust_config "" \
        --se_inic "${sePreFilterIC}_reg.nc"
      )
    else
      echo "Generating regional Frankengrid for HWRF"
      (set -x; ncl -n atm_to_cam.ncl 'datasource="HWRF"' \
        'dycore="'${DYCORE}'"' \
        numlevels=${numLevels} \
        YYYYMMDDHH=${yearstr}${monthstr}${daystr}${cyclestr} \
        'data_filename = "'${regional_src}'"' \
        'wgt_filename="'${TMPWGTFILE}'"' \
        'adjust_config=""' \
        'se_inic = "'${sePreFilterIC}_reg.nc'"'
      ) ; exit_status=$?
      check_ncl_exit "atm_to_cam.ncl" $exit_status
      set -e # Turn error checking back on
    fi

    echo "Overlay regional file on top of basefile"
    cp -v ${sePreFilterIC} ${sePreFilterIC}_base.nc
    (set -x; python overlay.py "${sePreFilterIC}" "${sePreFilterIC}_reg.nc" --maxLev 80. )

    echo "Cleaning up temporary ESMF files"
    rm -v $TMPWGTFILE
    rm -v hwrf_storm_scrip.nc
    rm -v ${sePreFilterIC}_reg.nc

  fi

fi #End debug if statement

############################### #### ###############################
##### ADD OR REMOVE VORTEX

if [ "$add_vortex" = true ] ; then
  cd $atm_to_cam_path/tcseed
  set +e
  echo "Adding or removing a TC from initial condition based on ${vortex_namelist}"

  echo "... finding fill parameters"
  if [ "$DO_PYTHON" = true ]; then
    (set -x; python find-tc-fill-params.py \
        --inic_file "${sePreFilterIC}" \
        --vortex_namelist ${vortex_namelist}
    )
  else
    (set -x; ncl -n find-tc-fill-params.ncl 'inic_file= "'${sePreFilterIC}'"' 'pthi = "'${vortex_namelist}'"' ) ; exit_status=$?
    check_ncl_exit "find-tc-fill-params.ncl" $exit_status
  fi

  echo "... seeding or unseeding TC"
  if [ "$DO_PYTHON" = true ]; then
    (set -x; python py-seed-tc-in-ncdata.py \
        --inic_file "${sePreFilterIC}" \
        --vortex_namelist ${vortex_namelist}
    )
  else
    (set -x; ncl -n seed-tc-in-ncdata.ncl 'seedfile = "'${sePreFilterIC}'"' 'pthi = "'${vortex_namelist}'"' ) ; exit_status=$?
    check_ncl_exit "seed-tc-in-ncdata.ncl" $exit_status
  fi

  set -e
fi

############################### #### ###############################
##### ADD WHITE NOISE PERTURBATIONS

if [ "$add_noise" = true ] ; then
  set +e
  echo "Adding white noise to initial condition"
  cd $atm_to_cam_path
  (set -x; ncl -n perturb_white_noise.ncl 'basFileName = "'${sePreFilterIC}'"' ) ; exit_status=$?
  check_ncl_exit "perturb_white_noise.ncl" $exit_status
  set -e
fi

############################### #### ###############################
##### ADD PERTURBATIONS

if [ "$add_perturbs" = true ] ; then
  echo "Adding perturbations"

  cd $atm_to_cam_path/perturb

  set +e

  ## Add perturbations to SST file
  sstFileIC_WPERT=${sstFileIC}_PERT.nc
  (set -x; ncl -n add_perturbations_to_sst.ncl 'BEFOREPERTFILE="'${sstFileIC}'"' \
     'AFTERPERTFILE = "'${sstFileIC_WPERT}'"' \
     'pthi="'${perturb_namelist}'"'
  ) ; exit_status=$?
  check_ncl_exit "add_perturbations_to_sst.ncl" $exit_status
  echo "SST perturbations added successfully"

  ## Add perturbations to ATM file
  sePreFilterIC_WPERT=${sePreFilterIC}_PERT.nc
  (set -x; ncl -n add_perturbations_to_cam.ncl 'BEFOREPERTFILE="'${sePreFilterIC}'"'  \
     'AFTERPERTFILE = "'${sePreFilterIC_WPERT}'"' \
     'gridfile = "'${modelgridfile}'"' \
     'MAPFILEPATH = "'${mapping_files_path}'"' \
     'pthi="'${perturb_namelist}'"'
  ) ; exit_status=$?
  check_ncl_exit "add_perturbations_to_cam.ncl" $exit_status
  echo "ATM NCL completed successfully"

  set -e # Turn error checking back on

  # move temp perturb files to overwrite unperturbed files
  mv ${sstFileIC_WPERT} ${sstFileIC}
  mv ${sePreFilterIC_WPERT} ${sePreFilterIC}
fi

cd $path_to_case

############################### CISM SETUP ###############################

if [ -f user_nl_cism ]; then
  sed -i '/dt_count/d' user_nl_cism
  echo "dt_count = 8" >> user_nl_cism
fi

############################### CLM SETUP ###############################

echo "................. configuring land and/or runoff"

## TEMPORARY CLM -> LAND FIX
## Check if clmstart exists but landstart doesn't, move if that is the case
if [[ ! -d ${landdir} ]] && [[ -d ${path_to_rundir}/clmstart/ ]] ; then
  echo "Moving ${path_to_rundir}/clmstart/ to ${landdir}"
  mv -v ${path_to_rundir}/clmstart/ ${landdir}
else
  echo "No need to modify land directory, either appropriately exists or will be created later"
fi
## END TEMP FIX

echo "Setting input land dataset"
# Clean up file to delete any special interp lines that may be needed later (but aren't needed for native init)
sed -i '/init_interp_fill_missing_with_natveg/d' user_nl_${lndName}
sed -i '/use_init_interp/d' user_nl_${lndName}

# Create a temp directory for now
RUNTMPDIR=$path_to_nc_files/tmp/
[ -z "$RUNTMPDIR" ] && { echo "RUNTMPDIR is not set. Exiting."; exit 1; }
# Create directory and make sure it was actually generated (to prevent user from not having write perms)
mkdir -vp $RUNTMPDIR
[ -d "$RUNTMPDIR" ] || { echo "Failed to create $RUNTMPDIR. Exiting."; exit 1; }
# Delete files in RUNTMPDIR
find "$RUNTMPDIR" -type f -exec rm -v {} +

# We want to check ${landdir} for land restart file with exact date match. If so, use that.
landrestartfile=$(find "${landdir}" -maxdepth 1 \( -type f -o -type l \) -name "${casename}.${lndSpecialName}.r.${yearstr}-${monthstr}-${daystr}-${cyclestrsec}.*" | head -n 1)

# Check to see if file exists on native SE land grid
if [ -n "${landrestartfile}" ] && [ -f "${landrestartfile}" ]; then
  echo "File ${landrestartfile} exists at exact time"
else
  # If we don't have an exact match from above, trying 00Z first, then check landrawdir
  echo "File ${landrestartfile} does not exist at exact time"
  landrestartfile=$(find "${landdir}" -maxdepth 1 -type f -name "${casename}.${lndSpecialName}.r.${yearstr}-${monthstr}-${daystr}-00000.*" | head -n 1)
  if [ -n "${landrestartfile}" ] && [ -f "${landrestartfile}" ]; then
    echo "File ${landrestartfile} exists at 00Z"
  else
    rawlandrestartfile=$(find "${landrawdir}" -maxdepth 1 -type f -name "*.${lndSpecialName}.r.${yearstr}-${monthstr}-${daystr}-${cyclestrsec}.*" | head -n 1)
    echo "rawlandrestartfile: ${rawlandrestartfile}"
    if [ -n "${rawlandrestartfile}" ] && [ -f "${rawlandrestartfile}" ]; then
      landrestartfile=${rawlandrestartfile}
    else
      echo "No LND restart file exists, setting landrestartfile to empty string."
      landrestartfile=
    fi
  fi
fi
echo "landrestartfile: ${landrestartfile}"

if [ -f "${landrestartfile}" ] && [[ "${landrestartfile}" != *.nc ]]; then
  echo "${landrestartfile} was found, but does not have an *.nc extension. Assuming compressed. Copying..."
  cp -v ${landrestartfile} $RUNTMPDIR
  try_uncompress $RUNTMPDIR/$(basename "${landrestartfile}")
  landrestartfile=$(find "${RUNTMPDIR}" -maxdepth 1 -type f -name "*.${lndSpecialName}.r.${yearstr}-${monthstr}-${daystr}-*.nc" | head -n 1)
  echo "Updated landrestartfile: ${landrestartfile}"
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
if [ $modelSystem -eq 0 ]; then   # CLM/CTSM
  sed -i '/check_finidat_pct_consistency/d' user_nl_${lndName}
  sed -i '/check_finidat_year_consistency/d' user_nl_${lndName}
  echo "check_finidat_pct_consistency = .false." >> user_nl_${lndName}
  echo "check_finidat_year_consistency = .false." >> user_nl_${lndName}
elif [ $modelSystem -eq 1 ]; then   # ELM
  sed -i '/check_finidat_fsurdat_consistency/d' user_nl_${lndName}
  echo "check_finidat_fsurdat_consistency = .false." >> user_nl_${lndName}
  # 2/25/24 CMZ added since ELM doesn't have use_init_interp support for rawlandrestartfile
  if [ -n "${rawlandrestartfile-}" ]; then # if rawlandrestartfile is SET *and* not empty...
    echo "WARNING USER_NL: We used a rawlandrestartfile, but ELM doesn't support interpolation, removing relevant lines"
    echo "WARNING USER_NL: If your model fails, check to make sure the rawlandrestartfile supports your target land grid"
    sed -i '/use_init_interp/d' user_nl_${lndName}
    sed -i '/init_interp_fill_missing_with_natveg/d' user_nl_${lndName}
  fi
else
  echo "Unknown modeling system set for modelSystem: $modelSystem"
  exit 1
fi

############################### ROF SETUP ###############################

if [ $do_runoff = true ]; then

  echo "USER_NL: Setting input ${rofName} dataset"

  # Delete any existing input data
  sed -i '/.*finidat_rtm/d' user_nl_${rofName}

  # We want to check ${landdir} for land restart files. If so, use those.
  rofrestartfile=$(find "${landdir}" -maxdepth 1 \( -type f -o -type l \) -name "${casename}.${rofSpecialName}.r.${yearstr}-${monthstr}-${daystr}-${cyclestrsec}.*" | head -n 1)
  echo $rofrestartfile

  # Check to see if file exists on native SE land grid
  if [ -n "${rofrestartfile}" ] && [ -f "${rofrestartfile}" ]; then
     echo "USER_NL: ${rofName} file exists at exact time"
  else
     echo "USER_NL: ${rofName} file does not exist at exact time"
     rofrestartfile=$(find "${landdir}" -maxdepth 1 -type f -name "${casename}.${rofSpecialName}.r.${yearstr}-${monthstr}-${daystr}-00000.*" | head -n 1)
     if [ -n "${rofrestartfile}" ] && [ -f "${rofrestartfile}" ]; then
       echo "USER_NL: ${rofName} file exists at 00Z"
    else
      rawrofrestartfile=$(find "${landrawdir}" -maxdepth 1 -type f -name "*.${rofSpecialName}.r.${yearstr}-${monthstr}-${daystr}-${cyclestrsec}.*" | head -n 1)
      echo "rawrofrestartfile: ${rawrofrestartfile}"
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
    cp -v ${rofrestartfile} $RUNTMPDIR
    try_uncompress $RUNTMPDIR/$(basename "${rofrestartfile}")
    rofrestartfile=$(find "${RUNTMPDIR}" -maxdepth 1 -type f -name "*.${rofSpecialName}.r.${yearstr}-${monthstr}-${daystr}-*.nc" | head -n 1)
    echo "Updated rofrestartfile: ${rofrestartfile}"
  fi

  ## Now modify user_nl_${rofName}
  if [ ${rofrestartfile} ] ; then
    echo "USER_NL: Adding ${rofrestartfile} to user_nl_${rofName}"
    echo "finidat_rtm='${rofrestartfile}'" >> user_nl_${rofName}
  else
    # Check to see if there is a raw ROF file in the landrawdir given by the user
    rawrofrestartfile=$(ls ${landrawdir}/*.${rofSpecialName}.r.${yearstr}-${monthstr}-${daystr}-${cyclestrsec}.nc || true)   # check for file, suppress failed ls error if true
    echo "rawrofrestartfile: ${rawrofrestartfile}"
    if [ ! -z ${rawrofrestartfile} ]; then   # if rawrofrestartfile string is NOT empty, add it.
      echo "WARNING USER_NL: Adding rawrofrestartfile for runoff, although no check to see if grid is consistent!"
      echo "USER_NL: Adding ${rawrofrestartfile} to user_nl_${rofName}"
      echo "finidat_rtm='${rawrofrestartfile}'" >> user_nl_${rofName}
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
xmlchange_verbose "REST_N" "9999"

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

cp user_nl_${atmName} user_nl_${atmName}.BAK

# Update env_run.xml with runtime parameters
xmlchange_verbose "RUN_STARTDATE" "$yearstr-$monthstr-$daystr"
xmlchange_verbose "START_TOD" "$cyclestrsec"
xmlchange_verbose "STOP_OPTION" "ndays"
xmlchange_verbose "STOP_N" "$numdays"

SEINIC=${sePreFilterIC}

############################### (IF) FILTER SETUP ###############################

if [ "$doFilter" = true ] ; then
  # If filtering, need to change these options

  xmlchange_verbose "STOP_OPTION" "$nhours"
  xmlchange_verbose "STOP_N" "$filterHourLength"

  sed -i '/.*ncdata/d' user_nl_${atmName}
  echo "ncdata='${SEINIC}'" >> user_nl_${atmName}

  # Do filter timestepping stability
  ATM_NCPL=192
  echo "ATM_NCPL: $ATM_NCPL  DTIME: $DTIME"
  if [[ "$DYCORE" == "se" ]]; then
    if [ "$use_nsplit" = true ]; then # if trad. SE nsplit timestepping
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

  if [ $debug = false ] ; then

    echo "Begin call to filter-run"
    run_CIME2 "$path_to_rundir" "$CIMEsubstring" "$CIMEbatchargs" true

    ## Run NCL filter
    cd $filter_path
    echo "Running filter"
    cp ${sePreFilterIC} ${sePostFilterIC}
    filtfile_name=${casename}.${atmName}.h0.$yearstr-$monthstr-$daystr-$cyclestrsec.nc
    set +e
    (set -x; ncl lowmemfilter.ncl \
     endhour=${filterHourLength} tcut=${filtTcut} \
    'filtfile_name = "'${path_to_rundir}'/'${filtfile_name}'"' \
    'writefile_name = "'${sePostFilterIC}'"' ) ; exit_status=$?
    check_ncl_exit "lowmemfilter.ncl" $exit_status
    set -e
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
    cp ${sePostFilterIC} ${FILTONLYDIR}/${casename}_FILTERED_${se_yearstr}-${se_monthstr}-${se_daystr}-${se_cyclestrsec}.nc
    cp ${sstFileIC} ${FILTONLYDIR}/${casename}_SST_1x1_${se_yearstr}-${se_monthstr}-${se_daystr}-${se_cyclestrsec}.nc
    mkdir ${FILTONLYDIR}/config_files
    mv ${path_to_nc_files}/*.gz ${FILTONLYDIR}/config_files
    mv ${path_to_nc_files}/*.nml ${FILTONLYDIR}/config_files
    mv ${path_to_nc_files}/*_in ${FILTONLYDIR}/config_files
    exit
  fi

  # Special things we have to do to "reset" model after filtering
  echo "Pushing model start back from $cyclestrsec to $se_cyclestrsec seconds..."
  cd $path_to_case
  xmlchange_verbose "RUN_STARTDATE" "$se_yearstr-$se_monthstr-$se_daystr"
  xmlchange_verbose "START_TOD" "$se_cyclestrsec"
  xmlchange_verbose "STOP_OPTION" "ndays"
  xmlchange_verbose "STOP_N" "$numdays"

  # Set new inic to point to post filter file
  SEINIC=${sePostFilterIC}
fi

############################### "ACTUAL" FORECAST RUN ###############################

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

# Delete existing ncdata and inject new one into user_nl_atm file
sed -i '/.*ncdata/d' user_nl_${atmName}
echo "ncdata='${SEINIC}'" >> user_nl_${atmName}

if [ $debug = false ] ; then
  echo "Begin call to forecast run"
  #run_CIME2 "$path_to_rundir" "$CIMEsubstring" "$CIMEbatchargs" true
  CIMESTATUS=1; CIMEITER=0
  while [ $CIMESTATUS != 0 ] ; do
    if [ "$CIMEITER" -ge "$CIMEMAXTRIES" ]; then
      echo "Exceeded the max tries: $CIMEMAXTRIES ... exiting"
      exit 1
    fi
    CIMEITER=$((CIMEITER+1))
    echo "Running CIME on CIMEITER $CIMEITER"
    run_CIME2 "$path_to_rundir" "$CIMEsubstring" "$CIMEbatchargs" false
  done
  echo "Returned status $CIMESTATUS"
fi

fi # end run model

cd $outputdir

# Generate folder structure and move NetCDF files
main_archive "$tmparchivecdir" "$atmName" "$lndName" "$rofName"

# Copy betacast configs to archive directory for posterity
cp -v $MACHINEFILE $NAMELISTFILE $OUTPUTSTREAMS $perturb_namelist $tmparchivecdir/betacast

# Copy user_nl* files to archive dir as well...
cp -v $path_to_case/user* $tmparchivecdir/nl_files

# Archive initial conditions?
if [ $archive_inic = true ]; then
  echo "BETACAST_USER: requesting initial condition files be archived"
  archive_inic "$tmparchivecdir" "$path_to_case" "$compress_history_nc" "$atmName" "$lndName" "$rofName" "$sstFileIC"
fi

# Archive nudging files generated by hindcasts
if [ $save_nudging_files = true ] ; then
  echo "BETACAST_USER: requesting nudging files be archived"
  archive_nudging "$tmparchivecdir" "$path_to_nc_files" "$compress_history_nc"
fi

## Move land files to new restart location
cd "$path_to_nc_files" || { echo "Failed to change directory to $path_to_nc_files"; exit 1; }
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
delete_leftovers "$path_to_nc_files" "$atmName" "$lndName" "$rofName"

## Delete run temp directory
rm -rfv $RUNTMPDIR

cd $path_to_nc_files
echo "Moving tmp archive directory to ARCHIVEDIR/YYYYMMDDHH"
if [ -d "${ARCHIVEDIR}/${yearstr}${monthstr}${daystr}${cyclestr}" ]
then
  echo "${ARCHIVEDIR}/${yearstr}${monthstr}${daystr}${cyclestr} already exists! Appending date string!"
  ARCHIVESUBDIR=${yearstr}${monthstr}${daystr}${cyclestr}_$(date +\%Y\%m\%d\%H\%M)
else
  ARCHIVESUBDIR=${yearstr}${monthstr}${daystr}${cyclestr}
fi
mv -v $tmparchivecdir ${ARCHIVEDIR}/${ARCHIVESUBDIR}

if [ "$dotracking" = true ] ; then
  echo "BETACAST_USER: requesting cyclones to be tracked."

  # Go to cyclone tracking folder...
  pushd "${SCRIPTPATH}/cyclone-tracking/" > /dev/null

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
      /bin/bash ./cyclone-tracking-driver.sh "${yearstr}${monthstr}${daystr}${cyclestr}" \
        "${casename}" \
        "${TCVITFILE}" \
        "${ATCFFILE}" \
        "${track_connectfile}" \
        "${path_to_rundir}" \
        "${track_sendhtml}" \
        "${track_hstream}" \
        "${track_stride}" \
        "${track_ATCFTECH}" \
        "${TE_SERIAL_DIR}"
    ) || { echo "Tracking code failed"; popd > /dev/null; exit 1; }
    cp -v "trajs.trajectories.txt.${casename}.png" "./fin-figs/trajs.trajectories.txt.${casename}.${yearstr}${monthstr}${daystr}${cyclestr}.png"
  else
    echo "No TCvitals file exists and/or no storms on said TCvitals file, no reason to run the tracking code"
  fi

  # Return to where we were...
  popd > /dev/null
fi

if [ "$sendplots" = true ]; then
  echo "BETACAST_USER: Sending plots to remote server!"
  upload_ncl_script="${SCRIPTPATH}/upload_ncl.sh"
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
    /bin/bash "${temp_upload_script}" "${nclPlotWeights}" "${outputdir}/${yearstr}${monthstr}${daystr}${cyclestr}" "${SCRIPTPATH}"
  )
  # Cleanup
  rm -f "${temp_upload_script}"
fi

# Compress model output streams
# Let's do this last so all the above scripts can operate on uncompressed files
if [ "$compress_history_nc" = true ]; then
  echo "BETACAST_USER: Requesting history be compressed"
  compress_history "${ARCHIVEDIR}/${ARCHIVESUBDIR}"
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
  cd ${SCRIPTPATH}
  #Remove top line from dates file
  remove_top_line_from_dates ${datesfile}

  AUTORESUB="yes"
  if [ $AUTORESUB == "yes" ]; then
    echo "*-*-*-* Automatically resubbing next date!"
    exec ./betacast.sh ${MACHINEFILE} ${NAMELISTFILE} ${OUTPUTSTREAMS}
  fi
fi

echo "DONE at $(date)"

exit 0
