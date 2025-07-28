#!/bin/bash

# -----------------------------------------------------------------------------
# wget_repeat
#
# Repeatedly attempts to download a file from a given URL using wget, retrying
# every 3 minutes upon failure. If the file cannot be successfully downloaded
# within a user-specified timeout, the function exits with error code 1.
#
# This function is safe to use in scripts with `set -e`, as it temporarily
# disables `errexit` during the retry loop and restores it afterward.
#
# Arguments:
#   $1 - URL of the file to download (required)
#   $2 - Maximum time allowed for retries, in seconds (optional; default: 10800s = 3 hours)
#
# Returns:
#   0 if the download succeeded within the timeout
#   1 if the timeout expired before a successful download
#
# Example usage:
#   wget_repeat "https://example.com/file.dat"           # default 3-hour timeout
#   wget_repeat "https://example.com/file.dat" 7200      # 2-hour timeout
#
# -----------------------------------------------------------------------------
wget_repeat() {
  local url="$1"
  local max_seconds="${2:-10800}"  # default 3 hours = 10800s
  local retry_delay=180            # 3 minutes between retries
  local start_time=$(date +%s)

  echo "Starting download retry loop for: $url"
  echo "Will give up after $((max_seconds / 60)) minutes"

  local error=1

  # Save and disable errexit (set -e) if it's on
  local errexit_was_set=0
  [[ $- == *e* ]] && errexit_was_set=1 && set +e

  while [ $error -ne 0 ]; do
    wget -nv --read-timeout=30 "$url"
    error=$?

    if [ $error -eq 0 ]; then
      echo "Download succeeded."
      break
    fi

    now=$(date +%s)
    elapsed=$((now - start_time))

    if [ $elapsed -ge $max_seconds ]; then
      echo "ERROR: Download failed after $elapsed seconds. Giving up."
      if [ "$errexit_was_set" -eq 1 ]; then set -e; fi
      return 1
    fi

    echo "Download failed. Waiting $retry_delay seconds and retrying..."
    sleep $retry_delay
  done

  # Restore errexit if it was set
  if [ "$errexit_was_set" -eq 1 ]; then set -e; fi

  return 0
}


get_gfs_atm() {
  echo "Getting GFS conditions"
  mkdir -p "$gfs_files_path"
  cd "$gfs_files_path" || exit
  LOCALGFSFILE="gfs_atm_${yearstr}${monthstr}${daystr}${cyclestr}.grib2"

  if [ ! -f "$LOCALGFSFILE" ]; then
    echo "Getting Atmo file"
    if [ "$islive" = true ] ; then
      gfsFTPPath="https://noaa-gfs-bdp-pds.s3.amazonaws.com/gfs.${yearstr}${monthstr}${daystr}/${cyclestr}/atmos/"
      gfsFTPFile="gfs.t${cyclestr}z.pgrb2.0p25.anl"
      rm -fv "$gfsFTPFile"
      echo "Attempting to download ${gfsFTPPath}${gfsFTPFile}"
      wget_repeat "${gfsFTPPath}${gfsFTPFile}"
    else
      rm -f gfs.t*pgrb2f00*
      gfsFTPPath="/glade/collections/rda/data/ds084.1/${yearstr}/${yearstr}${monthstr}${daystr}/"
      gfsFTPFile="gfs.0p25.${yearstr}${monthstr}${daystr}${cyclestr}.f000.grib2"
      cp "${gfsFTPPath}/${gfsFTPFile}" .
      echo "Attempting to copy ${gfsFTPPath}${gfsFTPFile}"
    fi
    mv -v "$gfsFTPFile" "$LOCALGFSFILE"
  else
    echo "${LOCALGFSFILE} already exists, skipping download"
  fi
  echo "... done with get_gfs_atm"
}


get_era_interim_atm() {
  echo "Using ERA-Interim forecast ICs"
  echo "Changing to ERA-Interim interpolation directory"
  mkdir -p "$era_files_path"
  cd "$era_files_path" || exit
  LOCALGFSFILE="ERA-Int_${yearstr}${monthstr}${daystr}${cyclestr}.nc"

  if [ ! -f "$LOCALGFSFILE" ]; then
    echo "Support broken for auto download ERA, please prestage"
    exit 1
  fi
  echo "... done with get_era_interim_atm"
}


get_cfsr_atm() {
  echo "Using CFSR ICs"
  echo "Changing to GFS interpolation directory since they are practically the same thing"
  mkdir -p "$gfs_files_path"
  cd "$gfs_files_path" || exit

  LOCALCFSRFILE="cfsr_atm_${yearstr}${monthstr}${daystr}${cyclestr}.grib2"
  if [ ! -f "$LOCALCFSRFILE" ]; then
    STCUTARR=(26 21 16 11 06 01)
    ENCUTARR=(99 25 20 15 10 05)
    index=0

    for FILEDAY in "${STCUTARR[@]}" ; do
      if [ "$daystr" -ge "$FILEDAY" ] ; then
        break
      fi
      index=$((index+1))
    done

    if [ "$index" -eq 0 ] ; then
      ENCUTARR[0]=$(date -d "$monthstr/1 + 1 month - 1 day" "+%d")
      echo "Last day of month ($monthstr) is ${ENCUTARR[0]}"
    fi

    echo "Getting file: ${CFSRFILENAME}"
    CFSRFILENAME="pgbhnl.gdas.${yearstr}${monthstr}${STCUTARR[$index]}-${yearstr}${monthstr}${ENCUTARR[$index]}.tar"
    if [[ $(hostname -s) = derecho* ]]; then
      cp "/glade/collections/rda/data/ds093.0/${yearstr}/${CFSRFILENAME}" .
    else
      wget -nv --load-cookies ~/.thecookies "http://rda.ucar.edu/data/ds093.0/${yearstr}/${CFSRFILENAME}"
    fi
    tar -xvf "$CFSRFILENAME"
    mv "pgbhnl.gdas.${yearstr}${monthstr}${daystr}${cyclestr}.grb2" "$LOCALCFSRFILE"
    rm pgbhnl.gdas.*
  fi
  echo "... done with get_cfsr_atm"
}


get_era5_atm() {
  echo "Using ERA5 forecast ICs"
  echo "Changing to ERA5 interpolation directory"
  mkdir -p "$era_files_path"
  cd "$era_files_path" || exit
  LOCALGFSFILE="ERA5_${yearstr}${monthstr}${daystr}${cyclestr}.nc"

  if [ ! -f "$LOCALGFSFILE" ]; then
    echo "Cannot find: ${era_files_path}/${LOCALGFSFILE}"
    if [[ "$MACHINEFILE" == *derecho* ]]; then
      echo "We are on derecho, so even though we lack a local file, we can use RDA"
      ERA5RDA=1
      RDADIR="/glade/campaign/collections/rda/data/d633000/"
    elif [[ "$MACHINEFILE" == *pm* ]]; then
      echo "We are on Cori, so even though we lack a local file, we can use RDA"
      ERA5RDA=1
      RDADIR="/global/cfs/projectdirs/m3522/cmip6/ERA5/"
    else
      echo "Support broken for auto download ERA, please prestage!"
      exit 1
    fi
  else
    echo "Found a local file --> using $PWD/$LOCALGFSFILE"
  fi
  echo "... done with get_era5_atm"
}


# Function to get GDAS SST data
get_gdas_sst() {
  SSTTYPE=GDAS
  echo "Getting SST data"
  mkdir -p "${sst_files_path}"
  cd "${sst_files_path}" || exit

  if [ "$islive" = true ]; then
    sstFTPPath="https://noaa-gfs-bdp-pds.s3.amazonaws.com/nsst.${yestyearstr}${yestmonthstr}${yestdaystr}/"
    sstFTPFile="rtgssthr_grb_0.5.grib2"
    rm -fv "${sstFTPFile}"
    echo "Attempting to download ${sstFTPPath}${sstFTPFile}"
    wget_repeat "${sstFTPPath}${sstFTPFile}"

    sstFile="gfs_sst_${yearstr}${monthstr}${daystr}${cyclestr}.grib2"
    mv -v "${sstFTPFile}" "${sstFile}"

    # Hack to convert grb2 to nc because ncl routines not working
    #sstFile="gfs_sst_${yearstr}${monthstr}${daystr}${cyclestr}.nc"
    #cdo -f nc copy "${sstFTPFile}" "${sstFile}" ; rm "${sstFTPFile}"

    iceFile="" # Do not need ice file since ice is stored in the SST file
  else
    echo "NCEP broke support for historical GDAS, use NOAAOI instead."
    exit 1
  fi
  echo "... done with get_gdas_sst"
}


# Function to handle unsupported ERA SST
get_erai_sst() {
  echo "ERA-I SST not quite supported yet..."
  exit 1
  echo "... done with get_erai_sst"
}


# Function to get NOAAOI SST data
get_noaaoi_sst() {
  SSTTYPE=NOAAOI
  echo "Using NOAAOI SSTs"
  mkdir -p "${sst_files_path}"
  cd "${sst_files_path}" || exit

  sstFile="sst.day.mean.${yearstr}.nc"
  if [ ! -f "${sst_files_path}/${sstFile}" ] || \
     { [ -f "${sst_files_path}/${sstFile}" ] && [ "$(ncdmnsz time "${sst_files_path}/${sstFile}")" -lt 365 ]; }
  then
    echo "NOAAOI SST file doesn't exist or has less than 365 timesteps, need to download"
    rm -f "${sst_files_path}/${sstFile}"
    sstFTPPath="ftp://ftp.cdc.noaa.gov/Datasets/noaa.oisst.v2.highres/"
    wget_repeat "${sstFTPPath}/${sstFile}"
  fi

  iceFile="icec.day.mean.${yearstr}.nc"
  if [ ! -f "${sst_files_path}/${iceFile}" ] || \
     { [ -f "${sst_files_path}/${iceFile}" ] && [ "$(ncdmnsz time "${sst_files_path}/${iceFile}")" -lt 365 ]; }
  then
    echo "NOAAOI ice file doesn't exist or has less than 365 timesteps, need to download"
    rm -f "${sst_files_path}/${iceFile}"
    sstFTPPath="ftp://ftp.cdc.noaa.gov/Datasets/noaa.oisst.v2.highres/"
    wget_repeat "${sstFTPPath}/${iceFile}"
  fi
  echo "... done with get_noaaoi_sst"
}
