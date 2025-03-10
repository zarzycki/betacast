#!/bin/bash

get_gfs_atm() {
  echo "Getting GFS conditions"
  mkdir -p "$gfs_files_path"
  cd "$gfs_files_path" || exit
  LOCALGFSFILE="gfs_atm_${yearstr}${monthstr}${daystr}${cyclestr}.grib2"

  if [ ! -f "$LOCALGFSFILE" ]; then
    echo "Getting Atmo file"
    if [ "$islive" = true ] ; then
      gfsFTPPath="https://ftpprd.ncep.noaa.gov/data/nccf/com/gfs/prod/gfs.${yearstr}${monthstr}${daystr}/${cyclestr}/atmos/"
      gfsFTPFile="gfs.t${cyclestr}z.pgrb2.0p25.anl"
      rm -fv "$gfsFTPFile"
      echo "Attempting to download ${gfsFTPPath}${gfsFTPFile}"

      error=1
      while [ $error != 0 ] ; do
        wget -nv --read-timeout=30 "$gfsFTPPath$gfsFTPFile"
        error=$?
        if [ $error -ne 0 ] ; then
          echo "Cannot get file, will wait 2 min and scrape again"
          sleep 120
        fi
      done
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
      RDADIR="/glade/campaign/collections/rda/data/ds633.0/"
    elif [[ "$MACHINEFILE" == *pm* ]]; then
      echo "We are on Cori, so even though we lack a local file, we can use RDA"
      ERA5RDA=1
      RDADIR="/global/cfs/projectdirs/m3522/cmip6/ERA5/"
    else
      echo "Support broken for auto download ERA, please prestage!"
      exit 1
    fi
  fi
}

# Function to get GDAS SST data
get_gdas_sst() {
  SSTTYPE=GDAS
  echo "Getting SST data"
  mkdir -p "${sst_files_path}"
  cd "${sst_files_path}" || exit

  if [ "$islive" = true ]; then
    sstFTPPath="https://ftpprd.ncep.noaa.gov/data/nccf/com/nsst/v1.2/nsst.${yestyearstr}${yestmonthstr}${yestdaystr}/"
    sstFTPFile="rtgssthr_grb_0.5.grib2"
    rm -fv "${sstFTPFile}"
    echo "Attempting to download ${sstFTPPath}${sstFTPFile}"

    error=1
    while [ $error != 0 ]; do
      wget -nv --read-timeout=30 "${sstFTPPath}${sstFTPFile}"
      error=$?
      if [ $error -ne 0 ]; then
        echo "Cannot get file, will wait 2 min and scrape again"
        sleep 120
      fi
    done

    #sstFile="gfs_sst_${yearstr}${monthstr}${daystr}${cyclestr}.grib2"
    #mv -v "${sstFTPFile}" "${sstFile}"

    # Hack to convert grb2 to nc because ncl routines not working
    sstFile="gfs_sst_${yearstr}${monthstr}${daystr}${cyclestr}.nc"
    cdo -f nc copy "${sstFTPFile}" "${sstFile}" ; rm "${sstFTPFile}"

    iceFile="" # Do not need ice file since ice is stored in the SST file
  else
    echo "NCEP broke support for historical GDAS, use NOAAOI instead."
    exit 1
  fi
}

# Function to handle unsupported ERA SST
get_erai_sst() {
  echo "ERA-I SST not quite supported yet..."
  exit 1
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
    error=1
    while [ $error != 0 ]; do
      wget -nv "${sstFTPPath}/${sstFile}"
      error=$?
      if [ $error -ne 0 ]; then
        echo "Cannot get file, will wait 2 min and scrape again"
        sleep 120
      fi
    done
  fi

  iceFile="icec.day.mean.${yearstr}.nc"
  if [ ! -f "${sst_files_path}/${iceFile}" ] || \
     { [ -f "${sst_files_path}/${iceFile}" ] && [ "$(ncdmnsz time "${sst_files_path}/${iceFile}")" -lt 365 ]; }
  then
    echo "NOAAOI ice file doesn't exist or has less than 365 timesteps, need to download"
    rm -f "${sst_files_path}/${iceFile}"
    sstFTPPath="ftp://ftp.cdc.noaa.gov/Datasets/noaa.oisst.v2.highres/"
    error=1
    while [ $error != 0 ]; do
      wget -nv "${sstFTPPath}/${iceFile}"
      error=$?
      if [ $error -ne 0 ]; then
        echo "Cannot get file, will wait 2 min and scrape again"
        sleep 120
      fi
    done
  fi
}
