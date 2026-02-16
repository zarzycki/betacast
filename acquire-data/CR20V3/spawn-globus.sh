NCAR_GLADE_GLOBUS="d33b3614-6d04-11e5-ba46-22000b92c6ec"
DOWNLOAD_DIR="/glade/derecho/scratch/zarzycki/globus"
YYYY=1862

#=================================================================================

# Defaults, only change if source data changes
NERSC_HPSS_GLOBUS="9cd89cfd-6d04-11e5-ba46-22000b92c6ec"

# Setting endpoints
DESTINATION_ENDPOINT=$NCAR_GLADE_GLOBUS
DESTINATION_DIRECTORY=$DOWNLOAD_DIR
SOURCE_ENDPOINT=$NERSC_HPSS_GLOBUS
SOURCE_DIRECTORY="/home/projects/incite11/www/20C_Reanalysis_version_3/everymember_anal_netcdf/subdaily/"

echo "Creating directory on destination endpoint..."
output=$(globus mkdir "$DESTINATION_ENDPOINT:$DESTINATION_DIRECTORY" 2>&1)
status=$?

if [ $status -ne 0 ]; then
  if echo "$output" | grep -q "MkdirFailed.Exists"; then
    echo "Directory already exists — continuing."
  else
    echo "Failed to create directory on destination endpoint."
    echo "$output"
    #exit 1
  fi
else
  echo "Directory created successfully on destination endpoint."
fi

echo "Getting list of possible files from source endpoint: ${SOURCE_ENDPOINT}"
globus ls --recursive "$SOURCE_ENDPOINT:$SOURCE_DIRECTORY" > tmp_files_from_hpss.txt

echo "Filtering list of files to only get relevant vars for Betacast"
bash _filter_globus.sh ${YYYY} tmp_files_from_hpss.txt tmp
rm -v tmp_files_from_hpss.txt

# Modify the batch request file for formatting purposes
awk -v base="$SOURCE_DIRECTORY" '{
    src=$0
    sub("^"base,"",src)
    split(src,a,"/")
    print src, a[length(a)]
}' tmp > files_${YYYY}.txt && rm -v tmp

echo "Submitting request"
globus transfer "$SOURCE_ENDPOINT:$SOURCE_DIRECTORY" "$DESTINATION_ENDPOINT:$DESTINATION_DIRECTORY" \
    --label "${YYYY} CR20V3" \
    --batch files_${YYYY}.txt

echo "Hooray!"