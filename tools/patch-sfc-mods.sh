#!/bin/bash

#Usage:
#~/betacast/tools/patch-sfc-mods.sh ~/betacast/ /glade/u/home/zarzycki/work/cam_20230623/ nuopc clm

BETACAST=$1
PATHTOCESM=$2
COUPLER=$3
COMPONENT_NAME=$4

declare -A component_types
component_types=(
    [clm]="lnd"
    [elm]="lnd"
    [mosart]="rof"
    [rtm]="rof"
)

## Here, we build a lookup table for the relevant component/coupler pair
declare -A component_patterns
component_patterns=(
    [clm]="*clm/*lnd_comp_$COUPLER.F90*"
    [elm]="*elm/*lnd_comp_$COUPLER.F90*"
    [mosart]="*mosart/*rof_comp_$COUPLER.F90*"
    [rtm]="*mosart/*rof_comp_$COUPLER.F90*"
)

## Get the clean component name (remove any underscores)
clean_component_name="${COMPONENT_NAME%%_*}"
echo $clean_component_name
if [[ $COMPONENT_NAME == *_* ]]; then
  # If underscore exists, extract special string
  special_string="_${COMPONENT_NAME#*_}"
else
  # Otherwise, we don't need a special string, so set to empty
  special_string=""
fi
echo $special_string

## If the user gave us a bad component key, stop.
if [[ -z ${component_patterns[$clean_component_name]} ]]; then
    echo "Invalid component key. Exiting."
    exit 1
fi

# Get this model's component type. This is "lnd", "rof", etc.
component_type=${component_types[$clean_component_name]}
echo "component_type: $component_type"

# Get the pattern from above
this_pattern=${component_patterns[$clean_component_name]}
echo "this_pattern: $this_pattern"

# Find the relevant file in the source tree using the pattern
echo "find ${PATHTOCESM} -type f -path "*/components/$this_pattern""
comp_restart_file=$(find ${PATHTOCESM} -type f -path "*/components/$this_pattern")
echo "comp_restart_file: $comp_restart_file"

# Check if the find command returned multiple results
count=$(echo "$comp_restart_file" | wc -l)
echo "count: $count"

# Error check if multiple files returned (need stricter pattern) or no files returned...
if [ "$count" -gt 1 ]; then
  echo "Multiple matches found. Exiting."
  echo "..." ; echo "comp_restart_file: $comp_restart_file" ; echo "..."
  echo "Need a stricter pattern..."
  exit 1
elif [ -z "$comp_restart_file" ]; then
  echo "No matches found. Exiting."
  exit 1
fi

echo "Found file: $comp_restart_file"

echo "Component: $COMPONENT_NAME"

filename="${component_type}_comp_$COUPLER.F90"
cp -v $comp_restart_file ./SourceMods/src.$clean_component_name
patch ./SourceMods/src.$clean_component_name/$filename < ${BETACAST}/patches/${filename%.F90}${special_string}.patch
