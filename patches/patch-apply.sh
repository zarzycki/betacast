#!/bin/bash

  local PATCHFLAGS=" "
  local PATCHDIR=${1}
  local filename=${2}
  local SUBPATH=${3}
  local CESMCOMPONENT=${4}
  echo "Patching ${filename}.F90..."
  cp ${CESMROOT}/${SUBPATH}/${filename}.F90 ./SourceMods/src.${CESMCOMPONENT}/
  patch ${PATCHFLAGS} ./SourceMods/src.${CESMCOMPONENT}/${filename}.F90 < ${PATCHDIR}/src.${CESMCOMPONENT}/${filename}.patch
  echo "----------------------------"

