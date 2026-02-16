#!/usr/bin/env bash
set -euo pipefail

YYYY=$1
INFILE=$2
OUTFILE=$3

# Standalone 2D vars
vars2D=(
  PRES
  TMP2m
)

# 3D vars + levels
vars=(
  HGT
  SPFH
  TMP
  UGRD
  VGRD
)

levels=(
  200 250 300 400 500
  600 650 700 750 800 850
  900 925 950 975 1000
)

# Build list of tags
tags=()
for v in "${vars2D[@]}"; do
  tags+=("$v")
done
for v in "${vars[@]}"; do
  for lev in "${levels[@]}"; do
    tags+=("${v}${lev}")
  done
done

# Build grep pattern: TAG/TAG_YYYY*.tar (but not .tar.idx)
pattern=$(printf '%s\n' "${tags[@]}" | sed "s|.*|^&/&_${YYYY}.*\\\\.tar\$|" | paste -sd'|')

grep -E "$pattern" "$INFILE" > "$OUTFILE"

echo "Extracted $(wc -l < "$OUTFILE") of ${#tags[@]} expected entries to $OUTFILE"

# Report missing tags
for tag in "${tags[@]}"; do
  if ! grep -q "^${tag}/" "$OUTFILE"; then
    echo "  MISSING: ${tag}"
  fi
done
