#!/bin/bash -l
#PBS -N casper-CR20V3-mean
#PBS -A P93300042
#PBS -l select=1:ncpus=1:mem=660GB
#PBS -l walltime=23:59:00
#PBS -q casper@casper-pbs
#PBS -j oe

TARDIR="/glade/derecho/scratch/zarzycki/globus/"
YEAR=1862

module load nco

cd $TARDIR/$YEAR

# Get unique variable prefixes (everything before _mem)
for prefix in $(ls *.${YEAR}_mem*.nc | grep -v daily | sed 's/_mem[0-9]*\.nc//' | sort -u); do
  outfile="${prefix}_mem999.nc"
  echo "========================================"
  echo "Processing: $prefix"
  echo "Output:   $outfile"
  echo "========================================"
  time ncea ${prefix}_mem*.nc "$outfile"
done
