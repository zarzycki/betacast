#!/bin/bash

module use /glade/p/ral/jntp/MET/MET_releases/modulefiles
module load met/10.0

tc_stat \
  -lookin MET_pairs.tcst \
  -job summary -by AMODEL,LEAD -lead 12,24,36,48,60,72,84,96,120 -event_equal FALSE \
  -column TK_ERR -column AMSLP-BMSLP -column AMAX_WIND-BMAX_WIND > MET_stats.tcst
