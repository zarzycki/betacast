#!/bin/bash

module use /glade/p/ral/jntp/MET/MET_releases/modulefiles
module load met/10.0

tc_pairs \
  -adeck ./dorian_sample/adeck.2019.al05 \
  -bdeck ./dorian_sample/bdeck.2019.al05 \
  -config TCPairsConfig \
  -out MET_pairs
