#!/usr/bin/env python
#
# USAGE:
# python getInterim.py YYYYMMDD HH
# python getInterim.py 20050825 00

import sys

thedate=str(sys.argv[1])
thetime=str(sys.argv[2])

# (C) Copyright 2012-2013 ECMWF.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
# In applying this licence, ECMWF does not waive the privileges and immunities 
# granted to it by virtue of its status as an ntergovernmental organisation nor
# does it submit to any jurisdiction.
#

from ecmwfapi import ECMWFDataServer

# To run this example, you need an API key 
# available from https://api.ecmwf.int/v1/key/
 
server = ECMWFDataServer()

#thedate="20050825"
#thetime="12"
  
server.retrieve({
          'stream'  : "oper",
          'dataset' : "interim",
          'levtype' : "pl",
          'levelist': "1/2/3/5/7/10/20/30/50/70/100/125/150/175/200/225/250/300/350/400/450/500/550/600/650/700/750/775/800/825/850/875/900/925/950/975/1000",
          'date'    : "%s"%(thedate),
          'time'    : "%s"%(thetime),
          'step'    : "0",
          'type'    : "an",
          'class'   : "ei",
          'grid'    : "0.75/0.75",
          'param'   : "129.128/130.128/131.128/132.128/133.128/135.128/157.128/246.128/247.128/248.128",
          'format'  : "netcdf",
          'target'  : "ERA-Int_ml_%s%s.nc"%(thedate,thetime)
          })
          
server.retrieve({
          'stream'  : "oper",
          'dataset' : "interim",
          'levtype' : "sfc",
          'date'    : "%s"%(thedate),
          'time'    : "%s"%(thetime),
          'step'    : "0",
          'type'    : "an",
          'class'   : "ei",
          'grid'    : "0.75/0.75",
          'param'   : "134.128/151.128/34.128",
          'format'  : "netcdf",
          'target'  : "ERA-Int_sfc_%s%s.nc"%(thedate,thetime)
          })
