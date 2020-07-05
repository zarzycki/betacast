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
          'levtype' : "ml",
          'levelist': "1/2/3/4/5/6/7/8/9/10/11/12/13/14/15/16/17/18/19/20/21/22/23/24/25/26/27/28/29/30/31/32/33/34/35/36/37/38/39/40/41/42/43/44/45/46/47/48/49/50/51/52/53/54/55/56/57/58/59/60",
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
          'param'   : "31.128/134.128/151.128/34.128",
          'format'  : "netcdf",
          'target'  : "ERA-Int_sfc_%s%s.nc"%(thedate,thetime)
          })
