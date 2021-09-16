#!/usr/bin/env python
#
# USAGE:
# python getInterim.py YYYYMMDD HH
# python getInterim.py 20050825 00

import sys 
import cdsapi

thedate=str(sys.argv[1])
theyear=thedate[0:4]
themonth=thedate[4:6]
theday=thedate[6:8]
thetime=str(sys.argv[2])

# (C) Copyright 2012-2013 ECMWF.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
# In applying this licence, ECMWF does not waive the privileges and immunities 
# granted to it by virtue of its status as an ntergovernmental organisation nor
# does it submit to any jurisdiction.
#


# To run this example, you need an API key 
# available from https://api.ecmwf.int/v1/key/
 

#thedate="20050825"
#thetime="12"

c = cdsapi.Client()

c.retrieve(
    'reanalysis-era5-single-levels',
    {
        'product_type': 'reanalysis',
        'format': 'netcdf',
        'variable': [
            'mean_sea_level_pressure', 'sea_ice_cover', 'sea_surface_temperature',
            'surface_pressure','2m_temperature','geopotential',
            ],
        'year': '%s'%(theyear),
        'month': '%s'%(themonth),
        'day': '%s'%(theday),
        'time': '%s:00'%(thetime),
    },
    'ERA5_sfc_%s%s.nc'%(thedate,thetime))
 

c.retrieve(
    'reanalysis-era5-pressure-levels',
    {
        'product_type': 'reanalysis',
        'format': 'netcdf',
        'variable': [
            'fraction_of_cloud_cover', 'geopotential', 'relative_humidity',
            'specific_cloud_ice_water_content', 'specific_cloud_liquid_water_content', 'specific_humidity',
            'temperature', 'u_component_of_wind', 'vertical_velocity',
            'v_component_of_wind',
            ],
        'pressure_level': [
            '1', '2', '3',
            '5', '7', '10',
            '20', '30', '50',
            '70', '100', '125',
            '150', '175', '200',
            '225', '250', '300',
            '350', '450', '500',
            '550', '600', '650',
            '700', '750', '775',
            '800', '825', '850',
            '875', '900', '925',
            '950', '975', '1000',
        ],
        'year': '%s'%(theyear),
        'month': '%s'%(themonth),
        'day': '%s'%(theday),
        'time': '%s:00'%(thetime),
    },
    'ERA5_ml_%s%s.nc'%(thedate,thetime))



