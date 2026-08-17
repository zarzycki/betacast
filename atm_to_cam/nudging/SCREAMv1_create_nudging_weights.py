'''
Created  by: Xingqiu 2023/06: This file read in the scream output netcdf file to create the nudging weights for regional nudging in scream
Modified by: zhang73 and bogensch 2024/10: The code has been modified to follow the habits of a E3SM/SCREAM user

USAGE:
Generate nudging weights file for your RRM, on the PG2 grid.

Example:
python3 SCREAMv1_create_nudging_weights.py -datafile USGS-gtopo30_CA_ne32_x32_v1_pg2_16xdel2.nc -nlev 128 -lat lat -lon lon -weightsfile your_weighting_file.nc

DETAILS:
 - datafile: Path to a file where the latitude and longitude coordinates are on the PG2 grid (NOT dynamics/GLL grid).
   The easist thing to do is to point to your TOPOGRAPHY file, which will satisy such a requirement.
 - nlev: Set the number of levels that your model run uses.  SCREAMv1 default is 128.
 - lat: The name of the latitude coordinate variable in your datafile.
 - lon: The name of the longitude coordinate variable in your datafile.
 - weightsfile: What you want your output file to be named.

You will also likely want to set the nudging window to fit your needs (see below):

'''
#######################################################
#SET THE NUDGING WINDOW
#---same as user_nl_eam passed to components/eam/src/physics/cam/nudging.F90.  See this file
#    for detailed comments for each of the following parameters or see the
#    components/eam/bld/namelist_files/namelist_definition.xml file.
#  Below are the select parameters that most users will be concerned with setting.  See the
#    create_nudging_weights function for some advanced settings.

Nudge_Hwin_lat0      = 37.2  # latitudinal center of window in degrees.
Nudge_Hwin_latWidth  = 0.0  # latitudinal width of window in degrees.
Nudge_Hwin_latDelta  = 1e-6   # latitudinal transition length of window in degrees.
Nudge_Hwin_lon0      = 240.6 # longitudinal center of window in degrees.
Nudge_Hwin_lonWidth  = 0.0  # longitudinal width of window in degrees.
Nudge_Hwin_lonDelta  = 1e-6   # longitudinal transition length of window in degrees.

# END USER DEFINED SETTINGS
########################################################

import argparse
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import numpy as np
import xarray as xr

def create_nudging_weights(opt):
   '''create nudging weights for scream regional nudging'''

   ds = xr.open_dataset(opt.datafile)

   lat  = ds[opt.lat]
   lon  = ds[opt.lon]
   nlev = int(opt.nlev)

   ncol=len(lat)
   print("Number of levels:", nlev)
   print("Number of columns:", ncol)

   alpha    = 1

   # Default horizontal nudging weights.
   Nudge_Hwin_lo        = 1.0   # LOW Coeffcient for Horizontal Window.
   Nudge_Hwin_hi        = 0.0   # HIGH Coeffcient for Horizontal Window.

   # Set nudging weights related to the vertical. Only advance users modify if you want
   #   to apply nudging to specific levels
   Nudge_Vwin_lo          =0.0  # LOW Coeffcient for Vertical Window.
   Nudge_Vwin_hi          =1.0  # HIGH Coeffcient for Vertical Window.
   Nudge_Vwin_Hindex      =92   # HI model index of transition (note, this is one level "down" from ndg midpoint (i.e., if set to 88, you really are nudging to 87
   Nudge_Vwin_Hdelta      =1.0  # HI transition length
   Nudge_Vwin_Lindex      =0.0  # LO model index of transition
   Nudge_Vwin_Ldelta      =0.1  # LO transition length

   levidx=np.arange(1,nlev+1)
   lev_lo=(levidx-Nudge_Vwin_Lindex)/Nudge_Vwin_Ldelta
   lev_hi=(Nudge_Vwin_Hindex-levidx)/Nudge_Vwin_Hdelta
   Wprof=((1.+np.tanh(lev_lo))/2.)*((1.+np.tanh(lev_hi))/2.)
   Vmax=np.max(Wprof)
   Vmin=np.min(Wprof)
   if Vmax <= Vmin:
       Vmax= max(Nudge_Vwin_lo,Nudge_Vwin_hi)
       Wprof[:] = Vmax
   else:
     Wprof[:] = (Wprof-Vmin)/(Vmax-Vmin)
     Wprof[:] = Nudge_Vwin_lo + Wprof*(Nudge_Vwin_hi-Nudge_Vwin_lo)

   # Assuming ntime is 1 since the new input doesn't have time dimension
   ntime = 1

   # Create Wprof as a DataArray with vertical levels (lev)
   Wprof = xr.DataArray(Wprof, dims=["lev"])
   # Now, create a 3D array of Wprof with dimensions (time, ncol, lev)
   # Since you now only have ncol, we broadcast Wprof over the ncol and time dimensions
   Wprof_3D = Wprof.expand_dims({"time": ntime, "ncol": ncol})

   Nudge_Hwin_lonWidthH=Nudge_Hwin_lonWidth/2.
   Nudge_Hwin_latWidthH=Nudge_Hwin_latWidth/2.

   lonp= 180.
   lon0=   0.
   lonn=-180.
   latp=  90.-Nudge_Hwin_lat0
   lat0=   0.
   latn= -90.-Nudge_Hwin_lat0

   Val1_p=(1.+np.tanh((Nudge_Hwin_lonWidthH+lonp)/Nudge_Hwin_lonDelta))/2.
   Val2_p=(1.+np.tanh((Nudge_Hwin_lonWidthH-lonp)/Nudge_Hwin_lonDelta))/2.
   Val3_p=(1.+np.tanh((Nudge_Hwin_latWidthH+latp)/Nudge_Hwin_latDelta))/2.
   Val4_p=(1.+np.tanh((Nudge_Hwin_latWidthH-latp)/Nudge_Hwin_latDelta))/2.

   Val1_0=(1.+np.tanh((Nudge_Hwin_lonWidthH+lon0)/Nudge_Hwin_lonDelta))/2.
   Val2_0=(1.+np.tanh((Nudge_Hwin_lonWidthH-lon0)/Nudge_Hwin_lonDelta))/2.
   Val3_0=(1.+np.tanh((Nudge_Hwin_latWidthH+lat0)/Nudge_Hwin_latDelta))/2.
   Val4_0=(1.+np.tanh((Nudge_Hwin_latWidthH-lat0)/Nudge_Hwin_latDelta))/2.

   Val1_n=(1.+np.tanh((Nudge_Hwin_lonWidthH+lonn)/Nudge_Hwin_lonDelta))/2.
   Val2_n=(1.+np.tanh((Nudge_Hwin_lonWidthH-lonn)/Nudge_Hwin_lonDelta))/2.
   Val3_n=(1.+np.tanh((Nudge_Hwin_latWidthH+latn)/Nudge_Hwin_latDelta))/2.
   Val4_n=(1.+np.tanh((Nudge_Hwin_latWidthH-latn)/Nudge_Hwin_latDelta))/2.

   Nudge_Hwin_max=     Val1_0*Val2_0*Val3_0*Val4_0
   Nudge_Hwin_min=min((Val1_p*Val2_p*Val3_n*Val4_n),
                      (Val1_p*Val2_p*Val3_p*Val4_p),
                      (Val1_n*Val2_n*Val3_n*Val4_n),
                      (Val1_n*Val2_n*Val3_p*Val4_p))

   phi_s                = Nudge_Hwin_lat0 - Nudge_Hwin_lonWidthH
   phi_n                = Nudge_Hwin_lat0 + Nudge_Hwin_lonWidthH
   delta_phi_ns         = Nudge_Hwin_latDelta
   phi_w                = Nudge_Hwin_lon0 - Nudge_Hwin_lonWidthH
   phi_e                = Nudge_Hwin_lon0 + Nudge_Hwin_lonWidthH
   delta_phi_we         = Nudge_Hwin_lonDelta


   phis_ratio = np.tanh((lat-phi_s)/delta_phi_ns)
   phin_ratio = np.tanh((lat-phi_n)/delta_phi_ns)
   phiw_ratio = np.tanh((lon-phi_w)/delta_phi_we)
   phie_ratio = np.tanh((lon-phi_e)/delta_phi_we)

   weight_lat_profile  = 0.25*alpha*(1.+phis_ratio)*(1-phin_ratio)
   weight_lon_profile  = 0.25*alpha*(1.+phiw_ratio)*(1-phie_ratio)

   Hcoef = weight_lat_profile * weight_lon_profile

   #---Scale the horizontal window coef for specified range of values.
   Hcoef=(Hcoef-Nudge_Hwin_min)/(Nudge_Hwin_max-Nudge_Hwin_min)
   Hcoef=(1.-Hcoef)*Nudge_Hwin_lo + Hcoef*Nudge_Hwin_hi
   nudging_weights = Wprof_3D * Hcoef

   nudging_weights_ = nudging_weights.to_dataset(name='nudging_weights')
   nudging_weights_['lat'] = lat
   nudging_weights_['lon'] = lon
   nudging_weights_.attrs['name']='nudging_weights'
   nudging_weights_.attrs['units']='None'
   nudging_weights_.attrs['long_name']='nudging weights'
   nudging_weights_.to_netcdf(opt.weightsfile)

   print("Created nudging weights file:",opt.weightsfile)

def main():

   parser = argparse.ArgumentParser()
   parser.add_argument('-datafile', default=None)
   parser.add_argument('-weightsfile', default='nudging_weights.nc')
   parser.add_argument('-nlev',default=128)
   parser.add_argument('-lat',default='lat')
   parser.add_argument('-lon',default='lon')

   opt = parser.parse_args()
   create_nudging_weights(opt)

if __name__ == '__main__':
   main()
