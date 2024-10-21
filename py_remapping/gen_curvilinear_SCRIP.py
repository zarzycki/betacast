import xarray as xr
from ESMF_regridding import curvilinear_to_SCRIP  # Make sure to import your custom functions

# Full path to the nc or GRIB file of the analysis
srcFileName = "/glade/u/home/zarzycki/work/2012_derecho/hrrr.t00z.wrfprsf00.grib2"
# Filename to write SCRIP grid to
srcGridName = "./hrrr_3km_scrip_tmp.nc"

# Open the source file
sfile = xr.open_dataset(srcFileName, engine="cfgrib")

# Extract the 2D latitude and longitude arrays
lat2d = sfile['gridlat_0'].values
lon2d = sfile['gridlon_0'].values

# Set options for the SCRIP file generation
Opt = {
    "ForceOverwrite": True,
    "PrintTimings": True,
    "Title": "RAP Grid",
    "Debug": True
    # Uncomment and modify the following line if you want to apply a mask
    # "GridMask": np.where(~np.isnan(sfile['some_variable'].values), 1, 0)
}

# Generate the SCRIP file
curvilinear_to_SCRIP(srcGridName, lat2d, lon2d, Opt)

print("SCRIP file generation completed.")

