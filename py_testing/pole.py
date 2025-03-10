import xarray as xr

# Load the files
file_ncl = 'ncl_cam_raw.nc'
file_py = 'py_cam_raw.nc'

# Open the datasets
ds_ncl = xr.open_dataset(file_ncl)
ds_py = xr.open_dataset(file_py)

# Get the div_cam variable from both datasets
div_cam_ncl = ds_ncl['div_cam']
div_cam_py = ds_py['div_cam']

# Extract the first and last latitude slices for the 0th lev
first_slice_ncl = div_cam_ncl.isel(lev=0).isel(lat=0)
last_slice_ncl = div_cam_ncl.isel(lev=0).isel(lat=-1)

first_slice_py = div_cam_py.isel(lev=0).isel(lat=0)
last_slice_py = div_cam_py.isel(lev=0).isel(lat=-1)

# Print the slices side by side
print("First latitude slice (NCL):")
print(first_slice_ncl.values)
print("First latitude slice (Python):")
print(first_slice_py.values)

print("\nLast latitude slice (NCL):")
print(last_slice_ncl.values)
print("Last latitude slice (Python):")
print(last_slice_py.values)
