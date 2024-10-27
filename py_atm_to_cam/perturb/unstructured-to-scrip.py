import xarray as xr
import os
import sys

module_path = os.path.abspath(os.path.join(os.path.dirname(__file__), '../..', 'py_remapping'))
if module_path not in sys.path:
    sys.path.append(module_path)
from ESMF_regridding import *


ds = xr.open_dataset("mpasa3-60-tclf004_gmted2010_modis_bedmachine_nc3000_Laplace0100_20240712.nc")

unstructured_to_scrip('test_mesh.nc', ds['lat'].values, ds['lon'].values, debug=True)

#     is_valid = validate_scrip_file('test_mesh.nc')
#     if not is_valid:
#         print("File may not be suitable for ESMF")

# CMZ note, ESMF is 1-based, so need to subtract one here
#ESMF_CELL=91
#test_specific_cell('test_mesh.nc', ESMF_CELL-1)

#test_grid_coverage('test_mesh.nc')


