import netCDF4 as nc
import numpy as np

def select_levels(input_file, selected_indices):

    num_selected_levels = len(selected_indices)

    output_file = f"L{num_selected_levels:02d}template_new.nc"

    with nc.Dataset(input_file, 'r') as src:
        with nc.Dataset(output_file, 'w') as dst:
            for name in src.ncattrs():
                dst.setncattr(name, src.getncattr(name))

            for name, dimension in src.dimensions.items():
                if name == 'lev':
                    dst.createDimension(name, len(selected_indices))
                elif name == 'ilev':
                    dst.createDimension(name, len(selected_indices) + 1)
                else:
                    dst.createDimension(name, len(dimension) if not dimension.isunlimited() else None)

            for name, variable in src.variables.items():
                if name in ['hyam', 'hybm', 'lev']:
                    data = variable[selected_indices]
                elif name in ['hyai', 'hybi', 'ilev']:
                    data = variable[selected_indices + [selected_indices[-1] + 1]]
                else:
                    data = variable[:]

                if name == 'lev':
                    print(data)

                dst_var = dst.createVariable(name, variable.datatype, variable.dimensions)

                dst_var.setncatts({k: variable.getncattr(k) for k in variable.ncattrs()})

                dst_var[:] = data

    print(f"Selected levels written to {output_file}")

# Example usage:
input_file = 'L26template.nc'
# Indices of the levels you grab from input_file
selected_indices = [0, 5, 10, 13, 17, 19, 21, 22, 23, 25]
select_levels(input_file, selected_indices)