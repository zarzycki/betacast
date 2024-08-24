import sys
import xarray as xr
import numpy as np

def calculate_mean_absolute_error(var1, var2, weights):
    return np.sum(np.abs(var1 - var2) * weights) / np.sum(weights)

def calculate_correlation(var1, var2, weights):
    var1_mean = np.sum(var1 * weights) / np.sum(weights)
    var2_mean = np.sum(var2 * weights) / np.sum(weights)
    covariance = np.sum(weights * (var1 - var1_mean) * (var2 - var2_mean)) / np.sum(weights)
    var1_variance = np.sum(weights * (var1 - var1_mean) ** 2) / np.sum(weights)
    var2_variance = np.sum(weights * (var2 - var2_mean) ** 2) / np.sum(weights)
    return covariance / np.sqrt(var1_variance * var2_variance)

def calculate_normalized_mean_bias(var1, var2, weights):
    mean_var1 = np.sum(var1 * weights) / np.sum(weights)
    if np.isnan(mean_var1) or mean_var1 == 0:
        print("   Warning: Mean of var1 is NaN or zero, NMB calculation skipped.")
        return np.nan
    return np.sum((var2 - var1) * weights) / (mean_var1 * np.sum(weights))

def calculate_normalized_rmse(var1, var2, weights):
    rmse = np.sqrt(np.sum(weights * (var1 - var2) ** 2) / np.sum(weights))
    mean_var1 = np.sum(var1 * weights) / np.sum(weights)
    if np.isnan(mean_var1) or mean_var1 == 0:
        print("   Warning: Mean of var1 is NaN or zero, NRMSE calculation skipped.")
        return np.nan
    return rmse / mean_var1

def check_same(file1, file2, variables_to_check=None):
    ds1 = xr.open_dataset(file1)
    ds2 = xr.open_dataset(file2)

    lat = ds1['lat'].values
    cos_lat = np.cos(np.radians(lat))

    common_vars = set(ds1.variables).intersection(ds2.variables)

    if variables_to_check:
        common_vars = common_vars.intersection(variables_to_check)

    for var_name in common_vars:
        data1 = ds1[var_name].values
        data2 = ds2[var_name].values

        # Skip comparison if data types are not compatible
        if data1.dtype != data2.dtype:
            print(f"Skipping variable {var_name} due to incompatible data types: {data1.dtype} vs {data2.dtype}")
            continue

        if np.issubdtype(data1.dtype, np.number):
            # Only perform calculations for numeric types
            if data1.shape == data2.shape:
                print(f"Variable: {var_name}")

                # Determine whether to expand cos_lat to match 2D or 3D data
                if len(data1.shape) == 3:
                    weights = cos_lat[np.newaxis, :, np.newaxis]  # Expand for 3D data
                elif len(data1.shape) == 2:
                    weights = cos_lat[:, np.newaxis]  # Expand for 2D data
                else:
                    print(f"Skipping variable {var_name} due to unexpected dimensions: {data1.shape}")
                    continue

                mae = calculate_mean_absolute_error(data1, data2, weights)
                correlation = calculate_correlation(data1, data2, weights)
                mean_value = np.sum(data1 * weights) / np.sum(weights)
                nmb = calculate_normalized_mean_bias(data1, data2, weights)
                nrmse = calculate_normalized_rmse(data1, data2, weights)

                print(f"   Mean Absolute Error: {mae}")
                print(f"   Correlation: {correlation}")
                print(f"   Mean Value: {mean_value}")
                print(f"   Normalized Mean Bias (NMB): {nmb}")
                print(f"   Normalized RMSE (NRMSE): {nrmse}")
            else:
                print(f"Variable: {var_name} has different shapes: {data1.shape} vs {data2.shape}")
        else:
            print(f"Skipping non-numeric variable {var_name} with dtype {data1.dtype}")

    ds1.close()
    ds2.close()

if __name__ == "__main__":
    if len(sys.argv) < 3 or len(sys.argv) > 4:
        print("Usage: python check_same.py FILE1.nc FILE2.nc [VAR1,VAR2,...]")
        sys.exit(1)

    file1 = sys.argv[1]
    file2 = sys.argv[2]

    variables_to_check = None
    if len(sys.argv) == 4:
        variables_to_check = set(sys.argv[3].split(','))

    check_same(file1, file2, variables_to_check)
