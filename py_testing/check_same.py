import sys
import xarray as xr
import numpy as np

def calculate_mean_absolute_error(var1, var2):
    return np.mean(np.abs(var1 - var2))

def calculate_correlation(var1, var2):
    return np.corrcoef(var1.ravel(), var2.ravel())[0, 1]

def calculate_normalized_mean_bias(var1, var2):
    mean_var1 = np.mean(var1)
    if np.isnan(mean_var1) or mean_var1 == 0:
        print("   Warning: Mean of var1 is NaN or zero, NMB calculation skipped.")
        return np.nan
    return np.mean(var2 - var1) / mean_var1

def calculate_normalized_rmse(var1, var2):
    rmse = np.sqrt(np.mean((var1 - var2) ** 2))
    mean_var1 = np.mean(var1)
    if np.isnan(mean_var1) or mean_var1 == 0:
        print("   Warning: Mean of var1 is NaN or zero, NRMSE calculation skipped.")
        return np.nan
    return rmse / mean_var1

def check_same(file1, file2, variables_to_check=None):
    ds1 = xr.open_dataset(file1)
    ds2 = xr.open_dataset(file2)

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

                mae = calculate_mean_absolute_error(data1, data2)
                correlation = calculate_correlation(data1, data2)
                mean_value = np.mean(data1)
                nmb = calculate_normalized_mean_bias(data1, data2)
                nrmse = calculate_normalized_rmse(data1, data2)

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
