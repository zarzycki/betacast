import sys
import xarray as xr
import numpy as np
import argparse

def calculate_mean_absolute_error(var1, var2):
    return np.nanmean(np.abs(var1 - var2))

def calculate_correlation(var1, var2):
    # Check for differing NaN locations
    nan_mismatch = np.isnan(var1) != np.isnan(var2)
    if np.any(nan_mismatch):
        print("   Warning: NaN locations differ between var1 and var2.")

    # Only calculate correlation where both var1 and var2 are not NaN
    valid_mask = ~np.isnan(var1) & ~np.isnan(var2)
    if valid_mask.sum() == 0:
        print("   Warning: No valid data points for correlation calculation.")
        return np.nan
    return np.corrcoef(var1[valid_mask].ravel(), var2[valid_mask].ravel())[0, 1]

def calculate_normalized_mean_bias(var1, var2):
    mean_var1 = np.nanmean(var1)
    if np.isnan(mean_var1) or mean_var1 == 0:
        print("   Warning: Mean of var1 is NaN or zero, NMB calculation skipped.")
        return np.nan
    return np.nanmean(var2 - var1) / mean_var1

def calculate_normalized_rmse(var1, var2):
    rmse = np.sqrt(np.nanmean((var1 - var2) ** 2))
    mean_var1 = np.nanmean(var1)
    if np.isnan(mean_var1) or mean_var1 == 0:
        print("   Warning: Mean of var1 is NaN or zero, NRMSE calculation skipped.")
        return np.nan
    return rmse / mean_var1

def check_same(file1, file2, variables_to_check=None):
    """Returns True if files match within acceptable tolerances, False otherwise."""
    ds1 = xr.open_dataset(file1)
    ds2 = xr.open_dataset(file2)

    common_vars = set(ds1.variables).intersection(ds2.variables)
    if variables_to_check:
        common_vars = common_vars.intersection(variables_to_check)

    all_vars_ok = True

    for var_name in common_vars:
        print("-----------------------------------------------------------------------")
        data1 = ds1[var_name].values
        data2 = ds2[var_name].values

        # Convert to float64 if one is float32 and the other is float64
        if (data1.dtype == np.float32 and data2.dtype == np.float64) or (data1.dtype == np.float64 and data2.dtype == np.float32):
            print(f"Converting both {var_name} to float64")
            data1 = data1.astype(np.float64)
            data2 = data2.astype(np.float64)

        # Skip comparison if data types are not compatible
        if data1.dtype != data2.dtype:
            print(f"Skipping variable {var_name} due to incompatible data types: {data1.dtype} vs {data2.dtype}")
            all_vars_ok = False
            continue

        if np.issubdtype(data1.dtype, np.number):
            # Only perform calculations for numeric types
            if data1.shape == data2.shape:
                print(f"Variable: {var_name}")

                mae = calculate_mean_absolute_error(data1, data2)
                correlation = calculate_correlation(data1, data2)
                mean_value = np.nanmean(data1)
                nmb = calculate_normalized_mean_bias(data1, data2)
                nrmse = calculate_normalized_rmse(data1, data2)

                # Check if this variable passes validation
                var_ok = True
                if correlation < 0.99:
                    print(f"******************** Correlation: {correlation}")
                    var_ok = False
                else:
                    print(f"   Correlation: {correlation}")

                print(f"   Mean Absolute Error: {mae}")
                print(f"   Mean Value: {mean_value}")
                print(f"   Normalized Mean Bias (NMB): {nmb}")

                if np.abs(nrmse) > 1.0:
                    print(f"******************** Normalized RMSE (NRMSE): {nrmse}")
                    var_ok = False
                else:
                    print(f"   Normalized RMSE (NRMSE): {nrmse}")

                all_vars_ok = all_vars_ok and var_ok

            else:
                print(f"Variable: {var_name} has different shapes: {data1.shape} vs {data2.shape}")
                all_vars_ok = False
        else:
            print(f"Skipping non-numeric variable {var_name} with dtype {data1.dtype}")

    ds1.close()
    ds2.close()
    return all_vars_ok

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Compare two NetCDF files')
    parser.add_argument('file1', help='First NetCDF file')
    parser.add_argument('file2', help='Second NetCDF file')
    parser.add_argument('variables', nargs='?', help='Comma-separated list of variables to compare')

    args = parser.parse_args()

    variables_to_check = None
    if args.variables:
        variables_to_check = set(args.variables.split(','))

    try:
        match = check_same(args.file1, args.file2, variables_to_check)
        # Exit with 0 if files match, 1 if they don't
        sys.exit(0 if match else 1)
    except Exception as e:
        print(f"Error: {str(e)}")
        sys.exit(2)  # Exit with 2 for errors