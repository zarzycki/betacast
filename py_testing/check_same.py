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
    if valid_mask.sum() <= 1:
        print("   Warning: Insufficient valid data points for correlation calculation.")
        return np.nan

    # Check if any variable has constant values (zero variance)
    var1_valid = var1[valid_mask]
    var2_valid = var2[valid_mask]

    var1_std = np.nanstd(var1_valid)
    var2_std = np.nanstd(var2_valid)

    if var1_std < 1e-10 or var2_std < 1e-10:
        print(f"   Note: At least one variable has constant values (std dev1: {var1_std:.1e}, std dev2: {var2_std:.1e}).")
        # Check if they're identical to each other
        if np.allclose(var1_valid, var2_valid, rtol=1e-10, atol=1e-10):
            print("   Note: Both variables have effectively identical values.")
            return 1.0  # Perfect correlation for identical values
        else:
            print("   Note: Variables have different values.")
            return np.nan  # Correlation is undefined for constant variables

    # Safe correlation calculation
    try:
        with np.errstate(invalid='ignore', divide='ignore'):
            corr = np.corrcoef(var1_valid.ravel(), var2_valid.ravel())[0, 1]
        return corr
    except Exception as e:
        print(f"   Warning: Correlation calculation failed: {str(e)}")
        return np.nan

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

def check_same(file1, file2, variables_to_check=None,
               corr_threshold_strict=0.999, corr_threshold_normal=0.99,
               nrmse_threshold_strict=0.1, nrmse_threshold_normal=1.0):
    """Returns True if files match within acceptable tolerances, False otherwise."""
    ds1 = xr.open_dataset(file1)
    ds2 = xr.open_dataset(file2)

    # Print file format information
    print(f"\nFile 1: {file1}")
    print(f"File 2: {file2}")

    # Get NetCDF format information (similar to ncdump -k)
    try:
        import netCDF4
        with netCDF4.Dataset(file1) as nc:
            nc_format1 = nc.file_format
        print(f"File 1 NetCDF format: {nc_format1}")
    except Exception as e:
        print(f"File 1 NetCDF format: Unable to determine - {str(e)}")

    try:
        import netCDF4
        with netCDF4.Dataset(file2) as nc:
            nc_format2 = nc.file_format
        print(f"File 2 NetCDF format: {nc_format2}")
    except Exception as e:
        print(f"File 2 NetCDF format: Unable to determine - {str(e)}")

    # Print dataset encoding information
    print("\nFile 1 encoding information:")
    if hasattr(ds1, 'encoding') and ds1.encoding:
        for key, value in ds1.encoding.items():
            print(f"   {key}: {value}")
    else:
        print("   No encoding information available")

    print("\nFile 2 encoding information:")
    if hasattr(ds2, 'encoding') and ds2.encoding:
        for key, value in ds2.encoding.items():
            print(f"   {key}: {value}")
    else:
        print("   No encoding information available")

    common_vars = set(ds1.variables).intersection(ds2.variables)
    if variables_to_check:
        print(f"\nChecking only these specific variables: {', '.join(sorted(variables_to_check))}")
        common_vars = common_vars.intersection(variables_to_check)
        if not common_vars:
            print("Warning: None of the specified variables were found in both files!")
    else:
        print("\nChecking all common variables between files")

    all_vars_ok = True

    # Threshold info
    print(f"\nUsing thresholds:")
    print(f"   Correlation - Strict variables: {corr_threshold_strict}, Normal variables: {corr_threshold_normal}")
    print(f"   NRMSE - Strict variables: {nrmse_threshold_strict}, Normal variables: {nrmse_threshold_normal}")

    for var_name in common_vars:

        print("-----------------------------------------------------------------------")
        print(f"----------- {var_name}")

        # Check data types
        dtype1 = ds1[var_name].dtype
        dtype2 = ds2[var_name].dtype

        if dtype1 != dtype2:
            print(f"WARNING: Variable {var_name} has different data types: {dtype1} vs {dtype2}")
            if (np.issubdtype(dtype1, np.integer) and np.issubdtype(dtype2, np.floating)) or \
               (np.issubdtype(dtype1, np.floating) and np.issubdtype(dtype2, np.integer)):
                print(f"   NOTE: Mixed integer and floating-point types!")

        data1 = ds1[var_name].values
        data2 = ds2[var_name].values

        # Check if data consists entirely of "special" values (zeros or nans)
        all_zeros1 = np.all(data1 == 0)
        all_zeros2 = np.all(data2 == 0)
        all_nans1 = np.all(np.isnan(data1)) if np.issubdtype(data1.dtype, np.floating) else False
        all_nans2 = np.all(np.isnan(data2)) if np.issubdtype(data2.dtype, np.floating) else False

        if all_zeros1 and all_zeros2:
            print("   Note: Both arrays consist entirely of zeros")
            # Skip numerical comparison for all-zero arrays
            print("   Correlation: 1.0 (identical zeros)")
            print("   Mean Absolute Error: 0.0")
            print("   Mean Value: 0.0")
            print("   Normalized Mean Bias (NMB): 0.0")
            print("   Normalized RMSE (NRMSE): 0.0")
            continue

        if all_nans1 and all_nans2:
            print("   Note: Both arrays consist entirely of NaN values")
            # Skip numerical comparison for all-NaN arrays
            print("   Correlation: 1.0 (identical NaNs)")
            print("   Mean Absolute Error: NaN")
            print("   Mean Value: NaN")
            print("   Normalized Mean Bias (NMB): NaN")
            print("   Normalized RMSE (NRMSE): NaN")
            continue

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

                strict_vars = {'PS', 'U', 'V', 'T', 'Q'}

                corr_threshold = corr_threshold_strict if var_name in strict_vars else corr_threshold_normal
                nrmse_threshold = nrmse_threshold_strict if var_name in strict_vars else nrmse_threshold_normal

                # Check if this variable passes validation
                var_ok = True
                if correlation < corr_threshold:
                   print(f"******************** Correlation: {correlation} (threshold: {corr_threshold})")
                   var_ok = False
                else:
                   print(f"   Correlation: {correlation}")

                print(f"   Mean Absolute Error: {mae}")
                print(f"   Mean Value: {mean_value}")
                print(f"   Normalized Mean Bias (NMB): {nmb}")

                if np.abs(nrmse) > nrmse_threshold:
                   print(f"******************** Normalized RMSE (NRMSE): {nrmse} (threshold: {nrmse_threshold})")
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
    parser.add_argument('file1',
                        help='First NetCDF file')
    parser.add_argument('file2',
                        help='Second NetCDF file')
    parser.add_argument('variables', nargs='?',
                        help='Comma-separated list of variables to compare')
    parser.add_argument('--corr-threshold-strict', type=float, default=0.999,
                        help='Correlation threshold for strict variables (PS, U, V, T, Q). Default: 0.999')
    parser.add_argument('--corr-threshold-normal', type=float, default=0.99,
                        help='Correlation threshold for normal variables. Default: 0.99')
    parser.add_argument('--nrmse-threshold-strict', type=float, default=0.1,
                        help='NRMSE threshold for strict variables (PS, U, V, T, Q). Default: 0.1')
    parser.add_argument('--nrmse-threshold-normal', type=float, default=1.0,
                        help='NRMSE threshold for normal variables. Default: 1.0')

    args = parser.parse_args()

    variables_to_check = None
    if args.variables:
        variables_to_check = set(args.variables.split(','))

    try:
        match = check_same(args.file1, args.file2, variables_to_check,
                  corr_threshold_strict=args.corr_threshold_strict,
                  corr_threshold_normal=args.corr_threshold_normal,
                  nrmse_threshold_strict=args.nrmse_threshold_strict,
                  nrmse_threshold_normal=args.nrmse_threshold_normal)
        # Exit with 0 if files match, 1 if they don't
        sys.exit(0 if match else 1)
    except Exception as e:
        print(f"Error: {str(e)}")
        sys.exit(2)  # Exit with 2 for errors
