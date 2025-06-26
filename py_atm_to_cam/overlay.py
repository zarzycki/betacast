import numpy as np
from netCDF4 import Dataset
from sklearn.neighbors import BallTree
import argparse
import time

start_time = time.time()

## Sigmoid weighting function to slowly transition from global to background grid
def sigmoid_weight(distance, threshold=100, scale=10):
    # normalize distance so that the sigmoid transitions around the threshold
    normalized_distance = (distance - threshold) / scale
    return 1 / (1 + np.exp(-normalized_distance))

def detect_dycore(dataset):
    """Detect if the dataset uses SE or SCREAM variable naming convention"""
    variables = list(dataset.variables.keys())

    # Check for SCREAM-specific variables
    if 'T_mid' in variables or 'qv' in variables or 'horiz_winds' in variables:
        return 'scream'
    # Check for SE-specific variables
    elif 'T' in variables or 'Q' in variables or 'U' in variables:
        return 'se'
    else:
        raise ValueError("Cannot determine dycore - no recognizable variable pattern found")

def load_variable(dataset, dycore, var_name):
    """Load a single variable accounting for dycore differences"""
    if dycore == 'scream':
        if var_name == 'PS':
            return dataset.variables['ps'][:]  # PS is 2D, no transpose needed
        elif var_name == 'T':
            return np.transpose(dataset.variables['T_mid'][:], (0, 2, 1))
        elif var_name == 'Q':
            return np.transpose(dataset.variables['qv'][:], (0, 2, 1))
        elif var_name == 'U':
            horiz_winds = dataset.variables['horiz_winds'][:]
            u_component = horiz_winds[:, :, 0, :]
            del horiz_winds  # Free memory immediately
            return np.transpose(u_component, (0, 2, 1))
        elif var_name == 'V':
            horiz_winds = dataset.variables['horiz_winds'][:]
            v_component = horiz_winds[:, :, 1, :]
            del horiz_winds  # Free memory immediately
            return np.transpose(v_component, (0, 2, 1))
    else:  # SE
        var_names_se = {'T': 'T', 'Q': 'Q', 'U': 'U', 'V': 'V', 'PS': 'PS'}
        return dataset.variables[var_names_se[var_name]][:]

def save_variable(dataset, dycore, var_name, data):
    """Save a single variable accounting for dycore differences"""
    if dycore == 'scream':
        if var_name == 'PS':
            dataset.variables['ps'][:] = data  # PS is 2D, no transpose needed
        elif var_name == 'T':
            dataset.variables['T_mid'][:] = np.transpose(data, (0, 2, 1))
        elif var_name == 'Q':
            dataset.variables['qv'][:] = np.transpose(data, (0, 2, 1))
        elif var_name == 'U':
            # Read existing horiz_winds, update U component, write back
            horiz_winds = dataset.variables['horiz_winds'][:]
            horiz_winds[:, :, 0, :] = np.transpose(data, (0, 2, 1))
            dataset.variables['horiz_winds'][:] = horiz_winds
            del horiz_winds
        elif var_name == 'V':
            # Read existing horiz_winds, update V component, write back
            horiz_winds = dataset.variables['horiz_winds'][:]
            horiz_winds[:, :, 1, :] = np.transpose(data, (0, 2, 1))
            dataset.variables['horiz_winds'][:] = horiz_winds
            del horiz_winds
    else:  # SE
        var_names_se = {'T': 'T', 'Q': 'Q', 'U': 'U', 'V': 'V', 'PS': 'PS'}
        dataset.variables[var_names_se[var_name]][:] = data

# Argument parser
parser = argparse.ArgumentParser(description="Process base_file and top_file")
parser.add_argument("base_file", type=str, help="The base file path")
parser.add_argument("top_file", type=str, help="The top file path")
parser.add_argument("--maxLev", type=float, default=100., help="The max level (optional)")
args = parser.parse_args()
print(f"Base_file is: {args.base_file}")
print(f"top_file is: {args.top_file}")
print(f"maxLev is: {args.maxLev}")

# Constants
R = 6371.0  # Radius of the earth in km

# Load the regional file
print(f"Loading {args.top_file}")
f2 = Dataset(args.top_file, "r")

# Detect dycore and get variable names
dycore_top = detect_dycore(f2)
print(f"Detected dycore for top file: {dycore_top}")

# Load PS first to calculate distances (needed for all variables)
PS_rgm = load_variable(f2, dycore_top, 'PS')

lat = f2.variables['lat'][:]
lon = f2.variables['lon'][:]
ncol = lat.size

print(f"Number of missing values in top_file is: {PS_rgm.mask.sum()}")

# Find distances in "valid" regional cells to nearest missing cell
# in other words, how many km are we from an edge

print("Finding distances for valid cells -> missing cells... ")

# Empty dists array of length ncol
dists = np.zeros_like(lat)

# Find indices of missing and non/missing data by using PS_rgm array (which is ncol)
not_missing_inds = np.where(~PS_rgm[0, :].mask)
yes_missing_inds = np.where(PS_rgm[0, :].mask)

# Prepare data
not_missing_points = np.stack((lat[not_missing_inds], lon[not_missing_inds]), axis=-1)
yes_missing_points = np.stack((lat[yes_missing_inds], lon[yes_missing_inds]), axis=-1)

# Build the tree (using haversine metric)
tree = BallTree(np.radians(yes_missing_points), metric='haversine')

# Query the tree for nearest neighbor
dist, _ = tree.query(np.radians(not_missing_points), k=1)  # dist is in radians

# Convert distance from radians to km
dist_kms = dist * R

# Fill the calculated distances in the appropriate indices
dists[not_missing_inds] = np.squeeze(dist_kms)

print("... done finding distances")

# This only overlays data in the troposphere because regional models are generally low-topped
levels = f2.variables['lev'][:]  # Assuming you have this level array in your file
mask = levels <= args.maxLev  # Create a mask where level is less or equal to 100mb

print(f"Loading {args.base_file}")
f1 = Dataset(args.base_file, "r+")

# Detect dycore for base file
dycore_base = detect_dycore(f1)
print(f"Detected dycore for base file: {dycore_base}")

# Check if dycores match
if dycore_top != dycore_base:
    raise ValueError(f"Dycore mismatch: top file is {dycore_top}, base file is {dycore_base}")

print("Tapering data with sigmoid function near boundary")
# Build weights matrix + expand for 3D and 2D
weights = 1 - sigmoid_weight(dists,threshold=250, scale=50)
#weights = 1 - sigmoid_weight(dists,threshold=1, scale=0.0001)
weights_2d = weights.reshape(1, -1)
weights_2d = np.where(PS_rgm.mask,1,weights_2d)
weights_3d = weights.reshape(1, 1, -1)
weights_3d = np.tile(weights_3d, (1, len(levels), 1))

print("Masking data external to regional domain")

# Process each variable separately to save memory
variables_to_process = ['T', 'Q', 'U', 'V', 'PS']

print(f"Writing regional information to {args.base_file}")
for var_name in variables_to_process:
    print(f"Processing {var_name}...")

    # Load regional variable
    var_rgm = load_variable(f2, dycore_top, var_name)

    # Load base variable
    var_base = load_variable(f1, dycore_base, var_name)

    if var_name == 'PS':
        # 2D variable processing
        # For valid rgm cells, taper the values by weighting them towards global background
        var_rgm = var_base * weights_2d + var_rgm * (1 - weights_2d)
        # Fully overwrite any masked values in the regional data with global data
        var_final = np.where(~var_rgm.mask, var_rgm, var_base)
    else:
        # 3D variable processing
        # Apply level mask first
        var_rgm.mask[:, mask, :] = True

        # Create weights for this variable
        weights_3d_var = np.where(var_rgm.mask, 1, weights_3d)

        # For valid rgm cells, taper the values by weighting them towards global background
        var_rgm = var_base * weights_3d_var + var_rgm * (1 - weights_3d_var)
        # Fully overwrite any masked values in the regional data with global data
        var_final = np.where(~var_rgm.mask, var_rgm, var_base)

    # Save the processed variable
    save_variable(f1, dycore_base, var_name, var_final)

    f1.sync()  # Force write to disk

    # Clean up memory
    del var_rgm, var_base, var_final
    if var_name != 'PS':
        del weights_3d_var

# For debugging, either overwrite existing dists/weights or create new vars
print("doing dists")
if 'dists' in f1.variables:
    f1.variables['dists'][:] = dists
else:
    ncol_dim = f1.dimensions['ncol']
    dists_var = f1.createVariable('dists', 'f8', (ncol_dim.name,))
    dists_var[:] = dists

print("doing weights")
if 'weights' in f1.variables:
    f1.variables['weights'][:] = weights
else:
    ncol_dim = f1.dimensions['ncol']
    weights_var = f1.createVariable('weights', 'f8', (ncol_dim.name,))
    weights_var[:] = weights

print("Processing complete, preparing to close files...")

print("Final sync of output file...")
f1.sync()

print("Closing output file (this may take a while for large files)...")
f1.close()

print("Closing input file...")
f2.close()

print("Files closed successfully")

end_time = time.time()
elapsed_time = end_time - start_time
print(f".. DONE. Elapsed time: {elapsed_time} seconds")