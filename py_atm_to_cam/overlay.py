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

T_rgm = f2.variables['T'][:]
Q_rgm = f2.variables['Q'][:]
U_rgm = f2.variables['U'][:]
V_rgm = f2.variables['V'][:]
PS_rgm = f2.variables['PS'][:]
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

# Append the mask to the corresponding dimension of your variables
T_rgm.mask[:, mask, :] = True
Q_rgm.mask[:, mask, :] = True
U_rgm.mask[:, mask, :] = True
V_rgm.mask[:, mask, :] = True

print(f"Loading {args.base_file}")
f1 = Dataset(args.base_file, "r+")

T_base = f1.variables['T'][:]
Q_base = f1.variables['Q'][:]
U_base = f1.variables['U'][:]
V_base = f1.variables['V'][:]
PS_base = f1.variables['PS'][:]

print("Tapering data with sigmoid function near boundary")
# Build weights matrix + expand for 3D and 2D
weights = 1 - sigmoid_weight(dists,threshold=250, scale=50)
#weights = 1 - sigmoid_weight(dists,threshold=1, scale=0.0001)
weights_2d = weights.reshape(1, -1)
weights_2d = np.where(PS_rgm.mask,1,weights_2d)
weights_3d = weights.reshape(1, 1, -1)
weights_3d = np.tile(weights_3d, (1, len(levels), 1))
weights_3d = np.where(T_rgm.mask,1,weights_3d)

# For valid rgm cells, taper the values by weighting them towards global background
# when they are near a domain boundary
T_rgm = T_base * weights_3d + T_rgm * (1 - weights_3d)
Q_rgm = Q_base * weights_3d + Q_rgm * (1 - weights_3d)
U_rgm = U_base * weights_3d + U_rgm * (1 - weights_3d)
V_rgm = V_base * weights_3d + V_rgm * (1 - weights_3d)
PS_rgm = PS_base * weights_2d + PS_rgm * (1 - weights_2d)

print("Masking data external to regional domain")
# Fully overwrite any masked values in the regional data with global data
T_base = np.where(~T_rgm.mask, T_rgm, T_base)
Q_base = np.where(~Q_rgm.mask, Q_rgm, Q_base)
U_base = np.where(~U_rgm.mask, U_rgm, U_base)
V_base = np.where(~V_rgm.mask, V_rgm, V_base)
PS_base = np.where(~PS_rgm.mask, PS_rgm, PS_base)

# Write the new variables back to the main dataset
print(f"Writing regional information to {args.base_file}")
f1.variables['T'][:] = T_base
f1.variables['Q'][:] = Q_base
f1.variables['U'][:] = U_base
f1.variables['V'][:] = V_base
f1.variables['PS'][:] = PS_base

# For debugging, either overwrite existing dists/weights or create new vars
if 'dists' in f1.variables:
  f1.variables['dists'][:] = dists
else:
  ncol_dim = f1.dimensions['ncol']
  dists_var = f1.createVariable('dists', 'f8', (ncol_dim.name,))
  dists_var[:] = dists

if 'weights' in f1.variables:
  f1.variables['weights'][:] = weights
else:
  ncol_dim = f1.dimensions['ncol']
  weights_var = f1.createVariable('weights', 'f8', (ncol_dim.name,))
  weights_var[:] = weights

# Cleanup
f1.close()
f2.close()

end_time = time.time()
elapsed_time = end_time - start_time
print(f".. DONE. Elapsed time: {elapsed_time} seconds")
