import numpy as np
from netCDF4 import Dataset
import sys
import datetime

# Usage: python perturb_white_noise.py /path/to/input/file.nc
if len(sys.argv) != 2:
    print("Usage: python perturb_white_noise.py /path/to/input/file.nc")
    sys.exit(1)

basFileName = sys.argv[1]

pertMag = 0.00001
print(f"We are perturbing: {basFileName}")

pertMagPerc = pertMag * 100
print(f"We are using perturbation magnitude of: {pertMagPerc}%")

# Open the NetCDF file
# Open with 'r+' which allows to read + overwrite vars (but not file)
basFile = Dataset(basFileName, 'r+')

################## Q var
# Figure out what Q var we want to perturb
potential_q_vars = ['Q','dpQ']
varq = None

for var in potential_q_vars:
    if var in basFile.variables:
        varq = var
        break

if varq is None:
    print("No Q var found!")
    basFile.close()
    sys.exit(1)

print(f"Using this Q var: {varq}")
q = basFile.variables[varq][:]

################## 2D
# Figure out what surface var we want to perturb
potential_2d_vars = ['PS', 'PS_dry', 'PSDRY']
var2d = None

for var in potential_2d_vars:
    if var in basFile.variables:
        var2d = var
        break

if var2d is None:
    print("No 2D var found!")
    basFile.close()
    sys.exit(1)

print(f"Using this 2D var: {var2d}")
ps = basFile.variables[var2d][:]

################## 3D
# Figure out what 3D var we want to perturb
potential_3d_vars = ['T', 'THETA', 'U', 'V']
var3d = None

for var in potential_3d_vars:
    if var in basFile.variables:
        var3d = var
        break

if var3d is None:
    print("No 3D var found!")
    basFile.close()
    sys.exit(1)

print(f"Using this 3D var: {var3d}")
t = basFile.variables[var3d][:]

##################

# Set the random seed
np.random.seed(int(datetime.datetime.now().timestamp()))

# Generate high and low end of perturbation range
low = 1.0 - pertMag
high = 1.0 + pertMag

# Perturb Q
print(f"orig Q max: {np.max(q)}")
q = q * np.random.uniform(low, high, q.shape)
print(f"new Q max: {np.max(q)}")

# Perturb the 3D var
print(f"orig {var3d} max: {np.max(t)}")
t = t * np.random.uniform(low, high, t.shape)
print(f"new {var3d} max: {np.max(t)}")

# Perturb PS
print(f"orig PS max: {np.max(ps)}")
ps = ps * np.random.uniform(low, high, ps.shape)
print(f"new PS max: {np.max(ps)}")

# Write perturbed variables back to basFile
basFile.variables[varq][:] = q
basFile.variables[var2d][:] = ps
basFile.variables[var3d][:] = t

basFile.close()

print("... done")
