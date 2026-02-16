"""
vertical_blend.py -- Vertically blend two atmospheric initial-condition files.

Usage:
    python vertical_blend.py <file_a> <file_b> [options]

Arguments:
    file_a          Base file (modified in place). Dominates below the blend level.
    file_b          Upper-atmosphere source (read-only). Dominates above the blend level.

Options:
    --blendLev      Pressure level (hPa) where the blend is centered. Default: 100
    --taperScale    Sigmoid scale in log-pressure space. Controls taper sharpness.
                    Smaller = sharper transition. Default: 0.5
                    Suggestions: 0.25 (tight), 0.5 (moderate), 1.0 (gradual)
    --vars          Comma-separated 3D variables to blend. Default: T,Q,U,V
                    PS is always skipped (surface field, no vertical dimension).
    --invert-blend  Flip the blend direction so file_b dominates below the blend
                    level and file_a dominates above it.

Notes:
    - Both files must be on the same grid with the same vertical levels.
    - The sigmoid operates in log-pressure space so the taper is symmetric on a
      log scale (physically more natural than linear pressure).
    - Supports SE and SCREAM dycore conventions (auto-detected).
    - Blend weights are stored in file_a as 'blend_weights' for diagnostics.
    - Per-level weights are printed at startup so you can verify the taper.
"""

import numpy as np
from netCDF4 import Dataset
import argparse
import time
import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'py_functions'))
from pyfuncs import detect_dycore_nc, load_variable_blending, save_variable_blending

start_time = time.time()

# Sigmoid weighting function for vertical tapering
def sigmoid_weight_vert(pressure, center, scale):
    """Returns 0 deep below center, 1 well above center (in log-pressure)."""
    # Work in log-pressure so the taper is symmetric on a log scale
    normalized = (np.log(center) - np.log(pressure)) / scale
    return 1 / (1 + np.exp(-normalized))

# Argument parser
parser = argparse.ArgumentParser(
    description="Vertically blend two atmospheric initial-condition files. "
                "Below the blend level file_a dominates; above it file_b dominates, "
                "with a smooth sigmoid taper in between. file_a is modified in place."
)
parser.add_argument("file_a", type=str,
                    help="Base file (modified in place). Dominates below blend level.")
parser.add_argument("file_b", type=str,
                    help="Upper-atmosphere source (read-only). Dominates above blend level.")
parser.add_argument("--blendLev", type=float, default=100.0,
                    help="Pressure level (hPa) at which the blend is centered (default: 100)")
parser.add_argument("--taperScale", type=float, default=0.5,
                    help="Sigmoid scale in log-pressure space. Smaller = sharper transition "
                         "(default: 0.5)")
parser.add_argument("--vars", type=str, default="T,Q,U,V",
                    help="Comma-separated list of variables to blend (default: T,Q,U,V). "
                         "PS is a surface field and is not blended.")
parser.add_argument("--invert-blend", action="store_true", default=False,
                    help="Flip blend direction: file_b dominates below blend level, "
                         "file_a dominates above.")
args = parser.parse_args()

var_list = [v.strip() for v in args.vars.split(',')]

print(f"file_a (base) : {args.file_a}")
print(f"file_b (upper): {args.file_b}")
print(f"blendLev      : {args.blendLev} hPa")
print(f"taperScale    : {args.taperScale}")
print(f"invert-blend  : {args.invert_blend}")
print(f"variables     : {var_list}")

# Open files and detect dycores
print(f"\nLoading {args.file_b}")
fb = Dataset(args.file_b, "r")
dycore_b = detect_dycore_nc(fb)
print(f"Detected dycore for file_b: {dycore_b}")

print(f"Loading {args.file_a}")
fa = Dataset(args.file_a, "r+")
dycore_a = detect_dycore_nc(fa)
print(f"Detected dycore for file_a: {dycore_a}")

if dycore_a != dycore_b:
    raise ValueError(f"Dycore mismatch: file_a is {dycore_a}, file_b is {dycore_b}")

# Build vertical weights from the level coordinate
levels = fa.variables['lev'][:]  # pressure levels in hPa
nlev = len(levels)

# w=0 deep below blendLev (file_a kept), w=1 well above blendLev (file_b used)
w_1d = sigmoid_weight_vert(levels, center=args.blendLev, scale=args.taperScale)

if args.invert_blend:
    print("\n** invert-blend is ON: flipping weights (file_b dominates below blend level) **")
    w_1d = 1 - w_1d

print(f"\nVertical weights (file_b fraction) at each level:")
for k in range(nlev):
    tag = ""
    if w_1d[k] < 0.01:
        tag = "  <- file_a"
    elif w_1d[k] > 0.99:
        tag = "  <- file_b"
    print(f"  lev {k:3d}  {levels[k]:10.4f} hPa  w={w_1d[k]:.4f}{tag}")

# Reshape for broadcasting: (1, nlev, 1) for 3D fields
w_3d = w_1d.reshape(1, nlev, 1)

# Blend each variable
print(f"\nBlending variables and writing to {args.file_a}")
for var_name in var_list:
    if var_name == 'PS':
        print(f"Skipping PS (surface field, no vertical dimension)")
        continue

    print(f"Processing {var_name}...")
    var_a = load_variable_blending(fa, dycore_a, var_name)
    var_b = load_variable_blending(fb, dycore_b, var_name)

    blended = var_a * (1 - w_3d) + var_b * w_3d
    save_variable_blending(fa, dycore_a, var_name, blended)
    fa.sync()

    del var_a, var_b, blended

# Store blend metadata for diagnostics
if 'blend_weights' in fa.variables:
    fa.variables['blend_weights'][:] = w_1d
else:
    lev_dim = fa.dimensions['lev']
    bw = fa.createVariable('blend_weights', 'f8', (lev_dim.name,))
    bw[:] = w_1d
    bw.long_name = "vertical blend weight (0=file_a, 1=file_b)"
    bw.blend_level_hPa = args.blendLev
    bw.taper_scale = args.taperScale

fa.sync()

# Clean up
print("\nClosing files...")
fa.close()
fb.close()

end_time = time.time()
print(f".. DONE. Elapsed time: {end_time - start_time:.1f} seconds")
