import numpy as np
import os
import subprocess
import time
import sys
import argparse
from datetime import datetime

# Betacast modules
module_path = os.path.abspath(os.path.join(os.path.dirname(__file__), '../..', 'py_functions'))
if module_path not in sys.path:
    sys.path.append(module_path)
from py_seedfuncs import keyword_values, gc_latlon

# Constants
minlat = 5.
maxlat = 15.
minlon = 0.0000001
maxlon = 359.99999
minrmw = 150000.
maxrmw = 450000.
mindp = 15.0
maxdp = 40.0
THRESHOLD = 15.

# Set up argument parser
parser = argparse.ArgumentParser(description='Process cyclone seeding parameters')
parser.add_argument('--filename', type=str, required=True, help='Path to the cyclones tempest file')
parser.add_argument('--pthi', type=str, required=True, help='Path to the namelist file')
args = parser.parse_args()

# Get command line arguments
filename = args.filename
pthi = args.pthi

# Check invert_vortex
invert_vortex = keyword_values(pthi, "invert_vortex", "bool")
if invert_vortex:
    print(f"invert_vortex is set to {invert_vortex}")
    print("We don't need to worry about where to seed vortices!")
    sys.exit(9)

if os.path.exists(filename):
    # Read existing tempest file
    with open(filename, 'r') as f:
        lines = f.readlines()
    nlines = len(lines)
    print(f"nlines: {nlines}")
    print("-------------------------------------")

    # Parse data
    data = np.genfromtxt(filename, usecols=(1, 2))  # lon, lat columns
    # Handle single line case
    if data.ndim == 1:
        data = data.reshape(1, -1)
    print(data)
    lons = data[:, 0]
    lats = data[:, 1]

    print("lats        lons")
    for lat, lon in zip(lats, lons):
        print(f"{lat:10.5f}     {lon:10.5f}")
    print("-------------------------------------")

    # Set random seed using current time
    np.random.seed(int(time.time_ns() % 1000000000))

    # Find suitable location
    DONE = False
    ITERMAX = 100
    ii = 0
    while not DONE:
        ii += 1
        thislat = np.random.uniform(minlat, maxlat)
        thislon = np.random.uniform(minlon, maxlon)
        print(f"lat/lon: {thislat} {thislon}")

        gcdist, _ = gc_latlon(thislat, thislon, lats, lons, 2, 2)
        min_dist = np.min(gcdist)
        print(f"minimum separation is: {min_dist}")

        if min_dist > THRESHOLD:
            DONE = True

        if ii >= ITERMAX:
            print(f"at iter {ii} could not find random point outside of THRESHOLD {THRESHOLD}")
            sys.exit(1)
else:
    # No existing TCs
    print("random_seed: no preexisting vortices!")
    thislat = np.random.uniform(minlat, maxlat)
    thislon = np.random.uniform(minlon, maxlon)

print(f"random seedlat: {thislat}")
print(f"random seedlon: {thislon}")

target_rmw = np.random.uniform(minrmw, maxrmw)
minp = np.random.uniform(mindp, maxdp)
print(f"random target_rmw: {target_rmw}")
print(f"random minp: {minp}")
inv_minp = -1 * minp * 100

print("updating namelist file")
# Use sed to update the namelist file
subprocess.run(["sed", "-i", f"s?.*psminlat=.*?psminlat={thislat}?", pthi])
subprocess.run(["sed", "-i", f"s?.*psminlon=.*?psminlon={thislon}?", pthi])
subprocess.run(["sed", "-i", f"s?.*target_rmw=.*?target_rmw={target_rmw}?", pthi])
subprocess.run(["sed", "-i", f"s?.*minp=.*?minp={inv_minp}?", pthi])

# Append new storm to tracking file
with open(filename, 'a') as f:
    f.write(f"     9999    {thislon}     {thislat}      9.99999e+11\n")

print("DONE figuring out seeds!")