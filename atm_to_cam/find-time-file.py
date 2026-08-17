import os
import sys
import numpy as np
import xarray as xr
from datetime import datetime, timedelta
import argparse
import cftime

def build_time_lookup(filelist):
    """
    Build a lookup table for filenames and time indices.

    Parameters:
    -----------
    filelist : list of str
        List of NetCDF file paths.

    Returns:
    --------
    dict
        A dictionary with keys 'times', 'files', 'units', and 'calendar'.
    """
    nfiles = len(filelist)
    ftmp = xr.open_dataset(filelist[0], decode_times=False)
    nftimes = ftmp.time.size
    timeArr = np.empty((nfiles, nftimes), dtype=ftmp.time.dtype)
    listArr = []

    for ii, filename in enumerate(filelist):
        print(f"Loading time for file #{ii} of {nfiles-1}: {filename}")
        ftmp = xr.open_dataset(filename, decode_times=False)
        nftimes = ftmp.time.size  # Adjust in case files at the end have fewer timesteps
        timeArr[ii, :nftimes] = ftmp.time.values
        listArr.append(filename)

    timeunits = ftmp.time.attrs['units']
    calendar = ftmp.time.attrs.get('calendar', 'standard')

    return {
        'times': timeArr,
        'files': listArr,
        'units': timeunits,
        'calendar': calendar
    }

def get_time_lookup(filelist, stashedfile):
    """
    Retrieve the time lookup from a stashed file.

    Parameters:
    -----------
    filelist : list of str
        List of NetCDF file paths.

    stashedfile : str
        Path to the file containing stashed time lookup data.

    Returns:
    --------
    dict
        A dictionary with keys 'times', 'files', 'units', and 'calendar'.
    """
    ds = xr.open_dataset(stashedfile, decode_times=False)
    timeArr = ds.timeArr.values
    ds.close()
    return {
        'times': timeArr,
        'files': filelist,
        'units': ds.timeArr.units,
        'calendar': ds.timeArr.calendar
    }

def get_file_and_time_from_lookup(timeArr, thisTime, eps=1.0e-6):
    """
    Find the file and time index closest to the specified time.

    Parameters:
    -----------
    timeArr : dict
        Dictionary returned by build_time_lookup or get_time_lookup.

    thisTime : float
        The desired time in the units of the time array.

    eps : float, optional
        Tolerance for time matching.

    Returns:
    --------
    tuple
        (filename, time index, error code)
    """
    time1d = timeArr['times'].flatten()
    timediff = np.abs(time1d - thisTime)
    min_diff = np.min(timediff)
    indices = np.unravel_index(np.argmin(timediff), timeArr['times'].shape)

    if min_diff > eps:
        print(f"WARNING: diff b/w best matched time in nc files & time from traj file = {min_diff}")
        print(f"WARNING: thisTime is: {thisTime}")
        return_error = int(min_diff)
    else:
        return_error = 0

    if np.any(np.isnan(indices)):
        trackindex = -1
        needed_file = ""
    else:
        fileix = indices[0]
        trackindex = indices[1]
        needed_file = timeArr['files'][fileix]

    return needed_file, trackindex, return_error

def find_time_in_files(DIR, YYYYMMDDHH, UQSTR):
    """
    Find the file containing the specified time.

    Parameters:
    -----------
    DIR : str
        Directory containing the NetCDF files.

    YYYYMMDDHH : str
        The desired time in the format YYYYMMDDHH.

    UQSTR : str
        A unique string to append to the output file.

    Returns:
    --------
    None
    """
    if not os.path.exists(DIR):
        print("find-time-file.py: user did not define a valid DIR, exiting")
        sys.exit(1)

    print(f"Finding file in: {DIR}")
    filelist = sorted([os.path.join(DIR, f) for f in os.listdir(DIR) if f.endswith(".nc")])

    if not filelist:
        print("No NetCDF files found in the directory, exiting")
        sys.exit(1)

    timeArr = build_time_lookup(filelist)

    # Convert YYYYMMDDHH to a datetime object
    yyyy = int(YYYYMMDDHH[:4])
    mm = int(YYYYMMDDHH[4:6])
    dd = int(YYYYMMDDHH[6:8])
    hh = int(YYYYMMDDHH[8:])

    print(f"Finding time: {yyyy} {mm} {dd} {hh} 0 0")

    # Convert to datetime
    desired_time = datetime(yyyy, mm, dd, hh)

    # Convert to a cftime object
    cftime_value = cftime.date2num(desired_time, units=timeArr['units'], calendar=timeArr['calendar'])

    # Cast to numpy float64
    desired_time_numeric = np.float64(cftime_value)

    # Find the file and time index
    theFile, trackindex, return_error = get_file_and_time_from_lookup(timeArr, desired_time_numeric)

    if theFile and os.path.exists(theFile):
        print(f"Found {YYYYMMDDHH} in file: {theFile}")
        with open(f"m2mfile.{UQSTR}", 'w') as f:
            f.write(theFile)
    else:
        print(f"find-time-file.py: Unable to correctly find file containing {YYYYMMDDHH}")
        sys.exit(1)

    # Exit cleanly
    sys.exit(0)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Find the file containing the specified time.")
    parser.add_argument("--DIR", type=str, required=True, help="Directory containing NetCDF files")
    parser.add_argument("--YYYYMMDDHH", type=str, required=True, help="Desired time in format YYYYMMDDHH")
    parser.add_argument("--UQSTR", type=str, required=True, help="Unique string for output file")

    args = parser.parse_args()

    find_time_in_files(args.DIR, args.YYYYMMDDHH, args.UQSTR)

