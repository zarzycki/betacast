import numpy as np
import re

def getTrajectories(filename,nVars,headerDelimStr,isUnstruc):
    """
    Read trajectory data from a TempestExtremes file format.

    Parameters:
    ----------
    filename (str): Path to input trajectory file
    nVars (int): Number of variables per trajectory point (-1 for auto-detection)
    headerDelimStr (str): String that marks header lines (e.g. "start")
    isUnstruc (bool): If True, adds an extra column for unstructured grid data

    Returns:
    ----------
    numtraj (int): Number of trajectories found
    maxNumPts (int): Maximum length of any trajectory
    ncols (int): Number of data columns per point
    prodata (np.ndarray): Array of shape (nvars, ntraj, maxpts) containing trajectory data
    """

    print("Getting trajectories from TempestExtremes file...")
    print("Running getTrajectories on '%s' with unstruc set to '%s'" % (filename, isUnstruc))
    print("nVars set to %d and headerDelimStr set to '%s'" % (nVars, headerDelimStr))

    # Using the newer with construct to close the file automatically.
    with open(filename) as f:
        data = f.readlines()

    # Find total number of trajectories and maximum length of trajectories
    numtraj=0
    numPts=[]

    for line in data:
        if headerDelimStr in line:
            # if header line, store number of points in given traj in numPts
            headArr = line.split()
            numtraj += 1
            numPts.append(int(headArr[1]))
        else:
            # if not a header line, and nVars = -1, find number of columns in data point
            if nVars < 0:
                nVars=len(line.split())

    maxNumPts = max(numPts) # Maximum length of ANY trajectory
    print("Found %d columns" % nVars)
    print("Found %d trajectories" % numtraj)

    # Initialize storm and line counter
    stormID=-1
    lineOfTraj=-1

    # Create array for data
    if isUnstruc:
        prodata = np.empty((nVars+1,numtraj,maxNumPts))
    else:
        prodata = np.empty((nVars,numtraj,maxNumPts))
    prodata[:] = np.nan

    for i, line in enumerate(data):
        if headerDelimStr in line:  # check if header string is satisfied
            stormID += 1      # increment storm
            lineOfTraj = 0    # reset trajectory line to zero
        else:
            ptArr = line.split()
            for jj in range(nVars):
                if isUnstruc:
                    prodata[jj+1,stormID,lineOfTraj]=ptArr[jj]
                else:
                    prodata[jj,stormID,lineOfTraj]=ptArr[jj]
            lineOfTraj += 1   # increment line

    # Make sure we return the correct size of the array
    if isUnstruc:
        ncols=nVars+1
    else:
        ncols=nVars

    print("... done reading data")
    return numtraj, maxNumPts, ncols, prodata


def writeTrajectories(filename, data, ntraj, npts, headerDelimStr="start"):
    """
    Write trajectory data to a file in TempestExtremes format.

    Parameters:
    ----------
    filename (str): Output file path
    data (np.ndarray): Data array of shape (nvars, ntraj, maxpts)
    ntraj (int): Number of trajectories
    npts (int): Max number of points
    headerDelimStr (str): Header delimiter string
    """

    print(f"Writing TempestExtremes file to {filename}...")
    with open(filename, 'w') as f:
        for storm in range(ntraj):
            # Find valid points for this storm
            valid_points = ~np.isnan(data[0,storm,:])
            num_valid = int(np.sum(valid_points))

            # Get time values from last 4 columns for header
            yyyy = int(data[-4,storm,0])
            mm = int(data[-3,storm,0])
            dd = int(data[-2,storm,0])
            hh = int(data[-1,storm,0])

            # Write header
            f.write(f"{headerDelimStr}\t{num_valid}\t{yyyy}\t{mm}\t{dd}\t{hh}\n")

            # Write data points
            for pt in range(num_valid):
                line_data = []
                # Process all columns except time columns
                for col in range(data.shape[0]-4):
                    val = data[col,storm,pt]
                    if col <= 1:  # First two columns are integers
                        line_data.append(f"{int(val)}")
                    else:
                        line_data.append(f"{val:e}")

                # Add time columns
                for timevar in [-4, -3, -2, -1]:
                    line_data.append(f"{int(data[timevar,storm,pt])}")

                f.write("\t" + "\t".join(line_data) + "\n")




def getNodes(filename, nVars, isUnstruc):
    """
    Read node data from a TempestExtremes nodes file format.

    Parameters:
    ----------
    filename (str): Path to input nodes file
    nVars (int): Number of variables per node (-1 for auto-detection)
    isUnstruc (bool): If True, adds extra columns for unstructured grid data

    Returns:
    ----------
    numnodetimes (int): Number of timesteps with nodes
    maxNumPts (int): Maximum number of nodes at any timestep
    prodata (np.ndarray): Array of shape (nvars+4/5, numnodetimes, maxpts) containing node data
    """

    print("Getting nodes from TempestExtremes file...")

    # Using the newer with construct to close the file automatically.
    with open(filename) as f:
        data = f.readlines()

    # Find total number of trajectories and maximum length of trajectories
    numnodetimes=0
    numPts=[]

    for line in data:
        if re.match(r'\w', line):
            # if header line, store number of points in given traj in numPts
            headArr = line.split()
            numnodetimes += 1
            numPts.append(int(headArr[3]))
        else:
            # if not a header line, and nVars = -1, find number of columns in data point
            if nVars < 0:
                nVars=len(line.split())

    maxNumPts = max(numPts) # Maximum length of ANY trajectory
    print("Found %d columns" % nVars)
    print("Found %d trajectories" % numnodetimes)
    print("Found %d maxNumPts" % maxNumPts)

    # Initialize storm and line counter
    stormID=-1
    lineOfTraj=-1

    # Create array for data
    if isUnstruc:
        prodata = np.empty((nVars+5,numnodetimes,maxNumPts))
    else:
        prodata = np.empty((nVars+4,numnodetimes,maxNumPts))

    prodata[:] = np.nan
    nextHeadLine=0

    for i, line in enumerate(data):
        if re.match(r'\w', line):  # check if header string is satisfied
            stormID += 1      # increment storm
            lineOfTraj = 0    # reset trajectory line to zero
            headArr = line.split()
            YYYY = int(headArr[0])
            MM = int(headArr[1])
            DD = int(headArr[2])
            HH = int(headArr[4])
        else:
            ptArr = line.split()
            for jj in range(nVars-1):
                if isUnstruc:
                    prodata[jj+1,stormID,lineOfTraj]=ptArr[jj]
                else:
                    prodata[jj,stormID,lineOfTraj]=ptArr[jj]
            if isUnstruc:
                prodata[nVars+1,stormID,lineOfTraj]=YYYY
                prodata[nVars+2,stormID,lineOfTraj]=MM
                prodata[nVars+3,stormID,lineOfTraj]=DD
                prodata[nVars+4,stormID,lineOfTraj]=HH
            else:
                prodata[nVars  ,stormID,lineOfTraj]=YYYY
                prodata[nVars+1,stormID,lineOfTraj]=MM
                prodata[nVars+2,stormID,lineOfTraj]=DD
                prodata[nVars+3,stormID,lineOfTraj]=HH
            lineOfTraj += 1   # increment line

    print("... done reading data")
    return numnodetimes, maxNumPts, prodata



def append_column(data, new_column, position=-5):
    """
    Add a new column to 3D numpy array at specified position.

    Parameters:
    ----------
    data (np.ndarray): Original TE 3D array of shape (nvars, ntraj, npts)
    new_column (np.ndarray): Column to add, must be shape (ntraj, npts)
    position (int): Position to insert column, counting from end. Default -5 puts column before last 4

    Returns:
    ----------
    new_data (np.ndarray): Array with new column inserted
    """
    # Check dimensions
    if new_column.shape != data.shape[1:]:
        raise ValueError(f"New column shape {new_column.shape} does not match required shape {data.shape[1:]}")

    # Create new array with space for additional column
    new_data = np.empty((data.shape[0]+1, data.shape[1], data.shape[2]))

    # Copy data before insertion point
    new_data[:position] = data[:position+1]

    # Insert new column
    new_data[position] = new_column

    # Copy remaining data after insertion point
    new_data[position+1:] = data[position+1:]

    return new_data