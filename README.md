# Betacast

### General workflow
1. Create a case directory with either a supported or new grid configuration and verify that it is stable/runs given arbitrary inputs.
2. Build an analysis/reanalysis -SE weight file.
3. Set up namelist files
3a. Edit machine file for your particular system
3b. Edit namelist file for your particular use case.
3c. Edit output streams file.
4. If running historical, set dates.txt
5. Submit driver script

The first step is to create a functional F compset. Broadly, this is active atmosphere, active land, and data ocean/ice. Other configurations may work (e.g., active runoff, wave, glacier models) but have not been tested.

### 1. Create case directory and test stability.

This step just requires a built and tested case of CESM. Any source mods should be made directly to this case. It is only necessary to show the model is stable for a few days.

**CESM example:**

```
$ ./case.setup
$ ### patch land mods for restart
$ ./case.build
$ ./case.submit
```

**E3SM example:**

```
$ cd ~/E3SM/cime/scripts
$ ./create_newcase --case ~/F-betacast-FC5AV1C --compset FC5AV1C --res ne30_ne30 --mach cori-knl --project m1637
$ cd ~/F-betacast-FC5AV1C
$ ./xmlchange CHARGE_ACCOUNT=m1637
$ ./xmlchange NTASKS=-8
$ ./xmlchange NTASKS_ESP=1
$ ./case.setup
$ ### patch land mods for restart
$ ./case.build
$ ./case.submit
```

Some notes:
1. In the "### patch" step, a small modification is made to the land model to enforce restart files to be printed every 12 hours. This is done since the land model is initialized via nudging with the atmosphere in this framework. This patch can be applied by copying CESM or E3SM's `lnd_comp_mct.F90` (from the land model source code) into `$CASEDIR/SourceMods/src.clm` (or equivalent) and running
`$ patch lnd_comp_mct.F90 < ${PATCHDIR}/lnd_comp_mct.patch`
over the top of the file, which injects the correct logic. This needs to be done before the `./case.build` step.

### 2. Generate an analysis/reanalysis to CAM weight file

This can be done with `${BETACAST}/remapping/gen_GFS_to_SE_weight_file.ncl`. This script requires two inputs that are directly modified in the script body, `srcGridFile` (a SCRIP grid file defining the GFS regular lat-lon grid) and `dstGridFile` (a SCRIP grid file defining the destination CAM grid).

Historically there have been two analysis grid sizes associated with publicly disseminated GFS/CFS/CFSR analyses, 0.5degree and 0.25degree.

### 3. Edit namelists

There is a rudimentary namelist capability I have written into Bash. Ideally, this would be in Python or something else but for now this is sufficent.

A bash namelist for this codebase consists of a text file and variables of the format

`VALUE = key`

Each VALUE line is read in, splitting on ` = `.

NOTEx1: There must be a space between and after the splitting `=`!
NOTEx2: If you want to pass an empty string in, you must define the key as `"___"` (three underscores) since spaces break the splitting.
NOTEx3: Namelist files may include comments by specifying `#` as the first character of a line.

### 3a. Edit machine file for your particular system

In `${BETACAST}/machine_files` there are sample files that define where folders and data files will be stored for your system. There are suggested configurations for Cheyenne and Cori-NERSC, but you may edit these for your workflow or copy/paste for a different system (i.e., university cluster).

| Namelist Variable | Description |
| --- | --- |
| path_to_case | Path to CESM "home" case directory |
| path_to_inputdata | Path to where re/analysis + model initial conditions/forcing data is stored |
| path_to_rundir | Path (top-level) to directory where CESM actively runs |
| sewxscriptsdir | Path to betacast repo (i.e., `${BETACAST}`) |

### 3b. Edit namelist file for your particular case

In `${BETACAST}/namelist_files` there are sample files that define the forecast configuration. This is the primary location where run settings are specified.

| Namelist Variable | Description |
| --- | --- |
| debug | Setting to 1 adds debugging options. Otherwise leave at 0 |
| islive | if 1 then pull GDAS/GFS from server in real-time, 0 is "hindcast" mode |
| runmodel | Unused, set to "true" |
| atmDataType | What ATM data we want to use? 1 = GFS ANL, 2 = ERA-I, 3 = CFSR |
| sstDataType | What SST data we want to use? 1 = GDAS, 2 = ERA, 3 = NOAAOI |
| numLevels | 72 -> E3SMv1, 32 -> CAM6, 30 -> CAM5, 26 -> CAM4 |
| numdays | How long for forecast to run (in days) |
| doFilter | Should we apply offline forward DFI? Generally "true" for SE, more diffusive cores can set to "false" |
| filterOnly | Exit code after the filter run if doFilter=true (useful for producing ncdata for ensembles) |
| numHoursSEStart | Centerpoint of filter duration (leave at 3)|
| filterHourLength | Filter duration (leave at 6)|
| filtTcut | Cut setting for filter (leave at 6) |
| add_perturbs | Add PGW perturbations for counterfactual runs (leave at false) |
| add_noise | Add white noise to ncdata for ensemble (leave at false) |
| land_spinup | Cycle land spinup only (unsupported currently, leave false) |
| gfs2seWeights | Path to file allowing for GFS -> ATM regridding |
| landrawdir | For CLM5, path to CLM restart files to check/interpolate from if native grid finidat does not exist |
| PROJECTID | Project ID for run submissions |
| FILTERWALLCLOCK | Wall clock time for filter run |
| FILTERQUEUE | Submission queue for filter run |
| RUNWALLCLOCK | Wall clock time for forecast run |
| RUNQUEUE | Submission queue for forecast run |
| usingCIME | Are we using CIME (set to "true" unless using a very old CESM tag or unsupported GCM) |
| DTIME | Physics timestep (in seconds) |
| FINERES | Finest resolution of SE grid |
| USERSTAB | Required dynamics timestep (in s), negative values try internal calculation, but use with caution |
| sendplots | Are we going to send live output to some external server? (generally false unless you are CMZ) |
| nclPlotWeights | Weights to go from unstructured -> lat/lon grid for plotting (generally false unless you are CMZ) |
| dotracking | Do online TC tracking and process to ATCF format? |

### 3c. Edit output streams

```
nhtfrq=-3
mfilt=1
fincl1='PS:I',U10:I','PRECT:I'
```

### 4. Dates file

If not running in real-time, dates to be simulated are passed into the script via a text file in the root betacast directory named `dates.{CASENAME}.txt`.

For example, let's say we are running ne30np4 forecasts with a casename of MY_NE30_BETACAST and want to run a simulation on 00Z August 1, 2nd, and 3rd, 2018. In the `${BETACAST}` directory we would add a text file named `dates.MY_NE30_BETACAST.txt` and add the following lines to the top of the file:

```
2018080100
2018080200
2018080300
```

... when islive is equal to 0 in the namelist, the code will extract the YYYY, MM, DD, and HH from the 1st line of this text file (deleting it in the process -- NOTE, if a run is terminated early

### 5. Run model

```$ ./betacast.sh machine_files/machine.cheyenne namelists/nl.conus30x8 output_streams/output.generic```
