# Betacast

**Betacast** (originally a portmanteau of "beta" and "forecast") is a software package designed to easily initialize CIME-ized modeling systems in "numerical weather prediction" mode. The code is capable of running in realtime or hindcast mode and possesses tools for generating balanced initial conditions and configuring the model appropriately. **Please cite the below manuscript and this repository if you use Betacast in your research.**

Reference: C. M. Zarzycki and C. Jablonowski. Experimental tropical cyclone forecasts using a variable-resolution global model. *Monthly Weather Review*, 1430 (10):0 4012--4037, 2015. 10.1175/MWR-D-15-0159.1.

### General workflow
1. Create a case directory with either a supported or new grid configuration and verify that it is stable/runs given arbitrary inputs.
2. Build an analysis/reanalysis grid to model grid weight file using ESMF (or TempestRemap).
3. Set up namelist files.
3a. Edit machine file for your particular system.
3b. Edit namelist file for your particular use case.
3c. Edit output streams file.
4. Set up data folder structure.
5. If running historical, set dates.txt
6. Submit driver script.

The first step is to create a functional F compset. Broadly, this is **active atmosphere, active land, active runoff (optional) and data ocean/ice**. Other configurations may work (e.g., active wave, glacier models) but have not been tested. B compsets (fully coupled) should work, but betacast does not initialize the ocean model at this time, so it would rely on the namelist default provided by the modeling system.

### 1. Create case directory and test stability.

This step just requires a built and tested case of CESM. Any source mods or other specific namelist settings should be applied directly to this case. It is only necessary to show the model is stable for a few hours.

**CESM example:**

```
$ cd ~/work/cesm-release/cime/scripts
$ ./create_newcase --case ~/F-betacast-F2000climo --compset F2000climo --res ne30_g16 --mach cheyenne --project UNSB0017 --run-unsupported
$ cd ~/F-betacast-F2000climo
$ ./case.setup
$ ### patch land mods for restart
$ ./case.build
$ ./case.submit
```

**E3SMv2 example:**

```
$ cd ~/E3SM/cime/scripts
$ ./create_newcase --case ~/F-betacast-FC5AV1C --compset F2010C5-CMIP6-LR --res ne30_ne30 --mach cori-knl --project m1637
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
over the top of the file, which injects the correct logic. A similar procedure is used if you want runoff model restart files using `rof_comp_mct.patch`. **This needs to be done before the `./case.build` step.**
2. For E3SM, current suggested compsets are: F2010C5-CMIP6-HR (ne120, VR) and F2010C5-CMIP6-LR (ne30). For SCREAM, these are FSCREAM-LR and FSCREAM-HR, respectively.
3. E3SMv2 (tags from approximately October 2020 onward) are only officially supported. E3SMv1 is effectively supported by choosing modelSystem = 0, although continual updates to support the evolution of EAM and ELM seperately from CAM and CLM/CTSM may eventually break this backwards compatibility.
4. For running with the runoff turned on, it may be necessary to create a new compset with MOSART (or other runoff model like RTM) instead of SROF (stub runoff). An example for E3SMv2 is:

	    <compset>
	      <alias>F2010C5-CMIP6-HR-ROF</alias>
	      <lname>2010_EAM%CMIP6-HR_ELM%SPBC_CICE%PRES_DOCN%DOM_MOSART_SGLC_SWAV</lname>
	    </compset>

5. Empty note.

### 2. Generate an analysis/reanalysis to CAM weight file

Betacast needs a file (ESMF format) that provides high-order weights to take the analysis data and horizontally remap it to the target grid.

This can be done with `${BETACAST}/remapping/gen_GFS_to_SE_weight_file.ncl`. This script requires two inputs that are directly modified in the script body, `srcGridFile` (a SCRIP grid file defining the analysis regular lat-lon grid) and `dstGridFile` (a SCRIP grid file defining the destination CAM grid).

Historically there have been two analysis grid sizes associated with publicly disseminated GFS/CFS/CFSR analyses, 0.5deg (CFSR and GFS pre-2017) and 0.25deg (GFS post-2017). ERA5 data from CDS is on a 0.25deg grid. The SCRIP files for these grids are located in `${BETACAST}/remapping/scrip/`.

### 3. Edit namelists

There is a rudimentary namelist capability betacast uses via Bash. Ideally, this would be in Python or something else but for now this is sufficent.

A bash namelist for this codebase consists of a text file and variables of the format

`VALUE = key`

Each VALUE line is read in, splitting on ` = `.

#### Important namelist notes!

1. There must be at least one ASCII space between and after the splitting `=`!
2. If you want to pass an empty string in, you must define the key as `"___"` (three underscores) since spaces break the splitting.
3. Namelist files may include comments by specifying `#` as the first character of a line.

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
| modelSystem | 0 = CESM + E3SMv1, 1 = E3SMv2 (defaults to 0 if empty or not included) |
| do_runoff | Include runoff model files. 0 = no, 1 = yes (defaults to 0 if empty or not included) |
| atmDataType | What ATM data we want to use? 1 = GFS ANL, 2 = ERA-I, 3 = CFSR, 4 = ERA5 |
| sstDataType | What SST data we want to use? 1 = GDAS, 2 = ERA, 3 = NOAAOI |
| numLevels | 72 -> E3SMv1, 32 -> CAM6, 30 -> CAM5, 26 -> CAM4 |
| numdays | How long for forecast to run (in days) |
| adjust_topo | Full path to a *model* topography file. If valid file/path, code will apply hydrostatic adjustment during atm initial condition step. Turn off by not including variable or setting to empty string. |
| adjust_flags | Hydrostatic adjustment options. Currently "a" (include TBOT adjustment) and "-" (PS adjustment only) are supported. Only applied with valid adjust_topo file. |
| doFilter | Should we apply offline forward DFI? Generally "false" for diffusive dycores and/or SE with hydrostatic adjustment. Set to "true" if using SE with no adjustment to minimize GW noise during first ~72 hours. |
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

In `output_streams`, you can generate a text file that specifies output streams that are to be appended to the model namelist. Some sample options for CESM are included in the repo.

An example of outputting surface pressure, 10-m wind, and precipitation rate every 3 hours, with one time per file (i.e., 8 files per day) is:

```
nhtfrq=-3
mfilt=1
fincl1='PS:I',U10:I','PRECT:I'
```

### 4. Set up data folder structure

Betacast requires a semi-permanent directory structure for which to download analysis data and write forcing data for CESM/E3SM. This does not have to be backed up but may be beneficial to not be truly on scratch space if simulations are to be carried out over a longer period of time.

Directories required are named in the machine namelist file (see 3a below). These can either be created by hand or created by running...
`$ ./tools/setup_data_dirs.sh ../machine_files/machine.MYMACHINEFILE`
from the top-level directory of betacast.

### 5. Dates file

If not running in real-time, dates to be simulated are passed into the script via a text file in the root betacast directory named `dates.{CASENAME}.txt`. Currently only 00Z and 12Z cycles are fully supported in order to minimize disk writes associated with restart files at intermediate cycles.

For example, let's say we are running ne30np4 forecasts with a casename of MY_NE30_BETACAST and want to run a simulation on 00Z August 1, 2nd, and 3rd, 2018. In the `${BETACAST}` directory we would add a text file named `dates.MY_NE30_BETACAST.txt` and add the following lines to the top of the file:

```
2018080100
2018080200
2018080300
2018080312
```

... when `islive` is equal to 0 in the namelist, the code will extract the YYYY, MM, DD, and HH from the 1st line of this text file. NOTE: this line will only be deleted upon successful forecast completion. If the code crashes mid-run, then the interrupted forecast remains at the head of the dates file.

### 5.1. Pre-stage atmospheric analysis data (optional)

Pre-staging the data consists of pulling the atmospheric data from an external server, renaming it, doing some basic file concatenation, variable arrangment, etc. and placing it into the relevant betacast folder. Doing so minimizes issues with attempting to download data from external servers at run-time, which can lead to crashes due to invalid credentials and other issues.

An example of pre-staging data is shown using ERA5. The pre-stage code in this case requires python3 and cdsapi.

```
# Example on Cheyenne of loading python + cdsapi (part of NCAR suite)
module load python
ncar_pylib

# Navigate to ECMWF data acquisition folder
cd ~/betacast/atm_to_cam/getECMWFdata

# Run prestage shell script. Two command line inputs are:
#   $1 = output directory for staged file
#   $2 = requested data in YYYYMMDDHH format.
./prestage-ERA5.sh /glade/work/${LOGNAME}/sewx/ECMWF/ 2011082512
```

### 6. Run model

```$ ./betacast.sh machine_files/machine.cheyenne namelists/nl.conus30x8 output_streams/output.generic```

## SOME NOTES

(WORK IN PROGRESS!)

#### Land model initialization.
It is really a much better idea to run

