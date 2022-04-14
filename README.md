# Betacast

**Betacast** (originally a portmanteau of "beta" and "forecast") is a software package designed to easily initialize CIME-ized modeling systems in "numerical weather prediction" mode. The code is capable of running in realtime or hindcast mode and possesses tools for generating balanced initial conditions and configuring the model appropriately. **Please cite the below manuscript and this repository if you use Betacast in your research.**

[![DOI](https://zenodo.org/badge/31787694.svg)](https://zenodo.org/badge/latestdoi/31787694)

Reference: C. M. Zarzycki and C. Jablonowski. Experimental tropical cyclone forecasts using a variable-resolution global model. *Monthly Weather Review*, 1430 (10):0 4012--4037, 2015. 10.1175/MWR-D-15-0159.1.

🔴 **IMPORTANT NOTE**: This README assumes some level of familiarity with CESM and/or E3SM. If a user has not used either model before, they are encouraged to view available tutorial materials before proceeding with Betacast.

### Workflow
1. Create a case directory with either a supported or new grid configuration and verify that it is stable/runs given arbitrary inputs.
2. Build an analysis/reanalysis grid to model grid weight file.
3. Set up namelist files.
    1. Edit/create machine file for your particular system.
    2. Create namelist file for your particular use case.
    3. Edit/create output streams file.
4. Set up Betacast data folder structure.
5. If running historical simulations, edit dates.txt
    1. Pre-stage atmospheric data.
6. Decide how to handle land initialization.
7. Run Betacast.

The first step is to create a functional F compset. Broadly, this is **active atmosphere, active land, active runoff (optional) and data ocean/ice**. Other configurations may work (e.g., active wave, glacier models) but have not been tested. B compsets (fully coupled) should work, but betacast does not initialize the ocean model at this time, so it would rely on the namelist default provided by the modeling system.

### 1. Create case directory and test stability.

This step just requires a built and tested case of CESM. Any source mods or other specific namelist settings should be applied directly to this case. It is only necessary to show the model is stable for a few hours.

🔴 **IMPORTANT NOTE**: The resolution defined in this step *must* be the resolution you are using for your experiments since CESM determines atmospheric resolution at build time. For example, if you set up 1deg CESM configuration and then use 0.25deg grids/initial conditions the model will crash! See the `--res` flag and associated CESM documentation for more information about supported resolutions.

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
`$ patch lnd_comp_mct.F90 < ${BETACAST}/patches/lnd_comp_mct.patch`
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

Betacast needs a file (ESMF format) that provides high-order weights to take the analysis data (e.g., ERA5, GFS) and horizontally remap it to the target grid (e.g., CAM, EAM).

This can be done with `${BETACAST}/remapping/gen_analysis_to_model_wgt_file.ncl`. This script requires **four** inputs that are directly modified in the script body.

- `dstGridName` a shortname describing the model grid (for naming purposes only).
- `dstGridFile` a full path to a file defining the destination model grid.
- `anlgrid` is the type of analysis and corresponding grid resolution (*three* are supported, see below).
- `wgtFileDir` is the directory where the weight file should be saved after being generated (this will dictate the path + file used in `gfs2seWeights` in the next section.

`dstGridFile` can be one of *three* formats. It can be a **SCRIP grid file** (contains variables like grid_corner_lat), an **ESMF grid file** (contains variables like nodeCoords), or an **SE/HOMME model output file** (contains dimension ncol). The script will automatically attempt to determine the type of file and create remapping weights accordingly.

Historically there have been two analysis grid sizes associated with publicly disseminated GFS/CFS/CFSR analyses, 0.5deg (CFSR and GFS pre-2017) and 0.25deg (GFS 2017-). ERA5 data from CDS is on a 0.25deg grid. The SCRIP files for these grids are located in `${BETACAST}/remapping/anl_scrip/`.

🔴 **IMPORTANT NOTE**: The CAM weight file needs to be the model grid read during initialization. This is particularly important to note for grids like FV (which has staggered winds) and SE/HOMME (which has dual grids for the dynamics and physics). In the case of SE/HOMME runs, the destination grid is defined by the *physics* grid.

### 3. Edit namelists

There is a rudimentary namelist capability betacast uses via Bash. Ideally, this would be in Python or something else but for now this is sufficent.

A bash namelist for this codebase consists of a text file and variables of the format

`VALUE = key`

Each VALUE line is read in, splitting on ` = `.

#### Important namelist notes!

1. There must be at least one ASCII space between and after the splitting `=`!
2. If you want to pass an empty string in, you must define the key as `"___"` (three underscores) since spaces break the splitting.
3. Namelist files may include comments by specifying `#` as the first character of a line.

### 3.1 Edit machine file for your particular system

In `${BETACAST}/machine_files` there are sample files that define where folders and data files will be stored for your system. There are suggested configurations for Cheyenne and Cori-NERSC, but you may edit these for your workflow or copy/paste for a different system (i.e., university cluster).

| Namelist Variable | Description |
| --- | --- |
| path_to_case | Path to CESM "home" case directory |
| path_to_inputdata | Path to where re/analysis + model initial conditions/forcing data is stored |
| path_to_rundir | Path (top-level) to directory where CESM actively runs |
| sewxscriptsdir | Path to betacast repo (i.e., `${BETACAST}`) |

### 3.2 Edit namelist file for your particular case

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
| adjust_topo | Full path to a *model* (i.e., bnd_topo) topography file. If a valid file/path, code will apply hydrostatic adjustment during atm initial condition step. Turn off by not including variable or setting to empty string. |
| adjust_flags | Hydrostatic adjustment options. Currently "a" (include TBOT adjustment) and "-" (PS adjustment only) are supported. Only applied with valid adjust_topo file. |
| doFilter | Should we apply offline forward DFI? Generally "false" for diffusive dycores and/or SE/HOMME with hydrostatic adjustment. Set to "true" if using SE with no adjustment (or unbalanced IC from another source) to minimize GW noise during first ~72 hours. |
| filterOnly | Exit code after the filter run if doFilter=true (useful for producing ncdata for ensembles) |
| numHoursSEStart | Centerpoint of filter duration (leave at 3), only used if doFilter |
| filterHourLength | Filter duration (leave at 6), only used if doFilter |
| filtTcut | Cut setting for filter (leave at 6), only used if doFilter |
| add_perturbs | Add PGW perturbations for counterfactual runs (leave at false) |
| add_noise | Add white noise to ncdata for ensemble (currently white noise is small, generally leave as false) |
| land_spinup | Cycle land spinup only (unsupported currently, leave false) |
| gfs2seWeights | Full path name of weights file for analysis -> model regridding (see previous section) |
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

### 3.3 Edit output streams

In `output_streams`, you can generate a text file that specifies output streams that are to be appended to the model namelist. Some sample options for CESM are included in the repo.

An example of outputting surface pressure, 10-m wind, and precipitation rate every 3 hours, with one time per file (i.e., 8 files per day) is:

```
nhtfrq=-3
mfilt=1
fincl1='PS:I',U10:I','PRECT:I'
```

### 4. Set up data folder structure

Betacast requires a 'permanent' directory structure for which to download analysis data and write forcing data for the model. This does not have to be backed up but may be beneficial to not be truly on scratch space if simulations are to be carried out over a longer period of time.

Directories required are named in the machine namelist file (see 3.1 above). These can either be created by hand or created by running...
`$ ./tools/setup_data_dirs.sh ../machine_files/machine.MYMACHINEFILE`
from the top-level directory of Betacast. **This only needs to be done once per user per system.**

### 5. Dates file

If not running in real-time, dates to be simulated are passed into the script via a text file in `${BETACAST}/dates/` in a file named `dates.${CASENAME}.txt`. Currently only 00Z and 12Z cycles are fully supported in order to minimize disk writes associated with restart files at intermediate cycles.

For example, let's say we are running ne30np4 forecasts with a casename of MY_NE30_BETACAST and want to run a simulation on 00Z August 1, 2nd, and 3rd, 2018. In the `${BETACAST}/dates/` directory we would add a text file named `dates.MY_NE30_BETACAST.txt` and add the following lines to the top of the file:

```
2018080100
2018080200
2018080300
2018080312
```

... when `islive` is equal to 0 in the namelist, the code will extract the YYYY, MM, DD, and HH from the 1st line of this text file. NOTE: this line will only be deleted upon successful forecast completion. If the code crashes mid-run, then the interrupted forecast remains at the head of the dates file.

**NOTE**: There is backwards compatibility such that the "dates" file can live in `${BETACAST}` root. This is not recommended because the root directory can get quite messy with multiple cases and/or ensembles. If `datestemplate` is specified in the namelist *AND* there is no existing dates file for a particular case, Betacast will copy the template file to `dates.${CASENAME}.txt`. This is useful for running an ensemble of simulations -- you would create one template file of all relevant dates and then each case would copy that file over as a starting point (and each subsuquent model run would use the case-specific file).

Workflow when `islive` is false is:

1. Look for `{BETACAST}/dates/dates.${CASENAME}.txt`.
2. If not **1**, look for `{BETACAST}/dates.${CASENAME}.txt`.
3. If not **2**, look for `datestemplate` in namelist *AND* that `datestemplate` exists. If yes, copy.
4. If not **1**, **2**, or **3**, exit.

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

### 6. Land initialization specification

(WORK IN PROGRESS!)

### 7. Run Betacast

```$ ./betacast.sh machine_files/machine.cheyenne namelists/nl.conus30x8 output_streams/output.generic```

This ends the betacast workflow.

## Generating DATM files

When nudging the land model to initialize CLM/ELM, we need 'forcing' files for DATM. While we can use some existing forcing files in the CESM/E3SM repo, it may be beneficial to initialize using ERA5 forcing. To do so, we need to generate the ERA5 forcing files and data streams for running with an I compset.

The process is pretty straightforward.

1. Download files from ERA5 repository.
2. Gen DATM files using this raw ERA5 data.
3. (optional) add climate deltas for counterfactual runs.
4. Add `user_datm_` files to your `I` case directory.

### 1. Download files

```
cd ${BETACAST}/land-spinup/gen_datm/get-era5
## Need to have cdsapi Python library loaded -- on NCAR see next line
ncar_pylib
## Edit ./driver-get-era5.sh for years + local location for download
nohup ./driver-get-era5.sh &
```

### 2. Generate DATM files

```
cd ${BETACAST}/land-spinup/gen_datm/get-era5
## Set years etc.
qsubcasper driver-gen-datm.sh
```

### 3. (optional) add climate deltas

NOTE: This code will create a duplicate directory of the DATM stream and add perturbations to every file! It would be smart to first set up a control DATM folder with the only forcing files needed (e.g., two years instead of twenty).

```
cd ${BETACAST}/land-spinup/gen_datm/add-perturbs
## Edit driver-perturb-datm.sh
## Edit add_perturbations_to_DATM.ncl
qsubcasper driver-perturb-datm.sh
```

### 4. Add `user_datm` files.

```
## Currently, we cheat and overwrite CRUNCEP with ERA5 until I learn how to create our own stream.
./xmlchange DATM_MODE=CLMCRUNCEPv7
cd ${MYCASEDIR}
cp user_datm .
## edit paths in user_datm files
## in user_nl_datm
## tintalgo = "coszen", "linear", "linear", "linear", "lower"
```

---

### Testing different offsets

For `coszen` it seems like ERA5 likes -10800s for an offset. This number matches the diurnal SWdown solar cycle when comparing a DATM run to an F compset run for an equivalent calendar day.

An example is shown below. This is 12Z averaged over past 6 hours (SWdown:A) on Jan 16th. Left is from an ne30 F compset run, right is from an I compset f09 run. Moving to higher negative offsets shifts the bullseye to the left (west) at a given time. Zero offset would shift the I compset SWdown to the right by 45deg.

/Screen Shot 2021-05-13 at 5.59.57 PM.png

For `linear`, the offset seems best to be 0. This creates a "centering" of the state variables (e.g., T, etc.) on the instantaneous synoptic field produced by ERA5. When doing an ncdiff of the model produced TS and comparing to ERA5, there is very little red/blue shift noted in the diff field compared to other offsets.

Some code:

```
ncks -d time,0,7,2 out.1992.08.nc tmp.nc
ncremap --alg_typ=neareststod --esmf_typ=neareststod -i tmp.nc -o tmp_regrid.nc -d ~/scratch/RoS-ICLM45-f09/run/RoS-ICLM45-f09.clm2.h0.1992-08-01-00000.nc
ncrename -v ssrd,SWdown tmp_regrid.nc
ncrename -v t2m,Tair tmp_regrid.nc
ncrename -v mtpr,RAIN tmp_regrid.nc
ncks -v SWdown,Tair,RAIN tmp_regrid.nc tmp_date.nc
ncrename -d latitude,lat -d longitude,lon tmp_date.nc
ncdiff tmp_date.nc ~/scratch/RoS-ICLM45-f09/run/RoS-ICLM45-f09.clm2.h0.1992-08-01-00000.nc diff.nc
```

---

### Some notes on user_nl_datm from the CLM documentation.

##### offset (in the stream file)

offset is the time offset in seconds to give to each stream of data. Normally it is NOT used because the time-stamps for data is set correctly for each stream of data. Note, the offset may NEED to be adjusted depending on the taxmode described above, or it may need to be adjusted to account for data that is time-stamped at the END of an interval rather than the middle or beginning of interval. The offset can is set in the stream file rather than on the stream namelist. For data with a taxmode method of coszen the time-stamp needs to be for the beginning of the interval, while for other data it should be the midpoint. The offset can be used to adjust the time-stamps to get the data to line up correctly. 
    
##### tintalgo

tintalgo is the time interpolation algorithm. For CLM we usually use one of three modes: coszen, nearest, or linear. We use coszen for solar data, nearest for precipitation data, and linear for everything else. If your data is half-hourly or hourly, nearest will work fine for everything. The coszen scaling is useful for longer periods (three hours or more) to try to get the solar to match the cosine of the solar zenith angle over that longer period of time. If you use linear for longer intervals, the solar will cut out at night-time anyway, and the straight line will be a poor approximation of the cosine of the solar zenith angle of actual solar data. nearest likewise would be bad for longer periods where it would be much higher than the actual values.

- Note: For coszen the time-stamps of the data should correspond to the beginning of the interval the data is measured for. Either make sure the time-stamps on the datafiles is set this way, or use the offset described above to set it. 
- Note: For nearest and linear the time-stamps of the data should correspond to the middle of the interval the data is measured for. Either make sure the time-stamps on the datafiles is set this way, or use the offset described above to set it. 

## Generating 'deltas' for counterfactual simulations

Todo
