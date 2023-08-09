# Betacast

**Betacast** (originally a portmanteau of "beta" and "forecast") is a software package designed to easily initialize CIME-ized modeling systems in "numerical weather prediction" mode. The code is capable of running in realtime or hindcast mode and possesses tools for generating balanced initial conditions and configuring the model appropriately. **Please cite both the below manuscript and this repository if you use Betacast in your research.**

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
    1. Pre-stage atmospheric data if necessary.
6. Decide how to handle land initialization.
7. Run Betacast.

The first step is to create a functional F compset. Broadly, this is **active atmosphere, active land, active runoff (optional) and data ocean/ice**. Other configurations may work (e.g., active wave, glacier models) but have not been tested. B compsets (fully coupled) should work, but betacast does not initialize the ocean model at this time, so it would rely on the namelist default provided by the modeling system.

### 1. Create case directory and test stability.

This step just requires a built and tested case of CESM. Any source mods or other specific namelist settings should be applied directly to this case. It is only necessary to show the model is stable for a few hours.

🔴 **IMPORTANT NOTE**: The resolution defined in this step *must* be the resolution you are using for your experiments since both CESM and E3SM determine atmospheric resolution at build time. For example, if you set up 1deg CESM configuration and then use 0.25deg grids/initial conditions the model will crash! See the `--res` flag and associated CESM documentation for more information about supported resolutions.

**CESM example:**

<!--
MODELROOT=/glade/u/home/zarzycki/work/cam_20230623/
BETACAST=~/betacast/
PROJECTID=P93300642
CASESDIR=~/tests-betacast/
-->

```
cd ${MODELROOT}/cime/scripts
./create_newcase --case ${CASESDIR}/F-betacast-F2000climo --compset F2000climo --res ne30_g16 --mach cheyenne --project ${PROJECTID} --run-unsupported
cd ${CASESDIR}/F-betacast-F2000climo
./case.setup
${BETACAST}/tools/patch-sfc-mods.sh ${BETACAST} ${MODELROOT} nuopc clm
./case.build
./case.submit
```

NOTE: The above uses the "nuopc" driver. Releases <=CESM2.2 will use "mct" by default.

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

1. In the "patch-sfc-mods" step, a small modification is made to the land model to enforce restart files to be printed every 12 hours. This is done since the land model is initialized via nudging with the atmosphere in this framework. A shell script `${BETACAST}/tools/patch-sfc-mods.sh` has been created to ease this application, which takes the Betacast directory, top-level model directory, driver (mct or nuopc) and model component (e.g., elm, clm, mosart, rtm) as inputs. These patches can be manually applied by copying CESM or E3SM's `lnd_comp_mct.F90` (from the land model source code) into `$CASEDIR/SourceMods/src.clm` (or equivalent) and running
`$ patch lnd_comp_mct.F90 < ${BETACAST}/patches/lnd_comp_mct.patch`
over the top of the file, which injects the correct logic. A similar procedure is used if you want runoff model restart files using `rof_comp_mct.patch`. **This needs to be done before the `./case.build` step.**
2. For E3SM, current suggested compsets are: F2010C5-CMIP6-HR (ne120, VR) and F2010C5-CMIP6-LR (ne30). For SCREAM, these are FSCREAM-LR and FSCREAM-HR, respectively.
3. E3SMv2 (tags from approximately October 2020 onward) are only officially supported. E3SMv1 is effectively supported by choosing modelSystem = 0, although continual updates to support the evolution of EAM and ELM seperately from CAM and CLM/CTSM may eventually break this backwards compatibility.
4. For running with the runoff turned on, it may be necessary to create a new compset with MOSART (or other runoff model like RTM) instead of SROF (stub runoff). An example for E3SMv2 is:

	    <compset>
	      <alias>F2010C5-CMIP6-HR-ROF</alias>
	      <lname>2010_EAM%CMIP6-HR_ELM%SPBC_CICE%PRES_DOCN%DOM_MOSART_SGLC_SWAV</lname>
	    </compset>

5. More notes to come?

### 2. Generate an analysis/reanalysis to CAM weight file

Betacast needs a file (ESMF format) that provides high-order weights to take the analysis data (e.g., ERA5, GFS) and horizontally remap it to the target grid (e.g., CAM, EAM).

This can be done with `${BETACAST}/remapping/gen_analysis_to_model_wgt_file.ncl`. This script requires **four** inputs that are directly modified in the script body.

- `dstGridName` a shortname describing the model grid (for naming purposes only).
- `dstGridFile` a full path to a file defining the destination model grid.
- `anlgrid` is the type of analysis and corresponding grid resolution (*three* are supported, see below).
- `wgtFileDir` is the directory where the weight file should be saved after being generated (this will dictate the path + file used in `anl2mdlWeights` in the next section.

`dstGridFile` can be one of *three* formats. It can be a **SCRIP grid file** (contains variables like grid_corner_lat), an **ESMF grid file** (contains variables like nodeCoords), or an **SE/HOMME model output file** (contains dimension ncol). The script will automatically attempt to determine the type of file and create remapping weights accordingly.

Historically there have been two analysis grid sizes associated with publicly disseminated GFS/CFS/CFSR analyses, 0.5deg (CFSR and GFS pre-2017) and 0.25deg (GFS 2017-). ERA5 data from CDS is on a 0.25deg grid. The SCRIP files for these grids are located in `${BETACAST}/remapping/anl_scrip/`.

🔴 **IMPORTANT NOTE**: The CAM weight file needs to be the model grid read during initialization. This is particularly important to note for grids like FV (which has staggered winds) and SE/HOMME (which has dual grids for the dynamics and physics). In the case of SE/HOMME runs, the destination grid is defined by the *physics* grid.

### 3. Edit namelists

There is a rudimentary namelist capability betacast uses via Bash. Ideally, this would be in Python or something else but for now this is sufficent.

A bash namelist for this codebase consists of a text file and variables of the format

`VALUE = key`

Each VALUE line is read in, splitting on ` = `.

#### Important namelist notes!

1. There must be at least one ASCII space between and after the splitting `=`! (As of 10/2022 this is no longer a requirement and `VALUE=key` is acceptable, although best practice remains to use `VALUE = key`)
2. If you want to pass an empty string in, you must define the key as `"___"` (three underscores) since spaces break the splitting.
3. Namelist files may include comments by specifying `#` as the first character of a line.

### 3.1 Edit machine file for your particular system

In `${BETACAST}/machine_files` there are sample files that define where folders and data files will be stored for your system. There are suggested configurations for Cheyenne and Cori-NERSC, but you may edit these for your workflow or copy/paste for a different system (i.e., university cluster).

| Namelist Variable | Description |
| --- | --- |
| path\_to\_case | Path to CESM "home" case directory |
| path\_to\_inputdata | Path to where re/analysis + model initial conditions/forcing data is stored |
| path\_to\_rundir | Path (top-level) to directory where CESM actively runs |
| sewxscriptsdir | Path to betacast repo (i.e., `${BETACAST}`) |

### 3.2 Edit namelist file for your particular case

In `${BETACAST}/namelist_files` there are sample files that define the forecast configuration. This is the primary location where run settings are specified. Betacast uses a general philosophy that 0 = false/no and 1 = true/yes.

| Namelist Variable | Description |
| --- | --- |
| ARCHIVEDIR | Top-level directory for archiving runs (ARCHIVEDIR/CASE/YYYYMMDDHH). Do not set for default archive in rundir. |
| debug | Setting to true (1) adds debugging options. Otherwise leave at false (0) |
| islive | if true (1) then pull GDAS/GFS from server in real-time, false (0) is "hindcast" mode |
| datestemplate | If islive = false, dates.XXX.txt file to copy if Betacast cannot find an existing dates file |
| runmodel | Unused, set to "true" |
| archive_inic | Add (NCO-compressed) initial conditions for component models to archive directory (0 = no (default), 1 = yes) |
| compress_history_nc | Use NCO lossless compression to compress history files (0 = no (default), 1 = yes) |
| tararchivedir | Should the archive folder be tarred?  (0 = no, 1 = yes (default)) |
| modelSystem | 0 = CESM + E3SMv1, 1 = E3SMv2+/SCREAM (defaults to 0 if empty or not included) |
| cime_coupler | Which driver to use? Can be "mct" or "nuopc". Default is "mct" if not specified. |
| do_runoff | Include runoff model files (false/true) (defaults to false if empty or not included) |
| atmDataType | What ATM data we want to use? 1 = GFS ANL, 2 = ERA-I, 3 = CFSR, 4 = ERA5 |
| sstDataType | What SST data we want to use? 1 = GDAS, 2 = ERA, 3 = NOAAOI |
| numLevels | 128 -> SCREAM, 72 -> E3SM, 58 -> CAM7, 32 -> CAM6, 30 -> CAM5, 26 -> CAM4 |
| numdays | How long for forecast to run (in days) |
| adjust_topo | Full path to a *model* (i.e., bnd_topo) topography file. If a valid file/path, code will apply hydrostatic adjustment during atm initial condition step. Turn off by not including variable or setting to empty string. |
| adjust_flags | Hydrostatic adjustment options. Currently "a" (include TBOT adjustment) and "-" (PS adjustment only) are supported. Only applied with valid adjust_topo file. |
| doFilter | Should we apply offline forward DFI? Generally "false" for diffusive dycores and/or SE/HOMME with hydrostatic adjustment. Set to "true" if using SE with no adjustment (or unbalanced IC from another source) to minimize GW noise during first ~72 hours. |
| filterOnly | Exit code after the filter run if doFilter=true (useful for producing ncdata for ensembles) |
| numHoursSEStart | Centerpoint of filter duration (leave at 3), only used if doFilter |
| filterHourLength | Filter duration (leave at 6), only used if doFilter |
| filtTcut | Cut setting for filter (leave at 6), only used if doFilter |
| add_perturbs | Add PGW perturbations for counterfactual runs? Leave at false generally. |
| perturb_namelist | Path to "perturbation" namelist for counterfactual climate simulations |
| add_noise | Add white noise to ncdata for ensemble (currently white noise is small, generally leave as false) |
| land_spinup | Cycle land spinup only (unsupported currently, leave false) |
| keep_land_restarts | 0 = delete land/rof restart files, 1 = archive land/rof restart files (possibly overwriting those in ${CASE}/run/landstart) |
| override_rest_check | If true, overrides internal check for SourceMods for lnd/rof restarts (default: false) |
| save_nudging_files | false (default) doesn't output initial condition files, true outputs and archives inithist files for use in future nudging runs |
| landrawdir | For CLM5, path to CLM restart files to check/interpolate from if native grid finidat does not exist |
| predict_docn | 0 = persist t=0 SST/ice fields for duration of simulation, 1 = superimpose initialization anomalies on time-varying climatology |
| anl2mdlWeights | Full path name of weights file for analysis -> model regridding (see previous section) |
| PROJECTID | Project ID for run submissions |
| FILTERWALLCLOCK | Wall clock time for filter run |
| FILTERQUEUE | Submission queue for filter run |
| RUNWALLCLOCK | Wall clock time for forecast run |
| RUNQUEUE | Submission queue for forecast run |
| usingCIME | Are we using CIME (set to "true" unless using a very old CESM tag or unsupported GCM) |
| DTIME | Physics timestep (in seconds) |
| FINERES | Finest resolution of SE grid |
| USERSTAB | Required dynamics timestep (in s), negative values try internal calculation, but use with caution |
| use_nsplit | If true, use the CESM/E3SMv1 SE/HOMME nsplit timestep logic, if false apply new E3SMv2 se_tstep parameter (equal to USERSTAB) |
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

[Land initialization link](#lnd_initial_conditions)

### 7. Run Betacast

```$ ./betacast.sh machine_files/machine.cheyenne namelists/nl.ne30.chey output_streams/output.generic```

This ends the betacast workflow.

---

## atm\_to\_cam

This tool takes a 3-D analysis field from supported reanalysis products (CFSR, ERA5, ERAI), NWP analyses (GFS, RAP), or CAM/EAM full states and maps the data correctly onto a target grid for use in initialization (or nudging) of CAM or EAM.

### General workflow

1. Generate a `wgt_filename` (see below) to provide information regarding how to map the analysis grid to the target (model) grid.
2. (Optional) If using an unsupported vertical grid, add new level template file.
3. (Optional) Pre-stage data if not on Cheyenne.
4. Set up command line options for atm_to_cam.ncl
5. Using the below examples as a template, run `atm_to_cam.ncl`.

### Command line options for atm_to_cam.ncl

| Namelist Variable | Type | Description | Default | Required? |
| --- | --- | --- | --- | --- |
| datasource | Str | CFSR, ERA5, ERAI, GFS, CAM, RAP | | Y |
| numlevels | Int | Number of vertical levels | | Y |
| YYYYMMDDHH | Int | Initialization date in YYYYMMDDHH format | | Y |
| dycore | Str | fv, se, mpas | se | N |
| RDADIR | Str | Base path to RDA folder | "" | N |
| data_filename | Str | Full path to file containing initial information | | Y |
| wgt_filename | Str | Full path to ESMF weight file from ANL -> MOD | | Y |
| mpas_as_cam | Bool | If true, write MPAS in CAM physics output (ncol, for nudging), if false write MPAS in MPAS-A format (nCells, for init) | false | N |
| compress_file | Bool | If true, will attempt NetCDF "chunking" compression within NCL | false | N |
| write_floats | Bool | If true, write outputs as single instead of double precision | false | N |
| add_cloud_vars | Bool | If true, add CLDICE and CLDLIQ to output file | true | N |
| adjust_config | String | String defining how to perform hydro adjustment: if string is not empty, will do config. If "a" also try and correct TBOT in addition to PS | "" | N |
| model_topo_file | String | If MPAS, an MPAS inic file, otherwise a file containing PHIS for FV or SE | "" | N |
| se_inic | String | Full path of file to write as output | | Y |
| mod_in_topo | Str | Full path to PHIS field from *downscaling* MOD for pressure surface calculation purposes | "" | N |
| mod_remap_file | Str | Full path to ESMF weight file that goes *downscaling* MOD -> ANL | "" | N |

### Examples

#### Regridding ERA5 Cheyenne RDA --> SE ne30

```
ncl -n atm_to_cam.ncl 'datasource="ERA5RDA"' \
  numlevels=32 \
  YYYYMMDDHH=2019120100 \
  'dycore="se"' \
  'data_filename = "/glade/collections/rda/data//ds633.0/e5.oper.invariant/197901/e5.oper.invariant.128_129_z.ll025sc.1979010100_1979010100.nc"' \
  'wgt_filename="/glade/u/home/$LOGNAME/betacast/remapping/map_gfs_0.25x0.25_TO_ne30_patc.nc"' \
  'RDADIR="/glade/collections/rda/data/ds633.0/"' \
  'model_topo_file="/glade/p/cesmdata/cseg/inputdata/atm/cam/topo/se/ne30np4_nc3000_Co060_Fi001_PF_nullRR_Nsw042_20171020.nc"' \
  'adjust_config=""' \
  compress_file=False \
  write_floats=True \
  add_cloud_vars=True \
  'se_inic = "/glade/scratch/$LOGNAME/my_se_initial_condition_file.nc"'
```

### Specific notes and settings:

#### Generating a weight file

`wgt_filename` contains an ESMF file that provides mapping weights to go from the analysis grid to the model grid. An example would be to go from ERA5 0.25x0.25deg to CAM-SE ne30np4.

This file can be generated by using the NCL regridding scripts inside `${BETACAST}/remapping/gen_analysis_to_model_wgt_file.ncl`. Here, you must have a SCRIP grid descriptor for your target mesh. A series of included grids are found in `${BETACAST}/remapping/anl_scrip` for common analyses.

#### Downscaling CAM data

For downscaling CAM data one must pass in:

- mod_in_topo
- mod_remap_file

`mod_remap_file` should be an ESMF file going from the CAM input data to the analysis grid of interest (unless the forcing data is high-res, using the ERA5 0.25deg grid seems to be safest -- e.g., ne120 -> 0.25x0.25). `mod_in_topo` should be the path to the bnd_topo file containing PHIS for the downscaling model grid (e.g., ne120).

`data_filename` must contain a snapshot of 3-D U, V, T, Q, and 2-D PS. The requested `YYYYMMDDHH` must lie within the bounds of the time dimension in `data_filename`.

#### MPAS data

For MPAS data, the input variable `model_topo_file` actually contains the full initial condition file generated by MPAS's `init_atmosphere` routine. If one sets `mpas_as_cam` to False, `model_topo_file` will be copied to `se_inic` and then the state fields overwritten by the script so that the `se_inic` file can be passed in as an `ncdata` file.

#### Hydrostatic adjustment

Currently the hydrostatic adjustment is only available when `dycore` is set to SE. The code is activated if a valid `model_topo_file` containing `PHIS` is passed in and `adjust_config` is not an empty string. If `adjust_config=a` that applies a special case where both PS and TBOT are adjusted. Otherwise, only PS is adjusted. All state fields are then reinterpolated based on the new hybrid pressure levels associated with the updated PS.

#### Using RDA data

On Cheyenne (and when mirrored on other systems) Betacast can go directly to RDA to get ERA5 data. To do this, point `RDADIR` to the datasets top-level directory (e.g., /glade/collections/rda/data/ds633.0/). In this case, this file will be used if the hydrostatic adjustment code is active where `data_filename` is the invariant field containing the surface Z field (e.g., `${RDADIR}/e5.oper.invariant/197901/e5.oper.invariant.128_129_z.ll025sc.1979010100_1979010100.nc)`.

--

<a name="lnd_initial_conditions"></a>
## Generating a CLM/ELM initial condition

Betacast uses analysis/observations to initialize the atmosphere, data ocean, and data ice models, but is not able to interpolate land surface initial conditions at this time. The preferred method of initializing the land surface for an initial betacast run is to 'spin up' the land surface model by forcing it with observed atmospheric fluxes, which allows the land surface to asymptote to a state in balance with observed atmospheric forcing. In CESM and E3SM, this is known as an I compset.

For Betacasts run in succession, Betacast will use a short forecast from a previous betacast (e.g., the 12Z Jan 21 betacast would use the +12 hr land forecast from the 00Z Jan 21 betacast) which is stored in the `landstart` subdirectory within `${CASE}/run`. However, for a brand new Betacast (or set of Betacasts) an initial land file must be generated to be passed into CLM/ELM as `finidat`. To do this, Betacast runs a data atmosphere model that forces the land model.

The main process is:

```
cd $BETACAST/land-spinup
./auto-script.sh $MODELSYSTEM $DOERA5 $DATEYYYYMMDD $NMONTHS $NCYCLES $ANOMYEAR $REFYEAR $MACHNAMELIST
```

where `auto-script` takes in eight command line inputs:

| Namelist Variable | Description |
| --- | --- |
| MODELSYSTEM | Modeling system to use (integer, 0 = CLM, 1 = ELM) |
| DOERA5 | Use [ERA5](#era5_data_files) to override CRU/NCEP? (integer, 0 = True, 1 = False) |
| DATEYYYYMMDD | Date a CLM/ELM restart file is needed in YYYYMMDD (00Z) |
| NMONTHS | Integer number of months to spinup (1-12) |
| NCYCLES | Integer number of cycles to spinup (>=1) |
| ANOMYEAR | Integer anomaly year (which year of anomalies to apply) |
| REFYEAR | Integer reference year to correct anomalies (if negative, use raw anomalies) |
| MACHNAMELIST | Namelist file including machine or experiment specific settings |

An example would be:

```
./auto-script.sh 1 0 19960113 12 3 2080 1996 nl.landspinup.cori
```

This would spinup ELM/E3SM (1) using ERA5 DATM (0) for a Betacast initialization on Jan 13 1996 (19960113) using 12 months of spinup and 3 cycles. The deltas are taken from the year 2080 and corrected relative to the year 1996 and the namelist settings are specified in `nl.landspinup.cori`.

The simplest land spinup, which should work natively on supported machines with standard CESM/E3SM data repositories (i.e., no ERA5, no anomalies) would be.

```
./auto-script.sh 1 1 19960113 12 3 -1 -1 nl.landspinup.cori
```

The namelist (ex: `nl.landspinup.cori`) includes the following settings, which are mainly model/system specific.

| Namelist Variable | Description |
| --- | --- |
| CIMEROOT | Path to root of E3SM/CESM |
| PATHTOCASE | Path where you want case generated |
| ICASENAME | Name for spinup case (this is a base name, the forecast date, nmonths, and anomaly year (if >0) are appended) |
| PROJECT | Project charging ID |
| MACHINE | Machine (CIME) |
| NNODES | Number of nodes you want to run on |
| RESOL | What resolution? Make sure land matches Betacast land resolution |
| RUNQUEUE | Queue to run in |
| WALLCLOCK | Queue to run in |
| NMONTHSSPIN | Integer number of months to spinup (1-12) |
| BETACAST | Absolute path to Betacast (only used if doERA5 = 0) |
| BETACAST\_DATM\_FORCING\_BASE | Path to DATM_FORCING files (only used if doERA5 = 0) |

So the general workflow is:

1. Ensure CLM/ELM is able to run (ex: do you have an `fsurdat` for the grid?).
2. Edit `$MACHNAMELIST` as needed.
3. (optional) Ensure [ERA5 DATM files](#era5_data_files) and/or anomalie files are available on your machine and the path is specified in `$MACHNAMELIST`.
4. Run `./auto-script.sh` specifying the above command line options. Wait for model to configure, build, and submit.
5. Pending successful completion of `auto-script`, stage initial files in `$CASENAME/landstart` for Betacast.

Invoking `./auto-script.sh` should configure, build, and submit the I compset. The model will then run based on the scheduler. Following the (hopefully successful) run, the spunup land restart (i.e., initial condition) file will be located in the I compset directory and will need to be copied to the `landstart` subfolder in the coupled Betacast directory to be use pulled for initialization. A simple bash snippet is below, although one could also copy these files to a common directory and symlink them, etc.

```
CASENAME=RoS-F2010C5-ne0conus30x8-001-PI
ICASENAME=RoS-ICLM45-ne0conus30x8-ERA5

# Define
LANDFILEDIR=/global/homes/c/czarzyck/scratch/e3sm_scratch/cori-knl/${CASENAME}/run/landstart/
ICASEDIR=/global/homes/c/czarzyck/scratch/e3sm_scratch/cori-knl/${ICASENAME}/run/

# Make land directory inside of Betacast (CASENAME) dir
mkdir -p ${LANDFILEDIR}

# cd to completed I compset directory and copy files
cd ${ICASEDIR}
cp -v ${ICASENAME}.elm.r.*.nc ${LANDFILEDIR}
cp -v ${ICASENAME}.mosart.r.*.nc ${LANDFILEDIR}

# cd to Betacast land directory and rename to match Betacast CASENAME
cd ${LANDFILEDIR}
rename ${ICASENAME} ${CASENAME} *.nc
```


## Running with the climate change 'deltas' code

Thermodynamic deltas (e.g., warming signal, climate change fingerprints, pseudo-global warming forcing) that are derived from long-term climate integrations can be applied over the top of a Betacast analysis to simulate how a particular event would conditionally evolve under a different climate. Mechanistically, Betacast treats this as if it was just a regular forecast, except after generating the initial conditions two seperate scripts are run to overlay these fingerprints over the top of the atmospheric and ocean initial conditions in a consistent fashion.

Some relevant reading:

- M. F. Wehner, C. M. Zarzycki, and C. Patricola. Estimating the human influence on tropical cyclone intensity as the climate changes. In *Hurricane Risk*, pp. 235-260, Springer Books, 2019. 10.1007/978-3-030-02402-4_12
- K. A. Reed, A. M. Stansfield, M. F. Wehner, and C. M. Zarzycki. Forecasted attribution of the human influence on Hurricane Florence. *Science Advances*, 6(1), eaaw9253, 2020. 10.1126/sciadv.aaw9253.

🔴 **IMPORTANT NOTE**: This code generates deltas for the atmosphere and ocean/ice boundary. It does *not* generate initial conditions for the land surface. These need to be generated with the offline spinup process. See instructions below.

To activate this code, the following bool must be selected in the primary namelist:

```
add_perturbs = true
```

and a seperate namelist containing information about how to apply the deltas must be generated and the full path to this secondary file be specified in the primary namelist:

```
perturb_namelist = /global/homes/c/czarzyck/betacast/namelists/perturb.sample.nl
```

`perturb_namelist` is a text file that is read by a routine that runs after atm_to_cam and sst_to_cam which contains the following variables (all are required, even if **False** or empty strings).

| Namelist Variable | Type | Description |
| --- | --- | --- |
| case | string | "Deltas" case (currently only CESMLENS supported) |
| basedir | string | Base directory where `case` deltas are stored |
| start_month | int | Path to where re/analysis + model initial conditions/forcing data is stored |
| end_month | int | Path (top-level) to directory where CESM actively runs |
| current_year | int | Reference year to calculate deltas from |
| comp_year | int | Target year of deltas |
| correct_sfc | bool | Shift entire T/Q vertical profile based on surface delta? (generally **False**) |
| plevs | bool | Are we on pressure levels? (if **False**, means hybrid model levels) |
| update_pressure | bool | Update PS based on PS deltas (generally **False**) |
| update_winds | bool | Update winds based on some wind deltas or thermal wind (generally **False**) |
| do_ps_corr | bool | Apply simple linear PS correction to attempt to minimize geostrophic shock? (generally **True**) |
| esmf_remap | bool | Use ESMF remapping instead of NCL internal (generally **True**) |
| keep_esmf | bool | Keep the ESMF files between calls to add_perturbations to CAM? (generally **True**) |
| smooth_deltas | bool | Smooth deltas using a 9-point smoother? (generally **False** unless big resolution mismatch) |
| smooth_delta_iter | int | If *smooth_deltas=True*, how many iterations to apply? (higher numbers mean more smoothing) |
| output_atm_diag | bool | Output a separate diagnostics file with deltas and other info? |
| extra_diags_atm | bool | Output additional diags? (currently just precipitable water) |
| adjust_ice | bool | Adjust ice fraction based on deltas + freezing/melting? (generally **True**) |
| output_sst_diag | bool | Output a separate diagnostics file with deltas and other info? |

A template namelist is in the repo under `$BETACAST/namelists/perturb.sample.nl`. It is strongly suggested to just copy this file and edit individual keys as desired to ensure all required variables are present on the file. Betacast should print information regarding how the deltas (and their successful application) which should be verified by the user.

---

## Miscellaneous notes regarding DATM

<a name="era5_data_files"></a>
### Generating ERA5 DATM files

When nudging the land model to initialize CLM/ELM, we need 'forcing' files for DATM. While we can use some existing forcing files in the CESM/E3SM repo, it may be beneficial to initialize using ERA5 forcing. To do so, we need to generate the ERA5 forcing files and data streams for running with an I compset.

The process is pretty straightforward.

1. Download files from ERA5 repository.
2. Gen DATM files using this raw ERA5 data.
3. (optional) add climate deltas for counterfactual runs -- although easier to do with `anomaly' streams instead.
4. Add `user_datm_` files to your `I` case directory.

#### 1. Download files

```
cd ${BETACAST}/land-spinup/gen_datm/get-era5
## Need to have cdsapi Python library loaded -- on NCAR see next line
ncar_pylib
## Edit ./driver-get-era5.sh for years + local location for download
nohup ./driver-get-era5.sh &
```

#### 2. Generate DATM files

```
cd ${BETACAST}/land-spinup/gen_datm/get-era5
## Set years etc.
qsubcasper driver-gen-datm.sh
```

#### 3. (optional) add climate deltas

🔴 **IMPORTANT NOTE**: This has been superceded by the ability to maintain a single deck of ERA5 files and then create seperate `anomaly' streams that are overlaid on ERA5. This code will remain (and is scientifically valid) but requires more disk space than the other method.

NOTE: This code will create a duplicate directory of the DATM stream and add perturbations to every file! It would be smart to first set up a control DATM folder with the only forcing files needed (e.g., two years instead of twenty).

```
cd ${BETACAST}/land-spinup/gen_datm/add-perturbs
## Edit driver-perturb-datm.sh
## Edit add_perturbations_to_DATM.ncl
qsubcasper driver-perturb-datm.sh
```

#### 4. Add `user_datm` files.

NOTE: See some commentary on [testing different offsets](#testing_different_offsets)

```
## Currently, we cheat and overwrite CRUNCEP with ERA5 until I learn how to create our own stream.
./xmlchange DATM_MODE=CLMCRUNCEPv7
cd ${MYCASEDIR}
cp user_datm .
## edit paths in user_datm files
## in user_nl_datm
## tintalgo = "coszen", "linear", "linear", "linear", "lower"
```


<a name="testing_different_offsets"></a>
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

### Some notes on user_nl_datm from the CLM documentation.

##### offset (in the stream file)

offset is the time offset in seconds to give to each stream of data. Normally it is NOT used because the time-stamps for data is set correctly for each stream of data. Note, the offset may NEED to be adjusted depending on the taxmode described above, or it may need to be adjusted to account for data that is time-stamped at the END of an interval rather than the middle or beginning of interval. The offset can is set in the stream file rather than on the stream namelist. For data with a taxmode method of coszen the time-stamp needs to be for the beginning of the interval, while for other data it should be the midpoint. The offset can be used to adjust the time-stamps to get the data to line up correctly.

##### tintalgo

tintalgo is the time interpolation algorithm. For CLM we usually use one of three modes: coszen, nearest, or linear. We use coszen for solar data, nearest for precipitation data, and linear for everything else. If your data is half-hourly or hourly, nearest will work fine for everything. The coszen scaling is useful for longer periods (three hours or more) to try to get the solar to match the cosine of the solar zenith angle over that longer period of time. If you use linear for longer intervals, the solar will cut out at night-time anyway, and the straight line will be a poor approximation of the cosine of the solar zenith angle of actual solar data. nearest likewise would be bad for longer periods where it would be much higher than the actual values.

- Note: For coszen the time-stamps of the data should correspond to the beginning of the interval the data is measured for. Either make sure the time-stamps on the datafiles is set this way, or use the offset described above to set it.
- Note: For nearest and linear the time-stamps of the data should correspond to the middle of the interval the data is measured for. Either make sure the time-stamps on the datafiles is set this way, or use the offset described above to set it.
