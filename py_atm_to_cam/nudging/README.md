# Betacast Nudging File Generator

A tool for generating nudging files for CESM and E3SM climate models from either reanalysis data (ERA5, CFSR) or existing CAM simulations.

## Overview

This tool creates atmospheric nudging files that can be used to constrain CESM/E3SM simulations toward either:

1. Reanalysis datasets (e.g., ERA5 or CFSR)
2. Other CAM/EAM simulations (downscaling or specified dynamics)

The tool handles multiple target grid types (SE/HOMME, MPAS, FV), performs necessary remapping, and supports various temporal frequencies (hourly, 6-hourly, etc.).

## Usage

The two command line arguments are:

1. Nudging file generation namelist.
2. Dates file. If dates file is not to be used, it is safest to pass NULL in.

### Basic Command Line Usage

For interactive or direct execution:

```bash
./gen-nudge.sh ndg.era5.nl NULL
```

### Batch Queue Submission

#### PBS (CASPER)

```bash
qsub -v NLFILE="ndg.era5.nl",input_dates_file="NULL" gen-nudge.sh
```

#### SLURM (NERSC)

```bash
sbatch gen-nudge.sh ndg.era5.nl NULL
```

## Configuration

The configuration is controlled through a bash namelist file that specifies various parameters. An example namelist is provided in the repository, although it is encouraged that users copy then and create their own versus overwriting, which reduces conflicts with the main repo.

### Basic Configuration Options

| Option | Description |
|--------|-------------|
| `STYR` | Start year (e.g., 2008) |
| `STMON` | Start month (e.g., "jan", "feb") |
| `STDAY` | Start day of month (e.g., 1, 31) |
| `ENYR` | End year (e.g., 2008) |
| `ENMON` | End month (e.g., "jan", "feb") |
| `ENDAY` | End day of month (e.g., 1, 31) |
| `HR_RES` | Time resolution in hours (1 = hourly, 6 = 4xdaily, 24 = daily) |
| `NDAYS_PER_DATE` | If specified dates mode (see below), how long is each forecast |

### Path Configuration

| Option | Description |
|--------|-------------|
| `BETACASTDIR` | Path to Betacast source code |
| `OUTDIR` | Top level output directory for nudging files |
| `RDADIR` | Path to RDA repo on NCAR or NERSC machines |

### Model Configuration

| Option | Description |
|--------|-------------|
| `DYCORE` | Target dynamical core (e.g., "se", "fv", "mpas") |
| `GRIDSTR` | Target grid identifier (e.g., "ne30pg3", "f09") |
| `NUMLEVS` | Number of vertical levels in target model |
| `BNDTOPO` | Path to target model topography file (or ncdata file for MPAS) |
| `WGTNAME` | Path to ESMF weight file for regridding (src->target) |
| `DESCSTR` | Description string for source data |

### Optional Options (fun!)

| Option | Description |
|--------|-------------|
| `dryrun` | If `true`, generate parallel command file without execution |
| `CAM_TO_CAM` | If `true`, interpolate from a 3D CAM state |
| `add_scream` | If `true`, perform SCREAM-specific post-processing |

### Date Modes

The code works in two modes.

#### Specified dates mode

If the user passes in a valid, non-empty dates file as the second command line argument to the bash script, the code will first pull all dates from the dates file. It will then use `HR_RES` and `NDAYS_PER_DATE` to create an array of nudgings needed for that forecast. All other date namelist options are ignored.

For example, if the dates file reads:

```bash
2005082800
2005082900
```

and `HR_RES` is 6 and `NDAYS_PER_DATE` is 2, the resulting array of nudging to generate is:

```bash
2005082800
2005082806
2005082812
2005082818
2005082900
2005082906
2005082912
2005082918
2005083000
2005083006
2005083012
2005083018
2005083100
```

This works by going down forecast-by-forecast (each line) and creating an array of required dates to nudge that forecast. The code will internally sort and remove duplicates.

**NOTE:** The special case of "just generate files for the dates in my list" can be achieved by setting `NDAYS_PER_DATE=0`. This may be useful if one just wants to create nudging files for specific user-defined dates only.

#### Auto-generated dates mode

If dates file (second command line arg) is non-existant (or some garbage string like NULL), the code will use the `STYR`, `STMON`, `STDAY`, `ENYR`, `ENMON`, `ENDAY`, and `HR_RES` variables to generate an array of files to be generated.

For example:

```bash
STYR=2021
ENYR=2021
STMON="jun"
STDAY=20
ENMON="jun"
ENDAY=28
HR_RES=6
```

This would process data from June 20-28, 2021 at 6-hour intervals (00Z, 06Z, 12Z, 18Z).

`NDAYS_PER_DATE` is ignored in this case.

### Model Source Selection

The script operates in one of two "source" modes. As noted above, the boolean `CAM_TO_CAM` governs this. If `true`, source from existing model sim; if `false`, source from reanalysis such as ERA5. The default (if not set) is `false`.

### CAM-to-CAM Specific Options

The following options are only needed when `CAM_TO_CAM=true`. They are unused when false.

| Option | Description |
|--------|-------------|
| `SUBNAME` | Name identifier of source CAM simulation |
| `BINLIST` | Path/pattern to source data files (supports wildcards) |
| `MODREMAPFILE` | ESMF weight file from source model to intermediate grid |
| `MODINTOPO` | Path to source model topography file |

### Special model processing

The default behavior is to produce nudging files compatible with CAM5+ and E3SMv0+. SCREAMv1+ use a different file format as of 12/25. Toggling `do_scream=true` will invoke extra steps (NCO required) to massage files into said format. See `SCREAM_CMDS` in the shell script.

## Other Features

### Grid Substitution

For ensemble or multi-grid processing, you can use placeholders in paths:

```bash
MODREMAPFILE=/path/to/map_ne0np4natlantic\$GRIDLOWER.ne30x4_TO_era5_0.25x0.25_patc.nc
```

The script will substitute `$GRIDLOWER` with the lowercase version of the grid identifier.

IMPORTANT: The placeholder must be defined earlier in the namelist file since Betacast reads from top to bottom!

## Output

Nudging files are created with the following naming convention:
```
ndg.[DESCSTR].[GRIDSTR].L[NUMLEVS].cam2.i.YYYY-MM-DD-SSSSS.nc
```

For example:
```
ndg.ERA5.ne30pg3.L58.cam2.i.2021-06-20-00000.nc
```

They are stored in `OUTDIR`.

## Nudging weight files and visualization

### CAM/E3SM

For CAM and E3SM, nudging is controlled by namelist settings (see the other Markdown file in this folder). There is a visualization script that reads the `user_nl_cam` file in the directory and creates a PNG with nudging diagnostics.

```
python Lookat_NudgeWindow.py --NLEV 32
```

### SCREAM

For SCREAM, nudging windows/weights are passed in via a standalone file containing a full 3D (lev, lat, lon) mesh of weights. This can be generated using a Python script that has been pulled from the EAMxx-scripts repository.

```
python SCREAMv1_create_nudging_weights.py -datafile /global/cfs/cdirs/e3sm/inputdata/atm/cam/topo/USGS-gtopo30_ne30np4pg2_x6t-SGH.c20210614.nc -nlev 128 -lat lat -lon lon -weightsfile your_weighting_file.nc
```

The datafile needs to contain the target physics grid (where the nudging is applied). For an npXpgY configuration, the topography file contains both meshes and is a logical choice.