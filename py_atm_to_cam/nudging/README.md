# Betacast Nudging File Generator

A tool for generating nudging files for CESM and E3SM climate models from either reanalysis data (ERA5, CFSR) or existing CAM simulations.

## Overview

This tool creates atmospheric nudging files that can be used to constrain CESM/E3SM simulations toward either:

1. Reanalysis datasets (e.g., ERA5 or CFSR)
2. Other CAM/EAM simulations (downscaling or specified dynamics)

The tool handles multiple target grid types (SE/HOMME, MPAS, FV), performs necessary remapping, and supports various temporal frequencies (hourly, 6-hourly, etc.).

## Usage

### Basic Command Line Usage

For interactive or direct execution:

```bash
./gen-nudge.sh ndg.era5.nl
```

### Queue Submission

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
| `dryrun` | If `true`, generate parallel command file without execution |
| `STYR` | Start year (e.g., 2008) |
| `STMON` | Start month (e.g., "jan", "feb") |
| `STDAY` | Start day of month (e.g., 1, 31) |
| `ENYR` | End year (e.g., 2008) |
| `ENMON` | End month (e.g., "jan", "feb") |
| `ENDAY` | End day of month (e.g., 1, 31) |
| `HR_RES` | Time resolution in hours (1 = hourly, 6 = 4xdaily, 24 = daily) |

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
| `BNDTOPO` | Path to target model topography file |
| `WGTNAME` | Path to ESMF weight file for regridding (src->target) |
| `DESCSTR` | Description string for source data |

### Mode Selection

The script operates in one of two modes:

| Option | Description |
|--------|-------------|
| `CAM_TO_CAM` | If `true`, source from existing model sim; if `false`, source from reanalysis |

### CAM-to-CAM Specific Options

The following options are only needed when `CAM_TO_CAM=true`. They are unused when false.

| Option | Description |
|--------|-------------|
| `SUBNAME` | Name identifier of source CAM simulation |
| `BINLIST` | Path/pattern to source data files (supports wildcards) |
| `MODREMAPFILE` | ESMF weight file from source model to intermediate grid |
| `MODINTOPO` | Path to source model topography file |

## Other Features

### Grid Substitution

For ensemble or multi-grid processing, you can use placeholders in paths:

```bash
MODREMAPFILE=/path/to/map_ne0np4natlantic\$GRIDLOWER.ne30x4_TO_era5_0.25x0.25_patc.nc
```

The script will substitute `$GRIDLOWER` with the lowercase version of the grid identifier.

IMPORTANT: The placeholder must be defined earlier in the namelist file since Betacast reads from top to bottom!

### Time Range Processing

The script intelligently processes date ranges and handles various calendar formats:

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