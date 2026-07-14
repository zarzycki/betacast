# Betacast Nudging File Generator

This tool leverages the atm_to_cam.py script to create atmospheric nudging files that can be used to constrain CESM/E3SM simulations toward either:

1. Reanalysis datasets (e.g., ERA5, CFSR, CR20V3, etc.)
2. Other CAM/EAM simulations (downscaling or specified dynamics)

The tool handles multiple target grid types (SE/HOMME, MPAS, FV), performs necessary remapping, and supports various temporal frequencies (hourly, 6-hourly, etc.).

## Prerequisites

Before running, you need:

1. **A namelist file** — (i.e., "the namelist") Copy one of the example namelists (e.g., `ndg.era5.nl`) and edit it for your case. Do not edit the examples directly to avoid conflicts with the main repo.

2. **A submit wrapper with your paths** — (i.e., "the submit wrapper") Copy the appropriate wrapper from `submit/` and update `BETACAST`, `OUTDIR`, `RDADIR`, and `NUMCORES` for your environment.

3. **An ESMF weight file (`WGTNAME`)** — This maps the source analysis grid to your target model grid. Generate one using `py_remapping/gen_analysis_to_model_wgt_file.py` if you don't already have one. Example for ERA5 → ne30np4:

   ```bash
   python py_remapping/gen_analysis_to_model_wgt_file.py \
     --ANLGRID era5_0.25x0.25 \
     --DSTGRIDNAME ne30np4 \
     --DSTGRIDFILE /path/to/grids/model_scrip/ne30np4_091226_pentagons.nc \
     --WGTFILEDIR /path/to/output/weights/
   ```

   SCRIP files for common analysis grids are in `grids/anl_scrip/`. SCRIP files for common model grids are in `grids/model_scrip/`. Valid `--ANLGRID` options include: `era5_0.25x0.25`, `era5_0.3gaus`, `era5_2deg`, `gfs_0.25x0.25`, `gfs_0.50x0.50`.

4. **A topography file (`BNDTOPO`)** — The target model topography (or MPAS `ncdata` file).

## Usage

The command line arguments are:

1. Nudging file generation namelist.
2. (Optional) Dates file. If dates file is not to be used, it is safest to pass NULL in.
3. (Optional) INDEX for grid substitution in ensemble workflows.

### Batch Queue Submission (Recommended)

In general, you probably want to submit to a batch queue to ensure memory usage isn't an issue and to leverage parallelism via GNU parallel.

Machine-specific submit wrappers live in the `submit/` directory. Each wrapper sets the scheduler preamble, loads conda, exports machine paths (`BETACAST`, `OUTDIR`, `RDADIR`, `NUMCORES`), and then calls `gen-nudge.sh`.

**Edit the wrapper for your machine before submitting**... at minimum update the path exports and account number.

#### PBS (Casper, Derecho)

PBS requires `-v` to pass arguments into the job:

```bash
qsub -v NLFILE="ndg.era5.nl" submit/casper.sh
qsub -v NLFILE="ndg.era5.nl",input_dates_file="dates.txt" submit/casper.sh
qsub -v NLFILE="ndg.era5.nl",input_dates_file="dates.txt",INDEX="001" submit/casper.sh
```

#### SLURM (Perlmutter)

SLURM passes trailing arguments directly:

```bash
sbatch submit/pm-cpu.sh ndg.era5.nl
sbatch submit/pm-cpu.sh ndg.era5.nl dates.txt
sbatch submit/pm-cpu.sh ndg.era5.nl dates.txt 001
```

### Interactive / Direct Execution

You can also run `gen-nudge.sh` directly. Just export the required environment variables first:

```bash
export BETACAST=/path/to/betacast
export OUTDIR=/scratch/ndg
export RDADIR=/path/to/era5/rda
export NUMCORES=12
./gen-nudge.sh ndg.era5.nl NULL
```

### Adding a New Machine

To add support for a new machine, create a new file in `submit/` following the pattern of the existing wrappers (e.g., `casper.sh` for PBS, `pm-cpu.sh` for SLURM). Set the appropriate scheduler headers, module loads, and path exports for your machine. If I've done things correctly, `gen-nudge.sh` is set up agnostic.

## Configuration

The configuration is controlled through the namelist that specifies various parameters. An example namelist is provided in the repository, although it is encouraged that users copy it and create their own versus overwriting. This will reduce conflicts with the main repo in the event a PR is needed.

### Path Configuration (in "the submit wrapper")

These are set by the submit wrapper (see `submit/` directory). If running interactively, export them as environment variables.

| Environment Variable | Description |
|--------|-------------|
| `BETACAST` | Path to Betacast source code (falls back to `BETACASTDIR` in namelist) |
| `OUTDIR` | Top level output directory for nudging files |
| `RDADIR` | Path to RDA data on NCAR or NERSC machines; defaults to ERA5 RDA path. If set in the namelist, the namelist value takes priority over the environment variable (useful for non-ERA5 datasets like CR20V3). |
| `NUMCORES` | Number of GNU parallel cores (default: 12) |

### Namelist Configuration Options

| Option | Description |
|--------|-------------|
| `STYR` | Start year (e.g., 2008) |
| `STMON` | Start month (e.g., "jan", "feb", "mar") |
| `STDAY` | Start day of month (e.g., 1, 31) |
| `ENYR` | End year (e.g., 2008) |
| `ENMON` | End month (e.g., "jan", "feb", "mar") |
| `ENDAY` | End day of month (e.g., 1, 31) |
| `HR_RES` | Time resolution in hours (1 = hourly, 6 = 4xdaily, 24 = daily) |
| `NDAYS_PER_DATE` | Only used in specified dates mode (see Date Modes below): number of forecast days per date |

NOTE: All dates are inclusive.

### Model Configuration (in "the namelist")

| Option | Description |
|--------|-------------|
| `DYCORE` | Target dynamical core (e.g., "se", "fv", "mpas") |
| `GRIDSTR` | Target grid identifier (e.g., "ne30pg3", "f09") |
| `NUMLEVS` | Number of vertical levels in target model |
| `BNDTOPO` | Path to target model topography file (or ncdata file for MPAS) |
| `WGTNAME` | Path to ESMF weight file for regridding (src->target); generate with `py_remapping/gen_analysis_to_model_wgt_file.py` |
| `DESCSTR` | Description string for source data |

### Optional Options (fun! in "the namelist")

| Option | Description |
|--------|-------------|
| `dryrun` | If `true`, generate parallel command file without execution |
| `SOURCE` | Data source: `ERA5` (default), `CR20V3`, or `CAM` |
| `add_scream` | If `true`, perform SCREAM-specific post-processing |

## Other stuff...

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

The `SOURCE` namelist variable controls where input data comes from. Valid options:

| Value | Description |
|-------|-------------|
| `ERA5` | ERA5 reanalysis from the NCAR/NERSC RDA (default) |
| `CR20V3` | 20th Century Reanalysis v3 from RDA; requires `RDADIR` to be set in the namelist |
| `CAM` | Existing CAM/EAM simulation (see CAM-to-CAM options below) |

If `SOURCE` is not set, it defaults to `ERA5`. The old `CAM_TO_CAM=true` boolean is still accepted for backward compatibility but is deprecated in favor of `SOURCE=CAM`.

### CAM-to-CAM Specific Options (in "the namelist")

The following options are only needed when `SOURCE=CAM`. They are unused otherwise.

| Option | Description |
|--------|-------------|
| `SUBNAME` | Name identifier of source CAM simulation |
| `BINLIST` | Path/pattern to source data files (supports wildcards) |
| `MODREMAPFILE` | ESMF weight file from source model to intermediate grid |
| `MODINTOPO` | Path to source model topography file |

### Special model processing

The default behavior is to produce nudging files compatible with CAM5+ and E3SMv0+. SCREAMv1+ use a different file format as of December 2025. Toggling `add_scream=true` will invoke extra steps (NCO required!) to massage files into said format. See `SCREAM_CMDS` in the shell script.

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

```
INDEX=004
python py_remapping/gen_analysis_to_model_wgt_file.py --ANLGRID era5_0.25x0.25 --DSTGRIDNAME TClandfall-${INDEX}_ne192x8_np4_scrip --DSTGRIDFILE /global/cfs/cdirs/m2637/E3SM_SCREAM_files/grids/scrip/TClandfall-${INDEX}_ne192x8_np4_scrip.nc --ANLGRIDPATH ../grids/anl_scrip/ --WGTFILEDIR /global/cfs/cdirs/m2637/betacast/sewx/mapping/
```

The datafile needs to contain the target physics grid (where the nudging is applied). For an npXpgY configuration, the topography file contains both meshes and is a logical choice.
