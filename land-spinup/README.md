# Land Spinup for Betacast

This directory contains tools to generate CLM/ELM land surface initial conditions for Betacast runs via offline land model spinup. The approach forces a prognostic land model with a data atmosphere (DATM) (an "I compset" in CIME terminology) to create a land state in balance with observed atmospheric forcing. You can think of this as a "poor person's data assimilation" for the land model.

## Background

Betacast initializes the atmosphere, data ocean, and data ice from reanalysis/analysis (e.g., ERA5, GFS), but cannot directly interpolate land surface initial conditions. The recommended approach is to spin up the land surface model offline:

- For **successive Betacasts**, the framework automatically uses the +12 hr land restart from the previous run (stored in `${CASE}/run/landstart/`).
- For a **new Betacast** (or a new set of hindcasts), an initial land state must be generated using the tools here.

## Quick Start

### 1. Run `auto-script.sh`

```bash
cd $BETACAST/land-spinup
./auto-script.sh MODELSYSTEM DATAFORCING DATEYYYYMMDD NMONTHS NCYCLES ANOMYEAR REFYEAR MACHNAMELIST
```

| Argument | Description | Required |
| --- | --- | --- |
| `MODELSYSTEM` | `0` = CLM (CESM), `1` = ELM (E3SM/E3SMv2+) | Y |
| `DATAFORCING` | `0` = ERA5, `1` = CRU/NCEP (default), `2` = CAM/E3SM model output, `3` = CR20v3 | Y |
| `DATEYYYYMMDD` | Target date for the restart file (YYYYMMDD or YYYYMMDDHH) | Y |
| `NMONTHS` | Number of months of spinup prior to `DATEYYYYMMDD`; use `0` to use the full available DATM period | Y |
| `NCYCLES` | Number of cycles over those months (≥1; total model time = NMONTHS × NCYCLES) | Y |
| `ANOMYEAR` | Year of climate anomalies/deltas to apply; use `-1` for no anomalies | Y |
| `REFYEAR` | Reference year to correct anomalies; use `-1` for raw anomalies (ignored if `ANOMYEAR < 0`) | Y |
| `MACHNAMELIST` | Path to machine/experiment namelist file | Y |

**Supported `DATAFORCING`:**

| Value | Dataset
| --- | ---
| `0` | ERA5
| `1` | CRU/NCEP (CRUNCEP, internal)
| `2` | CAM/E3SM model output
| `3` | 20th Century Reanalysis v3 (CR20v3)

### 2. Examples

**ERA5, CESM/CLM, historical (e.g., Hurricane Katrina):**

```bash
./auto-script.sh 0 0 20050828 24 1 -1 -1 nl.landspinup.derecho
```

Spins up CLM with ERA5 DATM for 24 months prior to August 28, 2005. No climate deltas.

**ERA5, E3SM/ELM, no anomalies:**

```bash
./auto-script.sh 1 0 19960113 12 1 -1 -1 nl.landspinup.pm-cpu
```

Spins up ELM with ERA5 DATM for 12 months prior to January 13, 1996. No climate deltas.

**ERA5, E3SM/ELM, with climate change deltas (PGW):**

```bash
./auto-script.sh 1 0 19960113 12 3 2080 1996 nl.landspinup.pm-cpu
```

Uses 12 months x 3 cycles = 36 months of spinup. Deltas from year 2080, corrected relative to 1996.

**CR20v3 forcing:**

```bash
./auto-script.sh 0 3 18620601 12 1 -1 -1 nl.landspinup.derecho
```

Spins up CLM with CR20V3 DATM for 12 months prior to June 6, 1862. No climate deltas.

**Simplest case (built-in CRU/NCEP forcing, no ERA5 needed):**

```bash
./auto-script.sh 1 1 19960113 12 3 -1 -1 nl.landspinup.pm-cpu
```

Spins up ELM with internal DATM files, 12 months x 3 cycles = 36 months of spinup

**NMONTHS vs. NCYCLES:** `NMONTHS` is the span of calendar months looped over; `NCYCLES` repeats that span. For example, `NMONTHS=24 NCYCLES=2` produces 48 months of spinup, cycling over the same 2-year window twice. Setting `NMONTHS=0` skips the spinup window calculation and uses the full available DATM period (useful for `BUILD_ONLY=true` to build a case for use with `drive-dates`).

### 3. Machine namelist

Copy and edit one of the existing `nl.landspinup.*` files for your system:

| Variable | Description | Required | Default |
| --- | --- | --- | --- |
| `CIMEROOT` | Path to CESM or E3SM root | Y | |
| `PATHTOCASE` | Directory where the I-compset case will be created | Y | |
| `ICASENAME` | Base name for the spinup case (date, spinup months, and anomaly year are appended automatically). May contain `RESSTRING` as a placeholder for `RESOL`. | Y | |
| `PROJECT` | HPC project/account charging ID | Y | |
| `MACHINE` | CIME machine name (e.g., `derecho`, `pm-cpu`) | Y | |
| `NNODES` | Number of nodes | Y | |
| `RESOL` | Model resolution (must match your Betacast run grid) | Y | |
| `RUNQUEUE` | Batch queue | Y | |
| `WALLCLOCK` | Wall clock time (e.g., `"06:00:00"`) | Y | |
| `BETACAST` | Absolute path to the Betacast repository | Y | |
| `BETACAST_DATM_FORCING_BASE` | Top-level directory of DATM forcing files. Must contain subdirectories `Solar/`, `Precip/`, and `TPQW/`. | Y* | |
| `BETACAST_DATM_ANOMALY_BASE` | Top-level directory of anomaly files. Must contain `ens_{FLDS,PRECT,QBOT,TBOT}_anom.nc`. | Y* | |
| `USER_FSURDAT` | Override default land surface dataset | | (model default) |
| `USER_FINIDAT` | Override default initial conditions | | cold start |
| `USER_ICOMPSET` | Override default I compset | | system default |
| `USER_JOB_PRIORITY` | Override job priority | | (queue default) |
| `BUILD_ONLY` | Set `true` to configure/build without submitting | | `false` |
| `FORCE_PURGE` | Set `true` to delete existing case + run dirs before recreating | | `false` |
| `FORCE_COLD` | Set `false` to skip forcing a cold start when no `USER_FINIDAT` is given | | `true` |
| `RUN_DIR_BASE` | Top-level scratch dir to purge when `FORCE_PURGE=true` | Y* | |
| `COUPLER` | `mct` or `nuopc` | | `mct` |

*`BETACAST_DATM_FORCING_BASE` is required when `DATAFORCING` = 0, 2, or 3. `BETACAST_DATM_ANOMALY_BASE` is required when `ANOMYEAR` >= 0. `RUN_DIR_BASE` is required when `FORCE_PURGE=true`.

### 4. Stage restart files for Betacast

After the I-compset run completes, copy the restart files to the `landstart/` subfolder of your Betacast case directory:

```bash
CASENAME=my-betacast-case
ICASENAME=my-spinup-case

LANDFILEDIR=/path/to/betacast/cases/${CASENAME}/run/landstart/
ICASEDIR=/path/to/scratch/${ICASENAME}/run/

mkdir -p ${LANDFILEDIR}
cd ${ICASEDIR}

# ELM restart (E3SM)
cp -v ${ICASENAME}.elm.r.*.nc    ${LANDFILEDIR}
cp -v ${ICASENAME}.mosart.r.*.nc ${LANDFILEDIR}   # optional, if runoff active

# CLM restart (CESM)
cp -v ${ICASENAME}.clm2.r.*.nc   ${LANDFILEDIR}
cp -v ${ICASENAME}.rtm.r.*.nc    ${LANDFILEDIR}   # optional, if runoff active

# Rename to match the Betacast CASENAME
cd ${LANDFILEDIR}
rename ${ICASENAME} ${CASENAME} *.nc
```

You may also store these in a generic `landrawdir` folder. See the top level Betacast directions.

## Generating DATM Forcing Files

### ERA5

ERA5 surface forcing must be downloaded and converted before it can be used as DATM input. Pre-staged ERA5 DATM files are available on Derecho/Glade and Perlmutter. Ask Colin for paths. To generate your own:

#### 1. Download ERA5 surface fields

```bash
cd ${BETACAST}/land-spinup/gen_datm/get-era5
module load conda ; conda activate npl   # (NCAR; adjust for your system)
# Edit driver-get-era5.sh for years and output path
nohup ./driver-get-era5.sh &
```

#### 2. Convert to streams

```bash
cd ${BETACAST}/land-spinup/gen_datm/gen-datm
# Edit driver-gen-datm.sh for years and paths
qsubcasper driver-gen-datm.sh   # (or equivalent scheduler command)
```

The output directory (`BETACAST_DATM_FORCING_BASE`) will contain three subdirectories:

```
Solar/    ← FSDS (shortwave down)
Precip/   ← PRECIPmms (precipitation)
TPQW/     ← TBOT, WIND, QBOT, PSRF, FLDS, ZBOT
```

### CR20v3

Generating CR20v3 DATM files follows the same structure:

```bash
cd ${BETACAST}/land-spinup/gen_datm/gen-datm
# Edit driver-gen-datm-cr20v3.sh
# Domain file: cr20v3-domain.nc (already present in this directory)
python gen-forcing-cr20v3.py
```

### CESM/E3SM model output as forcing

For counterfactual or model-to-model experiments, DATM files can be generated from model history output:

```bash
# Edit single-file-to-datm.sh (pay attention to MAPFILE and OUTDIRBASE)
${BETACAST}/land-spinup/datm-from-model/single-file-to-datm.sh

# Batch convert all files (uses GNU parallel if available)
${BETACAST}/land-spinup/datm-from-model/batch-model-to-datm.sh
```

Then invoke `auto-script.sh` with `DATAFORCING=2` and set `BETACAST_DATM_FORCING_BASE` to the output directory:
```bash
./auto-script.sh 1 2 19860101 12 1 -1 -1 nl.landspinup.pm-cpu
```

## Deltas/Anomalies

Anomaly files representing a climate change signal can be overlaid on ERA5 DATM. Four files are required in `BETACAST_DATM_ANOMALY_BASE`:

```
ens_FLDS_anom.nc
ens_PRECT_anom.nc
ens_QBOT_anom.nc
ens_TBOT_anom.nc
```

These are currently derived from CESM1 LENS ensemble mean anomalies but can be replaced with any equivalent files. When `REFYEAR > 0`, `normalize-datm-deltas.ncl` is called to normalize the anomalies relative to a specific reference year before applying them.

Set `ANOMYEAR` to the target future year and `REFYEAR` to the historical reference year:

```bash
./auto-script.sh 1 0 19960113 12 3 2080 1996 nl.landspinup.pm-cpu
```

> **Note:** Anomalies are drawn from CESM1 LENS data which starts in 1920. Setting `ANOMYEAR >= 0` requires that the spinup start date is no earlier than 1920.

## Bias Correction (`bias-correct/`)

In some situations, it may be useful to apply a spatially/temporally varying bias correction to the ERA5 DATM during spinup, for example, to correct a known warm or dry bias over a region of interest (Hi, Alan!). This reuses the DATM anomaly mechanism.

See `bias-correct/bias_instructions.md` for the full workflow. In brief:

1. Create bias-correction NetCDF files:
   - **Temperature, humidity, winds, pressure**: additive (e.g., `+1.0` = +1 K)
   - **Precipitation, SW/LW fluxes**: multiplicative ratio (e.g., `1.1` = +10%)
2. Run `auto-script.sh` with `BUILD_ONLY=true` to build but not submit.
3. Copy the stream files from `bias-correct/` into the case directory and edit them to point to your bias files.
4. Edit `user_nl_datm` in the case to reference the new anomaly streams.
5. Run `./case.submit`.

Example bias files and stream templates are provided in `bias-correct/`. Currently, this is a pretty "manual" process.

## Sequential Spinup for Many Dates (`drive-dates/`)

NOTE: This should be considered partially deprecated; see the next section.

For generating land restart files across many target dates in sequence (e.g., for a large hindcast ensemble), use the `drive-dates` workflow instead of running `auto-script.sh` repeatedly.

### Setup

First, use `auto-script.sh` with `BUILD_ONLY=True` and `NMONTHS=0` to create and build the case without running it:
```bash
./auto-script.sh 0 0 20020101 0 1 -1 -1 nl.landspinup.tclf.ERA5
```

Then run `drive-dates.sh` with a namelist that points to the built case and a file listing the target dates:

```bash
cd drive-dates/

# Run interactively (recommended; runs until all dates are complete)
nohup ./drive-dates.sh nl_dd.ERA5.ne30pg3 &

# Or run in the foreground
./drive-dates.sh nl_dd.ERA5.ne30pg3
```

`drive-dates.sh` is self-resubmitting: after each successful run it archives the restart, removes the completed date from the dates file, and calls itself again. It runs until the dates file is empty.

### Dates file

The `datesfile` is a chronologically ordered text file with one date per line (YYYYMMDDHH):

```
2003091612
2003091700
2003091712
2004081000
```

### `drive-dates` namelist variables

| Variable | Description | Required | Default |
| --- | --- | --- | --- |
| `CASESRC` | Directory containing the built I-compset case | Y | |
| `CASENAME` | Name of the case | Y | |
| `BASERUN` | Top-level scratch directory (parent of `CASENAME/run/`) | Y | |
| `BETACAST` | Path to Betacast | Y | |
| `datesfile` | Path to the dates file | Y | |
| `DIRSTASH` | Directory to archive restart files after each run | Y | |
| `WALLCLOCK` | Wallclock for standard runs | Y | |
| `RUNQUEUE` | Queue for standard runs | Y | |
| `RUNPRIORITY` | CIME job priority (`JOB_PRIORITY`) | | |
| `SHORTCLOCK` | Wallclock for short runs | | `$WALLCLOCK` |
| `SHORTQUEUE` | Queue for short runs | | `$RUNQUEUE` |
| `SHORTCUTOFFHRS` | Hours threshold between short and standard runs | | `120` |
| `CIMEsubstring` | Additional substring argument for CIME submission | | `""` |
| `CIMEbatchargs` | Additional batch arguments for CIME submission | | `""` |
| `CIMEMAXTRIES` | Max retries if a run fails before aborting | | `3` |

See `drive-dates/README.md` for full documentation.

## Linear Spinup with Inline Restart Writes (`inline-dates/`)

A (much better!) alternative to `drive-dates` for generating many restart files is the **inline** approach: configure and run a single continuous I-compset integration, but patch the land/runoff model so it writes restarts at user-specified times during the run. This can be more efficient than `drive-dates` for dense sets of dates over a short period.

### Workflow

1. Configure a case with `auto-script.sh` (with `BUILD_ONLY=true`).
2. Apply the Fortran patches in `inline-dates/SourceMods/` to the land and (optionally) runoff model, then clean/rebuild.
3. Create `restart_times.txt` in the run directory — a chronologically ordered ASCII file with one YYYYMMDDHH per line:
   ```
   2002080412
   2002080500
   2002090212
   ```
4. Run the model normally (set `RESUBMIT`, `CONTINUE_RUN`, etc. as needed).
5. After the run, archive using `inline-dates/compress-and-archive.sh`:
   - Compresses all but the last restart with zstd (Betacast decompresses at runtime).
   - Moves compressed files to a final archive directory.

See `inline-dates/README.md` for full details.