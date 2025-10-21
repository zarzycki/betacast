# Drive dates!

This script automates advancing and resubmitting data atmosphere CLM/ELM based on sequential dates defined in a namelist file. This is used to generate restart files for the land and runoff models that can be used as initial conditions for CESM/E3SM.

## Usage

```bash
./drive-dates.sh <NAMELISTFILE>
```

### Arguments

| Argument         | Description                                                                                                                                                            |
| ---------------- | ---------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `<NAMELISTFILE>` | Path to a bash-style namelist file defining key variables such as `CASESRC`, `CASENAME`, `BASERUN`, `BETACAST`, `datesfile`, `DIRSTASH`, `WALLCLOCK`, `RUNQUEUE`, etc. |

### Example Invocations

First you must configure an I-compset using `auto-script`. Set to the date you want to start on and set spinup months to 0. Set `BUILD_ONLY=True` so the script sets up everything but doesn't submit the initial run.

```
./auto-script.sh 0 0 20020101 0 1 -1 -1 nl.landspinup.tclf.ERA5 
```

The build your namelist (start with example or see below)

```bash
# Run interactively
./drive-dates.sh nl_dd.ERA5.ne128pg2

# Run in the background
nohup ./drive-dates.sh nl_dd.ERA5.ne128pg2 &
```

### Namelist Variables

| Variable        | Description                                                                                    | Default (if not set) | Required |
| --------------- | ---------------------------------------------------------------------------------------------- | -------------------- | -------- |
| `CASESRC`       | Base directory path where case source directories are stored. (e.g., `/global/homes/c/czarzyck/I-runs/`)      | —                    | ✅        |
| `CASENAME`      | Name of the case (e.g., `ERA5-ELMinic_2002010100_000_01`).              | —                    | ✅        |
| `BASERUN`       | Base directory for model run. Should not include `CASENAME/run`. (e.g., `/pscratch/sd/c/czarzyck/e3sm_scratch/pm-cpu/`)        | —                    | ✅        |
| `BETACAST`      | Root directory of the Betacast workflow repository. | —                    | ✅        |
| `datesfile`     | Path to the text file listing simulation dates for sequential runs.                                 | —                    | ✅        |
| `DIRSTASH`      | Directory where restart and log files are archived (copied/compressed) after runs.              | —                    | ✅        |
| `WALLCLOCK`     | Default long wallclock time for full-length runs.                                              | —                    | ✅        |
| `RUNQUEUE`      | Default long-queue name for full-length runs.                                                  | —                    | ✅        |
| `RUNPRIORITY`   | Optional CIME job priority to apply (`JOB_PRIORITY`).                                          | none                 | ⚪        |
| `SHORTCLOCK`    | Wallclock for short runs (`<=120 h`). Useful if debug/short queue available. Defaults to `$WALLCLOCK` if unset.                       | `$WALLCLOCK`         | ⚪        |
| `SHORTQUEUE`    | Queue for short runs (`<=120 h`). Useful if debug/short queue available. Defaults to `$RUNQUEUE` if unset.                            | `$RUNQUEUE`          | ⚪        |
| `CIMEsubstring` | Optional substring argument for CIME job submission.                                           | `""`                 | ⚪        |
| `CIMEbatchargs` | Optional batch arguments for CIME run submission.                                              | `""`                 | ⚪        |
| `CIMEMAXTRIES`  | Maximum number of retry attempts if CIME run fails.                                            | `3`                  | ⚪        |


### The "datesfile"

The datesfile is a **chronologically** ordered textfile containing the target dates for initialization over the integration period. For example:

```
2003091612
2003091700
2003091712
2004081000
2004081012
```


