# Inline dates

1. Configure case using auto-script.sh.
2. Patch land and runoff models, build (or rebuild if built in step 1).
3. Generates `restart_times.txt`, put in $SCRATCH/run/ directory.
4. Run model linearly (i.e., set RESUBMIT, CONTINUE\_RUN, etc. as desired).
5. Archive files using `compress-and-archive.sh`.

### `restart_times.txt`

This is just an ASCII file with chronological restart (i.e., land and runoff initial conditions) dates in YYYYMMDDHH format.

Example:

```
2002080412
2002080500
2002090212
2002090300
2002090312
2002090500
2002090512
2002090600
2002091200
2002091212
2002091300
```

### Archive files

In general, CLM and ELM initial conditions are large but very compressible. Since some model builds don't play nicely with chunked/compressed NetCDF, we just use zstd to compress them "offline" -- Betacast will uncompress at runtime.

`compress-and-archive.sh` will take all but the last restart file (needs to be saved in case it is needed for resubmit), zstd them using Betacast utils, and move them to an archival folder of your choice. This can be used in the Betacast namelist to draw initial conditions from later.

### Example auto-script namelist

```
CIMEROOT=/global/cfs/cdirs/m2637/E3SM_SCREAM_files/src/E3SM-20251110/
PATHTOCASE=/global/homes/c/czarzyck/I-runs/
ICASENAME=ERA5-ELMinic-conus256
PROJECT=m2637
MACHINE=pm-cpu
NNODES=8
RESOL=tc003lf192x8pg2_conustight256x8pg2
RUNQUEUE=regular
WALLCLOCK="03:00:00"
USER_FSURDAT="/global/homes/c/czarzyck/m2637/E3SM_SCREAM_files/fsurdat/surfdata_conus-tight_256x8_pg2_simyr2010_c240215.nc"
BETACAST=/global/homes/c/czarzyck/betacast/
BETACAST_DATM_FORCING_BASE=/pscratch/sd/c/czarzyck/DATM_FORCING/ERA5/
BUILD_ONLY=True
FORCE_PURGE=True
RUN_DIR_BASE=/pscratch/sd/c/czarzyck/e3sm_scratch/pm-cpu/
```

Run for "0" spinup months starting on the day of the linear initialization. Ideally, this is many months (or even a few years) beyond the first entry in `restart_times.txt`.
```
./auto-script.sh 1 0 20020101 0 1 -1 -1 nl.landspinup.pm-cpu.ERA5
```
