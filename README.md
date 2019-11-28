# Betacast

### General workflow
1. Create a case directory with either a supported or new grid configuration and verify that it is stable/runs given arbitrary inputs.
2. Build an analysis/reanalysis -SE weight file.
3. Edit machine file for your particular system
4. Edit namelist file for your particular use case.
5. Edit output streams file.
6. If running historical, set dates.txt

The first step is to create a functional F compset. Broadly, this is active atmosphere, active land, and data ocean/ice. Other configurations may work (e.g., active runoff, wave, glacier models) but have not been tested.

### 1. Create case directory and test stability.

This step just requires a built and tested case of CESM. Any source mods should be made directly to this case. It is only necessary to show the model is stable for a few days.

**CESM example:**

```
./case.setup
### patch land mods for restart
./case.build
./case.submit
```

**E3SM example:**

```
cd ~/E3SM/cime/scripts
./create_newcase --case ~/F-betacast-FC5AV1C --compset FC5AV1C --res ne30_ne30 --mach cori-knl --project m1637
cd ~/F-betacast-FC5AV1C
./xmlchange CHARGE_ACCOUNT=m1637
./xmlchange NTASKS=-8
./xmlchange NTASKS_ESP=1
./case.setup
### patch land mods for restart
./case.build
./case.submit
```

Some notes:
1. In the "patch" step, a small modification is made to the land model to enforce restart files to be printed every 12 hours. This is done since the land model is initialized via nudging with the atmosphere in this framework. This patch can be applied by copying CESM or E3SM's `lnd_comp_mct.F90` (from the land model source code) into `$CASEDIR/SourceMods/src.clm` (or equivalent) and running `patch lnd_comp_mct.F90 < ${PATCHDIR}/lnd_comp_mct.patch` over the top of the file, which injects the correct logic. This needs to be done before the `./case.build` step.

### 2. Generate an analysis/reanalysis to CAM weight file

This can be done with `${BETACAST}/remapping/gen_GFS_to_SE_weight_file.ncl`. This script requires two inputs that are directly modified in the script body, `srcGridFile` (a SCRIP grid file defining the GFS regular lat-lon grid) and `dstGridFile` (a SCRIP grid file defining the destination CAM grid).

Historically there have been two analysis grid sizes associated with publicly disseminated GFS/CFS/CFSR analyses, 0.5degree and 0.25degree.

### 3. Edit machine file for your particular system

In `${BETACAST}/machine_files` there are sample files that define where folders and data files will be stored for your system. There are suggested configurations for Cheyenne and Cori-NERSC, but you may edit these for your workflow or copy/paste for a differnt system (i.e., university cluster).

| Namelist Variable | Description |
| --- | --- |
| path_to_case | Path to CESM "home" case directory |
| path_to_inputdata | Path to where re/analysis + model initial conditions/forcing data is stored |
| path_to_rundir | Path (top-level) to directory where CESM actively runs |
| sewxscriptsdir | Path to betacast repo (i.e., `${BETACAST}`) |

### 4. Edit namelist file for your particular case

### 5. Edit output streams





