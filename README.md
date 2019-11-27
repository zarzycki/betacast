# Betacast

### General workflow
1. Create a case directory with either a supported or new grid configuration and verify that it is stable/runs given arbitrary inputs.
2. Build an analysis/reanalysis -SE weight file.
3. Edit machine file for your particular system
4. Edit namelist file for your particular use case.
5. Edit output streams file.
6. If running historical, set dates.txt

The first step is to create a functional F compset. Broadly, this is active atmosphere, active land, and data ocean/ice. Other configurations may work (e.g., active runoff, wave, glacier models) but have not been tested.

**E3SM:**

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


