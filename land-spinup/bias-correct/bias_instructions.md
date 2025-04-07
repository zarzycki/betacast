## Bias correction

There may be situations where it may be beneficial to apply a bias correction. For example:

- ERA5 DATM alone is too warm over a target region during the spinup period, producing low snow water equivalent for a snowmelt storyline.
- ERA5 DATM alone is biased dry for a specific event at the end of the spinup period, producing soil moisture that is biased low for a sequential flood storyline.

One can attempt to further improve the land surface initialization through "bias correction." To do this, we use the "anomaly" feature of DATM to add or subtract values to the target ERA5 meteorological fields (e.g., to "juice" or "dejuice" the observed meteorology used as forcing).

1. Create anomaly (bias correction) files for various fields (typically T and/or PRECT).
2. Run `land-spinup.sh` but set build-only to `true`.
3. Add `user_datm.streams.txt.Anomaly.Forcing.*` files to case directory. Edit to point to bias correction files.
4. Edit `user_nl_datm` to correctly point to the new anomaly data streams.
5. Run `./case.submit` to spawn land spinup run.

### Create anomaly files

These are NetCDF files that can be read by DATM that contain the anomalies. The can be spatially and temporally varying and contain corrections that are applied to the main DATM stream (e.g., ERA5). I have placed two example files for T and PRECT in `$BETACAST/land-spinup/bias-correct/`. T contains a temporally varying temperature adjustment over CONUS (0 elsewhere). PRECT contains a temporally and latitudinally varying precipitation adjustment in the eastern hemisphere that is postive in the tropics and negative in the polar regions. Both are just toy examples.

- For temperature this should be a linear correction. For example, if you want the May temperatures to be +1K warmer than ERA5, the anomaly file should be = 1K during the relevant May timestep(s).
- For precipitation this should be a "ratio" correction. For example, if you want the August precipitation to be 10% higher than ERA5, the anomaly file should be = 1.1 during the relevant August timestep(s).

The examples I have provided are just an annual cycle correction centered at month midpoints (the time unit is just Jan 16th, Feb 15th, etc. of the year 1). You may apply daily (or even hourly) adjustments if you'd like by modulating the time dimension, same with multi-year adjustments. I have also provided an example on a nominally 1x1 degree RLL grid (CESM LENS), but you could apply on a finer/coarser grid as long as you have follow the same file format and also have a relevant domain file that can be pointed to in `user_datm.streams.txt.Anomaly.*`. Refer to CTSM/ELM documentation for more details on DATM streams.

### Run `land-spinup.sh` in "build only" mode

Build a spinup case as usual:

```
./auto-script.sh 1 0 19830524 12 3 -1 -1 nl.landspinup.pm-cpu
```

except make sure:

```
BUILD_ONLY=True
```

is added to the machine namelist (e.g., `nl.landspinup.pm-cpu`). This will do everything *except* submit the job for execution.

### Add user_datm.streams.txt.Anomaly.* files

You will need to add XML streams to your case directory that contain user information to point to the anomaly files. I have included examples for T and PRECT in `$BETACAST/land-spinup/bias-correct/` you can copy and edit accordingly.

### Edit user\_nl\_datm

In user_nl_datm, you need to add the `anomaly_forcing`, update `tintalgo` for the interpolation for the forcings (typically linear), and append the datm.streams.txt.* from the previous step to `streams`.

For example:

```
tintalgo = "coszen", "nearest", "linear", "linear", "lower"

streams = "datm.streams.txt.CLMCRUNCEPv7.Solar 1982 1982 1983",
    "datm.streams.txt.CLMCRUNCEPv7.Precip 1982 1982 1983",
    "datm.streams.txt.CLMCRUNCEPv7.TPQW 1982 1982 1983",
    "datm.streams.txt.presaero.clim_2000 1 1 1",
    "datm.streams.txt.topo.observed 1 1 1"
```

becomes

```diff
+ anomaly_forcing = 'Anomaly.Forcing.Precip','Anomaly.Forcing.Temperature'
  
  tintalgo = "coszen", "nearest", "linear", "linear", "lower"
+ , "linear", "linear"
  
  streams = "datm.streams.txt.CLMCRUNCEPv7.Solar 1982 1982 1983",
      "datm.streams.txt.CLMCRUNCEPv7.Precip 1982 1982 1983",
      "datm.streams.txt.CLMCRUNCEPv7.TPQW 1982 1982 1983",
      "datm.streams.txt.presaero.clim_2000 1 1 1",
      "datm.streams.txt.topo.observed 1 1 1"
+ ,
+     "datm.streams.txt.Anomaly.Forcing.Precip 1 1 1",
+     "datm.streams.txt.Anomaly.Forcing.Temperature 1 1 1"
```

Here, the 1 1 1 corresponds to:

- `DATM_YR_ALIGN`: Simulation year corresponding to `DATM_YR_START`. A common usage is to set this to `RUN_STARTDATE`. With this setting, the forcing in the first year of the run will be the forcing of year `DATM_YR_START`. Another use case is to align the calendar of transient forcing with the model calendar. For example, setting `DATM_YR_ALIGN=DATM_YR_START` will lead to the forcing calendar being the same as the model calendar. The forcing for a given model year would be the forcing of the same year. This would be appropriate in transient runs where the model calendar is setup to span the same year range as the forcing data.
- `DATM_YR_START`: Starting year to loop data over
- `DATM_YR_END`: Ending year to loop data over

In other words:

```
    In datm_in, streams namelist input has the form
     streams = 'stream1.txt year_align year_first year_last ',
               'stream2.txt year_align year_first year_last ',
                ...
               'streamN.txt year_align year_first year_last '
```

(above information shamelessly stolen from [DATM documentation](https://escomp.github.io/CDEPS/versions/master/html/datm.html))

### Confirming bias correction is applied correctly

If you are like me, and perpetually petrified your simulation "runs" but gives "garbage" answers, you can test your bias correction implementation by running two configurations -- one with the bias correction, one without, and comparing the forcing files output from the land model.

E.g.,

```
# Vanilla spinup case
./auto-script.sh 1 0 19830524 12 3 -1 -1 nl.landspinup.pm-cpu
# Perform above workflow, and now run again (storing in a new casedir)
./auto-script.sh 1 0 19830524 12 3 -1 -1 nl.landspinup.pm-cpu
```

Compare the output from the two cases. For example:

```
# Concatenate "raw" ERA5 DATM run
cd $SCRATCH/PM-CPU_19830524_012_03/run/
ncrcat -O -v RAIN,SNOW,TBOT PM-CPU*.elm.h1.1981-* cat.nc
# Concatenate "bias-corrected" ERA5 DATM run
cd ../../PM-CPU_19830524_012_03_ANOM/run/
ncrcat -O -v RAIN,SNOW,TBOT PM-CPU*.elm.h1.1981-* cat.nc
# Diff the bias correction minus raw
ncdiff -O cat.nc ../../PM-CPU_19830524_012_03/run/cat.nc diff.nc
# Remap to usable grid if needed
ncremap -i diff.nc -o tmp_new.nc -m ~/m2637/betacast/sewx/maps/map_ne30_to_1x1glob_patch.nc 
ncview tmp_new.nc
# Do analysis here
```