## Generating a CLM initial condition for a CESM forecast

This README details the use of a CESM "I" compset (active land, data atmosphere) to "nudge" CLM into a state balanced to the atmosphere. This is done to generate a consistent CLM initial condition that can be used with a CAM initial data file derived from observations or reanalysis to produce an initialized simulation akin to a forecast run.

-------------------------------------------------------------------------------------------------------------------------

Since CLM has the capability to initialize by interpolating from other grids, this is optional (unless you want to run CLM4 in which case it's less optional)

1. Generate MPAS SCRIP file (ex: mp120a_scrip.nc)
2. Generate fsurdat file for relevant years (using SCRIP file)
3. Add grid to CESM via config_grids.xml in CIME (need req'd domain, mapping, runoff, etc. files as well.)
4. Follow the directions as below, replacing ne120 with your user-defined grid.

-------------------------------------------------------------------------------------------------------------------------

### The easy way (CLM5 and beyond)

Since CLM forcing for DATM is ~0.5deg and CLM5 has the capability to interpolate between grids during initialization there is not a tremendous amount of utility in running a very high-resolution CLM grid for spinup. The one benefit is some CLM land cover type data (fsurdat) is O(~10km). However, anecdotally, for most weather purposes, completing a 0.25deg initialization and then interpolating it to a CLM grid in forecast mode is adequate (personal communication with Ahmed Tawfik/Erik Kluzek/Dave Lawrence in CGD/TSS)

Preferably using the CESM tag that you would like to use for initialized forecasts…


* Ensure that there is a valid compset for CRU forcing with 2000 atmospheric conditions and the version of CLM you would like to use (CLM5-SP for forecasts with CESM2). This exists for CLM4 right now, but not CLM5 (as of cesm2_0_alpha07)
  * If a compset with CRU and CLM5-SP exists, ignore.
  * To add, open $CESMROOT/components/clm/cime_config/config_compsets.xml and add:

```
<compset>
<alias>I2000CLM5</alias>
<lname>2000_DATM%CRU_CLM50%SP_SICE_SOCN_RTM_SGLC_SWAV</lname>
</compset>
```

Now create a compset using the 0.25deg version of CLM5 (based on the ne120 SE grid).

`$ ./create_newcase -case ~/I-compsets/I2000CLM5_ne120 -res ne120_ne120 -compset I2000CLM5 -mach cheyenne --project 1111111111`

Navigate to your CASE directory, edit env_mach_pes.xml (optional), and then run `./case.setup`. The below PE configuration runs approximately 3.5 SYPD on Cheyenne:

After `./case.setup`, the user_nl_* files will be populated. NOTE: If you have added a new grid to CLM/CESM, you may get an error during this stage if you did not specify required datasets such as fsurdat and finidat in CLM’s default namelists. You can add them in user_nl_clm and re-run `./case.setup`.

The main modification is copying over updated user-defined DATM streams with updated CRU data (through the end of 2016 currently). On glade, these are currently found in /glade/u/home/zarzycki/cam_cesm_files/DATM-forcing/CLMCRUNCEP.20170901.fromTawfik. From the case dir, just run:

`cp /glade/u/home/zarzycki/cam_cesm_files/DATM-forcing/CLMCRUNCEP.20170901.fromTawfik/user_datm* .`

To be safe, you should also add end_restart=.true. in user_nl_cpl. This is not needed if you are careful in env_run.xml, but ensures the model will dump restart (initial) files at the termination of the run, regardless of how end_run is called.

Now all that is left is to modify env_run.xml. From personal discussions with Ahmed Tawfik, soil moisture and temperature are the two primary quantities of interest for short-term forecasts and do not require a great deal of spinup. It is helpful to nudge through multiple seasons, so we have agreed the most straightforward thing is to initialize DATM CLM on Jan 1st of the calendar year BEFORE the forecast date of interest.

Example forecast date of 5/6/2015 (we initialize DATM/CLM on 1/1/2014)

Settings for all forecasts:

```
./xmlchange STOP_N=3
./xmlchange STOP_OPTION='nyears'
./xmlchange DATM_CLMNCEP_YR_END=2016
./xmlchange REST_OPTION='end'
```

CESM will stop when it *first* hits *either* STOP_N and STOP_OPTION *or* STOP_DATE, so the STOP_N/STOP_OPTION needs to be beyond the forecast date of interest. DATM_CLMNCEP_YR_END just says that the end of the DATM stream is 2016.

Settings forecast-date specific:

```
./xmlchange DATM_CLMNCEP_YR_ALIGN=2014
./xmlchange DATM_CLMNCEP_YR_START=2014
./xmlchange RUN_STARTDATE=2014-01-01
./xmlchange STOP_DATE=20150506
```

RUN_STARTDATE is the initial date for the spinup run. The two ALIGNS are set to the CLM initialization date (Jan 1 of the calendar year prior to year CLM restart file is needed for forecast initialization). STOP_DATE is set to the date the forecast condition is needed. The model will stop and print a CLM restart file on this date.

Run the model using `./case.submit`. It should integrate for 1+ years, stopping on STOP_DATE and printing model restart files. The file *clm.r.* can be used as a CLM initial condition for coupled forecast runs.

To run a forecast with this CESM initial file, you will need to add the following to your user_nl_clm (in the coupled forecast case, NOT the DATM case). Note that use_init_interp is needed if finidat is on a different grid than the CLM grid used in the forecast compset (ex: using the ne120 grid here to generate an initial data file for an MPAS CLM grid, for example)

```
finidat="pathtomyfile/sample.ne120_ne120.clm2.r.2015-05-06-00000.nc"
use_init_interp=.true.
```

