# Nudging notes

`gen-nudge.sh` uses a bash namelist. If one is submitting to a queue, it should be passed in as an environmental variable named `NLFILE`.

```
qsubcasper -v NLFILE="ndg.era5.nl" gen-nudge.sh
```

If run on the command line, the code will look for the first command line option.

```
./gen-nudge.sh ndg.era5.nl 
```

SLURM can handle

```
sbatch gen-nudge.sh ndg.era5.nl 
```

| Namelist Variable | Description | CAM_TO_CAM |
| --- | --- | --- |
| CAM_TO_CAM | Are we running the "downscaling" code for CAM? If **false** -> re/analysis | |
| dryrun | If **true** run through the code and print commands file without actually running | |
| STYR | Start year (integer) of data to generate nudging for | |
| ENYR | End year (integer) of data to generate nudging for | |
| BETACASTDIR | Path to top-level Betacast source directory | |
| OUTDIR | Base output folder for nudging data to be output to | |
| DYCORE | What is our target dycore? Valid strings are "mpas" "se" and "fv" | |
| GRIDSTR | Name of the target grid, used for naming files and output directories | |
| BNDTOPO | bnd_topo file for *target* grid for hydrostatic correction | |
| WGTNAME | ESMF weight file that is regular lat-lon to target model grid | |
| NUMLEVS | Vertical levels in the model | |
| DESCSTR | String describing source data (e.g., "CAM5" "ERA5" "CFSR") | |
| HR_RES | Integer representing hours between nudging files (ex: 1 = hourly, 6 = 4xdaily) | |
| STMON | Short string representing start month to generate nudging for (e.g., "aug" "sep" "oct") | |
| ENMON | Short string representing end month to generate nudging for (e.g., "aug" "sep" "oct") | |
| STDAY | Integer representing start day to generate nudging for | |
| ENDAY | Integer representing end day to generate nudging for | |
| SUBNAME | Name of source CAM simulation in CAM_TO_CAM | x |
| BINLIST | Input file list (string w/ allowed wildcards) that contain 3-D model outputs to source from | x |
| MODREMAPFILE | ESMF weight file from source model grid -> intermediate RLL grid needed for WGTNAME | x |
| MODINTOPO | bnd_topo file for *source* grid | x |

--

`BINLIST` will be expanded if the nudging start and end dates straddle multiple files containing the 3-D fields (e.g., T, Q, U, V, PS, etc.). For example:

```
BINLIST="/glade/u/home/zarzycki/scratch/cam_to_cam/raw_files/storm_1354/CHEY.VR28.NATL.WAT.CAM5.4CLM5.0.dtime900.002.cam.h2.2013-08-*-00000.nc"
```

will load all August 2013 files and scan for required times in any of those files.

--

For Hyperion runs (or ensembles) an option is to pass a placeholder variable name in the input and then apply a parameter substitution in the code. For example.

```
MODREMAPFILE=/glade/u/home/zarzycki/betacast/remapping/map_ne0np4natlantic\$GRIDLOWER.ne30x4_TO_era5_0.25x0.25_patc.nc
```

in the namelist and then a modified nudging driver to include:

```
MODREMAPFILE="${MODREMAPFILE//\$GRIDLOWER/$GRIDLOWER}"
```

will substitute whatever is saved in `$GRIDLOWER` in the shell script.

# Nudging to ERA5 on Cheyenne

To nudge to hourly ERA5 data with either FV or SE (or any SE-like unstructured grid, like FV3) we can use the RDA dataset already mounted to Cheyenne in conjunction with the [Betacast](www.colinzarzycki.com) code.

- Generate ERA5 -> CAM grid map file
- Run batch script looping over atm_to_cam.ncl to generate nudging forcing files

### Generate ERA5 to map file

Betacast's (cloned to `${BETACAST}`) atm_to_cam requires an ESMF mapping file to regrid data from a host analysis grid to the CAM model grid. The easiest way to generate these is by using SCRIP grid files and NCL, Python, or the ESMF binaries.

Existing SCRIP files for ERA5 and CFSR (0.5) are in: `${BETACAST}/remapping/anl_scrip`

New host scrip files can be generated using `${BETACAST}/remapping/gen_reglatlon_SCRIP.ncl`

For CAM, supported SCRIP files can be dug up in the CESM data repo. The places I generally look for already processed SCRIP files are in:

- `/glade/p/cesmdata/inputdata/atm/cam/coords/`
- `/glade/p/cesmdata/inputdata/lnd/clm2/mappingdata/grids/`

After this, you can edit `${BETACAST}/remapping/gen_GFS_to_CAM_weight_file.ncl` to generate a weight file (ex: map_HOST_to_CAM_patc.nc.

### Run shell script

After generating a map file, you can run `${BETACAST}/atm_to_cam/nudging/gen-nudge.sh` script on Cheyenne's high memory nodes (see preamble). This job uses GNU parallel to spawn 36 jobs in parallel to use atm_to_cam to take data from the native ERA5 repo and create "nudging" (i.e., initial condition) files at each hourly time (e.g., 00Z, 01Z, etc.).

The only files that need to be set are `BNDTOPO` (CAM model topography for that grid) and `WGTNAME` (reanalysis to CAM weight file generated in the previous step).

NOTE: The default setting of `compress_file=True` in atm_to_cam uses NetCDF compression on the nudging files. This works fine on CESM2 nudging files but these files *cannot* be used for `ncdata` (initial conditions). If you want to use a file generated by atm_to_cam for initial conditions, you must set `compress_file=False`.
