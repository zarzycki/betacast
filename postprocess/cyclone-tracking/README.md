# Cyclone tracking in Betacast

### Overview

Two main steps:

1. Acquire TCVitals from NHC or other source.
2. Low-bar pass of TC tracks using TempestExtremes.
3. (Optional) Calculate statistics with MET-TC.

### Acquire TC Vitals file

First, we need to figure out exactly where the operational centers deemed observed TCs to be at the same time as the Betacast initialization. These can be gleaned from a couple different sources, with the below bash script attempting to figure out the correct source.

```
bash ./get-vitals.sh ${YYYYMMDDHH} ${TCVITFOLDER}
```

where: `${YYYYMMDDHH}` is the forecast cycle and `${TCVITFOLDER}` is the folder where the vitals file will be saved.

Example:

```
bash ./get-vitals.sh 2018082100 ./fin-tcvitals/
```

Will procure the TC Vitals file for Aug, 21. 2018 at 00Z and store it within a sub-folder named `fin-tcvitals`. Each init time will get a single file where there will be N rows consisting of the N observed cyclones globally at that time.

#### TC Vitals for older storms

In `cyclone-tracking/ebt-to-vits`, there is a script `convert-ebt-to-tcv.py` that will convert [Extended Best Track Data](https://rammb2.cira.colostate.edu/research/tropical-cyclones/tc_extended_best_track_dataset/) to TC Vitals. Importantly, this is not the historical Vitals archive, but really a surrogate used for tracking and verification, but actual Vitals files may have somewhat different values, particularly for storm center and wind radii.

To invoke, just pass the desired year to the script:

```
python convert-ebt-to-tcv.py 2003
```

### Track all cyclones and match to TC Vitals

Requirements: CAM or EAM data with PSL, UBOT, and VBOT on the files at <= 6-hourly time frequency. Code only tracks at 6-hourly frequency (to match WMO operations) so if higher time frequency is passed in (e.g., 3-hourly) the `TIMESTRIDE` argument below must be updated.

Then `./cyclone-tracking-driver.sh` can be invoked with the following format (11 (optional 12) command line arguments!).

```
bash cyclone-tracking-driver.sh \
	${YYYYMMDDHH} \
	${CASENAME} \
	${TCVITFILE} \
	${ATCFFILE} \
	${TE_CONNECTFILE} \
	${OUTPUTBASE} \
	${SENDHTML} \
	${HSTREAM} \
	${TIMESTRIDE} \
	${ATCFTECH} \
	${TE_NOMPI_DIR}
```

| Namelist Variable | Description |
| --- | --- |
| YYYYMMDDHH | Forecast cycle as 10-digit integer |
| CASENAME | Name of the Betacast case (ex: conus_30_x8_forecast) |
| TCVITFILE | Full path to TC Vitals file for matching |
| ATCFFILE | Full path to desired output ATCF file |
| TE_CONNECTFILE | Connectivity file required by TE |
| OUTPUTBASE | Path to Betacast case base (i.e., fullpath before date folder) |
| SENDHTML | Send files to Colin's server? Generally false. |
| HSTREAM | Which "h" stream for EAM/CAM contains the tracker variables required by TE? |
| TIMESTRIDE | Do I stride TE files? Should be 6 for hourly, 2 for 3-hourly, 1 for 6-hourly |
| ATCFTECH | Shortcode to use in ATCF files (ex: op centers using GFSO, HWRF, etc.) |
| TE_NOMPI_DIR | Top-level path to a serial version of TempestExtremes |
| ARCHIVEDIR | (optional) If set, archive candidates, tracks, and ATCF with nc files |

After this script is successful, an [ATCF](https://en.wikipedia.org/wiki/Automated_Tropical_Cyclone_Forecasting_System) file for this particular forecast cycle will be generated in `./fin-atcf`. If there is no ATCF for this particular cycle, a new one will be generated. If one already exists, any lines in the ATCF that match `ATCFTECH` will be purged, with new lines appended. This allows for serial appending of various ensemble members (e.g., CAM001, CAM002, CAM003, etc.).

#### Only track cyclones

The code is modularized so that one could run the tracking code only without the ATCF matching. This may be desirable if you want an operational track without a vitals file. To do so "offline" follow this.

```
bash run-te-tracking.sh \
  ${TE_CONNECTFILE} \
  ${OUTPUTBASE} \
  ${YYYYMMDDHH} \
  ${CASENAME} \
  ${HSTREAM} \
  ${TE_NOMPI_DIR} \
  ${TIMESTRIDE} \
  ${TRAJFILE_TO_OUT}
```

NOTE: `${CASENAME}` is really just a unique string here. `${TRAJFILE_TO_OUT}` is the resulting TE trajectory file.

An example in action is:

```
bash run-te-tracking.sh \
  "/glade/u/home/zarzycki/work/grids/connect/mpasa3-60-tclf009_connect.dat" \
  "/glade/derecho/scratch/zarzycki/SYNTH-tclf009-mp3a-L32/run/" \
  2003083100 \
  "MPAS" \
  "h3i" \
  "/glade/work/zarzycki/tempestextremes_noMPI/" \
  2 \
  traj.txt.MPAS
```

However, note that the code tracks *all* cyclonic storms (no warm core threshold, very lax PSL threshold, etc.) and is not well suited to tracking only TCs.

#### Generating a connectivity file (TE_CONNECTFILE) for TempestExtremes

For TE to work on an unstructured mesh (e.g., SE, HOMME, MPAS), you must generate a connectivity file using TE's included binary that tells how each cell is "attached" to its neighbors. Example:

```
cd ${TE_NOMPI_DIR}/bin/
./GenerateConnectivityFile --in_mesh $GRID.g --out_type CGLL --out_connect $GRID_connect_v2.dat
```

or

```
cd ${TE_NOMPI_DIR}/bin/
./GenerateConnectivityFile --in_mesh $GRID.scrip.nc --out_type FV --out_connect $GRID_connect_v2.dat
```

If you have a SCRIP file even for SE/HOMME grids, the latter is preferred for consistency.

### Calculate statistics with MET

An optional (but common) step is to calculate a variety of operational statistics with DTC's [Model Evaluation Tools](https://dtcenter.org/community-code/model-evaluation-tools-met). This is beyond the scope of this README right now, but the ATCF that comes out of this code can be used as any a-deck file would be for any model. See the above link for more details.

A sample set of scripts is provided in `./met-tc/`. After navigating to that folder, the basic workflow for calculating error statistics is:

1. Move ATCF files from `./fin-atcf/` and concatenate into a single a-deck file. This contains information about the forecasted tracks of TCs from Betacast.
2. (Optional) Acquire a-deck files and concatenate to model ATCF data.
3. Acquire b-deck files.
4. Gather TC Pairs using `./tc-pairs.sh`.
5. Calculate TC Stats using `./tc-stats.sh`.
6. Analyze resultant stats output with Python, etc.

#### Acquire a-deck files

A-deck files represent model forecasts of various TCs and can be acquired from a variety of sources. One example is the RAL repository. All a-deck files for a given hurricane season (e.g., 2021) can be downloaded via wget (O(100-200MB)).

```
wget -r --no-parent -nH --cut-dirs 3 -A 'aal*2021.dat' http://hurricanes.ral.ucar.edu/repository/data/adecks_open/2021/
```

These can be concatenated into one seasonal file (or even multiple seasons added to the same file).

```
mv 2021 adeck-2021
cd adeck-2021
cat aal*.dat > aal.2021.ALL
```

#### Acquire b-deck files

B-deck files include the 'best-track' information (i.e., observations) about TCs and can be also be acquired from a variety of sources. All b-deck files for a given hurricane season (e.g., 2021) can be downloaded via wget (O(1MB)).

```
wget -r --no-parent -nH --cut-dirs 3 -A 'bal*2021.dat' http://hurricanes.ral.ucar.edu/repository/data/bdecks_open/2021/
```

Likewise, these can be concatenated as follows:

```
mv 2021 bdeck-2021
cd bdeck-2021
cat bal*.dat > bal.2021.ALL
```

#### Using tc_pairs

This is a MET-TC tool. Generally, my workflow is to copy `TCPairsConfig_Default` from the MET repository for that particular version (e.g., v10.0). Then use that as my starting point to select storms, filter time windows, etc.

Running `./tc-pairs.sh` should generate a file that contains filtered and/or matched information from the a-deck and b-deck files that is required for calculating statistics.

NOTE: If you get an error like the below (experienced in v7.0)

```
WARNING:
WARNING: int ATCFTrackLine::read_line(LineDataFile * ldf) -> found fewer than the expected number of elements (6) in ATCF track line:
WARNING: AL, 05, 2021063018, 03, HWRF, 078,
WARNING:
```

This is *not* a warning and the malformed lines need to be deleted!

#### Using tc_stat

This is also a MET-TC tool. Here, I configure everything from the command line for basic statistics. See `./tc-stats.sh` for more information. The resulting file can be loaded in Python or Excel and one can filter by lead time, model, etc. to evaluate things like track error and intensity bias.
