# Cyclone tracking in Betacast

---

### Overview

Two main steps:

1. Acquire TCVitals from NHC or other source.
2. Low-bar pass of TC tracks using TempestExtremes.
3. (Optional) Perform statistics with MET-TC.

### Acquire TCVitals file

```
/bin/bash ./get-vitals.sh ${YYYYMMDDHH} ${TCVITFOLDER}
```

where: `${YYYYMMDDHH}` is the forecast cycle and `${TCVITFOLDER}` is the folder where the vitals file will be saved.

Example:

```
/bin/bash ./get-vitals.sh 2018082100 ./fin-tcvitals/
```

Will procure the TCVitals file for Aug, 21. 2018 at 00Z and store it within a sub-folder named `fin-tcvitals`.


### Track all cyclones and match to TCVitals

Requirements: CAM or EAM data with PSL, UBOT, and VBOT on the files at <= 6-hourly time frequency. Code only tracks at 6-hourly frequency (to match WMO operations) so if higher time frequency is passed in (e.g., 3-hourly) the `TIMESTRIDE` argument below must be updated.

Then `./drive-tracking.sh` can be invoked with the following format (11 command line arguments!).

```
/bin/bash ./drive-tracking.sh ${YYYYMMDDHH} \
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

After this script is successful, an ATCF file for this particular forecast cycle will be generated in `./fin-atcf`. If there is no ATCF for this particular cycle, a new one will be generated. If one already exists, any lines in the ATCF that match `ATCFTECH` will be purged, with new lines appended. This allows for serial appending of various ensemble members (e.g., CAM001, CAM002, CAM003, etc.).

### Calculate statistics with ATCF data

An optional (but common) step is to calculate a variety of operational statistics with DTC's [Model Evaluation Tools](https://dtcenter.org/community-code/model-evaluation-tools-met). This is beyond the scope of this README right now, but the ATCF that comes out of this code can be used as any a-deck file would be for any model. See the above link for more details. 
