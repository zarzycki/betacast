# CR20V3 Data Acquisition

Scripts for downloading 20th Century Reanalysis V3 (CR20V3) ensemble member data from NERSC HPSS via Globus. Each year of data is approximately 6 TB.

## Prerequisites

- Globus CLI (`globus`) (can be installed with conda)
- GNU parallel (`parallel`) (can also be installed with conda)

```
conda create -n globus -c conda-forge globus-cli parallel
```

## Step 1: Activate and setup Globus command line on destination machine

```
globus login
# Follow prompts
# Verify you are logged in with...
globus whoami
```

## Step 2: Download tar files via Globus

Edit `spawn-globus.sh`.

- NCAR_GLADE_GLOBUS is the endpoint for the destination locale (e.g., NCAR Glade)
- DOWNLOAD_DIR is the target download (top-level) directory on the destination locale
- YYYY is the target year you are requesting.

```bash
bash spawn-globus.sh
```

This lists available files, filters them to the needed variables/levels for the specified year (via `_filter_globus.sh`), and submits a Globus batch transfer.

**Variables downloaded:** PRES, TMP2m (2D); HGT, SPFH, TMP, UGRD, VGRD at up to 17 pressure levels (3D).

## Step 3: Extract and clean up tar files

After the files are downloaded, we need to untar them.

Edit `untar-and-rm.sh` and set `TARDIR` to match the Globus destination directory, then submit:

```
qsub untar-and-rm.sh
```

Extracts all `.tar` files in `TARDIR` using GNU parallel (8 concurrent jobs) and deletes each tar file after successful extraction.

Note you can set the number of CPUs to ~2 and run from the login nodes, I presume.

## Step 4 (Debugging + optional): Compute ensemble mean

Needs NCO installed.

Edit `create-mean-casper.sh` and set `TARDIR` and `YEAR`, then submit:

```
qsub create-mean-casper.sh
```

Computes the ensemble mean of all members for each variable/time using `ncea`. Output files are named `*_mem999.nc` (mem999 = ensemble mean).

This is just for debugging if you want to evaluate whether Betacast using the ensemble mean of the individual members matches the ensemble mean provided by NOAA (or in RDA).
