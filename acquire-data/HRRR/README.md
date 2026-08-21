# HRRR Data Download

`get-hrrr.py` downloads HRRR GRIB2 files from the NOAA HRRR AWS Open Data bucket ([registry.opendata.aws/noaa-hrrr-pds](https://registry.opendata.aws/noaa-hrrr-pds/)). The bucket is (as of Aug. 2026) public, so no AWS credentials are required — files are fetched anonymously over HTTPS with `requests`. This gets the analysis (f00).

### Usage

```bash
python get-hrrr.py --date YYYYMMDDHH --outdir OUTDIR [options]
```

#### Required arguments

- `--date` — cycle date/time in `YYYYMMDDHH` format, exactly 10 digits (e.g. `2024010100` for 2024-01-01 00Z).
- `--outdir` — directory to write the downloaded file(s) to (created if it doesn't exist).

#### Optional arguments

- `--product` — `wrfnat` (native/hybrid levels), `wrfprs` (pressure levels), `wrfsfc` (surface fields, incl. `APCP` precip accumulation), `both` (wrfnat+wrfprs, default), or `all` (wrfnat+wrfprs+wrfsfc).
- `--fxx` — forecast hour, `0`-`48` (default: `0`). For analyses this should always be `0`.
- `--region` — HRRR domain, e.g. `conus` or `alaska` (default: `conus`).
- `--skip-existing` — skip downloading a file if it already exists at the destination path (default: overwrite).

### Note on precipitation accumulation

`APCP` in `wrfsfc` is accumulated *since the cycle's forecast start*, so it is
~0 at `--fxx 0` (analysis). To get a nonzero precip accumulation you need
`--fxx` >= 1. Building a 6-hour accumulated precip field (as needed for
`difftc`) requires deciding how to compose that from HRRR's hourly cycles —
this is not yet settled (see `difftc/docs/proposal_context.md` §21) and this
script does not do that composition; it only fetches individual raw GRIB2
files per cycle/forecast-hour.

### Example

```bash
python get-hrrr.py --date 2024010100 --outdir $SCRATCH/sewx/hrrr
```

Downloads both the `wrfnat` and `wrfprs` f00 files for the 2024-01-01 00Z cycle from:

```
https://noaa-hrrr-bdp-pds.s3.amazonaws.com/hrrr.20240101/conus/hrrr.t00z.wrfnatf00.grib2
https://noaa-hrrr-bdp-pds.s3.amazonaws.com/hrrr.20240101/conus/hrrr.t00z.wrfprsf00.grib2
```

... into `$SCRATCH/sewx/hrrr/HRRR_wrfnat_f00_2024010100.grib2` and `$SCRATCH/sewx/hrrr/HRRR_wrfprs_f00_2024010100.grib2`.