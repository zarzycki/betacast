#!/usr/bin/env python
"""Download a HRRR GRIB2 file from the NOAA HRRR AWS Open Data bucket.

Example:
    python get-hrrr.py --date 2024010100 --outdir /tmp/hrrr
"""

import argparse
import os
import sys

import requests

BASE_URL = "https://noaa-hrrr-bdp-pds.s3.amazonaws.com"
VALID_PRODUCTS = ("wrfnat", "wrfprs")


def parse_args():
    parser = argparse.ArgumentParser(description="Download a HRRR GRIB2 file from AWS Open Data.")
    parser.add_argument("--date", type=str, required=True,
                        help="Cycle date/time in YYYYMMDDHH format (10 characters), e.g. 2024010100")
    parser.add_argument("--outdir", type=str, required=True,
                        help="Directory to write the downloaded GRIB2 file to")
    parser.add_argument("--product", type=str, choices=VALID_PRODUCTS + ("both",), default="both",
                        help="HRRR product: wrfnat (native levels), wrfprs (pressure levels), "
                             "or both (default: both)")
    parser.add_argument("--fxx", type=int, default=0,
                        help="Forecast hour, 0-48 (default: 0)")
    parser.add_argument("--region", type=str, default="conus",
                        help="HRRR region/domain, e.g. conus or alaska (default: conus)")
    return parser.parse_args()


def validate_date(date_str):
    if len(date_str) != 10 or not date_str.isdigit():
        sys.exit(f"ERROR: --date must be YYYYMMDDHH (10 digits), got '{date_str}'")
    yyyy, mm, dd, hh = date_str[0:4], date_str[4:6], date_str[6:8], date_str[8:10]
    if not (1 <= int(mm) <= 12):
        sys.exit(f"ERROR: invalid month in --date '{date_str}'")
    if not (1 <= int(dd) <= 31):
        sys.exit(f"ERROR: invalid day in --date '{date_str}'")
    if not (0 <= int(hh) <= 23):
        sys.exit(f"ERROR: invalid hour in --date '{date_str}'")
    return yyyy, mm, dd, hh


def build_url(yyyy, mm, dd, hh, region, product, fxx):
    remote_filename = f"hrrr.t{hh}z.{product}f{fxx:02d}.grib2"
    local_filename = f"HRRR_{product}_f{fxx:02d}_{yyyy}{mm}{dd}{hh}.grib2"
    return f"{BASE_URL}/hrrr.{yyyy}{mm}{dd}/{region}/{remote_filename}", local_filename


def download_file(url, dest_path):
    with requests.get(url, stream=True, timeout=60) as resp:
        if resp.status_code != 200:
            sys.exit(f"ERROR: failed to download {url} (HTTP {resp.status_code})")
        with open(dest_path, "wb") as f:
            for chunk in resp.iter_content(chunk_size=1024 * 1024):
                f.write(chunk)


def main():
    args = parse_args()

    yyyy, mm, dd, hh = validate_date(args.date)

    if not (0 <= args.fxx <= 48):
        sys.exit(f"ERROR: --fxx must be between 0 and 48, got {args.fxx}")

    os.makedirs(args.outdir, exist_ok=True)

    products = VALID_PRODUCTS if args.product == "both" else (args.product,)

    for product in products:
        url, filename = build_url(yyyy, mm, dd, hh, args.region, product, args.fxx)
        dest_path = os.path.join(args.outdir, filename)

        print(f"Downloading {url}")
        download_file(url, dest_path)
        print(f"Saved to {dest_path}")


if __name__ == "__main__":
    main()
