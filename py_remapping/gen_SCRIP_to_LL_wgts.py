# Usage examples:
#   Global grid, default patch interpolation:
#     python gen_SCRIP_to_LL_wgts.py --gridName philly128x8pg2 --outres 0.02x0.02 --srcGridDir /path/to/scrip/ --srcGridFile Philadelphia_TC_grid_v2_ne128x8_pg2_SCRIP.nc --dstGridDir /glade/u/home/zarzycki/scratch/
#
#   Regional grid with predefined domain:
#     python gen_SCRIP_to_LL_wgts.py --gridName philly128x8pg2 --outres 0.02x0.02 --srcGridDir /path/to/scrip/ --srcGridFile Philadelphia_TC_grid_v2_ne128x8_pg2_SCRIP.nc --dstGridDir /glade/u/home/zarzycki/scratch/ --regional --regdomain midatlantic
#
#   Different interpolation method:
#     python gen_SCRIP_to_LL_wgts.py --gridName philly128x8pg2 --outres 0.02x0.02 --srcGridDir /path/to/scrip/ --srcGridFile Philadelphia_TC_grid_v2_ne128x8_pg2_SCRIP.nc --dstGridDir /glade/u/home/zarzycki/scratch/ --InterpMethod bilinear
#
#   Custom domain via explicit corners:
#     python gen_SCRIP_to_LL_wgts.py --gridName philly128x8pg2 --outres 0.02x0.02 --srcGridDir /path/to/scrip/ --srcGridFile Philadelphia_TC_grid_v2_ne128x8_pg2_SCRIP.nc --dstGridDir /glade/u/home/zarzycki/scratch/ --regional --regdomain tclf --reg_minLLbnds 12.0 215.0 --reg_maxLLbnds 55.0 299.0
#
#   Custom domain via center + total widths:
#     python gen_SCRIP_to_LL_wgts.py --gridName philly128x8pg2 --outres 0.02x0.02 --srcGridDir /path/to/scrip/ --srcGridFile Philadelphia_TC_grid_v2_ne128x8_pg2_SCRIP.nc --dstGridDir /glade/u/home/zarzycki/scratch/ --regional --regdomain tclf --regcenlat 20.0 --regcenlon 25.0 --reglatwidth 20.0 --reglonwidth 24.0
#
#   A practical example taking a VR/RRM mesh to a high-res regional/urban RLL mesh
#   STRING="003"
#   python gen_SCRIP_to_LL_wgts.py \
#     --gridName "zclf${STRING}-mp3a" \
#     --outres 0.025x0.025 \
#     --srcGridDir "/glade/work/zarzycki/grids/scrip/" \
#     --srcGridFile "mpasa3-60-tclf${STRING}_scrip.nc" \
#     --dstGridDir "/glade/derecho/scratch/zarzycki/" \
#     --regional \
#     --regcenlat 28.4 \
#     --regcenlon -92.9 \
#     --reglatwidth 20.0 \
#     --reglonwidth 24.0

import os
import shutil
import logging
import argparse
import subprocess
from ESMF_regridding import latlon_to_SCRIP, esmf_regrid_gen_weights

# Set up logging
logging.basicConfig(level=logging.INFO)

# Arguments
parser = argparse.ArgumentParser(description="Generate SCRIP to lat-lon regridding weights")
parser.add_argument("--gridName", required=True, help="Name of the source grid (e.g. philly128x8pg2)")
parser.add_argument("--outres", required=True, help="Output resolution (e.g. 0.02x0.02)")
parser.add_argument("--srcGridDir", required=True, help="Directory containing the source SCRIP file")
parser.add_argument("--srcGridFile", required=True, help="Source SCRIP filename")
parser.add_argument("--dstGridDir", required=True, help="Directory for destination SCRIP and weight files")
parser.add_argument("--InterpMethod", default="patch",
                    choices=["bilinear", "patch", "conserve", "nearestdtos", "neareststod"],
                    help="Interpolation method (default: patch)")
parser.add_argument("--regional", action="store_true", default=False,
                    help="Generate a regional destination grid (default: False)")
parser.add_argument("--regdomain", help="Predefined regional domain name (e.g. midatlantic, asd, atlantic, tctest, florida, narr). For custom domains, omit or use any other name and provide bounds via --reg_minLLbnds/--reg_maxLLbnds or --regcenlat/--regcenlon/--reglatwidth/--reglonwidth")
parser.add_argument("--reg_minLLbnds", type=float, nargs=2, metavar=("LAT", "LON"),
                    help="Lower-left corner [lat lon] for custom regional domain")
parser.add_argument("--reg_maxLLbnds", type=float, nargs=2, metavar=("LAT", "LON"),
                    help="Upper-right corner [lat lon] for custom regional domain")
parser.add_argument("--regcenlat", type=float, help="Center latitude for custom regional domain")
parser.add_argument("--regcenlon", type=float, help="Center longitude for custom regional domain")
parser.add_argument("--reglatwidth", type=float, help="Total width in latitude (degrees) for custom regional domain")
parser.add_argument("--reglonwidth", type=float, help="Total width in longitude (degrees) for custom regional domain")
args = parser.parse_args()

# Setup
gridName = args.gridName
InterpMethod = args.InterpMethod
regional = args.regional
regdomain = args.regdomain if args.regdomain is not None else ""
outres = args.outres
srcGridDir = args.srcGridDir
srcGridFile = args.srcGridFile
dstGridDir = args.dstGridDir
# Auto-enable regional if any of the "reg____" flags were passed without --regional
reg_flags = [args.regdomain, args.reg_minLLbnds, args.reg_maxLLbnds,
             args.regcenlat, args.regcenlon, args.reglatwidth, args.reglonwidth]
if any(v is not None for v in reg_flags) and not regional:
    regional = True
    logging.warning("Regional flags were provided but --regional was not set; enabling regional mode. "
                    "Add --regional to suppress this warning.")

# Generate target RLL mesh name and final weight name
if regional:
    dstGridFile = f"{outres}_reg_SCRIP.nc"
    wgtFileName = f"map_{gridName}_to_{outres}reg_{InterpMethod}.nc"
else:
    dstGridFile = f"{outres}_SCRIP.nc"
    wgtFileName = f"map_{gridName}_to_{outres}glob_{InterpMethod}.nc"
wgtFileDir = dstGridDir

# Construct full paths
srcGridName = os.path.join(srcGridDir, srcGridFile)
dstGridName = os.path.join(dstGridDir, dstGridFile)

# Set options for generating the SCRIP file
Opt = {
    "ForceOverwrite": True,
    "PrintTimings": True,
    "Debug": True
}

# Set regional domain corners
if regional:
    rll_option1_args = [args.reg_minLLbnds, args.reg_maxLLbnds]
    rll_option2_args = [args.regcenlat, args.regcenlon, args.reglatwidth, args.reglonwidth]

    have_rll_opt1 = all(v is not None for v in rll_option1_args)
    have_partial_rll_opt1 = any(v is not None for v in rll_option1_args) and not have_rll_opt1

    have_rll_opt2 = all(v is not None for v in rll_option2_args)
    have_partial_rll_opt2 = any(v is not None for v in rll_option2_args) and not have_rll_opt2

    if have_partial_rll_opt1:
        parser.error("Must provide both --reg_minLLbnds and --reg_maxLLbnds together")
    if have_partial_rll_opt2:
        parser.error("Must provide all four of --regcenlat, --regcenlon, --reglatwidth, --reglonwidth together")

    if have_rll_opt2:
        if args.regdomain:
            logging.warning(f"--regdomain '{args.regdomain}' ignored because --regcenlat/--regcenlon/--reglatwidth/--reglonwidth were provided")
        cenlon = args.regcenlon + 360.0 if args.regcenlon < 0 else args.regcenlon
        Opt["LLCorner"] = [args.regcenlat - args.reglatwidth / 2.0, cenlon - args.reglonwidth / 2.0]
        Opt["URCorner"] = [args.regcenlat + args.reglatwidth / 2.0, cenlon + args.reglonwidth / 2.0]
    elif have_rll_opt1:
        if args.regdomain:
            logging.warning(f"--regdomain '{args.regdomain}' ignored because --reg_minLLbnds/--reg_maxLLbnds were provided")
        Opt["LLCorner"] = args.reg_minLLbnds
        Opt["URCorner"] = args.reg_maxLLbnds
    # Supported regional mesh domains
    elif regdomain == "asd":
        Opt["LLCorner"] = [10.0, 230.0]
        Opt["URCorner"] = [55.0, 299.0]
    elif regdomain == "atlantic":
        Opt["LLCorner"] = [5.0, 250.0]
        Opt["URCorner"] = [55.0, 355.0]
    elif regdomain == "tctest":
        Opt["LLCorner"] = [5.0, 250.0]
        Opt["URCorner"] = [55.0, 355.0]
    elif regdomain == "florida":
        Opt["LLCorner"] = [15.0, 271.0]
        Opt["URCorner"] = [35.0, 296.0]
    elif regdomain == "midatlantic":
        Opt["LLCorner"] = [31.0, 278.0]
        Opt["URCorner"] = [47.0, 293.0]
    elif regdomain == "narr":
        Opt["LLCorner"] = [12.0, 215.0]
        Opt["URCorner"] = [65.0, 298.0]
    else:
        parser.error("Regional grid requested but no supported domain specified. Use --regdomain, "
                     "--reg_minLLbnds/--reg_maxLLbnds, or --regcenlat/--regcenlon/--reglatwidth/--reglonwidth")

# Generate the regular SCRIP file
logging.info("Generating SCRIP file...")
latlon_to_SCRIP(dstGridName, outres, Opt)

# Reset options for generating weights
Opt = {
    "InterpMethod": InterpMethod,
    "ForceOverwrite": True,
    "PrintTimings": True,
    "SrcRegional": True
}

# Generate the weights file
logging.info("Generating weights file...")
# Note, if we were to use ncremap or TE, etc. to generate weights we'd add function here
esmf_regrid_gen_weights(srcGridName, dstGridName, os.path.join(wgtFileDir, wgtFileName), Opt)

# Compress the weight file using ncks if available
# This must be LOSSLESS compression using DEFLATE!! No bit grooming, etc.
ncks_exists = shutil.which("ncks") is not None
if ncks_exists:
    logging.info(f"NCO compressing: {os.path.join(wgtFileDir, wgtFileName)}")
    before_size = os.path.getsize(os.path.join(wgtFileDir, wgtFileName))
    subprocess.run(["ncks", "-O", "-4", "-L", "1", os.path.join(wgtFileDir, wgtFileName), os.path.join(wgtFileDir, wgtFileName)])
    after_size = os.path.getsize(os.path.join(wgtFileDir, wgtFileName))
    compression_percentage = 100.0 * (1 - after_size / before_size)
    logging.info(f"Compression percentage: {compression_percentage:.2f}%")
else:
    logging.warning(f"Cannot compress {os.path.join(wgtFileDir, wgtFileName)}, ncks not found")

pet_log = "PET0.RegridWeightGen.Log"
if os.path.exists(pet_log):
    os.remove(pet_log)
    logging.info(f"Removed {pet_log}")

logging.info("Completed SCRIP and weight file generation.")
