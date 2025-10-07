import os
import shutil
import logging
import subprocess
from ESMF_regridding import latlon_to_SCRIP, esmf_regrid_gen_weights

# Set up logging
logging.basicConfig(level=logging.INFO)

# User settings
gridName = "philly128x8pg2"
InterpMethod = "patch"  # bilinear, patch, conserve, nearestdtos, neareststod
regional = True
regdomain = "midatlantic"  # asd, atlantic, tctest, florida

# RLL GRID
outres = "0.02x0.02"  # Change this to desired resolution
srcGridDir = "/glade/u/home/zarzycki/work/grids/scrip/"
srcGridFile = "Philadelphia_TC_grid_v2_ne128x8_pg2_SCRIP.nc"

# Destination grid settings
dstGridDir = "/glade/u/home/zarzycki/scratch/"
if regional:
    dstGridFile = f"{outres}_reg_SCRIP.nc"
    wgtFileName = f"map_{gridName}_to_{outres}reg_{InterpMethod}.nc"
else:
    dstGridFile = f"{outres}_SCRIP.nc"
    wgtFileName = f"map_{gridName}_to_{outres}glob_{InterpMethod}.nc"
wgtFileDir = "/glade/u/home/zarzycki/scratch/"

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
    if regdomain == "asd":
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
esmf_regrid_gen_weights(srcGridName, dstGridName, os.path.join(wgtFileDir, wgtFileName), Opt)

# Compress the weight file using ncks if available
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

logging.info("Completed SCRIP and weight file generation.")
