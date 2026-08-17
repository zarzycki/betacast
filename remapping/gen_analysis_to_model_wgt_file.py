import os
import sys
import subprocess
import logging
import argparse
import shutil
import xarray as xr

from ESMF_regridding import esmf_regrid_gen_weights  # Import the function

# Betacast modules
module_path = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', 'py_functions'))
if module_path not in sys.path:
    sys.path.append(module_path)
import pyfuncs

logging.basicConfig(
    level=logging.INFO,  # Change to logging.DEBUG for more detailed logs
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
    handlers=[
        logging.StreamHandler()  # Ensure logs are sent to the console
    ]
)

# Is the Betacast path available to us?
BETACAST, PATHTOHERE = pyfuncs.get_betacast_path()
ANLSCRIPPATH = os.path.join(BETACAST, 'grids/anl_scrip/')

# Argument parser
parser = argparse.ArgumentParser(
    description=(
        "Generate ESMF patch-interpolation weight files that map an analysis/reanalysis grid "
        "to a model destination grid (or vice versa with --FLIP_MODEL_AND_ANALYSIS). "
        "The output weight file is written to WGTFILEDIR and optionally compressed with ncks.\n\n"
        "Supported analysis grids (--ANLGRID):\n"
        "  Global grids  : era5_0.25x0.25, era5_0.3gaus, era5_2deg,\n"
        "                  gfs_0.25x0.25, gfs_0.50x0.50\n"
        "  Regional grids: rap_13km, hrrr_3km, hwrf_storm\n\n"
        "Supported destination grid file formats (--DSTGRIDFILE):\n"
        "  SCRIP (.nc)  : variables grid_center_lat / grid_corner_lat\n"
        "  ESMF mesh    : variables nodeCoords / elementConn\n"
        "  (Exodus files are NOT supported; convert to a physics-grid SCRIP first)\n\n"
        "Example model SCRIP files are in betacast/grids/model_scrip/, e.g.:\n"
        "  ne30pg3_e3sm_scrip.nc, ne30pg2_e3sm_scrip.nc, ne120pg3_e3sm_scrip.nc,\n"
        "  ne30np4_091226_pentagons.nc, 0.9x1.25_c110307.nc, 1.9x2.5_c110308.nc"
    ),
    epilog=(
        "Examples:\n"
        "  # ERA5 0.25-deg → ne30pg3 (analysis-to-model, typical use case)\n"
        "  python gen_analysis_to_model_wgt_file.py \\\n"
        "      --ANLGRID era5_0.25x0.25 \\\n"
        "      --DSTGRIDNAME ne30pg3_e3sm \\\n"
        "      --DSTGRIDFILE /path/to/betacast/grids/model_scrip/ne30pg3_e3sm_scrip.nc \\\n"
        "      --WGTFILEDIR /scratch/weights/\n\n"
        "  # GFS 0.25-deg → FV 1-deg (analysis-to-model)\n"
        "  python gen_analysis_to_model_wgt_file.py \\\n"
        "      --ANLGRID gfs_0.25x0.25 \\\n"
        "      --DSTGRIDNAME fv_0.9x1.25 \\\n"
        "      --DSTGRIDFILE /path/to/betacast/grids/model_scrip/0.9x1.25_c110307.nc \\\n"
        "      --WGTFILEDIR ./\n\n"
        "  # ne30pg3 → ERA5 0.25-deg (model-to-analysis, for e.g. bias correction)\n"
        "  python gen_analysis_to_model_wgt_file.py \\\n"
        "      --ANLGRID era5_0.25x0.25 \\\n"
        "      --DSTGRIDNAME ne30pg3_e3sm \\\n"
        "      --DSTGRIDFILE /path/to/betacast/grids/model_scrip/ne30pg3_e3sm_scrip.nc \\\n"
        "      --WGTFILEDIR /scratch/weights/ \\\n"
        "      --FLIP_MODEL_AND_ANALYSIS"
    ),
    formatter_class=argparse.RawDescriptionHelpFormatter,
)
parser.add_argument(
    "--ANLGRID", type=str, default="era5_0.25x0.25",
    help=(
        "Analysis/reanalysis source grid name. "
        "Valid options: era5_0.25x0.25, era5_0.3gaus, era5_2deg, "
        "gfs_0.25x0.25, gfs_0.50x0.50, rap_13km, hrrr_3km, hwrf_storm. "
        "Regional grids (rap_13km, hrrr_3km, hwrf_storm) automatically set SrcRegional=True. "
        "Default: era5_0.25x0.25"
    ),
)
parser.add_argument(
    "--ANLGRIDPATH", type=str, default=ANLSCRIPPATH,
    help=(
        "Directory containing the analysis SCRIP files named <ANLGRID>_scrip.nc. "
        f"Default: {ANLSCRIPPATH}"
    ),
)
parser.add_argument(
    "--DSTGRIDNAME", type=str, default="Philadelphia_TC_grid_v2_ne128x8_pg2",
    help=(
        "Short label for the destination model grid, used in the output weight filename: "
        "map_<ANLGRID>_TO_<DSTGRIDNAME>_patc.nc. "
        "Default: Philadelphia_TC_grid_v2_ne128x8_pg2"
    ),
)
parser.add_argument(
    "--DSTGRIDFILE", type=str, required=True,
    help=(
        "Full path to the model destination grid file. "
        "Accepted formats: SCRIP (grid_center_lat/grid_corner_lat) or ESMF mesh (nodeCoords/elementConn). "
        "Pre-built SCRIP files for common grids are in betacast/grids/model_scrip/."
    ),
)
parser.add_argument(
    "--WGTFILEDIR", type=str, default="./",
    help="Directory where the output weight file will be written. Default: current directory",
)
parser.add_argument(
    "--FLIP_MODEL_AND_ANALYSIS", action="store_true",
    help=(
        "Reverse the mapping direction: model → analysis instead of analysis → model. "
        "Output filename becomes map_<DSTGRIDNAME>_TO_<ANLGRID>_patc.nc. "
        "Default: False (analysis → model)"
    ),
)

args = parser.parse_args()

# User settings from command line arguments
anlgrid = args.ANLGRID
anlgridpath = args.ANLGRIDPATH
dstGridName = args.DSTGRIDNAME
dstGridFile = args.DSTGRIDFILE
wgtFileDir = args.WGTFILEDIR
flip_model_and_analysis = args.FLIP_MODEL_AND_ANALYSIS

# Validate anlgrid
valid_anlgrids = ["era5_0.25x0.25", "era5_0.3gaus", "era5_2deg", "gfs_0.25x0.25", "gfs_0.50x0.50", "rap_13km", "hrrr_3km", "hwrf_storm"]
if anlgrid not in valid_anlgrids:
    logging.error(f"Unsupported analysis grid: {anlgrid}")
    exit()

srcGridFile = os.path.join(anlgridpath, f"{anlgrid}_scrip.nc")

# Check if the grid is regional
regionalGrid = anlgrid in ["rap_13km", "hrrr_3km", "hwrf_storm"]

# Determine destination grid type
f = xr.open_dataset(dstGridFile)
if "nodeCoords" in f.variables and "elementConn" in f.variables:
    dstType = "esmf"
elif "grid_center_lat" in f.variables and "grid_corner_lat" in f.variables:
    dstType = "scrip"
elif "eb_prop1" in f.variables:
    logging.error("Exodus files are not supported. Generate the physics grid and pass that in.")
    exit()
else:
    dstType = "model"
logging.info(f"Determined input grid type to be: {dstType}")

# Handle "model" type
if dstType == "model":
    exit()

InterpMethod = "patch"
shortInterpName = "patc" if InterpMethod == "patch" else InterpMethod

# Generate weights file
logging.info("Generating weights!")
if not flip_model_and_analysis:
    wgtFileName = f"map_{anlgrid}_TO_{dstGridName}_{shortInterpName}.nc"
    srcGridName = srcGridFile
    dstGridName = dstGridFile
else:
    wgtFileName = f"map_{dstGridName}_TO_{anlgrid}_{shortInterpName}.nc"
    srcGridName = dstGridFile
    dstGridName = srcGridFile

Opt = {
    "InterpMethod": InterpMethod,
    "ForceOverwrite": True,
    "PrintTimings": True,
    "NoPETLog": True,
    "RemovePETLog": True,
    "Debug": True,
    "NetCDFType": "netcdf4"
}

# Adjust options based on the grid type and whether it's a regional grid
if regionalGrid:
    if not flip_model_and_analysis:
        Opt["SrcRegional"] = True
    else:
        Opt["DstRegional"] = True

if dstType in ["model", "esmf"]:
    if not flip_model_and_analysis:
        Opt["DstESMF"] = True
        if dstType == "model":
            dstGridName = seGridName
    else:
        Opt["SrcESMF"] = True
        if dstType == "model":
            srcGridName = seGridName

# Run the regrid weight generation
success = esmf_regrid_gen_weights(srcGridName, dstGridName, os.path.join(wgtFileDir, wgtFileName), Opt)
if success is False:
    logging.error(f"Weight generation failed. Aborting.")
    sys.exit(1)

# Cleanup
if dstType == "model":
    logging.info("Removing online generated grid descriptor")
    os.remove(seGridName)

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

logging.info(f"Successfully generated: {os.path.join(wgtFileDir, wgtFileName)}")
logging.info("Done")