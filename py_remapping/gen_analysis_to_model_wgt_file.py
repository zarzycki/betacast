import os
import subprocess
import logging
from ESMF_regridding import esmf_regrid_gen_weights  # Import the function

# User settings
anlgrid = "era5_0.25x0.25"
anlgridpath = "./anl_scrip/"
dstGridName = "Philadelphia_TC_grid_v2_ne128x8_pg2"
dstGridFile = "/global/homes/c/czarzyck/m2637/E3SM_SCREAM_files/grids/scrip/Philadelphia_TC_grid_v2_ne128x8_pg2_SCRIP.nc"
wgtFileDir = "./"
flip_model_and_analysis = False

# Check if things are provided through environment variables (similar to command-line args)
anlgrid = os.getenv('ANLGRID', anlgrid)
anlgridpath = os.getenv('ANLGRIDPATH', anlgridpath)
dstGridName = os.getenv('DSTGRIDNAME', dstGridName)
dstGridFile = os.getenv('DSTGRIDFILE', dstGridFile)
wgtFileDir = os.getenv('WGTFILEDIR', wgtFileDir)
flip_model_and_analysis = os.getenv('FLIP_MODEL_AND_ANALYSIS', flip_model_and_analysis) == 'True'

# Validate anlgrid
valid_anlgrids = ["era5_0.25x0.25", "gfs_0.25x0.25", "gfs_0.50x0.50", "rap_13km", "hrrr_3km", "hwrf_storm"]
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
    logging.info("Generating grid file for model!")
    if "lat" not in f.variables or "lon" not in f.variables:
        logging.error("Destination file does not have lat/lon coords, exiting!")
        exit()

    lat = f["lat"].values
    lon = f["lon"].values

    # Generate SE grid
    seGridName = "grid_se.nc"
    Opt_se = {
        "ForceOverwrite": True,
        "PrintTimings": True,
        "Title": "SE Grid"
    }
    unstructured_to_ESMF(seGridName, lat, lon, Opt_se)
    f.close()

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
esmf_regrid_gen_weights(srcGridName, dstGridName, os.path.join(wgtFileDir, wgtFileName), Opt)

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

