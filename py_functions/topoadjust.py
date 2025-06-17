import numpy as np
import xarray as xr
import os
import vertremap
from constants import grav, Rd, gamma_s, p0, vert_interp_thresh, extrap_threshold
import pyfuncs
import packing
import logging
logger = logging.getLogger(__name__)


def load_model_orography(dycore, model_topo_file, dim_sePS):
    """Load model orography based on dycore type."""
    ttfile = xr.open_dataset(model_topo_file)
    if dycore == "se" or dycore == "scream":
        if "PHIS" in ttfile.variables and ttfile["PHIS"].shape == dim_sePS:
            return ttfile["PHIS"].values
        elif "PHIS_d" in ttfile.variables and ttfile["PHIS_d"].shape == dim_sePS:
            return ttfile["PHIS_d"].values
        elif "PHIS_gll" in ttfile.variables and ttfile["PHIS_gll"].shape == dim_sePS:
            return ttfile["PHIS_gll"].values
        else:
            raise ValueError(f"Cannot find valid PHIS or PHIS_d or PHIS_gll with dimensions that match ncol: {dim_sePS[0]}")
    else:
        return ttfile["PHIS"].values


def correct_state_variables(data_horiz, correct_or_not, tempadjustflag, dycore):

    vert_corrs = 0
    tcorriter = 0

    pyfuncs.print_min_max_dict(data_horiz)

    dim_sePS = data_horiz['ps'].shape
    ncol = dim_sePS[0]

    for kk in range(ncol):
        if not np.isnan(data_horiz['ts'][kk]):
            deltaPhi = data_horiz['phis'][kk] - data_horiz['phis_SE'][kk]
            if abs(deltaPhi) < 10.0:
                continue

            Tsfc_fv = data_horiz['ts'][kk] + gamma_s * (deltaPhi / grav)
            lwr_coef = 0.5
            Tlayermean = (lwr_coef * Tsfc_fv + (1. - lwr_coef) * data_horiz['ts'][kk]) if Tsfc_fv >= data_horiz['ts'][kk] else ((1. - lwr_coef) * Tsfc_fv + lwr_coef * data_horiz['ts'][kk])

            if Tlayermean < 255.0:
                Tlayermean = (255.0 + Tlayermean) / 2.
                tcorriter += 1
            elif Tlayermean > 290.5:
                Tlayermean = (290.5 + Tlayermean) / 2.
                tcorriter += 1

            ps_orig = data_horiz['ps'][kk]
            beta = np.exp(deltaPhi / Rd / Tlayermean)
            data_horiz['ps'][kk] *= beta

            correct_or_not[kk] = 1

            if tempadjustflag == "a" and dycore != "fv":
                nlev = len(data_horiz['hya'])
                data_horiz['t'][nlev - 1, kk] += Tsfc_fv - data_horiz['ts'][kk]

            if abs(ps_orig - data_horiz['ps'][kk]) > vert_interp_thresh:
                pm_orig = data_horiz['hya'] * p0 + data_horiz['hyb'] * ps_orig
                pm_corr = data_horiz['hya'] * p0 + data_horiz['hyb'] * data_horiz['ps'][kk]
                linlog = 2 if abs(ps_orig - data_horiz['ps'][kk]) > extrap_threshold else -2

                allowable_interp_vars = vertremap._VERT_REMAP_VARS

                for var in data_horiz:
                    if var in allowable_interp_vars:
                        interp_result = vertremap.int2p(pm_orig, data_horiz[var][:, kk], pm_corr, linlog)
                        data_horiz[var][:, kk] = np.where(np.isnan(interp_result), data_horiz[var][:, kk], interp_result)

                vert_corrs += 1

    return vert_corrs, tcorriter, correct_or_not


def topo_adjustment(data_horiz, dycore, model_topo_file, adjust_config):

    dim_sePS = data_horiz['ps'].shape

    correct_or_not = np.zeros_like(data_horiz['ps']).astype(np.float32)

    if dycore != "mpas" and model_topo_file and os.path.exists(model_topo_file):
        logging.info(f"TOPOADJUST: Performing hydrostatic correction for surface pressure using {model_topo_file}")

        tempadjustflag = adjust_config[0] if adjust_config else ""
        logging.info(f"TOPOADJUST: tempadjustflag set to: {tempadjustflag}")

        data_horiz['phis_SE'] = load_model_orography(dycore, model_topo_file, dim_sePS)

        if dycore == "fv":
            logging.info("TOPOADJUST_FV: unpacking fv vars")
            data_horiz['phis_SE'] = packing.latlon_to_ncol(data_horiz['phis_SE'])

        vert_corrs, tcorriter, correct_or_not = correct_state_variables(data_horiz, correct_or_not, tempadjustflag, dycore)

        logging.info(f"needed to correct {vert_corrs} vertical profiles for updated PS")
        logging.info(f"needed to correct {tcorriter} temps for being too cold or too hot")

    elif model_topo_file in [" ", "NULL", ""]:
        logging.info("Empty model topo file entered, not performing hydro adjustment")
        logging.info("continuing...")
    elif model_topo_file and not os.path.exists(model_topo_file):
        logging.info("model_topo_file passed in but cannot find file on Unix system")
        logging.info("if you do not want adjustment, specify NULL in the namelist")
        logging.info("exiting...")
        return
    else:
        logging.info("No model topo file passed into script, not performing hydro adjustment")
        logging.info("continuing...")

    return correct_or_not
