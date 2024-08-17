import numpy as np
import xarray as xr
import os
import glob
import horizremap
import vertremap
from constants import grav, Rd, gamma_s, p0, vert_interp_thresh, extrap_threshold
from numba import jit


def load_additional_fields(datasource, grb_file, RDADIR, yearstr, monthstr, daystr, cyclestr):
    """Load additional fields from reanalysis data based on datasource."""
    if datasource in ["GFS", "CFSR", "HWRF"]:
        topo_data = grb_file["HGT_P0_L1_GLL0"].values * grav
        sfct_data = grb_file["TMP_P0_L103_GLL0"].values if datasource == "HWRF" else grb_file["TMP_P0_L104_GLL0"].values
    elif datasource == "RAP":
        topo_data = grb_file["HGT_P0_L1_GLC0"].values * grav
        sfct_data = grb_file["TMP_P0_L1_GLC0"].values
    elif datasource in ["ERAI", "ERA5"]:
        topo_data = grb_file["phis"].isel(time=0).values
        sfct_data = grb_file["t2m"].isel(time=0).values
    elif datasource == "ERA5RDA":
        sf_dir = f"{RDADIR}/e5.oper.an.sfc/{yearstr}{monthstr}"
        topo_data = grb_file["Z"].isel(time=0).values
        rda_find_pattern = f"{sf_dir}/e5.oper.an.sfc.128_167_2t.ll025sc.{yearstr}{monthstr}0100_*.nc"
        rda_files = glob.glob(rda_find_pattern)
        if not rda_files:
            raise FileNotFoundError(f"No files found matching pattern: {rda_find_pattern}")
        rda_file = xr.open_mfdataset(rda_files, combine='by_coords')
        rda_time = rda_file["time"]
        rda_thistime = np.datetime64(f"{yearstr}-{monthstr}-{daystr}T{cyclestr[:2]}:00")
        sfct_data = rda_file["VAR_2T"].sel(time=rda_thistime, method="nearest").values
    elif datasource == "CAM":
        topo_data = grb_file["phis_rll"].values
        sfct_data = grb_file["tbot_mod"].values
    else:
        raise ValueError("No available reanalysis surface data for correction")
    return topo_data, sfct_data


def load_model_orography(dycore, model_topo_file, dim_sePS):
    """Load model orography based on dycore type."""
    ttfile = xr.open_dataset(model_topo_file)
    if dycore == "se":
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




@jit(nopython=True)
def correct_state_variables(ncol, t_fv, q_fv, u_fv, v_fv, cldliq_fv, cldice_fv, ps_fv,
                            topo_data_SE, topo_model_SE, sfct_data_SE, hya, hyb, lev,
                            tempadjustflag, dycore, add_cloud_vars):
    vert_corrs = 0
    tcorriter = 0
    correct_or_not = np.zeros_like(ps_fv)

    for kk in range(ncol):
        if not np.isnan(sfct_data_SE[kk]):
            deltaPhi = topo_data_SE[kk] - topo_model_SE[kk]
            if abs(deltaPhi) < 10.0:
                continue

            Tsfc_fv = sfct_data_SE[kk] + gamma_s * (deltaPhi / grav)
            lwr_coef = 0.5
            Tlayermean = (lwr_coef * Tsfc_fv + (1. - lwr_coef) * sfct_data_SE[kk]) if Tsfc_fv >= sfct_data_SE[kk] else ((1. - lwr_coef) * Tsfc_fv + lwr_coef * sfct_data_SE[kk])

            if Tlayermean < 255.0:
                Tlayermean = (255.0 + Tlayermean) / 2.
                tcorriter += 1
            elif Tlayermean > 290.5:
                Tlayermean = (290.5 + Tlayermean) / 2.
                tcorriter += 1

            ps_orig = ps_fv[kk]
            beta = np.exp(deltaPhi / Rd / Tlayermean)
            ps_fv[kk] *= beta

            correct_or_not[kk] = 1

            if tempadjustflag == "a" and dycore != "fv":
                nlev = len(lev)
                t_fv[nlev - 1, kk] += Tsfc_fv - sfct_data_SE[kk]

            if abs(ps_orig - ps_fv[kk]) > vert_interp_thresh:
                pm_orig = hya * p0 + hyb * ps_orig
                pm_corr = hya * p0 + hyb * ps_fv[kk]
                linlog = 2 if abs(ps_orig - ps_fv[kk]) > extrap_threshold else -2

                t_fv[:, kk] = np.where(np.isnan(vertremap.int2p(pm_orig, t_fv[:, kk], pm_corr, linlog)), t_fv[:, kk], vertremap.int2p(pm_orig, t_fv[:, kk], pm_corr, linlog))
                q_fv[:, kk] = np.where(np.isnan(vertremap.int2p(pm_orig, q_fv[:, kk], pm_corr, linlog)), q_fv[:, kk], vertremap.int2p(pm_orig, q_fv[:, kk], pm_corr, linlog))
                u_fv[:, kk] = np.where(np.isnan(vertremap.int2p(pm_orig, u_fv[:, kk], pm_corr, linlog)), u_fv[:, kk], vertremap.int2p(pm_orig, u_fv[:, kk], pm_corr, linlog))
                v_fv[:, kk] = np.where(np.isnan(vertremap.int2p(pm_orig, v_fv[:, kk], pm_corr, linlog)), v_fv[:, kk], vertremap.int2p(pm_orig, v_fv[:, kk], pm_corr, linlog))

                if add_cloud_vars:
                    cldice_fv[:, kk] = np.where(np.isnan(vertremap.int2p(pm_orig, cldice_fv[:, kk], pm_corr, linlog)), cldice_fv[:, kk], vertremap.int2p(pm_orig, cldice_fv[:, kk], pm_corr, linlog))
                    cldliq_fv[:, kk] = np.where(np.isnan(vertremap.int2p(pm_orig, cldliq_fv[:, kk], pm_corr, linlog)), cldliq_fv[:, kk], vertremap.int2p(pm_orig, cldliq_fv[:, kk], pm_corr, linlog))

                vert_corrs += 1

    return vert_corrs, tcorriter, correct_or_not


def topo_adjustment(ps_fv, t_fv, q_fv, u_fv, v_fv, cldliq_fv, cldice_fv, hya, hyb, dycore, model_topo_file, datasource, grb_file, lev, yearstr, monthstr, daystr, cyclestr, wgt_filename, adjust_config, RDADIR="", add_cloud_vars=False):
    dim_sePS = ps_fv.shape

    if dycore != "mpas" and model_topo_file and os.path.exists(model_topo_file):
        print(f"TOPOADJUST: Performing hydrostatic correction for surface pressure using {model_topo_file}")

        tempadjustflag = adjust_config[0] if adjust_config else ""
        print(f"TOPOADJUST: tempadjustflag set to: {tempadjustflag}")

        topo_data, sfct_data = load_additional_fields(datasource, grb_file, RDADIR, yearstr, monthstr, daystr, cyclestr)

        topo_data_SE, selat, selon = horizremap.remap_with_weights_wrapper(topo_data, wgt_filename)
        sfct_data_SE, selat, selon = horizremap.remap_with_weights_wrapper(sfct_data, wgt_filename)

        topo_model_SE = load_model_orography(dycore, model_topo_file, dim_sePS)

        if dycore == "fv":
            print("TOPOADJUST_FV: unpacking fv vars")
            topo_data_SE = topo_data_SE.flatten(order='F')
            sfct_data_SE = sfct_data_SE.flatten(order='F')
            topo_model_SE = topo_model_SE.flatten(order='F')
            ps_fv = latlon_to_ncol(ps_fv)
            t_fv = latlon_to_ncol(t_fv)
            q_fv = latlon_to_ncol(q_fv)
            u_fv = latlon_to_ncol(u_fv)
            v_fv = latlon_to_ncol(v_fv)
            cldice_fv = latlon_to_ncol(cldice_fv)
            cldliq_fv = latlon_to_ncol(cldliq_fv)
            ncol = dim_sePS[0] * dim_sePS[1]
        else:
            ncol = dim_sePS[0]

        vert_corrs, tcorriter, correct_or_not = correct_state_variables(
            ncol, t_fv, q_fv, u_fv, v_fv, cldliq_fv, cldice_fv, ps_fv,
            topo_data_SE, topo_model_SE, sfct_data_SE, hya, hyb, lev,
            tempadjustflag, dycore, add_cloud_vars)

        print(f"needed to correct {vert_corrs} vertical profiles for updated PS")
        print(f"needed to correct {tcorriter} temps for being too cold or too hot")

    elif model_topo_file in [" ", "NULL", ""]:
        print("Empty model topo file entered, not performing hydro adjustment")
        print("continuing...")
    elif model_topo_file and not os.path.exists(model_topo_file):
        print("model_topo_file passed in but cannot find file on Unix system")
        print("if you do not want adjustment, specify NULL in the namelist")
        print("exiting...")
        return
    else:
        print("No model topo file passed into script, not performing hydro adjustment")
        print("continuing...")

    return correct_or_not
