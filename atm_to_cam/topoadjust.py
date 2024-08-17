import numpy as np
import xarray as xr
import os
import glob
import horizremap
import vertremap
from constants import grav, Rd, gamma_s, p0

def topo_adjustment(ps_fv, t_fv, q_fv, u_fv, v_fv, cldliq_fv, cldice_fv, hya, hyb, dycore, model_topo_file, datasource, grb_file, lev, yearstr, monthstr, daystr, cyclestr, wgt_filename, adjust_config, RDADIR="", add_cloud_vars=False):
    sePS = ps_fv
    dim_sePS = sePS.shape

    if dycore != "mpas" and model_topo_file and os.path.exists(model_topo_file):
        print(f"TOPOADJUST: Performing hydrostatic correction for surface pressure using {model_topo_file}")

        # Set flags for config
        tempadjustflag = adjust_config[0] if adjust_config else ""
        print(f"TOPOADJUST: tempadjustflag set to: {tempadjustflag}")

        # Load additional fields from reanalysis
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
            print("TOPOADJUST: no available reanalysis surface data for correction... exiting...")
            return

        # Regrid to model grid
        topo_data_SE, selat, selon = horizremap.remap_with_weights_wrapper(topo_data, wgt_filename)
        sfct_data_SE, selat, selon = horizremap.remap_with_weights_wrapper(sfct_data, wgt_filename)

        # Load model orography
        ttfile = xr.open_dataset(model_topo_file)

        if dycore == "se":
            if "PHIS" in ttfile.variables and ttfile["PHIS"].shape == dim_sePS:
                print("TOPOADJUST_SE: Using PHIS for topo adjustment")
                topo_model_SE = ttfile["PHIS"].values
            elif "PHIS_d" in ttfile.variables and ttfile["PHIS_d"].shape == dim_sePS:
                print("TOPOADJUST_SE: Using PHIS_d for topo adjustment")
                topo_model_SE = ttfile["PHIS_d"].values
            elif "PHIS_gll" in ttfile.variables and ttfile["PHIS_gll"].shape == dim_sePS:
                print("TOPOADJUST_SE: Using PHIS_gll for topo adjustment")
                topo_model_SE = ttfile["PHIS_gll"].values
            else:
                print(f"TOPOADJUST_SE: Topo error on: {model_topo_file}")
                print(f"TOPOADJUST_SE: Cannot find valid PHIS or PHIS_d or PHIS_gll with dimensions that match ncol: {dim_sePS[0]}")
                return
        else:
            print("TOPOADJUST_nonSE: Using PHIS for topo adjustment")
            topo_model_SE = ttfile["PHIS"].values

        # Unpack FV
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

        # Check number of columns corrected
        vert_corrs = 0
        tcorriter = 0
        correct_or_not = np.zeros_like(topo_data_SE)

        for kk in range(ncol):
            if not np.isnan(sfct_data_SE[kk]):  # check if missing (useful for regional data)
                if kk % 10000 == 0:
                    percentage = 100. * kk / (ncol - 1)
                    print(f"Done with {percentage:.2f} % of correction loop")

                # Get difference in geopotential, skip if difference is small
                deltaPhi = topo_data_SE[kk] - topo_model_SE[kk]
                if abs(deltaPhi) < 10.0:
                    continue

                correct_or_not[kk] = 1

                # Estimate "model's" surface temperature
                Tsfc_fv = sfct_data_SE[kk] + gamma_s * (deltaPhi / grav)

                # Calculate layer mean temperature for use in hydrostatic
                lwr_coef = 0.5  # use 0.5 for straight average
                Tlayermean = (lwr_coef * Tsfc_fv + (1. - lwr_coef) * sfct_data_SE[kk]) if Tsfc_fv >= sfct_data_SE[kk] else ((1. - lwr_coef) * Tsfc_fv + lwr_coef * sfct_data_SE[kk])

                # Correct very warm and very cold layer means
                if Tlayermean < 255.0:
                    Tlayermean = (255.0 + Tlayermean) / 2.
                    tcorriter += 1
                elif Tlayermean > 290.5:
                    Tlayermean = (290.5 + Tlayermean) / 2.
                    tcorriter += 1

                # Correct ps_fv
                ps_orig = ps_fv[kk]
                beta = np.exp(deltaPhi / Rd / Tlayermean)
                ps_fv[kk] *= beta

                # Correct T by shifting Tbot "down" the same delta
                if tempadjustflag == "a" and dycore != "fv":
                    nlev = len(lev)
                    t_fv[nlev - 1, kk] += Tsfc_fv - sfct_data_SE[kk]

                # Correct other state variables
                vert_interp_thresh = 0.1  # ps corr diff (Pa) req. to interp vert profiles
                extrap_threshold = 5000.  # maximum ps corr diff (Pa) to allow extrapolation
                if abs(ps_orig - ps_fv[kk]) > vert_interp_thresh:
                    pm_orig = hya * p0 + hyb * ps_orig
                    pm_corr = hya * p0 + hyb * ps_fv[kk]
                    linlog = 2 if abs(ps_orig - ps_fv[kk]) > extrap_threshold else -2

                    if kk == 4520:
                        print(deltaPhi)
                        print(gamma_s)
                        print(topo_data_SE[kk])
                        print(topo_model_SE[kk])
                        print(sfct_data_SE[kk])
                        print("*")
                        print(t_fv[:, kk])
                        print(pm_orig)
                    t_fv[:, kk] = np.where(np.isnan(vertremap.int2p(pm_orig, t_fv[:, kk], pm_corr, linlog)), t_fv[:, kk], vertremap.int2p(pm_orig, t_fv[:, kk], pm_corr, linlog))
                    q_fv[:, kk] = np.where(np.isnan(vertremap.int2p(pm_orig, q_fv[:, kk], pm_corr, linlog)), q_fv[:, kk], vertremap.int2p(pm_orig, q_fv[:, kk], pm_corr, linlog))
                    u_fv[:, kk] = np.where(np.isnan(vertremap.int2p(pm_orig, u_fv[:, kk], pm_corr, linlog)), u_fv[:, kk], vertremap.int2p(pm_orig, u_fv[:, kk], pm_corr, linlog))
                    v_fv[:, kk] = np.where(np.isnan(vertremap.int2p(pm_orig, v_fv[:, kk], pm_corr, linlog)), v_fv[:, kk], vertremap.int2p(pm_orig, v_fv[:, kk], pm_corr, linlog))
                    if kk == 4520:
                        print(t_fv[:, kk])
                        print(pm_corr)
                        print("-------------------------------------")

                    if add_cloud_vars:
                        cldice_fv[:, kk] = np.where(np.isnan(vertremap.int2p(pm_orig, cldice_fv[:, kk], pm_corr, linlog)), cldice_fv[:, kk], vertremap.int2p(pm_orig, cldice_fv[:, kk], pm_corr, linlog))
                        cldliq_fv[:, kk] = np.where(np.isnan(vertremap.int2p(pm_orig, cldliq_fv[:, kk], pm_corr, linlog)), cldliq_fv[:, kk], vertremap.int2p(pm_orig, cldliq_fv[:, kk], pm_corr, linlog))

                    vert_corrs += 1

        print(f"needed to correct {vert_corrs} vertical profiles for updated PS")
        print(f"needed to correct {tcorriter} temps for being too cold or too hot")

        # Cleanup
        del beta, ps_orig, tcorriter, Tlayermean, Tsfc_fv, deltaPhi
        del ttfile, topo_model_SE, topo_data_SE, sfct_data_SE, topo_data, sfct_data

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
