# Running a nudged simulation

A "nudged" CAM simulation is just a free-running CAM simulation with additional namelist options. These options tell CAM to find analysis files that describe the target field (e.g., ERA5 reanalysis) and will "push" the model towards that state at each physics timestep.

Here, I assume nudging files are already generated. See below if you need to generate others.

Create a case as you usually would.

```
CESMROOT=/glade/work/zarzycki/cam_20230623/
CASEDIR=/glade/u/home/zarzycki/CPT/FHIST-ne30pg3-ndg-ERA5-Q24-N23-x101

$CESMROOT/cime/scripts/create_newcase --case "${CASEDIR}" --compset FHIST --res ne30pg3_ne30pg3_mg17 --mach cheyenne --pecount 432 --project P93300642 --run-unsupported

cd "${CASEDIR}"

## Specify start date
./xmlchange RUN_STARTDATE=2018-01-01

## Specify run settings
./xmlchange STOP_N=12
./xmlchange STOP_OPTION=nmonths
./xmlchange RESUBMIT=4

## Update SST file
./xmlchange SSTICE_YEAR_END=2022
./xmlchange SSTICE_DATA_FILENAME='$DIN_LOC_ROOT/atm/cam/sst/sst_HadOIBl_bc_1x1_1850_2022_c230628.nc'
```

To run over the 2018-2022 time period (initialized on 2018-01-01) use these user\_nl\_cam and user\_nl\_clm files. The output streams can be modified as desired.

#### user\_nl\_cam

```
ncdata='/glade/p/univ/upsu0032/ncdata/ERA5.ne30np4.L58.cam2.i.2018-01-01-00000.nc'

nhtfrq=0,-12,-3,1,1
mfilt=1,2,8,48,48

fincl2='Z3','PS','U','V','OMEGA','T','Q','Nudge_U','Nudge_V','Nudge_T','Nudge_Q','Target_U','Target_V','Target_T','Target_Q','PSL','U10','TS','TREFHT','TAUX','TAUY','QREFHT','USTAR','PBLH','SHFLX','LHFLX','PRECC','PRECL','PRECSC','PRECSL','FLDS','FLNS','FLNSC','FSNS','FSDS','FSDSC'

fincl3='PRECL','PRECC','PRECT','TMQ','FLUT'

fincl4='T','U','V','Q','PS','TS','TREFHT','U10','TAUX','TAUY','QREFHT','USTAR','PBLH','SHFLX','LHFLX','TAUBLJX','TAUBLJY','TAUGWX','TAUGWY','PRECC','PRECL','PRECSC','PRECSL','KVH','KVM','CLDLIQ','CLDICE','FLDS','FLNS','FLNSC','FSNS','FSDS','FSDSC','TTGWORO','VTGWORO','UTGWORO','Lscale','em','Kh_zt','Kh_zm','upwp','vpwp','wprtp','wpthlp','tau_zm'
fincl4lonlat='97.49w_36.61n','4.927e_51.971n','156.61w_71.32n','149.89w_70.50n','26.63e_67.37n','11.92e_78.92n','135.07w_60.71n','68.51w_63.74n','38.42w_72.58n','39.18e_88.75n'

fincl5='Nudge_U','Nudge_V','Nudge_T','Nudge_Q','Target_U','Target_V','Target_T','Target_Q','T','U','V','Q','PS','TS','TREFHT','U10','TAUX','TAUY','QREFHT','USTAR','PBLH','SHFLX','LHFLX','TAUBLJX','TAUBLJY','TAUGWX','TAUGWY','PRECC','PRECL','PRECSC','PRECSL','KVH','KVM','CLDLIQ','CLDICE','FLDS','FLNS','FLNSC','FSNS','FSDS','FSDSC','TTGWORO','VTGWORO','UTGWORO','Lscale','em','Kh_zt','Kh_zm','upwp','vpwp','wprtp','wpthlp','tau_zm','ICEFRAC','SNOWHICE'
!fincl5='T','U','V','Q','PS','TS','TREFHT','U10','TAUX','TAUY','QREFHT','USTAR','PBLH','SHFLX','LHFLX','TAUBLJX','TAUBLJY','TAUGWX','TAUGWY','PRECC','PRECL','PRECSC','PRECSL','KVH','KVM','CLDLIQ','CLDICE','FLDS','FLNS','FLNSC','FSNS','FSDS','FSDSC','TTGWORO','VTGWORO','UTGWORO','Lscale','em','Kh_zt','Kh_zm','upwp','vpwp','wprtp','wpthlp','tau_zm','ICEFRAC','SNOWHICE'
fincl5lonlat='0.001e:160e_81.5n:90n','340e:359.999e_81.5n:90n','270e:359.000e_87n:90n'

collect_column_output=.false.,.false.,.false.,.true.,.true.

avgflag_pertape='A','I','I','I','I'

!! ---- NUDGING START
Nudge_Model         = .true.
Nudge_Path          = '/glade/scratch/zarzycki/ndg/se_ne30pg3_L58/ERA5/'
Nudge_File_Template = 'ndg.ERA5.ne30pg3.L58.cam2.i.%y-%m-%d-%s.nc'
Nudge_Force_Opt     = 1
Nudge_TimeScale_Opt = 0
Nudge_Times_Per_Day = 8
Model_Times_Per_Day = 48
Nudge_Uprof         = 2
Nudge_Ucoef         = 1
Nudge_Vprof         = 2
Nudge_Vcoef         = 1
Nudge_Tprof         = 2
Nudge_Tcoef         = 1
Nudge_Qprof         = 0
Nudge_Qcoef         = 0
Nudge_PSprof        = 0
Nudge_PScoef        = 0
Nudge_Beg_Year      = 2017
Nudge_Beg_Month     = 1
Nudge_Beg_Day       = 1
Nudge_End_Year      = 2020
Nudge_End_Month     = 12
Nudge_End_Day       = 31
Nudge_Hwin_lat0     = 0.0
Nudge_Hwin_lon0     = 180.
Nudge_Hwin_latWidth = 9999.0
Nudge_Hwin_lonWidth = 9999.0
Nudge_Hwin_latDelta = 1.0
Nudge_Hwin_lonDelta = 1.0
Nudge_Hwin_Invert   = .false.
Nudge_Vwin_Hindex   = 39.
Nudge_Vwin_Hdelta   = 1.000
Nudge_Vwin_Lindex   = 0.
Nudge_Vwin_Ldelta   = 0.001
Nudge_Vwin_Invert   = .false.
!! ---- NUDGING END

!! Updated CMZ with cam_20230623
clubb_history=.true.
clubb_vars_zt = 'thlm', 'T_in_K', 'thvm', 'rtm', 'rcm', 'rfrzm', 'rvm', 'rel_humidity', 'um', 'vm', 'wm_zt', 'um_ref', 'vm_ref', 'ug', 'vg', 'cloud_frac', 'ice_supersat_frac', 'rcm_in_layer', 'rcm_in_cloud', 'cloud_cover', 'p_in_Pa', 'exner', 'rho_ds_zt', 'thv_ds_zt', 'Lscale', 'thlm_forcing', 'thlm_mc', 'rtm_forcing', 'rtm_mc', 'rvm_mc', 'rcm_mc', 'rcm_sd_mg_morr', 'thlm_mfl_min', 'thlm_mfl_max', 'thlm_enter_mfl', 'thlm_exit_mfl', 'thlm_old', 'thlm_without_ta', 'rtm_mfl_min', 'rtm_mfl_max', 'rtm_enter_mfl', 'rtm_exit_mfl', 'rtm_old', 'rtm_without_ta', 'wp3', 'wpup2', 'wpvp2', 'thlp3', 'rtp3', 'wpthlp2', 'wp2thlp', 'wprtp2', 'wp2rtp', 'Lscale_up', 'Lscale_down', 'Lscale_pert_1', 'Lscale_pert_2', 'tau_zt', 'invrs_tau_zt', 'Kh_zt', 'wp2thvp', 'wp2rcp', 'w_up_in_cloud', 'w_down_in_cloud', 'cld_updr_frac', 'cld_downdr_frac', 'wprtpthlp', 'rc_coef', 'sigma_sqd_w_zt', 'rho', 'Ncm', 'Nc_in_cloud', 'Nc_activated', 'Nccnm', 'Nim', 'snowslope', 'Nsm', 'Ngm', 'sed_rcm', 'rsat', 'rsati', 'rrm', 'rsm', 'rim', 'rgm', 'Nrm', 'm_vol_rad_rain', 'm_vol_rad_cloud', 'eff_rad_cloud', 'eff_rad_ice', 'eff_rad_snow', 'eff_rad_rain', 'eff_rad_graupel', 'precip_rate_zt', 'radht', 'radht_LW', 'radht_SW', 'diam', 'mass_ice_cryst', 'rcm_icedfs', 'u_T_cm', 'rtm_bt', 'rtm_ma', 'rtm_ta', 'rtm_mfl', 'rtm_tacl', 'rtm_cl', 'rtm_sdmp', 'rtm_pd', 'thlm_bt', 'thlm_ma', 'thlm_sdmp', 'thlm_ta', 'thlm_mfl', 'thlm_tacl', 'thlm_cl', 'wp3_bt', 'wp3_ma', 'wp3_ta', 'wp3_tp', 'wp3_ac', 'wp3_bp1', 'wp3_pr_tp', 'wp3_pr_turb', 'wp3_pr_dfsn', 'wp3_pr1', 'wp3_pr2', 'wp3_pr3', 'wp3_dp1', 'wp3_sdmp', 'wp3_cl', 'wp3_splat', 'rtp3_bt', 'rtp3_tp', 'rtp3_ac', 'rtp3_dp', 'thlp3_bt', 'thlp3_tp', 'thlp3_ac', 'thlp3_dp', 'rrm_bt', 'rrm_ma', 'rrm_sd', 'rrm_ts', 'rrm_sd_morr', 'rrm_ta', 'rrm_evap', 'rrm_auto', 'rrm_accr', 'rrm_evap_adj', 'rrm_src_adj', 'rrm_mc_nonadj', 'rrm_hf', 'rrm_wvhf', 'rrm_cl', 'rrm_mc', 'Nrm_bt', 'Nrm_ma', 'Nrm_sd', 'Nrm_ts', 'Nrm_ta', 'Nrm_evap', 'Nrm_auto', 'Nrm_evap_adj', 'Nrm_src_adj', 'Nrm_cl', 'Nrm_mc', 'rsm_bt', 'rsm_ma', 'rsm_sd', 'rsm_sd_morr', 'rsm_ta', 'rsm_mc', 'rsm_hf', 'rsm_wvhf', 'rsm_cl', 'Nsm_bt', 'Nsm_ma', 'Nsm_sd', 'Nsm_ta', 'Nsm_mc', 'Nsm_cl', 'rim_bt', 'rim_ma', 'rim_sd', 'rim_sd_mg_morr', 'rim_ta', 'rim_mc', 'rim_hf', 'rim_wvhf', 'rim_cl', 'rgm_bt', 'rgm_ma', 'rgm_sd', 'rgm_sd_morr', 'rgm_ta', 'rgm_mc', 'rgm_hf', 'rgm_wvhf', 'rgm_cl', 'Ngm_bt', 'Ngm_ma', 'Ngm_sd', 'Ngm_ta', 'Ngm_mc', 'Ngm_cl', 'Nim_bt', 'Nim_ma', 'Nim_sd', 'Nim_ta', 'Nim_mc', 'Nim_cl', 'Ncm_bt', 'Ncm_ma', 'Ncm_act', 'Ncm_ta', 'Ncm_mc', 'Ncm_cl', 'PSMLT', 'EVPMS', 'PRACS', 'EVPMG', 'PRACG', 'PGMLT', 'MNUCCC', 'PSACWS', 'PSACWI', 'QMULTS', 'QMULTG', 'PSACWG', 'PGSACW', 'PRD', 'PRCI', 'PRAI', 'QMULTR', 'QMULTRG', 'MNUCCD', 'PRACI', 'PRACIS', 'EPRD', 'MNUCCR', 'PIACR', 'PIACRS', 'PGRACS', 'PRDS', 'EPRDS', 'PSACR', 'PRDG', 'EPRDG', 'NGSTEN', 'NRSTEN', 'NISTEN', 'NSSTEN', 'NCSTEN', 'NPRC1', 'NRAGG', 'NPRACG', 'NSUBR', 'NSMLTR', 'NGMLTR', 'NPRACS', 'NNUCCR', 'NIACR', 'NIACRS', 'NGRACS', 'NSMLTS', 'NSAGG', 'NPRCI', 'NSCNG', 'NSUBS', 'PRC', 'PRA', 'PRE', 'PCC', 'NNUCCC', 'NPSACWS', 'NPRA', 'NPRC', 'NPSACWI', 'NPSACWG', 'NPRAI', 'NMULTS', 'NMULTG', 'NMULTR', 'NMULTRG', 'NNUCCD', 'NSUBI', 'NGMLTG', 'NSUBG', 'NACT', 'SIZEFIX_NR', 'SIZEFIX_NC', 'SIZEFIX_NI', 'SIZEFIX_NS', 'SIZEFIX_NG', 'NEGFIX_NR', 'NEGFIX_NC', 'NEGFIX_NI', 'NEGFIX_NS', 'NEGFIX_NG', 'NIM_MORR_CL', 'QC_INST', 'QR_INST', 'QI_INST', 'QS_INST', 'QG_INST', 'NC_INST', 'NR_INST', 'NI_INST', 'NS_INST', 'NG_INST', 'T_in_K_mc', 'w_KK_evap_covar_zt', 'rt_KK_evap_covar_zt', 'thl_KK_evap_covar_zt', 'w_KK_auto_covar_zt', 'rt_KK_auto_covar_zt', 'thl_KK_auto_covar_zt', 'w_KK_accr_covar_zt', 'rt_KK_accr_covar_zt', 'thl_KK_accr_covar_zt', 'rr_KK_mvr_covar_zt', 'Nr_KK_mvr_covar_zt', 'KK_mvr_variance_zt', 'vm_bt', 'vm_ma', 'vm_gf', 'vm_cf', 'vm_ta', 'vm_f', 'vm_sdmp', 'vm_ndg', 'vm_mfl', 'um_bt', 'um_ma', 'um_gf', 'um_cf', 'um_ta', 'um_f', 'um_sdmp', 'um_ndg', 'um_mfl', 'mixt_frac', 'w_1', 'w_2', 'varnce_w_1', 'varnce_w_2', 'thl_1', 'thl_2', 'varnce_thl_1', 'varnce_thl_2', 'rt_1', 'rt_2', 'varnce_rt_1', 'varnce_rt_2', 'rc_1', 'rc_2', 'rsatl_1', 'rsatl_2', 'cloud_frac_1', 'cloud_frac_2', 'chi_1', 'chi_2', 'stdev_chi_1', 'stdev_chi_2', 'chip2', 'stdev_eta_1', 'stdev_eta_2', 'covar_chi_eta_1', 'covar_chi_eta_2', 'corr_w_chi_1', 'corr_w_chi_2', 'corr_w_eta_1', 'corr_w_eta_2', 'corr_chi_eta_1', 'corr_chi_eta_2', 'corr_w_rt_1', 'corr_w_rt_2', 'corr_w_thl_1', 'corr_w_thl_2', 'corr_rt_thl_1', 'corr_rt_thl_2', 'crt_1', 'crt_2', 'cthl_1', 'cthl_2', 'F_w', 'F_rt', 'F_thl', 'min_F_w', 'max_F_w', 'min_F_rt', 'max_F_rt', 'min_F_thl', 'max_F_thl', 'coef_wprtp2_implicit', 'term_wprtp2_explicit', 'coef_wpthlp2_implicit', 'term_wpthlp2_explicit', 'coef_wprtpthlp_implicit', 'term_wprtpthlp_explicit', 'coef_wp2rtp_implicit', 'term_wp2rtp_explicit', 'coef_wp2thlp_implicit', 'term_wp2thlp_explicit', 'wp2_zt', 'thlp2_zt', 'wpthlp_zt', 'wprtp_zt', 'rtp2_zt', 'rtpthlp_zt', 'up2_zt', 'vp2_zt', 'upwp_zt', 'vpwp_zt', 'Skw_zt', 'Skthl_zt', 'Skrt_zt', 'rcm_supersat_adj', 'C11_Skw_fnc', 'chi', 'a3_coef_zt', 'wp3_on_wp2_zt', 'precip_frac', 'precip_frac_1', 'precip_frac_2', 'Ncnm', 'cloud_frac_refined', 'rcm_refined', 'hl_on_Cp_residual', 'qto_residual'
clubb_vars_zm = 'wp2', 'rtp2', 'thlp2', 'rtpthlp', 'wprtp', 'wpthlp', 'wp3_zm', 'thlp3_zm', 'rtp3_zm', 'wp2up2', 'wp2vp2', 'wp4', 'wpthvp', 'rtpthvp', 'thlpthvp', 'tau_zm', 'invrs_tau_zm', 'invrs_tau_xp2_zm', 'invrs_tau_wp2_zm', 'invrs_tau_wpxp_zm', 'invrs_tau_wp3_zm', 'invrs_tau_no_N2_zm', 'invrs_tau_bkgnd', 'invrs_tau_sfc', 'invrs_tau_shear', 'Kh_zm', 'wprcp', 'rc_coef_zm', 'thlprcp', 'rtprcp', 'rcp2', 'upwp', 'vpwp', 'upthlp', 'uprtp', 'vpthlp', 'vprtp', 'upthvp', 'uprcp', 'vpthvp', 'vprcp', 'rho_zm', 'sigma_sqd_w', 'rho_ds_zm', 'thv_ds_zm', 'em', 'shear', 'mean_w_up', 'mean_w_down', 'Frad', 'wpNcp', 'VNr', 'Vrr', 'VNc', 'Vrc', 'VNs', 'Vrs', 'Vrg', 'VNi', 'Vri', 'Vrrprrp', 'VNrpNrp', 'Vrrprrp_expcalc', 'VNrpNrp_expcalc', 'wp2_bt', 'wp2_ma', 'wp2_ta', 'wp2_ac', 'wp2_bp', 'wp2_pr1', 'wp2_pr2', 'wp2_pr3', 'wp2_pr_dfsn', 'wp2_dp1', 'wp2_dp2', 'wp2_sdmp', 'wp2_cl', 'wp2_pd', 'wp2_sf', 'wp2_splat', 'wprtp_bt', 'wprtp_ma', 'wprtp_ta', 'wprtp_tp', 'wprtp_ac', 'wprtp_bp', 'wprtp_pr1', 'wprtp_pr2', 'wprtp_pr3', 'wprtp_dp1', 'wprtp_mfl', 'wprtp_cl', 'wprtp_sicl', 'wprtp_pd', 'wprtp_forcing', 'wprtp_mc', 'wpthlp_bt', 'wpthlp_ma', 'wpthlp_ta', 'wpthlp_tp', 'wpthlp_ac', 'wpthlp_bp', 'wpthlp_pr1', 'wpthlp_pr2', 'wpthlp_pr3', 'wpthlp_dp1', 'wpthlp_mfl', 'wpthlp_cl', 'wpthlp_sicl', 'wpthlp_forcing', 'wpthlp_mc', 'upwp_bt', 'upwp_ma', 'upwp_ta', 'upwp_tp', 'upwp_ac', 'upwp_bp', 'upwp_pr1', 'upwp_pr2', 'upwp_pr3', 'upwp_pr4', 'upwp_dp1', 'upwp_mfl', 'upwp_cl', 'vpwp_bt', 'vpwp_ma', 'vpwp_ta', 'vpwp_tp', 'vpwp_ac', 'vpwp_bp', 'vpwp_pr1', 'vpwp_pr2', 'vpwp_pr3', 'vpwp_pr4', 'vpwp_dp1', 'vpwp_mfl', 'vpwp_cl', 'rtp2_bt', 'rtp2_ma', 'rtp2_ta', 'rtp2_tp', 'rtp2_dp1', 'rtp2_dp2', 'rtp2_cl', 'rtp2_pd', 'rtp2_sf', 'rtp2_forcing', 'rtp2_mc', 'thlp2_bt', 'thlp2_ma', 'thlp2_ta', 'thlp2_tp', 'thlp2_dp1', 'thlp2_dp2', 'thlp2_cl', 'thlp2_pd', 'thlp2_sf', 'thlp2_forcing', 'thlp2_mc', 'rtpthlp_bt', 'rtpthlp_ma', 'rtpthlp_ta', 'rtpthlp_tp1', 'rtpthlp_tp2', 'rtpthlp_dp1', 'rtpthlp_dp2', 'rtpthlp_cl', 'rtpthlp_sf', 'rtpthlp_forcing', 'rtpthlp_mc', 'up2', 'vp2', 'up2_bt', 'up2_ma', 'up2_ta', 'up2_tp', 'up2_dp1', 'up2_dp2', 'up2_pr1', 'up2_pr2', 'up2_sdmp', 'up2_cl', 'up2_pd', 'up2_sf', 'up2_splat', 'vp2_bt', 'vp2_ma', 'vp2_ta', 'vp2_tp', 'vp2_dp1', 'vp2_dp2', 'vp2_pr1', 'vp2_pr2', 'vp2_sdmp', 'vp2_cl', 'vp2_pd', 'vp2_sf', 'vp2_splat', 'wpthlp_enter_mfl', 'wpthlp_exit_mfl', 'wpthlp_mfl_min', 'wpthlp_mfl_max', 'wprtp_mfl_min', 'wprtp_mfl_max', 'wprtp_enter_mfl', 'wprtp_exit_mfl', 'wm_zm', 'cloud_frac_zm', 'ice_supersat_frac_zm', 'rcm_zm', 'rtm_zm', 'thlm_zm', 'w_1_zm', 'w_2_zm', 'varnce_w_1_zm', 'varnce_w_2_zm', 'mixt_frac_zm', 'Skw_velocity', 'gamma_Skw_fnc', 'C6rt_Skw_fnc', 'C6thl_Skw_fnc', 'C6_term', 'C7_Skw_fnc', 'C1_Skw_fnc', 'coef_wp4_implicit', 'bv_freq_sqd', 'sqrt_Ri_zm', 'Richardson_num', 'shear_sqd', 'a3_coef', 'wp3_on_wp2', 'Skw_zm', 'Skthl_zm', 'Skrt_zm', 'stability_correction', 'rtp2_from_chi', 'lh_rtp2_mc', 'lh_thlp2_mc', 'lh_wprtp_mc', 'lh_wpthlp_mc', 'lh_rtpthlp_mc'

ext_frc_specifier = 'H2O    -> /glade/p/cesmdata/cseg/inputdata/atm/cam/chem/emis/elev/H2OemissionCH4oxidationx2_3D_L70_1849-2101_CMIP6ensAvg_SSP3-7.0_c190403.nc',
         'num_a1 -> /glade/p/cesmdata/cseg/inputdata/atm/cam/chem/emis/emissions_ssp370/emissions-cmip6-ScenarioMIP_IAMC-AIM-ssp370-1-1_num_so4_a1_anthro-ene_vertical_mol_175001-210101_0.9x1.25_c20190222.nc',
         'num_a1 -> /glade/p/cesmdata/cseg/inputdata/atm/cam/chem/emis/CMIP6_emissions_1750_2015/emissions-cmip6_num_a1_so4_contvolcano_vertical_850-5000_0.9x1.25_c20170724.nc',
         'num_a2 -> /glade/p/cesmdata/cseg/inputdata/atm/cam/chem/emis/CMIP6_emissions_1750_2015/emissions-cmip6_num_a2_so4_contvolcano_vertical_850-5000_0.9x1.25_c20170724.nc',
         'SO2    -> /glade/p/cesmdata/cseg/inputdata/atm/cam/chem/emis/CMIP6_emissions_1750_2015/emissions-cmip6_SO2_contvolcano_vertical_850-5000_0.9x1.25_c20170724.nc',
         'so4_a1 -> /glade/p/cesmdata/cseg/inputdata/atm/cam/chem/emis/emissions_ssp370/emissions-cmip6-ScenarioMIP_IAMC-AIM-ssp370-1-1_so4_a1_anthro-ene_vertical_mol_175001-210101_0.9x1.25_c20190222.nc',
         'so4_a1 -> /glade/p/cesmdata/cseg/inputdata/atm/cam/chem/emis/CMIP6_emissions_1750_2015/emissions-cmip6_so4_a1_contvolcano_vertical_850-5000_0.9x1.25_c20170724.nc',
         'so4_a2 -> /glade/p/cesmdata/cseg/inputdata/atm/cam/chem/emis/CMIP6_emissions_1750_2015/emissions-cmip6_so4_a2_contvolcano_vertical_850-5000_0.9x1.25_c20170724.nc'
ext_frc_type = 'INTERP_MISSING_MONTHS'

srf_emis_specifier = 'bc_a4    -> /glade/p/cesmdata/cseg/inputdata/atm/cam/chem/emis/emissions_ssp370/emissions-cmip6-ScenarioMIP_IAMC-AIM-ssp370-1-1_bc_a4_anthro_surface_mol_175001-210101_0.9x1.25_c20190222.nc',
         'bc_a4    -> /glade/p/cesmdata/cseg/inputdata/atm/cam/chem/emis/emissions_ssp370/emissions-cmip6-ScenarioMIP_IAMC-AIM-ssp370-1-1_bc_a4_bb_surface_mol_175001-210101_0.9x1.25_c20190222.nc',
         'DMS      -> /glade/p/cesmdata/cseg/inputdata/atm/cam/chem/emis/emissions_ssp370/emissions-cmip6-ScenarioMIP_IAMC-AIM-ssp370-1-1_DMS_bb_surface_mol_175001-210101_0.9x1.25_c20190222.nc',
         'DMS      -> /glade/p/cesmdata/cseg/inputdata/atm/cam/chem/emis/emissions_ssp370/emissions-cmip6-SSP_DMS_other_surface_mol_175001-210101_0.9x1.25_c20190222.nc',
         'num_a1   -> /glade/p/cesmdata/cseg/inputdata/atm/cam/chem/emis/emissions_ssp370/emissions-cmip6-ScenarioMIP_IAMC-AIM-ssp370-1-1_num_so4_a1_bb_surface_mol_175001-210101_0.9x1.25_c20190222.nc',
         'num_a1   -> /glade/p/cesmdata/cseg/inputdata/atm/cam/chem/emis/emissions_ssp370/emissions-cmip6-ScenarioMIP_IAMC-AIM-ssp370-1-1_num_so4_a1_anthro-ag-ship_surface_mol_175001-210101_0.9x1.25_c20190222.nc',
         'num_a2   -> /glade/p/cesmdata/cseg/inputdata/atm/cam/chem/emis/emissions_ssp370/emissions-cmip6-ScenarioMIP_IAMC-AIM-ssp370-1-1_num_so4_a2_anthro-res_surface_mol_175001-210101_0.9x1.25_c20190222.nc',
         'num_a4   -> /glade/p/cesmdata/cseg/inputdata/atm/cam/chem/emis/emissions_ssp370/emissions-cmip6-ScenarioMIP_IAMC-AIM-ssp370-1-1_num_bc_a4_bb_surface_mol_175001-210101_0.9x1.25_c20190222.nc',
         'num_a4   -> /glade/p/cesmdata/cseg/inputdata/atm/cam/chem/emis/emissions_ssp370/emissions-cmip6-ScenarioMIP_IAMC-AIM-ssp370-1-1_num_bc_a4_anthro_surface_mol_175001-210101_0.9x1.25_c20190222.nc',
         'num_a4   -> /glade/p/cesmdata/cseg/inputdata/atm/cam/chem/emis/emissions_ssp370/emissions-cmip6-ScenarioMIP_IAMC-AIM-ssp370-1-1_num_pom_a4_anthro_surface_mol_175001-210101_0.9x1.25_c20190222.nc',
         'num_a4   -> /glade/p/cesmdata/cseg/inputdata/atm/cam/chem/emis/emissions_ssp370/emissions-cmip6-ScenarioMIP_IAMC-AIM-ssp370-1-1_num_pom_a4_bb_surface_mol_175001-210101_0.9x1.25_c20190222.nc',
         'pom_a4   -> /glade/p/cesmdata/cseg/inputdata/atm/cam/chem/emis/emissions_ssp370/emissions-cmip6-ScenarioMIP_IAMC-AIM-ssp370-1-1_pom_a4_anthro_surface_mol_175001-210101_0.9x1.25_c20190222.nc',
         'pom_a4   -> /glade/p/cesmdata/cseg/inputdata/atm/cam/chem/emis/emissions_ssp370/emissions-cmip6-ScenarioMIP_IAMC-AIM-ssp370-1-1_pom_a4_bb_surface_mol_175001-210101_0.9x1.25_c20190222.nc',
         'SO2      -> /glade/p/cesmdata/cseg/inputdata/atm/cam/chem/emis/emissions_ssp370/emissions-cmip6-ScenarioMIP_IAMC-AIM-ssp370-1-1_SO2_anthro-ag-ship-res_surface_mol_175001-210101_0.9x1.25_c20190222.nc',
         'SO2      -> /glade/p/cesmdata/cseg/inputdata/atm/cam/chem/emis/emissions_ssp370/emissions-cmip6-ScenarioMIP_IAMC-AIM-ssp370-1-1_SO2_anthro-ene_surface_mol_175001-210101_0.9x1.25_c20190222.nc',
         'SO2      -> /glade/p/cesmdata/cseg/inputdata/atm/cam/chem/emis/emissions_ssp370/emissions-cmip6-ScenarioMIP_IAMC-AIM-ssp370-1-1_SO2_bb_surface_mol_175001-210101_0.9x1.25_c20190222.nc',
         'so4_a1   -> /glade/p/cesmdata/cseg/inputdata/atm/cam/chem/emis/emissions_ssp370/emissions-cmip6-ScenarioMIP_IAMC-AIM-ssp370-1-1_so4_a1_anthro-ag-ship_surface_mol_175001-210101_0.9x1.25_c20190222.nc',
         'so4_a1   -> /glade/p/cesmdata/cseg/inputdata/atm/cam/chem/emis/emissions_ssp370/emissions-cmip6-ScenarioMIP_IAMC-AIM-ssp370-1-1_so4_a1_bb_surface_mol_175001-210101_0.9x1.25_c20190222.nc',
         'so4_a2   -> /glade/p/cesmdata/cseg/inputdata/atm/cam/chem/emis/emissions_ssp370/emissions-cmip6-ScenarioMIP_IAMC-AIM-ssp370-1-1_so4_a2_anthro-res_surface_mol_175001-210101_0.9x1.25_c20190222.nc',
         'SOAG     -> /glade/p/cesmdata/cseg/inputdata/atm/cam/chem/emis/emissions_ssp370/emissions-cmip6-ScenarioMIP_IAMC-AIM-ssp370-1-1_SOAGx1.5_anthro_surface_mol_175001-210101_0.9x1.25_c20200403.nc',
         'SOAG     -> /glade/p/cesmdata/cseg/inputdata/atm/cam/chem/emis/emissions_ssp370/emissions-cmip6-ScenarioMIP_IAMC-AIM-ssp370-1-1_SOAGx1.5_bb_surface_mol_175001-210101_0.9x1.25_c20200403.nc',
         'SOAG     -> /glade/p/cesmdata/cseg/inputdata/atm/cam/chem/emis/emissions_ssp/emissions-cmip6-SOAGx1.5_biogenic_surface_mol_175001-210101_0.9x1.25_c20190329.nc'
srf_emis_type = 'INTERP_MISSING_MONTHS'

prescribed_ozone_datapath = '/glade/p/cesmdata/cseg/inputdata/atm/cam/ozone_strataero'
prescribed_ozone_file = 'ozone_strataero_WACCM_L70_zm5day_18500101-21010201_CMIP6histEnsAvg_SSP370_c190403.nc'
prescribed_ozone_name = 'O3'
prescribed_ozone_type = 'SERIAL'

prescribed_strataero_datapath = '/glade/p/cesmdata/cseg/inputdata/atm/cam/ozone_strataero'
prescribed_strataero_file = 'ozone_strataero_WACCM_L70_zm5day_18500101-21010201_CMIP6histEnsAvg_SSP370_c190403.nc'
prescribed_strataero_type = 'SERIAL'
prescribed_strataero_use_chemtrop =  .true.

solar_irrad_data_file = '/glade/p/cesmdata/cseg/inputdata/atm/cam/solar/spectral_irradiance_Lean_1610-2140_ann_c100408.nc'

flbc_file = '/glade/p/cesmdata/cseg/inputdata/atm/waccm/lb/LBC_2014-2500_CMIP6_SSP370_0p5degLat_GlobAnnAvg_c190301.nc'
flbc_list = 'CO2','CH4','N2O','CFC11eq','CFC12'
flbc_type = 'SERIAL'
scenario_ghg = 'CHEM_LBC_FILE'

tracer_cnst_datapath = '/glade/p/cesmdata/cseg/inputdata/atm/cam/tracer_cnst'
tracer_cnst_file = 'tracer_cnst_halons_3D_L70_1849-2101_CMIP6ensAvg_SSP3-7.0_c190403.nc'
tracer_cnst_type = 'INTERP_MISSING_MONTHS'
```

#### user\_nl\_clm

```
finidat="/glade/p/univ/upsu0032/lnd-inic/FHIST-ne30-ndg-ERA5-Q24-N23-x001.clm2.r.2018-01-01-00000.nc"
init_interp_fill_missing_with_natveg = .true.
use_init_interp = .true.
```

Then setup the case, build the model, and run:

```
./case.setup

./case.build

./case.submit
```

A couple notes:

- A full description of the nudging namelist options is [HERE](https://ncar.github.io/CAM/doc/build/html/users_guide/physics-modifications-via-the-namelist.html)
- `Nudge_Vwin_Hindex` tells you the level index at which the nudging "stops" when counting from the top. However, `Nudge_Vwin_Hdelta` greater than 0.001 tapers the nudging a bit, such that XXXX.
- `Model_Times_Per_Day` should be set to `ATM_NCPL` from `env_run.xml`
- `Nudge_Times_Per_Day` should be set to 24 (hourly), 8 (3-hourly), 4 (6-hourly), 2 (12-hourly) or 1 (daily). Chris Kruse found that nudging throughout the atmosphere all the way to the surface with 24-hourly caused problems near orography. Nudging to a certain level seems less problematic, although if 3-hourly data works, that would seem to be preferred.
- `Nudge_?prof` should always be set to `2` unless you want nudging off (`0`) or just want nudging everywhere (`1`).
- The `Nudge_Beg_Year` and `Nudge_End_Year`, etc. options should cover the range of nudging files you have generated. You can also modify them to have the nudging only be performed during parts of a simulation (i.e., specify a range smaller than the available nudging data range). However, if you specify a range that is larger than the available nudging range *and* attempt to get the model to integrate in that "bad" part of the range, it will fail.

### Generating nudging profile

You can evaluate the nudging profile by using an NCL tool. You can modify the `user_nl_cam` file iteratively and then place back in the case directory.

```
cd $BETACAST/atm_to_cam/nudging
cp ~/PATHTOCASEDIR/user_nl_cam .
ncl Lookat_NudgeWindow.ncl NLEV=58
```

### Generating nudging files

The first step is to create nudging files from reanalysis. On Cheyenne this is easy using Betacast and the existing ERA5.

```
## edit ndg.era5.nl
qsub -v NLFILE="ndg.era5.nl" gen-nudge.sh
```
