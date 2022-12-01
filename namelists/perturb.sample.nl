! ***************** GENERAL SETTINGS
case         = "CESMLENS"
basedir      = "/global/cfs/cdirs/m2637/betacast/deltas/"
start_month  = 1
end_month    = 1
current_year = 1996
comp_year    = 2081

! ***************** ATM ONLY
correct_sfc       = False
plevs             = False
update_pressure   = False
update_winds      = False 
do_ps_corr        = True 
esmf_remap        = True 
keep_esmf         = True
smooth_deltas     = False 
smooth_delta_iter = 10
output_atm_diag   = True
extra_diags_atm   = False

! ***************** SST ONLY
adjust_ice        = True  
output_sst_diag   = True





