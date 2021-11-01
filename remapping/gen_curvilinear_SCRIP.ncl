load "$NCARG_ROOT/lib/ncarg/nclscripts/esmf/ESMF_regridding.ncl"

begin

  srcFileName = "/glade/u/home/zarzycki/work/2012_derecho/hrrr.t00z.wrfprsf00.grib2"
  srcGridName = "./anl_scrip/hrrr_3km_scrip_tmp.nc"
  
  sfile = addfile(srcFileName,"r")

  lat2d = sfile->gridlat_0
  lon2d = sfile->gridlon_0
      
  Opt                = True
  ;Opt@GridMask       = where(.not.ismissing(var),1,0)    ; use "Mask2D" in NCL V6.2.1 and earlier
  Opt@ForceOverwrite = True
  Opt@PrintTimings   = True
  Opt@Title          = "RAP Grid"
  Opt@Debug = True
  
  curvilinear_to_SCRIP(srcGridName,lat2d,lon2d,Opt)
    
end

