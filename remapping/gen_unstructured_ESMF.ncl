load "$NCARG_ROOT/lib/ncarg/nclscripts/esmf/ESMF_regridding.ncl"

begin

    InterpMethod= "patch"                      ; define interpolation method

;---Input file
    srcDirName  = "/glade/u/home/zarzycki/work/2012_derecho/"
    srcFileName = "ERA5.rrm.x4.ini.2012062906.nc"
    
    srcGridName = "./rrm.x4_ESMF.nc"

    sfile = addfile(srcDirName+"/"+srcFileName,"r")

    lat = sfile->lat
    lon = sfile->lon
        
    sfile = addfile(srcDirName+"/"+srcFileName,"r")
    Opt                = True
    Opt@ForceOverwrite = True
    Opt@PrintTimings   = True
    Opt@Debug = True

    unstructured_to_ESMF(srcGridName,lat,lon,Opt)
    
end

