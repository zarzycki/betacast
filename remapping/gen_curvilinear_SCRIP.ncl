load "$NCARG_ROOT/lib/ncarg/nclscripts/esmf/ESMF_regridding.ncl"

begin

    InterpMethod= "patch"                      ; define interpolation method

;---Input file
    srcDirName  = "./"
    srcFileName = "rap_130_20120629_0300_000.grb2"
    srcFilePath =  srcDirName + srcFileName
    srcGridName = "./anl_scrip/rap_13km_scrip.nc"
    
    sfile = addfile(srcFilePath,"r")

    lat2d = sfile->gridlat_0
    lon2d = sfile->gridlon_0
    
    var = sfile->MSLMA_P0_L101_GLC0
    
    sfile = addfile(srcFilePath,"r")
    Opt                = True
    ;Opt@GridMask       = where(.not.ismissing(var),1,0)    ; use "Mask2D" in NCL V6.2.1 and earlier
    Opt@ForceOverwrite = True
    Opt@PrintTimings   = True
    Opt@Title          = "RAP Grid"
Opt@Debug = True

    curvilinear_to_SCRIP(srcGridName,lat2d,lon2d,Opt)
    
end

