load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/esmf/ESMF_regridding.ncl"

begin

;---Input file
srcFileName = "ERA-Int_2005101800.nc"

;---Output (and input) files
srcGridName = "ERA-Int_SCRIP.nc"

;----------------------------------------------------------------------
; Convert ERA-Interim grid to a SCRIP convention file.
;----------------------------------------------------------------------
src_file = addfile(srcFileName,"r")
temp     = src_file->t

Opt                = True
Opt@ForceOverwrite = True
Opt@PrintTimings   = True
Opt@Title          = "NCEP Grid"

rectilinear_to_SCRIP(srcGridName,temp&latitude,temp&longitude,Opt)

;    Opt                = True
;    Opt@LLCorner       = (/ -60.d,   0.d/)
;    Opt@URCorner       = (/  60.d, 355.d/)
;    Opt@ForceOverwrite = True
;    Opt@PrintTimings   = True

;    latlon_to_SCRIP(srcGridName,"5x5",Opt)

delete(Opt)

;----------------------------------------------------------------------
; Get SE var-res SCRIP file
;----------------------------------------------------------------------

dstGridName = "/glade/p/work/zarzycki/grids/scrip/uniform_60.g_scrip.nc"

;----------------------------------------------------------------------
; Generate the weights that take you from the NCEP grid to a
; 5x5 degree grid.
;----------------------------------------------------------------------

wgtFileName = "ERA_to_uniform_60_patch.nc"

Opt                      = True
Opt@InterpMethod         = "patch"     ; default
Opt@ForceOverwrite       = True
Opt@PrintTimings         = True
Opt@Debug = True

ESMF_regrid_gen_weights(srcGridName,dstGridName,wgtFileName,Opt)

delete(Opt)
    
end
