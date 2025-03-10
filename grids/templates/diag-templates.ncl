load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/esmf/ESMF_regridding.ncl"

begin

top_template="./Z_L55template.nc"

tfile=addfile(top_template,"r")
thyam  = tfile->hyam
thybm  = tfile->hybm
thyai  = tfile->hyai
thybi  = tfile->hybi
tlev  = tfile->lev
tilev = tfile->ilev

tnlev = dimsizes(tlev)
tnint = dimsizes(tilev)

print("Number of input levels = "+tnlev)
print("Number of input interfaces = "+tnint)

ilev_thick = tilev(1:tnint-1) - tilev(0:tnint-2)

print(tilev(1:tnint-1)+"    "+tlev+"     "+ilev_thick)


end


