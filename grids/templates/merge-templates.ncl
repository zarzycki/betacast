load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/esmf/ESMF_regridding.ncl"

begin

keeplowestlevel=False

top_template="./L30template.nc"
bot_template="./L59template.nc"
splitlev=600.

tfile=addfile(top_template,"r")
thyam  = tfile->hyam
thybm  = tfile->hybm
thyai  = tfile->hyai
thybi  = tfile->hybi
tlev  = tfile->lev
tilev = tfile->ilev

bfile=addfile(bot_template,"r")
bhyam  = bfile->hyam
bhybm  = bfile->hybm
bhyai  = bfile->hyai
bhybi  = bfile->hybi
blev  = bfile->lev
bilev = bfile->ilev

tnlev = dimsizes(tlev)
tnint = dimsizes(tilev)
bnlev = dimsizes(blev)
bnint = dimsizes(bilev)

onlev = max((/bnlev,tnlev/))
onint = max((/bnint,tnint/))

print("Number of input levels = "+tnlev)
print("Number of input interfaces = "+tnint)
print("Number of input levels = "+bnlev)
print("Number of input interfaces = "+bnint)
print("Number of output levels = "+onlev)
print("Number of output interfaces = "+onint)

tilev = where(tilev .lt. splitlev,tilev,-1)
bilev = where(bilev .ge. splitlev,bilev,-1)

indtopi = ind(tilev .gt. 0)
indboti = ind(bilev .gt. 0)

print(indtopi)
print(indboti)

numtopi = dimsizes(indtopi)
numboti = dimsizes(indboti)

onint = numtopi + numboti
onlev = onint - 1

ohyam  = new(onlev,typeof(thyam))
ohybm  = new(onlev,typeof(thybm))
ohyai  = new(onint,typeof(thyai))
ohybi  = new(onint,typeof(thybi))
olev   = new(onlev,typeof(tlev))
oilev  = new(onint,typeof(tilev))

oilev(0:numtopi-1) = tilev(indtopi)
oilev(numtopi:onint-1) = bilev(indboti)

ohyai(0:numtopi-1) = thyai(indtopi)
ohyai(numtopi:onint-1) = bhyai(indboti)

ohybi(0:numtopi-1) = thybi(indtopi)
ohybi(numtopi:onint-1) = bhybi(indboti)

; now calculate midpoints, which are just literally the midpoints

do jj = 0,onlev-1
  ohyam(jj) = (ohyai(jj) + ohyai(jj+1)) / 2.
  ohybm(jj) = (ohybi(jj) + ohybi(jj+1)) / 2.
   olev(jj) = (oilev(jj) + oilev(jj+1)) / 2.
end do

olev!0="lev"
olev&lev=olev

copy_VarMeta(olev,ohyam)
copy_VarMeta(olev,ohybm)

;===================================================================
; write output template
;===================================================================
out_template="./_L"+onlev+"template.nc"
system("/bin/rm -f "+out_template)   ; remove any pre-existing file
ofile = addfile(out_template ,"c")  ; open output netCDF file

;fAtt               = True            ; assign file attributes
;fAtt@title         = "CAM vertical coord file "+onlev+" levels"
;fAtt@source_file   = ifile
;fAtt@creation_date = systemfunc ("date")
;fAtt@Conventions   = "None"
;fileattdef( ofile, fAtt )            ; copy file attributes

;===================================================================
; output variables directly; NCL will call appropriate functions
; to write the meta data associated with each variable
;===================================================================

ofile->hyam  = ohyam
ofile->hybm  = ohybm
ofile->hyai  = ohyai
ofile->hybi  = ohybi
ofile->lev   = olev
ofile->ilev  = todouble(tofloat(oilev))

end


