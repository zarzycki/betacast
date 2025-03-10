begin

top_template="./L56template.nc"

tfile=addfile(top_template,"r")
thyam  = tfile->hyam
thybm  = tfile->hybm
thyai  = tfile->hyai
thybi  = tfile->hybi
tlev  = tfile->lev
tilev = tfile->ilev

tnlev = dimsizes(tlev)
tnint = dimsizes(tilev)

keeplev = new(tnlev,"logical")
keepilev = new(tnint,"logical")

keepilev = True

;delLevs = (/tnint-2,tnint-3,tnint-5,tnint-6,tnint-8/)
delLevs = (/tnint-2/)

do ii = 0,tnint-1
  if (any(delLevs.eq.ii)) then
    keepilev(ii)=False
  end if
end do

ohyai = thyai(ind(keepilev))
ohybi = thybi(ind(keepilev))
oilev = tilev(ind(keepilev))

onint = dimsizes(ohyai)
onlev = onint-1

ohyam  = new(onlev,typeof(thyam))
ohybm  = new(onlev,typeof(thybm))
olev   = new(onlev,typeof(tlev))

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
out_template="./Z_L"+onlev+"template.nc"
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


