load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/esmf/ESMF_regridding.ncl"

begin

inttype="double" ; double, quadruple
keeplowestlevel=True

in_template="./L30template.nc"
ifile=addfile(in_template,"r")

ihyam  = ifile->hyam
ihybm  = ifile->hybm
ihyai  = ifile->hyai
ihybi  = ifile->hybi
ilev = ifile->lev
iilev = ifile->ilev

inlev = dimsizes(ilev)
inint = dimsizes(iilev)
onlev = 2*inlev
onint = onlev+1

if (keeplowestlevel) then
  onlev = onlev-1
  onint = onint-1
end if

print("Number of input levels = "+inlev)
print("Number of input interfaces = "+inint)
print("Number of output levels = "+onlev)
print("Number of output interfaces = "+onint)

ohyam  = new(onlev,typeof(ihyam))
ohybm  = new(onlev,typeof(ihybm))
ohyai  = new(onint,typeof(ihyai))
ohybi  = new(onint,typeof(ihybi))
olev   = new(onlev,typeof(ilev))
oilev  = new(onint,typeof(iilev))

do ii = 0,onint-1
  if (ii .eq. 0) then
    ohyai(0) = ihyai(0)
    ohybi(0) = ihybi(0)
    oilev(0) = iilev(0)
  else if (ii .eq. onint-1) then
    ohyai(ii) = ihyai(inint-1)
    ohybi(ii) = ihybi(inint-1)
    oilev(ii) = iilev(inint-1)
  else
    if (mod(ii,2) .eq. 0) then
      ohyai(ii) = ihyai(ii/2)
      ohybi(ii) = ihybi(ii/2)
      oilev(ii) = iilev(ii/2)
    else
      ohyai(ii) = ( ihyai((ii+1)/2) + ihyai((ii-1)/2) ) / 2.
      ohybi(ii) = ( ihybi((ii+1)/2) + ihybi((ii-1)/2) ) / 2.
      oilev(ii) = ( iilev((ii+1)/2) + iilev((ii-1)/2) ) / 2.
    end if
  end if
  end if
end do

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
out_template="./L"+onlev+"template.nc"
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
ofile->lev  = olev
ofile->ilev  = oilev

end


