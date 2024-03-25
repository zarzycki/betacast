begin

; add a pertMag of white noise to PS, Q, and some 3D state field
; usage: ncl perturb_white_noise.ncl 'basFileName="/glade/derecho/scratch/zarzycki/proj3_f45_4xco2.cam.r.0091-01-01-00000.nc"'

pertMag=0.001

basFile = addfile(basFileName,"w")
print("We are perturbing: "+basFileName)

pertMagPerc=pertMag*100.
print("We are using perturbation magnitude of: "+pertMagPerc+"%")

; get Q
q = basFile->Q
; get PS
ps = basFile->PS
if ( isfilevar(basFile, "T") ) then
  var3d="T"
else if ( isfilevar(basFile, "U") ) then
  var3d="U"
else
  print("no 3D var!")
  exit
end if
end if
; get 3D state var
print("Using this 3D var: "+var3d)
t = basFile->$var3d$

random_setallseed(tointeger(systemfunc("date +%s")),tointeger(systemfunc("date +%s"))+5)

; perturb T
low   = 1.0-pertMag
high  = 1.0+pertMag
dimT  = dimsizes(q)
unf   = random_uniform(low, high, dimT)

print("orig "+var3d+"max: "+max(t))
t     = t * unf
print("new "+var3d+"max: "+max(t))
unf   = random_uniform(low, high, dimT)
print("orig Qmax: "+max(q))
q     = q * unf
print("new Qmax: "+max(q))

delete(unf)

dimP  = dimsizes(ps)
unf   = random_uniform(low, high, dimP)
print("orig PSmax: "+max(ps))
ps     = ps * unf
print("new PSmax: "+max(ps))

basFile->PS=ps
basFile->$var3d$=t
basFile->Q=q

print("... done")
status_exit(9)

end
