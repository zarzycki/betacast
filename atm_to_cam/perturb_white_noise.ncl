begin

pertMag=0.002
;basFileName="BC_ne30_CAM5_climo.cam.r.0003-01-01-00000.nc"

basFile = addfile(basFileName,"w")

pertMagPerc=pertMag*100.
print("We are using perturbation magnitude of: "+pertMagPerc+"%")

ps = basFile->PS
t = basFile->T
q = basFile->Q

random_setallseed(tointeger(systemfunc("date +%s")),tointeger(systemfunc("date +%s"))+5)

; perturb T
low   = 1.0-pertMag
high  = 1.0+pertMag
dimT  = dimsizes(t)
unf   = random_uniform(low, high, dimT)
print("orig Tmax: "+max(t))
t     = t * unf
print("new Tmax: "+max(t))
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
basFile->T=t
basFile->Q=q

print("... done")
status_exit(9)

end
