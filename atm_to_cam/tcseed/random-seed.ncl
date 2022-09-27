load "../../functions/namelist.ncl"

begin

; pass in filename and cyclones_tempest.tmp

;filename="cyclones_tempest.tmp"
;pthi="seed.ape.nl"
minlat=5.
maxlat=20.
minlon=0.0000001
maxlon=359.99999
minrmw=150000.
maxrmw=450000.
mindp=15.0
maxdp=40.0
THRESHOLD=15.

invert_vortex  = keyword_values( pthi, "invert_vortex", "logical")
if (invert_vortex)
  print("invert_vortex is set to "+invert_vortex)
  print("We don't need to worry about where to seed vortices!")
  status_exit(9)
end if

;NROWS pass in via command line


nlines = toint(systemfunc("head -100 "+filename+" | wc -l"))
print("nlines: "+nlines)
print("-------------------------------------")

; nrows x ncol
data = asciiread(filename,(/nlines,4/),"float")

lats=data(0:nlines-1,2)
lons=data(0:nlines-1,1)

print("lats        lons")
print(lats+"     "+lons)

print("-------------------------------------")

rseed1 = toint(systemfunc("date +%s%N | cut -b11-22"))
rseed2 = toint(systemfunc("date +%s%N | cut -b13-22"))
random_setallseed(rseed1, rseed2)
    

DONE=False
ITERMAX=100
ii=0
do while(.not.DONE)
  ii = ii+1
  thislat = random_uniform(minlat,maxlat,1)
  thislon = random_uniform(minlon,maxlon,1)
  print("lat/lon: "+thislat+" "+thislon)
  gcdist = gc_latlon(thislat,thislon,lats,lons,2,2)
  print("minimum seperation is: "+min(gcdist))
  if (min(gcdist) .gt. THRESHOLD) then
    DONE=True
  end if
  if (ii .ge. ITERMAX) then
    print("at iter "+ii+"  could not find random point outside of THRESHOLD "+THRESHOLD)
    exit
  end if
end do

; add this storm to the filename of existing tracks
quote = inttochar(34)
system("echo "+quote+"     9999    "+thislon+"     "+thislat+"      9.99999e+11"+quote+" >> "+filename)

; update seed NL
;minp=950.
;target_rmw=-1.

target_rmw = random_uniform(minrmw,maxrmw,1)
minp = random_uniform(mindp,maxdp,1)
print("random target_rmw: "+target_rmw)
print("random minp: "+minp)
; flip dp since seed code only directly applies dp if negative
inv_minp = -1*minp*100.

print("updating namelist file")
system("sed -i "+quote+"s?.*psminlat=.*?psminlat="+thislat+"?"+quote+" "+pthi)
system("sed -i "+quote+"s?.*psminlon=.*?psminlon="+thislon+"?"+quote+" "+pthi)
system("sed -i "+quote+"s?.*target_rmw=.*?target_rmw="+target_rmw+"?"+quote+" "+pthi)
system("sed -i "+quote+"s?.*minp=.*?minp="+inv_minp+"?"+quote+" "+pthi)

status_exit(9)

end 
