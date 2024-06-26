begin

;inputres="180x360"
delimiter="x"

print("GEN_SST_DOMAIN: inputres: "+inputres)


parsed_str=str_split(inputres,delimiter)

nlat=toint(parsed_str(0))
nlon=toint(parsed_str(1))

dlat=180.0/todouble(nlat)
dlon=360.0/todouble(nlon)

print("GEN_SST_DOMAIN: nlat: "+nlat+"  nlon: "+nlon)
print("GEN_SST_DOMAIN: dlat: "+dlat+"  dlon: "+dlon)

dlatd2 = dlat/2.
dlond2 = dlon/2.

lat = todouble(fspan(-90.+dlatd2,90.-dlatd2,nlat))
lon = todouble(fspan(dlond2,360.-dlond2,nlon))

lat!0="nj"
lat@units="degrees_north"
lat@long_name="latitude"

lon!0="ni"
lon@units="degrees_east"
lon@long_name="longitude"

storedir="./domains/"
system("mkdir -vp "+storedir)
filename=storedir+"/domain.ocn."+nlat+delimiter+nlon+".nc"
system("/bin/rm -fv "+filename)   ; remove any pre-existing file
print("GEN_SST_DOMAIN: Writing: "+filename)
ncdf = addfile(filename,"c")  ; open output netCDF file

fAtt               = True            ; assign file attributes
fAtt@title         = "Betacast-generated DOCN domain file"
fAtt@resolution    =  nlat+"x"+nlon
fAtt@creation_date = systemfunc ("date")
fileattdef( ncdf, fAtt )            ; copy file attributes

ncdf->yc = lat
ncdf->xc = lon

; write out SCRIP grid

filename = str_sub_str(filename,"domain.ocn","scrip.ocn")
Opt=True
Opt@Debug=True
Opt@ForceOverwrite=True
rectilinear_to_SCRIP(filename,lat,lon,Opt)

print("GEN_SST_DOMAIN: done")

status_exit(9)

end
