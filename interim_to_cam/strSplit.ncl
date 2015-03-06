load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/esmf/ESMF_regridding.ncl"

begin

YYYYMMDDHH = "2014082112"
dtime_map = (/4,2,2,2/)
splitDate = str_split_by_length(YYYYMMDDHH,dtime_map)
yearstr=splitDate(0)
monthstr=splitDate(1)
daystr=splitDate(2)
cyclestr=splitDate(3)
print("Regridding analysis from: "+yearstr+" "+monthstr+" "+daystr+" "+cyclestr+"Z")

end
