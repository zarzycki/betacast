function buildTimeLookup(filelist)
; Pass in a filelist of CAM files, build a lookup table for filename and time indices
local ftmp, nfiles, nftimes, timeArr, listArr, ii
begin
  ; Here we build a "lookup" table to find what times are on what files...
  ; get list of h3 files
  ftmp = addfile(filelist(0),"r")
  nfiles = dimsizes(filelist)
  nftimes = dimsizes(ftmp->time)
  timeArr = new((/nfiles,nftimes/),typeof(ftmp->time))
  listArr = new((/nfiles/),"string")
  do ii = 0,nfiles-1
    print("loading time for file #"+ii+" of "+(nfiles-1)+"  "+filelist(ii))
    ftmp = addfile(filelist(ii),"r")
    nftimes=dimsizes(ftmp->time)   ; needed since files at end of series may have less than rest
    timeArr(ii,0:nftimes-1) = ftmp->time
    listArr(ii)=filelist(ii)
  end do
  timeArr!0="files"
  timeArr&files=listArr

  timeunits=ftmp->time@units
  opt=0
  opt@calendar=ftmp->time@calendar

  timeArr@units=timeunits
  timeArr@calendar=ftmp->time@calendar
  delete(ftmp)
  return(timeArr)
end

function getTimeLookup(filelist,stashedfile)
local f
begin

  f = addfile(stashedfile,"r")

  timeArr = f->timeArr
  timeArr&files=filelist

  delete(f)
  return(timeArr)

end

function getFileAndTimeFromLookup(timeArr,thisTime)
local timeArr, thisTime, time1d, dsizes_time, indices, fileix, trackindex, needed_file, timediff,eps
begin

  ; small number for checking time sync issues
  eps=1.0e-6

  ; look up time ix
  time1d      = ndtooned(timeArr)
  dsizes_time = dimsizes(timeArr)

  ; use this to get *exact match* for thisTime to some time in the time1d array
  ;timediff=0.0
  ;indices  = ind_resolve(ind(time1d.eq.thisTime),dsizes_time)

  ; here, find the index that is *nearest* to the thisTime
  ; which allows for slightly out of sync dates (e.g., 1 hour)
  timediff=abs(time1d-thisTime)
  indices  = ind_resolve(minind(timediff),dsizes_time)
  if (min(timediff) .gt. eps) then
    print("WARNING: diff b/w best matched time in nc files & time from traj file = "+min(timediff))
    print("WARNING: thisTime is: "+cd_calendar(thisTime,3))
    ;print(time1d+" "+timediff)
    print("---> this may mean your times are out of sync!")
    return_error = toint(min(timediff))
  else
    return_error = 0
  end if

  if (any(ismissing(indices))) then
    trackindex=-1
    needed_file=""
  else
    fileix     = indices(0,0)
    trackindex = indices(0,1)
    needed_file=timeArr&files(fileix)
  end if
  return( [/ needed_file,trackindex,return_error /] )
end

;----------------------------------------------------------------------------------------------------

begin

; given a directory of a bunch of parent files, find the single file containing the timestep needed
; for initializing a child simulation from a parent file
; return this as a single line ASCII file called m2mfile.UQSTR

;DIR="/Users/cmz5202/NetCDF/TGW/aux/"
;YYYYMMDDHH=1979091100
;UQSTR=uniqtime

if (.not.isvar("DIR")) then
  print("find-time-file.ncl: user did not define DIR, exiting")
  exit
end if
if (.not.isvar("YYYYMMDDHH")) then
  print("find-time-file.ncl: user did not define YYYYMMDDHH, exiting")
  exit
end if
if (.not.isvar("UQSTR")) then
  print("find-time-file.ncl: user did not define UQSTR, exiting")
  exit
end if

print("Finding file in: "+DIR)

FILES=systemfunc ("ls "+DIR+"*.nc")

timeArr=buildTimeLookup(FILES)

; ===== Getting date from YYYYMMDDHH
dtime_map = (/4,2,2,2/)
splitDate = str_split_by_length(tostring(YYYYMMDDHH),dtime_map)
yyyy=toint(splitDate(0))
mm=toint(splitDate(1))
dd=toint(splitDate(2))
hh=toint(splitDate(3))
mn=0
sc=0
print("Finding time: "+yyyy+" "+mm+" "+dd+" "+hh+" "+mn+" "+sc)

Opt=0
Opt@calendar=timeArr@calendar

desired_time = cd_inv_calendar(yyyy,mm,dd,hh,mn,sc,timeArr@units, Opt)

a = getFileAndTimeFromLookup(timeArr,desired_time)

theFile=a[0]

if .not.ismissing(theFile) .and. isfilepresent(theFile) then
  print("Found "+YYYYMMDDHH+" in file: "+a[0])
  asciiwrite ("m2mfile."+UQSTR, a[0])
else
  print("find-time-file.ncl: Unable to correctly find file containing "+YYYYMMDDHH)
  exit
end if

; exit cleanly for betacast
status_exit(9)

end
