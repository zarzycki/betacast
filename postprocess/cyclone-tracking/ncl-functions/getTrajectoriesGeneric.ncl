function getTrajectories(tfile,ncol,maxLines,headerDelimStr,tab_space,isHeader)

; tfile (string) -- name of trajectory output from Tempest/RECAP

; ncol (int) -- number of columns for each trajectory point

; maxLines (int) -- max # of lines allocated, but if set to negative integer, will try and
; determine internally. Only use hardcoded (e.g., 150) if you segfault or other OB errors

; headerDelimStr -- substring of separate trajectory headers (ex: start, Path), used to parse trajectories

; tab_space (string) -- delimiter of fields for each traj point. "tab" will convert to tab-delimited

; isHeader (boolean) -- does the trajectory file have a header? If yes, ditch first row of file.

begin

; read in trajectory data from ascii file
; get number of storms by counting number of lines with 'start'
nstorms_tot = stringtoint(systemfunc("grep '"+headerDelimStr+"' "+tfile+" | wc -l"))
print("GETTRAJECTORIES: Getting trajectories from: "+tfile)
print("GETTRAJECTORIES: Total number of storms = "+nstorms_tot)

tdata_tmp = asciiread(tfile,-1,"string")  ; read in trajectory data as strings
dims = dimsizes(tdata_tmp)    ; get dimensions
nrows = dims(0)               ; get number of rows

delete([/dims/])

; ====DO NOT TOUCH THIS==============
if (tab_space .eq. "tab") then
  tab_space="	 "
end if

if (isHeader) then
  startloopix=1
else
  startloopix=0
end if

; ===================================


; =============== FIGURE OUT MAX NUMBERS OF STORM TIMES ========

if (maxLines .le. 0)
  print("GETTRAJECTORIES: No max lines defined, trying to figure out!")
  randDate = systemfunc("date +%s%N")
  ML_filetmpFile="supertempcolinhi"+randDate+".txt"
  system("grep -n '"+headerDelimStr+"' "+tfile+" |cut -f1 -d: > "+ML_filetmpFile)
  ML_tmpStormRowsData = asciiread(ML_filetmpFile,-1,"integer")
  if (dimsizes(ML_tmpStormRowsData) .gt. 1) then
    ML_tmpStormRowsDiff = (ML_tmpStormRowsData(1:dimsizes(ML_tmpStormRowsData)-1)-ML_tmpStormRowsData(0:dimsizes(ML_tmpStormRowsData)-2))-1
  else
    ML_tmpStormRowsDiff = 0
  end if
  ML_finalStormSt=systemfunc("grep -n '"+headerDelimStr+"' "+tfile+" | tail -1 | cut -f1 -d:")
  ML_finalStormEn=systemfunc("wc -l < "+tfile)
  ML_lastStormSize = toint(ML_finalStormEn)-toint(ML_finalStormSt)
  maxLines = max((/max(ML_tmpStormRowsDiff),ML_lastStormSize/))
  system("rm "+ML_filetmpFile)
end if

print("GETTRAJECTORIES: Allocating... "+maxLines+" max lines per storm")

; ===============================================================

; create new array that will hold all data except rows containing start information
; add column to contain storm number
tdata = new((/nstorms_tot,ncol,maxLines/),float)

; loop over lines of trajectory data, convert string to float
nstorm = -1   ; assign storm number
n2 = 0       ; counter data without header lines

do n1 = startloopix,nrows-1  ; loop over all data with header lines
  headerCheck = str_match(tdata_tmp(n1),headerDelimStr)
  if (ismissing(headerCheck))
    headerLine=False
  else
    headerLine=True
  end if

  if (headerLine) then              ; header line
    ;print("Header line")
    nstorm = nstorm + 1                 ; set storm number
    n2 = 0                              ; reset "per indiv traj" counter
  else
    splitStr = str_split(tdata_tmp(n1),tab_space)
    do z = 0,ncol-1
      ;print(nstorm+" "+n2+" "+z+"   "+splitStr(z))
      ;print(" "+tdata_tmp(n1))
      tdata(nstorm,z,n2) = stringtofloat(splitStr(z))
    end do
    n2 = n2 + 1
      delete(splitStr)
  end if
end do

print("Done!")

delete([/tdata_tmp,n1,n2,nstorms_tot,nstorm,nrows/])

return(tdata)

end
