begin

; these are passed in as nl options
;basedir="/glade/p/univ/upsu0032/DATM_FORCING/ANOM_LENS/"
;current_year = 2019

print("basedir: "+basedir)
print("current_year: "+current_year)

; -------

start_year   = 1920
start_month  = 1
end_month    = 12
shouldermon  = 1
debug = False
do_compress = True

ens_vars = (/"FLDS","PRECT","QBOT","TBOT"/)
internal_vars = (/"FLDS","PRECT","Q","T"/)
operators=(/"ratio","ratio","diff","diff"/)

do ii = 0,dimsizes(ens_vars)-1

  newfile = basedir+"/ens_"+ens_vars(ii)+"_"+current_year+"ref_anom.nc"

  ; check if newfile exists .. if not, we need to build
  if (.not. isfilepresent(newfile))

    ; print some diagnostics...
    print("ANOM_DIAG: "+current_year+" "+ens_vars(ii)+" "+internal_vars(ii)+" "+operators(ii))

    ; load data
    f=addfile(basedir+"/ens_"+ens_vars(ii)+"_anom.nc","r")

    ; load the reference year
    refyr_stix=current_year*12-start_year*12+start_month-1-shouldermon
    refyr_enix=current_year*12-start_year*12+end_month-1
    print("ANOM_DIAG: refyr_stix: "+refyr_stix)
    print("ANOM_DIAG: refyr_enix: "+refyr_enix)
    ref_year  = f->$internal_vars(ii)$(refyr_stix:refyr_enix,:,:)

    printVarSummary(ref_year)

    ; average across time for the reference year
    ref_year_avg = dim_avg_n_Wrap(ref_year,0)

    ; load the existing anom for the whole datastream
    new_anom = f->$internal_vars(ii)$

    ; update the anom based on the ref year avg
    if (operators(ii) .eq. "diff") then
      new_anom = new_anom - conform(new_anom,ref_year_avg,(/1,2/))
    else
      new_anom = new_anom / conform(new_anom,ref_year_avg,(/1,2/))
    end if

    if (debug) then
      system("/bin/rm -f simple_"+ens_vars(ii)+".nc")   ; remove any pre-existing file
      ncdf = addfile("simple_"+ens_vars(ii)+".nc" ,"c")  ; open output netCDF file

      fAtt               = True            ; assign file attributes
      fAtt@creation_date = systemfunc ("date")
      fileattdef( ncdf, fAtt )            ; copy file attributes

      filedimdef(ncdf,"time",-1,True)

      ncdf->new_anom = new_anom
      ncdf->ref_year_removal = ref_year_removal
      ncdf->$internal_vars(ii)$ = f->$internal_vars(ii)$
    end if

    print("Writing "+newfile)

    system("/bin/rm -f "+newfile)   ; remove any pre-existing file
    ncdf = addfile(newfile ,"c")  ; open output netCDF file

    fAtt               = True            ; assign file attributes
    fAtt@creation_date = systemfunc ("date")
    fileattdef( ncdf, fAtt )            ; copy file attributes

    filedimdef(ncdf,"time",-1,True)

    ncdf->$internal_vars(ii)$ = new_anom
    ; lat, lon, time come over attached to the new anom var
    ;ncdf->lat = f->lat
    ;ncdf->lon = f->lon
    ;ncdf->time = f->time
    ; have to drag timebnd over
    ncdf->time_bnds=f->time_bnds

    if (do_compress) then
      print("... compressing")
      systemcall="ncks -4 -L 1 -O "+newfile+" "+newfile
      print(" "+systemcall)
      system(systemcall)
      print("... success!")
    end if

  else
    print(newfile+" already exists!")
  end if

end do

end
