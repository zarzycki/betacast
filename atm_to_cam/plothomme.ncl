load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"

begin
;******************************************************
;  Open a netCDF file containing the grid and data
;  from SE, formerly HOMME
;******************************************************

f2=addfile("/glade/work/zarzycki/sewx//INIC//MERGE.nc","r")

lat1d = ndtooned(f2->lat)          ; Pull off lat/lon variable and
lon1d = ndtooned(f2->lon)          ; convert both to 1D.
var  = ndtooned(f2->PS)      ; Read some data; convert to 1D.

print("max var: "+max(var))
print("min var: "+min(var))
;******************************************************
; create plot
;******************************************************
wks = gsn_open_wks("x11","homme")               ; open workstation
gsn_define_colormap(wks,"BlAqGrYeOrReVi200")   ; choose colormap

res                     = True         ; plot modifications desired
;res@gsnMaximize         = True         ; Maximize size of plot in frame
res@gsnSpreadColors     = True         ; Use full colormap, but start

res@sfXArray            = lon1d        ; Required to tell NCL where to
res@sfYArray            = lat1d        ; overlay data on globe.

res@cnFillOn            = True         ; Turn on contour fill
res@cnFillMode          = "RasterFill"   ; Style of fill. You can also
                                         ; use "CellFill" and "RasterFill"
res@cnLinesOn           = False        ; Turn off contour lines
res@cnLineLabelsOn      = False        ; Turn off contour line labels

res@lbLabelAutoStride   = True         ; Clean up labelbar labels.

map = gsn_csm_contour_map_ce(wks,var,res)

end
