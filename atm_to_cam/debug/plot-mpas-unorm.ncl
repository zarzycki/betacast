f = addfile("/glade/p/cesmdata/cseg/inputdata/atm/cam/inic/mpas/mpasa120.CFSR.L32.nc","r")
grid = addfile("/glade/p/cesmdata/cseg/inputdata/atm/cam/inic/mpas/mpasa120.CFSR.L32.nc","r")

RAD2DEG = get_r2d("double")      ; Radian to Degree

lonEdge = grid->lonEdge
latEdge = grid->latEdge
lonCell = grid->lonCell
latCell = grid->latCell

uNormal = f->u(0,:,0)

;---Start the graphics
wks = gsn_open_wks("x11","mpas_cmz")             ; send graphics to PNG file

lonCell = lonCell*RAD2DEG
latCell = latCell*RAD2DEG
lonEdge = lonEdge*RAD2DEG
latEdge = latEdge*RAD2DEG

cmap = read_colormap_file("WhiteBlueGreenYellowRed")

res                      = True              ; Plot mods desired.

res@gsnMaximize          = True              ; Maximize plot

res@cnFillOn             = True              ; color plot desired
res@cnFillMode           = "RasterFill"      ; turn raster on      
res@cnFillPalette        = cmap(48:208,:)       ; Don't use white
res@cnLinesOn            = False             ; turn off contour lines
res@cnLineLabelsOn       = False             ; turn off contour labels

;res@trGridType           = "TriangularMesh"  ; This is required to allow
                                           ; missing coordinates.

;res@cnLevelSelectionMode = "ManualLevels"
;res@cnMinLevelValF       = 55
;res@cnMaxLevelValF       = 100
;res@cnLevelSpacingF      = 2.5

;res@mpMinLonF = 50. ; select a subregion
;res@mpMaxLonF = 70.
;res@mpMinLatF = -25.
;res@mpMaxLatF = -10.

res@mpFillOn             = False

res@sfXArray             = lonEdge
res@sfYArray             = latEdge
res@gsnAddCyclic         = False

res@tiMainString         = "Surface pressure on MPAS grid (" + \
                         dimsizes(uNormal) + " cells)"

plot = gsn_csm_contour_map(wks,uNormal,res)  
