begin
	out_type = "png"

	;gridfile = "test_NewOrleans_TC_grid_ne32x8.g"

	if (.not.isfilepresent(gridfile)) then
		print((/"Can not find "+gridfile+"!"/))
		exit
	end if

  ; name for writing output file
	gridname = gridfile

	print((/"Plotting mesh from "+gridfile/))
	system("ln -s "+gridfile+" "+gridfile+".nc")

	pi4= atan(1.d);
	pi2 = acos(0.d)
	pi = pi2*2.d
	f = addfile(gridfile+".nc","r")
	num_el_blk=dimsizes(f->eb_status)

	if ((num_el_blk.ne.6).and.(num_el_blk.ne.1)) then
		print("Neither cube nor sphere, can't plot!")
		system("/bin/rm -f "+gridfile+".nc")
		exit
	end if

	system("mkdir -p ./pdfs")

	if (num_el_blk.eq.6) then
		IsCube = True
	else
		IsCube = False
	end if

	; cartesian coords:
	coordx = f->coord(0,:)
	coordy = f->coord(1,:)
	coordz = f->coord(2,:)

	; convert to lat/lon (project from cube to sphere!)
	if (IsCube) then
		tol = 1e-8
		print("Projecting coordinates onto sphere...")
		n_nodes = dimsizes(coordx)
		r = sqrt(coordx*coordx + coordy*coordy + coordz*coordz)
		lat = new(n_nodes,"double")
		lon = new(n_nodes,"double")
		face = new(n_nodes,"double")
		do i=0,n_nodes-1
			face(i) = -1
			if ((coordx(i).le.-pi4+tol).and.(face(i).eq.-1)) then
				face(i) = 2
				alpha = -coordy(i)
				beta = coordz(i)
			end if
			if ((coordx(i).ge.pi4-tol).and.(face(i).eq.-1)) then
				face(i) = 0
				alpha = coordy(i)
				beta = coordz(i)
			end if
			if ((coordy(i).le.-pi4+tol).and.(face(i).eq.-1)) then
				face(i) = 3
				alpha = coordx(i)
				beta = coordz(i)
			end if
			if ((coordy(i).ge.pi4-tol).and.(face(i).eq.-1)) then
				face(i) = 1
				alpha = -coordx(i)
				beta = coordz(i)
			end if
			if ((coordz(i).le.-pi4+tol).and.(face(i).eq.-1)) then
				face(i) = 5
				alpha = coordy(i)
				beta = coordx(i)
			end if
			if ((coordz(i).ge.pi4-tol).and.(face(i).eq.-1)) then
				face(i) = 4
				alpha = coordy(i)
				beta = -coordx(i)
			end if
			if (face(i).eq.-1) then
				print("Error: face not assigned to a node!")
				print((/i+" "+coordx(i)+" "+coordy(i)+" "+coordz(i)/))
				system("rm -f "+gridfile+".nc")
				exit
			end if
			if (face(i).lt.4) then
				lon(i) = alpha+pi2*face(i)
				lat(i) = atan(tan(beta)*cos(alpha))
			end if
			if (face(i).eq.4) then
				lon(i) = atan2(tan(alpha), -tan(beta))
				if (tan(alpha)*tan(alpha)+tan(beta)*tan(beta).lt.tol) then
					lat(i) = pi2
				else
					lat(i) = atan(1.d/sqrt(tan(alpha)*tan(alpha)+tan(beta)*tan(beta)))
				end if
			end if
			if (face(i).eq.5) then
				lon(i) = atan2(tan(alpha), tan(beta));-pi
				if (tan(alpha)*tan(alpha)+tan(beta)*tan(beta).lt.tol) then
					lat(i) = -pi2
				else
					lat(i) = -atan(1.d/sqrt(tan(alpha)*tan(alpha)+tan(beta)*tan(beta)))
				end if
			end if
		end do
		print("... done!")
	else
		lat = asin(coordz)
		lon = atan2(coordy,coordx)
	end if
	lat = lat*180.d/pi
	lon = lon*180.d/pi

	print("number of verticies = "+dimsizes(lon))
	print("lat min/max = "+min(lat)+" "+max(lat))

	dimnames = getvardims(f)
	dimsz = getfiledimsizes(f)
	do i=0,dimsizes(dimnames)-1
		if ( dimnames(i).eq."num_elem")
			print("dimname[i]="+dimnames(i)+" "+dimsz(i))
			nelem=dimsz(i)
		end if
	end do

	map = new ( (/nelem,4/), "integer", -99 )

	count=0
	fcnt = dimsizes(f->connect1)
	print("dimsizes(map1)="+fcnt)
	n=fcnt(0)
	map(count:(n-1),:) = (/ f->connect1 /)
	count=count+n

	if (IsCube)
		fcnt = dimsizes(f->connect2)
		print("dimsizes(map2)="+fcnt)
		n=fcnt(0)
		map(count:(count+n-1),:) = (/ f->connect2 /)
		count=count+n

		fcnt = dimsizes(f->connect3)
		print("dimsizes(map2)="+fcnt)
		n=fcnt(0)
		map(count:(count+n-1),:) = (/ f->connect3 /)
		count=count+n

		fcnt = dimsizes(f->connect4)
		print("dimsizes(map2)="+fcnt)
		n=fcnt(0)
		map(count:(count+n-1),:) = (/ f->connect4 /)
		count=count+n

		fcnt = dimsizes(f->connect5)
		print("dimsizes(map2)="+fcnt)
		n=fcnt(0)
		map(count:(count+n-1),:) = (/ f->connect5 /)
		count=count+n

		fcnt = dimsizes(f->connect6)
		print("dimsizes(map2)="+fcnt)
		n=fcnt(0)
		map(count:(count+n-1),:) = (/ f->connect6 /)
		count=count+n
	end if

	nelem=count
	print(dimsizes(map))
	print("number of elements = "+nelem)

	; Orthographic Projection
	wks = gsn_open_wks(out_type,"./pdfs/"+gridname+"_stereo")

	res = True
;	res@tiMainString = gridfile

	res@mpProjection      = "Orthographic"
        res@mpOutlineBoundarySets = "National"

; NEUS

;40.7128° N, 74.0060° W
	 res@mpCenterLatF      =  41.
	 res@mpCenterLonF      =  286.

	; USA
	; res@mpCenterLatF      =  45.
	; res@mpCenterLonF      =  260.


	; ATLANTIC
	 ;res@mpCenterLatF      =  10.
	 ;res@mpCenterLonF      =  330.

        ; EUROPE
         ;res@mpCenterLatF      =  45.
         ;res@mpCenterLonF      =  20.

	; Africa
	;res@mpCenterLatF      =  15.
	;res@mpCenterLonF      =  120.

        ; Pacfic
        ;res@mpCenterLatF      =  15.
        ;res@mpCenterLonF      =  100.

        ; res@mpCenterLatF      =  15.
        ; res@mpCenterLonF      =  200.

; RICH WESTPAC 30x4
     ; res@mpCenterLatF      =  10.
     ; res@mpCenterLonF      =  150.

; RICH INDIAN 30x4
     ; res@mpCenterLatF      =  0.
     ; res@mpCenterLonF      =  70.

; NORTH POLE
;      res@mpCenterLatF      =  90.
;      res@mpCenterLonF      =  0.

; SOUTH POLE
;      res@mpCenterLatF      =  -90.
;      res@mpCenterLonF      =  0.

; Australasia
	; res@mpCenterLatF      =  -27.
	; res@mpCenterLonF      =  147.

; New Zealand
;	 res@mpCenterLatF      =  -41.
;	 res@mpCenterLonF      =  175.

	res@vpXF      = 0.05
	res@vpYF      = 0.9
	res@vpWidthF  = 0.9
	res@vpHeightF = 0.8

	res@gsnDraw  = False       ; don't draw the plots now
	res@gsnFrame = False       ; or advance the frame

plot_outlines=True
if (plot_outlines) then
	res@mpOutlineOn            = True
	res@mpPerimOn              = False
	res@mpLandFillColor        = "tan"
	res@mpOceanFillColor       = "LightBlue"
	res@mpInlandWaterFillColor = "Blue"
	res@mpGreatCircleLinesOn = True
else
        res@mpOutlineOn            = False
        res@mpPerimOn              = False
        res@mpLandFillColor        = "White"
        res@mpOceanFillColor       = "White"
        res@mpInlandWaterFillColor = "White"
        res@mpGreatCircleLinesOn = True
end if

	gsn_define_colormap(wks,(/"white","black","tan","LightBlue","Blue"/))


		plot = gsn_csm_map(wks,res) ; create the plot

  draw(plot)

	; polygon resources
	res_p             = True

	res_p@gsLineThicknessF = 1.0
	res_p@gsLineColor   = "black"

	j1 = ispan(0,3,1)
	j2 = (j1+1)%4

	xlon = new ( (/5/), "double")
	xlat = new ( (/5/), "double")

	do i=0,nelem-1
		if ( mod(i,500).eq.0) then
			print ("i = "+i+"/"+(nelem-1))
		end if
		xlon(0:3) = lon(map(i,:)-1)
		xlat(0:3) = lat(map(i,:)-1)

		do j=0,2
			if ( abs(xlon(j+1)-xlon(j)) .gt. 180.0) then
				if (xlon(j+1) .gt. xlon(j) ) then
					xlon(j)=xlon(j)+360.
				else
					xlon(j+1)=xlon(j+1)+360.
				end if
			end if
		end do

		xlon(4)=xlon(0)
		xlat(4)=xlat(0)
		gsn_polyline(wks, plot, xlon,xlat,res_p)
	end do

;*****************************************
; plot min extrema over original data
;*****************************************

; polymarkers for cordex corners
  polyres               = True          ; poly marker mods desired
  polyres@gsMarkerIndex = 16            ; choose circle as polymarker
  polyres@gsMarkerSizeF = 8.0           ; select size to avoid streaking
  polyres@gsMarkerColor = (/"red"/)   ; choose color
  ;gsn_polymarker(wks,plot,-79.373333,43.741667,polyres)  ; draw polymarkers
  ;gsn_polymarker(wks,plot,-79.373333,47.741667,polyres)  ; draw polymarkers

  ;gsn_polymarker(wks,plot,1.92,71.84,polyres)  ; draw polymarkers
  ;gsn_polymarker(wks,plot,64.4,66.65,polyres)  ; draw polymarkers
  ;gsn_polymarker(wks,plot,338.23,42.36,polyres)  ; draw polymarkers
  ;gsn_polymarker(wks,plot,9.75,49.68,polyres)  ; draw polymarkers
  ;gsn_polymarker(wks,plot,44.77,46.72,polyres)  ; draw polymarkers
  ;gsn_polymarker(wks,plot,350.01,22.20,polyres)  ; draw polymarkers
  ;gsn_polymarker(wks,plot,12.48,27.34,polyres)  ; draw polymarkers
  ;gsn_polymarker(wks,plot,36.30,25.36,polyres)  ; draw polymarkers

if (0 .eq. 1) then
  gres                 = True
  gres@gsFillColor     = "red"
  gres@tfPolyDrawOrder = "Draw"  ; this can be used for polylines, polymarkers, or polygons
  gres@gsFillOpacityF=0.5
  gres@gsLineOpacityF=0.5
  ;polylat = (/40, 24, 24, 40/)   ; SE USA
  ;polylon = (/265, 265, 287, 287/)  ; SE USA
  ;polylat = (/-5., 25., 25., -5./)
  ;polylon = (/179., 179., 120., 120./)
  ; AUSTRALIA CORDEX (9)
  ;polylat = (/8.76, 3.87, 12.21, -18.03, -24.26, -27.90, -44.28, -52.36, -39.25/)
  ;polylon = (/110.19, 146.16, 182.02, 101.41, 147.63, 199.57, 89.25, 150.03, 206.57/)
  ; RICO / ATOMIC
  polylat = (/0.,35.,35.,0./)
  polylon = (/-110.,-110.,30.,30./)
  gsn_polygon(wks,plot,polylon,polylat,gres)

;TLC (110.19; 8.76)
;CNB (146.16; 3.87)
;TRC (182.02; 12.21)
;CWB (101.41; -18.03)
;CPD (147.63; -24.26)
;CEB (199.57; -27.90)
;BLC (89.25; -44.28)
;CSB (150.03; -52.36)
;BRC (206.57; -39.25)


end if

if (1 .eq. 0) then
  polyres               = True          ; poly marker mods desired
  polyres@gsMarkerIndex = 16            ; choose circle as polymarker
  polyres@gsMarkerSizeF = 10.0           ; select size to avoid streaking
  polyres@gsMarkerColor = (/"red"/)   ; choose color
  ; AUSTRALIA CORDEX (9)
  polylat = (/-25.,35.,35.,-25./)
  polylon = (/250.,250.,30.,30./)
  gsn_polymarker(wks,plot,polylon,polylat,polyres)

;TLC (110.19; 8.76)
;CNB (146.16; 3.87)
;TRC (182.02; 12.21)
;CWB (101.41; -18.03)
;CPD (147.63; -24.26)
;CEB (199.57; -27.90)
;BLC (89.25; -44.28)
;CSB (150.03; -52.36)
;BRC (206.57; -39.25)


end if

	if (out_type.eq."ps") then
		delete(wks)
		wks = gsn_open_wks(out_type,"./pdfs/"+gridname+"_stereo")
	else
		frame(wks)
	end if

	system("/bin/rm -f "+gridfile+".nc")

status_exit(9)

end

