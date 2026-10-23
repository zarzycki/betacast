function rainratetodbz(rainrate,units)

;========================================================
; This function takes in rainfall rate (here in mm/hr)
; and converts it to DbZ for plotting as simulated
; radar reflectivity
;
; rainrate = 2D (or 3D in time) array of integrated
; rainfall rate
; units, 0 = mm/hr, 1 = m/s, 2 cm/s
;
; Colin Zarzycki, zarzycki@umich.edu, 11/14/2012
;========================================================

begin

; If units of m/s, convert to mm/hr
if (units .eq. 1) then
    rainrate = rainrate*3600000
else if (units .eq. 2) then
    rainrate = rainrate*36000
else
end if
end if

dbz = rainrate
eps = 10^-25        ; very small number to prevent log from = 0
; Uses common NWS NEXRAD published conversion formula
;dbz = 10*log10(eps + (200 * rainrate^(1.6) ) )

; Z = 250 R^1.2
dbz = 10*log10(eps + (250*rainrate^(1.2)))

; if any NaN crop up in array
; set to fill value
if (any(isnan_ieee(dbz))) then
      value = -25000.000
      replace_ieeenan (dbz, value, 0)
      dbz@_FillValue = value
end if

dbz@units = "dBz"
dbz@long_name = "Simulated Radar Reflectivity"

return(dbz)

end

