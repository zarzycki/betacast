function tctestcase(cen_lon:numeric,cen_lat:numeric,dp:numeric,rp:numeric,lon[*]:numeric,lat[*]:numeric,p[*]:numeric,z[*]:numeric,zcoords:integer,psin[*]:numeric,uin[*]:numeric,vin[*]:numeric,Tin[*]:numeric,qin[*]:numeric)

;SUBROUTINE test5_tropical_cyclone (lon,lat,p,z,zcoords,u,v,w,t,phis,ps,rho,q)

begin

output = new(5,float)

;zcoords = 0 ;0 = pressure, 1 = z

;input longitude (radians) 
; 	latitude (radians) 
; 	pressure/height

; Constants
a	= 6371220.	        ; Earth's Radius (m)
Rd 	= 287.0	          ; Ideal gas const dry air (J kg^-1 K^1)
g	= 9.80616	          ; Gravity (m s^2)
omega   = 7.292115d-5 ; angular velocity 1/s
pi      = 4.*atan(1.)       ; pi
convert = 180./pi   ; conversion factor: radians to degrees
;rp         = 273000.  ; Radius for calculation of PS
;dp         = 2280.    ; Delta P for calculation of PS
zp         = 7000.    ; Height for calculation of P
q0         = 0.021    ; q at surface from Jordan
gamma_      = 0.007    ; lapse rate
Ts0        = 302.0   ; Surface temperature (SST)
p00        = 101500.  ; global mean surface pressure
p0         = 100000.  ; p for model level calculation
;cen_lat    = 20.      ; Center latitude of initial vortex
;cen_lon    = -40.     ; Center longitufe of initial vortex
zq1        = 3000.    ; Height 1 for q calculation
zq2        = 8000.    ; Height 2 for q calculation
exppr      = 1.5      ; Exponent for r dependence of p
exppz      = 2.       ; Exponent for z dependence of p
ztrop      = 15000.   ; Tropopause Height
qtrop      = 1.e-11     ; Tropopause specific humidity
rfpi       = 1000000. ; Radius within which to use fixed-point iter.
constTv    = 0.608    ; Constant for Virtual Temp Conversion
deltaz     = 2.e-5     ; Small number to ensure convergence in FPI
epsilon    = 1.e-25     ; Small number to aviod dividing by zero in wind calc
exponent = Rd*gamma_/g   ; exponent
T0    = Ts0*(1.+constTv*q0)   ; Surface temp
Ttrop = T0 - gamma_*ztrop     ; Tropopause temp
ptrop = p00*(Ttrop/T0)^(1./exponent); Tropopause pressure
      
f  = 2.*omega*sin(cen_lat/convert)           ; Coriolis parameter

; Great circle calculations
gr = a*acos(sin(cen_lat/convert)*sin(lat/convert) + \
           (cos(cen_lat/convert)*cos(lat/convert)*cos((lon-cen_lon)/convert)))          
;gr = gc_latlon(cen_lat,cen_lon,lat,lon,2,3)

; Define surface pressure as pressure in instead of constant  
p00 = tofloat(psin)
      
;-----------------------------------------------------------------------
;    initialize height field if provided pressure or pressure if provided z
;-----------------------------------------------------------------------

if (zcoords .eq. 1) then
; We are reading in Z coords and need to use hydrostatic to back out pressure

   height = z

   if (height .gt. ztrop) then
      p = ptrop*exp(-(g*(height-ztrop))/(Rd*Ttrop))
   else
      p = (p00-dp*exp(-(gr/rp)^exppr)*exp(-(height/zp)^exppz)) \
          * ((T0-gamma_*height)/T0)^(1/exponent)
   end if

else
; Reading in p coords



ps = -dp*exp(-1*(gr/rp)^exppr) + p00
      
;print("ps "+ps+"   p00 "+p00)

height = (T0/gamma_)*(1.-(p/ps)^exponent)

; If inside a certain distance of the center of the storm
; perform a Fixed-point iteration to calculate the height
; more accurately

     if (gr .lt. rfpi ) then
        ;print("We need to correct")
          zhere = height
          zherebefore = height
          
          do n = 1,1000
          
            fpidFdz =-exppz*zhere*dp*exp(-1*(gr/rp)^exppr)*exp(-1*(zhere/zp)^exppz)/(zp*zp)*((T0-gamma_*zhere)/T0)^(g/(Rd*gamma_)) \
            +g/(Rd*T0)*(p00-dp*exp(-1*(gr/rp)^exppr)*exp(-1*(zhere/zp)^exppz))*((T0-gamma_*zhere)/T0)^(g/(Rd*gamma_)-1.)
            
            fpiF = p-(p00-dp*exp(-1*(gr/rp)^exppr)*exp(-1*(zhere/zp)^exppz)) \
               *((T0-gamma_*zhere)/T0)^(g/(Rd*gamma_))
               
            zn = zhere - fpiF/fpidFdz
            
            if ( abs(zn-zhere)/abs(zn) .lt. deltaz) then
                ;print("We are done with "+n+" iterations")
                break
            end if
            
            zhere = zn
            
            if (n .eq. 1000) then
              print("FPI did not converge after 500 interations in q & T!!!")
            end if
          end do
          
          height = zn
   end if

end if
      
;-----------------------------------------------------------------------
;    initialize PS (surface pressure of moist air)
;-----------------------------------------------------------------------

      
;-----------------------------------------------------------------------
;    initialize U and V (wind components)
;-----------------------------------------------------------------------
      
      d1 = sin(cen_lat/convert)*cos(lat/convert) - cos(cen_lat/convert)*sin(lat/convert)*cos((lon/convert)-(cen_lon/convert))
      d2 = cos(cen_lat/convert)*sin((lon/convert)-(cen_lon/convert))
      d  = max( (/ epsilon, sqrt(d1^2. + d2^2.) /) )

      ufac = d1/d
      vfac = d2/d
      if (height .gt. ztrop) then
          u = 0.
          v = 0.
      else
          vt = (-f*gr/2+sqrt((f*gr/2)^(2) - (exppr*(gr/rp)^exppr)*(Rd*(T0-gamma_*height)) / (exppz*height*Rd*(T0-gamma_*height)/(g*zp^exppz) + 1. - p00/dp*exp((gr/rp)^exppr)*exp((height/zp)^exppz))))
          ;term1 = -f*gr/2
          ;term2 = f*f*gr*gr/4
          ;term3 = 1.5*(gr/rp)^(1.5)
          ;term4 = Rd*(T0-gamma_*height)
          ;term5 = 2*Rd*(T0-gamma_*height)*height/g/zp/zp
          ;term6 = (p00/dp)*exp((gr/rp)^(1.5))*exp((height/zp)^2)
          ;vt = term1 + sqrt(term2 - ((term3*term4)/(1+term5-term6)))
          ;print(term1)
          ;print(sqrt(term2 - ((term3*term4)/(1+term5-term6))))          
          v = vin + vfac*vt
          u = uin + ufac*vt
      end if
      
      ;-----------------------------------------------------------------------
      ;    set the vertical velocity to zero (only required for non-hydrostatic models)
      ;-----------------------------------------------------------------------
      
      w = 0.
      
      ;-----------------------------------------------------------------------
      ;    tracer q (specific humidity)
      ;-----------------------------------------------------------------------
      
      q = qin
      ; not modifying q
      ;if (height .gt. ztrop) then
      ;    q = qtrop
      ;else
      ;    q = q0*exp(-height/zq1)*exp(-1*(height/zq2)^exppz)
      ;end if
          
      ;-----------------------------------------------------------------------
      ;    initialize T (temperature)
      ;-----------------------------------------------------------------------
      
      ;if (height .gt. ztrop) then
      ;    t = Ttrop
      ;else
      ;    t = (T0-gamma_*height)/(1.+constTv*q)/(1.+exppz*Rd*(T0-gamma_*height)*height \
      ;        /(g*zp^exppz*(1.-p00/dp*exp((gr/rp)^exppr)*exp((height/zp)^exppz))))
      ;end if
      
      ; Uncomment this to debug with a uniform T
      ;Tin = 300.
      
      Tvin = Tin*(1+constTv*qin)
      
          t = (Tvin)/(1.+constTv*qin)/(1.+exppz*Rd*(Tvin)*height \
              /(g*zp^exppz*(1.-p00/dp*exp((gr/rp)^exppr)*exp((height/zp)^exppz))))

      
      ;-----------------------------------------------------------------------
      ;    initialize PHIS (surface geopotential)
      ;-----------------------------------------------------------------------
      
      phis = 0.  ; constant
      
      ;-----------------------------------------------------------------------
      ;    initialize RHO (density of moist air)
      ;-----------------------------------------------------------------------
      
      rho = p/(Rd*t*(1.+constTv*q))
      
      output(0) = tofloat(v)
      output(1) = tofloat(u)
      output(2) = tofloat(q)
      output(3) = tofloat(t)
      output(4) = tofloat(ps)
          
return(output)

end







function get_rp_from_dp_rmw(cen_lat:numeric,dp:numeric,target_rmw:numeric)

begin

;dp         = 5800.    ; Delta P for calculation of PS
;target_rmw = 50000.
;cen_lat    = 20.      ; Center latitude of initial vortex

rpi     = 200000.  ; First guess for calculation of PS
step    = 5000.    ; Initial step size (in m) for moving rpi
maxiter = 300      ; max number of iterations
err_stop = 111.    ; stop iterating if error between subsequent calls < err_stop

; Constants
a	= 6371220.	        ; Earth's Radius (m)
Rd 	= 287.0	          ; Ideal gas const dry air (J kg^-1 K^1)
g	= 9.80616	          ; Gravity (m s^2)
omega   = 7.292115d-5 ; angular velocity 1/s
pi      = 4.*atan(1.)       ; pi
convert = 180./pi   ; conversion factor: radians to degrees
zp         = 7000.    ; Height for calculation of P
q0         = 0.021    ; q at surface from Jordan
gamma_      = 0.007    ; lapse rate
Ts0        = 302.0   ; Surface temperature (SST)
p00        = 101500.  ; global mean surface pressure
p0         = 100000.  ; p for model level calculation
;cen_lat    = 20.      ; Center latitude of initial vortex
;cen_lon    = -40.     ; Center longitufe of initial vortex
zq1        = 3000.    ; Height 1 for q calculation
zq2        = 8000.    ; Height 2 for q calculation
exppr      = 1.5      ; Exponent for r dependence of p
exppz      = 2.       ; Exponent for z dependence of p
ztrop      = 15000.   ; Tropopause Height
qtrop      = 1.e-11     ; Tropopause specific humidity
rfpi       = 1000000. ; Radius within which to use fixed-point iter.
constTv    = 0.608    ; Constant for Virtual Temp Conversion
deltaz     = 2.e-5     ; Small number to ensure convergence in FPI
epsilon    = 1.e-25     ; Small number to aviod dividing by zero in wind calc
exponent = Rd*gamma_/g   ; exponent
T0    = Ts0*(1.+constTv*q0)   ; Surface temp
Ttrop = T0 - gamma_*ztrop     ; Tropopause temp
ptrop = p00*(Ttrop/T0)^(1./exponent); Tropopause pressure

; We can tweak these!

; calculate coriolis param
fc  = 2.*omega*sin(cen_lat/convert)           ; Coriolis parameter

do jj = 0,maxiter-1

  rp = rpi

  ; define initial arrays
  rr = fspan(0.,1000000.,10000+1)
  vt = rr
  vt = 0.

  ; find RMW
  do ii = 0,dimsizes(rr)-1
    r = rr(ii)
    T1 = - (fc * r) / 2.
    T2 =  (fc*fc*r*r / 4.)
    NUM = (3./2.)* ((r/rp)^(3./2.)) * T0 * Rd
    DEN = 1. - (p00/dp)*exp((r/rp)^(3./2.))
    vt(ii) = tofloat ( T1 + sqrt( T2 - (NUM / DEN) ) )
  end do

  vmax = max(vt)
  rmw = rr(maxind(vt))

  err_here = rmw - target_rmw
  
  print("iter: "+jj+"  rpi: "+rpi+"  rmw: "+rmw+"  target: "+target_rmw+"  err: "+abs(err_here))

  ; If we've gotten enough convergence, break out of loop, we're done
  if (abs(err_here) .lt. err_stop) then
    break
  end if
  
  ; check if error last step and error this step is different sign, constrict step
  if (jj .ne. 0) then
    if (erri * err_here .lt. 0) then
      step = step/2.
    end if
  end if
  
  ; if not break, update rpi
  if (err_here .ge. 0) then
    ; if predicted RMW > target RMW, make rp smaller
    rpi = rpi - step
  else
    ; otherwise, make rp bigger
    rpi = rpi + step
  end if
  
  ; update error
  erri = err_here
  
  if (jj .eq. maxiter-1)
    print("did not converge!")
    exit
  end if
  
end do

print("dp="+dp)
print("rp="+rpi)

;wks   = gsn_open_wks ("x11","xy")               ; send graphics to PNG file
;res                  = True                     ; plot mods desired
;res@tiMainString     = "Basic XY plot"          ; add title
;plot  = gsn_csm_xy (wks,rr,vt,res) ; create plot

return(rpi)

end