begin

; User-defined window cutoff value (hours)
tcut = 6

print("FILTER PROCESSING")
print("Using filter file: "+filtfile_name)
print("Will output post-filter file here: "+writefile_name)


filtfile = addfile(filtfile_name,"r")
writefile = addfile(writefile_name,"w")
vars = (/"PS","T","U","V","Q","CLDLIQ","CLDICE"/)

; Constants
pi = 3.14159
eps = 10^-16

dimT = dimsizes(filtfile->T)
nlev = dimT(1)
ncol = dimT(2)
time = filtfile->time
numtime = dimsizes(time)

h_k = new(numtime,double)
;theta_c = new(numtime,double)
omega = new(numtime,double)

M = numtime-1
N = M / 2
deltat = tofloat(endhour)/tofloat(M)  ; Timestep (hours)

tn = fspan(0, endhour, numtime)
k = ispan(-N, N, 1)

theta_c = 2*pi*deltat/tcut
eps = 0

print("-----------------")
print("Using a t_cut of: "+tcut)
print("Simulation is "+endhour+" hours long")
print("timestep is: "+deltat)

do i = 0,(dimsizes(k)-1)
  if k(i) .ne. 0 then
    omega(i) = sin(((eps+k(i))*pi)/(N+1))/(((eps+k(i))*pi)/(N+1)+eps)
    h_k(i) = (omega(i))*sin((eps+k(i))*theta_c)/((eps+k(i))*pi)
  else
    k_plus = 0.0001
    k_minus = -0.0001
    omega_plus = sin((k_plus*pi)/(N+1))/((k_plus*pi)/(N+1))
    h_k_plus = omega_plus*sin(k_plus*theta_c)/(k_plus*pi)
    omega_minus = sin((k_minus*pi)/(N+1))/((k_minus*pi)/(N+1))
    h_k_minus = omega_minus*sin(k_minus*theta_c)/(k_minus*pi)
    h_k(i) = (h_k_plus+h_k_minus)/2
  end if
end do

h_k = h_k * (1/sum(h_k))

print(h_k)

;print(tn)
;print(k)

;do i = 0,numtime-1
;  omega_k(i) = sin( (eps+k(i)*pi)/(M+1) ) / ( (eps+k(i)*pi)/(M+1) )
;  theta_c(i) = 2*pi*deltat/endhour
;  h_k(i) = omega_k(i) * ((sin (k(i)*theta_c(i))) / (eps+k(i)*pi))
;
;  omega_k(i) = sin( (eps+i*pi)/(N+1) ) / ( (eps+i*pi)/(N+1) )
;  theta_c(i) = 2*pi*deltat/endhour
;  h_k(i) = sin(i*theta_c(i))/(eps+i*pi)
;end do

; 6 hour cut
;h_k2 = (/-0.0000,-0.0025,-0.0076,-0.0134,-0.0167,-0.0133,0.0000,0.0246,0.0587,0.0970,0.1325, \
;    0.1576,0.1660,0.1576,0.1325,0.0970,0.0587,0.0246,0.0000,-0.0133,-0.0167,-0.0134,-0.0076,-0.0025,-0.0000/)

; 8 hour cut
;h_k = (/  -0.002137727887500,-0.004564236001700,-0.006031462869418,-0.005001378833234, \
;   0.000000000000000,0.009972620062365,0.025081193249744,0.044446573057312,0.066162990612804, \
;   0.087560955657867,0.105672778704921,0.117802518544982,0.122070351403714,0.117802518544982,0.105672778704921, \
;   0.087560955657867,0.066162990612804,0.044446573057312,0.025081193249744,0.009972620062365,0.000000000000000, \
;  -0.005001378833234,-0.006031462869418,-0.004564236001700,-0.002137727887500/)

; 12 hour cut

; h_k =(/   0.000000000000000, \
;    0.001397206131236, \
;    0.004660356786507, \
;    0.010098261183084, \
;    0.017783822199420, \
;    0.027505839019697, \
;    0.038759190500140, \
;    0.050778307285637, \
;    0.062611955069396, \
;    0.073230441921951, \
;    0.081650647952025, \
;    0.087060799548083, \
;    0.088926344805649, \
;    0.087060799548083, \
;    0.081650647952025, \
;    0.073230441921951, \
;    0.062611955069396, \
;    0.050778307285637, \
;    0.038759190500140, \
;    0.027505839019697, \
;    0.017783822199420, \
;    0.010098261183084, \
;    0.004660356786507, \
;    0.001397206131236, \
;    0.000000000000000/)

do z = 0,dimsizes(vars)-1

thisVar = vars(z)
print("Doing "+thisVar)

x = filtfile->$thisVar$

if (thisVar .eq. "PS") then
  x_filt = new((/1,ncol/),typeof(x))
else
  x_filt = new((/1,nlev,ncol/),typeof(x))
end if

i = 0
if (thisVar .eq. "PS") then
  do i = 0,ncol-1
    x_filt(0,i) = x(:,i)#totype(h_k(:),typeof(x_filt))
  end do
else
  do i = 0,ncol-1
    do j = 0,nlev-1
      x_filt(0,j,i) = x(:,j,i)#totype(h_k(:),typeof(x_filt))
    end do
  end do
end if

print("Writing "+thisVar+"...")

writefile->$thisVar$=(/x_filt/)

delete(x_filt)
delete(x)
delete(thisVar)
;delete(i)
;delete(j)

end do

print("done")

status_exit(9)

end
