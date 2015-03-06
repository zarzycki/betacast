load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"

begin

if (machineid .eq. 1) then
  filtfile_dir = "/glade/scratch/zarzycki/"+gridname+"/run/"
  writefile_dir = "/glade/p/work/zarzycki/sewx/INIC/"
  writefile_name = gridname+"_INIC_filter.nc"
else if (machineid .eq. 2) then
  print("Flux not configured yet")
  exit
else if (machineid .eq. 3) then
  filtfile_dir = "/home/zarzycki/"+gridname+"/run/"
  writefile_dir = "/home/zarzycki/sewx/INIC/"
  writefile_name = gridname+"_INIC_filter.nc"
else
  print("Invalid Machine ID specified")
  exit
end if
end if
end if

endhour = 6       ; End hour of run

filtfile = addfile(filtfile_dir+filtfile_name,"r")  
writefile = addfile(writefile_dir+writefile_name,"w")

PS = filtfile->PS
T = filtfile->T
U = filtfile->U
V = filtfile->V
Q = filtfile->Q
CLDLIQ = filtfile->CLDLIQ
CLDICE = filtfile->CLDICE

; Constants
pi = 3.14159
eps = 10^-16

dimPS = dimsizes(PS)
ncol = dimPS(1)
nlev = 30
time = filtfile->time
numtime = dimsizes(time)



h_k = new(numtime,double)
;theta_c = new(numtime,double)
omega = new(numtime,double)


M = numtime-1
N = M / 2
deltat = tofloat(endhour)/tofloat(M)  ; Timestep (hours)

; User-defined window cutoff value (hours)
tcut = 6

tn = fspan(0, endhour, numtime)
k = ispan(-N, N, 1)

theta_c = 2*pi*deltat/tcut
eps = 0

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

PS_filt = new((/1,ncol/),double)
T_filt = new((/1,nlev,ncol/),double)
U_filt = new((/1,nlev,ncol/),double)
V_filt = new((/1,nlev,ncol/),double)
Q_filt = new((/1,nlev,ncol/),double)
CLDLIQ_filt = new((/1,nlev,ncol/),double)
CLDICE_filt = new((/1,nlev,ncol/),double)

i = 0
do i = 0,ncol-1
  PS_filt(0,i) = PS(:,i)#h_k(:)
  do j = 0,nlev-1
    T_filt(0,j,i) = T(:,j,i)#h_k(:)
    U_filt(0,j,i) = U(:,j,i)#h_k(:)
    V_filt(0,j,i) = V(:,j,i)#h_k(:)
    Q_filt(0,j,i) = Q(:,j,i)#h_k(:)
    CLDLIQ_filt(0,j,i) = CLDLIQ(:,j,i)#h_k(:)
    CLDICE_filt(0,j,i) = CLDICE(:,j,i)#h_k(:)
  end do
end do

print("Writing file...")

  writefile->PS=(/PS_filt/)
  writefile->U=(/U_filt/)
  writefile->V=(/V_filt/)
  writefile->T=(/T_filt/)
  writefile->Q=(/Q_filt/)
  writefile->CLDLIQ=(/CLDLIQ_filt/)
  writefile->CLDICE=(/CLDICE_filt/)

  print("done")
  
end
