
clear all
close all
clc

%--------------------------------------------------------------------------
% CUSTOM PARAMETERS (can change below)

%--- Plate Parameters
% physical and elastic parameters
Lx      = 0.151 ;             %-- length along x [m]
Ly      = 0.08 ;             %-- length along y [m]
Lz      = 0.81e-3 ;           %-- thickness [m]
E       = 101e9 ;             %-- Young's mod [Pa]
rho     = 8765 ;              %-- density [kg/m^3]
nu      = 0.3 ;               %-- poisson's ratio

%-- damping parameters: T60 = 3*log(10)/(sig0+om^2*sig1)
sig0    = 0.5 ;
sig1    = 3e-9 ;

%-- input / output locations, FRACTIONS of [Lx Ly] (values must be >0 and <1) ;
in      = [0.54,0.78] ;
outL    = [0.57,0.75] ;
outR    = [0.56,0.65] ;

% elastic constants around the edges (this allows to set the various bcs)
BCs = [0, 0;
    1e15, 1e15;
    0, 0;
    0, 0];

% largest frequency required
% -- Note: for larger plates, a large maxFreq may require solving a huge eigenvalue problem ...
%--  best to test this out using small values, such as maxFreq = 2000, and to increase it slowly to assess performance
maxFreq   = 15000 ;


%--- Simulation Parameters
T        = 6 ;               % sim length [s]
fs       = 44100 ;           % sample rate [Hz]
AmpF     = 30 ;              % force amplitude [N]
twid     = 0.0006 ;          % forcing temporal width [s]
%----------------------------------

%--------------------------------------------------------------------------



%--------------------------------------------------------------------------
% DERIVED PARAMETERS (DON'T TOUCH BELOW)

Ts           = round(T*fs) ;
k            = 1 / fs ;
tv           = (0:Ts-1)*k ;       %-- time axis array
fv           = (0:Ts-1)*fs/Ts ;   %-- freq axis array

%-- eigenvalue problem
ldim                = [Lx Ly Lz] ;
D                   = E * Lz^3 / 12 / (1-nu^2);
h                   =  sqrt(sqrt(D/rho/Lz*16/(maxFreq*2*pi)^2)) ; %-- set according to largest freq
[Om,Q,Nx,Ny,biHarm,Dm] = magpie(rho,E,nu,ldim,h,BCs,[],'none') ;
fOm = 0 ;
Nmodes   = 0 ;
OmDsq    = 0 ; 
while fOm < maxFreq && Nmodes < (Nx+1)*(Ny+1) && OmDsq >= 0
    Nmodes = Nmodes + 1 ;
    fOm = Om(Nmodes)/2/pi;
    C   = sig0 + sig1*Om(Nmodes)^2 ;
    OmDsq = Om(Nmodes)^2 - C^2 ;
end

Nmodes = Nmodes - 1 ; 
fMax   = Om(Nmodes)/2/pi  %-- check if this is in the range of maxFreq

Om   = Om(1:Nmodes) ;
Q    = Q(:,1:Nmodes) ;
C    = sig0 + sig1*Om.^2 ;
OmD  = sqrt(Om.^2-C.^2) ;



%-- build input vector (spreading, lin interp)

Jin      = zeros((Nx+1)*(Ny+1),1) ;

tempinx  = in(1)*Nx ;
Min      = floor(tempinx) ;
alx      = tempinx - Min ;

tempiny  = in(2)*Ny ;
Nin      = floor(tempiny) ;
aly      = tempiny - Nin ;

Jin((Ny+1)*Min+Nin+1)      = alx*aly ;
Jin((Ny+1)*(Min+1)+Nin+1)  = (1-alx)*aly ;
Jin((Ny+1)*Min+Nin+2)      = alx*(1-aly) ;
Jin((Ny+1)*(Min+1)+Nin+2)  = (1-alx)*(1-aly) ;


Jin = Jin / h^2 / rho / Lz ;
Jin = pinv(Q)*Jin ;

%-- left output weights (in interp)
JoutL     = zeros(1,(Nx+1)*(Ny+1)) ;

tempoutx  = outL(1)*Nx ;
Mout      = floor(tempoutx) ;
alx       = tempoutx - Mout ;

tempouty  = outL(2)*Ny ;
Nout      = floor(tempouty) ;
aly       = tempouty - Nout ;

JoutL((Ny+1)*Mout+Nout+1)      = alx*aly ;
JoutL((Ny+1)*(Mout+1)+Nout+1)  = (1-alx)*aly ;
JoutL((Ny+1)*Mout+Nout+2)      = alx*(1-aly) ;
JoutL((Ny+1)*(Mout+1)+Nout+2)  = (1-alx)*(1-aly) ;

%-- right output weights (in interp)
JoutR     = zeros(1,(Nx+1)*(Ny+1)) ;

tempoutx  = outR(1)*Nx ;
Mout      = floor(tempoutx) ;
alx       = tempoutx - Mout ;

tempouty  = outR(2)*Ny ;
Nout      = floor(tempouty) ;
aly       = tempouty - Nout ;

JoutR((Ny+1)*Mout+Nout+1)      = alx*aly ;
JoutR((Ny+1)*(Mout+1)+Nout+1)  = (1-alx)*aly ;
JoutR((Ny+1)*Mout+Nout+2)      = alx*(1-aly) ;
JoutR((Ny+1)*(Mout+1)+Nout+2)  = (1-alx)*(1-aly) ;

%-- input forcing
Nfin  = floor(twid*fs) ;
fin   = zeros(Ts,1) ;
fin(1:Nfin) = 0.5*AmpF*(1 - cos(2*pi*(0:Nfin-1)/Nfin)) ;


%--- init
vm   = zeros(Nmodes,1) ;
v0   = zeros(Nmodes,1) ;
outL = zeros(Ts,1) ;
outR = zeros(Ts,1) ;
velL = zeros(Ts,1) ;
velR = zeros(Ts,1) ;
outLprev = 0 ;
outRprev = 0 ;

%----------------------------------
%-- main loop
for n = 1 : Ts

    vp      = 2*exp(-C*k).*cos(OmD*k).*v0 - exp(-2*C*k).*vm + k^2*Jin*fin(n) ;
    outLcur = JoutL * (Q*v0) ;
    outRcur = JoutR * (Q*v0) ;
    outL(n) = outLcur ;
    outR(n) = outRcur ;
    velL(n) = (outLcur-outLprev)/k ;
    velR(n) = (outRcur-outRprev)/k ;
    vm      = v0 ;
    v0      = vp ;
    outLprev = outLcur ; outRprev = outRcur ;

end

out = [outL,outR] ;
vel = [velL,velR] ;
%----------------------------------

soundsc(vel,fs)

subplot(2,1,1)
plot(tv,out/1e-3);
xlabel('t (s)') ; ylabel('w (mm)'); title('output displacement'); legend('left channel','right channel')
subplot(2,1,2)
plot(tv,vel)
xlabel('t (s)') ; ylabel('dw/dt (m/s)'); title('output velocity'); legend('left channel','right channel')


figure
subplot(2,1,1)
plot(fv,20*log10(abs(fft(out))));
xlabel('f (Hz)') ; ylabel('$|\hat w|$ (dB)','interpreter','latex'); title('fft displacement'); legend('left channel','right channel')
xlim([0,5000]) ;
subplot(2,1,2)
plot(fv,20*log10(abs(fft(vel))));
xlabel('f (Hz)') ; ylabel('$|\hat\frac{dw}{dt}|$ (dB)','interpreter','latex'); title('fft velocity'); legend('left channel','right channel')
xlim([0,5000]) ;
