function [u20_c,u21_c,u22_c,u10_c,u30_c,u40_c,u00_c,u31_c,u11_c] = D20_coeffs(Kx0,Rx0,h,D,nu)

u20_c = ((2*Kx0*Rx0*h^4 + 4*Kx0*D*h^3 - 12*D^2*nu^2 + 24*D^2*nu + 24*D^2)/(D*(2*D + Rx0*h)) - (8*D*nu)/(2*D + Rx0*h) - (8*(4*D + 4*D*nu))/(2*D + Rx0*h) + 20) ;

u21_c = ((8*(2*D - Rx0*h))/(2*D + Rx0*h) + (8*D^2*nu - 24*D^2)/(D*(2*D + Rx0*h)) - 8) ;

u22_c = ((2*D^2 + Rx0*h*D)/(D*(2*D + Rx0*h)) + 1) ;

u10_c = ((2*(4*D + 4*D*nu))/(2*D + Rx0*h) - (- 8*D^2*nu^2 + 16*D^2*nu + 8*D^2)/(D*(2*D + Rx0*h)) + (16*D*nu)/(2*D + Rx0*h) - 8) ;
u30_c = u10_c;

u40_c = ((- 2*D^2*nu^2 + 4*D^2*nu)/(D*(2*D + Rx0*h)) - (4*D*nu)/(2*D + Rx0*h) + 1) ;
u00_c = u40_c;

u31_c = (2 - (4*D^2*nu - 8*D^2)/(D*(2*D + Rx0*h)) - (2*(2*D - Rx0*h))/(2*D + Rx0*h)) ;
u11_c = u31_c;

end

