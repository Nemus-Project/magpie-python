function [u01_c,u11_c,u21_c,u00_c,u02_c,u03_c,u12_c,u10_c] = D01_coeffs(K0y,R0y,Rx0,h,D,nu)

u01_c  = ((4*D^2*nu^2 - 4*D^2 - 2*D*R0y*h + 2*D*Rx0*h + R0y*Rx0*h^2)/(4*D^2 - 4*D^2*nu^2 + 2*D*R0y*h + 2*D*Rx0*h + R0y*Rx0*h^2) - (8*(4*D + 4*D*nu))/(2*D + R0y*h) - (4*D*nu)/(2*D + R0y*h) + (2*K0y*Rx0*R0y^2*h^6 + 4*K0y*D*R0y^2*h^5 + 8*K0y*Rx0*D*R0y*h^5 - 8*K0y*D^2*R0y*h^4*nu^2 + 16*K0y*D^2*R0y*h^4 - 14*Rx0*D^2*R0y*h^2*nu^2 + 28*Rx0*D^2*R0y*h^2*nu + 24*Rx0*D^2*R0y*h^2 - 20*D^3*R0y*h*nu^2 + 40*D^3*R0y*h*nu + 48*D^3*R0y*h + 8*K0y*Rx0*D^2*h^4 - 16*K0y*D^3*h^3*nu^2 + 16*K0y*D^3*h^3 - 28*Rx0*D^3*h*nu^2 + 56*Rx0*D^3*h*nu + 48*Rx0*D^3*h + 40*D^4*nu^4 - 80*D^4*nu^3 - 136*D^4*nu^2 + 80*D^4*nu + 96*D^4)/(D*(2*D + R0y*h)*(4*D^2 - 4*D^2*nu^2 + 2*D*R0y*h + 2*D*Rx0*h + R0y*Rx0*h^2)) - (8*D*Rx0*h*nu)/(4*D^2 - 4*D^2*nu^2 + 2*D*R0y*h + 2*D*Rx0*h + R0y*Rx0*h^2) + 20) ;

u11_c  = ((8*(2*D - R0y*h))/(2*D + R0y*h) + (32*D^4*nu - 96*D^4 + 96*D^4*nu^2 - 32*D^4*nu^3 - 48*D^3*R0y*h - 48*D^3*Rx0*h + 16*D^3*R0y*h*nu + 16*D^3*Rx0*h*nu - 24*D^2*R0y*Rx0*h^2 + 8*D^2*R0y*Rx0*h^2*nu)/(D*(2*D + R0y*h)*(4*D^2 - 4*D^2*nu^2 + 2*D*R0y*h + 2*D*Rx0*h + R0y*Rx0*h^2)) - 8) ;

u21_c  = ((Rx0*D*R0y^2*h^3 + 2*D^2*R0y^2*h^2 + 4*Rx0*D^2*R0y*h^2 - 4*D^3*R0y*h*nu^2 + 8*D^3*R0y*h + 4*Rx0*D^3*h - 8*D^4*nu^2 + 8*D^4)/(D*(2*D + R0y*h)*(4*D^2 - 4*D^2*nu^2 + 2*D*R0y*h + 2*D*Rx0*h + R0y*Rx0*h^2)) + 1) ;

u00_c  = ((- 8*D^2*nu^2 + 4*R0y*h*D*nu + 8*D^2 + 4*R0y*h*D)/(4*D^2 - 4*D^2*nu^2 + 2*D*R0y*h + 2*D*Rx0*h + R0y*Rx0*h^2) + (2*(- 8*D^2*nu^2 + 4*Rx0*h*D*nu + 8*D^2 + 4*Rx0*h*D))/(4*D^2 - 4*D^2*nu^2 + 2*D*R0y*h + 2*D*Rx0*h + R0y*Rx0*h^2) + (16*D*nu)/(2*D + R0y*h) - (32*D^4*nu + 32*D^4 - 48*D^4*nu^2 - 32*D^4*nu^3 + 16*D^4*nu^4 + 16*D^3*R0y*h + 16*D^3*Rx0*h + 16*D^3*R0y*h*nu + 32*D^3*Rx0*h*nu + 8*D^2*R0y*Rx0*h^2 - 24*D^3*R0y*h*nu^2 + 8*D^3*R0y*h*nu^3 - 16*D^3*Rx0*h*nu^2 - 8*D^2*R0y*Rx0*h^2*nu^2 + 16*D^2*R0y*Rx0*h^2*nu)/(D*(2*D + R0y*h)*(4*D^2 - 4*D^2*nu^2 + 2*D*R0y*h + 2*D*Rx0*h + R0y*Rx0*h^2)) - 8) ;

u02_c  = ((2*(4*D + 4*D*nu))/(2*D + R0y*h) + (16*D*nu)/(2*D + R0y*h) - (64*D^4*nu + 32*D^4 - 64*D^4*nu^2 - 64*D^4*nu^3 + 32*D^4*nu^4 + 16*D^3*R0y*h + 16*D^3*Rx0*h + 32*D^3*R0y*h*nu + 32*D^3*Rx0*h*nu + 8*D^2*R0y*Rx0*h^2 - 16*D^3*R0y*h*nu^2 - 16*D^3*Rx0*h*nu^2 - 8*D^2*R0y*Rx0*h^2*nu^2 + 16*D^2*R0y*Rx0*h^2*nu)/(D*(2*D + R0y*h)*(4*D^2 - 4*D^2*nu^2 + 2*D*R0y*h + 2*D*Rx0*h + R0y*Rx0*h^2)) - 8) ;

u03_c  = ((16*D^4*nu - 8*D^4*nu^2 - 16*D^4*nu^3 + 8*D^4*nu^4 + 8*D^3*R0y*h*nu + 8*D^3*Rx0*h*nu - 4*D^3*R0y*h*nu^2 - 4*D^3*Rx0*h*nu^2 - 2*D^2*R0y*Rx0*h^2*nu^2 + 4*D^2*R0y*Rx0*h^2*nu)/(D*(2*D + R0y*h)*(4*D^2 - 4*D^2*nu^2 + 2*D*R0y*h + 2*D*Rx0*h + R0y*Rx0*h^2)) - (4*D*nu)/(2*D + R0y*h) + 1) ;

u12_c  = (2 - (16*D^4*nu - 32*D^4 + 32*D^4*nu^2 - 16*D^4*nu^3 - 16*D^3*R0y*h - 16*D^3*Rx0*h + 8*D^3*R0y*h*nu + 8*D^3*Rx0*h*nu - 8*D^2*R0y*Rx0*h^2 + 4*D^2*R0y*Rx0*h^2*nu)/(D*(2*D + R0y*h)*(4*D^2 - 4*D^2*nu^2 + 2*D*R0y*h + 2*D*Rx0*h + R0y*Rx0*h^2)) - (2*(2*D - R0y*h))/(2*D + R0y*h)) ;

u10_c  = ((2*(4*D^2*nu^2 - 4*D^2 + 2*D*R0y*h - 2*D*Rx0*h + R0y*Rx0*h^2))/(4*D^2 - 4*D^2*nu^2 + 2*D*R0y*h + 2*D*Rx0*h + R0y*Rx0*h^2) - (16*D^4*nu - 32*D^4 + 32*D^4*nu^2 - 16*D^4*nu^3 - 16*D^3*R0y*h - 16*D^3*Rx0*h + 8*D^3*R0y*h*nu + 8*D^3*Rx0*h*nu - 8*D^2*R0y*Rx0*h^2 + 16*D^3*R0y*h*nu^2 - 8*D^3*R0y*h*nu^3 + 4*D^2*R0y*Rx0*h^2*nu)/(D*(2*D + R0y*h)*(4*D^2 - 4*D^2*nu^2 + 2*D*R0y*h + 2*D*Rx0*h + R0y*Rx0*h^2)) - (4*D*R0y*h*nu)/(4*D^2 - 4*D^2*nu^2 + 2*D*R0y*h + 2*D*Rx0*h + R0y*Rx0*h^2) + 2) ;

end