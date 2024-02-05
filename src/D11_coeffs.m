function [u11_c,u12_c,u13_c,u10_c,u01_c,u21_c,u31_c,u22_c,u20_c,u00_c,u02_c] = D11_coeffs(R0y,Rx0,h,D,nu)

u11_c = (20 - (2*D - Rx0*h)/(2*D + Rx0*h) - (2*D - R0y*h)/(2*D + R0y*h)) ;

u12_c = - 8 ;

u13_c =  1 ;

u10_c = ((4*D + 4*D*nu)/(2*D + Rx0*h) - 8) ;

u01_c = ((4*D + 4*D*nu)/(2*D + R0y*h) - 8) ;

u21_c = - 8 ;

u31_c = 1 ;

u22_c = 2 ;

u20_c = (2 - (2*D*nu)/(2*D + Rx0*h)) ;

u00_c = (2 - (2*D*nu)/(2*D + Rx0*h) - (2*D*nu)/(2*D + R0y*h)) ;

u02_c = (2 - (2*D*nu)/(2*D + R0y*h)) ;

end
