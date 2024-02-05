function [u21_c,u22_c,u23_c,u20_c,u11_c,u31_c,u41_c,u01_c,u32_c,u30_c,u10_c,u12_c] = D21_coeffs(Rx0,h,D,nu)

u21_c = (20 - (2*D - Rx0*h)/(2*D + Rx0*h)) ;

u22_c = - 8 ;

u23_c = 1 ;

u20_c = + ((4*D + 4*D*nu)/(2*D + Rx0*h) - 8) ;

u11_c = - 8 ;

u31_c = - 8 ;

u41_c = 1 ;

u01_c = 1 ;

u32_c = 2 ;

u30_c = (2 - (2*D*nu)/(2*D + Rx0*h)) ;

u10_c = (2 - (2*D*nu)/(2*D + Rx0*h)) ;

u12_c = 2 ;

end




