function [u12_c,u13_c,u14_c,u11_c,u10_c,u02_c,u22_c,u32_c,u23_c,u21_c,u01_c,u03_c] = D12_coeffs(R0y,h,D,nu)

u12_c = (20 - (2*D - R0y*h)/(2*D + R0y*h)) ;

u13_c = - 8 ;

u14_c = 1 ;

u11_c = - 8 ;

u10_c = 1 ;

u02_c = ((4*D + 4*D*nu)/(2*D + R0y*h) - 8) ;

u22_c = - 8 ;

u32_c = 1 ;

u23_c = 2 ;

u21_c = 2 ;

u01_c = (2 - (2*D*nu)/(2*D + R0y*h)) ;

u03_c = (2 - (2*D*nu)/(2*D + R0y*h)) ;

end

