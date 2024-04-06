def D11_coeffs(R0y,Rx0,h,D,nu):
    return [(20 - (2*D - Rx0*h)/(2*D + Rx0*h) - (2*D - R0y*h)/(2*D + R0y*h)), - 8, 1,
    ((4*D + 4*D*nu)/(2*D + Rx0*h) - 8), ((4*D + 4*D*nu)/(2*D + R0y*h) - 8), - 8,
    1, 2, (2 - (2*D*nu)/(2*D + Rx0*h)), (2 - (2*D*nu)/(2*D + Rx0*h) - (2*D*nu)/(2*D + R0y*h)),
    (2 - (2*D*nu)/(2*D + R0y*h))];
