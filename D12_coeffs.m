function varargout = D12_coeffs(R0y,h,D,nu)
    out = [(20 - (2*D - R0y*h)/(2*D + R0y*h)), - 8, 1, - 8, 1, ((4*D + 4*D*nu)/(2*D + R0y*h) - 8), - 8, 1, 2, 2, (2 - (2*D*nu)/(2*D + R0y*h)), (2 - (2*D*nu)/(2*D + R0y*h))];
    varargout = num2cell(out);
end

