function [Q,Om,Nx,Ny,biHarm] = biharmonicEigs(rho,E,nu,Lx,Ly,Lz,h,K0y,R0y,Kx0,Rx0,KLy,RLy,KxL,RxL,Nmodes)
% BIHARMONICEIGS
%
%

    %%--- derived parameters (don't change here)
    D       = E * Lz^3 / 12 / (1-nu^2) ;
    Nx      = floor(Lx/h) ;
    Ny      = floor(Ly/h) ;
    %%----------------------------


    %% MATRIX BUILDER

    %%--- build matrix in blocks

    a0 = ones(Ny+1,1) ;
    a1 = ones(Ny,1) ;
    a2 = ones(Ny-1,1) ;

    [D00u00,D00u10,D00u20,D00u01,D00u02,D00u11] = D00_coeffs(K0y,R0y,Kx0,Rx0,h,D,nu) ;
    [D01u01,D01u11,D01u21,D01u00,D01u02,D01u03,D01u12,D01u10] = D01_coeffs(K0y,R0y,Rx0,h,D,nu) ;
    [D02u02,D02u12,D02u22,D02u01,D02u03,D02u04,D02u00,D02u13,D02u11] = D02_coeffs(K0y,R0y,h,D,nu) ;
    [D0Nu0N,D0Nu1N,D0Nu2N,D0Nu0Nm1,D0Nu0Nm2,D0Nu1Nm1] = D00_coeffs(K0y,R0y,KxL,RxL,h,D,nu) ;
    [D0Nm1u0Nm1,D0Nm1u1Nm1,D0Nm1u2Nm1,D0Nm1u0N,D0Nm1u0Nm2,D0Nm1u0Nm3,D0Nm1u1Nm2,D0Nm1u1N] = D01_coeffs(K0y,R0y,RxL,h,D,nu) ;
    %%
    %%-- Blk11
    D0 = D02u02*a0 ; D1 = D02u03*a1; D2 = D02u04*a2 ; Dm1 = D02u01*a1; Dm2 = D02u00*a2 ;

    Blk11               = sparse(diag(D0) + diag(D1,1) + diag(Dm1,-1) + diag(D2,2) + diag(Dm2,-2)) ;
    Blk11(1,1)          = D00u00 ;
    Blk11(1,2)          = D00u01 ;
    Blk11(1,3)          = D00u02 ;
    Blk11(2,1)          = D01u00 ;
    Blk11(2,2)          = D01u01 ;
    Blk11(2,3)          = D01u02 ;
    Blk11(2,4)          = D01u03 ;
    Blk11(end,end)      = D0Nu0N ;
    Blk11(end,end-1)    = D0Nu0Nm1 ;
    Blk11(end,end-2)    = D0Nu0Nm2 ;
    Blk11(end-1,end)    = D0Nm1u0N ;
    Blk11(end-1,end-1)  = D0Nm1u0Nm1 ;
    Blk11(end-1,end-2)  = D0Nm1u0Nm2 ;
    Blk11(end-1,end-3)  = D0Nm1u0Nm3 ;
    
    %%
    %%%-- Blk12
    D0 = D02u12*a0 ; D1 = D02u13*a1 ; Dm1 = D02u11*a1 ;

    Blk12               = sparse(diag(D0) + diag(D1,1) + diag(Dm1,-1))  ;
    Blk12(1,1)          = D00u10 ;
    Blk12(1,2)          = D00u11 ;
    Blk12(2,1)          = D01u10 ;
    Blk12(2,2)          = D01u11 ;
    Blk12(2,3)          = D01u12 ;
    Blk12(end,end)      = D0Nu1N ;
    Blk12(end,end-1)    = D0Nu1Nm1 ;
    Blk12(end-1,end)    = D0Nm1u1N ;
    Blk12(end-1,end-1)  = D0Nm1u1Nm1 ;
    Blk12(end-1,end-2)  = D0Nm1u1Nm2 ;

    %%
    %%-- Blk13
    D0 = D02u22*a0 ;

    Blk13               = sparse(diag(D0))   ;
    Blk13(1,1)          = D00u20 ;
    Blk13(2,2)          = D01u21 ;
    Blk13(end,end)      = D0Nu2N ;
    Blk13(end-1,end-1)  = D0Nm1u2Nm1 ;

    [D10u10, D10u20, D10u30, D10u00, D10u11, D10u12, D10u21, D10u01]                                                     = D10_coeffs(R0y,Kx0,Rx0,h,D,nu) ;
    [D11u11,D11u12,D11u13,D11u10,D11u01,D11u21,D11u31,D11u22,D11u20,D11u00,D11u02]                                       = D11_coeffs(R0y,Rx0,h,D,nu) ;
    [D12u12,D12u13,D12u14,D12u11,D12u10,D12u02,D12u22,D12u32,D12u23,D12u21,D12u01,D12u03]                                = D12_coeffs(R0y,h,D,nu) ;
    [D1Nu1N, D1Nu2N, D1Nu3N, D1Nu0N, D1Nu1Nm1, D1Nu1Nm2, D1Nu2Nm1, D1Nu0Nm1]                                             = D10_coeffs(R0y,KxL,RxL,h,D,nu) ;
    [D1Nm1u1Nm1,D1Nm1u1Nm2,D1Nm1u1Nm3,D1Nm1u1N,D1Nm1u0Nm1,D1Nm1u2Nm1,D1Nm1u3Nm1,D1Nm1u2Nm2,D1Nm1u2N,D1Nm1u0N,D1Nm1u0Nm2] = D11_coeffs(R0y,RxL,h,D,nu) ;
    %%
    %%-- Blk21
    D0 = D12u02*a0 ; D1 = D12u03*a1 ; Dm1 = D12u01*a1 ;

    Blk21               = sparse(diag(D0) + diag(D1,1) + diag(Dm1,-1))  ;
    Blk21(1,1)          = D10u00 ;
    Blk21(1,2)          = D10u01 ;
    Blk21(2,1)          = D11u00 ;
    Blk21(2,2)          = D11u01 ;
    Blk21(2,3)          = D11u02 ;
    Blk21(end,end)      = D1Nu0N ;
    Blk21(end,end-1)    = D1Nu0Nm1 ;
    Blk21(end-1,end)    = D1Nm1u0N ;
    Blk21(end-1,end-1)  = D1Nm1u0Nm1 ;
    Blk21(end-1,end-2)  = D1Nm1u0Nm2 ;

    %%
    %%-- Blk22
    D0 = D12u12*a0 ; D1 = D12u13*a1 ; D2 = D12u14*a2 ; Dm1 = D12u11*a1 ; Dm2 = D12u10*a2 ;

    Blk22               = sparse(diag(D0) + diag(D1,1) + diag(Dm1,-1) + diag(D2,2) + diag(Dm2,-2)) ;
    Blk22(1,1)          = D10u10 ;
    Blk22(1,2)          = D10u11 ;
    Blk22(1,3)          = D10u12 ;
    Blk22(2,1)          = D11u10 ;
    Blk22(2,2)          = D11u11 ;
    Blk22(2,3)          = D11u12 ;
    Blk22(2,4)          = D11u13 ;

    Blk22(end,end)      = D1Nu1N ;
    Blk22(end,end-1)    = D1Nu1Nm1 ;
    Blk22(end,end-2)    = D1Nu1Nm2 ;
    Blk22(end-1,end)    = D1Nm1u1N ;
    Blk22(end-1,end-1)  = D1Nm1u1Nm1 ;
    Blk22(end-1,end-2)  = D1Nm1u1Nm2 ;
    Blk22(end-1,end-3)  = D1Nm1u1Nm3 ;

    %%
    %%-- Blk23
    D0 = D12u22*a0 ; D1 = D12u23*a1 ; Dm1 = D12u21*a1 ;

    Blk23               = sparse(diag(D0) + diag(D1,1) + diag(Dm1,-1)) ;
    Blk23(1,1)          = D10u20 ;
    Blk23(1,2)          = D10u21 ;
    Blk23(2,1)          = D11u20 ;
    Blk23(2,2)          = D11u21 ;
    Blk23(2,3)          = D11u22 ;

    Blk23(end,end)      = D1Nu2N ;
    Blk23(end,end-1)    = D1Nu2Nm1 ;
    Blk23(end-1,end)    = D1Nm1u2N ;
    Blk23(end-1,end-1)  = D1Nm1u2Nm1 ;
    Blk23(end-1,end-2)  = D1Nm1u2Nm2 ;

    %%
    %%-- Blk24
    D0 = D12u32*a0 ;

    Blk24               = sparse(diag(D0))   ;
    Blk24(1,1)          = D10u30 ;
    Blk24(2,2)          = D11u31 ;
    Blk24(end,end)      = D1Nu3N ;
    Blk24(end-1,end-1)  = D1Nm1u3Nm1 ;

    [D20u20,D20u21,D20u22,D20u10,D20u30,D20u40,D20u00,D20u31,D20u11]                                                                = D20_coeffs(Kx0,Rx0,h,D,nu) ;
    [D21u21,D21u22,D21u23,D21u20,D21u11,D21u31,D21u41,D21u01,D21u32,D21u30,D21u10,D21u12]                                           = D21_coeffs(Rx0,h,D,nu) ;
    [D22u20,D22u11,D22u21,D22u31,D22u02,D22u12,D22u22,D22u32,D22u42,D22u13,D22u23,D22u33,D22u24]                                    = D22_coeffs ;
    [D2Nu2N,D2Nu2Nm1,D2Nu2Nm2,D2Nu1N,D2Nu3N,D2Nu4N,D2Nu0N,D2Nu3Nm1,D2Nu1Nm1]                                                        = D20_coeffs(KxL,RxL,h,D,nu) ;
    [D2Nm1u2Nm1,D2Nm1u2Nm2,D2Nm1u2Nm3,D2Nm1u2N,D2Nm1u1Nm1,D2Nm1u3Nm1,D2Nm1u4Nm1,D2Nm1u0Nm1,D2Nm1u3Nm2,D2Nm1u3N,D2Nm1u1N,D2Nm1u1Nm2] = D21_coeffs(RxL,h,D,nu) ;

    %%
    %%-- Blk31
    D0 = D22u02*a0 ;

    Blk31               = sparse(diag(D0))   ;
    Blk31(1,1)          = D20u00 ;
    Blk31(2,2)          = D21u01 ;
    Blk31(end,end)      = D2Nu0N ;
    Blk31(end-1,end-1)  = D2Nm1u0Nm1 ;

    %%
    %%-- Blk32
    D0 = D22u12*a0 ; D1 = D22u13*a1 ; Dm1 = D22u11*a1 ;

    Blk32               = sparse(diag(D0) + diag(D1,1) + diag(Dm1,-1)) ;
    Blk32(1,1)          = D20u10 ;
    Blk32(1,2)          = D20u11 ;
    Blk32(2,1)          = D21u10 ;
    Blk32(2,2)          = D21u11 ;
    Blk32(2,3)          = D21u12 ;

    Blk32(end,end)      = D2Nu1N  ;
    Blk32(end,end-1)    = D2Nu1Nm1 ;
    Blk32(end-1,end)    = D2Nm1u1N ;
    Blk32(end-1,end-1)  = D2Nm1u1Nm1 ;
    Blk32(end-1,end-2)  = D2Nm1u1Nm2 ;

    %%
    %%-- Blk33
    D0 = D22u22*a0 ; D1 = D22u23*a1; D2 = D22u24*a2 ; Dm1 = D22u21*a1; Dm2 = D22u20*a2 ;

    Blk33               = sparse(diag(D0) + diag(D1,1) + diag(Dm1,-1) + diag(D2,2) + diag(Dm2,-2)) ;
    Blk33(1,1)          = D20u20 ;
    Blk33(1,2)          = D20u21 ;
    Blk33(1,3)          = D20u22 ;
    Blk33(2,1)          = D21u20 ;
    Blk33(2,2)          = D21u21 ;
    Blk33(2,3)          = D21u22 ;
    Blk33(2,4)          = D21u23 ;

    Blk33(end,end)      = D2Nu2N ;
    Blk33(end,end-1)    = D2Nu2Nm1 ;
    Blk33(end,end-2)    = D2Nu2Nm2 ;
    Blk33(end-1,end)    = D2Nm1u2N ;
    Blk33(end-1,end-1)  = D2Nm1u2Nm1 ;
    Blk33(end-1,end-2)  = D2Nm1u2Nm2 ;
    Blk33(end-1,end-3)  = D2Nm1u2Nm3 ;

    %%
    %%-- Blk34
    D0 = D22u32*a0 ; D1 = D22u33*a1 ; Dm1 = D22u31*a1 ;

    Blk34               = sparse(diag(D0) + diag(D1,1) + diag(Dm1,-1)) ;
    Blk34(1,1)          =  D20u30 ;
    Blk34(1,2)          =  D20u31 ;
    Blk34(2,1)          =  D21u30 ;
    Blk34(2,2)          =  D21u31 ;
    Blk34(2,3)          =  D21u32 ;

    Blk34(end,end)      =  D2Nu3N;
    Blk34(end,end-1)    =  D2Nu3Nm1;
    Blk34(end-1,end)    =  D2Nm1u3N;
    Blk34(end-1,end-1)  =  D2Nm1u3Nm1;
    Blk34(end-1,end-2)  =  D2Nm1u3Nm2;

    %%
    %%-- Blk35
    D0 = D22u42*a0 ;

    Blk35               = sparse(diag(D0))   ;
    Blk35(1,1)          = D20u40 ;
    Blk35(2,2)          = D21u41 ;
    Blk35(end,end)      = D2Nu4N ;
    Blk35(end-1,end-1)  = D2Nm1u4Nm1 ;

    [D00u00,D00u10,D00u20,D00u01,D00u02,D00u11] = D00_coeffs(KLy,RLy,Kx0,Rx0,h,D,nu) ;
    [D01u01,D01u11,D01u21,D01u00,D01u02,D01u03,D01u12,D01u10] = D01_coeffs(KLy,RLy,Rx0,h,D,nu) ;
    [D02u02,D02u12,D02u22,D02u01,D02u03,D02u04,D02u00,D02u13,D02u11] = D02_coeffs(KLy,RLy,h,D,nu) ;
    [D0Nu0N,D0Nu1N,D0Nu2N,D0Nu0Nm1,D0Nu0Nm2,D0Nu1Nm1] = D00_coeffs(KLy,RLy,KxL,RxL,h,D,nu) ;
    [D0Nm1u0Nm1,D0Nm1u1Nm1,D0Nm1u2Nm1,D0Nm1u0N,D0Nm1u0Nm2,D0Nm1u0Nm3,D0Nm1u1Nm2,D0Nm1u1N] = D01_coeffs(KLy,RLy,RxL,h,D,nu) ;

    %%
    %%-- BlkMM
    D0 = D02u02*a0 ; D1 = D02u03*a1; D2 = D02u04*a2 ; Dm1 = D02u01*a1; Dm2 = D02u00*a2 ;

    BlkMM               = sparse(diag(D0) + diag(D1,1) + diag(Dm1,-1) + diag(D2,2) + diag(Dm2,-2)) ;
    BlkMM(1,1)          = D00u00 ;
    BlkMM(1,2)          = D00u01 ;
    BlkMM(1,3)          = D00u02 ;
    BlkMM(2,1)          = D01u00 ;
    BlkMM(2,2)          = D01u01 ;
    BlkMM(2,3)          = D01u02 ;
    BlkMM(2,4)          = D01u03 ;
    BlkMM(end,end)      = D0Nu0N ;
    BlkMM(end,end-1)    = D0Nu0Nm1 ;
    BlkMM(end,end-2)    = D0Nu0Nm2 ;
    BlkMM(end-1,end)    = D0Nm1u0N ;
    BlkMM(end-1,end-1)  = D0Nm1u0Nm1 ;
    BlkMM(end-1,end-2)  = D0Nm1u0Nm2 ;
    BlkMM(end-1,end-3)  = D0Nm1u0Nm3 ;

    %%
    %%-- BlkMMm1
    D0 = D02u12*a0 ; D1 = D02u13*a1 ; Dm1 = D02u11*a1 ;

    BlkMMm1               = sparse(diag(D0) + diag(D1,1) + diag(Dm1,-1))  ;
    BlkMMm1(1,1)          = D00u10 ;
    BlkMMm1(1,2)          = D00u11 ;
    BlkMMm1(2,1)          = D01u10 ;
    BlkMMm1(2,2)          = D01u11 ;
    BlkMMm1(2,3)          = D01u12 ;
    BlkMMm1(end,end)      = D0Nu1N ;
    BlkMMm1(end,end-1)    = D0Nu1Nm1 ;
    BlkMMm1(end-1,end)    = D0Nm1u1N ;
    BlkMMm1(end-1,end-1)  = D0Nm1u1Nm1 ;
    BlkMMm1(end-1,end-2)  = D0Nm1u1Nm2 ;

    %%
    %%-- BlkMMm2
    D0 = D02u22*a0 ;

    BlkMMm2               = sparse(diag(D0))   ;
    BlkMMm2(1,1)          = D00u20 ;
    BlkMMm2(2,2)          = D01u21 ;
    BlkMMm2(end,end)      = D0Nu2N ;
    BlkMMm2(end-1,end-1)  = D0Nm1u2Nm1 ;

    [D10u10, D10u20, D10u30, D10u00, D10u11, D10u12, D10u21, D10u01]                                                              = D10_coeffs(RLy,Kx0,Rx0,h,D,nu) ;
    [D11u11,D11u12,D11u13,D11u10,D11u01,D11u21,D11u31,D11u22,D11u20,D11u00,D11u02]                                                = D11_coeffs(RLy,Rx0,h,D,nu) ;
    [D12u12,D12u13,D12u14,D12u11,D12u10,D12u02,D12u22,D12u32,D12u23,D12u21,D12u01,D12u03]                                         = D12_coeffs(RLy,h,D,nu) ;
    [D1Nu1N, D1Nu2N, D1Nu3N, D1Nu0N, D1Nu1Nm1, D1Nu1Nm2, D1Nu2Nm1, D1Nu0Nm1]                                                      = D10_coeffs(RLy,KxL,RxL,h,D,nu) ;
    [D1Nm1u1Nm1,D1Nm1u1Nm2,D1Nm1u1Nm3,D1Nm1u1N,D1Nm1u0Nm1,D1Nm1u2Nm1,D1Nm1u3Nm1,D1Nm1u2Nm2,D1Nm1u2N,D1Nm1u0N,D1Nm1u0Nm2]          = D11_coeffs(RLy,RxL,h,D,nu) ;

    %%
    %%-- BlkMm1M
    D0 = D12u02*a0 ; D1 = D12u03*a1 ; Dm1 = D12u01*a1 ;

    BlkMm1M               = sparse(diag(D0) + diag(D1,1) + diag(Dm1,-1))  ;
    BlkMm1M(1,1)          = D10u00 ;
    BlkMm1M(1,2)          = D10u01 ;
    BlkMm1M(2,1)          = D11u00 ;
    BlkMm1M(2,2)          = D11u01 ;
    BlkMm1M(2,3)          = D11u02 ;
    BlkMm1M(end,end)      = D1Nu0N ;
    BlkMm1M(end,end-1)    = D1Nu0Nm1 ;
    BlkMm1M(end-1,end)    = D1Nm1u0N ;
    BlkMm1M(end-1,end-1)  = D1Nm1u0Nm1 ;
    BlkMm1M(end-1,end-2)  = D1Nm1u0Nm2 ;

    %%
    %%-- BlkMm1Mm1
    D0 = D12u12*a0 ; D1 = D12u13*a1 ; D2 = D12u14*a2 ; Dm1 = D12u11*a1 ; Dm2 = D12u10*a2 ;


    BlkMm1Mm1               = sparse(diag(D0) + diag(D1,1) + diag(Dm1,-1) + diag(D2,2) + diag(Dm2,-2)) ;
    BlkMm1Mm1(1,1)          = D10u10 ;
    BlkMm1Mm1(1,2)          = D10u11 ;
    BlkMm1Mm1(1,3)          = D10u12 ;
    BlkMm1Mm1(2,1)          = D11u10 ;
    BlkMm1Mm1(2,2)          = D11u11 ;
    BlkMm1Mm1(2,3)          = D11u12 ;
    BlkMm1Mm1(2,4)          = D11u13 ;

    BlkMm1Mm1(end,end)      = D1Nu1N ;
    BlkMm1Mm1(end,end-1)    = D1Nu1Nm1 ;
    BlkMm1Mm1(end,end-2)    = D1Nu1Nm2 ;
    BlkMm1Mm1(end-1,end)    = D1Nm1u1N ;
    BlkMm1Mm1(end-1,end-1)  = D1Nm1u1Nm1 ;
    BlkMm1Mm1(end-1,end-2)  = D1Nm1u1Nm2 ;
    BlkMm1Mm1(end-1,end-3)  = D1Nm1u1Nm3 ;

    %%
    %%-- BlkMm1Mm2
    D0 = D12u22*a0 ; D1 = D12u23*a1 ; Dm1 = D12u21*a1 ;

    BlkMm1Mm2               = sparse(diag(D0) + diag(D1,1) + diag(Dm1,-1)) ;
    BlkMm1Mm2(1,1)          = D10u20 ;
    BlkMm1Mm2(1,2)          = D10u21 ;
    BlkMm1Mm2(2,1)          = D11u20 ;
    BlkMm1Mm2(2,2)          = D11u21 ;
    BlkMm1Mm2(2,3)          = D11u22 ;

    BlkMm1Mm2(end,end)      = D1Nu2N ;
    BlkMm1Mm2(end,end-1)    = D1Nu2Nm1 ;
    BlkMm1Mm2(end-1,end)    = D1Nm1u2N ;
    BlkMm1Mm2(end-1,end-1)  = D1Nm1u2Nm1 ;
    BlkMm1Mm2(end-1,end-2)  = D1Nm1u2Nm2 ;

    %%
    %%-- BlkMm1Mm3
    D0 = D12u32*a0 ;

    BlkMm1Mm3               = sparse(diag(D0))   ;
    BlkMm1Mm3(1,1)          = D10u30 ;
    BlkMm1Mm3(2,2)          = D11u31 ;
    BlkMm1Mm3(end,end)      = D1Nu3N ;
    BlkMm1Mm3(end-1,end-1)  = D1Nm1u3Nm1 ;


    biHarm = sparse((Nx+1)*(Ny+1),(Nx+1)*(Ny+1)) ;


    for m = 3 : Nx - 1
        biHarm((Ny+1)*(m-1)+1 : (Ny+1)*m, (Ny+1)*(m-1)+1 : (Ny+1)*m)        = Blk33 ;
        biHarm((Ny+1)*(m-1)+1 : (Ny+1)*m, (Ny+1)*(m-2)+1 : (Ny+1)*(m-1))    = Blk32 ;
        biHarm((Ny+1)*(m-1)+1 : (Ny+1)*m, (Ny+1)*(m-3)+1 : (Ny+1)*(m-2))    = Blk31 ;
        biHarm((Ny+1)*(m-1)+1 : (Ny+1)*m, (Ny+1)*m+1 : (Ny+1)*(m+1))        = Blk34 ;
        biHarm((Ny+1)*(m-1)+1 : (Ny+1)*m, (Ny+1)*(m+1)+1 : (Ny+1)*(m+2))    = Blk35 ;
    end

    biHarm(1:Ny+1,1:Ny+1)           = Blk11 ;
    biHarm(1:Ny+1,Ny+2:2*Ny+2)      = Blk12 ;
    biHarm(1:Ny+1,2*Ny+3:3*Ny+3)    = Blk13 ;


    biHarm(Nx*(Ny+1)+1:(Ny+1)*(Nx+1),Nx*(Ny+1)+1:(Ny+1)*(Nx+1))         = BlkMM ;
    biHarm(Nx*(Ny+1)+1:(Ny+1)*(Nx+1),(Nx-1)*(Ny+1)+1:(Ny+1)*Nx)         = BlkMMm1 ;
    biHarm(Nx*(Ny+1)+1:(Ny+1)*(Nx+1),(Nx-2)*(Ny+1)+1:(Ny+1)*(Nx-1))     = BlkMMm2 ;


    biHarm(Ny+2:2*Ny+2,1:Ny+1)          = Blk21 ;
    biHarm(Ny+2:2*Ny+2,Ny+2:2*Ny+2)     = Blk22 ;
    biHarm(Ny+2:2*Ny+2,2*Ny+3:3*Ny+3)   = Blk23 ;
    biHarm(Ny+2:2*Ny+2,3*Ny+4:4*Ny+4)   = Blk24 ;

    biHarm((Ny+1)*(Nx-1)+1:Nx*(Ny+1),Nx*(Ny+1)+1:(Ny+1)*(Nx+1))         = BlkMm1M ;
    biHarm((Ny+1)*(Nx-1)+1:Nx*(Ny+1),(Nx-1)*(Ny+1)+1:(Ny+1)*(Nx))       = BlkMm1Mm1 ;
    biHarm((Ny+1)*(Nx-1)+1:Nx*(Ny+1),(Nx-2)*(Ny+1)+1:(Ny+1)*(Nx-1))     = BlkMm1Mm2 ;
    biHarm((Ny+1)*(Nx-1)+1:Nx*(Ny+1),(Nx-3)*(Ny+1)+1:(Ny+1)*(Nx-2))     = BlkMm1Mm3 ;

    biHarm = biHarm/h^4 ;

    %% EIGENVALUES

    [Q,Dm] = eigs(biHarm,Nmodes,'smallestabs') ;
    [~,indSort] = sort(diag((Dm))) ;
    Q = Q(:,indSort) ;


    freqs = sqrt(abs(diag(Dm)))*sqrt(D/rho/Lz)/2/pi ;
    Om    = 2*pi*freqs ;



    if Nmodes > 8
        colormap('parula') ;
        for m = 1 : 9
            mdShape = reshape(Q(:,m),[(Ny+1),(Nx+1)]) ;
            subplot(3,3,m)
            mesh(real(mdShape),(abs(mdShape)),'FaceColor','texturemap') ; view(2); axis equal; axis tight;
        end
    end
end
