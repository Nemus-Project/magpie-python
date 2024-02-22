import numpy as np
from Dxx_coeffs import *
import numbers
import numpy as np
from scipy.sparse import *


def bhmat(BCs: np.ndarray, Nxy: np.ndarray, h: float, Lz: float, E: float, nu: float):
    """

    :param BCs: density [kg/m^3]
    :param Nxy: Number of grad points [Nx Ny]
    :param h:   Grid spacing
    :param Lz:  plate thickness
    :param E:   Young's mod [Pa]
    :param nu:  poisson's ratio
    :return: biharmonic matrix size (Nx Ny)
    """
    ## validate
    assert BCs.shape == (4, 2)
    assert len(Nxy) == 2 and all([type(val) is int for val in Nxy]) and all([val > 0 for val in Nxy])
    assert h != 0.0
    assert Lz != 0.0
    assert E != 0.0
    assert nu != 0.0

    ## Unpack Variables
    K0y, Kx0, KLy, KxL, R0y, Rx0, RLy, RxL = BCs.flatten()
    Nx, Ny = Nxy

    D = E * (Lz ** 3) / 12 / (1 - (nu ** 2))
    ## MATRIX BUILDER

    ##--- build matrix in blocks
    a0 = np.ones((Ny + 1))
    a1 = np.ones((Ny))
    a2 = np.ones((Ny - 1))

    D00u00, D00u10, D00u20, D00u01, D00u02, D00u11 = D00_coeffs(K0y, R0y, Kx0, Rx0, h, D, nu)
    D01u01, D01u11, D01u21, D01u00, D01u02, D01u03, D01u12, D01u10 = D01_coeffs(K0y, R0y, Rx0, h, D, nu)
    D02u02, D02u12, D02u22, D02u01, D02u03, D02u04, D02u00, D02u13, D02u11 = D02_coeffs(K0y, R0y, h, D, nu)
    D0Nu0N, D0Nu1N, D0Nu2N, D0Nu0Nm1, D0Nu0Nm2, D0Nu1Nm1 = D00_coeffs(K0y, R0y, KxL, RxL, h, D, nu)
    D0Nm1u0Nm1, D0Nm1u1Nm1, D0Nm1u2Nm1, D0Nm1u0N, D0Nm1u0Nm2, D0Nm1u0Nm3, D0Nm1u1Nm2, D0Nm1u1N = D01_coeffs(K0y, R0y,
                                                                                                            RxL, h, D,
                                                                                                            nu)
    ##-- Blk11
    D0 = D02u02 * a0
    D1 = D02u03 * a1
    D2 = D02u04 * a2
    Dm1 = D02u01 * a1
    Dm2 = D02u00 * a2

    D0[[0, 1, -2, -1]] = [D00u00, D01u01, D0Nm1u0Nm1, D0Nu0N]
    D1[[0, 1, -1, -2]] = [D00u01, D01u02, D0Nu0Nm1, D0Nm1u0Nm2]
    D2[[0, 1]] = [D00u02, D01u03]
    Dm1[[0, -1]] = [D01u00, D0Nm1u0N]
    Dm1[[-1, -2]] = [D0Nu0Nm2, D0Nm1u0Nm3]
    Blk11 = diags(D0, 0) + diags(D1, 1) + diags(Dm1, -1) + diags(D2, 2) + diags(Dm2, -2)

    ###-- Blk12
    D0 = D02u12 * a0
    D1 = D02u13 * a1
    Dm1 = D02u11 * a1
    D0[[-1, -2, 0, 1]] = [D0Nu1N, D0Nm1u1Nm1, D00u10, D01u11]
    D1[[-1, -2, 0, 1]] = [D0Nu1Nm1, D0Nm1u1Nm2, D00u11, D01u12]
    Dm1[[-1, 0]] = [D0Nm1u1N, D01u10]

    Blk12 = diags(D0, 0) + diags(D1, 1) + diags(Dm1, -1)
    #
    #
    # -- Blk13
    D0 = D02u22 * a0
    D0[[0, 1, -1, -2]] = [D00u20, D01u21, D0Nu2N, D0Nm1u2Nm1]
    Blk13 = diags(D0, 0)

    [D10u10, D10u20, D10u30, D10u00, D10u11, D10u12, D10u21, D10u01] = D10_coeffs(R0y, Kx0, Rx0, h, D, nu)
    [D11u11, D11u12, D11u13, D11u10, D11u01, D11u21, D11u31, D11u22, D11u20, D11u00, D11u02] = D11_coeffs(R0y, Rx0, h,
                                                                                                          D, nu)
    [D12u12, D12u13, D12u14, D12u11, D12u10, D12u02, D12u22, D12u32, D12u23, D12u21, D12u01, D12u03] = D12_coeffs(R0y,
                                                                                                                  h, D,
                                                                                                                  nu)
    [D1Nu1N, D1Nu2N, D1Nu3N, D1Nu0N, D1Nu1Nm1, D1Nu1Nm2, D1Nu2Nm1, D1Nu0Nm1] = D10_coeffs(R0y, KxL, RxL, h, D, nu)
    [D1Nm1u1Nm1, D1Nm1u1Nm2, D1Nm1u1Nm3, D1Nm1u1N, D1Nm1u0Nm1, D1Nm1u2Nm1, D1Nm1u3Nm1, D1Nm1u2Nm2, D1Nm1u2N, D1Nm1u0N,
     D1Nm1u0Nm2] = D11_coeffs(R0y, RxL, h, D, nu)
    ##
    ##-- Blk21
    D0 = D12u02 * a0
    D1 = D12u03 * a1
    Dm1 = D12u01 * a1
    D0[[-1, -2, 0, 1]] = [D1Nu0N, D1Nm1u0Nm1, D10u00, D11u01]
    D1[[-1, -2, 0, 1]] = [D1Nu0Nm1, D1Nm1u0Nm2, D10u01, D11u02]
    Dm1[[-1, 0]] = [D1Nm1u0N, D11u00]

    Blk21 = diags(D0, 0) + diags(D1, 1) + diags(Dm1, -1)

    ##
    ##-- Blk22
    D0 = D12u12 * a0
    D1 = D12u13 * a1
    D2 = D12u14 * a2
    Dm1 = D12u11 * a1
    Dm2 = D12u10 * a2
    D0[[0, 1, -1, -2]] = [D10u10, D11u11, D1Nu1N, D1Nm1u1Nm1]
    D1[[0, 1, -1, -2]] = [D10u11, D11u12, D1Nu1Nm1, D1Nm1u1Nm2]
    D2[[0, 1]] = [D10u12, D11u13]
    Dm1[[0, -1]] = [D11u10, D1Nm1u1N]
    Dm1[[-1, -2]] = [D1Nu1Nm2, D1Nm1u1Nm3]
    Blk22 = diags(D0, 0) + diags(D1, 1) + diags(Dm1, -1) + diags(D2, 2) + diags(Dm2, -2)

    ##-- Blk23
    D0 = D12u22 * a0
    D1 = D12u23 * a1
    Dm1 = D12u21 * a1
    D0[[-1, -2, 0, 1]] = [D1Nu2N, D1Nm1u2Nm1, D10u20, D11u21]
    D1[[-1, -2, 0, 1]] = [D1Nu2Nm1, D1Nm1u2Nm2, D10u21, D11u22]
    Dm1[[-1, 0]] = [D1Nm1u2N, D11u20]
    Blk23 = diags(D0, 0) + diags(D1, 1) + diags(Dm1, -1)
    #
    #
    # -- Blk24
    D0 = D12u32 * a0
    #
    D0[[0, 1, -1, -2]] = [D10u30, D11u31, D1Nu3N, D1Nm1u3Nm1]
    Blk24 = diags(D0)
    #
    [D20u20, D20u21, D20u22, D20u10, D20u30, D20u40, D20u00, D20u31, D20u11] = D20_coeffs(Kx0, Rx0, h, D, nu)
    [D21u21, D21u22, D21u23, D21u20, D21u11, D21u31, D21u41, D21u01, D21u32, D21u30, D21u10, D21u12] = D21_coeffs(Rx0,
                                                                                                                  h, D,
                                                                                                                  nu)
    [D22u20, D22u11, D22u21, D22u31, D22u02, D22u12, D22u22, D22u32, D22u42, D22u13, D22u23, D22u33,
     D22u24] = D22_coeffs()
    [D2Nu2N, D2Nu2Nm1, D2Nu2Nm2, D2Nu1N, D2Nu3N, D2Nu4N, D2Nu0N, D2Nu3Nm1, D2Nu1Nm1] = D20_coeffs(KxL, RxL, h, D, nu)
    [D2Nm1u2Nm1, D2Nm1u2Nm2, D2Nm1u2Nm3, D2Nm1u2N, D2Nm1u1Nm1, D2Nm1u3Nm1, D2Nm1u4Nm1, D2Nm1u0Nm1, D2Nm1u3Nm2, D2Nm1u3N,
     D2Nm1u1N, D2Nm1u1Nm2] = D21_coeffs(RxL, h, D, nu)

    ##-- Blk31
    D0 = D22u02 * a0
    D0[[0, 1, -1, -2]] = [D20u00, D21u01, D2Nu0N, D2Nm1u0Nm1]
    Blk31 = diags(D0, 0)

    # -- Blk32
    D0 = D22u12 * a0
    D1 = D22u13 * a1
    Dm1 = D22u11 * a1
    D0[[-1, -2, 0, 1]] = [D2Nu1N, D2Nm1u1Nm1, D20u10, D21u11]
    D1[[-1, -2, 0, 1]] = [D2Nu1Nm1, D2Nm1u1Nm2, D20u11, D21u12]
    Dm1[[-1, 0]] = [D2Nm1u1N, D21u10]

    Blk32 = diags(D0, 0) + diags(D1, 1) + diags(Dm1, -1)

    ##-- Blk33
    D0 = D22u22 * a0
    D1 = D22u23 * a1
    D2 = D22u24 * a2
    Dm1 = D22u21 * a1
    Dm2 = D22u20 * a2
    D0[[-1, -2, 0, 1]] = [D2Nu2N, D2Nm1u2Nm1, D20u20, D21u21]
    D1[[-1, -2, 0, 1]] = [D2Nu2Nm1, D2Nm1u2Nm2, D20u21, D21u22]
    D2[[0, 1]] = [D20u22, D21u23]
    Dm1[[-1, 0]] = [D2Nm1u2N, D21u20]
    Dm1[[-1, -2]] = [D2Nu2Nm2, D2Nm1u2Nm3]
    Blk33 = diags(D0, 0) + diags(D1, 1) + diags(Dm1, -1) + diags(D2, 2) + diags(Dm2, -2)

    ##-- Blk34
    D0 = D22u32 * a0
    D1 = D22u33 * a1
    Dm1 = D22u31 * a1
    D0[[-1, -2, 0, 1]] = [D2Nu3N, D2Nm1u3Nm1, D20u30, D21u31]
    D1[[-1, -2, 0, 1]] = [D2Nu3Nm1, D2Nm1u3Nm2, D20u31, D21u32]
    Dm1[[-1, 0]] = [D2Nm1u3N, D21u30]

    Blk34 = diags(D0, 0) + diags(D1, 1) + diags(Dm1, -1)

    # -- Blk35
    D0 = D22u42 * a0
    D0[[-1, -2, 0, 1]] = [D2Nu4N, D2Nm1u4Nm1, D20u40, D21u41]
    Blk35 = diags(D0)
    #
    [D00u00, D00u10, D00u20, D00u01, D00u02, D00u11] = D00_coeffs(KLy, RLy, Kx0, Rx0, h, D, nu)
    [D01u01, D01u11, D01u21, D01u00, D01u02, D01u03, D01u12, D01u10] = D01_coeffs(KLy, RLy, Rx0, h, D, nu)
    [D02u02, D02u12, D02u22, D02u01, D02u03, D02u04, D02u00, D02u13, D02u11] = D02_coeffs(KLy, RLy, h, D, nu)
    [D0Nu0N, D0Nu1N, D0Nu2N, D0Nu0Nm1, D0Nu0Nm2, D0Nu1Nm1] = D00_coeffs(KLy, RLy, KxL, RxL, h, D, nu)
    [D0Nm1u0Nm1, D0Nm1u1Nm1, D0Nm1u2Nm1, D0Nm1u0N, D0Nm1u0Nm2, D0Nm1u0Nm3, D0Nm1u1Nm2, D0Nm1u1N] = D01_coeffs(KLy, RLy,
                                                                                                              RxL, h, D,
                                                                                                              nu)

    # -- BlkMM
    D0 = D02u02 * a0
    D1 = D02u03 * a1
    D2 = D02u04 * a2
    Dm1 = D02u01 * a1
    Dm2 = D02u00 * a2
    #
    D0[[-1, -2, 0, 1]] = [D0Nu0N, D0Nm1u0Nm1, D00u00, D01u01]
    D1[[-1, -2, 0, 1]] = [D0Nu0Nm1, D0Nm1u0Nm2, D00u01, D01u02]
    D2[[0, 1]] = [D00u02, D01u03]
    Dm1[[-1, 0]] = [D0Nm1u0N, D01u00]
    Dm1[[-1, -2]] = [D0Nu0Nm2, D0Nm1u0Nm3]

    BlkMM = diags(D0, 0) + diags(D1, 1) + diags(Dm1, -1) + diags(D2, 2) + diags(Dm2, -2)

    # -- BlkMMm1
    D0 = D02u12 * a0
    D1 = D02u13 * a1
    Dm1 = D02u11 * a1
    D0[[-1, -2, 0, 1]] = [D0Nu1N, D0Nm1u1Nm1, D00u10, D01u11]
    D1[[-1, -2, 0, 1]] = [D0Nu1Nm1, D0Nm1u1Nm2, D00u11, D01u12]
    Dm1[[-1, 0]] = [D0Nm1u1N, D01u10]

    BlkMMm1 = diags(D0, 0) + diags(D1, 1) + diags(Dm1, -1)

    ##-- BlkMMm2
    D0 = D02u22 * a0
    #
    D0[[0, 1, -1, -2]] = [D00u20, D01u21, D0Nu2N, D0Nm1u2Nm1]
    BlkMMm2 = diags(D0)
    #
    [D10u10, D10u20, D10u30, D10u00, D10u11, D10u12, D10u21, D10u01] = D10_coeffs(RLy, Kx0, Rx0, h, D, nu)
    [D11u11, D11u12, D11u13, D11u10, D11u01, D11u21, D11u31, D11u22, D11u20, D11u00, D11u02] = D11_coeffs(RLy, Rx0, h,
                                                                                                          D, nu)
    [D12u12, D12u13, D12u14, D12u11, D12u10, D12u02, D12u22, D12u32, D12u23, D12u21, D12u01, D12u03] = D12_coeffs(RLy,
                                                                                                                  h, D,
                                                                                                                  nu)
    [D1Nu1N, D1Nu2N, D1Nu3N, D1Nu0N, D1Nu1Nm1, D1Nu1Nm2, D1Nu2Nm1, D1Nu0Nm1] = D10_coeffs(RLy, KxL, RxL, h, D, nu)
    [D1Nm1u1Nm1, D1Nm1u1Nm2, D1Nm1u1Nm3, D1Nm1u1N, D1Nm1u0Nm1, D1Nm1u2Nm1, D1Nm1u3Nm1, D1Nm1u2Nm2, D1Nm1u2N, D1Nm1u0N,
     D1Nm1u0Nm2] = D11_coeffs(RLy, RxL, h, D, nu)

    ##-- BlkMm1M
    D0 = D12u02 * a0
    D1 = D12u03 * a1
    Dm1 = D12u01 * a1
    #
    D0[[-1, -2, 0, 1]] = [D1Nu0N, D1Nm1u0Nm1, D10u00, D11u01]
    D1[[-1, -2, 0, 1]] = [D1Nu0Nm1, D1Nm1u0Nm2, D10u01, D11u02]
    Dm1[[-1, 0]] = [D1Nm1u0N, D11u00]
    BlkMm1M = diags(D0, 0) + diags(D1, 1) + diags(Dm1, -1)

    ##-- BlkMm1Mm1
    D0 = D12u12 * a0
    D1 = D12u13 * a1
    D2 = D12u14 * a2
    Dm1 = D12u11 * a1
    Dm2 = D12u10 * a2

    D0[[-1, -2, 0, 1]] = [D1Nu1N, D1Nm1u1Nm1, D10u10, D11u11]
    D1[[-1, -2, 0, 1]] = [D1Nu1Nm1, D1Nm1u1Nm2, D10u11, D11u12]
    D2[[0, 1]] = [D10u12, D11u13]
    Dm1[[-1, 0]] = [D1Nm1u1N, D11u10]
    Dm1[[-1, -2]] = [D1Nu1Nm2, D1Nm1u1Nm3]

    #
    BlkMm1Mm1 = diags(D0) + diags(D1, 1) + diags(Dm1, -1) + diags(D2, 2) + diags(Dm2, -2)
    #
    # -- BlkMm1Mm2
    D0 = D12u22 * a0
    D1 = D12u23 * a1
    Dm1 = D12u21 * a1
    #
    D0[[-1, -2, 0, 1]] = [D1Nu2N, D1Nm1u2Nm1, D10u20, D11u21]
    D1[[-1, -2, 0, 1]] = [D1Nu2Nm1, D1Nm1u2Nm2, D10u21, D11u22]
    Dm1[[-1, 0]] = [D1Nm1u2N, D11u20]

    BlkMm1Mm2 = diags(D0) + diags(D1, 1) + diags(Dm1, -1)
    
    # -- BlkMm1Mm3
    D0 = D12u32 * a0

    D0[[0, 1, -1, -2]] = [D10u30, D11u31, D1Nu3N, D1Nm1u3Nm1]
    BlkMm1Mm3 = diags(D0)

    # biHarm = sparse((Nx+1)*(Ny+1),(Nx+1)*(Ny+1))
    #
    #
    # for m in np.r_[1 : Nx - 1]:
    #     biHarm[(Ny+1)*(m-1)+1 : (Ny+1)*m, (Ny+1)*(m-1)+1 : (Ny+1)*])        = Blk33
    #     biHarm[(Ny+1)*(m-1)+1 : (Ny+1)*m, (Ny+1)*(m-2)+1 : (Ny+1)*(m-1)]    = Blk32
    #     biHarm[(Ny+1)*(m-1)+1 : (Ny+1)*m, (Ny+1)*(m-3)+1 : (Ny+1)*(m-2)]    = Blk31
    #     biHarm[(Ny+1)*(m-1)+1 : (Ny+1)*m, (Ny+1)*m+1 : (Ny+1)*(m+1)]        = Blk34
    #     biHarm[(Ny+1)*(m-1)+1 : (Ny+1)*m, (Ny+1)*(m+1)+1 : (Ny+1)*(m+2)]    = Blk35
    #
    #
    # biHarm[0:Ny+1,0:Ny+1]           = Blk11
    # biHarm[0:Ny+1,Ny+2:2*Ny+2]      = Blk12
    # biHarm[0:Ny+1,2*Ny+3:3*Ny+3]    = Blk13
    #
    #
    # biHarm[Nx*(Ny+1)+1:(Ny+1)*(Nx+1),Nx*(Ny+1)+1:(Ny+1)*[Nx+1)] ]        = BlkMM
    # biHarm[Nx*(Ny+1)+1:(Ny+1)*(Nx+1),(Nx-1)*(Ny+1)+1:[Ny+1)*Nx]]         = BlkMMm1
    # biHarm[Nx*(Ny+1)+1:(Ny+1)*(Nx+1),(Nx-2)*(Ny+1)+1:(Ny+1)*[Nx-1)] ]    = BlkMMm2
    #
    #
    # biHarm[Ny+2:2*Ny+2,1:Ny+1]          = Blk21
    # biHarm[Ny+2:2*Ny+2,Ny+2:2*Ny+2]     = Blk22
    # biHarm[Ny+2:2*Ny+2,2*Ny+3:3*Ny+3]   = Blk23
    # biHarm[Ny+2:2*Ny+2,3*Ny+4:4*Ny+4]   = Blk24
    #
    # biHarm[(Ny+1)*(Nx-1)+1:Nx*(Ny+1),Nx*(Ny+1)+1:(Ny+1)*[Nx+1)]]         = BlkMm1M
    # biHarm[(Ny+1)*(Nx-1)+1:Nx*(Ny+1),(Nx-1)*(Ny+1)+1:(Ny+1)*[Nx)]]       = BlkMm1Mm1
    # biHarm[(Ny+1)*(Nx-1)+1:Nx*(Ny+1),(Nx-2)*(Ny+1)+1:(Ny+1)*[Nx-1)]]     = BlkMm1Mm2
    # biHarm[(Ny+1)*(Nx-1)+1:Nx*(Ny+1),(Nx-3)*(Ny+1)+1:(Ny+1)*[Nx-2)]]     = BlkMm1Mm3
    #
    # biHarm = biHarm / (h ** 4)
    # return biHarm


def main():
    Lx = 1.10
    Ly = 0.8
    Lz = 5e-3
    ldim = [Lx, Ly, Lz]  # -- plate dimensions [x, y, z] in metres
    E = 9.0e+9  # -- Young's mod [Pa]
    rho = 8765  # -- density [kg/m^3]
    nu = 0.3  # -- poisson's ratio
    Nmodes = 25  # -- number of modes to compute
    h = np.sqrt(Lx * Ly) * 0.01  # -- grid spacing
    BCs = np.ones((4, 2)) * 1e15  # -- elastic constants around the edges

    D = E * (Lz ** 3) / 12 / (1 - (nu ** 2))
    Nx = int(np.floor(Lx / h))
    Ny = int(np.floor(Ly / h))

    biHarm = bhmat(BCs, [Nx, Ny], h, Lz, E, nu)


if __name__ == '__main__':
    main()
