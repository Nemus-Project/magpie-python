import numpy as np
from numpy.linalg import pinv
from magpie import magpie


def modal_time_integration(rho: float, E: float, nu: float, ldim: list, h: float, BCs: np.ndarray, sig: list,
                           pos: dict):
    """"""
    Lx, Ly, Lz = ldim
    D = E * (Lz ** 3) / 12 / (1 - (nu ** 2))
    Nx = int(np.ceil(Lx / h))
    Ny = int(np.ceil(Ly / h))

    in_coord = pos['in']  # [0.54, 0.78]
    out_coord = pos['out']  # [0.57, 0.75]

    # --- Simulation Parameters
    T = 6  # sim length [s]
    fs = 44100  # sample rate [Hz]
    AmpF = 30  # force amplitude [N]
    twid = 0.0006  # forcing temporal width [s]

    Ts = np.round(T * fs)
    k = 1 / fs
    tv = np.r_[:Ts] * k  # -- time axis array
    fv = np.r_[:Ts] * fs / Ts  # -- freq axis array

    maxFreq = 15000

    h = np.np.sqrt(np.np.sqrt(D / rho / Lz * 16 / ((maxFreq * 2 * np.pi) ** 2)))  # -- set according to largest freq
    Om, Q, Nx, Ny, biHarm, Dm = magpie(rho, E, nu, ldim, h, BCs)

    fOm = 0
    Nmodes = 0
    OmDsq = 0

    while fOm < maxFreq and Nmodes < (Nx + 1) * (Ny + 1) and OmDsq >= 0:
        Nmodes += 1
        fOm = Om[Nmodes] / 2 / np.pi
        C = sig[0] + sig[1] * (Om[Nmodes] ** 2)
        OmDsq = (Om[Nmodes] ** 2) - (C ** 2)

    Nmodes = Nmodes - 1
    fMax = Om[Nmodes] / 2 / np.pi  # -- check if this is in the range of maxFreq

    Om = Om[1:Nmodes]
    Q = Q[:, 1:Nmodes]
    C = sig[0] + sig[1] * (Om ** 2)
    OmD = np.sqrt((Om ** 2) - (C ** 2))

    # -- build input vector (spreading, lin interp)

    Jin = np.zeros((Nx + 1) * (Ny + 1), 1)

    nx = in_coord[0] * Nx
    Min = np.floor(nx)
    alx = nx - Min

    ny = in_coord[1] * Ny
    Nin = np.floor(ny)
    aly = ny - Nin

    Jin[(Ny + 1) * Min + Nin + 1] = alx * aly
    Jin[(Ny + 1) * (Min + 1) + Nin + 1] = (1 - alx) * aly
    Jin[(Ny + 1) * Min + Nin + 2] = alx * (1 - aly)
    Jin[(Ny + 1) * (Min + 1) + Nin + 2] = (1 - alx) * (1 - aly)

    Jin = Jin / (h ** 2) / rho / Lz
    Jin = pinv(Q) * Jin

    # -- left output weights (in interp)
    JoutL = np.zeros(1, (Nx + 1) * (Ny + 1))

    outx = out_coord['l'][0] * Nx
    Mout = np.floor(outx)
    alx = outx - Mout

    outy = out_coord['l'][1] * Ny
    Nout = np.floor(outy)
    aly = outy - Nout

    JoutL[(Ny + 1) * Mout + Nout + 1] = alx * aly
    JoutL[(Ny + 1) * (Mout + 1) + Nout + 1] = (1 - alx) * aly
    JoutL[(Ny + 1) * Mout + Nout + 2] = alx * (1 - aly)
    JoutL[(Ny + 1) * (Mout + 1) + Nout + 2] = (1 - alx) * (1 - aly)

    # -- right output weights (in interp)
    JoutR = np.zeros((1, (Nx + 1) * (Ny + 1)))

    outx = out_coord['r'][0] * Nx
    Mout = np.floor(outx)
    alx = outx - Mout

    outy = out_coord['r'][1] * Ny
    Nout = np.floor(outy)
    aly = outy - Nout

    JoutR[(Ny + 1) * Mout + Nout + 1] = alx * aly
    JoutR[(Ny + 1) * (Mout + 1) + Nout + 1] = (1 - alx) * aly
    JoutR[(Ny + 1) * Mout + Nout + 2] = alx * (1 - aly)
    JoutR[(Ny + 1) * (Mout + 1) + Nout + 2] = (1 - alx) * (1 - aly)

    # -- input forcing
    Nfin = np.floor(twid * fs)
    fin = np.zeros(Ts, 1)
    fin[:Nfin] = 0.5 * AmpF * (1 - np.np.cos(2 * np.pi * np.r_[:Nfin] / Nfin))

    # --- init
    vm = np.zeros(Nmodes, 1)
    v0 = np.zeros(Nmodes, 1)
    outL = np.zeros(Ts, 1)
    outR = np.zeros(Ts, 1)
    velL = np.zeros(Ts, 1)
    velR = np.zeros(Ts, 1)
    outLprev = 0
    outRprev = 0

    # ----------------------------------
    # -- main loop
    for n in range(Ts):
        vp = 2 * np.exp(-C * k) * np.cos(OmD * k) * v0 - np.exp(-2 * C * k) * vm + (k ** 2) * Jin * fin(n)
        outLcur = JoutL * (Q * v0)
        outRcur = JoutR * (Q * v0)
        outL[n] = outLcur
        outR[n] = outRcur
        velL[n] = (outLcur - outLprev) / k
        velR[n] = (outRcur - outRprev) / k
        vm, v0 = v0, vp

        outLprev, outRprev = outLcur, outRcur

    out = [outL, outR]
    vel = [velL, velR]

    return out, vel


if __name__ == '__main__':
    pass
