import numpy as np
from numpy.linalg import pinv
from magpie import magpie


def modal_time_integration(rho: float, E: float, nu: float, ldim: list, BCs: np.ndarray, sig: list, maxFreq: float,
                           pos: dict, T: float = 6, fs: float = 44100, AmpF: float = 30, twid: float = 0.0006):
    """
    
    :param rho:
    :param E:
    :param nu:
    :param ldim:
    :param maxFreq:
    :param BCs:
    :param sig:
    :param pos:
    :param T:
    :param fs:
    :param AmpF:
    :param twid:
    :return:
    """

    Lx, Ly, Lz = ldim
    D = E * (Lz ** 3) / 12 / (1 - (nu ** 2))
    in_coord = pos['in']
    out_coord = {k: pos[k] for k in ('l', 'r')}

    # --- Simulation Parameters

    Ts = np.round(T * fs)
    k = 1 / fs
    tv = np.r_[:Ts] * k  # -- time axis array
    fv = np.r_[:Ts] * fs / Ts  # -- freq axis array

    h = np.sqrt(np.sqrt(D / rho / Lz * 16 / ((maxFreq * 2 * np.pi) ** 2)))  # -- set according to largest freq

    Om, Q, N, biHarm = magpie(rho, E, nu, ldim, h, BCs)

    fOm = 0
    Nmodes = 0
    OmDsq = 0

    Nx = N['x']
    Ny = N['y']

    while fOm < maxFreq and Nmodes < (Nx + 1) * (Ny + 1) and OmDsq >= 0:
        Nmodes += 1
        fOm = Om[Nmodes] / 2 / np.pi
        C = sig[0] + sig[1] * (Om[Nmodes] ** 2)
        OmDsq = (Om[Nmodes] ** 2) - (C ** 2)

    # Nmodes = Nmodes - 1
    fMax = Om[Nmodes] / 2 / np.pi  # -- check if this is in the range of maxFreq

    Om = Om[:Nmodes]
    Q = Q[:, :Nmodes].real
    C = (sig[0] + sig[1] * (Om ** 2))
    OmD = np.sqrt((Om ** 2) - (C ** 2))

    # -- build input vector (spreading, lin interp)

    Jin = np.zeros(((Nx) * (Ny), 1))

    nx = in_coord[0] * Nx
    Min = np.floor(nx)
    alx = nx - Min

    ny = in_coord[1] * Ny
    Nin = np.floor(ny)
    aly = ny - Nin

    Jin[int((Ny) * Min + Nin + 1)] = alx * aly
    Jin[int((Ny) * (Min + 1) + Nin + 1)] = (1 - alx) * aly
    Jin[int((Ny) * Min + Nin + 2)] = alx * (1 - aly)
    Jin[int((Ny) * (Min + 1) + Nin + 2)] = (1 - alx) * (1 - aly)

    Jin = Jin / (h ** 2) / rho / Lz
    Jin = (pinv(Q).real @ Jin).flatten()

    # -- left output weights (in interp)
    JoutL = np.zeros(((Nx) * (Ny)))

    outx = out_coord['l'][0] * Nx
    Mout = np.floor(outx)
    alx = outx - Mout

    outy = out_coord['l'][1] * Ny
    Nout = np.floor(outy)
    aly = outy - Nout

    JoutL[int((Ny) * Mout + Nout + 1)] = alx * aly
    JoutL[int((Ny) * (Mout + 1) + Nout + 1)] = (1 - alx) * aly
    JoutL[int((Ny) * Mout + Nout + 2)] = alx * (1 - aly)
    JoutL[int((Ny) * (Mout + 1) + Nout + 2)] = (1 - alx) * (1 - aly)

    # -- right output weights (in interp)
    JoutR = np.zeros(((Nx ) * (Ny )))

    outx = out_coord['r'][0] * Nx
    Mout = np.floor(outx)
    alx = outx - Mout

    outy = out_coord['r'][1] * Ny
    Nout = np.floor(outy)
    aly = outy - Nout

    JoutR[int((Ny) * Mout + Nout + 1)] = alx * aly
    JoutR[int((Ny) * (Mout + 1) + Nout + 1)] = (1 - alx) * aly
    JoutR[int((Ny) * Mout + Nout + 2)] = alx * (1 - aly)
    JoutR[int((Ny) * (Mout + 1) + Nout + 2)] = (1 - alx) * (1 - aly)

    # -- input forcing
    Nfin = int(np.floor(twid * fs))
    fin = np.zeros((Ts))
    fin[:Nfin] = 0.5 * AmpF * (1 - np.cos(2 * np.pi * np.r_[:Nfin] / Nfin))

    # --- init
    vm = np.zeros((Nmodes))
    v0 = np.zeros((Nmodes))
    outL = np.zeros((Ts))
    outR = np.zeros((Ts))
    velL = np.zeros((Ts))
    velR = np.zeros((Ts))
    outLprev = 0
    outRprev = 0

    # ----------------------------------
    # -- main loop
    for n in range(Ts):
        vp = 2 * np.exp(-C * k) * np.cos(OmD * k) * v0 - np.exp(-2 * C * k) * vm + (k ** 2) * Jin * fin[n]
        outLcur = JoutL @ (Q @ v0)
        outRcur = JoutR @ (Q @ v0)
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
    rho = 8765  # -- density [kg/m^3]
    E = 101e9  # -- Young's mod [Pa]
    nu = 0.3  # -- poisson's ratio

    ldim = [0.151, 0.08, 0.81e-3]

    # elastic constants around the edges (this allows to set the various bcs)
    BCs = np.zeros((4, 2)) * 1e15  # -- elastic constants around the edges
    BCs[1, :] = 1e15

    sig = [5e-3, 3e-9]  # -- damping parameters: T60 = 3*log(10)./(sig0+Om.^2*sig1)

    maxFreq = 15000.0  # max frequency to consider in hz

    # -- input / output locations, FRACTIONS of [Lx Ly] (values must be >0 and <1)
    pos = {
        'in': [0.54, 0.78],
        'l': [0.57, 0.75],
        'r': [0.56, 0.65]
    }

    modal_time_integration(rho, E, nu, ldim, BCs, sig, maxFreq, pos)
