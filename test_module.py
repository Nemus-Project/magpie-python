import magpie as mg
import numpy as np
import sounddevice as sd

if __name__ == '__main__':
    BCs = np.ones((4, 2)) * 1e15  # -- elastic constants around the edges

    rho = 7820
    E = 200e9
    nu = 0.3
    h = 0.01

    [Om, Q, N, biharm] = mg.magpie(rho, E, nu, [1, 0.8, 5e-3], 0.01, BCs, 5)

    Lx = 1.51
    Ly = 0.8
    Lz = 0.81e-3
    ldim = [Lx, Ly, Ly]

    sig = [5e-3, 3e-9]  # -- damping parameters: T60 = 3*log(10)./(sig0+Om.^2*sig1)

    maxFreq = 15000.0  # max frequency to consider in hz

    # -- input / output locations, FRACTIONS of [Lx Ly] (values must be >0 and <1)
    pos = {
        'in': [0.54, 0.78],
        'l': [0.57, 0.75],
        'r': [0.56, 0.65]
    }

    simulation_time = 1.0

    audio, _ = mg.modal_time_integration(rho, E, nu, ldim, BCs, sig, maxFreq, pos, T=simulation_time)
    norm_gain = np.abs(audio).max()
    audio /= norm_gain
    sd.play(audio, 44100)


    ExpFreq = [73.2, 148, 376, 431, 559, 910]  # -- these are measured from a plate
    h = np.sqrt(Lx * Ly) * 0.01  # -- grid spacing [m]
    E = mg.youngcalc(rho, ldim, h, BCs, ExpFreq, 3)

    print(E)
    Nx = 30
    Ny = 30
    biHarm = mg.bhmat(BCs, [Nx, Ny], h, Lz, E, nu)
    print(type(biHarm))