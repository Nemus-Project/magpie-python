import numpy as np

try:
    from magpie import magpie
except ImportError:
    from .magpie import magpie

def youngcalc(rho: float, ldim: list, h: float, BCs: np.ndarray, ExpFreq: list, Ntrain: int, should_plot: bool = False):
    """
    Estimate Young's modulus (E) of an experimental plate starting from a batch
    of experimentally measured frequencies, leveraging MAGPIE
    
    :param rho: the experimental plate density
    :type rho: float
    :param ldim:  a 3 x 1 array containing the Lx Ly Lz dimensions
    :type ldim: list
    :param h: the grid spacing of the FD scheme
    :type h: float
    :param BCs: a 4 x 2 array containing the rigidities of the boundary supports of the experimental plate
    :type BCs: np.array
    :param ExpFreq: an array contaning the first Nmodes measured modal frequencies
    :type ExpFreq: np.array
    :param Ntrain: an integer. Must be <= Nmodes. It is the number of training modes out of the available batch
    :type Ntrain: int
    :return: The estimated young's modulus for the given plate
    :rtype: float

    :Example:

    .. code:: python
        :linenos:

        ExpFreq = [73.2 148 376 431 559 910]   #-- these are measured from a plate
        rho     = 8765             #-- density [kg/(m ** 3)]
        Lx      = 0.1
        Ly      = 0.08
        Lz      = 0.00081
        BCs = np.array([[0,0],
                        [1e15,1e15],
                        [0,0],
                        [0,0],])

        ldim    = [Lx Ly Lz]
        h       = np.sqrt(Lx*Ly)*0.01   #-- grid spacing [m]

        E  = youngcalc(rho, ldim, h, BCs, ExpFreq, 3)


    """

    Nmodes = len(ExpFreq)

    if (Ntrain > Nmodes):
        raise Exception('Ntrain should be less than total experimental Freqs (i.e. Ntrain < len(ExpFreq))')

    TrainFreq = np.array(ExpFreq[0:Ntrain])
    TestFreq = [] if Ntrain == Nmodes else ExpFreq[Ntrain + 1:]

    # -- zero parameters
    E0 = 2e11  # -- Young's modulus [Pa] (just a number here, results should not change if this changes)
    nu = 0.3  # -- poisson's ratio (average value for metals)
    Lx, Ly, Lz = ldim

    A = Lx * Ly  # -- area
    D0 = E0 * (Lz ** 3) / 12 / (1 - (nu ** 2))  # -- zero-rigidity

    # -- run magpie and get non-dimensional freqs
    Om, _, _, _ = magpie(rho, E0, nu, ldim, h, BCs, Ntrain, "none")
    OmNDim = (Om / np.sqrt(D0)) * np.sqrt(rho * (A ** 2) * Lz)
    OmNDimsq = np.reshape(OmNDim ** 2,(-1,1)).T

    # -- least-square (LS) optimisation
    psi = np.reshape(((TrainFreq * 2 * np.pi) ** 2) * rho * Lz * (A ** 2),(-1,1))
    DLS = ((OmNDimsq @ psi) / (OmNDimsq @ OmNDimsq.T))[0,0]
    ELS = DLS / ((Lz ** 3) / 12 / (1 - (nu ** 2)))

    # -- launch a numerical simulation to get the frequencies of the numerical model
    # -- using the estimated Youngs Mod
    # -- and compare against the experimental freqs
    NumOm, _, _, _ = magpie(rho, ELS, nu, ldim, h, BCs, Nmodes, "none")

    return ELS

if __name__ == '__main__':
    ExpFreq = [73.2, 148, 376, 431, 559, 910]  # -- these are measured from a plate
    rho = 8765  # -- density [kg/(m ** 3)]
    Lx = 0.1
    Ly = 0.08
    Lz = 0.00081
    ldim = [Lx, Ly, Lz]

    BCs = np.zeros((4, 2)) * 1e15  # -- elastic constants around the edges
    BCs[1, :] = 1e15

    h = np.sqrt(Lx * Ly) * 0.01  # -- grid spacing [m]

    E = youngcalc(rho, ldim, h, BCs, ExpFreq, 3)
