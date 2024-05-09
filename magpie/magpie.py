import numpy as np
from scipy.sparse import *
from scipy.sparse.linalg import eigs
from scipy.sparse import eye
import matplotlib
from matplotlib import pyplot as plt

try:
    from bhmat import bhmat
except ImportError:
    from .bhmat import bhmat
    


def magpie(rho: float, E: float, nu: float, ldim: list, h: float, BCs: np.ndarray, Nm: int = 0, plot_type: str = None,
           base_mode: float = 0.0):
    """

    :param rho: density [kg/m^3]
    :param E: Young's mod [Pa]
    :param nu: poisson's ratio
    :param ldim: plate dimensions in meters [Lx, Ly, Lz], where Lz is thickness
    :param h: grid spacing
    :param BCs:
    :param Nm: Number of modes, if 0 maximum number of modes are calculated
    :param plot_type: style to plot mode shapes 'chladni' or '3D'
    :return: [Om, Q, {'x': Nx, 'y': Ny}, biharm]
    """
    ## Validate
    assert rho is not None
    assert E is not None
    assert nu is not None
    assert len(ldim) == 3
    assert h is not None
    assert BCs.shape == (4, 2)
    assert Nm is not None and Nm >= 0
    assert plot_type in ["chladni", "3D", "none"] or plot_type is None
    ##----------------------------
    Lx, Ly, Lz = ldim
    D = E * (Lz ** 3) / 12 / (1 - (nu ** 2))
    Nx = int(np.ceil(Lx / h))
    Ny = int(np.ceil(Ly / h))
    Nmodes = ((Nx * Ny) - 2) if Nm == 0 else Nm
    ##----------------------------
    ## Build BiHarmonic
    biharm = bhmat(BCs, [Nx, Ny], h, Lz, E, nu)

    ## EIGENVALUES
    [Dm, Q] = eigs(biharm, k=Nmodes, sigma=base_mode, which='LR')

    Om = np.sqrt(abs(Dm)) * np.sqrt(D / rho / Lz)
    hz = Om / (2 * np.pi)
    indSort = np.argsort(Dm)
    
    Om = Om[indSort]
    Q = Q[:, indSort]

    X = np.arange(0, Ny)
    Y = np.arange(0, Nx)
    X, Y = np.meshgrid(X, Y)

    sq = int(np.ceil(np.sqrt(Nmodes)))

    if plot_type == 'chladni':
        for m in range(Nmodes):
            ax = plt.subplot(sq, sq, m + 1)
            Z = abs(np.reshape(Q[:, m], [Nx, Ny]))
            chladni = plt.pcolormesh(Z.T, cmap='copper_r', shading='gouraud')
            ax.set_axis_off()
            chladni.set_clim(0.000, 0.002)

        plt.show()

    elif plot_type == '3D':

        fig, axes = plt.subplots(sq, sq, subplot_kw={"projection": "3d"})

        for m in range(Nmodes):
            Z = np.reshape(Q[:, m], [Nx, Ny])
            axes[m // sq][m % sq].plot_surface(X, Y, Z, rstride=1, cstride=1, cmap='hot', linewidth=0,
                                               antialiased=False)
        plt.show()

    return [Om, Q, {'x': Nx, 'y': Ny}, biharm]


def main():
    Lx = 1.10
    Ly = 0.8
    Lz = 5e-3
    ldim = [Lx, Ly, Lz]  # -- plate dimensions [x, y, z] in metres
    E = 9.0e+9  # -- Young's mod [Pa]
    rho = 8765  # -- density [kg/m^3]
    nu = 0.3  # -- poisson's ratio
    Nmodes = 4  # -- number of modes to compute
    h = np.sqrt(Lx * Ly) * 0.01  # -- grid spacing
    BCs = np.zeros((4, 2)) * 1e15  # -- elastic constants around the edges
    BCs[0, :] = 1e15

    matplotlib.use('macosx')
    return magpie(rho, E, nu, ldim, h, BCs, Nm=Nmodes)


class Magpie():
    """"""
    Om: np.ndarray
    "Angular Modal Frequencies"
    Q: np.ndarray
    "Eigenvector for each mode"
    N: dict
    "Dictionary Number of grid points {'x':Nx, 'y':Ny}"
    biharm: coo_matrix
    "biharmonic of the described system"
    rho: float
    "density [kg/m^3]"
    E: float
    "Young's mod [Pa]"
    nu: float
    "poisson's ratio"
    ldim: list
    "plate dimensions in meters [Lx, Ly, Lz], where Lz is thickness"
    D: float
    ""
    h: float
    "spacing between grid points of the biharmonic"
    BCs: np.ndarray
    ""
    Nmodes: int = 0,
    ""

    @classmethod
    def __init__(cls, rho: float, E: float, nu: float, ldim: list, h: float, BCs: np.ndarray):
        """

        :param rho: density [kg/m^3]
        :param E: Young's mod [Pa]
        :param nu: poisson's ratio
        :param ldim: plate dimensions in meters [Lx, Ly, Lz], where Lz is thickness
        :param h: grid spacing
        :param BCs
        """
        assert rho is not None
        assert E is not None
        assert nu is not None
        assert len(ldim) == 3
        assert h is not None
        assert BCs.shape == (4, 2)

        cls.rho = rho
        cls.E = E
        cls.nu = nu
        cls.ldim = ldim
        cls.h = h
        cls.BCs = BCs

        ##----------------------------
        Lx, Ly, Lz = ldim
        cls.D = E * (Lz ** 3) / 12 / (1 - (nu ** 2))
        ##----------------------------

    @classmethod
    def magpie(cls, Nm: int = 0,
               plot_type: str = None, base_mode: float = 0.0):
        """
        :param Nm: Number of modes, if 0 maximum number of modes are calculated
        :param plot_type: style to plot mode shapes 'chladni' or '3D'
        :return:
        """
        assert Nm is not None and Nm >= 0
        assert plot_type in ["chladni", "3D", "none"] or plot_type is None
        ##----------------------------
        Lx, Ly, Lz = cls.ldim
        Nx = int(np.ceil(Lx / cls.h))
        Ny = int(np.ceil(Ly / cls.h))
        cls.N = {'x': Nx, 'y': Ny, }
        Nmodes = ((Nx * Ny) - 2) if Nm == 0 else Nm
        cls.Nmodes = Nmodes
        ##----------------------------
        ## Build BiHarmonic
        biharm = bhmat(cls.BCs, [Nx, Ny], cls.h, Lz, cls.E, cls.nu)

        ## EIGENVALUES
        [Dm, Q] = eigs(biharm, k=Nmodes, sigma=base_mode, which='LR')

        cls.Om = np.sqrt(abs(Dm)) * np.sqrt(cls.D / cls.rho / Lz)
        indSort = np.argsort(Dm)
        cls.Q = Q[:, indSort]
        if plot_type == 'chladni':
            cls.plot_chladni()
        elif plot_type == '3D':
            cls.plot_3d()

        return [cls.Om, Q, {'x': Nx, 'y': Ny}, biharm]

    @classmethod
    def plot_chladni(cls):
        Nx = cls.N['x']
        Ny = cls.N['y']
        X = np.arange(0, Ny)
        Y = np.arange(0, Nx)
        X, Y = np.meshgrid(X, Y)

        sq = int(np.ceil(np.sqrt(cls.Nmodes)))

        for m in range(cls.Nmodes):
            ax = plt.subplot(sq, sq, m + 1)
            Z = abs(np.reshape(cls.Q[:, m], [Nx, Ny]))
            chladni = plt.pcolormesh(Z.T, cmap='copper_r', shading='gouraud')
            ax.set_axis_off()
            chladni.set_clim(0.000, 0.002)

        plt.show()

    @classmethod
    def plot_3d(cls):
        Nx = cls.N['x']
        Ny = cls.N['y']
        Q = cls.Q

        X = np.arange(0, Ny)
        Y = np.arange(0, Nx)
        X, Y = np.meshgrid(X, Y)

        sq = int(np.ceil(np.sqrt(cls.Nmodes)))

        fig, axes = plt.subplots(sq, sq, subplot_kw={"projection": "3d"})

        for m in range(cls.Nmodes):
            Z = np.reshape(Q[:, m], [Nx, Ny])
            axes[m // sq][m % sq].plot_surface(X, Y, Z, rstride=1, cstride=1, cmap='hot', linewidth=0,
                                               antialiased=False)
        plt.show()

    @classmethod
    def estimate_youngs(cls, ExpFreq: list):
        """
        get Young's modulus (E) of an experimental plate starting from a batch of experimentally measured frequencies, leveraging MAGPIE

        Example usage:

        ```
        ExpFreq = [73.2 148 376 431 559 910]   #-- these are measured from a plate
        rho     = 8765             #-- density [kg/(m ** 3)]
        Lx      = 0.1
        Ly      = 0.08
        Lz      = 0.00081

        BCs = [0,    0
               1e15, 1e15
               0,    0
               0,    0]

        ldim    = [Lx Ly Lz]
        h       = np.sqrt(Lx*Ly)*0.01   #-- grid spacing [m]


        Magpie(rho,E,nu,ldim,h,BCs)
        E  = Magpie.estimate_youngs(ExpFreq,3)
        ```

        :param rho: the experimental plate density
        :param ldim:  a 3X1 array containing the Lx Ly Lz dimensions
        :param h: the grid spacing of the FD scheme
        :param BCs: a 4X2 array containing the rigidities of the boundary supports of the experimental plate
        :param ExpFreq: an array contaning the first Nmodes measured modal frequencies
        :param Ntrain: an integer. Must be <= Nmodes. It is the number of training modes out of the available batch
        :return:
        """

        Ntrain = len(ExpFreq)
        TrainFreq = np.array(ExpFreq)

        # -- zero parameters
        Lx, Ly, Lz = cls.ldim

        A = Lx * Ly  # -- area
        D0 = cls.E * (Lz ** 3) / 12 / (1 - (cls.nu ** 2))  # -- zero-rigidity

        # -- run magpie and get non-dimensional freqs
        Om, _, _, _ = magpie(Nm=Ntrain)
        OmNDim = (Om / np.sqrt(D0)) * np.sqrt(cls.rho * (A ** 2) * Lz)
        OmNDimsq = np.reshape(OmNDim ** 2, (-1, 1)).T

        # -- least-square (LS) optimisation
        psi = np.reshape(((TrainFreq * 2 * np.pi) ** 2) * cls.rho * Lz * (A ** 2), (-1, 1))
        DLS = ((OmNDimsq @ psi) / (OmNDimsq @ OmNDimsq.T))[0, 0]
        ELS = DLS / ((Lz ** 3) / 12 / (1 - (cls.nu ** 2)))

        # -- launch a numerical simulation to get the frequencies of the numerical model
        # -- using the estimated Youngs Mod
        # -- and compare against the experimental freqs
        cls.E = ELS
        cls.magpie(Nm=Ntrain)

        # if Ntrain < Nmodes
        #
        #     subplot(2, 2, 1)
        #     Y = [TrainFreq, NumFreq(1:Ntrain)];
        #     X = 1:Ntrain;
        #     bar(X, Y)
        #     xlabel('Mode Number');
        #     ylabel('f (Hz)')
        #     legend('Exp', 'Num')
        #     title('Training Set')
        #
        #     subplot(2, 2, 3)
        #     errTrain = (1 - Y(:, 2). / Y(:, 1))*100;
        #     bar(X, errTrain)
        #     xlabel('Mode Number');
        #     ylabel('rel err (%)')
        #
        #     subplot(2, 2, 2)
        #     Y = [TestFreq, NumFreq(Ntrain + 1:end)];
        #     X = Ntrain + 1:Nmodes;
        #     bar(X, Y)
        #     xlabel('Mode Number');
        #     ylabel('f (Hz)')
        #     legend('Exp', 'Num')
        #     title('Testing Set')
        #
        #     subplot(2, 2, 4)
        #     errTest = (1 - Y(:, 2). / Y(:, 1))*100;
        #     bar(X, errTest)
        #     xlabel('Mode Number');
        #     ylabel('rel err (%)')
        #
        # else
        #
        #     subplot(2, 1, 1)
        #     Y = [TrainFreq, NumFreq(1:Ntrain)];
        #     X = 1:Ntrain;
        #     bar(X, Y)
        #     xlabel('Mode Number');
        #     ylabel('f (Hz)')
        #     legend('Exp', 'Num')
        #     title('Training Set')
        #
        #     subplot(2, 1, 2)
        #     errTrain = (1 - Y(:, 2). / Y(:, 1))*100;
        #     bar(X, errTrain)
        #     xlabel('Mode Number');
        #     ylabel('rel err (%)')
        #
        # end

        return ELS


def benchmark_eigs():
    import timeit
    BCs = np.zeros((4, 2)) * 1e15  # -- elastic constants around the edges
    BCs[0, :] = 1e15
    Lz = 1e-3
    nu = 0.3
    E = 1e9
    h = 0.01

    nbpts = 30

    results = []

    import csv

    for n in range(2, 100):
        Nx = int(10 * n)
        Ny = int(1000)
        h = 0.01

        biHarm = bhmat(BCs, [Nx, Ny], h, Lz, E, nu)

        avg_time10 = timeit.repeat(lambda: eigs(biHarm, k=k0, sigma=0.0, which="LR"),
                                   number=1, )
        print("10", end=" ")
        avg_time20 = timeit.repeat(lambda: eigs(biHarm, k=20, sigma=0.0, which="LR"),
                                   number=1)
        print("20", end=" ")
        avg_time50 = timeit.repeat(lambda: eigs(biHarm, k=50, sigma=0.0, which="LR"),
                                   number=1)
        print("30", end=" ")
        result = [Nx * Ny, min(avg_time10), min(avg_time20), min(avg_time50)]

        results.append(result)

        print(n, result)
        with open('python-eigs-timings.csv', 'a') as f:
            write = csv.writer(f)
            write.writerow(result)

    import csv

    with open('python-eigs-timings.csv', 'w') as f:
        fields = ['points', '10 modes', '20 modes', '50 modes']
        write = csv.writer(f)

        write.writerow(fields)
        write.writerows(results)


if __name__ == '__main__':
    # benchmark_eigs()

    BCs = np.zeros((4, 2))  # -- elastic constants around the edges
    Lz = 1e-3
    nu = 0.3
    E = 1e9
    h = 0.01
    Nx = int(250)
    Ny = int(250)
    h = 0.01
    k = 20
    rho = 7000
    D = E * (Lz ** 3) / 12 / (1 - (nu ** 2))

    eig_to_freq = lambda w: np.sqrt(np.abs(np.real(w))) * np.sqrt(D / rho / Lz)

    import timeit

    formats = ['bsr', 'coo', 'csc', 'csr', 'dia', 'dok', 'lil']

    for format in formats:
        biHarm = bhmat(BCs, [Nx, Ny], h, Lz, E, nu, format=format)
        print("Sigma: 0.0, LM " + format, end=': ')
        t = timeit.repeat(lambda: eigsh(biHarm,  # biharmonic
                                        k=k,  # number of eigen values
                                        M=None,  # We don't have
                                        sigma=5.566171537907070e-06,  # 0.01?
                                        which='LM',  # ‘LM’ | ‘SM’ | ‘LR’ | ‘SR’
                                        mode='normal',
                                        v0=None,  # nope
                                        ncv=None,  # nope
                                        maxiter=None,  # to test
                                        tol=0,  # to test
                                        return_eigenvectors=False,
                                        Minv=None,  # nope
                                        OPinv=None),  # nope,
                          number=1, repeat=10)
        print(min(t), "seconds")