import numpy as np
from bhmat import bhmat
from scipy.sparse.linalg import eigs
from scipy.sparse import eye
import matplotlib
from matplotlib import pyplot as plt


def magpie(rho: float, E: float, nu: float, ldim: list, h: float, BCs: np.ndarray, Nmodes: int, plot_type: str=None, base_mode: float=0.0):
    """

    :param rho:
    :param E:
    :param nu:
    :param ldim:
    :param h:
    :param BCs:
    :param Nmodes:
    :param plot_type:
    :return:
    """
    ##----------------------------
    Lx, Ly, Lz = ldim
    D = E * (Lz ** 3) / 12 / (1 - (nu ** 2))
    Nx = int(np.ceil(Lx / h))
    Ny = int(np.ceil(Ly / h))
    ##----------------------------
    ## Build BiHarmonic
    biharm = bhmat(BCs, [Nx, Ny], h, Lz, E, nu)

    ## EIGENVALUES
    [Dm, Q] = eigs(biharm, Nmodes, sigma=base_mode, which='LR')

    Om = np.sqrt(abs(Dm)) * np.sqrt(D / rho / Lz) 
    hz = Om / (2 * np.pi)
    indSort = np.argsort(Dm)
    Q = Q[:, indSort]

    X = np.arange(0, Ny)
    Y = np.arange(0, Nx)
    X, Y = np.meshgrid(X, Y)

    sq = int(np.ceil(np.sqrt(Nmodes)))

    if plot_type == 'chladni':
        for m in range(Nmodes):
            ax = plt.subplot(sq, sq, m + 1)
            Z = abs(np.reshape(Q[:, m], [Nx, Ny]))
            chladni = plt.pcolormesh( Z.T, cmap='copper_r', shading='gouraud')
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

    return [Q, Om, {'x': Nx,'y': Ny}, biharm]

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
    BCs[0,:] = 1e15

    matplotlib.use('macosx')
    return magpie(rho, E, nu, ldim, h, BCs, Nmodes, 'chladni')

if __name__ == '__main__':

    Q, Om, N, biharm = main()
    # Nmodes = Q.shape[1]
    # sq = int(np.ceil(np.sqrt(Nmodes)))
    #
    # X = np.arange(0, N['y'])
    # Y = np.arange(0, N['x'])
    # X, Y = np.meshgrid(X, Y)
    # for m in range(Nmodes):
    #     plt.subplot(sq, sq, m+1)
    #     Z = np.real(np.reshape(Q[:, m], [N['x'], N['y']]))
    #     plt.pcolormesh(X,Y, Z, cmap='copper', shading='gouraud')
    # plt.show()
