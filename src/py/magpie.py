import numpy as np
from bhmat import bhmat
from scipy.sparse.linalg import eigs
import matplotlib
from matplotlib import pyplot as plt
matplotlib.use('macosx')

def magpie(rho: float, E: float, nu: float, ldim: list, h: float, BCs: np.ndarray, Nmodes: int, plot_type: str):
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

    [Dm, Q] = eigs(biharm, Nmodes, which='SM')

    Om = np.sqrt(abs(Dm)) * np.sqrt(D / rho / Lz) / 2 / np.pi
    hz = Om / (2 * np.pi)
    indSort = np.argsort(Dm)
    Q = Q[:, indSort]



    m = 1;
    X = np.arange(0, Ny)
    Y = np.arange(0, Nx)
    X, Y = np.meshgrid(X, Y)

    sq = int(np.ceil(np.sqrt(Nmodes)))
    fig , axes = plt.subplots(sq, sq, sharey=True, subplot_kw={"projection": "3d"})

    for m in range(Nmodes):
        Z = np.reshape(Q[:, m], [Nx, Ny])
        axes[m//sq][m%sq].plot_surface(X, Y, Z, rstride=1, cstride=1, cmap='hot', linewidth=0, antialiased=False)

    plt.show()
    # if plot_type == 'chladni':
    #     subs = ceil(sqrt(Nmodes));
    #
    #     colormap('copper');
    #     cmp = colormap;
    #     cmp = flipud(cmp);
    #     colormap(cmp);
    #
    #     for m in np.r_[1: Nmodes]:
    #         mdShape = reshape(Q(:, m), [(Ny + 1), (Nx + 1)]);
    #         subplot(subs, subs, m)
    #         mesh(3e3 * real(mdShape), (abs(mdShape)), 'FaceColor', 'texturemap');
    #         view(2);
    #         axis
    #         equal;
    #         axis
    #         tight;
    #         axis
    #         off;
    #         clim([0.00005 0.002]);
    #
    #
    # elif plot_type == '3D':
    #
    #     subs = ceil(sqrt(Nmodes));
    #     colormap('parula');
    #
    #     for m in np.r_[1: Nmodes]:
    #         mdShape = reshape(Q(:, m), [(Ny + 1), (Nx + 1)]);
    #         subplot(subs, subs, m)
    #         mesh(3000 * (mdShape), (abs(mdShape)), 'FaceColor', 'texturemap');

def main():
    Lx = 1.10
    Ly = 0.8
    Lz = 5e-3
    ldim = [Lx, Ly, Lz]  # -- plate dimensions [x, y, z] in metres
    E = 9.0e+9  # -- Young's mod [Pa]
    rho = 8765  # -- density [kg/m^3]
    nu = 0.3  # -- poisson's ratio
    Nmodes = 16  # -- number of modes to compute
    h = np.sqrt(Lx * Ly) * 0.05  # -- grid spacing
    BCs = np.zeros((4, 2)) * 1e15  # -- elastic constants around the edges
    BCs[0,:] = 1e15

    magpie(rho, E, nu, ldim, h, BCs, Nmodes, 'chladni')

if __name__ == '__main__':
    main()
