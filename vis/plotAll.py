import sys
import math
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import rayUtil as ru
import rayPlot as rp

def intensity(t, X, U, pars):

    F = np.zeros(X.shape[0])

    r, theta, phi = ru.getSphericalCoords(X, pars)

    surface = np.fabs(theta-0.5*np.pi) < 0.1

    F[surface] = 1.0 + 0.5*np.cos(phi)
    return F


def plotAll(filename):

    fig, ax = plt.subplots(2,3, figsize=(16,9))

    tracks = ru.loadTracks(filename)
    map = ru.loadMap(filename)
    pars = ru.loadPars(filename)
    
    for track in tracks:
        t = track.t
        x0 = track.x[:,0]
        x1 = track.x[:,1]
        x2 = track.x[:,2]
        x3 = track.x[:,3]
        u0 = track.u[:,0]
        u3 = track.u[:,3]

        ax[0,0].plot(x1, x2, 'k', alpha=0.3)
        ax[0,1].plot(x3, x2, 'k', alpha=0.3)
        ax[1,0].plot(x1, x3, 'k', alpha=0.3)
        ax[0,2].plot(t, u0, 'k', alpha=0.3)
        ax[1,2].plot(t, u3, 'k', alpha=0.3)

    ax[0,0].set_xlabel(r'$x^1$')
    ax[0,0].set_ylabel(r'$x^2$')
    ax[0,1].set_xlabel(r'$x^3$')
    ax[0,1].set_ylabel(r'$x^2$')
    ax[1,0].set_xlabel(r'$x^1$')
    ax[1,0].set_ylabel(r'$x^3$')
    ax[0,2].set_xlabel(r'$t$')
    ax[0,2].set_ylabel(r'$u_0$')
    ax[1,2].set_xlabel(r'$t$')
    ax[1,2].set_ylabel(r'$u_3$')


    t1 = map.T[:,0]
    X1 = map.X[:,0,:]
    U1 = map.U[:,0,:]

    I = rp.intensity_face(t1, X1, U1, pars, R=10000.0)
    #I = intensity(t1, X1, U1, pars)

    print(I.min(), I.max())
    rp.plotMapBasic(ax[1,1], I, map, pars)

    fig.tight_layout()
    fig.savefig("plot.png")
    plt.close(fig)


if __name__ == "__main__":

    if(len(sys.argv) < 2):
        print("Need a map file")

    plotAll(sys.argv[1])
