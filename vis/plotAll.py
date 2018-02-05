import sys
import math
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import rayUtil as ru

def intensity_face(X, U, pars, R=30.0):

    x0 = X[:,0]
    x1 = X[:,1]
    x2 = X[:,2]
    x3 = X[:,3]

    x, y, z = ru.getCartesianCoords(X, pars)
    r, theta, phi = ru.getSphericalCoords(X, pars)

    F = np.zeros(x.shape)

    offsurface = np.fabs(z) > 1.0e-1

    F[x*x+y*y < R*R] = 1.0

    left_eye = (x+R/np.sqrt(8))**2 + (y-R/np.sqrt(8))**2 < (0.25*R)**2
    right_eye = (np.fabs(x-R/np.sqrt(8)) < 0.25*R) * (
                    np.fabs(y-R/np.sqrt(8)) < 0.1*R)
    mouth = (np.fabs(r*np.sin(theta) - 0.6*R) < 0.1*R) * (y<0)

    F[left_eye] = 0.5
    F[right_eye] = 0.5
    F[mouth] = 0.5

    F[offsurface] = 0.0

    return F


def intensity(X, U, pars):

    F = np.zeros(X.shape[0])

    r, theta, phi = ru.getSphericalCoords(X, pars)

    surface = np.fabs(theta-0.5*np.pi) < 0.1

    F[surface] = 1.0 + 0.5*np.cos(phi)
    return F

def plotMap(ax, filename):

    map = ru.loadMap(filename)
    pars = ru.loadPars(filename)

    t1 = map.T[:,0]
    X1 = map.X[:,0,:]
    U1 = map.U[:,0,:]

    F = intensity_face(X1, U1, pars)
    #F = intensity(X1, U1, pars)

    print(F.min(), F.max())

    ax.tricontourf(map.phiC, map.thetaC, F, 256, cmap=mpl.cm.inferno)
    ax.set_xlim(map.phiC.max(), map.phiC.min())
    ax.set_ylim(map.thetaC.max(), map.thetaC.min())
    ax.set_aspect('equal')

def plotAll(filename):

    fig, ax = plt.subplots(2,2, figsize=(12,9))

    tracks = ru.loadTracks(filename)
    
    for track in tracks:
        x0 = track.x[:,0]
        x1 = track.x[:,1]
        x2 = track.x[:,2]
        x3 = track.x[:,3]

        ax[0,0].plot(x1, x2, 'k', alpha=0.3)
        ax[0,1].plot(x3, x2, 'k', alpha=0.3)
        ax[1,0].plot(x1, x3, 'k', alpha=0.3)

    ax[0,0].set_xlabel(r'$x^1$')
    ax[0,0].set_ylabel(r'$x^2$')
    ax[0,1].set_xlabel(r'$x^3$')
    ax[0,1].set_ylabel(r'$x^2$')
    ax[1,0].set_xlabel(r'$x^1$')
    ax[1,0].set_ylabel(r'$x^3$')

    plotMap(ax[1,1], filename)

    fig.tight_layout()
    fig.savefig("plot.png")
    plt.close(fig)


if __name__ == "__main__":

    if(len(sys.argv) < 2):
        print("Need a map file")

    plotAll(sys.argv[1])
