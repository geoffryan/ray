import sys
import math
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import h5py as h5

def kep_u(r, M, a):
    return 0

def intensity_disk_T(X, U, T=1.0, M=1.0):

    x0 = X[:,0]
    x1 = X[:,1]
    x2 = X[:,2]
    x3 = X[:,3]

    F = np.zeros(x1.shape)
    r = x1
    phi = x3

    surface = np.fabs(x2-0.5*np.pi) < 0.1

    F[surface] = 1.0 + 0.5*np.cos(x3[surface])
    return F

def intensity_face(X, U, R=30.0):

    x0 = X[:,0]
    x1 = X[:,1]
    x2 = X[:,2]
    x3 = X[:,3]

    F = np.zeros(x1.shape)

    r = x1
    theta = x2
    phi = x3
    x = r*np.sin(theta)*np.cos(phi)
    y = r*np.sin(theta)*np.sin(phi)
    z = r*np.cos(theta)
    offsurface = np.fabs(theta-0.5*np.pi) > 1.0e-6

    #x = x1
    #y = x2
    #z = x3
    #r = np.sqrt(x*x + y*y + z*z)
    #theta = np.arccos(z/r)
    #phi = np.arctan2(y, x)
    #offsurface = np.fabs(z) > 1.0e-1

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


def intensity(X, U):

    F = np.zeros(X.shape[0])

    x0 = X[:,0]
    x1 = X[:,1]
    x2 = X[:,2]
    x3 = X[:,3]

    surface = np.fabs(x2-0.5*np.pi) < 0.1

    F[surface] = 1.0 + 0.5*np.cos(x3[surface])
    return F

def plotMap(ax, mapFilename):

    f = h5.File(mapFilename, "r")
    skyLoc = f["Map/thC"][...]
    t = f["Map/t"][...]
    X0 = f["Map/x0"][...]
    U0 = f["Map/u0"][...]
    X1 = f["Map/x1"][...]
    U1 = f["Map/u1"][...]
    f.close()

    thC = skyLoc[:,0]
    phC = skyLoc[:,1]

    F = intensity_face(X1, U1)
    #F = intensity(X1, U1)

    print(F.min(), F.max())

    ax.tricontourf(phC, thC, F, 256, cmap=mpl.cm.inferno)
    ax.set_xlim(phC.max(), phC.min())
    ax.set_ylim(thC.max(), thC.min())
    ax.set_aspect('equal')

def plotMapNice(mapFilename):

    f = h5.File(mapFilename, "r")
    skyLoc = f["Map/thC"][...]
    t = f["Map/t"][...]
    X0 = f["Map/x0"][...]
    U0 = f["Map/u0"][...]
    X1 = f["Map/x1"][...]
    U1 = f["Map/u1"][...]
    f.close()

    thC = skyLoc[:,0]
    phC = skyLoc[:,1]

    F = intensity(X1, U1)

    lat = 0.5*np.pi - thC
    lon = phC.copy()
    while (lon > np.pi).any():
        lon[lon>np.pi] -= 2*np.pi
    while (lon < -np.pi).any():
        lon[lon<-np.pi] += 2*np.pi

    #lon *= 180.0/np.pi
    #lat *= 180.0/np.pi

    print(lat.min(), lat.max())
    print(lon.min(), lon.max())
    print(F.min(), F.max())

    fig = plt.figure()
    ax = fig.add_subplot(111, projection='hammer')

    numt = len(thC)

    ax.tricontourf(lat, lon, F, 256, cmap=mpl.cm.inferno)

    fig.savefig('map.png')
    plt.close(fig)

def plotAll(mapFilename, trackFilenames):

    fig, ax = plt.subplots(2,2, figsize=(12,9))

    for trackFilename in trackFilenames:
        t, x0, x1, x2, x3 = np.loadtxt(trackFilename, usecols=[0,1,2,3,4],
                                        unpack=True)

        ax[0,0].plot(x1, x2, 'k', alpha=0.3)
        ax[0,1].plot(x3, x2, 'k', alpha=0.3)
        ax[1,0].plot(x1, x3, 'k', alpha=0.3)
        #x = x1*np.sin(x2)*np.cos(x3)
        #y = x1*np.sin(x2)*np.sin(x3)
        #z = x1*np.cos(x2)
        #ax[0,0].plot(x, y, 'k', alpha=0.3)
        #ax[0,1].plot(z, y, 'k', alpha=0.3)
        #ax[1,0].plot(x, z, 'k', alpha=0.3)

    ax[0,0].set_xlabel(r'$x^1$')
    ax[0,0].set_ylabel(r'$x^2$')
    ax[0,1].set_xlabel(r'$x^3$')
    ax[0,1].set_ylabel(r'$x^2$')
    ax[1,0].set_xlabel(r'$x^1$')
    ax[1,0].set_ylabel(r'$x^3$')

    plotMap(ax[1,1], mapFilename)

    fig.tight_layout()
    fig.savefig("plot.png")
    plt.close(fig)

if __name__ == "__main__":

    if(len(sys.argv) < 2):
        print("Need a map and tracks")

    plotAll(sys.argv[1], sys.argv[2:])
    #plotMapNice(sys.argv[1])
