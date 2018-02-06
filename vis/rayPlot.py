import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import rayUtil as ru

def plotMapBasic(ax, I, map, pars):


    if pars['MeshType'] == 0:
        Nt = pars['N1']
        Np = pars['N2']
        tmin = pars['X1a']
        tmax = pars['X1b']
        pmin = pars['X2a']
        pmax = pars['X2b']
        theta = np.linspace(tmin, tmax, Nt+1)
        phi = np.linspace(pmin, pmax, Np+1)
        pp, tt = np.meshgrid(phi, theta)
        II = I.reshape((Nt, Np))
        ax.pcolormesh(pp, tt, II, cmap=mpl.cm.inferno)
        ax.set_xlim(pmax, pmin)
        ax.set_ylim(tmax, tmin)
    else:
        ax.tricontourf(map.phiC, map.thetaC, I, 256, cmap=mpl.cm.inferno)
        ax.set_xlim(map.phiC.max(), map.phiC.min())
        ax.set_ylim(map.thetaC.max(), map.thetaC.min())
    ax.set_aspect('equal')


def plotMap(ax, I, map, pars, axes=True, label=True, labelSize=12, 
                tickLabelSize=12, cmap=mpl.cm.inferno, colorbar=False, 
                colorbarLabel=None, colorbarLabelSize=12):

    if pars['MeshType'] == 0:
        Nt = pars['N1']
        Np = pars['N2']
        tmin = pars['X1a']
        tmax = pars['X1b']
        pmin = pars['X2a']
        pmax = pars['X2b']
        theta = np.linspace(tmin, tmax, Nt+1)
        phi = np.linspace(pmin, pmax, Np+1)
        pp, tt = np.meshgrid(phi, theta)
        II = I.reshape((Nt, Np))
        C = ax.pcolormesh(pp, tt, II, cmap=mpl.cm.inferno)
        ax.set_xlim(pmax, pmin)
        ax.set_ylim(tmax, tmin)
    else:
        C = ax.tricontourf(map.phiC, map.thetaC, I, 256, cmap=mpl.cm.inferno)
        ax.set_xlim(map.phiC.max(), map.phiC.min())
        ax.set_ylim(map.thetaC.max(), map.thetaC.min())
    ax.set_aspect('equal')

    if axes:
        if label:
            if label is True:
                ax.set_xlabel(r"$\phi_C$", fontsize=labelSize)
                ax.set_ylabel(r"$\theta_C$", fontsize=labelSize)
            else:
                ax.set_xlabel(label[0], fontsize=labelSize)
                ax.set_ylabel(label[1], fontsize=labelSize)
    else:
        ax.axis('off')

    if colorbar:
        fig = ax.get_figure()
        cb = fig.colorbar(C, ax=ax)
        if colorbarLabel is not None:
            cb.set_label(colorbarLabel, fontsize=colorbarLabelSize)


def intensity_face(t, X, U, pars, R=30.0, t0=0.0, tnow=0.0, center=[0,0,0],
                    v=[0,0,0], fastlight=True):

    x0 = X[:,0]
    x1 = X[:,1]
    x2 = X[:,2]
    x3 = X[:,3]

    x, y, z = ru.getCartesianCoords(X, pars)
    if fastlight:
        x -= center[0] + v[0]*(tnow-t0)
        y -= center[1] + v[1]*(tnow-t0)
        z -= center[2] + v[2]*(tnow-t0)
    else:
        tave = x0.mean()
        x -= center[0] + v[0]*(tnow-t0+x0-tave)
        y -= center[1] + v[1]*(tnow-t0+x0-tave)
        z -= center[2] + v[2]*(tnow-t0+x0-tave)

    r = np.sqrt(x*x + y*y)
    #r, theta, phi = ru.getSphericalCoords(X, pars)

    F = np.zeros(x.shape)

    offsurface = np.fabs(z) > 1.0e-1

    F[r*r < R*R] = 1.0

    left_eye = (x+R/np.sqrt(8))**2 + (y-R/np.sqrt(8))**2 < (0.25*R)**2
    right_eye = (np.fabs(x-R/np.sqrt(8)) < 0.25*R) * (
                    np.fabs(y-R/np.sqrt(8)) < 0.1*R)
    mouth = (np.fabs(r - 0.6*R) < 0.1*R) * (y<0)

    F[left_eye] = 0.5
    F[right_eye] = 0.5
    F[mouth] = 0.5

    F[offsurface] = 0.0

    return F

