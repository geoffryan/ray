import h5py as h5
import numpy as np

sphMetrics = [1,2,3]
cartMetrics = [0,4,5]

class Track:

    id = -1
    t = None
    x = None
    u = None

    def __init__(self, id, T, X, U):
        self.id = id
        self.t = T
        self.x = X
        self.u = U

class Map:

    thetaC = None
    phiC = None
    T0 = None
    X0 = None
    U0 = None
    T = None
    X = None
    U = None
    nhits = None

    def __init__(self, thetaC, phiC, T0, X0, U0, T, X, U, nhits):
        self.thetaC = thetaC
        self.phiC = phiC
        self.T0 = T0
        self.X0 = X0
        self.U0 = U0
        self.T = T
        self.X = X
        self.U = U
        self.nhits = nhits


def loadPars(filename):

    pars = {}
    f = h5.File(filename, "r")
    parG = f['Pars']

    for key in parG:
        val = parG[key]
        pars[key] = parG[key][0]

    return pars

def loadMap(filename):

    f = h5.File(filename, "r")
    map = f['Map']
    T0 = map['t0'][...]
    X0 = map['x0'][...]
    U0 = map['u0'][...]
    T = map['t'][...]
    X = map['x'][...]
    U = map['u'][...]
    skyLoc = map['thC'][...]
    nhits = map['nhits'][...]
    f.close()

    thetaC = skyLoc[:,0]
    phiC = skyLoc[:,1]

    map = Map(thetaC, phiC, T0, X0, U0, T, X, U, nhits)

    return map

def loadTracks(filename):
    
    tracks = []

    f = h5.File(filename, "r")
    trackG = f['Tracks']

    for key in trackG:
        id = int(key[7:])
        t = trackG[key]['t'][...]
        x = trackG[key]['x'][...]
        u = trackG[key]['u'][...]

        tracks.append(Track(id, t, x, u))

    f.close()

    return tracks

def getCartesianCoords(X, pars):

    if pars['Metric'] in cartMetrics:
        return X[...,1].copy(), X[...,2].copy(), X[...,3].copy()
    else:
        r = X[...,1]
        cost = np.cos(X[...,2])
        sint = np.sin(X[...,2])
        cosp = np.cos(X[...,3])
        sinp = np.sin(X[...,3])

        return r*sint*cosp, r*sint*sinp, r*cost

def getSphericalCoords(X, pars):

    if pars['Metric'] in sphMetrics:
        return X[...,1].copy(), X[...,2].copy(), X[...,3].copy()
    else:
        x = X[...,1]
        y = X[...,2]
        z = X[...,3]
        r = np.sqrt(x*x + y*y)
        R = np.sqrt(r*r + z*z)
        phi = np.arctan2(y, x)
        theta = np.arccos(z/R)

        return R, theta, phi

