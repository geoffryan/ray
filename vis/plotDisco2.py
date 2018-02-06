import sys
import math
from itertools import imap, izip
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import h5py as h5

class Map:
    t = None
    thC = None
    X0 = None
    U0 = None
    X1 = None
    U1 = None

    def __init__(self, mapfile):
        f = h5.File(mapfile, "r")
        self.t = f["Map/t"][...]
        self.thC = f["Map/thC"][...]
        self.X0 = f["Map/x0"][...]
        self.U0 = f["Map/u0"][...]
        self.X1 = f["Map/x1"][...]
        self.U1 = f["Map/u1"][...]
        f.close()

    def scale(self, M):
        self.t *= M
        self.X0[:,0] *= M
        self.X0[:,1] *= M
        self.X1[:,0] *= M
        self.X1[:,1] *= M
        self.U0[:,2] *= M
        self.U0[:,3] *= M
        self.U1[:,2] *= M
        self.U1[:,3] *= M

def readCheckpoint(filename, nq=-1):
    """
    Return the data from a checkpoint file.  Remove duplicate entries.
    """

    f = h5.File(filename, "r")
    Data = f['Data'][...]
    gravMass = f['GravMass'][...]
    t = f['T'][0]
    f.close()

    #Standard way of removing duplicates (thnks to dfm and dh)
    Data = np.array(dict(izip(imap(hash, imap(tuple,Data[:,0:3])), Data)).values())
    
    #That step messed up the data, sort it by z, r, phi
    sort_inds = np.lexsort((Data[:,0], Data[:,1], Data[:,2]))
    Data = Data[sort_inds]

    #read in coords
    piph = np.array(Data[:,0])
    r = np.array(Data[:,1])
    z = np.array(Data[:,2])

    #read in data
    rho = np.array(Data[:,3])
    P = np.array(Data[:,4])
    vr = np.array(Data[:,5])
    vp = np.array(Data[:,6])
    vz = np.array(Data[:,7])
    q = []

    if nq < 1 or nq > Data.shape[1]:
        nq = Data.shape[1]-3
        
    for i in xrange(nq-5):
        try:
            q.append(np.array(Data[:,8+i]))
        except:
            q.append(np.zeros(r.shape))
    q = tuple(q)

    #calculate cell volumes
    z_vals = np.unique(z)
    r_vals = np.unique(r)
    
    phi = np.ones(len(r))
    dV = np.ones(len(r))
    dr = np.ones(len(r))
    dphi = np.ones(len(r))

    #dz
    if len(z_vals) > 1:
        dV[:] = (z_vals[len(z_vals)-1] - z_vals[0]) / (len(z_vals) - 1.0)

    #dr
    dr_vals = np.ones(r_vals.shape)
    Rm = 0.5*(r_vals[-2]+r_vals[-1])
    Rp = 1.5*r_vals[-1] - 0.5*r_vals[-2]
    for i in range(len(r_vals)-1, -1, -1):
        Rm = 2*r_vals[i] - Rp
        dr_vals[i] = Rp - Rm
        Rp = Rm

    for i in range(len(r_vals)):
        dV[r==r_vals[i]] *= dr_vals[i]
        dr[r==r_vals[i]] = dr_vals[i]

    #dphi
    for i in range(len(r_vals)):
        inds = (r==r_vals[i])
        my_phi = piph[inds]
        dp = np.zeros(my_phi.shape)
        dp[1:] = my_phi[1:] - my_phi[:-1]
        dp[0] = my_phi[0] - my_phi[-1]
        dp[dp>2*math.pi] -= 2*math.pi
        dp[dp<0] += 2*math.pi
        #dp = 2*math.pi/len(my_phi)
        dphi[inds] = dp
        phi[inds] = piph[inds] - 0.5*dp
        dV[inds] *= r[inds] * dphi[inds]

    return t, r, phi, z, rho, P, vr, vp, vz, dV, q, piph, dphi, gravMass

def imageDiskFastLight(mapFile, checkpoints, nu, M=1.0):

    map = Map(mapFile)

    map.scale(M)

    for i, checkpoint in enumerate(checkpoints):
        dat = readCheckpoint(checkpoint)
        t, flux = makeImage(map, dat, i, nu, M=M)

def makeImage(map, dat, id, nu, M=1.0):

    t = dat[0]
    r = dat[1]
    phi = dat[2]
    rho = dat[3]
    P = dat[4]
    vr = dat[5]
    vp = dat[6]
    piph = dat[7]
    dphi = dat[8]

    R = np.unique(r)
    rjph = np.empty(len(R)+1)
    rjph[0] = 0.0
    rjph[1:-1] = 0.5*(R[1:]+R[:-1])
    rjph[-1] = R[-1] + R[-1]-rjph[-2]

    thC = map.thC
    X = map.X1
    U = map.U1
    N = thC.shape[0]

    T = np.empty(N)
    zp1 = np.empty(N)

    for n in range(N):
        j = np.searchsorted(rjph, X[n,1]) - 1
        ind = (r==R[j])
        piph_loc = piph[ind]
        sind = np.argsort(piph_loc)
        piph_loc = piph_loc[sind]
        #dphi_loc = dphi[ind]
        dphi_loc = piph_loc - np.roll(piph_loc, 1)
        print(piph_loc)
        while (dphi_loc > 2*np.pi).any():
            dphi_loc[dphi_loc>2*np.pi] -= 2*np.pi
        while (dphi_loc < 0.0).any():
            dphi_loc[dphi_loc<0.0] += 2*np.pi
        phi_loc = piph_loc - 0.5*dphi_loc
        nphi = piph_loc.shape[0]
        i=-1
        for ii in range(nphi):
            diff = X[n,3] - phi_loc[ii]
            while diff > np.pi:
                diff -= 2*np.pi
            while diff < -np.pi:
                diff += 2*np.pi
            print("   {0} {1} {2}".format(diff, phi_loc[ii], dphi_loc[ii]))
            if diff <= 0.5*dphi_loc[ii] and diff >= -0.5*dphi_loc[ii]:
                i = ii
                break
        print(j, i, X[n])

        T[n] = P[ind][i] / rho[ind][i]
        if rho[ind][i] == 0.0:
            T[n] = 0.0
        u = calc_u(R[j], phi_loc[i], vr[ind][i], vp[ind][i], M)
        zp1[n] = -U[n,0]*u[0]-U[n,1]*u[1]-U[n,3]*u[2]

    print((T.min(), T.max()))

    Inu = (nu*nu*nu)[None,:] / (np.exp(-nu[None,:]/T[:,None])-1.0)

    Inu[-np.isfinite(Inu)] = 0.0

    for i in range(len(nu)):
        print(i)
        fig, ax = plt.subplots(1,1, figsize=(12,9))
        ax.tricontourf(thC[:,1], thC[:,0], Inu[:,i], 256, cmap=mpl.cm.inferno)
        ax.set_xlim(thC[:,1].max(), thC[:,1].min())
        ax.set_ylim(thC[:,0].max(), thC[:,0].min())
        fig.savefig("image_{0:03d}_{1:02d}.png".format(id, i))
        plt.close(fig)

    return t, Inu.sum()

def calc_u(r, phi, vr, vp, M):

    try:
        u0 = 1.0/math.sqrt(1.0-2*M/r-4*M/r*vr - (1.0+2*M/r)*vr*vr - r*r*vp*vp)
        ur = u0*vr
        up = u0*vp
    except:
        u0 = 1.0
        ur = 0.0
        up = 0.0

    return np.array([u0, ur, up])


if __name__ == "__main__":

    mapFile = sys.argv[1]
    checkpoints = sys.argv[2:]

    M = 1.0/30.0
    nu = np.array([0.01, 0.03, 0.1, 0.3, 1.0])

    imageDiskFastLight(mapFile, checkpoints, nu, M)
