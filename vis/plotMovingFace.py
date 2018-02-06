import sys
import numpy as np
import matplotlib.pyplot as plt
import rayUtil as ru
import rayPlot as rp

def plotFaces(filename, ts, v):

    map = ru.loadMap(filename)
    pars = ru.loadPars(filename)

    t = map.T[:,0]
    X1 = map.X[:,0,:]
    U1 = map.U[:,0,:]

    name = ".".join(filename.split(".")[:-1])
    dpi = 100
    figw = 4.0
    figh = figw * (pars['X1b']-pars['X1a'])/(pars['X2b']-pars['X2a'])

    if int(figh*dpi) % 2 != 0:
        figh = float(int(figh*dpi)+1) / dpi

    for i, ti in enumerate(ts):
        print("frame {0:d} ({1:d})".format(i, len(ts)))

        I = rp.intensity_face(t, X1, U1, pars, R=30.0, tnow=ti, v=v,
                                fastlight=False)

        plotname = "plot_{0:s}_face_{1:04d}.png".format(name, i)

        fig = plt.figure(figsize=(figw,figh))
        ax = fig.add_axes([0,0,1,1])
        rp.plotMap(ax, I, map, pars, axes=False, colorbar=False)

        fig.savefig(plotname)
        plt.close(fig)


if __name__ == "__main__":

    if len(sys.argv) < 2:
        print("Please give me Ray hdf5 files!")

    t = np.linspace(-60, 60, 121)
    v = np.array([1.0,0,0])

    for filename in sys.argv[1:]:
        plotFaces(filename, t, v)


