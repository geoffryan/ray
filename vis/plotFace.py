import sys
import numpy as np
import matplotlib.pyplot as plt
import rayUtil as ru
import rayPlot as rp

def plotFace(filename, imageOnly=False):

    map = ru.loadMap(filename)
    pars = ru.loadPars(filename)

    t = map.T[:,0]
    X1 = map.X[:,0,:]
    U1 = map.U[:,0,:]

    I = rp.intensity_face(t, X1, U1, pars)

    name = ".".join(filename.split(".")[:-1])
    plotname = "plot_{0:s}_face.png".format(name)

    figw = 12.0

    if imageOnly:
        figh = (figw * (pars['X1b']-pars['X1a']))/(pars['X2b']-pars['X2a'])
        fig = plt.figure(figsize=(figw,figh))
        ax = fig.add_axes([0,0,1,1])
        rp.plotMap(ax, I, map, pars, axes=False, colorbar=False)
    else:
        figh = 0.75 * figw
        fig, ax = plt.subplots(1,1,figsize=(figw,figh))
        rp.plotMap(ax, I, map, pars, colorbar=True)
        fig.tight_layout()

    fig.savefig(plotname)
    plt.close(fig)



if __name__ == "__main__":

    if len(sys.argv) < 2:
        print("Please give me Ray hdf5 files!")

    for filename in sys.argv[1:]:
        plotFace(filename, False)


