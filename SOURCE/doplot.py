import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D


def sp(geo,res):
    fig = plt.figure()
    ax = fig.add_subplot(111,projection='3d')

    #ax.scatter(s[::100,0],s[::100,1],s[::100,2])
    ax.scatter(geo.uphor[:,0],geo.uphor[:,1],geo.uphor[:,2],color='GREEN',s=100)
    ax.scatter(geo.lphor[:,0],geo.lphor[:,1],geo.lphor[:,2],color='black')

    sp = np.copy(res.hitpoint)
    mask = np.where(res.hitpoint[:,0] != 0)
    sc = ax.scatter(sp[mask,0],sp[mask,1],sp[mask,2],c=res.E_birth,alpha=0.5)
    cbar = plt.colorbar(sc)


    foils = np.shape(geo.foil)[0]
    for i in range(0,foils):
        ax.scatter(geo.foil[i,:,0],geo.foil[i,:,1],geo.foil[i,:,2],color='red')

