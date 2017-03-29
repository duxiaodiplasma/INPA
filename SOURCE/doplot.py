import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

def orbit(geo,res):
    fig = plt.figure()
    ax = fig.add_subplot(111,projection='3d')

    foils = np.shape(geo.foil)[0]
    for i in range(0,foils):
        ax.plot(geo.foil[i,:,0],geo.foil[i,:,1],geo.foil[i,:,2],color='red')
    ax.scatter(geo.uphor[:,0],geo.uphor[:,1],geo.uphor[:,2],color='GREEN',s=100)
    ax.scatter(geo.lphor[:,0],geo.lphor[:,1],geo.lphor[:,2],color='black')
    ax.scatter(res.sol[::10,0],res.sol[::10,1],res.sol[::10,2])
    return 'PLOT ORBIT'


def sp(geo,res):
    fig = plt.figure()
    ax = fig.add_subplot(111,projection='3d')

    ax.scatter(geo.uphor[:,0],geo.uphor[:,1],geo.uphor[:,2],color='GREEN',s=100)
    ax.scatter(geo.lphor[:,0],geo.lphor[:,1],geo.lphor[:,2],color='black')

    sp = np.copy(res.hitpoint)
    mask = np.where(res.hitpoint[:,0] != 0)
    #sc = ax.scatter(sp[mask,0],sp[mask,1],sp[mask,2],c=res.R_birth[:,0],alpha=0.5)
    sc = ax.scatter(sp[mask,0],sp[mask,1],sp[mask,2],c=res.E_birth,alpha=0.5)
    cbar = plt.colorbar(sc)


    foils = np.shape(geo.foil)[0]
    for i in range(0,foils):
        ax.plot(geo.foil[i,:,0],geo.foil[i,:,1],geo.foil[i,:,2],color='red')

    return 'PLOT STRIKING POINTS'

def reso(ini,res):
    fig = plt.figure()
    levels = np.linspace(0,1,res.segs)
    cs = plt.contour(ini.r,ini.z,np.sqrt(np.transpose(ini.psirz)),levels,vmin=0,vmax=1)
    cb = plt.colorbar(cs, shrink=0.8, extend='both',ticks=levels[0::2])


    plt.clabel(cs, levels[0::2], inline=1, fontsize=10)
    for i in np.arange(ini.mc)[::10]:
        plt.plot(res.rzseg[i,:,0], res.rzseg[i,:,1])


    #f, axarr = plt.subplots(2)
    #levels = np.linspace(0,1,21)
    #cs = axarr[0].contour(rr,zz,np.sqrt(np.transpose(psirz)),levels)
    ##cb = colorbar(cs, shrink=0.8, extend='both',ticks=levels[0::2])
    #
    #
    ##plt.clabel(cs, levels[0::2], inline=1, fontsize=10)
    #for i in np.arange(300):
    #    axarr[0].plot(rz_seg[i,:,0], rz_seg[i,:,1])
    #    axarr[1].plot(rz_seg[i,:,0],pitch[i,:]/np.pi*180)
