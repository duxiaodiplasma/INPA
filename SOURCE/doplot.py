import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

def orbit3d(geo,res):
    fig = plt.figure()
    ax = fig.add_subplot(111,projection='3d')

    foils = np.shape(geo.foil)[0]
    for i in range(0,foils):
        ax.plot(geo.foil[i,:,0],geo.foil[i,:,1],geo.foil[i,:,2],color='red')
    ax.scatter(geo.uphor[:,0],geo.uphor[:,1],geo.uphor[:,2],color='GREEN',s=res.bins)
    ax.scatter(geo.lphor[:,0],geo.lphor[:,1],geo.lphor[:,2],color='black')
    ax.scatter(res.sol[::10,0],res.sol[::10,1],res.sol[::10,2])
    return 'PLOT ORBIT'

def sp2d(geo,res):
    plt.clf()

    #sc = plt.scatter(res.rot_hitpoint[:,0],res.rot_hitpoint[:,2],c=res.E_birth,alpha=0.5)
    sc = plt.scatter(res.rot_hitpoint[:,0],res.rot_hitpoint[:,2],c=res.rho,alpha=0.5)
    plt.plot(geo.rot_foil[0,:,0],geo.rot_foil[0,:,2],alpha=0.5)
    plt.plot(geo.rot_lphor[:,0],geo.rot_lphor[:,2],alpha=0.5)
    cb = plt.colorbar(sc)
    return 'PLOT STRIKING POINTS'


def sp3d(geo,res):
    plt.clf()
    fig = plt.figure()
    ax = fig.add_subplot(111,projection='3d')

    ax.scatter(geo.uphor[:,0],geo.uphor[:,1],geo.uphor[:,2],color='GREEN',s=res.bins)
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

def dr_resolution(ini,res):
    fig = plt.figure(1)
    levels=[0,0.2,0.4,0.6,0.8,1.0]
    rho = res.drho
    rhomin = res.drhomin
    rhomax = res.drhomax

    sc = plt.contour(rhomin[1][0:res.bins[0]],rhomin[2][0:res.bins[1]],
                    rhomin[0].T,levels=levels,linewidth=40)

    sc = plt.contour(rho[1][0:res.bins[0]],rho[2][0:res.bins[1]],
                    rho[0].T,levels=levels,linewidth=40)

    sc = plt.contour(rhomax[1][0:res.bins[0]],rhomax[2][0:res.bins[1]],
                    rhomax[0].T,levels=levels,linewidth=4)
    cb=plt.colorbar(sc)
    return 'SPATIAL RESOLUTION'

def dE_resolution(ini,res):
    fig = plt.figure(2)
    dE = res.dE
    dEmin = res.dEmin
    dEmax = res.dEmax
    levels = [25,40,55,70]
    sc = plt.contour(dEmin[1][0:res.bins[0]],dEmin[2][0:res.bins[1]],
                    dEmin[0].T,levels=levels,linewidth=40)

    sc = plt.contour(dE[1][0:res.bins[0]],dE[2][0:res.bins[1]],
                    dE[0].T,levels=levels,linewidth=40)

    sc = plt.contour(dEmax[1][0:res.bins[0]],dEmax[2][0:res.bins[1]],
                    dEmax[0].T,levels=levels,linewidth=4)
    cb=plt.colorbar(sc)
    return 'ENERGY RESOLUTION'

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
