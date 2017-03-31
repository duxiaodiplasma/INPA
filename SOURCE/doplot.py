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
#    plt.clf()

    #sc = plt.scatter(res.rot_hitpoint[:,0],res.rot_hitpoint[:,2],c=res.E_birth,alpha=0.5)
    sc = plt.scatter(res.rot_hitpoint[:,0],res.rot_hitpoint[:,1],c=res.rho,alpha=0.5)
    rot_foil = np.vstack((geo.rot_foil[0,:,:],geo.rot_foil[0,0,:]))
    rot_pinhole = np.vstack((geo.rot_pinhole,geo.rot_pinhole[0,:]))
    rot_lphor = np.vstack((geo.rot_lphor,geo.rot_lphor[0,:]))
    plt.plot(rot_foil[:,0],rot_foil[:,1],alpha=0.5)
    plt.plot(rot_lphor[:,0],rot_lphor[:,1],alpha=0.5)
    plt.plot(rot_pinhole[:,0],rot_pinhole[:,1],alpha=0.5)
    cb = plt.colorbar(sc)
    plt.xlabel('x[m]')
    plt.ylabel('y[m]')
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

def dr_resolution(ini,geo,res):
    fig = plt.figure(1)
    levels=[0.05,0.15,0.35,0.65,0.95]
    rho = res.drho
    rhomin = res.drhomin
    rhomax = res.drhomax

    rot_foil = np.vstack((geo.rot_foil[0,:,:],geo.rot_foil[0,0,:]))
    rot_pinhole = np.vstack((geo.rot_pinhole,geo.rot_pinhole[0,:]))
    rot_lphor = np.vstack((geo.rot_lphor,geo.rot_lphor[0,:]))
    plt.plot(rot_foil[:,0],rot_foil[:,1],alpha=0.5)
    plt.plot(rot_lphor[:,0],rot_lphor[:,1],alpha=0.5)
    plt.plot(rot_pinhole[:,0],rot_pinhole[:,1],alpha=0.5)

    sc = plt.contour(rhomin[1][0:res.bins[0]],rhomin[2][0:res.bins[1]],
                    rhomin[0].T,levels=levels,linewidth=40)

    #sc = plt.contour(rho[1][0:res.bins[0]],rho[2][0:res.bins[1]],
    #                rho[0].T,levels=levels,linewidth=40)

    sc = plt.contour(rhomax[1][0:res.bins[0]],rhomax[2][0:res.bins[1]],
                    rhomax[0].T,levels=levels,linewidth=4)
    cb=plt.colorbar(sc)
    plt.xlabel('x[m]')
    plt.ylabel('y[m]')
    return 'SPATIAL RESOLUTION'

def dE_resolution(ini,geo,res):
    fig = plt.figure(2)
    dE = res.dE
    dEmin = res.dEmin
    dEmax = res.dEmax
    levels = [25,40,55,70]

    rot_foil = np.vstack((geo.rot_foil[0,:,:],geo.rot_foil[0,0,:]))
    rot_pinhole = np.vstack((geo.rot_pinhole,geo.rot_pinhole[0,:]))
    rot_lphor = np.vstack((geo.rot_lphor,geo.rot_lphor[0,:]))
    plt.plot(rot_foil[:,0],rot_foil[:,1],alpha=0.5)
    plt.plot(rot_lphor[:,0],rot_lphor[:,1],alpha=0.5)
    plt.plot(rot_pinhole[:,0],rot_pinhole[:,1],alpha=0.5)

    sc = plt.contour(dEmin[1][0:res.bins[0]],dEmin[2][0:res.bins[1]],
                    dEmin[0].T,levels=levels,linewidth=40)

    #sc = plt.contour(dE[1][0:res.bins[0]],dE[2][0:res.bins[1]],
    #                dE[0].T,levels=levels,linewidth=40)

    sc = plt.contour(dEmax[1][0:res.bins[0]],dEmax[2][0:res.bins[1]],
                    dEmax[0].T,levels=levels,linewidth=4)
    cb=plt.colorbar(sc)
    plt.xlabel('x[m]')
    plt.ylabel('y[m]')
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
    #    axarr[1].plot(rz_seg[i,:,0],pitch[i,:]/np.pi*180

def grid(ini,geo,res):
    fig = plt.figure(1)
    levels=[0.05,0.15,0.35,0.65,0.95]
    rho = res.drho
    dE = res.dE



    rot_foil = np.vstack((geo.rot_foil[0,:,:],geo.rot_foil[0,0,:]))
    rot_pinhole = np.vstack((geo.rot_pinhole,geo.rot_pinhole[0,:]))
    rot_lphor = np.vstack((geo.rot_lphor,geo.rot_lphor[0,:]))
    plt.plot(rot_foil[:,0],rot_foil[:,1],alpha=0.5)
    plt.plot(rot_lphor[:,0],rot_lphor[:,1],alpha=0.5)
    plt.plot(rot_pinhole[:,0],rot_pinhole[:,1],alpha=0.5)

    sc1 = plt.contour(rho[1][0:res.bins[0]],rho[2][0:res.bins[1]],
                    rho[0].T,levels=[0.1,0.3,0.5,0.7,0.9],linewidth=40)

    sc2 = plt.contour(dE[1][0:res.bins[0]],dE[2][0:res.bins[1]],
                    dE[0].T,levels=[30,45,60,75],linewidth=40)
    cb1=plt.colorbar(sc1)
    cb2=plt.colorbar(sc2)
    plt.xlabel('x[m]')
    plt.ylabel('y[m]')
    return 'SPATIAL RESOLUTION'
