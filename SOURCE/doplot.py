import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.gridspec as gridspec

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

    sc = plt.scatter(res.rot_hitpoint[:,0],res.rot_hitpoint[:,1],c=res.E_birth,alpha=0.5)
    #sc = plt.scatter(res.rot_hitpoint[:,0],res.rot_hitpoint[:,1],c=res.rho,alpha=0.5)
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
    levels=[0.1,0.3,0.5,0.7,0.9]
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
    levels = [30,45,60,75]

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

def grid(ini,geo,res):
    fig = plt.figure(1)
    levels=[0.05,0.15,0.35,0.65,0.95]
    rho = res.drho
    dE = res.dE



    for i in range(0,7):
        rot_foil = np.vstack((geo.rot_foil[i,:,:],geo.rot_foil[i,0,:]))
        plt.plot(rot_foil[:,0],rot_foil[:,1],alpha=0.5)

    rot_pinhole = np.vstack((geo.rot_pinhole,geo.rot_pinhole[0,:]))
    rot_lphor = np.vstack((geo.rot_lphor,geo.rot_lphor[0,:]))
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

def all(ini,geo,res):
    plt.rcParams['figure.figsize'] = (8,8)
    fig = plt.figure(figsize=(8,8))
    gs = gridspec.GridSpec(100,100)
    ax0 = plt.subplot(gs[5:35,10:65])
    ax1 = plt.subplot(gs[35:65,10:65])
    ax2 = plt.subplot(gs[65:95,10:65])

    for i in range(0,7):
        rot_foil = np.vstack((geo.rot_foil[i,:,:],geo.rot_foil[i,0,:]))
        ax0.plot(rot_foil[:,0],rot_foil[:,1],alpha=0.5)

    rot_pinhole = np.vstack((geo.rot_pinhole,geo.rot_pinhole[0,:]))
    rot_lphor = np.vstack((geo.rot_lphor,geo.rot_lphor[0,:]))
    ax0.plot(rot_lphor[:,0],rot_lphor[:,1],alpha=0.5)
    ax0.plot(rot_pinhole[:,0],rot_pinhole[:,1],alpha=0.5)

    for i in range(0,7):
        rot_foil = np.vstack((geo.rot_foil[i,:,:],geo.rot_foil[i,0,:]))
        ax1.plot(rot_foil[:,0],rot_foil[:,1],alpha=0.5)
    ax1.plot(rot_lphor[:,0],rot_lphor[:,1],alpha=0.5)
    ax1.plot(rot_pinhole[:,0],rot_pinhole[:,1],alpha=0.5)

    for i in range(0,7):
        rot_foil = np.vstack((geo.rot_foil[i,:,:],geo.rot_foil[i,0,:]))
        ax2.plot(rot_foil[:,0],rot_foil[:,1],alpha=0.5)
    ax2.plot(rot_lphor[:,0],rot_lphor[:,1],alpha=0.5)
    ax2.plot(rot_pinhole[:,0],rot_pinhole[:,1],alpha=0.5)

    rho = res.drho
    rhomin = res.drhomin
    rhomax = res.drhomax
    drhostd = res.rhostd
    # 80% confidnece
    cl = 1.28

    for i in range(0,7):
        rot_foil = np.vstack((geo.rot_foil[i,:,:],geo.rot_foil[i,0,:]))
        ax0.plot(rot_foil[:,0],rot_foil[:,1],alpha=0.5)
    rot_pinhole = np.vstack((geo.rot_pinhole,geo.rot_pinhole[0,:]))
    rot_lphor = np.vstack((geo.rot_lphor,geo.rot_lphor[0,:]))


    levels_rho = [0.1,0.3,0.5,0.7,0.9]
    sc = ax0.contour(rhomin[1][0:res.bins[0]],rhomin[2][0:res.bins[1]],
                    rho[0].T-cl*drhostd[0].T,levels=levels_rho,linewidth=40)
    #plt.colorbar(sc,use_gridspec=True)
    ax0.contour(rhomax[1][0:res.bins[0]],rhomax[2][0:res.bins[1]],
                    rho[0].T+cl*drhostd[0].T,levels=levels_rho,linewidth=4)

    dE = res.dE
    dEmin = res.dEmin
    dEmax = res.dEmax
    dEstd = res.Estd
    levels_E = [25,40,55,70]
    sc1=ax1.contour(dEmin[1][0:res.bins[0]],dEmin[2][0:res.bins[1]],
                    dE[0].T-cl*dEstd[0].T,levels=levels_E,linewidth=40)

    ax1.contour(dEmax[1][0:res.bins[0]],dEmax[2][0:res.bins[1]],
                    dE[0].T+cl*dEstd[0].T,levels=levels_E,linewidth=40)

    sc21=ax2.contour(rho[1][0:res.bins[0]],rho[2][0:res.bins[1]],
                    rho[0].T,levels=levels_rho,linewidth=40)

    sc22=ax2.contour(dE[1][0:res.bins[0]],dE[2][0:res.bins[1]],
                    dE[0].T,levels=levels_E,linewidth=40)
    #plt.clabel(sc22, fmt = '%2.1d', fontsize=10) #contour line labels
    #plt.clabel(sc21, fmt = '%3.1f', fontsize=10) #contour line labels


    ax2.set_xlabel('x[m]')
    ax0.set_ylabel('y[m]')
    ax1.set_ylabel('y[m]')
    ax2.set_ylabel('y[m]')

    # and transform them after to get the ABSOLUTE POSITION AND DIMENSIONS
    x0, y0, width, height = [1.2, 0.5, 0.5, 0.1]
    fig.tight_layout()
    import matplotlib
    Bbox = matplotlib.transforms.Bbox.from_bounds(x0, y0, width, height)
    trans = ax0.transAxes + fig.transFigure.inverted()
    l, b, w, h = matplotlib.transforms.TransformedBbox(Bbox, trans).bounds
    # Now just create the axes and the colorbar
    cbaxes = fig.add_axes([l, b, w, h])
    cbar = plt.colorbar(sc, cax=cbaxes, orientation='horizontal')
    cbar.ax.tick_params(labelsize=9)

    trans = ax1.transAxes + fig.transFigure.inverted()
    l, b, w, h = matplotlib.transforms.TransformedBbox(Bbox, trans).bounds
    # Now just create the axes and the colorbar
    cbaxes = fig.add_axes([l, b, w, h])
    cbar = plt.colorbar(sc1, cax=cbaxes, orientation='horizontal')
    cbar.ax.tick_params(labelsize=9)

    import time
    comment0 = 'INPASIM '
    comment1 =  time.strftime("%c")+' BY X.D. Du'
    comment2 = 'Filename: ' + np.str(ini.fn)
    comment3 = 'Comment: ' + np.str(ini.comment)
    comment4 =  np.str(ini.gfile)
    comment5 = 'confidence level: ' + np.str(cl*100)+'% of standard deviation'
    #comment5 = '$\\tilde{S}$ = ' + np.str(np.abs(a.S[0]))
    #comment6 = 'wavelength from ' + np.str(wl[0]/10.)+' to '+np.str(wl[1]/10.) + ' [nm]'
    plt.text(1.02, 0.8, comment0, ha='left', va='center', transform=ax2.transAxes,fontsize=9.0)
    plt.text(1.02, 0.7, comment1, ha='left', va='center', transform=ax2.transAxes,fontsize=9.0)
    plt.text(1.02, 0.5, comment2, ha='left', va='center', transform=ax2.transAxes,fontsize=9.0)
    plt.text(1.02, 0.4,comment3, ha='left', va='center', transform=ax2.transAxes,fontsize=9.0)
    plt.text(1.02, 0.25, comment4, ha='left', va='center', transform=ax2.transAxes,fontsize=8.0)
    plt.text(1.02, 0.15,comment5, ha='left', va='center', transform=ax2.transAxes,fontsize=8.0)
    #plt.text(1.05, 0.1, comment6, ha='left', va='center', transform=ax2.transAxes,fontsize=10.0)
    fig.savefig(ini.fn+'.ps',dpi=fig.dpi)
    import os
    return 'file save in: '+os.getcwd()+'/'+ini.fn+'.ps'


