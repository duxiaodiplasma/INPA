import numpy as np
import fbw
import lib
import bfield
import scipy.stats

def main(ini,geo,res):
    geo.nv = geo.uphoreq[0:3]
    geo.targetnv = [0,1,0]
    axis = np.cross(geo.nv, geo.targetnv)
    angle = lib.angle_between(geo.nv,geo.targetnv)

    # rotation matrix
    matrix = lib.rotation_matrix(axis,angle)

    #
    res.rot_hitpoint = np.zeros((ini.mc,3))
    for i in range(0,ini.mc):
        res.rot_hitpoint[i,:] = np.dot(matrix,res.hitpoint[i,:])

    # rotating phorsphor, pinhole
    geo.rot_lphor = np.zeros((4,3))
    geo.rot_uphor = np.zeros((4,3))
    geo.rot_pinhole = np.zeros((4,3))
    for i in range(0,4):
        geo.rot_lphor[i,:] = np.dot(matrix,geo.lphor[i,:])
        geo.rot_uphor[i,:] = np.dot(matrix,geo.uphor[i,:])
        geo.rot_pinhole[i,:] = np.dot(matrix,geo.pinhole[i,:])

    geo.rot_foil = np.zeros((1,4,3))
    for i in range(0,4):
        geo.rot_foil[0,i,:] = np.dot(matrix,geo.foil[0,i,:])

    rhomax = scipy.stats.binned_statistic_2d(res.rot_hitpoint[:,0],res.rot_hitpoint[:,2],res.rho,bins=100,statistic=np.max)
    rhomin = scipy.stats.binned_statistic_2d(res.rot_hitpoint[:,0],res.rot_hitpoint[:,2],res.rho,bins=100,statistic=np.min)
    rho = scipy.stats.binned_statistic_2d(res.rot_hitpoint[:,0],res.rot_hitpoint[:,2],res.rho,bins=100,statistic=np.mean)














    return res



# -- split line --
# plot results
#fig = plt.figure()
#levels = np.linspace(0,1,21)
#cs = plt.contour(rr,zz,np.sqrt(np.transpose(psirz)),levels)
#cb = plt.colorbar(cs, shrink=0.8, extend='both',ticks=levels[0::2])
#
#
##plt.clabel(cs, levels[0::2], inline=1, fontsize=10)
#for i in np.arange(100):
#    plt.plot(rz_seg[i,:,0], rz_seg[i,:,1])
#

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
#
#axarr[1].set_xlabel('R[m]')
#axarr[1].set_ylabel('pitch angle [degree]')
##plt.xlim(1.1,2.5)
##plt.ylim(-1.6,1.6)


