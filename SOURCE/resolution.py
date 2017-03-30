import numpy as np
import fbw
import lib
import bfield


def main(ini,res):
    # segments of each birth sightline
    segs = 21

    # segments between crossing points of NB
    rz_seg = np.zeros((ini.mc,segs,3))

    #??
    vsl = np.zeros((ini.mc,segs,3))

    for i in np.arange(ini.mc):
        n0_p = res.birthsl[i,0,:]
        n0_f = res.birthsl[i,1,:]
        rz_seg[i,:,:] = fbw.main(n0_p,n0_f,segs)

    # calculate local psi
    psi_seg = np.zeros((ini.mc,segs))

    # rho resolution
    rho = np.zeros((ini.mc,segs))
    rhoavg = np.zeros(ini.mc)
    drho = np.zeros(ini.mc)

    # pitch resolution
    dpitch = np.zeros(ini.mc)

    # ??
    R_mid_avg = np.zeros(ini.mc)
    # b vector
    vb = np.zeros((ini.mc,segs,3))
    # ??
    pitch = np.zeros((ini.mc,segs))
    pitchavg = np.zeros(ini.mc)

    for i in np.arange(ini.mc):
        for j in np.arange(segs):
            psi_seg[i,j] = ini.fpsi(rz_seg[i,j,1],rz_seg[i,j,0])
            rho[i,j] = np.sqrt(psi_seg[i,j])
            vb[i,j,:] = bfield.brzt(rz_seg[i,j,0],rz_seg[i,j,1],rz_seg[i,j,2],ini)
            pitch[i,j] = lib.angle_between(vb[i,j,:], res.birthsl[i,1,:]-res.birthsl[i,0,:])

        R_mid_avg[i] = np.mean(rz_seg[i,:,0])
        rhoavg[i] = np.mean(rho[i,:])
        pitchavg[i] = np.mean(pitch[i,:])
        drho[i] = np.amax(np.sqrt(psi_seg[i,:]))-np.amin(np.sqrt(psi_seg[i,:]))
        dpitch[i] = np.amax(pitch[i,:]) - np.amin(pitch[i,:])

    res.rho = rhoavg
    res.rho_eb = drho

    res.pitch = pitchavg
    res.pitch_eb = drho

    res.rzseg = rz_seg
    res.segs = segs
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


