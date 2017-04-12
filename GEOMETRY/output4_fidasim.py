"""
;+
; NAME:
;     output4_fidasim
; PURPOSE:
      calculate the geometry of INPA for input of the FIDASIM
; EXPLANATION:
;
; CALLING SEQUENCE:
#;
; OPTIONAL INPUT:
;
; OPTIONAL KEYWORD INPUT:
;
; OPTIONAL KEYWORD OUTPUTS:

; PROCEDURE:

; EXAMPLE:
; RESTRICTIONS:
;
; PROCEDURES USED:
;
; REVISION HISTORY:
"""

import numpy as np
import lib

def find_corners(origin,v1,v2):
    """
       v1
    p1 ->- p2
    |      |v2
    p4 -<- p3
    """
    p1 = origin
    p2 = origin + v1
    p3 = p2 + v2
    p4 = p3 - v1
    return p1,p2,p3,p4

def main3(fn,geo,w=2,h=10):
    """
    OUTPUT FOR FIDASIM IN CDR PHASE
    """

    # number of foil
    nfoil = np.shape(geo.foil)[0]

    # total channels
    nchan = nfoil*w*h
    a_cent = np.zeros((nchan,3))
    a_tedge = np.zeros((nchan,3))
    a_redge = np.zeros((nchan,3))

    # pinhole for each channel
    for i in range(0,nchan):
        p1 = geo.pinhole[0,:]
        p2 = geo.pinhole[1,:]
        p3 = geo.pinhole[2,:]
        p4 = geo.pinhole[3,:]

        a_cent[i,:] = (p1+p2+p3+p4)/4.
        a_tedge[i,:] = (p1+p2)/2.
        a_redge[i,:] = (p2+p3)/2.

    grids = np.zeros((nfoil,h,w,4,3))
    origin = np.zeros(3)
    # produce grids on each foil
    for i in range(0,nfoil):
        f1 = geo.foil[i,0,:]
        f2 = geo.foil[i,1,:]
        f3 = geo.foil[i,2,:]
        f4 = geo.foil[i,3,:]
        v1 = lib.unit_vector(f2-f1)*(np.linalg.norm(f2-f1)/w)
        v2 = lib.unit_vector(f4-f1)*(np.linalg.norm(f4-f1)/h)
        for j in range(0,h):
            # starting point
            origin[0] = np.linspace(f1[0],f4[0],h+1)[j]
            origin[1] = np.linspace(f1[1],f4[1],h+1)[j]
            origin[2] = np.linspace(f1[2],f4[2],h+1)[j]

            for k in range(0,w):
                grids[i,j,k,0,:], \
                grids[i,j,k,1,:], \
                grids[i,j,k,2,:], \
                grids[i,j,k,3,:], \
                = find_corners(origin,v1,v2)
                origin = origin + v1
                print(grids[i,j,k,:,:])

    # check the grids
    import matplotlib.pylab as plt
    from mpl_toolkits.mplot3d import Axes3D
    fig = plt.figure()
    fig.clf()
    ax = fig.add_subplot(111, projection='3d')
    for i in range(0,nfoil):
        for j in range(0,h):
            for k in range(0,w):
                tmp = np.vstack((grids[i,j,k,:,:],grids[i,j,k,0,:]))
                ax.plot(tmp[:,0],tmp[:,1],tmp[:,2])
                ax.set_xlabel('x[m]')
                ax.set_ylabel('y[m]')
                ax.set_zlabel('z[m]')

    # translate to FIDASIM language
    # dummy foil df
    # looking through pinhole at the detector
    # d_tedge: detector top edge
    # d_redge: detector right edge
    d_cent = np.zeros((nfoil,h,w,3))
    d_tedge = np.zeros((nfoil,h,w,3))
    d_redge = np.zeros((nfoil,h,w,3))
    for i in range(0,nfoil):
        for j in range(0,h):
            for k in range(0,w):
                df1 = grids[i,j,k,0,:]
                df2 = grids[i,j,k,1,:]
                df3 = grids[i,j,k,2,:]
                df4 = grids[i,j,k,3,:]
                d_cent[i,j,k,:] = (df1+df2+df3+df4)/4.
                d_tedge[i,j,k,:] = (df1+df2)/2.
                d_redge[i,j,k,:] = (df2+df3)/2.

    d_cent = np.reshape(d_cent,(nfoil*w*h,3))
    d_tedge = np.reshape(d_tedge,(nfoil*w*h,3))
    d_redge = np.reshape(d_redge,(nfoil*w*h,3))

    np.savetxt(fn,np.concatenate((a_cent,a_tedge,a_redge,d_cent,d_tedge,d_redge)))
    print('INPA GEOMETRY IS SAVED IN: '+ fn)


def main2(fn,geo):

    a_cent = np.zeros((7,3))
    a_tedge = np.zeros((7,3))
    a_redge = np.zeros((7,3))

    for i in range(0,7):
        p1 = geo.pinhole[0,:]
        p2 = geo.pinhole[1,:]
        p3 = geo.pinhole[2,:]
        p4 = geo.pinhole[3,:]

        a_cent[i,:] = (p1+p2+p3+p4)/4.
        a_tedge[i,:] = (p1+p2)/2.
        a_redge[i,:] = (p2+p3)/2.

    d_cent = np.zeros((7,3))
    d_tedge = np.zeros((7,3))
    d_redge = np.zeros((7,3))

    for i in range(0,7):
        f1 = geo.foil[i,0,:]
        f2 = geo.foil[i,1,:]
        f3 = geo.foil[i,2,:]
        f4 = geo.foil[i,3,:]

        d_cent[i,:] = (f1+f2+f3+f4)/4.
        d_tedge[i,:] = (f1+f2)/2.
        d_redge[i,:] = (f2+f3)/2.

    np.savetxt(fn,np.concatenate((a_cent,a_tedge,a_redge,d_cent,d_tedge,d_redge)))
    print('INPA GEOMETRY IS SAVED IN: '+ fn)
    return



def main(fn,geo,w,h):

    p1 = geo.pinhole[0,:]
    p2 = geo.pinhole[1,:]
    p3 = geo.pinhole[2,:]
    p4 = geo.pinhole[3,:]

    f1 = geo.foil[0,0,:]
    f2 = geo.foil[0,1,:]
    f3 = geo.foil[0,2,:]
    f4 = geo.foil[0,3,:]

    a_cent = np.zeros((w*h,3))
    a_tedge = np.zeros((w*h,3))
    a_redge = np.zeros((w*h,3))
    for k in np.arange(w*h):
       a_cent[k,:]  = (p1+p2+p3+p4)/4.
       a_tedge[k,:] = (p1+p2)/2.
       # seems wrong...
       a_redge[k,:] = (p1+p4)/2.
       #a_redge[k,:] = (p2+p3)/2.


    redge = np.zeros((h,3))
    redge[:,0] = np.linspace(f2[0],f3[0],2*h+1)[1::2]
    redge[:,1] = np.linspace(f2[1],f3[1],2*h+1)[1::2]
    redge[:,2] = np.linspace(f2[2],f3[2],2*h+1)[1::2]

    ledge = np.zeros((h,3))
    ledge[:,0] = np.linspace(f1[0],f4[0],2*h+1)[1::2]
    ledge[:,1] = np.linspace(f1[1],f4[1],2*h+1)[1::2]
    ledge[:,2] = np.linspace(f1[2],f4[2],2*h+1)[1::2]

    d_cent = np.zeros((w,h,3))
    for i in np.arange(h):
        d_cent[:,i,0] = np.linspace(ledge[i,0],redge[i,0],2*w+1)[1::2]
        d_cent[:,i,1] = np.linspace(ledge[i,1],redge[i,1],2*w+1)[1::2]
        d_cent[:,i,2] = np.linspace(ledge[i,2],redge[i,2],2*w+1)[1::2]
    d_cent = np.reshape(d_cent,(w*h,3))

    d_redge = np.zeros((w,h,3))
    for i in np.arange(h):
        d_redge[:,i,0] = np.linspace(ledge[i,0],redge[i,0],2*w+1)[0::2][0:w]
        d_redge[:,i,1] = np.linspace(ledge[i,1],redge[i,1],2*w+1)[0::2][0:w]
        d_redge[:,i,2] = np.linspace(ledge[i,2],redge[i,2],2*w+1)[0::2][0:w]
        #d_redge[:,i,0] = np.linspace(ledge[i,0],redge[i,0],2*w+1)[2::2][0:w]
        #d_redge[:,i,1] = np.linspace(ledge[i,1],redge[i,1],2*w+1)[2::2][0:w]
        #d_redge[:,i,2] = np.linspace(ledge[i,2],redge[i,2],2*w+1)[2::2][0:w]
    d_redge = np.reshape(d_redge,(w*h,3))

    # --
    redge = np.zeros((h,3))
    redge[:,0] = np.linspace(f2[0],f3[0],h+1)[0:h]
    redge[:,1] = np.linspace(f2[1],f3[1],h+1)[0:h]
    redge[:,2] = np.linspace(f2[2],f3[2],h+1)[0:h]

    ledge = np.zeros((h,3))
    ledge[:,0] = np.linspace(f1[0],f4[0],h+1)[0:h]
    ledge[:,1] = np.linspace(f1[1],f4[1],h+1)[0:h]
    ledge[:,2] = np.linspace(f1[2],f4[2],h+1)[0:h]

    d_tedge = np.zeros((w,h,3))
    for i in np.arange(h):
        d_tedge[:,i,0] = np.linspace(ledge[i,0],redge[i,0],2*w+1)[1::2]
        d_tedge[:,i,1] = np.linspace(ledge[i,1],redge[i,1],2*w+1)[1::2]
        d_tedge[:,i,2] = np.linspace(ledge[i,2],redge[i,2],2*w+1)[1::2]
    d_tedge = np.reshape(d_tedge,(w*h,3))

    np.savetxt(fn,
               np.concatenate((a_cent,a_tedge,a_redge,d_cent,d_tedge,d_redge)))

##    CHECK THE GEOMETRY
    foil = np.zeros((4,3))
    foil[0,:] = f1
    foil[1,:] = f2
    foil[2,:] = f3
    foil[3,:] = f4
    import matplotlib.pylab as plt
    from mpl_toolkits.mplot3d import Axes3D
    fig = plt.figure()
    fig.clf()
    ax = fig.add_subplot(111, projection='3d')

    ax.scatter(foil[:,0],foil[:,1],foil[:,2],color='green')
    ax.scatter(d_cent[:,0],d_cent[:,1],d_cent[:,2])
    ax.scatter(d_redge[:,0],d_redge[:,1],d_redge[:,2],color='red')
    ax.scatter(d_tedge[:,0],d_tedge[:,1],d_tedge[:,2],color='black')
#
    print('sucessfully generate geometry of INPA for FIDASIM')

#
    grids = np.zeros((w,h,4,3))
    g_ledge = np.zeros((h+1,3))
    g_redge = np.zeros((h+1,3))
    g_tedge = np.zeros((w+1,h+1,3))

    g_ledge[:,0] = np.linspace(f1[0],f4[0],h+1)
    g_ledge[:,1] = np.linspace(f1[1],f4[1],h+1)
    g_ledge[:,2] = np.linspace(f1[2],f4[2],h+1)

    g_redge[:,0] = np.linspace(f2[0],f3[0],h+1)
    g_redge[:,1] = np.linspace(f2[1],f3[1],h+1)
    g_redge[:,2] = np.linspace(f2[2],f3[2],h+1)

    for i in np.arange(h+1):
        g_tedge[:,i,0] = np.linspace(g_ledge[i,0],g_redge[i,0],w+1)
        g_tedge[:,i,1] = np.linspace(g_ledge[i,1],g_redge[i,1],w+1)
        g_tedge[:,i,2] = np.linspace(g_ledge[i,2],g_redge[i,2],w+1)

    for i in np.arange(w):
        for j in np.arange(h):
            grids[i,j,0,:] = g_tedge[i,j,:]
            grids[i,j,1,:] = g_tedge[i,j+1,:]
            grids[i,j,2,:] = g_tedge[i+1,j+1,:]
            grids[i,j,3,:] = g_tedge[i+1,j,:]


    return grids

#
