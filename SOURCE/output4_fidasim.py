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

def input(fn,p1,p2,p3,p4,f1,f2,f3,f4,w,h):

# for testing the subroutine
#aaa=0.5
#if aaa < 1 :
#
#    fn = './output/inpa.geometry'
#    p1 = np.array([-2.173863, 0.58766122, -0.71087029])
#    p2 = np.array([-2.17645119,  0.57800196,-0.71087029])
#    p3 = np.array([-2.17280322,  0.57702449, -0.72012971])
#    p4 = np.array([-2.17021503,0.58668375,-0.72012971])
#
#    f2,f3,f4,f1 = [-2.29053266,0.41084848, -0.81155492],\
#                  [-2.28636801,0.40973257, -0.82212581],\
#                  [-2.25311062, 0.60285913, -0.77833584],\
#                  [-2.25705344, 0.60391561,-0.76832801]
#
#    w = 20
#    h = 10 
#

    a_cent = np.zeros((w*h,3))
    a_tedge = np.zeros((w*h,3))
    a_redge = np.zeros((w*h,3))
    for k in np.arange(w*h):
       a_cent[k,:]  = (p1+p2+p3+p4)/4.
       a_tedge[k,:] = (p1+p2)/2.
       a_redge[k,:] = (p1+p4)/2.
    
    
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
#    foil = np.zeros((4,3))
#    foil[0,:] = f1
#    foil[1,:] = f2
#    foil[2,:] = f3
#    foil[3,:] = f4
#    import matplotlib.pylab as plt
#    from mpl_toolkits.mplot3d import Axes3D
#    fig = plt.figure()
#    fig.clf()
#    ax = fig.add_subplot(111, projection='3d')
#   
#    ax.scatter(foil[:,0],foil[:,1],foil[:,2],color='green')
#    ax.scatter(d_cent[:,0],d_cent[:,1],d_cent[:,2])
#    ax.scatter(d_redge[:,0],d_redge[:,1],d_redge[:,2],color='red')
#    ax.scatter(d_tedge[:,0],d_tedge[:,1],d_tedge[:,2],color='black')
#
#    print('sucessfully generate geometry of INPA for FIDASIM')

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
    
    grids = np.reshape(grids,(w*h,4,3))
    np.savez('./output/grids',grids = grids) 
 
    return grids

#
