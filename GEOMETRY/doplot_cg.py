import sys
import numpy as np

"""local lib"""
import matplotlib.pylab as plt
from mpl_toolkits.mplot3d import Axes3D
import h5py as h5

import geometry

fn_grid = '/home/duxiaodi/fidasim/output/159243/20170412/159243.inpa801_neutrals.h5'
fgrid = h5.File(fn_grid)

nx = fgrid['grid']['nx'].value
ny = fgrid['grid']['ny'].value
nz = fgrid['grid']['nz'].value

# cm -> m
x = fgrid['grid']['x_grid'].value/100.
y = fgrid['grid']['y_grid'].value/100.
z = fgrid['grid']['z_grid'].value/100.

plt.clf()
fig = plt.figure()
ax = fig.add_subplot(111,projection='3d')
ax.scatter(x[0::2],y[0::2],z[0::2],color='black',alpha=0.04)

geo =  geometry.main()

ax.scatter(geo.pinhole[:,0],geo.pinhole[:,1],geo.pinhole[:,2])

d_cent = np.zeros((7,3))
a_cent = np.zeros((7,3))
for i in range(0,7):
    #foil[i,:,:] = np.vstack((geo.foil[i,:,:],geo.foil[i,0,:]))
    d_cent[i,:] = (geo.foil[i,0,:]+geo.foil[i,1,:]+geo.foil[i,2,:]+geo.foil[i,3,:])/4.
    a_cent[i,:] = (geo.pinhole[0,:]+geo.pinhole[1,:]+geo.pinhole[2,:]+geo.pinhole[3,:])/4.

dummy = np.zeros((7,3))
for i in range(0,7):
    dummy[i,:] = (a_cent[i,:] - d_cent[i,:])*30+d_cent[i,:]
    ax.plot([d_cent[i,0],a_cent[i,0],dummy[i,0]],
            [d_cent[i,1],a_cent[i,1],dummy[i,1]],
            [d_cent[i,2],a_cent[i,2],dummy[i,2]])
ax.set_xlabel('x')
ax.set_ylabel('y')
ax.set_zlabel('z')


