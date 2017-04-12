import numpy as np
import scipy
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from scipy.interpolate import griddata


def main():

   '''visulazation of the results'''
   geometry = np.load('./output/geometry.npz')
   foil = geometry['foil']
   pinhole = geometry['pinhole']

   output = np.load('./output/inpa_grid.npz')
   sp = output['strike_pos']
   e_ini = output['E_ini']
   r_ini = output['R_ini']

   num = 100000
   point = np.zeros((num,2))
   x = sp[:,1]/np.cos(15/180*np.pi)
   z = sp[:,2]/np.cos(15/180*np.pi)

   xfit = np.linspace(0.35,0.6,100)

   R = np.arange(2.2,1.4,-0.1)
   E = np.arange(20,80,10)
   fit = np.zeros((len(R),2))
   grid_R = np.zeros((len(R),len(xfit)))
   grid_E = np.zeros((len(E),len(xfit)))

   for i in np.arange(len(R)):
       mask = np.where(((r_ini[:,0]>R[i]-0.002) & (r_ini[:,0]<R[i]+0.002)))
       fit[i,:] = np.polyfit(x[mask],z[mask],1)
       fit_fn = np.poly1d(fit[i,:])
       grid_R[i,:] = fit_fn(xfit)

   fit = np.zeros((len(E),6))
   for i in np.arange(len(E)):
       mask = np.where(((e_ini>E[i]-1) & (e_ini<E[i]+1)))
       fit[i,:] = np.polyfit(x[mask],z[mask],5)
       fit_fn = np.poly1d(fit[i,:])
       grid_E[i,:] = fit_fn(xfit)

#   fig = plt.figure(1)
#   fig.clf()
#   for i in np.arange(len(R)):
#       plt.plot(xfit,grid_R[i,:],color='black')
#   for i in np.arange(len(E)):
#       plt.plot(xfit,grid_E[i,:],color='black')
#

   return xfit,grid_R,  grid_E, R,E














