import h5py as h5
import numpy as np

def main(fn):
    f = h5.File(fn)
    flux_all =  f['flux'].value
    energy = f['energy'].value
    pitch = f['pitch'].value
    
    prob = np.zeros((200,50))
    for i in np.arange(50):
        prob[:,i] = flux_all[:,i]/(np.sum(flux_all[:,i])+1)
    
    return prob, energy

