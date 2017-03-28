import h5py as h5
import numpy as np

def main(grids,prob,E,E_fidasim):

    index_E = np.argmin(np.abs(E-E_fidasim))

    index = np.random.choice(np.arange(200), p=prob[:,index_E])

    f1 = grids[index,0,:]
    f2 = grids[index,1,:]
    f3 = grids[index,2,:]
    f4 = grids[index,3,:]

    return f1, f2, f3, f4

    

