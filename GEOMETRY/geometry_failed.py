import sys
sys.path.insert(0,'/home/duxiaodi/inpa/SOURCE')
import creatobj
import numpy as np


def main():
    # lower phorsphor
    lphor = np.array([
    [-2241.9,    600.1,   -793.2],
    [-2212.4,    603.3,   -833.35],
    [-2202.6,    701.7,   -818.3],
    [-2232.14,   698.5,   -778.1]
                    ])

    # pinhole
    pinhole = np.array([
    [-2177,  708.9, -714.9],
    [-2174,  709.2, -719],
    [-2171.5,705.1, -717.4],
    [-2174.5,704.8, -713.4]
                    ])

    # stripping foil
    foil = np.zeros((1,4,3))

    foil[0,:,:] = np.array([
    [-2253.8,    598.9,   -777.08],
    [-2247.9,    599.5,   -785.1],
    [-2238.1,    697.9,   -770.1],
    [-2244,      697.2,   -762.0]
                            ])

    geo = creatobj.geometry(pinhole/1e3,foil/1e3,lphor/1e3,lphor/1e3)

    return geo
