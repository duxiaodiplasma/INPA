import sys
sys.path.insert(0,'/home/duxiaodi/inpa/SOURCE')
import creatobj
import numpy as np


def main():
    # lower phorsphor
    lphor = np.array([
    [-2201,  592, -764],
    [-2172,  596, -804],
    [-2186,  654, -810],
    [-2216,  651, -770]
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
    [-2225,  639, -753],
    [-2219,  640, -761],
    [-2233,  699, -767],
    [-2239,  698, -759]
                            ])

    geo = creatobj.geometry(pinhole/1e3,foil/1e3,lphor/1e3,lphor/1e3)

    return geo
