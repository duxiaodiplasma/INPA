import sys
sys.path.insert(0,'/home/duxiaodi/inpa/SOURCE')
import creatobj
import numpy as np


def main():
    # lower phorsphor
    lphor = np.array([
    [-2216,  592, -770],
    [-2187,  592, -811],
    [-2180,  662, -806],
    [-2209,  662, -766]
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
    [-2206,  702, -746],
    [-2212,  702, -738]
                            ])

    geo = creatobj.geometry(pinhole/1e3,foil/1e3,lphor/1e3,lphor/1e3)

    return geo
