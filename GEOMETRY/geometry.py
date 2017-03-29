import sys
sys.path.insert(0,'/home/duxiaodi/inpa/SOURCE')
import creatobj
import numpy as np


def main():
    # lower phorsphor
    lphor = np.array([
    [-2223,  593, -780],
    [-2193,  596, -820],
    [-2181,  658, -806],
    [-2210,  655, -765]
                    ])

    ## pinhole
    #pinhole = np.array([
    #[-2177,  708.9, -714.9],
    #[-2174,  709.2, -719],
    #[-2171.5,705.1, -717.4],
    #[-2174.5,704.8, -713.4]
    #                ])

    # 0.2mm*0.2mm pinhole
    pinhole = np.array([
    [-2174.4,  706.1, -715.1],
    [-2173.2,  706.2, -716.7],
    [-2174.2, 707.8, -717.3],
    [-2175.4, 707.7, -715.7]
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
