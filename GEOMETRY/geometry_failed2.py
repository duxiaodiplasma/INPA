import sys
sys.path.insert(0,'/home/duxiaodi/inpa/SOURCE')
import creatobj
import numpy as np


def main():
    # lower phorsphor
    lphor = np.array([
    [-2212,   621.5, -767.4],
    [-2182.4, 624.7, -807.6],
    [-2232.9, 705.4, -838.2],
    [-2262.4, 702.2, -798.1]
                    ])

    # upper phorsphor
    uphor = np.array([
    [-2253.8, 606.4, -738.3],
    [-2225,   621.7, -745.2],
    [-2275.5, 702.4, -775.9],
    [-2304.3, 687.1, -769.0]
                    ])

    # pinhole
    pinhole = np.array([
    [-2177,  708.9, -714.9],
    [-2174,  709.2, -719],
    [-2171.5,705.1, -717.4],
    [-2174.5,704.8, -713.4]
                    ])

    # stripping foil
    foil = np.zeros((10,4,3))

    foil[0,:,:] = np.array([
    [-2272.7160078,699.1895561,-783.4911767],
    [-2269.4082698,693.8987660,-781.4791583],
    [-2263.9345410,694.4913348,-788.9196986],
    [-2267.2422790,699.7821249,-790.9317169]
                            ])

    foil[1,:,:] = np.array([
    [-2267.6668528,691.1133363,-780.4198958],
    [-2264.3591148,685.8225462,-778.4078775],
    [-2258.3860600,686.8810933,-785.4439714],
    [-2262.1931240,691.7059051,-787.8604361]
                            ])

    foil[2,:,:] = np.array([
    [-2262.6176979,683.0371165,-777.3486150],
    [-2259.3099599,677.7463264,-775.3365966],
    [-2253.8362311,678.3388952,-782.7771369],
    [-2257.1439691,683.6296853,-784.7891552]
                            ])

    foil[3,:,:] = np.array([
    [-2257.5685429,674.9608967,-774.2773341],
    [-2254.2608050,669.6701067,-772.2653158],
    [-2248.2877501,670.7286537,-779.3014097],
    [-2252.0948141,675.5534655,-781.7178744]
                            ])

    foil[4,:,:] = np.array([
    [-2252.5193880,666.8846769,-771.2060533],
    [-2249.2116500,661.5938869,-769.1940349],
    [-2243.7379212,662.1864557,-776.6345752],
    [-2247.0456592,667.4772457,-778.6465935]
                            ])

    foil[5,:,:] = np.array([
    [-2247.4702330,658.8084571,-768.1347724],
    [-2244.1624951,653.5176671,-766.1227541],
    [-2238.1894403,654.5762141,-773.1588480],
    [-2241.4971782,659.8670042,-775.1708663]
                            ])

    foil[6,:,:] = np.array([
    [-2242.4210781,650.7322374,-765.0634916],
    [-2239.1133401,645.4414473,-763.0514732],
    [-2233.1402853,646.4999944,-770.0875672],
    [-2236.9473493,651.3248062,-772.5040318]
                            ])

    foil[7,:,:] = np.array([
    [-2237.3719231,642.6560176,-761.9922107],
    [-2234.0641852,637.3652275,-759.9801924],
    [-2228.5904564,637.9577963,-767.4207327],
    [-2231.8981943,643.2485864,-769.4327510]
                            ])

    foil[8,:,:] = np.array([
    [-2232.3227682,634.5797978,-758.9209299],
    [-2229.0150302,629.2890077,-756.9089115],
    [-2223.0419754,630.3475548,-763.9450055],
    [-2226.8490394,635.1723666,-766.3614701]
                            ])

    foil[9,:,:] = np.array([
    [-2227.2736132,626.5035780,-755.8496490],
    [-2223.9658753,621.2127879,-753.8376307],
    [-2217.9928205,622.2713350,-760.8737246],
    [-2221.7998845,627.0961468,-763.2901893]
                            ])

    geo = creatobj.geometry(pinhole/1e3,foil/1e3,uphor/1e3,lphor/1e3)

    return geo
