import numpy as np
import lib
class geometry(object):
    """
    Attributes:
        phorphor: phorsphor 4 cornors' location
                 dimension (4,3)
        foil : foil cornors location
                  (10,4,3)
        pinhole: pinhole coordinate
                  dimension (4,3)
    """

    def __init__(self,pinhole,foil,uphor,lphor):
        """Return INPA GEMOETRY"""
        self.uphor = uphor
        self.lphor = lphor
        self.foil = foil
        self.pinhole = pinhole

        A = self.uphor[0,:]
        B = self.uphor[1,:]
        C = self.uphor[2,:]
        v1 = lib.unit_vector(B-A)
        v2 = lib.unit_vector(C-A)
        # normal vector
        nv = np.cross(v1,v2)
        k = -np.dot(nv,A)
        self.uphoreq = np.concatenate((nv,[k]))

        A = self.lphor[0,:]
        B = self.lphor[1,:]
        C = self.lphor[2,:]
        v1 = lib.unit_vector(B-A)
        v2 = lib.unit_vector(C-A)
        # normal vector
        nv = np.cross(v1,v2)
        k = -np.dot(nv,A)
        self.lphoreq = np.concatenate((nv,[k]))


class ini(object):
    """
    Attributes:
        mc: monte carlo numbers
                 dimension (?)
        tstep : time step
                 dimension (1)
        steps: total steps for tracking
                 dimension (1)
        Emin: minimium energy
        Emax: maximum energy
    """
    def __init__(self,mc,tstep,steps,Emin,Emax,gfile):
        self.mc = mc
        self.tstep = tstep
        self.steps = steps
        self.Emin = Emin
        self.Emax = Emax
        self.gfile = gfile

class result(object):
    def __init__(self,R_birth,P_birth,E_birth):
        self.R_birth = R_birth
        self.P_birth = P_birth
        self.E_birth = E_birth


