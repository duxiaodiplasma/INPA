import sys
sys.path.insert(0,'/home/duxiaodi/inpa/SOURCE')
sys.path.insert(0,'/home/duxiaodi/inpa/GEOMETRY')
sys.path.insert(0,'/home/duxiaodi/inpa/PLOT')
import creatobj
import bfield
import numpy as np
import inpa_odepara
import doplot
import geometry
import resolution
from joblib import Parallel, delayed


fpath = '/home/duxiaodi/inpa/GFILE/159243C04_801.geq'
#  -- INPUT --     mc,  tstep, steps, Eini, Emax, gfile
ini = creatobj.ini(5000, 1e-11, 10000,  20,  80,  fpath)
ini.charge = 1.6*(1e-19) # D charge
ini.mass = 2*1.67*(1e-27)# D mass
ini.cores = 10 # CPU CORES

# -- EQUILIBRIUM --
ini.fr,ini.fz,ini.fpsi,ini.equ,ini.r,ini.z,ini.psirz \
= bfield.prepare(ini.gfile)

# -- INPA GEOMETRY --
geo = geometry.main()

# --- RNG --
ini = inpa_odepara.RNG(ini,geo)

# -- PREPARE OUTPUT --
res = creatobj.result(ini.R_birth,-1.0*np.cos(ini.P_birth),ini.E_birth)
res.hitpoint = np.zeros((ini.mc,3))
res.birth_sl = np.zeros((ini.mc,2,3))

# ---------------------- #
#  parallel calculation  #
# ---------------------- #
if __name__ == '__main__':
   output = Parallel(n_jobs=ini.cores,verbose=5)(delayed(inpa_odepara.PARA_mc)(ini,geo,i) \
            for i in range(0,ini.mc))

for i in range(0,ini.mc):
    res.hitpoint[i,:] = output[i]
    res.birth_sl[i,0,:] = ini.n0_p[i,:]
    res.birth_sl[i,1,:] = ini.n0_f[i,:]





