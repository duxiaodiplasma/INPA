import sys
sys.path.insert(0,'/home/duxiaodi/inpa/SOURCE')
sys.path.insert(0,'/home/duxiaodi/inpa/GEOMETRY')
sys.path.insert(0,'/home/duxiaodi/inpa/PLOT')
import creatobj
import bfield
import numpy as np
import inpa
import inpa_ode
import doplot
import geometry

fpath = '/home/duxiaodi/inpa/GFILE/159243C04_801.geq'
#  -- INPUT --     mc,  tstep, steps, Eini, Emax, gfile
ini = creatobj.ini(3000, 1e-11, 10000,  20,  80,  fpath)

# -- EQUILIBRIUM --
ini.fr,ini.fz,ini.fpsi,ini.equ,ini.r,ini.z,ini.psirz \
= bfield.prepare(ini.gfile)

# -- INPA GEOMETRY --
geo = geometry.main()

# -- calculate striking point --
res = inpa_ode.main(geo,ini)

# calculuate the resolution due to the finite beam width effect
res = resolution.main(ini,res)







