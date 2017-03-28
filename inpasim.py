import sys
sys.path.insert(0,'/home/duxiaodi/inpa/SOURCE')
sys.path.insert(0,'/home/duxiaodi/inpa/GEOMETRY')
sys.path.insert(0,'/home/duxiaodi/inpa/PLOT')
import creatobj
import bfield
import numpy as np
import inpa
import inpa_ode
import inpa_odepara
import doplot

# -- test various geometry --
import geometry
import geometry_failed3
import geometry_upgrade
import geometry_upgrade2
import geometry_upgrade3

fpath = '/home/duxiaodi/inpa/GFILE/159243C04_801.geq'
#  -- INPUT --     mc,  tstep, steps, Eini, Emax, gfile
ini = creatobj.ini(5000, 1e-11, 10000,  20,  80,  fpath)
ini.charge = 1.6*(1e-19) # D charge
ini.mass = 2*1.67*(1e-27)# D mass

# -- EQUILIBRIUM --
ini.fr,ini.fz,ini.equ = bfield.prepare(ini.gfile)

# -- INPA GEOMETRY --
geo = geometry_failed3.main()

# -- CALL MAIN --
res = inpa_ode.main(geo,ini)





