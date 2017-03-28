import sys
sys.path.insert(0,'/home/duxiaodi/inpa/SOURCE')
sys.path.insert(0,'/home/duxiaodi/inpa/GEOMETRY')
sys.path.insert(0,'/home/duxiaodi/inpa/PLOT')
import creatobj
import bfield
import numpy as np
import inpa
import doplot

# -- test various geometry --
import geometry
import geometry_failed3
import geometry_upgrade
import geometry_upgrade2
import geometry_upgrade3

fpath = '/home/duxiaodi/inpa/INPUT/159243C04_801.geq'
#  -- INPUT --     mc,  tstep, steps, Eini, Emax, gfile
ini = creatobj.ini(100, 1e-11, 10000,  20,  80,  fpath)

# -- EQUILIBRIUM --
ini.fr,ini.fz,ini.equ = bfield.prepare(ini.gfile)

# -- INPA GEOMETRY --
geo = geometry_upgrade.main()

# -- CALL MAIN --
res = inpa.main(geo,ini)

#



