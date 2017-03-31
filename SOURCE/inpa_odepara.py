"""
+
 NAME:
     INPA SIMULATION
 PURPOSE:
     This Monte Carlo code calculates the orbit from carbon foil to the scintillator
 EXPLANATION:
     The neutral are randomly generated from plasma following the rules:
     (1) The energy of the neutral is randomly generated from 0keV to 80keV
     (2) The direction of neutral flying is randomly vector between one point
     in aperture and the other in carbon foil.
     (3) The initial velocity [vx,vy,vz] is then determined.
     (4) The larmor motion of the particle is then traced until it hits
         the scintillator, i.e.,vz(t0)*vz(t0+dt)<0.

 CALLING SEQUENCE:
       def geometry() in design.py is depenedent file

 OPTIONAL INPUT:

 OPTIONAL KEYWORD INPUT:

 OPTIONAL KEYWORD OUTPUTS:

 PROCEDURE:

 EXAMPLE:
      run inpa.py
 RESTRICTIONS:

 PROCEDURES USED:

 REVISION HISTORY:
     first usable version is inpa_v2.py
   	by DU,Xiaodi on 08-09-2016

     change the design function purely 3D
        design.geometry_3d(la, lb, lc)
   	by DU,Xiaodi on 08-16-2016

     add a routine to generate the input file of NPA GEOMETRY for FIDASIM
        output4_fidasim.input(p1,p2,p3,p4,f1,f2,f3,f4,w,h)
   	by DU,Xiaodi on 11-3-2016

"""
import sys
sys.path.insert(0,'/home/duxiaodi/inpa/SOURCE')
sys.path.insert(0,'/home/duxiaodi/inpa/GEOMETRY')
sys.path.insert(0,'/home/duxiaodi/inpa/PLOT')
import creatobj
import lib
import geqdsk
import bfield
import specif
import output4_fidasim
import read_hd5
import particle_launcher
import numpy as np
from scipy.integrate import odeint

'''particle location at t+dt'''
def track(s_t0, v, tstep):
    s_t1 = np.zeros(3)
    s_t1[0] = s_t0[0] + v[0]*tstep
    s_t1[1] = s_t0[1] + v[1]*tstep
    s_t1[2] = s_t0[2] + v[2]*tstep
    return s_t1


'''velocity at t+dt'''
def velo(v_t0, a, tstep):
    v_t1 = np.zeros(3)
    v_t1[0] = v_t0[0] + a[0]*tstep
    v_t1[1] = v_t0[1] + a[1]*tstep
    v_t1[2] = v_t0[2] + a[2]*tstep
    return v_t1

'''unit vector'''
def unit_vector(vector):
    # Returns the unit vector of the vector
    return vector / np.linalg.norm(vector)

'''generate one random point between p0 and p1 in the line'''
def random_3d(p0,p1):
    p3 = p0 + np.random.random(1)*(p1-p0)
    return p3


'''initial position on scintillator'''
def initials(geo,ini):

    p1 = geo.pinhole[0,:]
    p2 = geo.pinhole[1,:]
    p3 = geo.pinhole[2,:]
    p4 = geo.pinhole[3,:]

    ifoil = np.random.randint(geo.foil.shape[0])
    f1 = geo.foil[ifoil,0,:]
    f2 = geo.foil[ifoil,1,:]
    f3 = geo.foil[ifoil,2,:]
    f4 = geo.foil[ifoil,3,:]

    # random position at pinhole
    n0_p_12 = random_3d(p1,p2)
    n0_p_34 = random_3d(p3,p4)
    # neutral on pinhole plane
    n0_p = random_3d(n0_p_12,n0_p_34)

    # random position at foil
    n0_f_12 = random_3d(f1,f2)
    n0_f_34 = random_3d(f3,f4)
    # neutral on foil plane
    n0_f = random_3d(n0_f_12,n0_f_34)

    # pitch angle and energy monitored
    r_loc, pitch_loc = specif.main(n0_p, n0_f,ini)

    # transfer xyz to rzphi
    lr, lz, lphi = lib.xyz_to_rzphi(n0_f[0],n0_f[1],n0_f[2])

    # local magnetic field at the foil
    bt = bfield.brzt(lr,lz,lphi,ini)

    return n0_p, n0_f, bt, r_loc, pitch_loc


'''initial velocity of neutrals'''
def initial_v(E_energy,n0_p,n0_f):
    # initial velocity
    mass = 2*1.67*(1e-27) #deuteriam mass
    abs_v = np.sqrt(2*E_energy*1e3*(1.6e-19)/mass)

    v_ini = np.zeros(3)
    v_ini[0] = n0_f[0]-n0_p[0]
    v_ini[1] = n0_f[1]-n0_p[1]
    v_ini[2] = n0_f[2]-n0_p[2]

    v_ini = abs_v*unit_vector(v_ini)

    return v_ini

def ode_equ(s,t,fr,fz,rc,bc,charge,mass):
    x,y,z = s[0],s[1],s[2]
    vx,vy,vz = s[3],s[4],s[5]

    # 20170327 calculate local bt

    # xyz to rzphi
    lr = np.sqrt(x**2+y**2)
    lz = z
    phi = -1.0*np.arctan(x/y)
                      # note that the +Bt in DIII-D is C.W. direction
                      # multiply -1 to change the coordinate from LHS to RHS
                      # phi rotates also in C.C.W direction
                      # the same coordinate with TRIP3D and M3d-c1

    # local magnetic field strength
    l_dpsir = fr(lz,lr)
    l_dpsiz = fz(lz,lr)

    br = -1.0*-(1/lr)*l_dpsiz
    bz = -1.0* (1/lr)*l_dpsir # -1.0 for unify the defination
                  # of the plus current direction in DIII-D,$
                  # which is CCW direction for plus
    bt = np.asarray(rc*bc/lr)

    # rz coordinate to xyz coordinate
    bx = -1.0*bt*np.cos(phi)-br*np.sin(phi)
    by = -1.0*bt*np.sin(phi)+br*np.cos(phi)
    b = np.concatenate([bx,by,bz])

    ax,ay,az = charge/mass*np.cross([vx,vy,vz],b)

    return np.array([vx,vy,vz,ax,ay,az])

# -- RANDOM NUMBER GENERATOR FOR PARALLEL CALCULATION --
def RNG(ini,geo):
    ini.n0_p = np.zeros((ini.mc,3))
    ini.n0_f = np.zeros((ini.mc,3))
    ini.B_birth = np.zeros((ini.mc,3))
    ini.R_birth = np.zeros((ini.mc,2))
    ini.P_birth = np.zeros(ini.mc)
    ini.E_birth = np.random.uniform(ini.Emin,ini.Emax,size=ini.mc)

    for i in range(0,ini.mc):
        ini.n0_p[i,:], ini.n0_f[i,:], ini.B_birth[i,:], ini.R_birth[i], ini.P_birth[i] \
         = initials(geo,ini)
    return ini

# -- PARALLEL CALCULATION -- #
def PARA_mc(ini,geo,i):

    # initail conditions for orbit tracing
    n0_p, n0_f, E_ini \
 	= ini.n0_p[i,:], ini.n0_f[i,:], ini.E_birth[i]

    # initial velocity in xyz coordinate
    v_ini = initial_v(E_ini,n0_p,n0_f)

    # -----------------------
    #     solve ode         |
    # -----------------------

    ini_y0 = np.concatenate((n0_f,v_ini))

    # time slices for integrator
    tstep = np.copy(ini.tstep)
    steps = np.copy(ini.steps)
    tsol = np.arange(0,tstep*steps,tstep)

    #ODE INTEGRATOR
    sol = odeint(ode_equ,ini_y0,tsol,
          args=(ini.fr,ini.fz,ini.equ['rcentr'],ini.equ['bcentr'],ini.charge,ini.mass))

    # determine hitting points
    tmplower =  np.sign(
       geo.lphoreq[0]*sol[:,0] +
       geo.lphoreq[1]*sol[:,1] +
       geo.lphoreq[2]*sol[:,2] +
       geo.lphoreq[3]
       )
    u_or_l = np.where(np.diff(tmplower))[0]
    if len(u_or_l) > 1:
       hit_index = u_or_l[1]
       hit_point = sol[hit_index,0:3]
    # if not hit lower one than check upper one
    else:
       tmpupper =  np.sign(
       geo.uphoreq[0]*sol[:,0] +
       geo.uphoreq[1]*sol[:,1] +
       geo.uphoreq[2]*sol[:,2] +
       geo.uphoreq[3]
       )
       hit_index = np.where(np.diff(tmpupper))[0][1]
       hit_point = sol[hit_index,0:3]

    return hit_point

