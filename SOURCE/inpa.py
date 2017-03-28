"""
;+
; NAME:
;     INPA SIMULATION
; PURPOSE:
      This Monte Carlo code calculates the orbit from carbon foil to the scintillator
; EXPLANATION:
      The neutral are randomly generated from plasma following the rules:
      (1) The energy of the neutral is randomly generated from 0keV to 80keV
      (2) The direction of neutral flying is randomly vector between one point
	  in aperture and the other in carbon foil.
      (3) The initial velocity [vx,vy,vz] is then determined.
      (4) The larmor motion of the particle is then traced until it hits
          the scintillator, i.e.,vz(t0)*vz(t0+dt)<0.
;
; CALLING SEQUENCE:
        def geometry() in design.py is depenedent file
#;
; OPTIONAL INPUT:
;
; OPTIONAL KEYWORD INPUT:
;
; OPTIONAL KEYWORD OUTPUTS:

; PROCEDURE:

; EXAMPLE:
       run inpa.py
; RESTRICTIONS:
;
; PROCEDURES USED:
;
; REVISION HISTORY:
      first usable version is inpa_v2.py
		by DU,Xiaodi on 08-09-2016

      change the design function purely 3D
	     design.geometry_3d(la, lb, lc)
		by DU,Xiaodi on 08-16-2016

      add a routine to generate the input file of NPA GEOMETRY for FIDASIM
	     output4_fidasim.input(p1,p2,p3,p4,f1,f2,f3,f4,w,h)
		by DU,Xiaodi on 11-3-2016

"""
#import sys
#sys.path.insert(0,'/home/duxiaodi/inpa/SOURCE')
#sys.path.insert(0,'/home/duxiaodi/inpa/GEOMETRY')
#sys.path.insert(0,'/home/duxiaodi/inpa/PLOT')
import creatobj
import lib
import geqdsk
import bfield
import specif
import output4_fidasim
import read_hd5
import particle_launcher
import numpy as np
import scipy

'''particle location at t+dt'''
def track(s_t0, v, tstep):
    s_t1 = np.zeros(3)
    s_t1[0] = s_t0[0] + v[0]*tstep
    s_t1[1] = s_t0[1] + v[1]*tstep
    s_t1[2] = s_t0[2] + v[2]*tstep
    return s_t1

'''accerlation by lorentz force at t+dt'''
def acce(v, Bt):
    a = np.zeros(3)
    charge = 1.6*(1e-19)
    mass_p = 1.67*(1e-27)
    mass_d = 2*mass_p
    a = charge/mass_d*np.cross(v,Bt)
    return a

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
def initial_pos(geo,ini):

    p1 = geo.pinhole[0,:]
    p2 = geo.pinhole[1,:]
    p3 = geo.pinhole[2,:]
    p4 = geo.pinhole[3,:]

    #ifoil = np.random.randint(10)
    ifoil = 0
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

'''
    initial energy of neutrals
'''
def initial_E(ini):
    # initial energy [keV]
    Eth  = ini.Emin
    Enbi = ini.Emax
    E_energy =Eth + (Enbi-Eth)*np.random.random(1) #keV

    return E_energy


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


'''--------------'''
'''MAIN Program'''
'''--------------'''
def main(geo,ini):
    # --- INPUT ZONE ---
    # number of monte carlo particles
    mc_num = ini.mc
    tstep = ini.tstep
    steps = ini.steps

    # inital pitch
    P_ini = np.zeros(mc_num)

    # inital energy
    E_ini = np.zeros(mc_num)

    # record everthing in history of one particle
    s_hist = np.zeros((steps,3))
    s_histx = np.zeros((mc_num,steps,3))
    a_hist = np.zeros((steps,3))
    v_hist = np.zeros((steps,3))
    bt_hist = np.zeros((steps,3))

    # birth sightline
    birth_sl = np.zeros((mc_num,2,3))

    # position of striking points
    strike_pos = np.zeros((mc_num,3))

    # light emission on scintillator
    light = np.zeros(mc_num)

    # inital R position
    R_ini = np.zeros((mc_num,2))

    for mc in np.arange(mc_num):

       # initial energy
       E_energy = initial_E(ini)

       # generate initail values for orbit tracing
       n0_p, n0_f, bt, r_loc, pitch_loc \
    	 = initial_pos(geo,ini)

       # initial velocity in xyz coordinate
       v_ini = initial_v(E_energy,n0_p,n0_f)

       # birth place at the pinhole and foil
       birth_sl[mc,0,:] = np.array(n0_p)
       birth_sl[mc,1,:] = np.array(n0_f)

       # birth place of R, pitch, and Energy
       R_ini[mc,:] = r_loc
       P_ini[mc] = pitch_loc
       E_ini[mc] = E_energy

       # -----------------------
       #     MAIN CALCULATION  |
       # -----------------------

       # trace orbit from foil
       # s is the position in XYZ coordinate
       s = n0_f
       # v is the initial velocity in XYZ coordinate
       v = v_ini

       # trace the Larmor motion
       for i in np.arange(steps):

           # calculate the orbit
           s = track(s,v,tstep)

           # 20170327 calculate local bt
           sr,sz,sphi = lib.xyz_to_rzphi(s[0],s[1],s[2])
           bt = bfield.brzt(sr,sz,sphi,ini)

           a = acce(v, bt)
           v = velo(v, a, tstep)
           s_hist[i,:] = s
           s_histx[mc,i,:] = s
           a_hist[i,:] = a
           v_hist[i,:] = v
           bt_hist[i,:] = bt

           # Judge whether hit the sintillator or not
           if i > 900:
              hitornot1 = \
              geo.lphoreq[0]*s_hist[i-1,0]+ \
              geo.lphoreq[1]*s_hist[i-1,1]+ \
              geo.lphoreq[2]*s_hist[i-1,2]+ \
              geo.lphoreq[3]

              hitornot2 = \
              geo.lphoreq[0]*s_hist[i,0]+ \
              geo.lphoreq[1]*s_hist[i,1]+ \
              geo.lphoreq[2]*s_hist[i,2]+ \
              geo.lphoreq[3]

              # if two neighbour points are locating at
              # two sides of the foil plane, it will be
              # considered as that one collision occurs
              if hitornot1*hitornot2<0:
                 strike_pos[mc,:]  = s

#             if switch_real_img == 1:
#                light[mc] = E_energy

                 break;

    # structure for output
    res = creatobj.result(R_ini,-1.0*np.cos(P_ini/180.*np.pi),E_ini)
    res.hitpoint = strike_pos
    res.birthsl = birth_sl
    return res



