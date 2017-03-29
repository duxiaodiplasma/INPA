"""
;+
; NAME:
;     bfield
; PURPOSE:
      calcualte magnetic field strength from G-file
; EXPLANATION:
;
; CALLING SEQUENCE:
	geqdsk.py   (read gfile)
	util_tokamak.py	(tools for reading gfile)
#;
; OPTIONAL INPUT:

; OPTIONAL KEYWORD INPUT:
;
; OPTIONAL KEYWORD OUTPUTS:

; PROCEDURE:

; EXAMPLE:

; RESTRICTIONS:
;
; PROCEDURES USED:
;
; FIRST USABLE VERSION
	10-06-2016 by X.D.DU
; REVISION HISTORY:
"""

import geqdsk
import lib
import scipy
import numpy as np
from scipy import interpolate
import sys

def prepare(fn):

   # read gfile
   res = geqdsk.read(fn)

   nr = res['nw']
   nz = res['nh']
   r = res['r']
   z = res['z']
   psirz = res['psirz']
   rc = res['rcentr']
   bc = res['bcentr']

   # first deravitive of the psi in RZ direction
   dpsi_dr = np.zeros([nr,nz])
   dpsi_dz = np.zeros([nr,nz])

   for i in range(1, nr-1):
       for j in range(1, nz-1):
           dpsi_dr[i,j] = (psirz[i+1,j]-psirz[i-1,j])/(r[i+1,j]-r[i-1,j])
           dpsi_dz[i,j] = (psirz[i,j+1]-psirz[i,j-1])/(z[i,j+1]-z[i,j-1])

   # 2D interpolation for dpsirz
   rr = np.zeros(nr)
   zz = np.zeros(nz)
   for i in range(0, nr):
       rr[i] = r[i, 0]
   for i in range(0, nz):
       zz[i] = z[0, i]

   fr = interpolate.interp2d(zz, rr, dpsi_dr, kind='cubic')
   fz = interpolate.interp2d(zz, rr, dpsi_dz, kind='cubic')

   sibry = res['sibry']
   simag = res['simag']
   psirz = np.array((psirz-simag)/(sibry-simag))
   fpsi = interpolate.interp2d(zz, rr, psirz, kind='cubic')
   return fr,fz,fpsi,res,rr,zz,psirz


def brzt(lr,lz,phi,ini):
   fr = ini.fr
   fz = ini.fz
   rc = ini.equ['rcentr']
   bc = ini.equ['bcentr']

   # local magnetic field strength
   l_dpsir = fr(lz,lr)
   l_dpsiz = fz(lz,lr)

   br = -1.0*-(1/lr)*l_dpsiz
   bz = -1.0* (1/lr)*l_dpsir # -1.0 for unify the defination
			     # of the plus current direction in DIII-D,$
			     # which is CCW direction for plus
   bt = np.asarray(rc*bc/lr)
   # print(br,bz,bt)

   # rz coordinate to xyz coordinate
   bx = -1.0*bt*np.cos(phi)-br*np.sin(phi)
   by = -1.0*bt*np.sin(phi)+br*np.cos(phi)
   bz = np.copy(bz)

   return np.concatenate([bx, by, bz])


