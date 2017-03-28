"""
;+
; NAME:
;     INPA Slightline Design
; PURPOSE:
      This code is to design the slighlines for the INPA
; EXPLANATION:
;
; CALLING SEQUENCE:
;
; INPUT:
      (1) dimension of the aperture
      (2) center of the aperture
      (3) wanted region to monitor in plasma
      (4) NBI center beamline
;
; OPTIONAL INPUT:
;
; OPTIONAL KEYWORD INPUT:
;
; OUTPUT:
      4 pinhole corners
      4 foil corners
      sightlines

; OPTIONAL OUTPUTS:
      pitch angle of each sightline

; PROCEDURE:

; EXAMPLE:
       design.geometry(lx,ly,lz)
; RESTRICTIONS:
;
; PROCEDURES USED:
;
; REVISION HISTORY:
(1) first usable version is design.py
                by DU,Xiaodi on 08-09-2016

(2) Moidify the code for arbitary use
    in 3D cartisian coordinate
                by DU,Xiaodi on 08-16-2016
"""

import bfield
import lib
import numpy as np


def geometry_3d(la,lb,lc):

# dimension of the pinhole (half)
#la = 0.005/1. # height
#lb = 0.005/1.  # width
#lc = 0.1  # pinhole to foil
#debug_code = 1
#if debug_code > 0:
# --- input zone ---
   # source location
   bp1 = lib.inchtometer(np.array([-78.9, 310.1]))

   # beam crowwover point
   bp2 = lib.inchtometer(np.array([-55.3, 92.1]))

   # Plasma outboard boundary
   Rmidout = np.array(2.279)

   # Plasma inboard boundary
   Rmidin  = np.array(1.107)

   # magnetic axis
   Rm = np.array(1.734)

   npa = 'inpa_rm1'
  # npa = 'inpa_rm2'
   #npa = 'inpa_r0'

      #v2 0  -0.06
      #v3 -15  -0.06
      #v4 -30  -0.06
      #v5 -30  -0.5

   if npa == 'inpa_rm2':
      # define surface of pinhole
      # (share the same normal vector with foil)
      phi = 15
      surf_p1 = np.array([-1.8*np.cos(phi/180.*np.pi),1.8*np.sin(phi/180.*np.pi),-1.172])
      surf_p2 = np.array([-2.014*np.cos(phi/180.*np.pi),\
          2.014*np.sin(phi/180.*np.pi),-1.049])

      tang_vec = np.append(lib.rotate_vec([-np.cos(phi/180.*np.pi), np.sin(phi/180.*np.pi)],np.pi/2),0)

      # normal vector of suface aperture
      surf_normal = np.cross(surf_p1-surf_p2,tang_vec)

      #unit vectors
      v1 = lib.unit_vector(surf_p1-surf_p2)
      v2 = lib.unit_vector(tang_vec)
      v3 = lib.unit_vector(surf_normal)

      # pinhole center position
      pinhole = (surf_p1+surf_p2)/2.

   if npa == 'inpa_rm1':
      # define surface of pinhole
      # (share the same normal vector with foil)
      phi = 15
      surf_p1 = np.array([-2.136*np.cos(phi/180.*np.pi), 2.136*np.sin(phi/180.*np.pi),-0.995])
      surf_p2 = np.array([-2.364*np.cos(phi/180.*np.pi), 2.364*np.sin(phi/180.*np.pi),-0.436])

      tang_vec = np.append(lib.rotate_vec([-np.cos(phi/180.*np.pi), np.sin(phi/180.*np.pi)],np.pi/2),0)

      # normal vector of suface aperture
      surf_normal = np.cross(surf_p1-surf_p2,tang_vec)

      #unit vectors
      v1 = lib.unit_vector(surf_p1-surf_p2)
      v2 = lib.unit_vector(tang_vec)
      v3 = lib.unit_vector(surf_normal)

      # pinhole center position
      pinhole = (surf_p1+surf_p2)/2.

      # rotate this pinhole
      pinhole[0:2] = lib.rotate_vec(np.array(pinhole[0:2]),-0/180.*np.pi)
#      pinhole[2] = pinhole[2]-0.06 #v2,v3,v4,v5



   if npa == 'inpa_r0':
      phi = 15

      v3 = np.array([-Rmidout*np.cos(phi/180.*np.pi), Rmidout*np.sin(phi/180.*np.pi),0])
      v3 = lib.unit_vector(v3)

      # counter clockwise 90 degree
      v2 = lib.unit_vector(np.append(lib.rotate_vec(v3[0:2], np.pi/2),0))
      v1 = np.array([0,0,-1])

      # pinhole center position
      pinhole = np.array([-2.3271,0.410331,0])

# part I ends here

# -------------------
   # from detector see pinhole
   # four corner of the pinhole (clockwise)
   #  1 - - - 2
   #  |       |
   #  4 - - - 3
   pinhole_p = np.zeros((4,3))
   pinhole_p[0,:] = pinhole - la*v1 - lb*v2
   pinhole_p[1,:] = pinhole - la*v1 + lb*v2
   pinhole_p[2,:] = pinhole + la*v1 + lb*v2
   pinhole_p[3,:] = pinhole + la*v1 - lb*v2

   # center of the foil
   foil = pinhole + lc*v3

   # general equation for foil plane
   f_plane = np.zeros(4)
   f_plane[0:3] = v3
   f_plane[3] = -np.dot(v3,foil)

   # segments of the neutral source (nbi)
   num = 100
   source = np.zeros((num,3))
   radius_mid = np.zeros(num)

   # each sightline
   foil_p = np.zeros((num,4,3))
   pitch = np.zeros(num)
   vb = np.zeros((num,3))

   foil_p_center = np.zeros(3)
   dis_foil_beam = np.zeros(num)
   s_angle = np.zeros(num)

   j=0
   for i in np.linspace(Rmidout, Rmidin+0.1,num):
          # get the cross point A between nbi and bt
          cs1,cs2 = lib.line_to_circle(bp1,bp2,i)
          cs1 = np.append(cs1,0)
          source[j,:] = cs1
          radius_mid[j] = i

          # define a line
          foil_p[j,0,:] = lib.isect_line_plane_v3(cs1,pinhole_p[0],foil,v3,epsilon=1e-6)
          foil_p[j,1,:] = lib.isect_line_plane_v3(cs1,pinhole_p[1],foil,v3,epsilon=1e-6)
          foil_p[j,2,:] = lib.isect_line_plane_v3(cs1,pinhole_p[2],foil,v3,epsilon=1e-6)
          foil_p[j,3,:] = lib.isect_line_plane_v3(cs1,pinhole_p[3],foil,v3,epsilon=1e-6)

          # calculate the solid angle
          foil_p_center = (foil_p[j,0,:]+foil_p[j,1,:]+foil_p[j,2,:]+foil_p[j,3,:])/4.
          dis_foil_beam[j] = np.sqrt((foil_p_center[0]-cs1[0])**2+\
               			     (foil_p_center[1]-cs1[1])**2+\
	         		     (foil_p_center[2]-cs1[2])**2)
	  s_area = np.linalg.norm(np.cross(foil_p[j,0,:]-foil_p[j,1,:], foil_p[j,2,:]-foil_p[j,1,:]))
          s_angle[j] = s_area/(4*np.pi*dis_foil_beam[j]**2)

          # observed pitch angle
          r_loc = cs1
          lr = np.sqrt(r_loc[0]**2+r_loc[1]**2)
          lz = r_loc[2]
          lphi = -1.0*np.arctan(r_loc[0]/r_loc[1])
          vb[j,:] = bfield.brzt(lr,lz,lphi)

          vsl = lib.unit_vector(pinhole-cs1)
          pitch[j] = lib.angle_between(vb[j,:],vsl)/np.pi*180

          j = j + 1

   foil_bry = np.asarray([foil_p[num-1,0,:],foil_p[0,1,:],foil_p[0,2,:],foil_p[num-1,3,:]])
   geometry = [pinhole_p, foil_bry]
   bp = np.zeros((3,3))
   bp[0,:] = np.append(bp1,0)
   bp[1,:] = np.append(bp2,0)
   bp[2,:] = bp[1,:] +2.5*lib.unit_vector(bp[1,:]-bp[0,:])


   np.savez('./output/geometry',pinhole = geometry[0],foil = geometry[1])
   np.savez('./output/sightline',source=source,pinhole_p= pinhole_p,foil_p=foil_p
,bp = bp,vb = vb,s_angle=s_angle,s_area=s_area)

#
   return pinhole_p[0,:], \
          pinhole_p[1,:], \
          pinhole_p[2,:], \
          pinhole_p[3,:], \
          foil_p[num-1,0,:],\
          foil_p[0,1,:],  \
          foil_p[0,2,:],  \
          foil_p[num-1,3,:],\
          f_plane
