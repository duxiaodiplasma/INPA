"""
;+
; NAME:
  FBW.py (finite beam width)

; PURPOSE:
  to calcualte the crossing point between the sightlines 
  and beam. Include the beam width of 330L in DIII-D tokamak

; EXPLANATION:
;
; CALLING SEQUENCE:
;
; INPUT:
  input the two points on the sightline
;
; OPTIONAL INPUT:
;
; OPTIONAL KEYWORD INPUT:
;
; OUTPUT:
  crossing point with the beam
  transfer these points to RZ plane

; OPTIONAL OUTPUTS:

; PROCEDURE:
  (1) calcualte the crossing points with four beam planes
  (2) judge the points with the beam or not
      note that here to judge the point in triangular or not
	source: http://blackpawn.com/texts/pointinpoly/default.html


; EXAMPLE:
  rz_seg = fbw(p1,p2)
; RESTRICTIONS:
;
; PROCEDURES USED:
;
; REVISION HISTORY:
"""

import numpy as np
import bfield
import lib

def fbw(p0, p1):

#    p0 = np.array([-2.17180301,  0.5827061 , -0.71889309])
#    p1 = np.array([-2.27512821,  0.49102035, -0.79716202])
    
    # source location
    bp1 = lib.inchtometer(np.array([-78.9, 310.1, 0]))
    
    # beam crowwover point
    bp2 = lib.inchtometer(np.array([-55.3, 92.1, 0]))
        
    # unit vector of center of beamline
    vec_bp2 = lib.unit_vector(bp2-bp1)
    
    # half width of the beamline [m]
    awidw = 8.85*(1e-2) 
    
    # half height of the beamline [m]
    awidh = 24.*(1e-2) 
    
    # upper point of bp2
    bp2_upper = np.array(bp2)
    bp2_upper[2] = bp2[2] + awidh
    
    bp2_lower = np.array(bp2)
    bp2_lower[2] = bp2[2] - awidh
    
    # unit vector of bp2_upper and bp2
    vec_bp2_upper = lib.unit_vector(bp2_upper - bp2)
    vec_bp2_lower = -1.0*vec_bp2_upper
    
    # unit vector to the beam right
    # Right: viewing from source point to crossover point
    vec_bp2_right = np.cross(vec_bp2, vec_bp2_upper)
    vec_bp2_left = -1.0 * vec_bp2_right
    
    # 1------2
    # |plasma|
    # 4------3
    
    bp21 = bp2 + awidw*vec_bp2_left  + awidh*vec_bp2_upper
    bp22 = bp2 + awidw*vec_bp2_right + awidh*vec_bp2_upper
    bp23 = bp2 + awidw*vec_bp2_right + awidh*vec_bp2_lower
    bp24 = bp2 + awidw*vec_bp2_left  + awidh*vec_bp2_lower
   
    print bp21,bp22,bp23,bp24
    
    # 5------6
    # |plasma|
    # 8------7
    
    bp25 = bp21 + 3.*vec_bp2
    bp26 = bp22 + 3.*vec_bp2
    bp27 = bp23 + 3.*vec_bp2
    bp28 = bp24 + 3.*vec_bp2
    
    print bp25,bp26,bp27,bp28

    # find the intersection point
    x1562 = lib.isect_line_plane_v3(p0,p1,bp21,vec_bp2_upper)
    x4873 = lib.isect_line_plane_v3(p0,p1,bp24,vec_bp2_upper)
    x2673 = lib.isect_line_plane_v3(p0,p1,bp22,vec_bp2_left )
    x1584 = lib.isect_line_plane_v3(p0,p1,bp21,vec_bp2_left )
    
    # check point  
    xx = np.zeros((2,3))
    i = 0 
    if lib.point_in_tri(x1562,bp21,bp22,bp25)==1 or \
       lib.point_in_tri(x1562,bp22,bp25,bp26)==1:
       xx[i,:] = x1562
       i = i + 1
    if lib.point_in_tri(x4873,bp24,bp23,bp28)==1 or \
       lib.point_in_tri(x4873,bp27,bp23,bp28)==1:
       xx[i,:] = x4873
       i = i + 1
    if lib.point_in_tri(x2673,bp22,bp23,bp26)==1 or \
       lib.point_in_tri(x2673,bp27,bp23,bp26)==1:
       xx[i,:] = x2673
       i = i + 1
    if lib.point_in_tri(x1584,bp21,bp25,bp28)==1 or \
       lib.point_in_tri(x1584,bp21,bp28,bp24)==1:
       xx[i,:] = x1584
    
    # transfer to RZ coordinate
    segs = 21
    vec_xx = lib.unit_vector(xx[1,:]-xx[0,:])
    dis_xx = np.linspace(0, np.linalg.norm(xx[1,:]-xx[0,:]), segs)
    
    xx_seg = np.zeros((segs,3))
    rz_seg = np.zeros((segs,3))
    
    for i in np.arange(segs):
        xx_seg[i,:] = xx[0,:]+dis_xx[i]*vec_xx
	rz_seg[i,:] = lib.xyz_to_rzphi(xx_seg[i,0], xx_seg[i,1], xx_seg[i,2]) 
#   
    return rz_seg















