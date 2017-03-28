import numpy as np
import lib
import bfield
#from mpl_toolkits.mplot3d import Axes3D
#from matplotlib import cm
#import matplotlib.pyplot as plt

def main(p0,p1,ini):

# --- input zone ---
   # source location
   bp1 = lib.inchtometer(np.array([-78.9, 310.1, 0]))

   # beam crowwover point
   bp2 = lib.inchtometer(np.array([-55.3, 92.1, 0]))

   # define normal vector of beam plane
   bp_normal = np.cross(bp2 - bp1, np.asarray([0.,0.,1.]))

   # p1 at foil, p0 at aperture,
   r_loc = lib.isect_line_plane_v3(p0,p1,bp1,bp_normal,epsilon=1e-6)

   # b_field at r_loc
   lr,lz,lphi = lib.xyz_to_rzphi(r_loc[0],r_loc[1],r_loc[2])
   b_loc = bfield.brzt(lr,lz,lphi,ini)

   # pitch angle between sightline and magnetic field
   pitch_loc = lib.angle_between(p1-p0,b_loc)/np.pi*180.

   return np.asarray([lr,lz]), pitch_loc




