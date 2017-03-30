import numpy as np
import math


''' judge the point P is in the triangular in ABC '''
''' if being inside, return 1, if being outside, return 0'''
def point_in_tri(P, A, B, C):

    # Compute vectors
    v0 = C - A
    v1 = B - A
    v2 = P - A
    # compute dot products
    dot00 = np.dot(v0, v0)
    dot01 = np.dot(v0, v1)
    dot02 = np.dot(v0, v2)
    dot11 = np.dot(v1, v1)
    dot12 = np.dot(v1, v2)

    # Compute barycentric coordinates
    invDenom = 1. / (dot00 * dot11 - dot01 * dot01)
    u = (dot11 * dot02 - dot01 * dot12) * invDenom
    v = (dot00 * dot12 - dot01 * dot02) * invDenom

    # Check if point is in triangle
    if u >= 0 and v >= 0 and u + v < 1:
       return 1
    else:
       return 0


''' transfer X,Y,Z coordinate to R, Z, Phi coordinate'''
def xyz_to_rzphi(x,y,z):

    r = np.sqrt(x**2+y**2)
    phi = -1.0*np.arctan(x/y)
                      # note that the +Bt in DIII-D is C.W. direction
                      # multiply -1 to change the coordinate from LHS to RHS
                      # phi rotates also in C.C.W direction
                      # the same coordinate with TRIP3D and M3d-c1
    return r, z, phi


'''inch to meter'''
def inchtometer(inch):
    meter = 0.0254*inch
    return meter


'''intersections of any line and cycle centered at (0,0)'''
def line_to_circle(p1,p2,r0):
    k = (p2[1] - p1[1]) / (p2[0]- p1[0])
    b = p2[1] - k*p2[0]
    #solution
    xs1 = (-2*k*b - np.sqrt((2*k*b)**2 - 4*(k**2+1)*(b**2-r0**2)))/(2*(k**2+1))
    xs2 = (-2*k*b + np.sqrt((2*k*b)**2 - 4*(k**2+1)*(b**2-r0**2)))/(2*(k**2+1))
    ys1 = k*xs1 + b
    ys2 = k*xs2 + b
    return [xs1,ys1],[xs2,ys2]

'''rotate one vector at angle of theta'''
def rotate_vec(vec,theta):
    # rotating angle
    rotate = np.array([[np.cos(theta), -np.sin(theta)], [np.sin(theta),  np.cos(theta)]])
    vec_new = np.dot(rotate,vec)
    return vec_new

'''unit vector'''
def unit_vector(vector):
    # Returns the unit vector of the vector.  """
    return vector / np.linalg.norm(vector)

'''calculate angle between two vectors'''
def angle_between(v1, v2):
    # Returns the angle in radians between vectors 'v1' and 'v2'::
    v1_u = unit_vector(v1)
    v2_u = unit_vector(v2)
    return np.arccos(np.clip(np.dot(v1_u, v2_u), -1.0, 1.0))

'''define a line'''
def line(p1, p2):
    A = (p1[1] - p2[1])
    B = (p2[0] - p1[0])
    C = (p1[0]*p2[1] - p2[0]*p1[1])
    return A, B, -C

'''intersection of two lines'''
def intersection(L1, L2):
    D  = L1[0] * L2[1] - L1[1] * L2[0]
    Dx = L1[2] * L2[1] - L1[1] * L2[2]
    Dy = L1[0] * L2[2] - L1[2] * L2[0]
    if D != 0:
        x = Dx / D
        y = Dy / D
        return x,y
    else:
        return False

''' intersection function '''
def isect_line_plane_v3(p0, p1, p_co, p_no, epsilon=1e-6):
    """
    p0, p1: define the line
    p_co, p_no: define the plane:
        p_co is a point on the plane (plane coordinate).
        p_no is a normal vector defining the plane direction; does not need to be normalized.

    return a Vector or None (when the intersection can't be found).
    """

    u = sub_v3v3(p1, p0)
    dot = dot_v3v3(p_no, u)

    if abs(dot) > epsilon:
        # the factor of the point between p0 -> p1 (0 - 1)
        # if 'fac' is between (0 - 1) the point intersects with the segment.
        # otherwise:
        #  < 0.0: behind p0.
        #  > 1.0: infront of p1.
        w = sub_v3v3(p0, p_co)
        fac = -dot_v3v3(p_no, w) / dot
        u = mul_v3_fl(u, fac)
        return add_v3v3(p0, u)
    else:
        # The segment is parallel to plane
        return None

# generic math functions for intersection point 3D

def add_v3v3(v0, v1):
    return (
        v0[0] + v1[0],
        v0[1] + v1[1],
        v0[2] + v1[2],
        )

def sub_v3v3(v0, v1):
    return (
        v0[0] - v1[0],
        v0[1] - v1[1],
        v0[2] - v1[2],
        )

def dot_v3v3(v0, v1):
    return (
        (v0[0] * v1[0]) +
        (v0[1] * v1[1]) +
        (v0[2] * v1[2])
        )

def len_squared_v3(v0):
    return dot_v3v3(v0, v0)

def mul_v3_fl(v0, f):
    return (
        v0[0] * f,
        v0[1] * f,
        v0[2] * f,
        )

def rotation_matrix(axis, theta):
    """
    Return the rotation matrix associated with counterclockwise rotation about
    the given axis by theta radians.
    """
    axis = np.asarray(axis)
    axis = axis/math.sqrt(np.dot(axis, axis))
    a = math.cos(theta/2.0)
    b, c, d = -axis*math.sin(theta/2.0)
    aa, bb, cc, dd = a*a, b*b, c*c, d*d
    bc, ad, ac, ab, bd, cd = b*c, a*d, a*c, a*b, b*d, c*d
    return np.array([[aa+bb-cc-dd, 2*(bc+ad), 2*(bd-ac)],
                     [2*(bc-ad), aa+cc-bb-dd, 2*(cd+ab)],
                     [2*(bd+ac), 2*(cd-ab), aa+dd-bb-cc]])

