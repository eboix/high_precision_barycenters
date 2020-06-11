# This power diagram code is adapted from Devert Alexandre's code: https://gist.github.com/marmakoide/45d5389252683ae09c2df49d0548a627#file-laguerre-voronoi-2d-py.
# It is modified in the following ways:
# a) compatibility with python3
# b) some minor performance improvements are made
# c) infinite rays in the power diagram are ignored, because barycenter_min_oracle.py adds "virtual" points to the diagram to ensure that infinite rays do not matter for our code.
# d) the power diagrams are truncated to lie in a box

import itertools
import numpy as np
from scipy.spatial import ConvexHull
from sklearn.preprocessing import normalize

# ------------------ Code to truncate segments to box -------------------------------
def pt_in_box(pt, bounding_box):
    # Returns whether pt is in bounding_box
    x,y = pt
    (xmin, ymin), (xmax, ymax) = bounding_box
    return (xmin <= x and xmax >= x and ymin <= y and ymax >= y)

def line_intersection_with_vertical_line(pt1, slope, x):
    # Returns value of line given by pt1 and slope, at given x-coordinate
    if pt1[0] == x:
        return pt1[1], 0
    if slope[0] == 0:
        return None, None

    t = (x - pt1[0]) / slope[0]
    dy = t * slope[1]
    return pt1[1] + dy, t
        
    
def line_intersection_with_horizontal_line(pt1, slope, y):
    # Returns value of line given by pt1 and slope, at given y coordinate
    if pt1[1] == y:
        return pt1[0], 0
    if slope[1] == 0:
        return None, None

    t = (y - pt1[1]) / slope[1]
    dx = t * slope[0]
    return pt1[0] + dx, t

def get_segment_box_intersection_points(pt1, pt2, bounding_box):
    # Helper method: returns possible intersection points of line segment (pt1, pt2) with bounding_box
    xmin = bounding_box[0][0]
    ymin = bounding_box[0][1]
    xmax = bounding_box[1][0]
    ymax = bounding_box[1][1]

    slope = pt2 - pt1

    y_inter_xmin, xmint = line_intersection_with_vertical_line(pt1, slope, xmin)
    x_inter_ymin, ymint = line_intersection_with_horizontal_line(pt1, slope, ymin)
    y_inter_xmax, xmaxt = line_intersection_with_vertical_line(pt1, slope, xmax)
    x_inter_ymax, ymaxt = line_intersection_with_horizontal_line(pt1, slope, ymax)

    inter_pts = []
    curr_inter_pt = 0
    if not xmint is None:
        if ymin <= y_inter_xmin and y_inter_xmin <= ymax and 0 <= xmint and xmint <= 1:
           inter_pts.append((xmin, y_inter_xmin))
    if not xmaxt is None:
        if ymin <= y_inter_xmax and y_inter_xmax <= ymax and 0 <= xmaxt and xmaxt <= 1:
            inter_pts.append((xmax, y_inter_xmax))
    if not ymint is None:
        if xmin < x_inter_ymin and x_inter_ymin < xmax and 0 <= ymint and ymint <= 1:
            inter_pts.append((x_inter_ymin, ymin))
    if not ymaxt is None:
        if xmin < x_inter_ymax and x_inter_ymax < xmax and 0 <= ymaxt and ymaxt <= 1:
            inter_pts.append((x_inter_ymax, ymax))

    return np.asarray(inter_pts)

def truncate_finite_segment_to_box(pt1, pt2, bounding_box):
    # Return the finite segment (pt1, pt2) truncated to fit in bounding_box, or
    # return None, None if (pt1, pt2) is not in bounding_box.

    xmin = bounding_box[0][0]
    ymin = bounding_box[0][1]
    xmax = bounding_box[1][0]
    ymax = bounding_box[1][1]

    # First break into the 3 possible cases -- to avoid errors due to rounding.
    inbox1 = pt_in_box(pt1, bounding_box)
    inbox2 = pt_in_box(pt2, bounding_box)
    if inbox1 and inbox2:
        # Both endpoints are in the box
        return pt1, pt2
    elif not (inbox1 or inbox2):
        # Neither endpoint is in the box
        inter_pts = get_segment_box_intersection_points(pt1, pt2, bounding_box)
        if len(inter_pts) == 0:
            return None, None
        assert(len(inter_pts) == 2)
        return inter_pts[0], inter_pts[1]
    else:
        # Exactly one endpoint is in the box

        # Without loss of generality, pt1 is in box, and pt2 is not
        if inbox2:
            pt1, pt2 = pt2, pt1
        assert(pt_in_box(pt1, bounding_box))
        assert(not pt_in_box(pt2, bounding_box))

        # Compute the ONE intersection of the line segment with the box.
        slope = pt2 - pt1
        y_inter_xmin, xmint = line_intersection_with_vertical_line(pt1, slope, xmin)
        x_inter_ymin, ymint = line_intersection_with_horizontal_line(pt1, slope, ymin)
        y_inter_xmax, xmaxt = line_intersection_with_vertical_line(pt1, slope, xmax)
        x_inter_ymax, ymaxt = line_intersection_with_horizontal_line(pt1, slope, ymax)

        currintert = 2
        inter_pt = None, None
        if xmint is not None and xmint > 0 and xmint < currintert:
            currintert = xmint
            inter_pt = [xmin, y_inter_xmin]
        if ymint is not None and ymint > 0 and ymint < currintert:
            currintert = ymint
            inter_pt = [x_inter_ymin, ymin]
        if xmaxt is not None and xmaxt > 0 and xmaxt < currintert:
            currintert = xmaxt
            inter_pt = [xmax, y_inter_xmax]
        if ymaxt is not None and ymaxt > 0 and ymaxt < currintert:
            currintert = ymaxt
            inter_pt = [x_inter_ymax, ymax]

        if inter_pt[0] is None:
            return None, None
        assert(inter_pt[0] is not None)
        assert(currintert <= 1 + 1e-10)
        inter_pt = np.asarray(inter_pt)

        # Restore the original orientation of the line segment, before returning it
        pt1ret, pt2ret = pt1, inter_pt
        if inbox2:
            pt1ret, pt2ret = pt2ret, pt1ret

        return pt1ret, pt2ret

# --- Compute Delaunay triangulation --------------------------------------------------

def norm2(X):
    return np.sqrt(np.sum(X ** 2))

def normalized(X):
    return X / norm2(X)

def get_triangle_normal(A, B, C):
    return normalized(np.cross(A, B) + np.cross(B, C) + np.cross(C, A))


def get_power_circumcenter(A, B, C):
    N = get_triangle_normal(A, B, C)
    return (-.5 / N[2]) * N[:2]


def is_ccw_triangle_old(A, B, C):
    M = np.concatenate([np.stack([A, B, C]), np.ones((3, 1))], axis = 1)
    return np.linalg.det(M) > 0

def is_ccw_triangle(A,B,C):
    return (C[1]-A[1]) * (B[0]-A[0]) > (B[1]-A[1]) * (C[0]-A[0])

# Return true if line segments AB and CD intersect
def segments_intersect(l1,l2):
    A,B = l1
    C,D = l2
    return ccw(A,C,D) != ccw(B,C,D) and ccw(A,B,C) != ccw(A,B,D)


def get_power_triangulation(S, R):

    # Compute the lifted weighted points
    S_norm = np.sum(S ** 2, axis = 1) - R ** 2
    S_lifted = np.concatenate([S, S_norm[:,None]], axis = 1)

    # Special case for 3 points
    if S.shape[0] == 3:
        if is_ccw_triangle(S[0], S[1], S[2]):
            return [[0, 1, 2]], np.array([get_power_circumcenter(*S_lifted)])
        else:
            return [[0, 2, 1]], np.array([get_power_circumcenter(*S_lifted)])

    # Compute the convex hull of the lifted weighted points.
    # (~ 0.9 seconds on my machine for 100,000 points)
    hull = ConvexHull(S_lifted)

    # Extract the Delaunay triangulation from the lower hull.
    # This step is a little too slow for comfort. (~ 1.4 seconds on 100,000 points)
    tri_list = tuple([a, b, c] if is_ccw_triangle(S[a], S[b], S[c]) else [a, c, b]  for (a, b, c), eq in zip(hull.simplices, hull.equations) if eq[2] < 0) # Note that I set eq[2] < 0 instead of eq[2] <= 0, to avoid simplices that when projected to the plane are degenerate. This doesn't break anything because we are ignoring infinite rays in the power diagram.
    
    # Compute the Voronoi points. These are the power circumcenters of the triangles.
    # (~ 0.6 seconds on 100,000 points)
    ## Old code was:
    # other_V = [get_power_circumcenter(*S_lifted[tri]) for tri in tri_list]
    # print(np.all(np.abs(other_V - V) < 1e-10))

    num_tri = len(tri_list)
    A = np.zeros((num_tri, 3))
    B = np.zeros((num_tri, 3))
    C = np.zeros((num_tri, 3))

    for i, tri in enumerate(tri_list):
        A[i,:] = S_lifted[tri[0],:]
        B[i,:] = S_lifted[tri[1],:]
        C[i,:] = S_lifted[tri[2],:]
    N = np.cross(A, B) + np.cross(B, C) + np.cross(C, A)
    N = normalize(N)
    V = -N[:,:2] * 0.5

    tri_list = np.asarray(tri_list)

    V[:,0] = np.divide(V[:,0], N[:,2])
    V[:,1] = np.divide(V[:,1], N[:,2])

    return tri_list, V



# --- Compute Voronoi cells ---------------------------------------------------

'''
Compute arrangement of segments of Voronoi cells. Not necessary to store face information,
since the nearest-neighbor queries will be solved later with a KD-tree.

The segments are truncated to fit in the bounding box.
Infinite segments are ignored (unlike in the original code on which this is based)
'''
def get_power_diagram_truncated_from_power_triangulation(S, V, tri_list, bounding_box):

    segment_list = [] # (start_vtx, end_vtx)

    # Keep track of which edges separate which triangles
    edge_map = { }
    for i, tri in enumerate(tri_list):
        for edge in itertools.combinations(tri, 2):
            edge = tuple(sorted(edge))
            if edge in edge_map:
                edge_map[edge].append(i)
            else:
                edge_map[edge] = [i]

    for i, (a, b, c) in enumerate(tri_list):
        # For each edge of the triangle
        for u, v, w in ((a, b, c), (b, c, a), (c, a, b)):
            edge = tuple(sorted((u, v)))
            if len(edge_map[edge]) == 2:
            # Finite Voronoi edge
                j, k = edge_map[edge]
                if k == i:
                    j, k = k, j
                    continue
                if np.all(V[j] == V[k]):
                    continue

                pt1, pt2 = truncate_finite_segment_to_box(V[j], V[k], bounding_box)
                if not pt1 is None:
                    segment_list.append((pt1, pt2))
            else: 
                # Infinite Voronoi edge -- ignore
                pass

    return segment_list


def get_power_diagram_truncated(S, R, bounding_box):
    tri_list, V = get_power_triangulation(S, R)
    power_segments = get_power_diagram_truncated_from_power_triangulation(S, V, tri_list, bounding_box)
    return power_segments


##def display_triang(S, R, tri_list, V, voronoi_segments=None):
##    # Setup
##    fig, ax = plt.subplots()
##    plt.axis('equal')
##    plt.axis('off')    

##    # Set min/max display size, as Matplotlib does it wrong
##    min_corner = np.amin(S, axis = 0) - np.max(R)
##    max_corner = np.amax(S, axis = 0) + np.max(R)
##    plt.xlim((min_corner[0], max_corner[0]))
##    plt.ylim((min_corner[1], max_corner[1]))

##    # Plot the samples
##    for Si, Ri in zip(S, R):
##        Ri = Ri * 1e8 + 0.01
##        ax.add_artist(plt.Circle(Si, Ri, fill = True, alpha = .4, lw = 0., color = '#8080f0', zorder = 1))

##    # Plot the power triangulation
##    edge_set = frozenset(tuple(sorted(edge)) for tri in tri_list for edge in itertools.combinations(tri, 2))
##    line_list = LineCollection([(S[i], S[j]) for i, j in edge_set], lw = 1., colors = '.9')
##    line_list.set_zorder(0)
##    ax.add_collection(line_list)

##    if voronoi_segments is not None:
####        voronoi_segments = [x for x in voronoi_segments if np.any(x[0] != x[1])]
####        print(voronoi_segments)
####        for degseg in voronoi_segments:
####            ax.add_artist(plt.Circle(degseg[0], 0.001, fill=True))
##        seg_list = LineCollection(voronoi_segments)
##        ax.add_collection(seg_list)

##    for i, s in enumerate(S):
##        plt.text(s[0],s[1],str(i))

###    for v in V:
###        ax.add_artist(plt.Circle(v,0.1))

##    # Job done
##    plt.show()

