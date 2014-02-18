
from math import *
from itertools import *

EPS = 1e-9

def mag((x, y)):
    return sqrt(x**2+y**2)

def cross((px, py), (qx, qy), (rx, ry)):
    return (rx-qx)*(py-qy)-(ry-qy)*(px-qx)

# returns true if r is on the left side of the line pq
def ccw(p, q, r):
    return cross(p, q, r) > 0

def angle((ax, ay), (bx, by), (cx, cy)):
    def clamp(x):
        if x > 1.0:
            return 1.0
        elif x < -1.0:
            return -1.0
        else:
            return x
    ux=bx-ax
    uy=by-ay
    vx=cx-ax
    vy=cy-ay
    return acos(clamp((ux*vx+uy*vy)/
                   sqrt((ux*ux+uy*uy)*(vx*vx+vy*vy))))

def theta((ax, ay)):
    return atan2(ay, ax)

def in_polygon(p, ps):
    n = len(ps)
    summ = 0
    for i in xrange(n):
        if (cross(p, ps[i], ps[(i+1)%n]) < 0):
            summ -= angle(p, ps[i], ps[(i+1)%n])
        else:
            summ += angle(p, ps[i], ps[(i+1)%n])
    return abs(summ-2*pi) < EPS or abs(summ+2*pi) < EPS

def circumcenter((ax, ay), (bx, by), (cx, cy)):
    d = 2*(ax*(by-cy)+bx*(cy-ay)+cx*(ay-by))
    ux = ((ax**2+ay**2)*(by-cy)+(bx**2+by**2)*(cy-ay)+(cx**2+cy**2)*(ay-by))/d
    uy = ((ax**2+ay**2)*(cx-bx)+(bx**2+by**2)*(ax-cx)+(cx**2+cy**2)*(bx-ax))/d
    return (ux, uy)

def trans_vec((ax, ay), (bx, by)):
    return (ax-bx, ay-by)

def dist(a, b):
    return mag(trans_vec(a, b))

def translate((ax, ay), (bx, by)):
    return (ax+bx, ay+by)

def scale(alpha, (ax, ay)):
    return (alpha*ax, alpha*ay)

def unit(a):
    return scale(1/mag(a), a)

def midpoint(a, b):
    return translate(scale(0.5, trans_vec(a, b)), b)

def in_circle(p, cent, rad):
    dist = mag(trans_vec(p, cent))
    return rad-EPS > dist

def in_triangle(a, b, c):
    return in_polygon(circumcenter(a, b, c), (a, b, c))

def poly_vectors(ps):
    n = len(ps)
    for i in xrange(n):
        yield trans_vec(ps[(i+1)%n], ps[i])

def poly_lines(ps):
    return izip(ps, poly_vectors(ps))

def are_parallel(a, b):
    return abs(angle(a, (0, 0), b)) < EPS

def mat_inv((a, b, c, d)):
    det = a*d-b*c
    return (d/det, -b/det, -c/det, a/det)

def mat_vec_mult((a, b, c, d), (e, f)):
    return (a*e+b*f, c*e+d*f)

def rotate(theta, a):
    return mat_vec_mult((cos(theta), -sin(theta), sin(theta), cos(theta)), a)

def bisector(a, b):
    return (midpoint(a, b), unit(rotate(pi/2, trans_vec(a, b))))

# p+t*r = q+u*s
def intersect(((px, py), (rx, ry)), ((qx, qy), (sx, sy))):
    (t, u) = mat_vec_mult(mat_inv((rx, -sx, ry, -sy)), (qx-px, qy-py))
    return (t, u)
