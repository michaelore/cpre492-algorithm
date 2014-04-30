#!/usr/bin/env python2
# CprE 491 May14-31 approximation prototype
# Author: Michael Ore

from __future__ import division
import scipy as sp
import matplotlib.pyplot as plt
from scipy import stats
from scipy import spatial
from scipy.spatial import kdtree
from scipy import random
import random
import sys

### Basic geometry functions

def cross((px, py), (qx, qy), (rx, ry)):
    return (rx-qx)*(py-qy)-(ry-qy)*(px-qx)

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
    return sp.arccos((ux*vx+uy*vy)/
                   sp.sqrt((ux*ux+uy*uy)*(vx*vx+vy*vy)))

def mag((x, y)):
    return x**2+y**2

def in_polygon(p, ps):
    n = ps.shape[0]
    summ = 0
    for i in xrange(n):
        if (cross(p, ps[i], ps[(i+1)%n]) < 0):
            summ -= angle(p, ps[i], ps[(i+1)%n])
        else:
            summ += angle(p, ps[i], ps[(i+1)%n])
    return abs(summ-2*sp.pi) < 1e-9 or abs(summ+2*sp.pi) < 1e-9

def bounding_box(ps):
    xmin = sys.float_info.max
    xmax = -sys.float_info.max
    ymin = sys.float_info.max
    ymax = -sys.float_info.max
    for (x, y) in ps:
        xmin = min(xmin, x)
        xmax = max(xmax, x)
        ymin = min(ymin, y)
        ymax = max(ymax, y)
    return ((xmin, xmax), (ymin, ymax))

def rot_mat(angle):
    return sp.array([[sp.cos(angle), -sp.sin(angle)], [sp.sin(angle), sp.cos(angle)]])

def rotate(ps, angle):
    rot = rot_mat(angle)
    return rot.dot(ps)

def circumcenter((ax, ay), (bx, by), (cx, cy)):
    d = 2*(ax*(by-cy)+bx*(cy-ay)+cx*(ay-by))
    ux = ((ax**2+ay**2)*(by-cy)+(bx**2+by**2)*(cy-ay)+(cx**2+cy**2)*(ay-by))/d
    uy = ((ax**2+ay**2)*(cx-bx)+(bx**2+by**2)*(ax-cx)+(cx**2+cy**2)*(bx-ax))/d
    return (ux, uy)

def in_circle((px, py), (ax, ay), (bx, by), (cx, cy)):
    eps = 1e-9
    (centx, centy) = circumcenter((ax, ay), (bx, by), (cx, cy))
    d = mag((ax-centx, ay-centy))
    return d-eps > mag((px-centx, py-centy))

# A crude polygon for Florida's border
florida_border = sp.array([[30.95, -87.6],
                           [30.98, -85],
                           [30.72, -84.9],
                           [30.53, -81.45],
                           [25.27, -80.02],
                           [25.28, -80.98],
                           [27.85, -82.83],
                           [28.9, -82.73],
                           [29.73, -83.55],
                           [29.7, -85.27],
                           [30.32, -87.32]])

f = open(sys.argv[1], "r")

for line in f:
    parsed = line.split('\t')
    lat = float(parsed[0])
    lon = float(parsed[1])
    if in_polygon((lat, lon), florida_border):
        print(line.strip())
