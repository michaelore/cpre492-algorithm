#!/usr/bin/env python

from __future__ import division
import sys
from math import *
from itertools import *
import matplotlib.pyplot as plt

import geometry as geo
from geometry import EPS

import localmaxchecker as lmc

raw = sys.stdin.readlines()

nrows = len(raw)

with open(sys.argv[1]) as f:
    raw_border = f.readlines()

def line_to_latlon(strng):
    return tuple(map(float, strng.split('\t')))

def transpose((x, y)):
    return (y, x)

data = [transpose(line_to_latlon(row)) for row in raw]

border = [transpose(line_to_latlon(row)) for row in raw_border]

nborder = len(border)

t3_circle_groups = [[] for i in xrange(nrows)]

for icircle in combinations(xrange(nrows), 3):
    (p1, p2, p3) = (data[icircle[0]], data[icircle[1]], data[icircle[2]])
    center = geo.circumcenter(p1, p2, p3)
    rad = geo.dist(center, data[icircle[0]])
    if not geo.in_polygon(center, border):
        continue
    count = 0
    intervals = []
    for ip in xrange(nrows):
        vec = geo.trans_vec(data[ip], center)
        d = geo.mag(vec)
        if d < rad-EPS:
            count += 1
        elif d < rad+EPS:
            theta = geo.theta(vec)
            intervals.append((theta-pi/2, theta+pi/2, 1))
    for i in xrange(lmc.removability(intervals)):
        t3_circle_groups[count+i].append((center, rad))

print "Category 3 Circles:"
for i in xrange(len(t3_circle_groups)):
    print i, len(t3_circle_groups[i])

t2_circle_groups = [[] for i in xrange(nrows)]

for iline in combinations(xrange(nrows), 2):
    bisec = geo.bisector(data[iline[0]], data[iline[1]])
    for border_line in geo.poly_lines(border):
        if geo.are_parallel(bisec[1], border_line[1]):
            continue
        (t, u) = geo.intersect(bisec, border_line)
        center = geo.translate(border_line[0], geo.scale(u, border_line[1]))
        rad = geo.dist(center, data[iline[0]])
        if u <= 0 or u >= 1:
            continue
        count = 0
        intervals = []
        for ip in xrange(nrows):
            vec = geo.trans_vec(data[ip], center)
            d = geo.mag(vec)
            if d < rad-EPS:
                count += 1
            elif d < rad+EPS:
                theta = geo.theta(vec)
                intervals.append((theta-pi/2, theta+pi/2, 1))
        limit = len(intervals)
        limstart = geo.theta(border_line[1])
        limend = limstart+pi
        intervals.append((limstart, limend, limit+1))
        for i in xrange(lmc.removability(intervals)):
            t2_circle_groups[count+i].append((center, rad))

print "Category 2 Circles:"
for i in xrange(len(t2_circle_groups)):
    print i, len(t2_circle_groups[i])

t1_circle_groups = [[] for i in xrange(nrows)]

for ipoint in xrange(nrows):
    for iborder in xrange(nborder):
        center = border[iborder]
        rad = geo.dist(center, data[ipoint])
        count = 0
        intervals = []
        for ip in xrange(nrows):
            vec = geo.trans_vec(data[ip], center)
            d = geo.mag(vec)
            if d < rad-EPS:
                count += 1
            elif d < rad+EPS:
                theta = geo.theta(vec)
                intervals.append((theta-pi/2, theta+pi/2, 1))
        backwardvec = geo.trans_vec(border[(iborder-1)%nborder], center)
        forwardvec = geo.trans_vec(border[(iborder+1)%nborder], center)
        limit = len(intervals)
        limstart = geo.theta(forwardvec)
        limend = geo.theta(backwardvec)
        intervals.append((limstart, limend, limit+1))
        for i in xrange(lmc.removability(intervals)):
            t1_circle_groups[count+i].append((center, rad))

print "Category 1 Circles:"
for i in xrange(len(t1_circle_groups)):
    print i, len(t1_circle_groups[i])

k = int(floor(0.01*nrows))

def pltCircle((center, rad)):
    return plt.Circle(center,rad,color='b',fill=False)

plt.figure(figsize=(8,12))

plt.scatter(map(lambda (x, y): x, data), map(lambda (x, y): y, data), facecolors='none', edgecolors='r')

plt.plot(map(lambda (x, y): x, border)+[border[0][0]], map(lambda (x, y): y, border)+[border[0][1]], 'k--')

fig = plt.gcf()
for circ in t3_circle_groups[k] + t2_circle_groups[k] + t1_circle_groups[k]:
    fig.gca().add_artist(pltCircle(circ))

plt.show()
