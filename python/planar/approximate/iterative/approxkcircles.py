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

### Optimization functions
# Conceptually, we are optimizing a function that takes an (x, y) coordinate
# as input and outputs the distance from that point to the kth-closest data
# point; the "kth-closest distance function". Maximizing this function is the
# same as finding maximal k-circles.
#
# The algorithm is inspired by gradient descent and simulated annealing.
# Ordinary gradient descent won't converge because the gradient is non-zero at
# the extrema and behavior at boundaries is undefined. Here, the gradient is
# used to determine direction only and the magnitude steadily decreases. If the
# boundary is hit, the direction becomes semi-random.

# kth finds the distance to the kth-closest point, and its index in the data
def kth_closest(kd, k, p):
    (d, i) = kd.query(p, k+1)
    #print "%20d, %s, %10d" % (hash(str(i[:-1])), str(i[-8:-1]), i[-1])
    #print "%20d, %10d, %20d" % (hash(str(i[:-1])), i[-1], hash(str(sorted(i[:-1]))))
    print "%20d, %10d, %10d" % (hash(str(sorted(i[:-2]))), i[-2], i[-1])
    kthd = d[-2]
    kthi = i[-2]
    return (kthd, kthi)

# This is the gradient of the kth-closest distance function; a unit vector
# pointing away from the kth-closest point.
def gradient(kd, k, p):
    (kthd, kthi) = kth_closest(kd, k, p)
    kth = kd.data[kthi]
    return (1/kthd)*(p-kth)

def insert(kd, graph, i):
    r = str(sorted(i[:-3]))
    newcon = False
    if (r, i[-3]) not in graph:
        graph[(r, i[-3])] = set([(i[-2], i[-1])])
        #print "new connection"
        newcon = True
    else:
        s = graph[(r, i[-3])]
        if (i[-2], i[-1]) not in s:
            s.add((i[-2], i[-1]))
            #print "new connection"
            newcon = True
    if newcon:
        if (r, i[-2]) in graph:
            #print "i[-2] is in the graph"
            s = graph[(r, i[-2])]
            if (i[-1], i[-3]) in s or (i[-3], i[-1]) in s:
                #print "i[-2] matches"
                if (r, i[-1]) in graph:
                    #print "i[-1] is in the graph"
                    s = graph[(r, i[-1])]
                    if (i[-2], i[-3]) in s or (i[-3], i[-2]) in s:
                        #print "i[-1] matches"
                        (k3, k2, k1) = (kd.data[i[-3]], kd.data[i[-2]], kd.data[i[-1]])
                        center = circumcenter(k3, k2, k1)
                        if in_polygon(center, sp.array((k3, k2, k1))):
                            print "inside"
                            return (k3, k2, k1)
                        else:
                            print "not inside"
    return False

# One iteration of "gradient ascent". alpha is a scalar for the gradient. If
# the result would be beyond the border polygon, alpha is reduced and a new
# direction is selected randomly.
def gradient_ascent_once(kd, border, k, p, alpha, graph):
    (d, i) = kd.query(p, k+2)
    #print "%20d, %10d, %10d, %10d" % (hash(str(sorted(i[:-3]))), i[-3], i[-2], i[-1])
    kthd = d[-3]
    kthi = i[-3]
    kth = kd.data[kthi]
    grad = (1/kthd)*(p-kth)
    dispersion = random.gauss(0, 0.1)
    grad = rotate(grad, dispersion)
    result = p+alpha*grad
    maxi = 100
    while (not in_polygon(result, border)) and maxi > 0:
        print "bump"
        #return ((0, 0), True)
        alpha *= 0.8
        angle = sp.random.uniform(0, 2*sp.pi)
        grad = rotate(grad, angle)
        result = p+alpha*grad
        maxi -= 1
    return (result, insert(kd, graph, i))

# Gradient descent for "its" iterations. The "alpha"s for each iteration
# decrease as with the harmonic series. This should result in convergence only
# when zeroing in on a local maximum. Returns the results of all iterations.
def gradient_ascent(kd, border, k, p, alpha, its):
    graph = dict()
    early = False
    for i in xrange(1, its):
        (p, early) = gradient_ascent_once(kd, border, k, p, alpha*(0.8)**i, graph)
        if early:
            print "early", i
            return (p, early)
    print "not early"
    return ((0, 0), False)


# Generates n random points inside the border
def gen_points(border, n):
    ((xmin, xmax), (ymin, ymax)) = bounding_box(border)
    print ((xmin, xmax), (ymin, ymax))
    points = sp.empty((n, 2))
    for i in xrange(n):
        p = sp.array((sp.random.uniform(xmin, xmax), sp.random.uniform(ymin, ymax)))
        while (not in_polygon(p, border)):
            p = sp.array((sp.random.uniform(xmin, xmax), sp.random.uniform(ymin, ymax)))
    	points[i] = p
    return points

# Apply gradient_ascent to n random points
def gen_circles(kd, border, k, n, alpha, its):
    appl = lambda p: gradient_ascent(kd, border, k, p, alpha, its)
    #return sp.apply_along_axis(appl, 1, gen_points(border, n))
    return [appl(p) for p in gen_points(border, n)]

### Data input and processing

# Radius of Earth, in km
R = 6371

# Functions for converting coordinates from degrees to km, under Gall-Peters
# projection. Not strictly necessary, but I think it makes sense to use this
# or another equal-area projection given that we are judging circles by their
# area.
#gallf_x = sp.vectorize(lambda lam: R*sp.pi*lam/180.0/sp.sqrt(2))
#gallf_y = sp.vectorize(lambda phi: R*sp.sqrt(2)*sp.sin(sp.pi*phi/180))
gallf_x = sp.vectorize(lambda lam: lam)
gallf_y = sp.vectorize(lambda phi: phi)

# Reading data as tab-seperated-values
fldata = sp.genfromtxt('FL.uniquesites.txt', delimiter='\t')

# Apply Gall-Peters projection to input data
fldata_x = gallf_x(fldata[:,1])
fldata_y = gallf_y(fldata[:,0])
fldata_gall = sp.transpose(sp.vstack((fldata_x, fldata_y)))

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

# Apply Gall-Peters projection to the border
florida_border_x = gallf_x(florida_border[:,1])
florida_border_y = gallf_y(florida_border[:,0])
florida_border_gall = sp.transpose(sp.vstack((florida_border_x, florida_border_y)))

# The file had some data from all over the US, so here we restrict the data to
# a rough box around Florida.
((xmin, xmax), (ymin, ymax)) = bounding_box(florida_border_gall)
xmin -= 1
xmax += 1
ymin -= 1
ymax += 1
florida_gall = sp.array(filter(lambda (x, y): x >= xmin and x <= xmax and y >= ymin and y <= ymax, fldata_gall))

# Convex Hull as a default location constraint
#hull = spatial.ConvexHull(florida_gall)

# KDTree for k-nearest-neighbor queries
kd = kdtree.KDTree(florida_gall)

# Pick k to be 0.1% of the size of the data-set
k = int(0.001*florida_gall.shape[0])
#k = 5
results = gen_circles(kd, florida_border_gall, k, 100, 1, 50)

### Display code

# Plot a circle
def plt_circle(kd, k, p):
    (kthd, kthi) = kth_closest(kd, k, p)
    return plt.Circle(p,kthd,color='b',fill=False)

# Set size of figure
plt.figure(figsize=(8,12))

# Plot original data
plt.scatter(florida_gall[:,0], florida_gall[:,1], facecolors='none', edgecolors='r')

# Plot border
plt.plot(florida_border_x, florida_border_y, 'k--')

# Plot circles
fig = plt.gcf()
for (junk, circ) in results:
    if circ and circ[0][0] != 0:
        center = circumcenter(circ[0], circ[1], circ[2])
        plt.plot(circ[0][0], circ[0][1], 'o')
        plt.plot(circ[1][0], circ[1][1], 'o')
        plt.plot(circ[2][0], circ[2][1], 'o')
        plt.plot(center[0], center[1], 'o')
        fig.gca().add_artist(plt_circle(kd, k, center))

# We're done! Display the plot
plt.show()

"""
n = florida_gall.shape[0]
nsamples = 100
sample_count = sp.empty((nsamples))
sample_radius = sp.empty((nsamples))
cursor = 0
while cursor < nsamples:
    (i, j, k) = (sp.random.randint(0, n), sp.random.randint(0, n), sp.random.randint(0, n))
    if (i == j or i == k or j == k):
        continue
    (pi, pj, pk) = (florida_gall[i], florida_gall[j], florida_gall[k])
    cent = circumcenter(pi, pj, pk)
    if not in_polygon(cent, sp.array((pi, pj, pk))):
        continue
    if not in_polygon(cent, florida_border_gall):
        continue
    count = 0
    for t in xrange(n):
        if (in_circle(florida_gall[t], pi, pj, pk)):
            count += 1
        if count > 1000:
            break
    if count > 1000:
        continue
    sample_count[cursor] = count
    sample_radius[cursor] = mag((pi[0]-cent[0], pi[1]-cent[1]))
    print cursor
    cursor += 1

plt.figure(figsize=(12,12))
plt.scatter(sample_count, sample_radius)
plt.show()
"""
