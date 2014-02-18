
import itertools
from math import *

TAU = 2*pi

class AngleWell:
    def __init__(self, intervals=[]):
        self.angles = []
        self.base = 0
        self.add_intervals(intervals)

    def add_intervals(self, intervals):
        for (start, end, weight) in intervals:
            start = start % TAU
            end = end % TAU
            self.angles += [(start, 0, weight), (end, weight, 0)]
            if end < start:
                self.base += weight
        newangles = []
        for (angle, angles) in itertools.groupby(sorted(self.angles), lambda(a, l, r): a):
            leftsum = 0
            rightsum = 0
            for (a, l, r) in angles:
                leftsum += l
                rightsum += r
            newangles.append((angle, leftsum, rightsum))
        self.angles = newangles

    def removability(self):
        minval = self.base
        cumsum = self.base
        for (angle, left, right) in self.angles:
            cumsum -= left
            minval = min(minval, cumsum)
            cumsum += right
        return minval

def removability(intervals):
    return AngleWell(intervals).removability()
