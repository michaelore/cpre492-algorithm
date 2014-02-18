#!/usr/bin/env python

import sys
from random import *

raw = sys.stdin.readlines()

def line_to_latlon(strng):
    return tuple(map(float, strng.split('\t')))

data = [line_to_latlon(row) for row in raw]

for (lat, lon) in data:
    print "%f\t%f" % (lat+gauss(0, 1e-1), lon+gauss(0, 1e-1))
