#!/usr/bin/env python2

import sys
from random import *

N = int(sys.argv[1])

for i in range(N):
    lat = (random()-0.5)*180
    lon = (random()-0.5)*360
    print("%s\t%s" % (lat, lon))
