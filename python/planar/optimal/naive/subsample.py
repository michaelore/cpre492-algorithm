#!/usr/bin/env python

import sys
import random

nitems = int(sys.argv[1])
seed = sys.argv[2]
random.seed(seed)
data = sys.stdin.readlines()
for line in random.sample(data, nitems):
    sys.stdout.write(line)
