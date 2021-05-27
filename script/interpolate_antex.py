#! /usr/bin/python

import numpy as np
from scipy import interpolate
import matplotlib.pyplot as plt

## x-axis is zanith angle in degrees
x_start = 0
x_stop = 90e0
x_step = 5e0
x_range_inclusive = True
x_numpts = int(np.floor((x_stop - x_start) / x_step + 1))
x, s = np.linspace(x_start, x_stop, num=x_numpts, endpoint=True, retstep=True)
assert(int(s) == int(x_step))
