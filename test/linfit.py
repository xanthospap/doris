#! /usr/bin/python

import numpy as np
import matplotlib.pyplot as plt
import sys
import os

if len(sys.argv) < 2:
    print('Usage linfit [FILE]', file=sys.stderr)
    sys.exit(1)

if not os.path.isfile(sys.argv[1]):
    print('ERROR File {:} does not exist!'.format(sys.argv[1]), file=sys.stderr)
    sys.exit(1)

xl = []
yl = []
with open(sys.argv[1], 'r') as fin:
    for line in fin.readlines():
        l = line.split()
        x, y = [float(x) for x in line.split()[0:2]]
        xl.append(x)
        yl.append(y)

coefs = np.polyfit(np.array(xl), np.array(yl), 1)
print('Model: y = {:.4f} + {:.4f} * x'.format(coefs[1], coefs[0]))

poly = np.poly1d(coefs)
new_x = np.linspace(xl[0], xl[-1])
new_y = poly(new_x)
plt.plot(xl, yl, "o", new_x, new_y)
plt.show()