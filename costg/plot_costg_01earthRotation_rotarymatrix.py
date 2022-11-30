#! /usr/bin/python

import datetime
import sys, os
import math
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import argparse
import requests
from scipy.interpolate import interp1d
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.ticker import MultipleLocator
import fileinput

def component(key, lofd): return [ dct[key] for dct in lofd ]
def units(_): return '[-]'

##
## Read file from STDIN, expect 
## Mjd  xx, xy, xz, yx, yy, yz, zx, zy, zz [-]
##
dct = []
for line in fileinput.input():
    ## parse file/input
    if line[0] != '#':
        l = line.split()
        try:
            dct.append({'mjd': float(l[0]), 'xx':float(l[1]), 'xy':float(l[2]), 
                'xz':float(l[3]), 'yx': float(l[4]), 'yy':float(l[5]), 
                'yz':float(l[6]), 'zx':float(l[7]), 'zy':float(l[8]), 'zz':float(l[9])})
        except:
            print("Failed to parse line: [{:}]; skipped!".format(line.strip()), 
                file=sys.stderr)

## plot data/diffs
i = 0; j = 0
fig, axs = plt.subplots(3, 3, figsize=(10, 6), constrained_layout=True)
t = component('mjd', dct)
for y in ['xx', 'xy', 'xz', 'yx', 'yy', 'yz', 'zx', 'zy', 'zz']:
    r,c = i%3, j//3
    i+=1
    j+=1
    x = component(y, dct)
    axs[r,c].scatter(t,x,s=1,color='black')
    axs[r,c].set_title(y + ' in ' + units(y))

## Rotate date labels automatically
fig.suptitle('COST-G Benchmark Diffs\n', fontsize=16)
plt.show()
plt.savefig('01earthRotation_rotaryMatrix.png')
