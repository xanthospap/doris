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
from scipy import stats

def component(key, lofd): return [ dct[key] for dct in lofd ]
def units(_): return '[m/sec^2]'

##
## Read file from STDIN, expect 
## Mjd xacc yacc zacc [m/sec^2]
##
dct = []
for line in fileinput.input():
    ## parse file/input
    if line[0] != '#':
        l = line.split()
        try:
            dct.append({'mjd': float(l[0]), 'xacc':float(l[1]), 'yacc':float(l[2]), 
                'zacc':float(l[3])})
        except:
            print("Failed to parse line: [{:}]; skipped!".format(line.strip()), 
                file=sys.stderr)

## plot data/diffs
i = 0; j = 0
fig, axs = plt.subplots(3, 1, figsize=(10, 6), constrained_layout=True)
t = component('mjd', dct)
for y in ['xacc', 'yacc', 'zacc']:
    r,c = i%3, j
    i+=1
    x = component(y, dct)
    sts = stats.describe(x)
    _stats = ' {:+.1e} +- {:.1e}, |max|: {:.1e}'.format(sts.mean, math.sqrt(sts.variance), max(abs(sts.minmax[0]), abs(sts.minmax[1])))
    print(_stats)
    axs[r].scatter(t,x,s=1,color='black',label=_stats)
    axs[r].set_title(y + ' in ' + units(y) + _stats)

## Rotate date labels automatically
fig.suptitle('COST-G Benchmark Diffs\n', fontsize=16)
plt.savefig('02gravityfield_itrf.png')
plt.show()
