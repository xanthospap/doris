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
def units(key):
    if key == 'mjd': return 'days'
    if key in ['xp', 'yp', 'X', 'Y', 'cio', 'tio']: return 'arcsec'
    return 'sec'

##
## Read file from STDIN, expect 
## Mjd    xp('')    yp('') dut1 (sec)  lod (sec)    X ('')    Y ('')  CIO ('')  TIO ('')
##
dct = []
for line in fileinput.input():
    ## parse file/input
    if line[0] != '#':
        l = line.split()
        try:
            dct.append({'mjd': float(l[0]), 'xp':float(l[1]), 'yp':float(l[2]), 
                'dut1':float(l[3]), 'lod': float(l[4]), 'X':float(l[5]), 
                'Y':float(l[6]), 'cio':float(l[7]), 'tio':float(l[8])})
        except:
            print("Failed to parse line: [{:}]; skipped!".format(line.strip()), 
                file=sys.stderr)

## plot data/diffs
i = 0; j = 0
fig, axs = plt.subplots(4, 2, figsize=(10, 6), constrained_layout=True)
t = component('mjd', dct)
for y in ['xp', 'yp', 'X', 'Y', 'dut1', 'lod', 'cio', 'tio']:
    r,c = i%4, j//4
    i+=1
    j+=1
    x = component(y, dct)
    axs[r,c].scatter(t,x,s=1,color='black')
    axs[r,c].set_title(y + ' in ' + units(y))

## Rotate date labels automatically
fig.suptitle('COST-G Benchmark Diffs\n', fontsize=16)
plt.savefig('01earthRotation_interpolatedEOP.png')
plt.show()
