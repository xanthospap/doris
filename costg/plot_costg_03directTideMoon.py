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

# scale for x, y and z accelerations (m/s2)
scale = 1e15

## use TeX to render Matplotlib text
plt.rcParams['text.usetex'] = True

use_hours_of_day = True
if use_hours_of_day:
  try:
    import julian
  except:
    print('Failed to import julian module (see https://pypi.org/project/julian/)', file=sys.stderr)
    print('Either install it, or set the \"use_hours_of_day\" option to false', file=sys.stderr)
    sys.exit(1)

def component(key, lofd, fac=None): 
  return [ dct[key] for dct in lofd ] if not fac else [ dct[key]*fac for dct in lofd ]
def units(_): 
  return r"$m/s^2$"
def handle_date(mjd): return julian.from_jd(mjd, fmt='mjd') if use_hours_of_day else mjd
def tex_ytitle(cmp):
  if cmp[0] == 'x': return r"$\ddot{x}$ [" + units(cmp) + r"] $\times 10^{-15}$"
  elif cmp[0] == 'y': return r"$\ddot{y}$ [" + units(cmp) + r"] $\times 10^{-15}$"
  else: return r"$\ddot{z}$ [" + units(cmp) + r"] $\times 10^{-15}$"

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
            dct.append({'mjd': handle_date(float(l[0])), 'xacc':float(l[1]), 'yacc':float(l[2]), 
                'zacc':float(l[3])})
        except:
            print("Failed to parse line: [{:}]; skipped!".format(line.strip()), 
                file=sys.stderr)

## plot data/diffs
i = 0; j = 0
fig, axs = plt.subplots(3, 1, figsize=(10, 6), sharex=True, sharey=True, constrained_layout=True)
t = component('mjd', dct)
for y in ['xacc', 'yacc', 'zacc']:
    r,c = i%3, j
    i+=1
    x = component(y, dct, scale)
    sts = stats.describe(x)
    print('Component: {:} #pts {:} | Min {:+.2e} Max {:+.2e} Mean {:.2e} Std. Deviation {:.2e}'.format(y[0],sts.nobs, sts.minmax[0], sts.minmax[1], sts.mean, math.sqrt(sts.variance)))
    axs[r].grid(axis='y', color='0.85')
    axs[r].scatter(t,x,s=1,color='black')
    axs[r].set_ylabel(tex_ytitle(y), fontsize=14)
    # use formatters to specify major and minor ticks
    axs[r].xaxis.set_major_formatter(mdates.DateFormatter("%H:%M"))
    axs[r].xaxis.set_minor_locator(mdates.HourLocator())
axs[-1].set_xlabel("Date {:}".format(t[0].strftime("%Y/%m/%d")), fontsize=12)

## Rotate date labels automatically
fig.suptitle('COST-G Benchmark Diffs\nMoon Tide Acceleration Comparisson (CRF)', fontsize=16)
#fig.suptitle('Moon Tide Acceleration (CRF)', fontsize=16)
plt.savefig('03directMoonTide.png')
plt.show()
