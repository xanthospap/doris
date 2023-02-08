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

# scale for x, y and z acclerations (m/s2)
scale = 1e0

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

def units(cmp):
  return r"mas" if cmp in ["xp", "yp"] else r"msec"
  #return r"$m/s^2$"

def handle_date(mjd): return julian.from_jd(mjd, fmt='mjd') if use_hours_of_day else mjd

def tex_ytitle(cmp):
  if cmp == 'xp': return r"$x_p$ " + units(cmp)
  elif cmp == 'yp': return r"$y_p$ " + units(cmp)
  elif cmp == 'dut1': return r"$\delta UT1$ " + units(cmp)
  elif cmp == 'lod': return r"$\delta LOD$ " + units(cmp)
  else: raise RuntimeError()

##
## Read file from STDIN, expect 
## [INTRP || RDATA] mjd 59218.996874999997 xp +0.06209981 yp +0.30840237 dut1 -0.17436115 lod +0.00015511 dx +0.00016087 dy +0.00006598
##
intrp_dct = []
rdata_dct = []
for line in fileinput.input():
    ## parse file/input
    if line[0] != '#' :
        l = line.split()
        s = l[0]+' '.join([l[i] for i in range(1,len(l)) if i%2!=0])
        if s in ["[INTRP]mjd xp yp dut1 lod dx dy", "[RDATA]mjd xp yp dut1 lod dx dy"]:
            try:
                if l[0] == "[INTRP]":
                    intrp_dct.append({'mjd': handle_date(float(l[2])), 'xp':float(l[4]), 'yp':float(l[6]), 
                        'dut1':float(l[8]), 'lod':float(l[10]), 'dX':float(l[12]),'dY':float(l[14]),})
                else:
                    rdata_dct.append({'mjd': handle_date(float(l[2])), 'xp':float(l[4]), 'yp':float(l[6]), 
                        'dut1':float(l[8]), 'lod':float(l[10]), 'dX':float(l[12]),'dY':float(l[14]),})
            except:
                print("Failed to parse line: [{:}]; skipped!".format(line.strip()), 
                    file=sys.stderr)

## plot data/diffs
i = 0; j = 0
fig, axs = plt.subplots(4, 1, figsize=(10, 6), sharex=True, sharey=False, constrained_layout=True)
for y in ['xp', 'yp', 'dut1', 'lod']:
    r,c = i%4, j
    i+=1
    scale = 1e3
    axs[r].scatter(component('mjd', intrp_dct),component(y, intrp_dct, scale),s=1,color='black')
    axs[r].scatter(component('mjd', rdata_dct),component(y, rdata_dct, scale),s=30,color='red',marker='x')
    axs[r].set_ylabel(tex_ytitle(y), fontsize=14)
    # use formatters to specify major and minor ticks
    axs[r].xaxis.set_major_formatter(mdates.DateFormatter("%d/%m %H:%M"))
    #axs[r].xaxis.set_minor_locator(mdates.HourLocator())
axs[-1].set_xlabel("Date (UTC), year = 2021", fontsize=12)

## Rotate date labels automatically
fig.suptitle('Earth Orientation Parameters Interpolation', fontsize=16)
plt.savefig('eop_interpolation.png')
plt.show()
