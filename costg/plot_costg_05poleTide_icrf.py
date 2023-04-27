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
scale = 1e9
def unitsstr(): return r'$[m/s^2]\times{:.1e}$'.format(1e0/scale) if scale != 1e0 else r'$[m/s^2]$'
use_hours_of_day = True

## use TeX to render Matplotlib text
plt.style.use('seaborn-v0_8-pastel')
tex_fonts = {
    # Use LaTeX to write all text
    "text.usetex": True,
    "font.family": "serif",
    # Use 10pt font in plots, to match 10pt font in document
    "axes.labelsize": 10,
    "font.size": 10,
    # Make the legend/label fonts a little smaller
    "legend.fontsize": 7,
    "xtick.labelsize": 9,
    "ytick.labelsize": 9
}
plt.rcParams.update(tex_fonts)

if use_hours_of_day:
  try:
    import julian
  except:
    print('Failed to import julian module (see https://pypi.org/project/julian/)', file=sys.stderr)
    print('Either install it, or set the \"use_hours_of_day\" option to false', file=sys.stderr)
    sys.exit(1)

def handle_date(mjd): return julian.from_jd(mjd, fmt='mjd') if use_hours_of_day else mjd
def get_diffs(dct, component):
    t = []; y= [];
    k1,k2 = list(zip(['rax','ray','raz'],['ax','ay','az']))[component]
    for entry in dct:
        t.append(entry['mjd'])
        y.append(entry[k1] - entry[k2])
    return t,y
def get_diffs_norm(dct):
    t = []; y= [];
    for entry in dct:
        t.append(entry['mjd'])
        n = 0e0
        for component in range(3):
            k1,k2 = list(zip(['rax','ray','raz'],['ax','ay','az']))[component]
            n += (entry[k1] - entry[k2])*(entry[k1] - entry[k2])
        y.append(math.sqrt(n))
    return t,y

dct = []
for line in fileinput.input():
    ## parse file/input
    if line[0] != '#':
        l = line.split()
        try:
            dct.append({'mjd': handle_date(float(l[0])), 'rax':float(l[1]), 'ray':float(l[2]), 
                'raz':float(l[3]), 'ax':float(l[4]), 'ay':float(l[5]), 'az':float(l[6])})
        except:
            print("Failed to parse line: [{:}]; skipped!".format(line.strip()), 
                file=sys.stderr)

## plot data/diffs
fig, axs = plt.subplots(4, 1, figsize=(10, 6), sharex=True, sharey=True, constrained_layout=True)
for i in range(3):
    t,y = get_diffs(dct, i)
    y = [yy * scale for yy in y]
    sts = stats.describe(y)
    axs[i].scatter(t,y,s=1,color='black')
    axs[i].axhline(sts.mean, color='r', linestyle='--', linewidth=1)
    axs[i].text(t[0],sts.minmax[0],r'{:+.1e} $\pm$ {:+.3e}'.format(sts.mean, math.sqrt(sts.variance)),color='red')
    axs[i].text(t[0],sts.minmax[1],r'min:{:+.1e} max:{:+.1e}'.format(sts.minmax[0], sts.minmax[1]),color='red')
    axs[i].set_ylabel([r'$\delta\ddot{x}$', r'$\delta\ddot{y}$', r'$\delta\ddot{z}$'][i]+unitsstr())
    axs[i].xaxis.set_major_formatter(mdates.DateFormatter("%H:%M"))
    axs[i].xaxis.set_minor_locator(mdates.HourLocator())

# plot norms
t,y = get_diffs_norm(dct)
y = [yy * scale for yy in y]
sts = stats.describe(y)
axs[3].scatter(t,y,s=1,color='black')
axs[3].axhline(sts.mean, color='r', linestyle='--', linewidth=1)
axs[3].text(t[0],-sts.minmax[1],r'{:+.1e} $\pm$ {:+.3e}'.format(sts.mean, math.sqrt(sts.variance)),color='red')
axs[3].text(t[0],sts.minmax[1],r'min:{:+.1e} max:{:+.1e}'.format(sts.minmax[0], sts.minmax[1]),color='red')
axs[3].set_ylabel(r'$\delta\ddot{a}$'+unitsstr())
axs[3].xaxis.set_major_formatter(mdates.DateFormatter("%H:%M"))
axs[3].xaxis.set_minor_locator(mdates.HourLocator())

axs[-1].set_xlabel("Date {:}".format(t[0].strftime("%Y/%m/%d")), fontsize=12)

## Rotate date labels automatically
fig.suptitle('COST-G Benchmark Diffs\nSolid Earth Pole Tide Acceleration Comparisson', fontsize=16)
plt.savefig('05poleTide_icrf.pdf')
plt.show()
