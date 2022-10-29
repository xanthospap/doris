#! /usr/bin/python

import datetime
import sys
import math
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.ticker import MultipleLocator

def t2sec(ts, t0):
    return [(t-t0).total_seconds() for t in ts]

def parse(fn):
    dct = {}
    with open(fn, 'r') as fin:
        for line in fin.readlines():
            l = line.split()
            try:
                t = datetime.datetime.strptime(' '.join([l[0],l[1][0:-3]]), '%Y-%m-%d %H:%M:%S.%f')
            except:
                t = datetime.datetime.strptime(' '.join([l[0],l[1]]), '%Y-%m-%d %H:%M:%S.%f')
            x,y,z,vx,vy,vz,mjd = [float(k) for k in l[2:]]
            dct[t] = [x,y,z,vx,vy,vz,mjd]
    return dct

def makeDiffs(dct1, dct2):
    dct = {}
    for t,vals in dct1.items():
        if t not in dct2:
            print('## Warning! failed to find record for date: {:}'.format(t))
        else:
            vals2 = dct2[t]
            if vals[-1] != vals2[-1]: print('Ops! difference is {:.15e}'.format(abs(vals[-1] - vals2[-1])))
            assert(abs(vals[-1] - vals2[-1]) < 1e-9)
            diffs = [ abs(x[0]-x[1]) for x in zip(vals,vals2)]
            dct[t] = diffs
    return dct

def colAsArray(dct,col,fac=1e0):
    return [ vals[col]*fac for _,vals in dct.items() ]

def whichCol(component, posvel):
    return component + int(posvel==1)*3

def whichTitle(component, posvel, fac=1e0):
    title = ' [{:}m]'.format('k' if fac==1e-3 else '')
    if posvel == 1:
        title = ' Velocity [{:}m/sec]'.format('k' if fac == 1e-3 else '')
    return ['X','Y','Z'][component] + title

def reduce(dates, dct):
    dctcp = {}
    for date in dates:
        if date in dct:
            dctcp[date] = dct[date]
        else:
            raise RuntimeError('ERROR WTF date')
    return dctcp

def ColStatistics(dct,col,fac=1e0):
    array = colAsArray(dct,col,fac)
    sum = 0e0
    sumSquared = 0e0
    max = -1e99
    min = 1e99
    for value in array:
        sum += value
        sumSquared += value *value
        if value > max: max = value
        if value < min: min = value
    n = float(len(array))
    return min, max, sum / n, math.sqrt(sumSquared / n)

dct1 = parse(sys.argv[1])
dct2 = parse(sys.argv[2])
diffs = makeDiffs(dct1,dct2) if len(dct2) > len(dct1) else makeDiffs(dct2,dct1)

t0 = datetime.datetime.max
for t in dct1:
  if t < t0: t0 = t

start_sec = 0

fig, axs = plt.subplots(3, 2, figsize=(10, 6), constrained_layout=True)
t = [ ti for ti in diffs ]
for component in range(3):
    for pv in range(2):
        fac = 1e-3 if pv == 0 else 1e0 ## noop .... just meters now
        fac = 1e0
        index = whichCol(component, pv)
        
        axs[component, pv].scatter(t2sec(t,t0),colAsArray(diffs,index,fac),s=1,color='black')
        
        axs[component, pv].set_title(whichTitle(component, pv, fac))
        _, end = axs[component, pv].get_xlim()
        axs[component, pv].xaxis.set_minor_locator(MultipleLocator(3600))
        axs[component, pv].xaxis.set_ticks(np.arange(start_sec, end, 3600*2))
        
        minv, maxv, average, rms = ColStatistics(diffs,index,fac)
        text_x = t[-1] - (t[-1] - t[0]) / 10e0
        text_y = maxv - (maxv-minv) / 10e0
        axs[component, pv].text(text_x, text_y, 'Average: {:+.4f}, Rms: {:.4f}'.format(average, rms), style='italic',
        bbox={'facecolor': 'red', 'alpha': 0.5, 'pad': 10})

        axs[component, pv].grid(True, 'both', 'x')

# plt.tight_layout()

if len(sys.argv) > 3:
    subtitle = ' '.join(sys.argv[3:])
else:
    subtitle = ''
fig.suptitle('Sp3 - Integrator Diffs\n{}'.format(subtitle), fontsize=16)

plt.savefig(sys.argv[2] + '.png', transparent=True)
plt.show()
