#! /usr/bin/python

import datetime
import sys
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages


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

def colAsArray(dct,col):
    return [ vals[col] for _,vals in dct.items() ]

def whichCol(component, posvel):
    return component + int(posvel==1)*3

def whichTitle(component, posvel):
    title = ' [m]'
    if posvel == 1: title = ' Velocity [m/sec]'
    return ['X','Y','Z'][component] + title

def reduce(dates, dct):
    dctcp = {}
    for date in dates:
        if date in dct:
            dctcp[date] = dct[date]
        else:
            raise RuntimeError('ERROR WTF date')
    return dctcp

def dump2Table(dct1, dct2):
  print('{:22s} {:15s} {:15s} {:15s} {:15s} {:15s} {:15s} {:15s} {:15s} {:15s} {:15s} {:15s} {:15s} {:15s} {:15s} {:15s} {:15s} {:15s} {:15s} {:15s} {:15s}'.format('Epoch', 'X_sp3', 'Y_sp3', 'Z_sp3', 'Vx_sp3', 'Vy_sp3', 'Vz_sp3', 'X_integrator', 'Y_integrator', 'Z_integrator', 'Vx_integrator', 'Vy_integrator', 'Vz_integrator', 'Mjd', 'SecInDay', 'DX', 'DY', 'DZ', 'DVx', 'DVy', 'DVz'))

  t0 = datetime.datetime.max
  for t in dct1:
    if t < t0: t0 = t

  for t,vals in dct1.items():
    if t in dct2:
      vals2 = dct2[t]
      if vals[-1] == vals2[-1]:
        #assert(vals[-1] == vals2[-1])
        print('{:} {:+15.3f} {:+15.3f} {:+15.3f} {:+15.3f} {:+15.3f} {:+15.3f} {:+15.3f} {:+15.3f} {:+15.3f} {:+15.3f} {:+15.3f} {:+15.3f} {:+15.6f} {:+15.1f} {:+15.1f} {:+15.1f} {:+15.1f} {:+15.1f} {:+15.1f} {:+15.1f}'.format(t.strftime('%Y-%m-%dT%H:%M:%S'), *vals[:-1], *vals2, (t-t0).total_seconds(), *[x[0]-x[1] for x in zip(vals[0:6],vals2[0:6])]))

dct1 = parse(sys.argv[1])
dct2 = parse(sys.argv[2])
#dump2Table(dct1, dct2)

diffs = makeDiffs(dct1,dct2) if len(dct2) > len(dct1) else makeDiffs(dct2,dct1)

t0 = datetime.datetime.max
for t in dct1:
  if t < t0: t0 = t

book = sys.argv[2] + '.pdf'
start_sec = 0

with PdfPages(book) as pdf:

    fig, axs = plt.subplots(3, 2, figsize=(10, 6), constrained_layout=True)
    t = [ ti for ti in diffs ]
    for component in range(3):
      for pv in range(2):
          index = whichCol(component, pv)
          axs[component, pv].scatter(t2sec(t,t0),colAsArray(diffs,index),s=1,color='black')
          axs[component, pv].set_title(whichTitle(component, pv))
          _, end = axs[component, pv].get_xlim()
          axs[component, pv].xaxis.set_ticks(np.arange(start_sec, end, 3600))

    plt.tight_layout()
    fig.suptitle('Sp3 - Integrator Diffs')
    #plt.savefig('diffs.png', bbox_inches='tight')
    pdf.savefig()
    plt.close()

    fig, axs = plt.subplots(3, 2, figsize=(10, 6), constrained_layout=True)
    t = [ ti for ti in diffs ]
    dct1 = reduce(t, dct1)
    for component in range(3):
      for pv in range(2):
          index = whichCol(component, pv)
          axs[component, pv].scatter(t2sec(t,t0),colAsArray(dct1,index), s=1,color='black')
          axs[component, pv].set_title(whichTitle(component, pv))
          _, end = axs[component, pv].get_xlim()
          axs[component, pv].xaxis.set_ticks(np.arange(start_sec, end, 3600))
    plt.tight_layout()
    fig.suptitle('Sp3 State')
    #plt.savefig('dct1.png', bbox_inches='tight')
    pdf.savefig()
    plt.close()

    fig, axs = plt.subplots(3, 2, figsize=(10, 6), constrained_layout=True)
    t = [ ti for ti in diffs ]
    dct2 = reduce(t, dct2)
    for component in range(3):
      for pv in range(2):
          index = whichCol(component, pv)
          axs[component, pv].scatter(t2sec(t,t0),colAsArray(dct2,index), s=1,color='black')
          axs[component, pv].set_title(whichTitle(component, pv))
          _, end = axs[component, pv].get_xlim()
          axs[component, pv].xaxis.set_ticks(np.arange(start_sec, end, 3600))
    plt.tight_layout()
    fig.suptitle('Integrator State')
    #plt.savefig('dct2.png', bbox_inches='tight')
    pdf.savefig()
    plt.close()

    fig, axs = plt.subplots(3, 2, figsize=(10, 6), constrained_layout=True)
    t = [ ti for ti in diffs ]
    for component in range(3):
      for pv in range(2):
          index = whichCol(component, pv)
          axs[component, pv].scatter(t2sec(t,t0),colAsArray(dct1,index),s=1.2,color='black')
          axs[component, pv].scatter(t2sec(t,t0),colAsArray(dct2,index),s=1,color='red')
          axs[component, pv].set_title(whichTitle(component, pv))
          _, end = axs[component, pv].get_xlim()
          axs[component, pv].xaxis.set_ticks(np.arange(start_sec, end, 3600))
    plt.tight_layout()
    fig.suptitle('Sp2 and Integrator State Overlay')
    # plt.savefig('foo3.png', bbox_inches='tight')
    pdf.savefig()
    plt.close()

    d = pdf.infodict()
    d['Title'] = '{}'.format(sys.argv[2])
    d['Author'] = 'xanthos'
    d['Subject'] = 'testing orbit integrator'
    d['Keywords'] = 'PdfPages multipage keywords author title subject'
    d['CreationDate'] = datetime.datetime(2022, 7, 10)
    d['ModDate'] = datetime.datetime.today()
