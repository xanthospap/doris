#! /usr/bin/python
import datetime
import sys
import math
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import argparse
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.ticker import MultipleLocator

## Assuming
## %Y-%m-%d %H:%M:%S.%f TAI Beacon Arc x y z Vx Vy Vz Df/f Lw C Res El Mjd
## 0        1           2   3      4   5     8        11   12   14  15 16
def parse(fn):
    dct = {}
    with open(fn, 'r') as fin:
        for line in fin.readlines():
            l = line.split()
            try:
                t = datetime.datetime.strptime(' '.join([l[0],l[1][0:-3]]), '%Y-%m-%d %H:%M:%S.%f')
            except:
                t = datetime.datetime.strptime(' '.join([l[0],l[1]]), '%Y-%m-%d %H:%M:%S.%f')
            x,y,z,vx,vy,vz,dff,lw,c,res,el,mjd = [float(k) for k in l[5:]]
            beacon = l[3]
            arc_nr = int(l[4])
            dct[t] = [beacon,arc_nr,x,y,z,vx,vy,vz,dff,lw,c,res,el,mjd]
    return dct

## Assuming
## %Y-%m-%d %H:%M:%S.%f x y z Vx Vy Vz Mjd
## 0        1           2     5        8
def parse_reference(fn):
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

def colAsArray(dct,col,fac=1e0,site=None):
    if site:
      return [ (t,vals[col]*fac) for t,vals in dct.items() if vals[0]==site]
    else:
      return [ vals[col]*fac for _,vals in dct.items() ]

def reduce(dates, dct):
    dctcp = {}
    for date in dates:
        if date in dct:
            dctcp[date] = dct[date]
        else:
            raise RuntimeError('ERROR WTF date')
    return dctcp

def ColStatistics(dct,col,fac=1e0):
    if dct:
      array = colAsArray(dct,col,fac)
    else:
      array = col
    sum = 0e0
    sumSquared = 0e0
    max = -1e99
    min = 1e99
    for value in array:
        sum += value
        sumSquared += value * value
        if value > max: max = value
        if value < min: min = value
    n = float(len(array))
    return min, max, sum / n, math.sqrt(sumSquared / n)

def remove_outliers(dct, col, fac=1e0):
  min, max, sumDn, sigma = ColStatistics(dct,col,fac)
  xs=[];ys=[]
  for t,vals in dct.items():
    if abs(vals[col]*fac) < 3e0 * sigma:
      xs.append(t)
      ys.append(vals[col]*fac)
  min, max, sumDn, sigma = ColStatistics(None,ys,fac)
  xsnew=[];ysnew=[]
  for entry in zip(xs,ys):
    if abs(entry[1]*fac) < 3e0 * sigma:
      xsnew.append(entry[0])
      ysnew.append(entry[1]*fac)
  return xsnew,ysnew

class myFormatter(argparse.ArgumentDefaultsHelpFormatter,
                  argparse.RawTextHelpFormatter):
    pass

parser = argparse.ArgumentParser(
    formatter_class=myFormatter,
    description=
    'Plot Script for DORIS POD tests',
    epilog=('''National Technical University of Athens,
    Dionysos Satellite Observatory\n
    Send bug reports to:
    Xanthos Papanikolaou, xanthos@mail.ntua.gr
    Aug, 2022'''))

parser.add_argument(
    'input',
    metavar='INPUT',
    help='File holding results')

parser.add_argument(
    '-r',
    '--ref-state',
    metavar='REF_STATE',
    dest='ref_state',
    default=None,
    required=False,
    help='File holding reference state')

parser.add_argument(
    '-s',
    '--sites',
    metavar='SITES',
    dest='sites',
    nargs='+',
    default=None,
    required=False,
    help='File holding reference state')

parser.add_argument(
  '--plot-residuals',
  action='store_true',
  dest='plot_res')

def plot_residuals(fn):
  residuals_index = 11
  elevation_index = residuals_index + 1
  dct = parse(fn)
  t0 = datetime.datetime.max
  for t in dct:
    if t < t0: t0 = t

  t = [ ti for ti in dct ]
  fig, ax = plt.subplots(2,2)
  fac = 1e0
  ax[0,0].scatter(t,colAsArray(dct,residuals_index,fac),s=1,color='black')
  ax[0,0].set_title('Residuals [m]')
  ax[0,0].xaxis.set_major_formatter(mdates.DateFormatter("%H:%M"))
  ax[0,0].xaxis.set_major_locator(mdates.HourLocator(interval=3))
  ax[0,0].xaxis.set_minor_locator(mdates.HourLocator(interval=1))
  ax[0,0].grid(True, 'both', 'x')
  # get statistics and remove outliers
  ts,rs=remove_outliers(dct, residuals_index)
  ax[1,0].scatter(ts,rs,s=1,color='black')
  ax[1,0].set_title('Residuals [m]')
  ax[1,0].xaxis.set_major_formatter(mdates.DateFormatter("%H:%M"))
  ax[1,0].xaxis.set_major_locator(mdates.HourLocator(interval=3))
  ax[1,0].xaxis.set_minor_locator(mdates.HourLocator(interval=1))
  ax[1,0].grid(True, 'both', 'x')
  fig.autofmt_xdate()
  #
  # counts, bins = np.histogram(colAsArray(dct,residuals_index,fac))
  n_bins = 20
  ax[0,1].hist(colAsArray(dct,residuals_index,fac), bins=n_bins)
  # elevation Vs res
  ax[1,1].scatter(colAsArray(dct,elevation_index,fac),colAsArray(dct,residuals_index,fac),s=1,color='black')
  #
  fig.suptitle('DORIS residuals @ {:}\n'.format(t0.strftime('%Y-%m-%d')), fontsize=16)
  plt.show()

def plot_site(fn, site):
  site = site.upper()
  dct = parse(fn)
  tdff = colAsArray(dct,8,1e0,site)
  tlwz = colAsArray(dct,9,1e0,site)
  t1 = [ x[0] for x in tdff ]
  dff = [ x[1] for x  in tdff ]
  t2 = [ x[0] for x in tlwz ]
  lwz = [ x[1] for x  in tlwz ]
  t0 = min(t1[0],t2[0])
  if tdff == []:
      print("No data for site: {:}".format(site), file=sys.stderr)
  fig, ax = plt.subplots(2,1)
  fac = 1e0
  ax[0].scatter(t1,dff,s=1,color='black')
  ax[0].set_title('Df/f [-]')
  ax[0].xaxis.set_major_formatter(mdates.DateFormatter("%H:%M"))
  ax[0].xaxis.set_major_locator(mdates.HourLocator(interval=3))
  ax[0].xaxis.set_minor_locator(mdates.HourLocator(interval=1))
  ax[0].grid(True, 'both', 'x')
  ax[1].scatter(t2,lwz,s=1,color='black')
  ax[1].set_title('Lw Zenith [m]')
  ax[1].xaxis.set_major_formatter(mdates.DateFormatter("%H:%M"))
  ax[1].xaxis.set_major_locator(mdates.HourLocator(interval=3))
  ax[1].xaxis.set_minor_locator(mdates.HourLocator(interval=1))
  ax[1].grid(True, 'both', 'x')
  fig.autofmt_xdate()
  fig.suptitle('Site {:}@{:}\n'.format(site, t0.strftime('%Y-%m-%d')), fontsize=16)
  plt.show()

def plot_state_diffs(fnref, fntest):
  def whichCol(component, posvel):
      return component + int(posvel==1)*3
  def whichTitle(component, posvel, fac=1e0):
      title = ' [{:}m]'.format('k' if fac==1e-3 else '')
      if posvel == 1:
          title = ' Velocity [{:}m/sec]'.format('k' if fac == 1e-3 else '')
      return ['X','Y','Z'][component] + title
  def make_state_diffs(ref_dct, dct):
      resdct = {}
      for t,vals in ref_dct.items():
          if t not in dct:
              print('## Warning! failed to find record for date: {:}'.format(t))
          else:
              vals2 = dct[t][2:8] + [dct[t][-1]]
              if vals[-1] != vals2[-1]: print('Ops! difference is {:.15e}'.format(abs(vals[-1] - vals2[-1])))
              assert(abs(vals[-1] - vals2[-1]) < 1e-9)
              diffs = [ x[0]-x[1] for x in zip(vals,vals2)]
              resdct[t] = diffs
      return resdct

  dct1 = parse_reference(fnref)
  dct2 = parse(fntest)
  diffs = make_state_diffs(dct1,dct2)

  t0 = datetime.datetime.max
  for t in dct1:
    if t < t0: t0 = t

  fig, axs = plt.subplots(3, 2, figsize=(10, 6), constrained_layout=True)
  t = [ ti for ti in diffs ]
  ## x-, y-, z-components ...
  for component in range(3):
      ## first position, then velocity ...
      for pv in range(2):
          # fac = 1e-3 if pv == 0 else 1e0 ## noop .... just meters now
          fac = 1e0
          index = whichCol(component, pv)
          #axs[component, pv].scatter(t2sec(t,t0),colAsArray(diffs,index,fac),s=1,color='black')
          axs[component, pv].scatter(t,colAsArray(diffs,index,fac),s=1,color='black')
          axs[component, pv].set_title(whichTitle(component, pv, fac))
          _, end = axs[component, pv].get_xlim()
          #axs[component, pv].xaxis.set_minor_locator(MultipleLocator(3600))
          #axs[component, pv].xaxis.set_ticks(np.arange(start_sec, end, 3600*2))
          axs[component, pv].xaxis.set_major_formatter(mdates.DateFormatter("%H:%M"))
          axs[component, pv].xaxis.set_major_locator(mdates.HourLocator(interval=3))
          axs[component, pv].xaxis.set_minor_locator(mdates.HourLocator(interval=1))
          minv, maxv, average, rms = ColStatistics(diffs,index,fac)
          text_x = t[-1] - (t[-1] - t[0]) / 10e0
          text_y = maxv - (maxv-minv) / 10e0
          axs[component, pv].text(text_x, text_y, 'Average: {:+.4f}, Rms: {:.4f}'.format(average, rms), style='italic',
          bbox={'facecolor': 'red', 'alpha': 0.5, 'pad': 10})
          axs[component, pv].grid(True, 'both', 'x')
  ## Rotate date labels automatically
  fig.autofmt_xdate()
  fig.suptitle('Sp3 - Integrator Diffs. at {:}\n'.format(t0.strftime('%Y-%m-%d')), fontsize=16)
  plt.show()

if __name__ == '__main__':

    args = parser.parse_args()

    if args.ref_state:
      plot_state_diffs(args.ref_state, args.input)

    if args.plot_res:
      plot_residuals(args.input)

    if args.sites:
        for s in args.sites:
            plot_site(args.input, s)
