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
from scipy import stats
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.ticker import MultipleLocator
import julian

## VMF3 DORIS-site-grid URI
VMF3_DG = "https://vmf.geo.tuwien.ac.at/trop_products/DORIS/V3GR/V3GR_OP/daily/"
def daily_vmf3_fn(t):
    url = VMF3_DG + '/{:4d}/'.format(t.year)
    fn = '{:4d}{:}.v3gr_d'.format(t.year,t.strftime("%j"))
    return url, fn

def fetch_vmf3_doris_grid(t):
    url,fn = daily_vmf3_fn(t)
    r = requests.get(url+fn, allow_redirects=True)
    open(fn, 'wb').write(r.content)
    return

def get_daily_vmf3_for(site, t0, t1, fn):
    site = site.upper()
    mjd0 = julian.to_jd(t0, fmt='mjd')
    mjd1 = julian.to_jd(t1, fmt='mjd')
    t = []
    lwz = []
    with open(fn, 'r') as fin:
        for line in fin.readlines():
            if line[0:4] == site:
                l = line.split()
                tc = julian.from_jd(float(l[1]), fmt='mjd')
                if (tc-t0).total_seconds()/3600e0 <= 18 and (t1-tc).total_seconds()/3600e0 <= 18:
                    t.append(julian.from_jd(float(l[1]), fmt='mjd'))
                    lwz.append(float(l[5]))
    return t,lwz

def interpolate(t, v, t0, t1):
    mjd = [ julian.to_jd(ti, fmt='mjd') for ti in t ]
    f = interp1d(mjd, v)
    newx = np.linspace(mjd[0], mjd[-1], num=50, endpoint=True)
    return [ julian.from_jd(ti, fmt='mjd') for ti in newx ], f(newx)

def trim_excess_epochs(dct, max_hours):
    if max_hours == sys.float_info.max: return dct
    t0 = min(dct)
    new_dct = {}
    for t,vals in dct.items():
        if t - t0 <= datetime.timedelta(hours=max_hours):
            new_dct[t] = vals
    return new_dct

## Assuming
## %Y-%m-%d %H:%M:%S.%f x y z Vx Vy Vz
## 0        1           2 
def parse(fn, max_hours):
    dct = {}
    with open(fn, 'r') as fin:
        for line in fin.readlines():
            if line[0] != '#': 
              l = line.split()
              line_is_valid = False
              try:
                  t = datetime.datetime.strptime(' '.join([l[0],l[1][0:-3]]), '%Y-%m-%d %H:%M:%S.%f')
                  line_is_valid = True
              except:
                  try:
                    t = datetime.datetime.strptime(' '.join([l[0],l[1]]), '%Y-%m-%d %H:%M:%S.%f')
                    line_is_valid = True
                  except:
                    pass
              if line_is_valid:
                #x,y,z,vx,vy,vz,dff,lw,c,res,el = [float(k) for k in l[5:16]]
                #beacon = l[3]
                #arc_nr = int(l[4])
                #dct[t] = [beacon,arc_nr,x,y,z,vx,vy,vz,dff,lw,c,res,el]
                x,y,z,vx,vy,vz = [float(k) for k in l[2:]]
                dct[t] = [x,y,z,vx,vy,vz]
    return trim_excess_epochs(dct, max_hours)

## Assuming
## %Y-%m-%d %H:%M:%S.%f x y z Vx Vy Vz
## 0        1           2     5       
def parse_reference(fn, max_hours):
    dct = {}
    with open(fn, 'r') as fin:
        for line in fin.readlines():
            l = line.split()
            try:
                t = datetime.datetime.strptime(' '.join([l[0],l[1][0:-3]]), '%Y-%m-%d %H:%M:%S.%f')
            except:
                t = datetime.datetime.strptime(' '.join([l[0],l[1]]), '%Y-%m-%d %H:%M:%S.%f')
            x,y,z,vx,vy,vz = [float(k) for k in l[2:]]
            dct[t] = [x,y,z,vx,vy,vz]
    return trim_excess_epochs(dct, max_hours)

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
    return stats.describe(array)

def remove_outliers(dct, col, fac=1e0):
  #min, max, sumDn, sigma = ColStatistics(dct,col,fac)
  vstats = ColStatistics(dct,col,fac)
  sigma = math.sqrt(vstats.variance)
  xs=[];ys=[]
  for t,vals in dct.items():
    if abs(vals[col]*fac) < 3e0 * sigma:
      xs.append(t)
      ys.append(vals[col]*fac)
  vstats = ColStatistics(None,ys,fac)
  sigma = math.sqrt(vstats.variance)
  xsnew=[];ysnew=[]
  for entry in zip(xs,ys):
    if abs(entry[1]*fac) < 3e0 * sigma:
      xsnew.append(entry[0])
      ysnew.append(entry[1]*fac)
  return xsnew,ysnew

# [SHDW ] Changing shadow status from 0.000 to 0.013 at 59827.871805994
def shadow_passes(fn):
    passes = []
    with open(fn) as fin:
        cmin = 1e22
        cmax = 1e-22
        for line in fin.readlines():
            if line[0:7] == "[SHDW ]":
                l = line.split()
                t = float(l[-1])
                if cmax != 1e-22 and t > cmax + 0.05:
                    #passes.append([cmin,cmax])
                    passes.append([julian.from_jd(cmin, fmt='mjd'), julian.from_jd(cmax, fmt='mjd')])
                    cmin = 1e22
                    cmax = 1e-22
                if t < cmin: cmin = t
                if t > cmax: cmax = t
    passes.append([julian.from_jd(cmin, fmt='mjd'), julian.from_jd(cmax, fmt='mjd')])
    return passes

def parse_residuals(fn):
    saa_passes = []
    dct = {}
    with open(fn) as fin:
        cmin = datetime.datetime(2050,1,1)
        cmax = datetime.datetime(1950,1,1)
        for line in fin.readlines():
            if line[0:5] == "[RES]":
                l = line.split()
                line_is_valid = False
                try:
                    t = datetime.datetime.strptime(' '.join([l[1],l[2][0:-3]]), '%Y-%m-%d %H:%M:%S.%f')
                    line_is_valid = True
                except:
                    try:
                        t = datetime.datetime.strptime(' '.join([l[1],l[2]]), '%Y-%m-%d %H:%M:%S.%f')
                        line_is_valid = True
                    except:
                        pass
                if line_is_valid:
                    if cmax != datetime.datetime(1950,1,1) and t > cmax + datetime.timedelta(0, 5*60e0):
                        saa_passes.append([cmin,cmax])
                        cmin = datetime.datetime(2050,1,1)
                        cmax = datetime.datetime(1950,1,1)
                    if t < cmin: cmin = t
                    if t > cmax: cmax = t
                dct[t] = {'beacon': l[3], 'res': float(l[4]), 'df': float(l[9]), 'cd': float(l[10]), 'cr': float(l[11][:-1]), 'el': float(l[12])}
    saa_passes.append([cmin,cmax])
    return saa_passes, dct

def parse_block(fn):
    dct = {}
    with open(fn) as fin:
        for line in fin.readlines():
            if line[0:5] == "[BLC]":
                l = line.split()
                line_is_valid = False
                try:
                    t = datetime.datetime.strptime(' '.join([l[1],l[2][0:-3]]), '%Y-%m-%d %H:%M:%S.%f')
                    line_is_valid = True
                except:
                     try:
                        t = datetime.datetime.strptime(' '.join([l[1],l[2]]), '%Y-%m-%d %H:%M:%S.%f')
                        line_is_valid = True
                     except:
                        pass
                if line_is_valid:
                    dct[t] = {'obs': int(l[4]), 'lowele': int(l[6]), 'saa': int(l[8]), 'oc': int(l[10]), 'filter': int(l[12])}
    return dct

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
    '-o',
    '--save-plot',
    metavar='SAVE_PLOT_AS',
    dest='save_as',
    default=None,
    required=False,
    help='Save plot as')

parser.add_argument(
    '-c',
    '--shadow',
    metavar='SHADOW_PASSES',
    dest='shadow_pass_fn',
    default=None,
    required=False,
    help='File to extract shadow passes')

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
    '--max-hours',
    metavar='MAX_HOURS',
    dest='max_hours',
    type=float,
    default=sys.float_info.max,
    required=False,
    help='Plot results up to MAX HOURS')

parser.add_argument(
  '--plot-residuals',
  action='store_true',
  dest='plot_res')

def plot_residuals(fn):
  passes, dct = parse_residuals(fn)
  t0 = datetime.datetime(2050,1,1)
  for t in dct:
    if t < t0: t0 = t

  t = [ ti for ti in dct ]
  fig, ax = plt.subplots(2,1)
  fac = 1e0
  ax[0].scatter(t,colAsArray(dct,'res',fac),s=1,color='black')
  ax[0].set_title('Residuals [m]')
  ax[0].xaxis.set_major_formatter(mdates.DateFormatter("%H:%M"))
  ax[0].xaxis.set_major_locator(mdates.HourLocator(interval=4))
  ax[0].xaxis.set_minor_locator(mdates.HourLocator(interval=2))
  ax[0].grid(True, 'both', 'x')
  # 
  dct = parse_block(fn)
  for t in dct: 
      if t < t0: t0 = t
  t = [ ti for ti in dct ]
  ax[1].plot(t,colAsArray(dct,'obs',fac),color='black')
  ax[1].plot(t,colAsArray(dct,'oc',fac),color='red')
  ax[1].plot(t,colAsArray(dct,'filter',fac),color='green')
  ax[1].set_title('Beacons per Block')
  ax[1].xaxis.set_major_formatter(mdates.DateFormatter("%H:%M"))
  ax[1].xaxis.set_major_locator(mdates.HourLocator(interval=4))
  ax[1].xaxis.set_minor_locator(mdates.HourLocator(interval=2))
  ax[1].grid(True, 'both', 'x')
  fig.autofmt_xdate()
  #
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
  tmin = min(t1[0],t2[0])
  tmax = max(t1[-1], t2[-1])
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

  _,vmf3fn = daily_vmf3_fn(tmin)
  if not os.path.isfile(vmf3fn):
      fetch_vmf3_doris_grid(tmin)
  tv,lwzv = get_daily_vmf3_for(site, tmin, tmax, vmf3fn)
  tv,lwzv = interpolate(tv, lwzv, tmin, tmax)

  ax[1].scatter(t2,lwz,s=1,color='black')
  ax[1].scatter(tv,lwzv,s=.5,color='red')
  ax[1].set_title('Lw Zenith [m]')
  ax[1].xaxis.set_major_formatter(mdates.DateFormatter("%H:%M"))
  ax[1].xaxis.set_major_locator(mdates.HourLocator(interval=3))
  ax[1].xaxis.set_minor_locator(mdates.HourLocator(interval=1))
  ax[1].grid(True, 'both', 'x')
  fig.autofmt_xdate()
  fig.suptitle('Site {:}@{:}\n'.format(site, tmin.strftime('%Y-%m-%d')), fontsize=16)
  plt.show()

def plot_state_diffs(fnref, fntest, max_hours, save_as, shadow_passes_fn):
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
              vals2 = dct[t]
              diffs = [ x[0]-x[1] for x in zip(vals,vals2)]
              resdct[t] = diffs
      return resdct

  try:
    dct1 = parse_reference(fnref, max_hours)
    print('[DEBUG] File {:} parsed as reference'.format(fnref))
  except:
    dct1 = parse(fnref, max_hours)
    print('[DEBUG] File {:} parsed as data'.format(fnref))
  dct2 = parse(fntest, max_hours)
  diffs = make_state_diffs(dct1,dct2)

  if shadow_passes_fn:
    shd_passes = shadow_passes(shadow_passes_fn)
    try:
        saa_passes,_ = parse_residuals(shadow_passes_fn)
    except:
        print('[DEBUG] Tried but failed to extract saa passes from file {:}'.format(fn), file=sys.stderr)
        saa_passes = []
  else:
    shd_passes = []
    saa_passes = []

  t0 = datetime.datetime.max
  for t in dct1:
    if t < t0: t0 = t

  fig, axs = plt.subplots(3, 2, figsize=(10, 6), constrained_layout=True)
  t = [ ti for ti in diffs ]
  ## x-, y-, z-components ...
  for component in range(3):
      ## first position, then velocity ...
      for pv in range(2):
          fac = 1e0
          index = whichCol(component, pv)
          axs[component, pv].scatter(t,colAsArray(diffs,index,fac),s=1,color='black')
          axs[component, pv].set_title(whichTitle(component, pv, fac))
          _, end = axs[component, pv].get_xlim()
          axs[component, pv].xaxis.set_major_formatter(mdates.DateFormatter("%H:%M"))
          axs[component, pv].xaxis.set_major_locator(mdates.HourLocator(interval=3))
          axs[component, pv].xaxis.set_minor_locator(mdates.HourLocator(interval=1))
          dsts = ColStatistics(diffs,index,fac)
          text_x = t[0] #t[-1] - (t[-1] - t[0]) / 10e0
          text_y = dsts.minmax[0] # dsts.minmax[1] - (dsts.minmax[1]-dsts.minmax[0]) / 10e0
          if pv == 0:
            axs[component, pv].text(text_x, text_y, 'Mean: {:+.1f} +/- {:.1f}'.format(dsts.mean, math.sqrt(dsts.variance)), style='italic', bbox={'facecolor': 'red', 'alpha': 0.5, 'pad': 10})
          else:
            axs[component, pv].text(text_x, text_y, 'Mean: {:+.4f} +/- {:.4f}'.format(dsts.mean, math.sqrt(dsts.variance)), style='italic', bbox={'facecolor': 'red', 'alpha': 0.5, 'pad': 10})
          axs[component, pv].grid(True, 'both', 'x')
          if shd_passes != []:
              for intrv in shd_passes:
                axs[component, pv].axvspan(intrv[0], intrv[1], alpha=0.2, color='red')
          if saa_passes != []:
              for intrv in saa_passes:
                axs[component, pv].axvspan(intrv[0], intrv[1], alpha=0.1, color='green')
  ## Rotate date labels automatically
  fig.autofmt_xdate()
  fig.suptitle('Sp3 - Integrator Diffs. at {:}\n'.format(t0.strftime('%Y-%m-%d')), fontsize=16)
  if save_as:
      print('Saving figure to {:}'.format(save_as))
      plt.savefig(save_as)
  plt.show()

if __name__ == '__main__':

    args = parser.parse_args()

    if args.ref_state:
      plot_state_diffs(args.ref_state, args.input, args.max_hours, args.save_as, args.shadow_pass_fn)

    if args.plot_res:
      plot_residuals(args.shadow_pass_fn)

    if args.sites:
        for s in args.sites:
            plot_site(args.input, s)
