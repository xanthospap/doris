#! /usr/bin/python3

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
from matplotlib.ticker import PercentFormatter
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
def parse_orbit(fn):
    def strdct2dct(line):
      dct = {}
      for l in line.split():
        if ':' not in l:
          dct['site'] = l.strip()
        else:
          try:
            dct[l.split(':')[0]] = float(l.split(':')[1])
          except:
            dct[l.split(':')[0]] = l.split(':')[1]
      return dct
    
    dct = {}
    revolutions = []
    with open(fn, 'r') as fin:
        for line in fin.readlines():
            if line[0:5] == '[ORB]': 
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
                if t in dct:
                  dct[t].append(strdct2dct(' '.join(l[3:])))
                else:
                  dct[t] = [strdct2dct(' '.join(l[3:]))]
            
            elif line[0:5] == "[REV]":
                try:
                    t = datetime.datetime.strptime(' '.join([l[1],l[2][0:-3]]), '%Y-%m-%d %H:%M:%S.%f')
                    revolutions.append(t)
                except:
                    try:
                        t = datetime.datetime.strptime(' '.join([l[1],l[2]]), '%Y-%m-%d %H:%M:%S.%f')
                        revolutions.append(t)
                    except:
                        pass
    return revolutions, dct

## Assuming
## %Y-%m-%d %H:%M:%S.%f x y z Vx Vy Vz
## 0        1           2     5       
def parse_reference(fn):
    dct = {}
    with open(fn, 'r') as fin:
        for line in fin.readlines():
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
              labels=['x','y','z','vx','vy','vz']
              values=[float(k) for k in l[2:]]
              dct[t] = [dict(zip(labels, values))]
    return dct

def colAsArray(dct,col,fac=1e0,site=None):
    x=[];y=[];sy=[];
    ## if the requested column has a variance entry, it should be col_var
    varcol = col + '_var'
    for epoch,entries in dct.items():
      ## note that entries is a list of dictionaries ...
      for entry in entries:
        if col in entry:
          x.append(epoch)
          y.append(entry[col])
          if varcol in entry:
            sy.append(math.sqrt(entry[varcol]))
    return x,y,sy

def filter_dict_site(dct, site):
  newd = {}
  for t,entry in dct.items():
    for d in entry:
      if 'site' in d and site.upper() == d['site']:
        newd[t] = [d]
  return newd

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
    def strdct2dct(line):
      dct = {}
      for l in line.split():
        if ':' not in l:
          dct['site'] = l.strip()
        else:
          try:
            dct[l.split(':')[0]] = float(l.split(':')[1])
          except:
            dct[l.split(':')[0]] = l.split(':')[1]
      return dct
    ## Note that we ca have multiple entries for a given time-stamp. Hence, 
    ## a simple dictionary (i.e. key-value pairs) would be ambiguous. In this
    ## case, we will hold multiple values for each key, so that we can handle
    ## multiple entries for one epoch. The dictionary to return thus, will be:
    ## { t1: [{'beacon':...,'res',...,}, {'beacon':...,'res',...,}, ..., {'beacon':...,'res',...,}],
    ## { t2: [{'beacon':...,'res',...,}, {'beacon':...,'res',...,}, ..., {'beacon':...,'res',...,}],
    ## { tn: [{'beacon':...,'res',...,}, {'beacon':...,'res',...,}, ..., {'beacon':...,'res',...,}]}
    saa_passes = []
    revolutions = []
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
                    # [RES] 2022-09-05 00:01:51.854741246 HEMB -66.429022 [-66.429022 29979056.940859] (+1888.591 -1955.020 0.000e+00 2.000000 1.000000) 10.205
                    # 0     1          2                  3     4         5           6                 7          8        9         10       11        12
                    if cmax != datetime.datetime(1950,1,1) and t > cmax + datetime.timedelta(0, 5*60e0):
                        saa_passes.append([cmin,cmax])
                        cmin = datetime.datetime(2050,1,1)
                        cmax = datetime.datetime(1950,1,1)
                    if t < cmin: cmin = t
                    if t > cmax: cmax = t
                if t in dct:
                  dct[t].append(strdct2dct(' '.join(l[3:])))
                else:
                  dct[t] = [strdct2dct(' '.join(l[3:]))]
            
            elif line[0:5] == "[REV]":
                try:
                    t = datetime.datetime.strptime(' '.join([l[1],l[2][0:-3]]), '%Y-%m-%d %H:%M:%S.%f')
                    revolutions.append(t)
                except:
                    try:
                        t = datetime.datetime.strptime(' '.join([l[1],l[2]]), '%Y-%m-%d %H:%M:%S.%f')
                        revolutions.append(t)
                    except:
                        pass
              
    saa_passes.append([cmin,cmax])
    return saa_passes, revolutions, dct

def parse_block(fn):
    def strdct2dct(line):
      dct = {}
      for l in line.split():
        if ':' not in l:
          dct['site'] = l.strip()
        else:
          try:
            dct[l.split(':')[0]] = float(l.split(':')[1])
          except:
            dct[l.split(':')[0]] = l.split(':')[1]
      return dct
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
                  if t in dct:
                    dct[t].append(strdct2dct(' '.join(l[3:])))
                  else:
                    dct[t] = [strdct2dct(' '.join(l[3:]))]
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

parser.add_argument(
  '--plot-dynamics',
  action='store_true',
  dest='plot_dynamic_params')

def plot_residuals(fn, saveas=None):
  passes, revs, dct = parse_residuals(fn)

  ## Residuals Vs Time
  print("Plotting residuals Vs time ...", end='')
  fig, ax = plt.subplots()
  fac = 1e0
  t,y,_ = colAsArray(dct,'res',fac)
  ax.scatter(t,y,s=1,color='black')
  if revs != []:
    for r in revs:
      ax.axvline(x = r, color = 'green', label = 'Rev')
  ax.xaxis.set_major_formatter(mdates.DateFormatter("%H:%M"))
  ax.xaxis.set_major_locator(mdates.HourLocator(interval=4))
  ax.xaxis.set_minor_locator(mdates.HourLocator(interval=2))
  ax.grid(True, 'both', 'x')
  ax.set(xlabel='Time', ylabel='Residual [m/s]', title='Residuals')
  fig.suptitle('DORIS residuals @ {:}\n'.format(min(t).strftime('%Y-%m-%d')), fontsize=16)
  fig.autofmt_xdate()
  stats = ColStatistics(None,y,fac=1e0)
  ax.text(min(t), stats.minmax[0], r'#obs {:}, mean {:+.1f} $\pm$ {:.1f}'.format(stats.nobs, stats.mean, math.sqrt(stats.variance)), style='italic', bbox={'facecolor': 'red', 'alpha': 0.2, 'pad': 10})
  if saveas:
    fig.savefig(saveas + "-res.png")
  plt.show()
  print(" done")
  
  # Beacons per Block
  print("Plotting beacons per block ...", end='')
  dct = parse_block(fn)
  fac = 1
  f, (a0, a1) = plt.subplots(1, 2, gridspec_kw={'width_ratios': [3, 1]})
  t,y1,_ = colAsArray(dct,'obs',fac)
  t,y2,_ = colAsArray(dct,'ufilter',fac)
  a0.plot(t, y1, 'o-', color='blue', linewidth=1, markersize=2.5, label="beacons per block")
  a0.plot(t, y2, 'x-', color='black', linewidth=1, markersize=2.5, label="beacons used")
  a0.set(title='Beacons per Block', ylabel="Num Beacons", xlabel="Time")
  a0.legend()
  a0.xaxis.set_major_formatter(mdates.DateFormatter("%H:%M"))
  a0.xaxis.set_major_locator(mdates.HourLocator(interval=4))
  a0.xaxis.set_minor_locator(mdates.HourLocator(interval=2))
  a0.grid(True, 'both', 'x')
  a1.hist(y2, bins=np.arange(10)-.5, density=True, ls='dashed', fc=(.1,.1,.8,0.2), edgecolor='blue', label="beacons used")
  a1.hist(y1, bins=np.arange(10)-.5, density=True, ls='dashed', fc=(.8,0.01,.1,0.3), edgecolor='red', label="beacons per block")
  a1.set_xticks([0,1,2,3,4,5,6,7,8])
  a1.legend()
  if saveas:
    f.savefig(saveas + "-beacons.png")
  plt.show()
  print("done")

def obs_per_site(fn):
  _, _, dct = parse_residuals(fn)
  sites = []; nobs = []; mres = [];
  for ti,lofd in dct.items():
    for d in lofd:
      site = d['site']
      if site not in sites:
        sites.append(site)
        nobs.append(1)
        mres.append(d['res'])
      else:
        index = sites.index(site)
        nobs[index] += 1
        mres[index] += d['res']
  for i in range(len(mres)):
    mres[i] = mres[i] / nobs[i]
  fig, ax = plt.subplots()
  ax.bar(range(len(sites)), nobs)
  ax.set_xticklabels(sites)
  ax.grid(True, 'both', 'y')
  #ax.set(xlabel='Time', ylabel='Residual [m/s]', title='Residuals')
  #fig.suptitle('DORIS residuals @ {:}\n'.format(min(t).strftime('%Y-%m-%d')), fontsize=16)
  plt.show()

def plot_site(fn, site):
  saa_passes, revs, dct = parse_residuals(fn)
  dct = filter_dict_site(dct, site)
  fig, ax = plt.subplots(2,1)
  fac = 1e0

  t,y,err = colAsArray(dct,'res',fac)
  t1,ypred,_= colAsArray(dct,'prediction',fac)
  if revs != []:
    for r in revs:
      ax[0].axvline(x=r,color='green',label='Rev')
  print('{:} vs {:} vs {:}'.format(len(t), len(y), len(err)))
  ax[0].errorbar(t,y,yerr=err,fmt='none', ecolor='red', elinewidth=.1, alpha=.1)
  ax[0].scatter(t,y,s=1,color='black')
  ax[0].scatter(t1,ypred,s=1,color='blue')
  ax[0].set_title('Residuals [m/s]')
  ax[0].xaxis.set_major_formatter(mdates.DateFormatter("%H:%M"))
  ax[0].xaxis.set_major_locator(mdates.HourLocator(interval=3))
  ax[0].xaxis.set_minor_locator(mdates.HourLocator(interval=1))
  ax[0].grid(True, 'both', 'x')
  
  t,y,_ = colAsArray(dct,'df',fac)
  if revs != []:
    for r in revs:
      ax[1].axvline(x=r,color='green',label='Rev')
  ax[1].scatter(t,y,s=1,color='black')
  ax[1].set_title('Df/f [-]')
  ax[1].xaxis.set_major_formatter(mdates.DateFormatter("%H:%M"))
  ax[1].xaxis.set_major_locator(mdates.HourLocator(interval=3))
  ax[1].xaxis.set_minor_locator(mdates.HourLocator(interval=1))
  ax[1].grid(True, 'both', 'x')
  
  fig.autofmt_xdate()
  fig.suptitle('Site {:}@{:}\n'.format(site, t[0].strftime('%Y-%m-%d')), fontsize=16)
  plt.show()

  """
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
  """

def plot_dynamic_params(fn, saveas=None, shadow_passes_fn=None):
  saa_passes, revs, dct = parse_residuals(fn)
  if shadow_passes_fn:
    shd_passes = shadow_passes(shadow_passes_fn)
  else:
    shd_passes = []

  fig, ax = plt.subplots(2,1)
  fac = 1e0
  
  t,y,err = colAsArray(dct,'cd',fac)
  ax[0].errorbar(t,y,yerr=err,fmt='none', ecolor='red', elinewidth=.1, alpha=.1)
  ax[0].scatter(t,y,s=2,color='black')
  ax[0].set_title(r'Drag Coefficient $C_d$')
  if revs != []:
    for r in revs:
      ax[0].axvline(x = r, color = 'green', label = 'Rev')
  ax[0].xaxis.set_major_formatter(mdates.DateFormatter("%H:%M"))
  ax[0].xaxis.set_major_locator(mdates.HourLocator(interval=3))
  ax[0].xaxis.set_minor_locator(mdates.HourLocator(interval=1))
  ax[0].grid(True, 'both', 'x')
  
  t,y,err = colAsArray(dct,'cr',fac)
  ax[1].scatter(t,y,s=1,color='black')
  ax[1].errorbar(t,y,yerr=err,fmt='none', ecolor='red', elinewidth=.1, alpha=.1)
  if revs != []:
    for r in revs:
      ax[1].axvline(x = r, color = 'green', label = 'Rev')
  ax[1].set_title(r'Srp Coefficient $C_r$')
  ax[1].xaxis.set_major_formatter(mdates.DateFormatter("%H:%M"))
  ax[1].xaxis.set_major_locator(mdates.HourLocator(interval=3))
  ax[1].xaxis.set_minor_locator(mdates.HourLocator(interval=1))
  ax[1].grid(True, 'both', 'x')
  
  fig.autofmt_xdate()
  t0 = min(t)
  fig.suptitle('Dynamic Parameters {:}\n'.format(t0.strftime('%Y-%m-%d')), fontsize=16)
  plt.show()

def plot_state_diffs(fnref, fntest, save_as):
  def whichCol(component, posvel):
      labels = ['x','y','z','vx','vy','vz']
      return labels[component + int(posvel==1)*3]
  def whichTitle(component, posvel, fac=1e0):
      title = ' [{:}m]'.format('k' if fac==1e-3 else '')
      if posvel == 1:
          title = ' Velocity [{:}m/sec]'.format('k' if fac == 1e-3 else '')
      return ['X','Y','Z'][component] + title
  def make_state_diffs(ref_dct, dct):
      resdct = {}
      for t,d1 in ref_dct.items():
          if t not in dct:
              print('## Warning! failed to find record for date: {:}'.format(t))
          else:
              if len(dct[t]) != 1:
                print('## Warning! more than one estimates found for block')
              d1 = d1[0]
              d2 = dct[t][-1]
              resdct[t] = [{}]
              for key,val2 in d2.items():
                if key in d1:
                  resdct[t][0][key] = val2 - d1[key]
                else:
                  resdct[t][0][key] = val2
      return resdct

  dct1 = parse_reference(fnref)
  revs, dct2 = parse_orbit(fntest)
  diffs = make_state_diffs(dct1,dct2)
  shd_passes = shadow_passes(fntest)

  fig, axs = plt.subplots(3, 2, sharey='col', figsize=(10, 6), constrained_layout=True)
  ## x-, y-, z-components ...
  for component in range(3):
      ## first position, then velocity ...
      for pv in range(2):
          fac = 1e0
          key = whichCol(component, pv)
          if revs != []:
            for r in revs:
              axs[component, pv].axvline(x=r,color='green',label='Rev')
          t,y,err = colAsArray(diffs,key,fac)
          #if (component+pv==0): axs[component, pv].errorbar(t,y,yerr=err, fmt='none', ecolor='red', elinewidth=.1, alpha=.1)
          axs[component, pv].scatter(t,y,s=1,color='black')
          axs[component, pv].set_title(whichTitle(component, pv, fac))
          axs[component, pv].xaxis.set_major_formatter(mdates.DateFormatter("%H:%M"))
          axs[component, pv].xaxis.set_major_locator(mdates.HourLocator(interval=3))
          axs[component, pv].xaxis.set_minor_locator(mdates.HourLocator(interval=1))
          axs[component, pv].yaxis.set_minor_locator(MultipleLocator(1))
          dsts = ColStatistics(None,y,fac)
          if pv == 0:
            axs[component, pv].text(t[0], dsts.minmax[0], r'Mean: {:+.1f} $\pm$ {:.1f}'.format(dsts.mean, math.sqrt(dsts.variance)), style='italic', bbox={'facecolor': 'red', 'alpha': 0.5, 'pad': 10})
          else:
            axs[component, pv].text(t[0], dsts.minmax[0], r'Mean: {:+.4f} $\pm$ {:.4f}'.format(dsts.mean, math.sqrt(dsts.variance)), style='italic', bbox={'facecolor': 'red', 'alpha': 0.5, 'pad': 10})
          axs[component, pv].grid(True, 'both', 'x', color='grey', linewidth=.1)
          axs[component, pv].grid(True, 'both', 'y', color='green', linewidth=.2)
          if shd_passes != []:
              for intrv in shd_passes:
                axs[component, pv].axvspan(intrv[0], intrv[1], alpha=0.2, color='red')
          #if saa_passes != []:
          #    for intrv in saa_passes:
          #      axs[component, pv].axvspan(intrv[0], intrv[1], alpha=0.1, color='green')
  ## Rotate date labels automatically
  fig.autofmt_xdate()
  fig.suptitle('State Differences Estimated - Reference\n{:}'.format(t[0].strftime('%Y-%m-%d')), fontsize=16)
  if save_as:
      print('Saving figure to {:}'.format(save_as))
      plt.savefig(save_as)
  plt.show()

if __name__ == '__main__':

    args = parser.parse_args()

    if args.ref_state:
      plot_state_diffs(args.ref_state, args.input, args.save_as)

    if args.plot_res:
      plot_residuals(args.input)

    if args.plot_dynamic_params:
      plot_dynamic_params(args.input)

    if args.sites:
        for s in args.sites:
            plot_site(args.input, s)

    obs_per_site(args.input)
