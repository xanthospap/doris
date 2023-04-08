#! /usr/bin/python3

import datetime
import sys, os
import math
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import argparse
import requests
import copy
from scipy.interpolate import interp1d
from scipy import stats
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.ticker import MultipleLocator
from matplotlib.ticker import PercentFormatter
import scipy.signal as signal
import julian

width_pt = 418.25368

# Using seaborn's style
plt.style.use('seaborn')

tex_fonts = {
    # Use LaTeX to write all text
    "text.usetex": True,
    "font.family": "serif",
    # Use 10pt font in plots, to match 10pt font in document
    "axes.labelsize": 10,
    "font.size": 10,
    # Make the legend/label fonts a little smaller
    "legend.fontsize": 8,
    "xtick.labelsize": 8,
    "ytick.labelsize": 8
}
tex_fonts_small = {
    # Use LaTeX to write all text
    "text.usetex": True,
    "font.family": "serif",
    # Use 10pt font in plots, to match 10pt font in document
    "axes.labelsize": 8,
    "font.size": 6,
    # Make the legend/label fonts a little smaller
    "legend.fontsize": 11,
    "xtick.labelsize": 6,
    "ytick.labelsize": 6
}

plt.rcParams.update(tex_fonts)
def set_size(width, fraction=1, subplots=(1, 1)):
    """Set figure dimensions to avoid scaling in LaTeX.
       see https://jwalton.info/Embed-Publication-Matplotlib-Latex/

    Parameters
    ----------
    width: float or string
            Document width in points, or string of predined document type
    fraction: float, optional
            Fraction of the width which you wish the figure to occupy
    subplots: array-like, optional
            The number of rows and columns of subplots.
    Returns
    -------
    fig_dim: tuple
            Dimensions of figure in inches
    """
    # Width of figure (in pts)
    fig_width_pt = width_pt * fraction
    # Convert from pt to inches
    inches_per_pt = 1 / 72.27
    # Golden ratio to set aesthetic figure height
    # https://disq.us/p/2940ij3
    golden_ratio = (5**.5 - 1) / 2
    # Figure width in inches
    fig_width_in = fig_width_pt * inches_per_pt
    # Figure height in inches
    fig_height_in = fig_width_in * golden_ratio * (subplots[0] / subplots[1])
    return (fig_width_in, fig_height_in)

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

def get_beacon_type(snxfn, t):
    def secday2hms(sdstr):
        return str(datetime.timedelta(seconds=float(sdstr)))
    dct = {}
    with open(snxfn, 'r') as fin:
        line = fin.readline()
        # %=SNX 2.02 IDS 22:329:00000 IDS 93:003:00000 22:002:00000 D 01872 2 X V
        l = line.split()
        snxend = datetime.datetime.strptime(':'.join([l[6][0:6], secday2hms(l[6][7:])]), '%y:%j:%H:%M:%S')
        while line and line.strip() != '+SITE/ANTENNA':
            line = fin.readline()
        line = fin.readline()
        assert (line.strip() == '*Code PT SOLN T Data_start__ Data_end____ Description_________ S/N__')
        line = fin.readline()
        while line and line.strip() != '-SITE/ANTENNA':
            # ADEA  A    1 D 93:003:00000 98:084:11545              ALCATEL -----
            l = line.split()
            beacon = l[0]
            start = datetime.datetime.strptime(':'.join([l[4][0:6], secday2hms(l[4][7:])]), '%y:%j:%H:%M:%S')
            stop  = datetime.datetime.strptime(':'.join([l[5][0:6], secday2hms(l[5][7:])]), '%y:%j:%H:%M:%S')
            if (snxend-stop).total_seconds() < 3600: stop = datetime.datetime.now()
            if t >= start and t < stop:
                dct[beacon] = l[6]
            line = fin.readline()
    return dct

def remove_outliers(x,y,max_remove_percent):
  s = stats.describe(y)
  std = math.sqrt(s.variance)
  fac = 3e0; it = 0; minworks=1e25; rem_indexes = [];
  keepon = True
  while (keepon):
    print('[DEBUG] Outlier detection with factor: {:.1f}'.format(fac))
    ouliers = 0
    new_indexes = []
    for i in range(len(y)):
      if abs(y[i] - s.mean) > fac * std:
        new_indexes.append(i)
    percent = len(new_indexes) * 100e0 / len(y)
    if percent < max_remove_percent/2e0:
      print('[DEBUG] Factor {:.1f} only removes about {:.2f}% as outliers; trying for smaller factor'.format(fac, percent))
      fac -= .2
      rem_indexes = copy.deepcopy(new_indexes)
    elif percent < max_remove_percent:
      print('[DEBUG] Removing about {:.2f}% as outliers; factor = {:.1f}'.format(percent, fac))
      rem_indexes = copy.deepcopy(new_indexes)
      break
    else:
      if rem_indexes != []:
          break
      fac = fac + .1
      print('[DEBUG] Factor {:.1f} only removes about {:.2f}% as outliers; trying for larger factor'.format(fac, percent))
    it += 1
    assert( it < 100 )
  xn = [ x[i] for i in range(len(x)) if i not in rem_indexes ] 
  yn = [ y[i] for i in range(len(y)) if i not in rem_indexes ] 
  return xn,yn

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

def parse_orbit_eci(fn):
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
    with open(fn, 'r') as fin:
        for line in fin.readlines():
            if line[0:5] == '[ECI]': 
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

## Assuming
## %Y-%m-%d %H:%M:%S.%f x y z Vx Vy Vz
## 0        1           2     5       
def parse_reference(fn):
    dct = {}
    with open(fn, 'r') as fin:
        for line in fin.readlines():
            if line[0] != '#' and line[0:5] != '[ECI]':
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

def parse_reference_eci(fn):
    dct = {}
    with open(fn, 'r') as fin:
        for line in fin.readlines():
            if line[0:5] == '[ECI]':
              l = line[5:].split()
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

def colAsArray(dct,coly,fac=1e0,colx='t',site_list=[]):
    x=[];y=[];sy=[];
    ## if the requested column has a variance entry, it should be col_var
    varcol = coly + '_var'
    for epoch,entries in dct.items():
      ## note that entries is a list of dictionaries ...
      for entry in entries:
        if coly in entry:
          if site_list==[] or entry['site'] in site_list:
              x.append(epoch if colx=='t' else entry[colx])
              y.append(entry[coly])
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
                    passes.append([julian.from_jd(cmin, fmt='mjd'), julian.from_jd(cmax, fmt='mjd')])
                    cmin = 1e22
                    cmax = 1e-22
                if t < cmin: cmin = t
                if t > cmax: cmax = t
    if cmax != 1e-22: passes.append([julian.from_jd(cmin, fmt='mjd'), julian.from_jd(cmax, fmt='mjd')])
    return passes

def parse_acceleration(fn):
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
            if line[0:5] == "[ACC]":
                l = line.split()
                line_is_valid = False
                try:
                    t = julian.from_jd(float(l[1]), fmt='mjd')
                    line_is_valid = True
                except:
                    pass
                if line_is_valid:
                    if t in dct:
                      dct[t].append(strdct2dct(' '.join(l[2:])))
                    else:
                      dct[t] = [strdct2dct(' '.join(l[2:]))]
              
    return dct

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
    '-x',
    '--snx',
    metavar='SINEX_DPOD',
    dest='snx',
    default=None,
    required=False,
    help='Sinex file to read antenna types from')

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
  '--plot-acceleration',
  action='store_true',
  dest='plot_acc')

parser.add_argument(
  '--plot-periodogram',
  action='store_true',
  dest='plot_per')

parser.add_argument(
  '--use-eci',
  action='store_true',
  dest='use_eci')

parser.add_argument(
  '--plot-dynamics',
  action='store_true',
  dest='plot_dynamic_params')

def plot_residuals(fn, plot_regression=False, saveas=None, snx=None, plot_range=(-.20,.20), regression_range=(-.25,.25)):
  passes, revs, dct = parse_residuals(fn)
  Regression = plot_regression

  def reducexy(x,y,miny,maxy):
    xx = []; yy = []
    for i in range(len(y)):
      if y[i] >= miny and y[i] <= maxy:
        xx.append(x[i])
        yy.append(y[i])
    return np.array(xx),np.array(yy)
  
  print("Plotting residuals Vs elevation ...", end='')
  fig, ax = plt.subplots(1, 1, figsize=set_size(width_pt - width_pt/3))
  fac = 1e0
  t,y,_ = colAsArray(dct,'res',fac)
  el,y,_ = colAsArray(dct,'res',fac,'el')
  ax.set_ylim([plot_range[0], plot_range[1]])
  ax.set_xlim([0, 90])
  ax.scatter(el,y,s=1,color='black')
  ax.grid(True, 'both', 'x')
  if Regression:
    xx, yy = reducexy(el, y, regression_range[0], regression_range[1])
    b = np.average(yy)
    ax.axhline(b, color='r', linestyle='--', linewidth=1)
    ax.text(1,b+.01,"{:+.2f} [mm]".format(b*1e3),color='red')
  ax.set(xlabel=r'Elevation Angle ($^\circ$)', ylabel='Residual [m/s]')
  fig.suptitle('DORIS residuals at {:}\n'.format(min(t).strftime('%Y-%m-%d')), fontsize=16)
  if saveas:
      full_name = saveas + 'resVsEle.pdf'
      print('Saving figure to {:}'.format(full_name))
      plt.savefig(full_name, format='pdf', bbox_inches='tight')
  plt.show()
  print(" done")
  
  if snx:
      fig, ax = plt.subplots(1, 1, figsize=set_size(width_pt - width_pt/3))
      fac = 1e0
      t,y,_ = colAsArray(dct,'res',fac)
      antennae = get_beacon_type(snx, t[0])
      print('Number of beacons matched in snx: {:}'.format(len(antennae)))
      ## unique antenna types
      antenna_types = list(dict.fromkeys([ v for _,v in antennae.items() ]))
      colors = ['blue', 'red', 'green', 'black']
      ci = 0
      for at in antenna_types:
        print('plotting residuals for antenna type: [{:}] color={:}'.format(at,colors[ci]))
        el,y,_ = colAsArray(dct,'res',fac,'el',[b for b,t in antennae.items() if t==at])
        ax.scatter(el,y,s=2,facecolor=colors[ci], label=at, lw=0,alpha=.5)
        ci += 1
      ax.set_ylim([plot_range[0], plot_range[1]])
      ax.set_xlim([0, 90])
      ax.grid(True, 'both', 'x')
      ax.legend()
      #if Regression:
      #  xx, yy = reducexy(el, y, regression_range[0], regression_range[1])
      #  b = np.average(yy)
      #  ax.axhline(b, color='r', linestyle='--', linewidth=1)
      #  ax.text(1,b+.01,"{:+.2f} [mm]".format(b*1e3),color='red')
      ax.set(xlabel=r'Elevation Angle ($^\circ$)', ylabel='Residual [m/s]')
      fig.suptitle('DORIS residuals at {:}\n'.format(min(t).strftime('%Y-%m-%d')), fontsize=16)
      if saveas:
          full_name = saveas + 'resVsEleVsAntenna.pdf'
          print('Saving figure to {:}'.format(full_name))
          plt.savefig(full_name, format='pdf', bbox_inches='tight')
      plt.show()
  
  ## Residuals Vs Time
  print("Plotting residuals Vs time ...", end='')
  fig, ax = plt.subplots(1, 1, figsize=set_size(width_pt - width_pt/3))
  fac = 1e0
  t,y,_ = colAsArray(dct,'res',fac)
  #t,y = remove_outliers(t,y,5)
  ax.set_ylim([plot_range[0], plot_range[1]])
  ax.scatter(t,y,s=1,color='black')
  if revs != []:
    for r in revs:
      ax.axvline(x = r, color = 'green', label = 'Rev')
  ax.xaxis.set_major_formatter(mdates.DateFormatter("%H:%M"))
  ax.xaxis.set_major_locator(mdates.HourLocator(interval=4))
  ax.xaxis.set_minor_locator(mdates.HourLocator(interval=2))
  ax.grid(True, 'both', 'x')
  ax.set(xlabel='Time', ylabel='Residual [m/s]')
  fig.suptitle('DORIS residuals at {:}\n'.format(min(t).strftime('%Y-%m-%d')), fontsize=16)
  if Regression:
    # obtain m (slope) and b(intercept) of linear regression line
    xx, yy = reducexy(t,y,regression_range[0], regression_range[1])
    #x = mdates.date2num(xx)
    #m, b = np.polyfit(x, yy, 1)
    #print(m,b)
    # add linear regression line to scatterplot
    #ax.plot(xx, m*x+b, color='red', linewidth=1)
    b = np.average(yy)
    ax.axhline(b, color='r', linestyle='--', linewidth=1)
    ax.text(t[0],b+.01,"{:+.2f} [mm]".format(b*1e3),color='red')
  fig.autofmt_xdate()
  if saveas:
      full_name = saveas + 'resVsTime.pdf'
      print('Saving figure to {:}'.format(full_name))
      plt.savefig(full_name, format='pdf', bbox_inches='tight')
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

def obs_per_site(fn, saveas=None):
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
  fig, ax = plt.subplots(1, 1, figsize=set_size(width_pt/3))
  ax.bar(range(len(sites)), nobs)
  ax.set_xticklabels(sites)
  plt.xticks(range(len(mres)))
  plt.xticks(rotation=90, ha='right')
  ax.grid(True, 'both', 'y')
  ax.set(xlabel='Beacon/Site', ylabel='Num. Obs', title='Observations per Beacon')
  #fig.suptitle('DORIS residuals @ {:}\n'.format(min(t).strftime('%Y-%m-%d')), fontsize=16)
  if saveas:
      full_name = saveas + 'obsBeacons.pdf'
      print('Saving figure to {:}'.format(full_name))
      plt.savefig(full_name, format='pdf', bbox_inches='tight')
  plt.show()

def plot_site(fn, site):
  saa_passes, revs, dct = parse_residuals(fn)
  dct = filter_dict_site(dct, site)
  fig, ax = plt.subplots(2,1)
  fac = 1e0

  t,y,err = colAsArray(dct,'res',fac)
  t1,ypred,_= colAsArray(dct,'res_prediction',fac)
  for i in range(len(t1)):
      print("{:} {:.3f} {:.3f} {:.3f} {:.3f}".format(t[i], y[i],ypred[i],y[i]-ypred[i], err[i]))
  if revs != []:
    for r in revs:
      ax[0].axvline(x=r,color='green',label='Rev')
  #ax[0].errorbar(t,y,yerr=err,fmt='none', ecolor='red', elinewidth=.1, alpha=.1)
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

  fig, ax = plt.subplots(2, 1, sharey='col', figsize=set_size(width_pt, subplots=(2,1)))
  fac = 1e0
  
  t,y,err = colAsArray(dct,'cd',fac)
  ax[0].errorbar(t,y,yerr=err,fmt='none', ecolor='red', elinewidth=.1, alpha=.2, capsize=1)
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
  if saveas:
      full_name = saveas + 'dynamicParams.pdf'
      plt.savefig(full_name, format='pdf', bbox_inches='tight')
  plt.show()

def plot_state_diffs(fnref, fntest, plot_regression=False, plot_details=False, saveas=None):
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

  Regression = plot_regression
  Details = plot_details

  fig, axs = plt.subplots(3, 2, sharey='col', figsize=set_size(width_pt, subplots=(3,2)))
  ## x-, y-, z-components ...
  for component in range(3):
      ## first position, then velocity ...
      for pv in range(2):
          fac = 1e0
          key = whichCol(component, pv)
          if Details and revs != []:
              for r in revs: axs[component, pv].axvline(x=r,color='green',label='Rev')
          t,y,err = colAsArray(diffs,key,fac)
          axs[component, pv].scatter(t,y,s=.3,color='black')
          axs[component, pv].set_title(whichTitle(component, pv, fac))
          axs[component, pv].xaxis.set_major_formatter(mdates.DateFormatter("%H:%M"))
          axs[component, pv].xaxis.set_major_locator(mdates.HourLocator(interval=3))
          axs[component, pv].xaxis.set_minor_locator(mdates.HourLocator(interval=1))
          if pv == 0: 
              axs[component, pv].yaxis.set_major_locator(MultipleLocator(1))
              axs[component, pv].yaxis.set_minor_locator(MultipleLocator(1))
              axs[component, pv].set_ylim(-3.6,3.6)
          else: 
              axs[component, pv].yaxis.set_major_locator(MultipleLocator(.001))
              axs[component, pv].yaxis.set_minor_locator(MultipleLocator(.0005))
          if Regression:
            # obtain m (slope) and b(intercept) of linear regression line
            x = mdates.date2num(t)
            m, b = np.polyfit(x, y, 1)
            # add linear regression line to scatterplot
            axs[component, pv].plot(t, m*x+b)
          dsts = ColStatistics(None,y,fac)
          if not Details:
              axs[component, pv].grid(True, 'both', 'x', color='grey', linewidth=.1)
              axs[component, pv].grid(True, 'both', 'y', color='green', linewidth=.2)
          if Details and shd_passes != []:
              for intrv in shd_passes: axs[component, pv].axvspan(intrv[0], intrv[1], alpha=0.2, color='red')
          print('Component {:}/{:}: mean {:+.6f} +/- {:.9f} Max: {:.6f} Min: {:.6f}'.format(component, pv, dsts.mean, math.sqrt(dsts.variance), dsts.minmax[1], dsts.minmax[0]))
  ## Rotate date labels automatically
  fig.autofmt_xdate()
  fig.suptitle('State Differences Estimated - Reference\n{:}'.format(t[0].strftime('%Y-%m-%d')), fontsize=16)
  if saveas:
      full_name = saveas + 'statediffs.pdf'
      print('Saving figure to {:}'.format(full_name))
      plt.savefig(full_name, format='pdf', bbox_inches='tight')
  plt.show()

def plot_state_diffs_eci(fnref, fntest, plot_regression=False, plot_details=False, saveas=None):
  def whichCol(component, posvel):
      labels = ['x','y','z','vx','vy','vz']
      return labels[component + int(posvel==1)*3]
  def whichTitle(component, posvel, fac=1e0):
      labels = [r'$x$',r'$y$',r'$z$',r'$v_x$',r'$v_y$',r'$v_z$']
      title = r'{:} [m]'.format(labels[component])
      if posvel!=0:title = r'{:} [m]'.format(labels[component+3])
      return title
  def whichTitleRTN(component, posvel, fac=1e0):
      title = ' [{:}m]'.format('k' if fac==1e-3 else '')
      if posvel == 1:
          title = ' Velocity [{:}m/sec]'.format('k' if fac == 1e-3 else '')
      return ['Radial','In-Track','Cross-track'][component] + title
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
  def rotMatrix(pos,vel):
    p = np.array(pos)
    v = np.array(vel) 
    T = np.identity(3) 
    r1 = p / np.sqrt(np.sum(p**2))
    T[0] = r1 ## first row
    r2 = np.cross(p,v)
    T[1] = r2 / np.sqrt(np.sum(r2**2))
    r3 = np.cross(T[1], T[0])
    T[2] = r3 / np.sqrt(np.sum(r3**2))
    np.set_printoptions(precision=3)
    # print(T)
    return T
  def dcar2drac(dct1,dct2):
    resdct = {}
    for t,d1 in dct1.items():
      if t not in dct2:
          print('## Warning! failed to find record for date: {:}'.format(t))
      else:
          if len(dct2[t]) != 1:
            print('## Warning! more than one estimates found for block')
          d1 = d1[0]
          d2 = dct2[t][-1]
          refpos = np.array([d1['x'], d1['y'], d1['z']])
          refvel = np.array([d1['vx'], d1['vy'], d1['vz']])
          pos2 = np.array([d2['x'], d2['y'], d2['z']])
          rtn = rotMatrix(refpos,refvel).dot(refpos-pos2)
          resdct[t] = [{'x':rtn[0], 'y':rtn[1], 'z':rtn[2]}]
    return resdct

  dct1 = parse_reference_eci(fnref)
  dct2 = parse_orbit_eci(fntest)
  diffs = make_state_diffs(dct1,dct2)
  shd_passes = shadow_passes(fntest)

  Regression = plot_regression
  Details = plot_details

  fig, axs = plt.subplots(3, 2, sharey='col', figsize=set_size(width_pt, subplots=(3,2)))
  ## x-, y-, z-components ...
  for component in range(3):
      ## first position, then velocity ...
      for pv in range(2):
          fac = 1e0
          key = whichCol(component, pv)
          t,y,err = colAsArray(diffs,key,fac)
          axs[component, pv].scatter(t,y,s=.3,color='black')
          axs[component, pv].set_title(whichTitle(component, pv, fac))
          axs[component, pv].xaxis.set_major_formatter(mdates.DateFormatter("%H:%M"))
          axs[component, pv].xaxis.set_major_locator(mdates.HourLocator(interval=3))
          axs[component, pv].xaxis.set_minor_locator(mdates.HourLocator(interval=1))
          if pv == 0: 
              axs[component, pv].yaxis.set_major_locator(MultipleLocator(1))
              axs[component, pv].yaxis.set_minor_locator(MultipleLocator(1))
              axs[component, pv].set_ylim(-3.6,3.6)
          else: 
              axs[component, pv].yaxis.set_major_locator(MultipleLocator(.001))
              axs[component, pv].yaxis.set_minor_locator(MultipleLocator(.0005))
          if Regression:
            # obtain m (slope) and b(intercept) of linear regression line
            x = mdates.date2num(t)
            m, b = np.polyfit(x, y, 1)
            # add linear regression line to scatterplot
            if (abs(m)>.5):
                axs[component, pv].plot(t, m*x+b, color='red', linewidth=1)
                axs[component, pv].text(t[0], -3, 'slope: {:+.2f}'.format(m), color='red')
          dsts = ColStatistics(None,y,fac)
          if not Details:
              axs[component, pv].grid(True, 'both', 'x', color='grey', linewidth=.1)
              axs[component, pv].grid(True, 'both', 'y', color='green', linewidth=.2)
          if Details and shd_passes != []:
              for intrv in shd_passes: axs[component, pv].axvspan(intrv[0], intrv[1], alpha=0.2, color='red')
          print('Component {:}/{:}: mean {:+.6f} +/- {:.9f} Max: {:.6f} Min: {:.6f}'.format(component, pv, dsts.mean, math.sqrt(dsts.variance), dsts.minmax[1], dsts.minmax[0]))
  ## Rotate date labels automatically
  fig.autofmt_xdate()
  fig.suptitle('State Differences Estimated - Reference (ECI)\n{:}'.format(t[0].strftime('%Y-%m-%d')), fontsize=16)
  if saveas:
      full_name = saveas + 'stateDiffsEci.pdf'
      print('Saving figure to {:}'.format(full_name))
      plt.savefig(full_name, format='pdf', bbox_inches='tight')
  plt.show()
  
  diffs = dcar2drac(dct1,dct2)
  fig, axs = plt.subplots(3, 1, sharey='col', figsize=set_size(width_pt-width_pt/3, subplots=(3,1)))
  ## x-, y-, z-components ...
  for component in range(3):
      fac = 1e0
      key = whichCol(component, 0)
      t,y,_= colAsArray(diffs,key,fac)
      axs[component].scatter(t,y,s=.3,color='black')
      axs[component].set_title(whichTitleRTN(component, 0, fac))
      axs[component].xaxis.set_major_formatter(mdates.DateFormatter("%H:%M"))
      axs[component].xaxis.set_major_locator(mdates.HourLocator(interval=3))
      axs[component].xaxis.set_minor_locator(mdates.HourLocator(interval=1))
      axs[component].yaxis.set_major_locator(MultipleLocator(1))
      axs[component].yaxis.set_minor_locator(MultipleLocator(1))
      axs[component].set_ylim(-4,4)
      #if Regression:
      #  # obtain m (slope) and b(intercept) of linear regression line
      #  x = mdates.date2num(t)
      #  m, b = np.polyfit(x, y, 1)
      #  # add linear regression line to scatterplot
      #  axs[component].plot(t, m*x+b, color='red', linewidth=1)
      #  if (abs(m)>.5):
      #      axs[component].text(t[0], -3, 'slope: {:+.2f}'.format(m), color='red')
      dsts = ColStatistics(None,y,fac)
      print('Component {:}/{:}: mean {:+.6f} +/- {:.9f} Max: {:.6f} Min: {:.6f}'.format(component, 0, dsts.mean, math.sqrt(dsts.variance), dsts.minmax[1], dsts.minmax[0]))
  fig.autofmt_xdate()
  fig.suptitle('State Differences Estimated - Reference (RTN)\n{:}'.format(t[0].strftime('%Y-%m-%d')), fontsize=16)
  if saveas:
      full_name = saveas + 'stateDiffsRtn.pdf'
      print('Saving figure to {:}'.format(full_name))
      plt.savefig(full_name, format='pdf', bbox_inches='tight')
  plt.show()

def plot_state(fntest, plot_regression=False, plot_details=False, saveas=None):
  def whichCol(component, posvel):
      labels = ['x','y','z','vx','vy','vz']
      return labels[component + int(posvel==1)*3]
  def whichTitle(component, posvel, fac=1e0):
      title = ' [{:}m]'.format('k' if fac==1e-3 else '')
      if posvel == 1:
          title = ' Velocity [{:}m/sec]'.format('k' if fac == 1e-3 else '')
      return ['X','Y','Z'][component] + title

  revs, dct2 = parse_orbit(fntest)
  diffs = dct2
  shd_passes = shadow_passes(fntest)

  Regression = plot_regression
  Details = plot_details

  fig, axs = plt.subplots(3, 2, sharey='col', figsize=set_size(width_pt, subplots=(3,2)))
  ## x-, y-, z-components ...
  for component in range(3):
      ## first position, then velocity ...
      for pv in range(2):
          fac = 1e0
          key = whichCol(component, pv)
          if Details and revs != []:
              for r in revs: axs[component, pv].axvline(x=r,color='green',label='Rev')
          t,y,err = colAsArray(diffs,key,fac)
          axs[component, pv].scatter(t,y,s=.3,color='black')
          axs[component, pv].set_title(whichTitle(component, pv, fac))
          axs[component, pv].xaxis.set_major_formatter(mdates.DateFormatter("%H:%M"))
          axs[component, pv].xaxis.set_major_locator(mdates.HourLocator(interval=3))
          axs[component, pv].xaxis.set_minor_locator(mdates.HourLocator(interval=1))
          if Regression:
            # obtain m (slope) and b(intercept) of linear regression line
            x = mdates.date2num(t)
            m, b = np.polyfit(x, y, 1)
            # add linear regression line to scatterplot
            axs[component, pv].plot(t, m*x+b)
          dsts = ColStatistics(None,y,fac)
          if not Details:
              axs[component, pv].grid(True, 'both', 'x', color='grey', linewidth=.1)
              axs[component, pv].grid(True, 'both', 'y', color='green', linewidth=.2)
          if Details and shd_passes != []:
              for intrv in shd_passes: axs[component, pv].axvspan(intrv[0], intrv[1], alpha=0.2, color='red')
  ## Rotate date labels automatically
  fig.autofmt_xdate()
  fig.suptitle('Estimated ECEF state at {:}'.format(t[0].strftime('%Y-%m-%d')), fontsize=16)
  if saveas:
      full_name = saveas + 'state.pdf'
      print('Saving figure to {:}'.format(full_name))
      plt.savefig(full_name, format='pdf', bbox_inches='tight')
  plt.show()

def plot_accelerations(fn, saveas=None):
    def filterbyt(t,y,mindt):
        tt = [t[0]]
        yy = [y[0]]
        for i in range(len(t)):
            if (t[i] - tt[-1]).total_seconds() < mindt:
                yy[-1] = (yy[-1] + y[i])/2.
            else:
                tt.append(t[i])
                yy.append(y[i])
        return tt,yy

    tstart = datetime.datetime(2050,1,1)
    dct = parse_acceleration(fn)
    fig, ax = plt.subplots(1, 1, figsize=set_size(width_pt))
    ax.set_yscale('log')
    ax.set_ylabel(r'Accelerations [$m/s^2$]')
    ax.xaxis.set_major_formatter(mdates.DateFormatter("%H:%M"))
    ax.xaxis.set_major_locator(mdates.HourLocator(interval=3))
    ax.xaxis.set_minor_locator(mdates.HourLocator(interval=1))
    colors = ['tab:blue', 'tab:orange', 'tab:green', 'tab:red', 'tab:purple', 'tab:brown', 'tab:pink', 'tab:gray', 'tab:olive', 'tab:cyan']
    labels = [r'gravity $n>1$', 'Moon', 'Sun', 'earth tide', 'ocean tides', 'srp', 'drag']
    ci = 0
    ctlist = []
    
    for key in ['gm', 'moon', 'sun', 'setide', 'octide', 'srp', 'drag']:
    #for key in ['gm', 'moon', 'sun', 'setide', 'octide', 'srp']:
        dtsec = 60. if key != 'srp' else 1.
        t,y1,_ = colAsArray(dct,key,1e0)
        print('Original data size {:}'.format(len(y1)), end='')
        t1,y1 = filterbyt(t,y1,dtsec)
        print('reduced to {:}'.format(len(y1)))
        if key == 'gm':
            t,y,_ = colAsArray(dct,'grav',1e0)
            t,y = filterbyt(t,y,dtsec)
            y2 = [ abs(y[0] - yy[1]) for yy in zip(y,y1) ]
            ax.scatter(t,y2,s=1,color=colors[ci], label=labels[ci])
            ctlist.append(colors[ci])
            ax.scatter(t,y1,s=1,color='black',label=r'gravity $n=0$')
            ctlist.append('black')
            if t[0] < tstart: tstart = t[0]
        else:
            ax.scatter(t1,y1,s=1,color=colors[ci],label=labels[ci])
            ctlist.append(colors[ci])
        ci += 1

    #ax.set_ylim([1e-14,1e2])
    
    fig.autofmt_xdate()
    ax.tick_params(axis="x",direction="in", pad=-15)
    for tick in ax.xaxis.get_majorticklabels(): tick.set_horizontalalignment("left")
    # Shrink current axis's height by 10% on the bottom
    box = ax.get_position()
    ax.set_position([box.x0, box.y0 + box.height * 0.1,
                 box.width, box.height * 0.9])
    # Put a legend below current axis
    leg = ax.legend(loc='upper center', bbox_to_anchor=(0.5, -0.005), fancybox=True, shadow=True, ncol=4)
    for h, t in zip(leg.legendHandles, leg.get_texts()): t.set_color(h.get_facecolor()[0])
    for h in leg.legendHandles: h._sizes = [30]
        
    fig.suptitle('Computed accelerations at {:}'.format(tstart.strftime('%Y-%m-%d')), fontsize=16)
    if saveas:
        full_name = saveas + 'accelerations.pdf'
        print('Saving figure to {:}'.format(full_name))
        plt.savefig(full_name, format='pdf', bbox_inches='tight')
    plt.show()

def plot_state_diff_periodograms(fnref, fntest, use_eci=False, saveas=None):
  def omega2period(w): return 2. * math.pi / w
  def period2omega(T): return 2. * math.pi / T
  def whichCol(component, posvel):
      labels = ['x','y','z','vx','vy','vz']
      return labels[component + int(posvel==1)*3]
  def whichTitle(component, posvel, fac=1e0):
    r = [r'$x$', r'$y$', r'$z$']
    v = [r'$v_x$', r'$v_y$', r'$v_z$']
    return r'Component {:}'.format(r[component] if posvel==0 else v[component])
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
  def filter_peaks(peaks, w, pgram):
    wvals = [ w[i] for i in peaks ]
    pvals = [ pgram[i] for i in peaks ]
    pi = [ x for x,_ in reversed(sorted(zip(peaks,pvals), key=lambda pair: pair[1])) ]
    return pi[0:3] if len(pi)>=3 else pi
  def remove_trend(t,y):
    x = mdates.date2num(t)
    m, b = np.polyfit(x, y, 1)
    yy = []
    for i in range(len(y)):
        yy.append(y[i]-(m*x[i]+b))
    return yy

  if use_eci:
      dct1 = parse_reference_eci(fnref)
      dct2 = parse_orbit_eci(fntest)
  else:
      dct1 = parse_reference(fnref)
      revs, dct2 = parse_orbit(fntest)
  diffs = make_state_diffs(dct1,dct2)
  
  nout = 50000
  Tmin =  6. * 3600    # 6 hours [sec]
  Tmax =  .5 * 3600    # .1 of hour [sec]
  Wmax = period2omega(Tmax)
  Wmin = period2omega(Tmin)
  assert(Wmin < Wmax)
  w = np.linspace(Wmin, Wmax, nout)
  
  fig, axs = plt.subplots(3, 2, sharex=True, sharey=True, figsize=set_size(width_pt, subplots=(3,2)))
  ## x-, y-, z-components ...
  for component in range(3):
      ## first position, then velocity ...
      for pv in range(2):
          fac = 1e0
          key = whichCol(component, pv)
          t,y,_ = colAsArray(diffs,key,fac)
          y = remove_trend(t,y)
          sec = []
          t0 = t[0]
          for i in range(len(t)): sec.append((t[i]-t0).total_seconds())
          pgram = signal.lombscargle(sec, y, w, normalize=True)
          axs[component, pv].plot(w, pgram)
          if component == 2: axs[component, pv].set_xlabel('Angular frequency [rad/s]')
          if component == 1: axs[component, pv].set_ylabel('Normalized amplitude')
          peaks = signal.find_peaks(pgram, height=.005)[0] ## peak indexes
          peaks = filter_peaks(peaks, w, pgram)
          marks = [ w[i] for i in peaks ]
          axs[component, pv].vlines(x=marks, ymin=0, ymax=1, color = 'r', linestyle='--', linewidth=1)
          textystart = .4
          for i, x in enumerate(marks):
            axs[component, pv].text(x, textystart, "{:.1f}[hr]".format(omega2period(x)/3600), rotation=0, verticalalignment='center', fontsize=10)
            textystart += .2
          axs[component, pv].set_title(whichTitle(component, pv))
  fig.suptitle('Estimated - Reference at {:}\nLomb-Scargle Periodogram'.format(t[0].strftime('%Y-%m-%d')), fontsize=12)
  if saveas:
      full_name = saveas + 'periodogram.pdf'
      print('Saving figure to {:}'.format(full_name))
      plt.savefig(full_name, format='pdf', bbox_inches='tight')
  plt.show()

if __name__ == '__main__':

    args = parser.parse_args()

    if args.plot_acc:
        plot_accelerations(args.input, args.save_as)

    if args.ref_state:
      if args.use_eci:
        plot_state_diffs_eci(args.ref_state, args.input, True, False, args.save_as)
      else:
        plot_state_diffs(args.ref_state, args.input, False, True, args.save_as)
      #plot_state(args.input, False, False, args.save_as)
      if args.plot_per:
        plot_state_diff_periodograms(args.ref_state, args.input, args.use_eci, args.save_as)

    if args.plot_res:
      plot_residuals(args.input, True, args.save_as, args.snx)

    if args.plot_dynamic_params:
      plot_dynamic_params(args.input, args.save_as)

    if args.sites:
        for s in args.sites:
            plot_site(args.input, s)

    if 1==2:
        obs_per_site(args.input, args.save_as)
