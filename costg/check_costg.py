#! /usr/bin/python

import os, sys, datetime
import argparse
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import argparse
from scipy.interpolate import interp1d
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.ticker import MultipleLocator
import fileinput
from scipy import stats
import julian
import subprocess
import math

plt.style.use('seaborn-v0_8-pastel')
tex_fonts = {
    "text.usetex": True,
    "font.family": "serif",
    "axes.labelsize": 12,
    "font.size": 12,
    "legend.fontsize": 10,
    "xtick.labelsize": 10,
    "ytick.labelsize": 10
}
plt.rcParams.update(tex_fonts)

def unitsstr(scale): return r'$[m/s^2]\times{:.1e}$'.format(1e0/scale) if scale != 1e0 else r'$[m/s^2]$'
def handle_date(mjd): return julian.from_jd(mjd, fmt='mjd')
def get_diffs(dct, component,raw=False):
    t = []; y= [];
    k1,k2 = list(zip(['rax','ray','raz'],['ax','ay','az']))[component]
    for entry in dct:
        t.append(entry['mjd'])
        y.append(entry[k1] - entry[k2]) if not raw else y.append(entry[k2])
    return t,y
def get_diffs_norm(dct,raw=False):
    t = []; y= [];
    for entry in dct:
        t.append(entry['mjd'])
        n = 0e0
        for component in range(3):
            k1,k2 = list(zip(['rax','ray','raz'],['ax','ay','az']))[component]
            if not raw:
                n += (entry[k1] - entry[k2])*(entry[k1] - entry[k2])
            else:
                n += entry[k2]*entry[k2]
        y.append(math.sqrt(n))
    return t,y
def parse(fn, force=None):
    dct = []
    with open(fn, 'r') as fin:
        for line in fin.readlines():
            ## parse file/input
            if line[0] != '#':
                l = line.split()
                if not force or line.startswith(force):
                    ofset = 1 if force else 0
                    try:
                        dct.append({'mjd': handle_date(float(l[0+ofset])), 'rax':float(l[1+ofset]), 'ray':float(l[2+ofset]), 
                            'raz':float(l[3+ofset]), 'ax':float(l[4+ofset]), 'ay':float(l[5+ofset]), 'az':float(l[6+ofset])})
                    except:
                        #print("Failed to parse line: [{:}]; skipped!".format(line.strip()), 
                        #    file=sys.stderr)
                        pass
    return dct
def plot2(fn, scale, title, saveas, forceid=[''], raw=False):
    for f in forceid:
        print('|{:}|'.format(' '.join([title,f])))
        forcestr = None if f==' ' else f
        dct = parse(fn, forcestr)
        fig, axs = plt.subplots(4, 1, figsize=(10, 6), sharex=True, sharey=True, constrained_layout=True)
        # cartesian components ...
        for i in range(3):
            t,y = get_diffs(dct, i, raw)
            y = [yy * scale for yy in y]
            sts = stats.describe(y)
            axs[i].scatter(t,y,s=1,color='black')
            axs[i].axhline(sts.mean, color='r', linestyle='--', linewidth=1)
            axs[i].text(t[0],sts.minmax[0],r'{:+.1e} $\pm$ {:+.3e}'.format(sts.mean, math.sqrt(sts.variance)),color='red')
            axs[i].text(t[0],sts.minmax[1],r'min:{:+.1e} max:{:+.1e}'.format(sts.minmax[0], sts.minmax[1]),color='red')
            axs[i].set_ylabel([r'$\delta\ddot{x}$', r'$\delta\ddot{y}$', r'$\delta\ddot{z}$'][i]+unitsstr(scale))
            axs[i].xaxis.set_major_formatter(mdates.DateFormatter("%H:%M"))
            axs[i].xaxis.set_minor_locator(mdates.HourLocator())
            # print statistics in MarkDown format (as a single line)
            print('|{:}|{:+.3e}|{:+.3e}|{:+.3e}|{:+.3e}|{:+.3e}|'.format(['x','y','z'][i], sts.mean, math.sqrt(sts.variance), sts.minmax[0], sts.minmax[1], scale))
        # plot norm ...
        t,y = get_diffs_norm(dct,raw)
        y = [yy * scale for yy in y]
        sts = stats.describe(y)
        axs[3].scatter(t,y,s=1,color='black')
        axs[3].axhline(sts.mean, color='r', linestyle='--', linewidth=1)
        axs[3].text(t[0],-sts.minmax[1],r'{:+.1e} $\pm$ {:+.3e}'.format(sts.mean, math.sqrt(sts.variance)),color='red')
        axs[3].text(t[0],sts.minmax[1],r'min:{:+.1e} max:{:+.1e}'.format(sts.minmax[0], sts.minmax[1]),color='red')
        axs[3].set_ylabel(r'$\delta\ddot{a}$'+unitsstr(scale))
        axs[3].xaxis.set_major_formatter(mdates.DateFormatter("%H:%M"))
        axs[3].xaxis.set_minor_locator(mdates.HourLocator())
        print('|{:}|{:+.3e}|{:+.3e}|{:+.3e}|{:+.3e}|{:+.3e}|'.format('a', sts.mean, math.sqrt(sts.variance), sts.minmax[0], sts.minmax[1], scale))
        # set x-axis lable ...
        axs[-1].set_xlabel("Date {:}".format(t[0].strftime("%Y/%m/%d")))
        ## title and save ...
        line1_title = 'COST-G Benchmark Diffs' if not raw else 'Acceleration Reults'
        fig.suptitle('{:}\n{:}'.format(line1_title, title+('' if not forcestr else ' '+forcestr)))
        savefn = saveas if not forcestr else saveas.replace('.png','-'+forcestr)+'.png'
        plt.savefig(savefn)
        # plt.show()

def plot(fn, scale, title, saveas, forceid=[''], raw=False):
    for f in forceid:
        print('|{:}|'.format(' '.join([title,f])))
        forcestr = None if f==' ' else f
        dct = parse(fn, forcestr)
        fig, axs = plt.subplots(3, 1, figsize=(10, 6), sharex=True, sharey=True, constrained_layout=True)
        # cartesian components ...
        for i in range(3):
            t,y = get_diffs(dct, i, raw)
            t=t[1:-1]
            y=y[1:-1]
            y = [yy * scale for yy in y]
            sts = stats.describe(y)
            axs[i].scatter(t,y,s=1,color='black')
            axs[i].axhline(sts.mean, color='r', linestyle='--', linewidth=1)
            axs[i].set_ylabel([r'$\delta\ddot{x} [m/s^2]$', r'$\delta\ddot{y} [m/s^2]$', r'$\delta\ddot{z} [m/s^2]$'][i])
            axs[i].xaxis.set_major_formatter(mdates.DateFormatter("%H:%M"))
            axs[i].xaxis.set_minor_locator(mdates.HourLocator())
            # print statistics in MarkDown format (as a single line)
            print('|{:}|{:+.3e}|{:+.3e}|{:+.3e}|{:+.3e}|{:+.3e}|'.format(['x','y','z'][i], sts.mean, math.sqrt(sts.variance), sts.minmax[0], sts.minmax[1], scale))
        # set x-axis lable ...
        axs[-1].set_xlabel("Date {:}".format(t[0].strftime("%Y/%m/%d")))
        ## title and save ...
        line1_title = 'COST-G Benchmark Diffs' if not raw else 'Acceleration Reults'
        fig.suptitle('{:}\n{:}, Scale {:.1e}'.format(line1_title, title+('' if not forcestr else ' '+forcestr), scale))
        savefn = saveas if not forcestr else saveas.replace('.png','-'+forcestr)+'.png'
        plt.savefig(savefn)
        # plt.show()

class myFormatter(argparse.ArgumentDefaultsHelpFormatter,
                  argparse.RawTextHelpFormatter):
    pass

parser = argparse.ArgumentParser(
    formatter_class=myFormatter,
    description=
    'Run validatation programs against COST-G benchmark tests',
    epilog=('''National Technical University of Athens,
    Dionysos Satellite Observatory\n
    Send bug reports to:
    Xanthos Papanikolaou, xanthos@mail.ntua.gr
    May, 2023'''))

parser.add_argument(
    '--costg-dir',
    metavar='COSTG_DIR',
    dest='costg_dir',
    default='../costG',
    required=False,
    help='Top-level directory of the COSTG benchamrk test.')

parser.add_argument(
    '--progs-dir',
    metavar='PROGS_DIR',
    dest='progs_dir',
    default=os.path.abspath(os.getcwd()),
    required=False,
    help='Directory with COSTG test executables')

parser.add_argument(
    '--plots-dir',
    metavar='PLOTS_DIR',
    dest='plots_dir',
    required=False,
    default='figures',
    help='Directory to save plots at.')

parser.add_argument(
    '--jpl-dir',
    metavar='JPL_DIR',
    dest='jpl_dir',
    required=False,
    default='../data/jpl',
    help='Directory to find jpl products, e.g. DE ephemerides.')

parser.add_argument(
    '--factor',
    metavar='FACTOR',
    dest='factor',
    type=float,
    required=False,
    default=1e12,
    help='Factor to be used when plotting, e.g. 1e9')

parser.add_argument(
    '--verbose',
    action='store_true',
    dest='verbose',
    help='Verbose mode on')

parser.add_argument(
    '--raw',
    action='store_true',
    dest='plot_raw',
    help='Plot \'raw\' results, i.e. not differences but computed acceleration results')

if __name__ == '__main__':

    ## parse cmd
    args = parser.parse_args()

    verboseprint = print if args.verbose else lambda *a, **k: None

    # executables including path
    nopath_prog_list = ['check-gravity-field.out', 'check-third-body.out', 'check-solid-earth-tide.out', 'check-pole-tide.out', 'check-ocean-pole-tide.out', 'check-relativistic.out', 'check-fes14b-ocean-tide.out', 'check-aod1b-atmosphericTides.out', 'check-aod1b-atmosphericTides-s1.out']
    prog_list = [ os.path.join(args.progs_dir, x) for x in nopath_prog_list ]

    # plot file names (output)
    ext = '.png'
    plot_fns = [ os.path.join(args.plots_dir, x.replace('check-','').replace('.out', ext)) for x in nopath_prog_list ]
    plot_titles = [ 'Gravity Field n=(2..180)', 'Third Body', 'Solid Earth Tide', 'Solid Earth Pole Tide', 'Ocean Pole Tide n=(2..180)', 'Relativistic Correction', 'FES2014 Major Waves n=(2..180)', 'AOD1B RL06 all waves n=(2..180)', 'AOD1B RL06 S1 wave n=(2..180)']

    # path to costg models
    mdlp = os.path.join(args.costg_dir, 'models')
    stlp = os.path.join(args.costg_dir, 'satellite')
    pltp = os.path.join(args.progs_dir, 'ploters')
    outp = os.path.join(args.plots_dir, '')
    jplp = os.path.join(args.jpl_dir, '')
    for p in [mdlp, stlp, pltp, outp, jplp]:
        if not os.path.isdir(p):
            print('ERROR. Cannot find directory {:}'.format(p), file=sys.stderr)
            sys.exit(1)
    
    # args list
    args_list = [
        { 'args': [os.path.join(mdlp, 'EIGEN6-C4.gfc'), os.path.join(stlp, '02gravityfield_itrf.txt'), os.path.join(stlp, '00orbit_itrf.txt')] },
        { 'args': [os.path.join(stlp, '03directTideMoon_icrf.txt'), os.path.join(stlp, '03directTideSun_icrf.txt'), os.path.join(stlp, '00orbit_icrf.txt'), os.path.join(jplp,'de421.bsp'), os.path.join(jplp, 'naif0012.tls')] },
        { 'args': [os.path.join(stlp, '04solidEarthTide_icrf.txt'), os.path.join(mdlp, 'eopc04_14_IAU2000.62-now'), os.path.join(stlp, '00orbit_icrf.txt'), os.path.join(jplp,'de421.bsp'), os.path.join(jplp, 'naif0012.tls'), os.path.join(stlp, '01earthRotation_quaternion.txt') ] },
        { 'args': [os.path.join(stlp, '05poleTide_icrf.txt'), os.path.join(mdlp, 'eopc04_14_IAU2000.62-now'), os.path.join(stlp, '00orbit_icrf.txt'), os.path.join(stlp, '01earthRotation_quaternion.txt') ] },
        { 'args': [os.path.join(stlp, '06oceanPoleTide_icrf.txt'), os.path.join(mdlp, 'eopc04_14_IAU2000.62-now'), os.path.join(stlp, '00orbit_icrf.txt'), os.path.join(stlp, '01earthRotation_quaternion.txt'), os.path.join(mdlp, 'desaiscopolecoef.txt') ] },
        { 'args': [os.path.join(stlp, '07relativistic_icrf.txt'), os.path.join(stlp, '00orbit_icrf.txt'), os.path.join(jplp,'de421.bsp'), os.path.join(jplp, 'naif0012.tls') ] },
        { 'args': [os.path.join(stlp, '11oceanTide_fes2014b_34major_icrf.txt'), os.path.join(mdlp, 'eopc04_14_IAU2000.62-now'), os.path.join(stlp, '00orbit_icrf.txt'), os.path.join(mdlp, 'FES2014b_oceanTide/oceanTide_FES2014b.potential.iers.txt'), os.path.join(stlp, '01earthRotation_quaternion.txt') ] },
        { 'args': [os.path.join(stlp, '09aod1b_atmosphericTides_icrf.txt'), os.path.join(mdlp, 'eopc04_14_IAU2000.62-now'), os.path.join(stlp, '00orbit_icrf.txt'), os.path.join(mdlp, 'AOD1B_tides_groops_gfc/atmosTides_AOD1BRL06.potential.iers.txt'), os.path.join(stlp, '01earthRotation_quaternion.txt') ] },
        { 'args': [os.path.join(stlp, '09aod1b_atmosphericTides_S1_icrf.txt'), os.path.join(mdlp, 'eopc04_14_IAU2000.62-now'), os.path.join(stlp, '00orbit_icrf.txt'), os.path.join(mdlp, 'AOD1B_tides_groops_gfc/atmosTides_AOD1BRL06.potential.iers.txt'), os.path.join(stlp, '01earthRotation_quaternion.txt') ] }
    ]
    for d in args_list:
        files = [ x for x in d['args'] ]
        for f in files:
            if not os.path.isfile(f):
                print('ERROR Failed to located input file {:}'.format(f), file=sys.stderr)
                sys.exit(1)

# go ahead and run executables
temp_fn = ".tmp"
for i,prog in enumerate(prog_list):
    cmd = [prog] + [ x for x in args_list[i]['args'] ]
    verboseprint('Running command: [{:}]'.format(' '.join(cmd)))
    forces = ['']
    ftmp = open(temp_fn, "w")
    try:
        result = subprocess.run(cmd, stdout=ftmp, check=True)
    except:
        print('ERROR Command failed! Aborting test suite ...')
        sys.exit(1)
    if os.path.basename(prog) == 'check-third-body.out':
        forces = ['[SUN]', '[MOON]']
    plot(temp_fn, args.factor, plot_titles[i], plot_fns[i], forces, args.plot_raw)
    ftmp.close()
    os.remove(".tmp")
