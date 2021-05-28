#! /usr/bin/python

import sys
import argparse
import numpy as np
from scipy.interpolate import interp1d
import matplotlib.pyplot as plt

def parse_antex(antex, antenna_type):
    with open(antex, 'r') as fin:
        line = fin.readline()
        if line[60:].strip() != 'ANTEX VERSION / SYST':
            print('[ERROR] Expected first line of ANTEX, found:\n\'{:}\''.format(line.strip()))
            raise RuntimeError('Error reading ANTEX file')
        dct = {}
        while line:
            if line[60:].strip() == 'TYPE / SERIAL NO':
                cantenna = line[0:20].strip()
                print('[D] Checking antenna type \'{:}\' against \'{:}\''.format(cantenna.lower(), antenna_type.lower()))
                if (cantenna.lower() == antenna_type.lower()):
                    line = fin.readline()
                    assert(line[60:].strip() == 'METH / BY / # / DATE')
                    line = fin.readline()
                    assert(line[60:].strip() == 'DAZI')
                    dazi = float(line[0:20].strip())
                    line = fin.readline()
                    assert(line[60:].strip() == 'ZEN1 / ZEN2 / DZEN')
                    zen1, zen2, dzen = float(line[2:8]), float(line[8:14]), float(line[14:20])
                    dct['zen1'] = zen1
                    dct['zen2'] = zen2
                    dct['dzen'] = dzen
                    line = fin.readline()
                    assert(line[60:].strip() == '# OF FREQUENCIES')
                    num_freqs = int(line[0:20])
                    line = fin.readline()
                    assert(line[60:].strip() == 'SINEX CODE')
                    cur_freq = 0
                    while cur_freq < num_freqs:
                        print('Collecting corrections for freq#{:}'.format(cur_freq))
                        line = fin.readline()
                        while line[60:].strip() == 'COMMENT': line = fin.readline()
                        assert(line[60:].strip() == 'START OF FREQUENCY')
                        freq = line[0:20].strip()
                        dct[freq]={}
                        line = fin.readline()
                        assert(line[60:].strip() == 'NORTH / EAST / UP')
                        n, e, u = float(line[0:10]), float(line[10:20]), float(line[20:30])
                        dct[freq]['neu'] = [n, e, u]
                        line = fin.readline()
                        assert(line[3:8] == 'NOAZI')
                        noazi = [ float(line[(i+1)*8:(i+2)*8]) for i in range(int((zen2 - zen1)//dzen)+1)]
                        dct[freq]['noazi'] = noazi
                        if dazi > 0:
                            az = 0
                            while az <= 360e0:
                                line = fin.readline()
                                az += dazi
                        line = fin.readline()
                        assert(line[60:].strip() == 'END OF FREQUENCY')
                        cur_freq+=1
                    return dct
            line = fin.readline()
        return None

class myFormatter(argparse.ArgumentDefaultsHelpFormatter,
                  argparse.RawTextHelpFormatter):
    pass

parser = argparse.ArgumentParser(
    formatter_class=myFormatter,
    description=
    'Plot DORIS phase center variations for an antenna, given an antex file.',
    epilog=('''National Technical University of Athens,
    Dionysos Satellite Observatory\n
    Send bug reports to:
    Xanthos Papanikolaou, xanthos@mail.ntua.gr
    Dimitris Anastasiou,danast@mail.ntua.gr
    May, 2021'''))

parser.add_argument(
    '-x',
    '--antex-file',
    metavar='ANTEX_FILE',
    dest='antex_file',
    required=True,
    help=
    'Provide an ANTEX file with records for the given antenna type.'
)
parser.add_argument(
    '-a',
    '--antenna_type',
    metavar='ANTENNA_TYPE',
    dest='antenna_type',
    required=False,
    default='alcatel',
    help=
    'Choose the antenna type.'
)

parser.add_argument(
    '-f',
    '--frequency',
    metavar='FREQUENCY',
    dest='frequency',
    required=False,
    default='D01',
    help=
    'Frequency to plot as recorded in the ANTEX file (e.g. \'D01\', \'R01\', etc)'
)

if __name__ == '__main__':

    args = parser.parse_args()
    dct = parse_antex(args.antex_file, args.antenna_type)

    if not dct:
        print('Failed to find antex block for antenna \'{:}\''.format(args.antenna_type))
        sys.exit(2)
    if args.frequency not in dct or 'noazi' not in dct[args.frequency]:
        print('Failed find antex block for frequency \'{:}\''.format(args.frequency))
        sys.exit(3)

    ## x-axis is zanith angle in degrees
    x_start = dct['zen1']
    x_stop = dct['zen2']
    x_step = dct['dzen']
    x_range_inclusive = True

    ## x array (based on x-axis)
    x_numpts = int(np.floor((x_stop - x_start) / x_step + 1))
    x, s = np.linspace(x_start, x_stop, num=x_numpts, endpoint=True, retstep=True)
    assert(int(s) == int(x_step))

    ## y array (phase law corrections on x-axis points)
    y = np.array(dct[args.frequency]['noazi'])

    ## validate equal sizes (x-y arrays)
    assert(x.size == y.size)

    ## interpolation types
    intrp_type = ['linear', 'nearest', 'nearest-up', 'zero', 'slinear', 'quadratic', 'cubic', 'previous', 'next']

    ## prepare plot
    fig = plt.figure()
    gs = fig.add_gridspec(3, 3, hspace=0e0, wspace=0e0)
    ax = gs.subplots(sharex='col', sharey='row')
    fig.suptitle('Antex NOAZI records for {:}::{:}'.format(args.antenna_type, args.frequency))

    for i, tp in enumerate(intrp_type):
        yi = i//3
        xi = i - yi*3

        ## interpolate
        f = interp1d(x, y, intrp_type[i], assume_sorted=True)

        xnew = np.linspace(x_start+1, x_stop, num=x_numpts*20, endpoint=True, retstep=False)
        ynew = f(xnew)

        ax[xi, yi].plot(x, y, 'o', xnew, ynew, '-')
        ax[xi, yi].title.set_text(intrp_type[i])

    for ax in ax.flat: ax.label_outer()
    fig.show()
    fig.savefig("foo.pdf", bbox_inches='tight')
