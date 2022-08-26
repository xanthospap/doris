#! /usr/bin/python3

import datetime
import matplotlib
import argparse
import sys, os, re

class myFormatter(argparse.ArgumentDefaultsHelpFormatter,
                  argparse.RawTextHelpFormatter):
    pass

parser = argparse.ArgumentParser(
    formatter_class=myFormatter,
    description=
    'Scatter Plot using matplotlib',
    epilog=('''National Technical University of Athens,
    Dionysos Satellite Observatory\n
    Send bug reports to:
    Xanthos Papanikolaou, xanthos@mail.ntua.gr
    Jul, 2022'''))

parser.add_argument(
    '-f',
    '--files',
    metavar='INPUT_FILES',
    dest='infiles',
    nargs = '+',
    default=None,
    required=False,
    help='Input files')

parser.add_argument(
    '-c',
    '--columns',
    metavar='COLUMNS',
    dest='columns',
    nargs = '+',
    default=None,
    required=True,
    help='Columns to plot. Either 2 or 2 * INPUT_FILES')

parser.add_argument(
    '--time',
    action='store_true',
    dest='x_is_time',
    help='Signal that the x column actually contains datetime data'
)

def parse_computation_string(fstr):
    """ Expect e.g. $1+2*$2+$3/$14 (for one file), or
                    0.5*$0:1+3.14/$1:11*$0:2
        which will result in the strings:
                    fcl[0][1]+2*fcl[0][2]+fcl[0][3]/fcl[0][14] and
                    0.5*fcl[0][1]+3.14/fcl[1][11]*fcl[0][2]
        respectively

        WARNING in the fcl 2D array, the file index is 0-offset but the column
                index is 1-offset, aka fcl[1][11] is the 10th column of the 1st
                file
    """
    print("Parsing {:}".format(fstr))
    if ":" not in fstr:
        return re.sub(r"\$([0-9]+)", r"fcl[0][\1]", fstr)
    for i in range(10):
        fstr = re.sub(r"\${:}:([0-9]+)".format(i), r"fcl[{:}][\1]".format(i), fstr)
    return fstr

def mapfile(fn, datetime_column=-1):
    with open(fn, 'r') as fin:
        if datetime_column<0:
            return [ line for line in fin.readlines() if not line.startswith('#') ]
        else:
            lines = [line for line in fin.readlines() if not line.startswith('#')]
            ## guess datetime format
            dtformats = ['%Y-%m-%d %H:%M:%S.%f', '%Y-%m-%d %H:%M:%S']
            testl = lines[0].split()
            limit_second_resoltion = 0
            for i in range(len(dtformats)):
                try:
                    t = datetime.datetime.strptime(' '.join([testl[datetime_column], testl[datetime_column+1]]), dtformats[i])
                    break
                except ValueError as v:
                    ## well, maybe second resolution is too big!
                    if i == 0 and len(v.args) > 0 and v.args[0].startswith('unconverted data remains: '):
                        nl = ' '.join([testl[datetime_column], testl[datetime_column+1]])[:-(len(v.args[0]) - 26)]
                        t = datetime.datetime.strptime(nl, dtformats[i])
                        limit_second_resoltion = (len(v.args[0]) - 26)
                        print("WARNING! Limmiting second resolution for file {:}".format(fn))
                        break
                except:
                    if i == len(dtformats-1):
                        errmsg = "ERROR. Failed to parse datetime string in file {:} from columns {:} to {:}".format(fn,datetime_column,datetime_column+1)
                        raise RuntimeError(errmsg)
            newlines = []
            for line in lines:
                l = line.split()
                t = datetime.datetime.strptime(
                    ' '.join([testl[datetime_column], testl[datetime_column+1]])[:-limit_second_resoltion], dtformats[i])
                newlines.append(
                    ' '.join(l[0:datetime_column] + [t] + l[datetime_column+1:]))
            return newlines

if __name__ == '__main__':

    args = parser.parse_args()

    #for fn in args.infiles:
    #    if not os.path.isdir(fn):
    #        print('ERROR. Failed to find file: {:}'.format(fn), file=sys.stderr)
    #        sys.exit(1)
    
    for c in args.columns:
        eval_str = parse_computation_string(c)
        print(eval_str)
    
    for fn in args.infiles:
        print(mapfile(fn,0))

