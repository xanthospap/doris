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

def lines2column_dictionary(lines):
    """ Give a list of strings (lines), parse them to a dictionary of columns.
        The resulting dictionary will have as keys the column number/indexes
        and values lists with the columns
        E.g.
        lines = ['abc def foo bar', 'koko baz ttt aaa'] will return
        {0:[abc,koko], 1:[def,baz], 2:[foo,ttt], 3:[bar,aaa]}
    """ 
    ## number of columns in each non-comment line
    numcolumns = [len(l.split()) for l in lines if not l.startswith('#')]
    ## all lines should have the same number of columns
    assert(all([nc==numcolumns[0] for nc in numcolumns]))
    ## a dictionary of n lists, where n is the number of columns. keys are 
    # column indexes
    numcolumns = numcolumns[0]
    colsarray = {}
    for i in range(numcolumns): colsarray[i] = []
    ## append column arrays, one at a time
    dummy_it = 0
    for line in lines:
        if not line.startswith('#'):
            l = line.split()
            for i in range(numcolumns):
                if i==0 and dummy_it < 3:
                    dummy_it += 1
                colsarray[i].append(l[i])
    return colsarray

def column_dictionary2time(coldict,datetime_column):
    """ Given a dictionay (coldict) which is the results of a call to 
        lines2column_dictionary, this function will try to concatenate columns 
        datetime_column and datetime_column+1 and parse them as a single 
        datetime instance.
        The function will return coldict where column with key=datetime_column
        will be an array of datetime.datetime instances. All other arrays of
        the dictionary will be kept as are at input.
    """
    dtformats = ['%Y-%m-%d %H:%M:%S.%f', '%Y-%m-%d %H:%M:%S']
    test_column = ' '.join([coldict[datetime_column][0],coldict[datetime_column+1][0]])
    clip_at = 0
    for format_idx in range(len(dtformats)):
        try:
            t = datetime.datetime.strptime(test_column, dtformats[format_idx])
            break
        except ValueError as v:
            if format_idx == 0 and len(v.args) > 0 and v.args[0].startswith('unconverted data remains: '):
                test_column = ' '.join([coldict[datetime_column][0],coldict[datetime_column+1][0][:-(len(v.args[0]) - 26)]])
                t = datetime.datetime.strptime(test_column, dtformats[format_idx])
                clip_at = len(v.args[0]) - 26
                print("WARNING! Truncating fractional seconds for Python to parse!", file=sys.stderr)
                break
        except:
            if format_idx == len(dtformats-1):
                raise RuntimeError()
    datetime_array = [datetime.datetime.strptime(' '.join([dc[0],dc[1][:-clip_at]]), dtformats[format_idx]) for dc in zip(coldict[datetime_column], coldict[datetime_column+1])]
    coldict[datetime_column] = datetime_array
    return coldict

def drop_cols(coldict,cols_to_keep):
    for key in coldict:
        if key not in cols_to_keep:
            coldict[key]=[]
    return coldict

def mapfile(fn, datetime_column=-1):
    with open(fn, 'r') as fin:
        d = lines2column_dictionary(fin.readlines())
        if datetime_column >= 0:
            d = column_dictionary2time(d,datetime_column)
    return d

if __name__ == '__main__':

    args = parser.parse_args()

    eval_strings = []
    for c in args.columns:
        eval_strings.append(parse_computation_string(c))
    assert(not len(eval_strings)%2)
    
    for fn in args.infiles:
        if not os.path.isdir(fn):
            print('ERROR. Failed to find file: {:}'.format(fn), file=sys.stderr)
            sys.exit(1)
        d = mapfile(fn, 