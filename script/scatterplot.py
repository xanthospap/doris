#! /usr/bin/python3

import datetime
import matplotlib.pyplot as plt
import numpy as np
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
    '-w',
    '--where',
    metavar='WHERE',
    dest='where',
    nargs = '+',
    default=None,
    required=False,
    help='where clause, e.g. \"-w \'$0:4\'=DION')

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
    """ Given a list of strings (lines), parse them to a dictionary of columns.
        The resulting dictionary will have as keys the column number/indexes
        and values lists with the columns
        E.g.
        lines = ['abc def foo bar', 'koko baz ttt aaa'] will return
        {1:[abc,koko], 2:[def,baz], 3:[foo,ttt], 4:[bar,aaa]}

        !! Warning !! Column indexes (keys in the returned dictionary) are 
                      1-oofset (not 0-offset)
    """ 
    ## number of columns in each non-comment line
    numcolumns = [len(l.split()) for l in lines if not l.startswith('#')]
    ## all lines should have the same number of columns
    assert(all([nc==numcolumns[0] for nc in numcolumns]))
    ## a dictionary of n lists, where n is the number of columns. keys are 
    # column indexes
    numcolumns = numcolumns[0]
    colsarray = {}
    for i in range(numcolumns): colsarray[i+1] = []
    ## append column arrays, one at a time
    dummy_it = 0
    for line in lines:
        if not line.startswith('#'):
            l = line.split()
            for i in range(numcolumns):
                if i==0 and dummy_it < 3:
                    dummy_it += 1
                colsarray[i+1].append(l[i])
    ## lists to numpy arrays if possible
    for k in colsarray:
      try:
        float(colsarray[k][0])
        colsarray[k] = np.asarray(colsarray[k],dtype=np.float64)
      except:
        pass
    return colsarray

def column_dictionary2time(coldict,datetime_column):
    """ Given a dictionay (coldict) which is the results of a call to 
        lines2column_dictionary, this function will try to concatenate columns 
        datetime_column and datetime_column+1 and parse them as a single 
        datetime instance.
        The function will return coldict where column with key=datetime_column
        will be an array of datetime.datetime instances. All other arrays of
        the dictionary will be kept as are at input.

        !! datetime_column should be 1-offset (not 0-offset)!!
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

#def apply_where_clause(coldict, clause_dict):
#  """ clause_dict = {0:'DION', 4:'repro22'}
#  """
#  for key,val in clause_dict:
#      if key in coldict:


def drop_cols(fncount, cols_dict, eval_str_list):
  colset = set()
  ## check the elements of the eval string, and make a list of the columns
  ## we are actually going to need
  for ele in eval_str_list:
    g = re.findall(r"(?=fcl\[{:}\]\[([0-9]+)\])".format(fncount), ele)
    if g and g != []:
      for c in g: colset.add(int(c))
  ## new dictionay with columns to keep
  newd = {}
  for c in colset:
    newd[c] = cols_dict[c]
  return newd

def mapfile(fn, datetime_column=-1):
    with open(fn, 'r') as fin:
        d = lines2column_dictionary(fin.readlines())
        if datetime_column >= 0:
            d = column_dictionary2time(d,datetime_column)
    return d

def datetime_start_column(estr):
  """ Assuming that the 'estr' string is of type:
    fcl[int1][int2] where int1 is the filename counter and int2 is the
    1-offste column index (i.e. a string parsed from 'parse_computation_string',
    return int1 and int2
  """
  g = re.match(r"fcl\[([0-9]+)\]\[([0-9]+)\]", estr)
  if not len(g.groups()) == 2:
    raise RuntimeError('ERROR Failed to parse x-time column')
  else:
    return int(g.group(1)), int(g.group(2))

if __name__ == '__main__':

    args = parser.parse_args()

    ## parse eval strings/columns
    eval_strings = []
    for c in args.columns:
        eval_strings.append(parse_computation_string(c))
    #assert(not len(eval_strings)%2)
    #print(eval_strings)

    # one list entry per file
    fcl = []
  
    ## for each input file, parse columns and keep the ones we need  
    for count,fn in enumerate(args.infiles):
        if not os.path.isfile(fn):
            print('ERROR. Failed to find file: {:}'.format(fn), file=sys.stderr)
            sys.exit(1)
        
        ## datetime column(s) ?
        if not args.x_is_time:
          datetime_column = -1
        else:
          if len(eval_strings) == 2:
            fnc, datetime_column = datetime_start_column(eval_strings[0])
          else:
            fnc, datetime_column = datetime_start_column(eval_strings[count*2])
          assert(fnc == count)

        ## parse the file (time as datetime)
        d = mapfile(fn, datetime_column)
        d = drop_cols(count, d, eval_strings)
        fcl.append(d)
        # print(d)

    ## make the scatterplot
    nsubplots = len(eval_strings) // 2
    fig, axs = plt.subplots(nsubplots, 1, figsize=(10, 6), constrained_layout=True)
    for nplt in range(nsubplots):
      xdata = eval(eval_strings[nplt*2])
      #print(xdata)
      # ydata_str = re.sub(r"(fcl\[[0-9]+\]\[[0-9]+\])", r"np.asarray(\1)", eval_strings[nplt*2+1])
      ydata_str = eval_strings[nplt*2+1]
      # print(ydata_str)
      ydata = eval(ydata_str)
      #print(ydata)
      #ydata_str = eval_strings[nsubplots*2+1]
      if nsubplots == 1:
        axs.scatter(xdata,ydata,color='black')
      else:
        axs[nplt].scatter(xdata,ydata,color='black')
      # axs[nplt].set_title('{:}'.format())
    plt.show()
  
