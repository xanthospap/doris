#!/usr/bin/python3
#-*- coding: utf-8 -*-

##
## 
##

import datetime
import argparse
import sys, os

def snx_stamp2dt(snx_str):
  ydoy, secday = snx_str[0:6], snx_str[7:]
  return datetime.datetime.strptime(ydoy, "%y:%j") + datetime.timedelta(days=0, seconds=int(secday))

def sn_solnestimate_block(fn):
  def goto_solnestimate_block(snx):
    line = snx.readline()
    while line and line.strip() != "+SOLUTION/ESTIMATE":
      line = snx.readline()
    if not line:
      print('Failed to find start of block \'+SOLUTION/ESTIMATE\' in SINEX file', file=sys.sdterr)
      raise RuntimeError()
    return 
  
  data = {}
  with open(fn, 'r') as snx:
    goto_solnestimate_block(snx)
    line = snx.readline()
    assert(line.strip() == "*INDEX TYPE__ CODE PT SOLN _REF_EPOCH__ UNIT S __ESTIMATED VALUE____ _STD_DEV___")
    while True:
      line = snx.readline()
      if line.strip() == "-SOLUTION/ESTIMATE": break
      l = line.split()
      index, type, code, pt, soln, ref_epoch, unit, s, estimated_val, std_dev = l
      ref_epoch = snx_stamp2dt(ref_epoch)
      estimated_val = float(estimated_val)
      std_dev = float(std_dev)
      if type[0:3] == "STA":
        assert(type[-1] in ["X","Y","Z"] and unit == "m")
      elif type[0:3] == "VEL":
        assert(type[-1] in ["X","Y","Z"] and unit == "m/y")
      else:
        print("Unknown solution type in SINEX: [{:}]".format(type), file=sys.stderr)
        raise RuntimeError()
      component = type[-1]
      d = {'index': index, 'type': type, 'pt': pt, 'soln': soln, 'ref_epoch': ref_epoch, 's': s, 'estimated_val': estimated_val, 'std_dev': std_dev}
      if code in data:
        data[code][type] = d
      else:
        data[code]={type: d}
    return data

class myFormatter(argparse.ArgumentDefaultsHelpFormatter,
                  argparse.RawTextHelpFormatter):
    pass

parser = argparse.ArgumentParser(
    formatter_class=myFormatter,
    description=
    'Read-in, parse and dump station information from a list of IDS log files',
    epilog=('''National Technical University of Athens,
    Dionysos Satellite Observatory\n
    Send bug reports to:
    Xanthos Papanikolaou, xanthos@mail.ntua.gr
    Aug, 2022'''))

parser.add_argument(
    '-o',
    '--output',
    metavar='OUTPUT',
    dest='save_as',
    default=None,
    required=False,
    help='Save the beacon information file using this file(name); can include path.')

parser.add_argument(
    '-s',
    '--snx',
    metavar='SINEX',
    dest='sinex',
    required=True,
    help='SINEX file to extract station coordinates from')

parser.add_argument(
    '-d',
    '--date',
    metavar='DATE',
    dest='date',
    required=True,
    type=lambda s: datetime.datetime.strptime(s, '%Y-%m-%d'),
    help='Extrapolate coordinate estimates to this date')

if __name__ == '__main__':

    args = parser.parse_args()

    snxdata = sn_solnestimate_block(args.sinex)

    ## for every site in the SINEX file
    for site, dct in snxdata.items():
      xyz = []
      for component in ['X', 'Y', 'Z']:
        dt = (args.date - dct['STA'+component]['ref_epoch']).days / 365.25e0
        xyz.append(dct['STA'+component]['estimated_val'] + dt * dct['VEL'+component]['estimated_val'])
        assert(dct['STA'+component]['ref_epoch'] == dct['VEL'+component]['ref_epoch'])
      print('{:24s} {:16.3f}{:16.3f}{:16.3f}'.format(site, xyz[0], xyz[1], xyz[2]))
