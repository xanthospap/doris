#!/usr/bin/python3
#-*- coding: utf-8 -*-

##
## Read-in, parse and dump station information from a list of IDS log files
##

import datetime
import argparse
import sys, os

def snx_stamp2dt(sxn_str):
  return datetime.datetime.strptime(snx_str, "%y:%j:%S")

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
    '-l',
    '--log-dir',
    metavar='LOG_DIR',
    dest='log_dir',
    required=False,
    default=os.getcwd(),
    help='Directory where the IDS log files are placed. Log files are expected to have the \'.LOG\' extension')

if __name__ == '__main__':

    args = parser.parse_args()

    if not os.path.isdir(args.log_dir):
        print('ERROR. Failed to locate directory: {:}'.format(args.log_dir), file=sys.stderr)
        sys.exit(1)
    
    if not args.save_as:
        fout = sys.stdout
    else:
        fout = open(args.save_as, "w")

    print("{:4s} {:19s} {:19s} {:9s}".format("Site", "From", "To", "Height"), file=fout)

    for fn in os.listdir(args.log_dir):
        if os.path.splitext(fn)[1] in [".log", ".LOG"]:
            idslog.dump_log_antenna_info(os.path.join(args.log_dir,fn),fout)
