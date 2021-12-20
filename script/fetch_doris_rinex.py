#! /usr/bin/python

## ------------------------------------------------------- ##
##  https://ids-doris.org/about-doris-rinex-format.html    ##
## ------------------------------------------------------- ##

import sys
import os
import datetime
import argparse
import shutil
import requests
from ftplib import FTP_TLS

## General Configuration - needs to be moved to a system-wide dir -
cddis_url = { 'domain': 'gdc.cddis.eosdis.nasa.gov', 'root_dir': '/pub/doris/data/' }

analysis_center = { 'grg': 'CNES/GRGS', 'gsc': 'NASA/GSFC', 'lca': 'LEGOS-CLS', 'ssa': 'SSALTO'}
all_acn = [ k for k in analysis_center ]
all_acn += [ analysis_center[k] for k in analysis_center ]

satellite_names = {
    'cs2' : 'Cryosat-2',
    'en1' : 'Envisat-1',
    'h2a' : 'HY-2A',
    'ja1' : 'Jason-1',
    'ja2' : 'Jason-2',
    'ja3' : 'Jason-3',
    's3a' : 'Sentinel-3A',
    's3b' : 'Sentinel-3B',
    's6a' : 'Sentinel-6A',
    'sp2' : 'SPOT-2',
    'sp3' : 'SPOT-3',
    'sp4' : 'SPOT-4',
    'sp5' : 'SPOT-5',
    'srl' : 'Saral',
    'top' : 'TOPEX/Poseidon'}
sats = [ k for k in satellite_names ]
sats += [ satellite_names[k] for k in satellite_names ]
##

def fetch_file(target_dir, target_fn, local_dir, local_fn=None):
    local_fn = target_fn if local_fn is None else local_fn
    local_fn = os.path.join(local_dir, local_fn)

    url = '{:}://{:}'.format('https', cddis_url['domain'])
    print('Downloading remote file {:} to {:} ...'.format(target_fn, local_fn), end='')

    error = 0
    try:
        ftps = FTP_TLS(host = cddis_url['domain'])
        ftps.login(user='anonymous', passwd='xanthos@mail.ntua.gr')
        ftps.prot_p()
        ftps.cwd(target_dir)
        ftps.retrbinary("RETR " + target_fn, open(local_fn, 'wb').write)
        print('done!')
    except:
        print('ERROR Failed to download file!')
        if os.path.isfile(target_fn): os.remove(target_fn)
        error = 1

    return error

def make_target_filename(sat_name, dt, version_nr=1):
    sss = sat_name
    rx  = 'rx'
    yy  = dt.strftime('%y')
    ddd = dt.strftime('%j')
    lll = '{:03d}'.format(version_nr)
    z   = 'Z'
    return sss + rx + yy + ddd + '.' + lll + '.' + z

def make_target_path(sat_name, dt):
    sss = sat_name
    yyyy = dt.strftime('%Y')
    return '{:}{:}/{:}/'.format(cddis_url['root_dir'], sss, yyyy)

class myFormatter(argparse.ArgumentDefaultsHelpFormatter,
                  argparse.RawTextHelpFormatter):
    pass


parser = argparse.ArgumentParser(
    formatter_class=myFormatter,
    description=
    'Download DORIS RINEX data file(s)',
    epilog=('''National Technical University of Athens,
    Dionysos Satellite Observatory\n
    Send bug reports to:
    Xanthos Papanikolaou, xanthos@mail.ntua.gr
    Dimitris Anastasiou,danast@mail.ntua.gr
    May, 2021'''))

parser.add_argument('-y',
                    '--year',
                    metavar='YEAR',
                    dest='year',
                    type=int,
                    required=False,
                    help='The year of date.')

parser.add_argument('-d',
                    '--doy',
                    metavar='DOY',
                    dest='doy',
                    type=int,
                    required=False,
                    help='The day-of-year (doy) of date.')

#parser.add_argument(
#    '-o',
#    '--output',
#    metavar='OUTPUT',
#    dest='save_as',
#    required=False,
#    help='Save the downloaded file using this file(name); can include path.')

parser.add_argument(
    '-O',
    '--output-dir',
    metavar='OUTPUT_DIR',
    dest='save_dir',
    required=False,
    default=os.getcwd(),
    help='Save the downloaded file under the given directory name.')

parser.add_argument(
    '-s',
    '--satellite',
    choices=sats,
    metavar='SATELLITE',
    dest='satellite',
    required=True,
    help=
    'Choose Satellite'
)

#parser.add_argument(
#    '-c',
#    '--data-center',
#    choices=['cddis', 'ign'],
#    metavar='DATA_CENTER',
#    dest='data_center',
#    default = 'cddis',
#    required=False,
#    help=
#    'Choose Data Center.'
#)

parser.add_argument('-l',
                    '--rinex-version-nr',
                    metavar='RINEX_VERSION_NUMBER',
                    dest='version_nr',
                    type=int,
                    required=False,
                    default=1,
                    help='Version number (starting with 001 for the initial version) when the file is replaced.')

if __name__ == '__main__':

    args = parser.parse_args()

    if (args.year is not None and args.doy is None) or (args.doy is not None and args.year is None):
        print('[ERROR] Need to specify both Year and DayOfYear', file=sys.stderr)
        sys.exit(1)

    if args.year is not None:
        t = datetime.datetime.strptime('{:4d}-{:03d}'.format(args.year, args.doy), '%Y-%j')
    else:
        t = datetime.datetime.now()

    if not os.path.isdir(args.save_dir):
        print('ERROR! Failed to find directory {:}'.format(args.save_dir), file=sys.stderr)
        sys.exit(1)
    
    if args.satellite not in sats:
        print('ERROR! Invalid satellite name; choose one from:', file=sys.stderr)
        for shortn, longn in satellite_names.items():
            print("\t{:10s} aka {:15s}".format(shortn, longn), file=sys.stderr)
        sys.exit(1)
    elif args.satellite in satellite_names:
        satellite = args.satellite
    else:
        for shortn, longn in satellite_names.items():
            if args.satellite == longn:
                satellite = shortn
                break

    target_rinex = make_target_filename(satellite, t, args.version_nr)
    target_path = make_target_path(satellite, t)
    sys.exit(fetch_file(target_path, target_rinex, args.save_dir))
