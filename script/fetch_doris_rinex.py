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

satellite_abbreviation_dct = { 
    'jason-2': 'ja2', 
    'cryosat-2': 'cs2', 
    'hy-2a': 'h2a', 
    'saral': 'srl', 
    'jason-3': 'ja3', 
    'sentinel-3a': 's3a', 
    'sentinel-3b': 's3b', 
    'hy-2c':'h2c',
    'sentinel-6a': 's6a'}
sats = [k for k in satellite_abbreviation_dct]
cddis_url = { 'domain': 'gdc.cddis.eosdis.nasa.gov', 'root_dir': '/pub/doris/data/' }
cddis_credentials = { 'username':'xanthos', 'password':'Xanthos1984'}

def fetch_file(target_dir, target_fn, local=None):
    local = target_fn if local is None else local

    url = '{:}://{:}'.format('https', cddis_url['domain'])

    print('Downloading remote file {:} ...'.format(url+target_dir+target_fn))

    error = 0
    try:
        ftps = FTP_TLS(host = cddis_url['domain'])
        ftps.login(user='anonymous', passwd='xanthos@mail.ntua.gr')
        ftps.prot_p()
        ftps.cwd(target_dir)
        ftps.retrbinary("RETR " + target_fn, open(target_fn, 'wb').write)
    except:
        print('ERROR Failed to download file!')
        if os.path.isfile(target_fn): os.remove(target_fn)
        error = 1

    return error

def make_target_filename(sat_name, dt, version_nr=1):
    sss = satellite_abbreviation_dct[sat_name.lower()] if len(sat_name)>3 else sat_name
    rx  = 'rx'
    yy  = dt.strftime('%y')
    ddd = dt.strftime('%j')
    lll = '{:03d}'.format(version_nr)
    z   = 'Z'
    return sss + rx + yy + ddd + '.' + lll + '.' + z

def make_target_path(sat_name, dt):
    sss = satellite_abbreviation_dct[sat_name.lower()] if len(sat_name)>3 else sat_name
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

parser.add_argument(
    '-o',
    '--output',
    metavar='OUTPUT',
    dest='save_as',
    required=False,
    help='Save the downloaded file using this file(name); can include path.')

parser.add_argument(
    '-O',
    '--output-dir',
    metavar='OUTPUT_DIR',
    dest='save_dir',
    required=False,
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
        print('[ERROR] Need to specify both Year and DayOfYear')
        sys.exit(1)

    if args.year is not None:
        t = datetime.datetime.strptime('{:4d}-{:03d}'.format(args.year, args.doy), '%Y-%j')
    else:
        t = datetime.datetime.now()

    input_dct = {}
    if args.save_as:
        input_dct['save_as'] = args.save_as
    if args.save_dir:
        input_dct['save_dir'] = args.save_dir

    target_rinex = make_target_filename(args.satellite, t, args.version_nr)
    target_path = make_target_path(args.satellite, t)
    sys.exit(fetch_file(target_path, target_rinex))
