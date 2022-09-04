#!/usr/bin/python3

##
## fetch space weather data from CelesTrak database, in CSV format
##

import os
import sys
import datetime
import argparse
import requests
from requests.auth import HTTPBasicAuth
from urllib.parse import urlparse


class myFormatter(argparse.ArgumentDefaultsHelpFormatter,
                  argparse.RawTextHelpFormatter):
    pass


parser = argparse.ArgumentParser(
    formatter_class=myFormatter,
    description='Download Body and Solar Panel Qauternion files from CDDIS\nSee https://cddis.nasa.gov/Data_and_Derived_Products/DORIS/DORIS_quaternion_files.html',
    epilog=('''National Technical University of Athens,
    Dionysos Satellite Observatory\n
    Send bug reports to:
    Xanthos Papanikolaou, xanthos@mail.ntua.gr
    Jul, 2022'''))

parser.add_argument(
    '-s',
    '--satellite',
    metavar='SATELLITE',
    dest='satellite',
    default=None,
    required=True,
    choices=['ja1', 'ja2', 'ja3'],
    help='Target satellite (using a 3-character id')

parser.add_argument(
    '-d',
    '--date',
    metavar='DATE',
    dest='date',
    default=None,
    required=True,
    help='Target date in YYYY-MM-DD format')

parser.add_argument(
    '-O',
    '--output-dir',
    metavar='OUTPUT_DIR',
    dest='save_dir',
    required=False,
    default=os.getcwd(),
    help='Save the downloaded file under the given directory name.')

cddis_url = 'https://cddis.nasa.gov/archive/doris/ancillary/quaternions/ja3/2021/'

if __name__ == '__main__':

    args = parser.parse_args()

    if not os.path.isdir(args.save_dir):
        print('ERROR. Failed to locate directory: {:}'.format(
            args.save_dir), file=sys.stderr)
        sys.exit(1)
    
    try:
        t = datetime.datetime.strptime(args.date, '%Y-%m-%d')
    except:
        print('ERROR. Failed to parse date string from \'{:}\''.format(args.date), file=sys.stderr)
        sys.exit(1)

    cddis_url_append = '{:}/{:}'.format(args.satellite, t.year)
    remote_dir = os.path.join(cddis_url, cddis_url_append)

    ## get yearly directory listing
    md5sums = os.path.join(remote_dir, "MD5SUMS")
    print('target file: {:}'.format(md5sums))
    r = requests.get(md5sums, auth=HTTPBasicAuth('xanthos', 'Xanthos1984'))
    if not r.ok:
        print('ERROR Failed fetching remote file {:}. Error code {:}'.format(md5sums, r.status_code), file=sys.stderr)
        sys.exit(1)
    # Opens a local file of same name as remote file for writing to
    with open('foo', 'wb') as fd:
        for chunk in r.iter_content(chunk_size=1000):
            fd.write(chunk)

    # Closes local file
    fd.close()


