#!/usr/bin/python3

##
## fetch space weather data from CelesTrak database, in CSV format
##

import os, sys
import argparse
import requests
from urllib.parse import urlparse

class myFormatter(argparse.ArgumentDefaultsHelpFormatter,
                  argparse.RawTextHelpFormatter):
    pass

parser = argparse.ArgumentParser(
    formatter_class=myFormatter,
    description=
    'Download CelesTrak Space Weather Data in CSV format.\nSee https://celestrak.org/SpaceData/',
    epilog=('''National Technical University of Athens,
    Dionysos Satellite Observatory\n
    Send bug reports to:
    Xanthos Papanikolaou, xanthos@mail.ntua.gr
    Jul, 2022'''))

parser.add_argument(
    '-o',
    '--output',
    metavar='OUTPUT',
    dest='save_as',
    default=None,
    required=False,
    help='Save the downloaded file using this file(name); can include path.')

parser.add_argument(
    '-O',
    '--output-dir',
    metavar='OUTPUT_DIR',
    dest='save_dir',
    required=False,
    default=os.getcwd(),
    help='Save the downloaded file under the given directory name.')

parser.add_argument(
    '--csv-all',
    action='store_true',
    dest='csv_all',
    help='Download the \'historic\' CSV file, containing data from 1957 to now. By default, the script will download the latest CSV file, containing data for the last 5 years.'
)

remote_csv_all = "https://celestrak.org/SpaceData/SW-All.csv"
remote_csv_5yr = "https://celestrak.org/SpaceData/SW-Last5Years.csv"

if __name__ == '__main__':

    args = parser.parse_args()

    if not os.path.isdir(args.save_dir):
        print('ERROR. Failed to locate directory: {:}'.format(args.save_dir), file=sys.stderr)
        sys.exit(1)

    remote = remote_csv_all if args.csv_all else remote_csv_5yr

    local_fn = os.path.basename(remote) if not args.save_as else args.save_as
    local = os.path.join(args.save_dir, local_fn)

    data = requests.get(remote, allow_redirects=True)
    open(local, 'wb').write(data.content)

    sys.exit(0)
