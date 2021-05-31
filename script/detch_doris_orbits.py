#! /usr/bin/python

##
##  https://cddis.nasa.gov/Data_and_Derived_Products/DORIS/doris_idsorbit.html
##

import sys
import os
import datetime
import argparse
import urllib.request
from pybern.products.downloaders.retrieve import ftp_retrieve
from pybern.products.downloaders.ftplist import ftp_dirlist

URL = 'https://cddis.nasa.gov/archive/doris/products/orbits/ssa/'

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
    'sp2' : 'SPOT-2',
    'sp3' : 'SPOT-3',
    'sp4' : 'SPOT-4',
    'sp5' : 'SPOT-5',
    'srl' : 'Saral',
    'top' : 'TOPEX/Poseidon'}
all_sat_names = [ k for k in satellite_names ]
all_sat_names += [ satellite_names[k] for k in satellite_names ]

def make_target_filename(acn='ssa', sss='cs2', dt=datetime.datetime.now(), vv=1, has_doris=True, has_gps=False, has_slr=False, lll=1):
    """ cccsssVV.bXXDDD.eYYEEE.dgs.sp3.LLL.Z
        ccc : analysis_center
        sss : satellite name
        VV  : version number
        XX  : last two digits of year of first position
        DDD : three-digit day of year of first position
        YY  : last two digits of year of last position
        EEE : three-digit day of year of last position
        dgs : d = "D" or "_" if DORIS data are used or not
            : g = "G" or "_" if GPS data are used or not
            : s = "S" or "_" if SLR data are used or not
        LLL : version number (starting with 001) when the file is replaced
    """
    if acn not in analysis_center:
        for k,v in analysis_center.items():
            if acn == v:
                acn = k
                break
    if acn not in analysis_center:
        print('[ERROR] Invalid analysis center \'{:}\''.format(acn))
        raise RuntimeError('Invalid analysis center')

    if sss not in satellite_names:
        for k,v in satellite_names.items():
            if sss == v:
                sss=k
                break
    if sss not in satellite_names:
        print('[ERROR] Invalid satellite id')
        raise RuntimeError('Invalid satellite id')

    vv  = '{:02d}'.format(vv)
    lll = '{:03d}'.format(lll)

    dgs = ['___']
    if has_doris:
        dgs = ['D__']
    else:
        dgs.append('D__')
    
    if has_gps:
        for d in dgs:
            d[1] = 'G'
    else:
        tt = []
        for d in dgs:
            new_d = d[0] + 'G' + d[2]
            tt.append(new_d)
        dgs += [ d for d in tt if d not in dgs ]

    if has_slr:
        for d in dgs:
            d[2] = 'S'
    else:
        tt = []
        for d in dgs:
            new_d = d[0] + d[1] + 'S'
            tt.append(new_d)
        dgs += [ d for d in tt if d not in dgs ]

    xx = dt.strftime('%y')
    ddd = dt.strftime('%j')

    return [ '{:}{:}{:}{:}{:}YYEEE{:}{:}.Z'.format(acn, sss, vv, xx, ddd, dgs_, lll) for dgs_ in dgs ]

def make_target_path(acn='ssa', sss='cs2'):
    if acn not in analysis_center:
        for k,v in analysis_center:
            if acn == v:
                acn = k
                break
    if acn not in analysis_center:
        print('[ERROR] Invalid analysis center \'{:}\''.format(acn))
        raise RuntimeError('Invalid analysis center')

    if sss not in satellite_names:
        for k,v in satellite_names:
            if sss == v:
                sss=k
                break
    if sss not in satellite_names:
        print('[ERROR] Invalid satellite id')
        raise RuntimeError('Invalid satellite id')

    return '{:}{:}/{:}/'.format(URL, acn, sss)

class myFormatter(argparse.ArgumentDefaultsHelpFormatter,
                  argparse.RawTextHelpFormatter):
    pass


parser = argparse.ArgumentParser(
    formatter_class=myFormatter,
    description=
    'Download DORIS Sp3c/d data file(s)',
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
                    required=True,
                    help='The year of date.')

parser.add_argument('-d',
                    '--doy',
                    metavar='DOY',
                    dest='doy',
                    type=int,
                    required=True,
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
    choices=all_sat_names,
    metavar='SATELLITE',
    dest='satellite',
    required=True,
    help=
    'Choose Satellite'
)

parser.add_argument(
    '-c',
    '--data-center',
    choices=all_acn,
    metavar='DATA_CENTER',
    dest='data_center',
    default = 'ssa',
    required=False,
    help=
    'Choose Data Center.'
)

parser.add_argument('-l',
                    '--sp3-version-nr',
                    metavar='SP3_VERSION_NUMBER',
                    dest='version_nr',
                    type=int,
                    required=False,
                    default=1,
                    help='Version number (starting with 001 for the initial version) when the file is replaced.')

if __name__ == '__main__':

    args = parser.parse_args()

    t = datetime.datetime.strptime('{:}-{:03d}'.format(args.year, args.doy), '%Y-%j')

    input_dct = {}
    if args.save_as:
        input_dct['save_as'] = args.save_as
    if args.save_dir:
        input_dct['save_dir'] = args.save_dir

    vv = 1

    target_sp3 = make_target_filename(args.data_center, args.satellite, t, 1, True, False, False, args.version_nr)
    target_path = make_target_path(args.data_center, args.satellite)
    print('target_path \'{:}\''.format(target_path))
    for i in target_sp3: print('\ttarget: \'{:}\''.format(i))
    remote_files = ftp_dirlist(
    """
    status = 200
    for sp3 in target_sp3:
        url = '{:}{:}'.format(target_path, sp3)
        try :
            status, target, saveas = ftp_retrieve(url, **input_dct)
            break
        except RuntimeError:
            print('Failed to download target file \'{:}\''.format(sp3))
    sys.exit(status)
    """
