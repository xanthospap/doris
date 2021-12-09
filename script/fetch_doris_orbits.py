#! /usr/bin/python

##
##  https://cddis.nasa.gov/Data_and_Derived_Products/DORIS/doris_idsorbit.html
##

import sys
import os
import re
import datetime
import argparse
from ftplib import FTP_TLS

analysis_center = { 'grg': 'CNES/GRGS', 'gsc': 'NASA/GSFC', 'lca': 'LEGOS-CLS', 'ssa': 'SSALTO'}
all_acn = [ k for k in analysis_center ]
all_acn += [ analysis_center[k] for k in analysis_center ]
cddis_url = { 'domain': 'gdc.cddis.eosdis.nasa.gov', 'root_dir': '/pub/doris/data/' }
cddis_credentials = { 'username':'xanthos', 'password':'Xanthos1984'}

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
all_sat_names = [ k for k in satellite_names ]
all_sat_names += [ satellite_names[k] for k in satellite_names ]

def fetch_file(target_dir, target_fn, local_dir, local_fn=None):
    local_fn = target_fn if local_fn is None else local_fn
    local_fn = os.path.join(local_dir, local_fn)

    url = '{:}://{:}'.format('https', cddis_url['domain'])
    # print('Downloading remote file {:} ...'.format(url+target_dir+target_fn))

    error = 0
    try:
        ftps = FTP_TLS(host = cddis_url['domain'])
        ftps.login(user='anonymous', passwd='xanthos@mail.ntua.gr')
        ftps.prot_p()
        ftps.cwd(target_dir)
        ftps.retrbinary("RETR " + target_fn, open(local_fn, 'wb').write)
    except:
        print('ERROR Failed to download file!')
        if os.path.isfile(target_fn): os.remove(target_fn)
        error = 1

    return error

def listFD(target_dir):
    ftps = FTP_TLS(host = 'gdc.cddis.eosdis.nasa.gov')
    ftps.login(user='anonymous', passwd="xanthos@mail.ntua.gr")
    ftps.prot_p()
    ftps.cwd(target_dir)
    listing = []
    ftps.dir(listing.append)
    """ Answer is a list of strings like:
    -rw-rw-r--    1 ftp      ftp          7304 Dec 04 14:20 SHA512SUMS
    -rw-rw-rw-    1 ftp      ftp        465685 Apr 12  2019 ssasrl20.b13073.e13080.D__.sp3.001.Z
    -rw-rw-rw-    1 ftp      ftp        462369 Apr 12  2019 ssasrl20.b13080.e13087.D__.sp3.001.Z
    """
    return [ row.split()[8] for row in listing ]

def match_orbit_file(ftp_file_list, pattern_list, dt):
    """ e.g pattern[0] = 'ssasrl01.bXXDDD.eYYEEE.D__.sp3.001.Z'
    """
    matched_files = []
    for pattern in pattern_list:
        ## compile regex to match against
        pattern = re.sub(r'\.bXXDDD\.', '.b([0-9]{5}).', pattern)
        pattern = re.sub(r'\.eYYEEE\.', '.e([0-9]{5}).', pattern)
        pattern = re.sub(r'\.', '\.', pattern)
        rgx = re.compile(pattern)
        # print('>> regex to match against \'{:}\''.format(pattern))
     
        for remote_fn in ftp_file_list:
            m = re.match(rgx, remote_fn)
            if m:
                # print('>> matched remote file {:}'.format(remote_fn))
                start = datetime.datetime.strptime(m.group(1), '%y%j')
                stop  = datetime.datetime.strptime(m.group(2), '%y%j')
                if dt >= start and dt < stop:
                    matched_files.append(remote_fn)
        
        if matched_files != []: return matched_files
    return matched_files


def make_target_filename(acn='ssa', sss='cs2', vv='[012]{2}', has_doris=True, has_gps=False, has_slr=False, lll=-1):
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

    if lll < 0:
        lll = '[0-9]{3}'
    else:
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

    return [ '{:}{:}{:}.bXXDDD.eYYEEE.{:}.sp3.{:}.Z'.format(acn, sss, vv, dgs_, lll) for dgs_ in dgs ]

def make_target_dir(acn='ssa', sss='cs2'):
    if acn not in analysis_center:
        for k,v in analysis_center:
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

    return '{:}/{:}/{:}/'.format('/doris/products/orbits', acn, sss)

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
    default=os.getcwd(),
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
                    default=-1,
                    help='Version number (starting with 001 for the initial version) when the file is replaced. -1 Means any version.')

if __name__ == '__main__':

    args = parser.parse_args()

    t = datetime.datetime.strptime('{:}-{:03d}'.format(args.year, args.doy), '%Y-%j')

    ## handle version
    if args.version_nr < 0:
        version = '[0-9]{2}'
    else:
        version = '{:02d}'.format(args.version_nr)
    
    if not os.path.isdir(args.save_dir):
        print('ERROR! Failed to find directory {:}'.format(args.save_dir))
        sys.exit(1)

    ## make a list of target sp3 filenames
    target_sp3 = make_target_filename(args.data_center, args.satellite, version, True, False, False, 1)

    ## make the directory path
    target_dir = make_target_dir(args.data_center, args.satellite)

    ## get the remote directory listing
    ftp_file_list = listFD(target_dir)

    ## match the target sp3
    matched_remote_fns = match_orbit_file(ftp_file_list, target_sp3, t)
    #for remote_fn in matched_remote_fns:
    #    print('>> File to download: {:}'.format(remote_fn))a

    if len(matched_remote_fns) > 1:
        print('WARNING: More than one orbit files match the description:')
        for rfn in matched_remote_fns: print('\t Matched file: {:}'.format(rfn))

    sys.exit(fetch_file(target_dir, matched_remote_fns[0], args.save_dir))
