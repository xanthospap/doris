#!/usr/bin/python3

##
## Fetch IERS Bulletin B/C04 file(s) for a given date. 
## Note that an IERS Bulletin B/C04 file cantains final values for up to one
## month (prior to publication) and preliminery values for up to one month 
## after this. This script will download the respective file for the given 
## date. If the date is within N days of the start of preliminery values, then
## it will also download the next (following) file concatenate the two to a
## new, merged Bulletin B/C04 file, thus containing final values for more than 
## one month!
##

import os, sys, datetime
import argparse
import requests
from urllib.parse import urlparse

class myFormatter(argparse.ArgumentDefaultsHelpFormatter,
                  argparse.RawTextHelpFormatter):
    pass

parser = argparse.ArgumentParser(
    formatter_class=myFormatter,
    description=
    'Download IERS Bulletin B/C04 file(s) for a given date and if needed create a new, merged one.',
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
    '-N',
    '--day-limits',
    metavar='DAY_LIMIT',
    dest='day_limit',
    required=False,
    type=int,
    default=3,
    help='The target date and -N days and +N days to target date must be included in the final file')

parser.add_argument(
    '-d',
    '--date',
    metavar='DATE',
    dest='date',
    required=True,
    help='The target date')


bc04_url = "https://hpiers.obspm.fr/iers/bul/bulb_new/bulletinb."
bc04_nr0 = 252
bc04_first_date = datetime.datetime.strptime("20081212", "%Y%m%d")

## Give a Bulletin B/C04 file, open it and extract the time interval 
## (start/stop datetimes) of the final data block.
def c04_get_final_data(c04fn):
    final_data = []
    error = 0
    
    with open(c04fn, 'r') as fin:
        line = fin.readline()
        
        lines_matched = 0
        while line.strip() != 'Mean formal error':
            line = fin.readline()
            if not line:
                error = 1
                break
            if line.strip() == "1 - DAILY FINAL VALUES OF x, y, UT1-UTC, dX, dY":
              lines_matched += 1
            elif ' '.join(line.strip().split()) == "DATE MJD x y UT1-UTC dX dY x err y err UT1 err dX err dY err":
              lines_matched += 1
            elif ' '.join(line.strip().split()) == "(0 h UTC) mas mas ms mas mas mas mas ms mas mas":
              lines_matched += 1
        
        if lines_matched != 3:
          print('ERROR Failed to validate lines in Bulletin B file {:}'.format(c04fn), file=sys.stderr)
        
        if not error:
            line = fin.readline() # empty ...
            line = fin.readline() # Mean formal error
            if not line.lstrip().startswith('Mean formal error') :
                print('ERROR. Expected \'Mean formal error\' not found in file {:}'.format(c04fn), file=sys.stderr)
                error = 1
            line = fin.readline() # empty ...

        if not error:
            ## read actual final data
            line = fin.readline()
            while line.strip() != '':
                ## int values
                int_list = [int(d) for d in line.split()[0:4]]
                ## floating values, should be 10
                float_list = [ float(f) for f in line.split()[4:]
                assert(len(float_list) == 10)
                ## append collected data
                final_data.append(int_list + float_list)
                ## next line
                line = fin.readline()
    
    if not error:
        return dates_in_file[0], dates_in_file[-1]
    else:
        raise RuntimeError('Failed parsing Bulletin B/C04 file')

## Give a Bulletin B/C04 file, open it and extract the EOP data for the 
## (final data block.
def c04_date_span(c04fn):
    error = 0
    dates_in_file = []
    
    with open(c04fn, 'r') as fin:
        line = fin.readline()
        
        while line.strip() != 'Final values':
            line = fin.readline()
            if not line:
                error = 1
                print('ERROR. Failed finding \'Final values\' in file {:}'.format(c04fn), file=sys.stderr)
                break
        
        if not error:
            line = fin.readline() # empty ...
            line = fin.readline() # Mean formal error
            if not line.lstrip().startswith('Mean formal error') :
                print('ERROR. Expected \'Mean formal error\' not found in file {:}'.format(c04fn), file=sys.stderr)
                error = 1
            line = fin.readline() # empty ...

        if not error:
            ## read actual final data
            line = fin.readline()
            while line.strip() != '':
                y,m,d = [int(d) for d in line.split()[0:3]]
                d = '{:}-{:02d}-{:02d}'.format(y,m,d)
                dates_in_file.append(datetime.datetime.strptime(d,'%Y-%m-%d'))
                line = fin.readline()
    
    if not error:
        return dates_in_file[0], dates_in_file[-1]
    else:
        raise RuntimeError('Failed parsing Bulletin B/C04 file')

## cast a dictionary of type {filename1:
##                           {'t0':datetime_start, 't':datetime_stop},
##                           {filename2:
##                           {'t0':datetime_start, 't':datetime_stop}, ...
## to  a sorted list of filenames, base on the value 't0'
def c04_dict_to_sorted_list(c04_dict):
    lst = []
    ct = datetime.datetime.max()
    for fn,vals in c04_dict.items():
        if vals[0] < ct:
            lst = [fn] + lst
            ct = vals[0]
        else:
            lst = lst.append(fn)
            if len(lst) == 1: ct = vals[0]
    return lst

def merge_c04_files(c04_dict):
    if len(c04_dict) == 0: return None
    elif len(c04_dict) == 1: return list(c04_dict.keys())[0]
    else:
      ## filenames in chronological order
      c04_list = c04_dict_to_sorted_list(c04_dict)
      final_data_extra = [ c04_get_final_data(fn) for f in c04_list[1:]]

if __name__ == '__main__':

    args = parser.parse_args()

    t = datetime.datetime.strptime(args.date, '%Y-%m-%d')
    t_start = t - datetime.timedelta(days=args.day_limit)
    t_end   = t + datetime.timedelta(days=args.day_limit)

    # covered span(s)
    prior_covered = False
    post_covered = False
    
    ## files to merge
    files_to_merge = {}

    # guess the number of the file
    dummy_it = 0
    guess = bc04_nr0 + (t_start - bc04_first_date).days // 31
    while not prior_covered or not post_covered and dummy_it < 1000:

        ## download remote file
        remote = bc04_url + '{:}'.format(guess)
        local  = os.path.basename(remote)
        data = requests.get(remote, allow_redirects=True)
        open(local, 'wb').write(data.content)

        ## get its time span
        d0,dn = c04_date_span(local)

        ## will we be needing it?
        file_is_of_use = False

        ## check if it overlaps the requested interval
        if not prior_covered and d0 < t_start and t_start < dn:
            prior_covered = True
            file_is_of_use = True

        if not post_covered and d0 < t_end and t_end < dn:
            post_covered = True
            file_is_of_use = True

        ## if not to be used, remove local
        if file_is_of_use:
            files_to_merge[local] = {'t0': d0, 'tn': dn}
        else:
            os.remove(local)

        ## next guess
        if d0 > t_start: guess -= 1
        elif dn < t_end: guess += 1
        else:
           raise RuntimeError('WTF??')

        dummy_it += 1
