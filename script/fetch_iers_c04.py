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
## Update 08/10/2022 - Add values for LOD
##

import os, sys, datetime
import argparse
import requests
import copy
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
    '-n',
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

parser.add_argument(
    '--verbose',
    action='store_true',
    dest='verbose',
    help='Verbose mode on')

bc04_url = "https://hpiers.obspm.fr/iers/bul/bulb_new/bulletinb."
bc04_nr0 = 252
bc04_first_date = datetime.datetime.strptime("20081212", "%Y%m%d")
header = """#1 - DAILY FINAL VALUES OF x, y, UT1-UTC, dX, dY
# Angular unit is milliarcsecond (mas), time unit is millisecond (ms). 
# Upgraded solution from March 1 2017 - consistent with ITRF 2014.
# 
#      DATE     MJD       x        y       UT1-UTC     dX        dY       x err     y err    UT1 err   dX err    dY err   LOD             sigma           OMEGA           sigma
#   (0 h UTC)            mas      mas        ms        mas       mas       mas       mas        ms       mas       mas      ms            ms              mas/ms          mas/ms
"""

## Given a Bulletin B/C04 file, open it and extract the EOP data for the 
## final data block, from subsection 3: 
## 3 - EARTH ANGULAR VELOCITY : DAILY FINAL VALUES OF LOD, OMEGA AT 0hUTC
def c04_get_final_data_03(c04fn):
    final_data = []
    error = 0
    
    with open(c04fn, 'r') as fin:
        line = fin.readline()

        lines_to_match = [
          '3 - EARTH ANGULAR VELOCITY : DAILY FINAL VALUES OF LOD, OMEGA AT 0hUTC',
          'LOD : Excess of the Length of day - 86400 s TAI',
          'OMEGA : Earth angular velocity',
          'DATE MJD LOD sigma OMEGA sigma'
          # '(0 h UTC) ms ms mas/ms mas/ms'
        ]
        
        lines_matched = 0;
        while not line.strip().lower().startswith('(0 h utc)              ms       ms        mas/ms          mas/ms'):
            if not line:
                error = 1
                break
            if ' '.join(line.strip().split()) in lines_to_match:
              lines_matched += 1
            line = fin.readline()
        
        if lines_matched != len(lines_to_match) or error:
          errmsg = 'Failed to extract data block 3 from Bulletin B/C04 file {:} (error #{:})'.format(c04fn, error)
          raise RuntimeError(errmsg)
        
        line = fin.readline() # empty ...
        assert(line.strip() == '')

        keys = ["year", "month", "day", "mjd", "lod", "lod_sigma", "omega", "omega_sigma"]
        if not error:
            ## read actual final data
            line = fin.readline()
            while line.strip() != '':
                d = {}
                ## int values (year, month, day_of_month, mjd)
                int_list = [int(d) for d in line.split()[0:4]]
                ## floating values, should be 4
                float_list = [ float(f) for f in line.split()[4:]]
                assert(len(float_list) == 4)
                ## make dictionary
                for k,v in zip(keys, int_list+float_list): d[k]=v
                ## append collected data
                final_data.append(d)
                ## next line
                line = fin.readline()
 
    if not error:
        return final_data
    else:
        errmsg = 'Failed to extract data block 3 from Bulletin B/C04 file {:} (error #{:})'.format(c04fn, error)
        raise RuntimeError(errmsg)


## Given a Bulletin B/C04 file, open it and extract the EOP data for the 
## final data block, from subsection 1: 
## "1 - DAILY FINAL VALUES OF x, y, UT1-UTC, dX, dY"
def c04_get_final_data_01(c04fn):
    final_data = []
    error = 0
    
    with open(c04fn, 'r') as fin:
        line = fin.readline()

        lines_to_match = [
          '1 - DAILY FINAL VALUES OF x, y, UT1-UTC, dX, dY',
          'DATE MJD x y UT1-UTC dX dY x err y err UT1 err dX err dY err',
          '(0 h UTC) mas mas ms mas mas mas mas ms mas mas'
        ]
        
        lines_matched = 0;
        while not line.strip().lower().startswith('mean formal error'):
            if not line:
                error = 1
                break
            if ' '.join(line.strip().split()) in lines_to_match:
              lines_matched += 1
            line = fin.readline()
        
        if lines_matched != len(lines_to_match) or error:
          errmsg = 'Failed to extract data block 1 from Bulletin B/C04 file {:} (error #{:})'.format(c04fn, error)
          raise RuntimeError(errmsg)
        
        line = fin.readline() # empty ...
        assert(line.strip() == '')

        keys = ["year", "month", "day", "mjd", "x", "y", "ut1-utc", "dX", "dY", "x_sigma", "y_sigma", "ut1-utc_sigma", "dX_sigma", "dY_sigma"]
        if not error:
            ## read actual final data
            line = fin.readline()
            while line.strip() != '':
                d = {}
                ## int values (year, month, day_of_month, mjd)
                int_list = [int(d) for d in line.split()[0:4]]
                ## floating values, should be 10
                float_list = [ float(f) for f in line.split()[4:]]
                assert(len(float_list) == 10)
                ## make dictionary
                for k,v in zip(keys, int_list+float_list): d[k]=v
                ## append collected data
                final_data.append(d)
                ## next line
                line = fin.readline()
 
    if not error:
        return final_data
    else:
        errmsg = 'Failed to extract data block 1 from Bulletin B/C04 file {:} (error #{:})'.format(c04fn, error)
        raise RuntimeError(errmsg)

def get_dlist_record_for(mjd, dlist):
  for d in dlist:
    if d['mjd'] == mjd: return copy.deepcopy(d)
  return None

def c04_get_final_data(c04fn):
  eop_data = []
  ## first collect block 1 data
  block1_eop = c04_get_final_data_01(c04fn)
  ## collect block 3 data
  block3_eop = c04_get_final_data_03(c04fn)
  ## merge based on datetime (mjd)
  for d1 in block1_eop:
    d3 = get_dlist_record_for(d1['mjd'], block3_eop)
    ## drop duplicate records from d3 (e.g. year, month, etc ...)
    for k in d1.keys(): d3.pop(k, None)
    ## append key/value pair and add to result list
    for k, v in d3.items(): d1[k] = v
    eop_data.append(d1)
  return eop_data

## Given a Bulletin B/C04 file, open it and extract the time interval 
## (start/stop datetimes) of the final data block.
def c04_date_span(c04fn):
  ## get the datetimes of the final data in the given file
  dates = [ datetime.datetime.strptime('{:}-{:02d}-{:02d}'.format(d['year'],d['month'],d['day']), '%Y-%m-%d') for d in c04_get_final_data(c04fn) ]
  ## assert they are sorted
  assert(sorted(dates) == dates)
  ## return first and last
  return dates[0], dates[-1]

## cast a dictionary of type {filename1:
##                           {'t0':datetime_start, 't':datetime_stop},
##                           {filename2:
##                           {'t0':datetime_start, 't':datetime_stop}, ...
## to  a sorted list of filenames, base on the value 't0'
def c04_dict_to_sorted_list(c04_dict):
    lst = []
    ct = datetime.datetime.max
    for fn,vals in c04_dict.items():
        if vals['t0'] < ct:
            lst = [fn] + lst
            ct = vals['t0']
        else:
            lst.append(fn)
            if len(lst) == 1: ct = vals['t0']
    return lst

## write out a list (of dicts) of IERS Bulletin B/C04 final data
def c04_cat(outfn, *c04fns):
  #with open(outfn, 'w') as fout:
  fout = open(outfn, 'w') if outfn else sys.stdout

  print(header, end='', file=fout)
  for i,cfn in enumerate(c04fns):
    fdata = c04_get_final_data(cfn)
    for j,line in enumerate(fdata):
      wlst = [ '{:4d}{:>4d}{:>4d}{:>8d}'.format(line['year'], line['month'], line['day'], line['mjd']) ]
      print(*(wlst), end='', file=fout)
      
      keys = ["x", "y", "ut1-utc", "dX", "dY", "x_sigma", "y_sigma", "ut1-utc_sigma", "dX_sigma", "dY_sigma"]
      for k in keys:
        print('{:+10.4f}'.format(line[k]), end='', file=fout)
      print(' ', end='', file=fout)
      
      keys = ["lod", "lod_sigma", "omega", "omega_sigma"]
      for k in keys:
        print('{:+15.12f} '.format(line[k]), end='', file=fout)
      
      if i==len(c04fns)-1 and j==len(fdata)-1:
        pass
      else:
        print('', file=fout) ## aka newline 
      # flst = [ '{:8.3f}{:9.3f}{:10.4f}{:9.3f}{:7.3f}{:9.3f}{:9.3f}{:10.4f}{:7.3f}{:7.3f}'.format(*line[4:])]

  if fout is not sys.stdout:
    fout.close()

def merge_c04_files(outfn, c04_dict):
  c04_list = c04_dict_to_sorted_list(c04_dict)
  c04_cat(outfn, *c04_list)

if __name__ == '__main__':

    ## parse cmd
    args = parser.parse_args()

    ## verbose print if requested
    vprint = print if args.verbose else lambda *a, **k: None

    ## datetime stamps we are going to use
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
        vprint('Downloaded file: {:}'.format(remote))

        ## get its time span
        d0,dn = c04_date_span(local)
        vprint('\tTime span: {:} to {:}'.format(d0, dn))

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
        elif prior_covered and post_covered:
          pass
        else:
           raise RuntimeError('WTF??')

        dummy_it += 1

    ## what file are we writing in?
    outfile = None if not args.save_as else args.save_as

    ## merge files to one
    merge_c04_files(outfile, files_to_merge)

    ## delete downloaded files
    #for fn in files_to_merge: os.remove(fn)

    ## all done!
    sys.exit(0)
