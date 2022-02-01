#! /usr/bin/python

## 
##  Simple python script to check the data section of a icgem (aka '.gfc')
##+ file, and replace any FORTRAN-formatted double numerical to a C/C++
##+ based onw. Aka, change any numerical value of the type: '123.123[dD]-10'
##+ to '123.123e-10.
##  The transformation will only take place in the data section, aka the part 
##+ of the file after the string 'end_of_head' up untill EOF
##

import sys
import re

if len(sys.argv) != 2:
  print('Usage: {:} [icgem-file]'.format(sys.argv[0]))
  sys.exit(1)

with open(sys.argv[1], encoding="utf8", errors='ignore') as fin:

  ## skip the header segment ...
  while True: 
    line = fin.readline()
    print('{:}'.format(line.strip()))
    if line.lstrip().startswith('end_of_head'): break
  
  ## data segment ...
  line = fin.readline()
  while line:
    line = re.sub(r"([+-]?[0-9]?\.[0-9]+)([Dd])([+-]?[0-9]+)", r"\1e\3", line)
    print('{:}'.format(line.strip()))
    line = fin.readline()
