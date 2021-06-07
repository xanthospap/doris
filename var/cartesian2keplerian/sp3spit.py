#! /usr/bin/python

import sys

if len(sys.argv) != 2:
  print('Usage: {:} [SP3_FILE]'.format(sys.argv[0]));
  sys.exit(1)

with open(sys.argv[1], 'r') as sp3:
  for line in sp3.readlines():
    if line[0] in 'PV' and line[1] in 'LGRES' and line[2].isdigit() and line[3].isdigit():
      if line[0] == 'P':
        x, y, z = [ float(b)*1e3 for b in line[4:46].split() ]
        sv_pos = line[1:4]
        position_ok = True
      elif line[0] == 'V':
        vx, vy, vz = [ float(b)*1e-1 for b in line[4:46].split() ]
        sv_vel = line[1:4]
        velocity_ok = True
        if position_ok and sv_vel == sv_pos:
          print('{:+15.3f} {:+15.3f} {:+15.3f} {:+12.6f} {:+12.6f} {:+12.6f}'.format(x, y, z, vx, vy, vz))
      else:
        position_ok = False
        velocity_ok = False
