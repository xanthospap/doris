#! /usr/bin/python

import sys
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
from datetime import datetime
import julian

##
## sys.argv[1] -> an IERS C04 EOP file
## sys.argv[2] -> regularized values for EOPS as 
## 'mjd float dut1 float lod float'
## (file can be easily constructed using complete_eop_example.cpp)
##

tstart = datetime.strptime('2021-12-28', '%Y-%m-%d')
tend = datetime.strptime('2022-01-03', '%Y-%m-%d')

dts = []
dxps = []
dyps = []
dut1s = []
dlods = []
with open(sys.argv[1], 'r') as fin:
  for line in fin.readlines():
    if len(line)>100 and line[0] != ' ' and line[0] != '#':
      l = line.split()
      t = datetime.strptime('{:} {:02d} {:02d}'.format(int(l[0]), int(l[1]), int(l[2])), '%Y %m %d')
      xp, yp = float(l[4]), float(l[5])
      dut1, lod = float(l[6]), float(l[7])
      if t >= tstart and t< tend:
        dts.append(t)
        dut1s.append(dut1)
        dlods.append(lod)
        dxps.append(xp)
        dyps.append(yp)

dts2 = []
dut1s2 = []
dlods2 = []
dxps2 = []
dyps2 = []
with open(sys.argv[2], 'r') as fin:
  for line in fin.readlines():
    l = line.split()
    t = julian.from_jd(float(l[1]), fmt='mjd')
    if t >= tstart and t< tend:
        dts2.append(t)
        dxps2.append(float(l[3]))
        dyps2.append(float(l[5]))
        dut1s2.append(float(l[7]))
        dlods2.append(float(l[9]))

###################################################################### ERPs
## plot
f, axarr = plt.subplots(2, sharex=True)
axarr[0].scatter(dts2, dut1s2, s=0.5, c='red', alpha=0.8)
axarr[0].scatter(dts, dut1s, s=50, c='black', marker='x')
axarr[0].set_ylabel('UT1-UTC [seconds]')
axarr[0].set_title('IERS C04 EOPs & Interpolated')
axarr[1].scatter(dts2, dlods2, s=0.5, c='red', alpha=0.8)
axarr[1].scatter(dts, dlods, s=50, c='black', marker='x')
axarr[1].set_ylabel('LOD [seconds/day]')

## date tick-marks
myFmt = mdates.DateFormatter("%m/%d:%H")
axarr[0].xaxis.set_major_formatter(myFmt)
axarr[0].xaxis.set_major_locator(mdates.HourLocator(interval=12))
axarr[0].xaxis.set_minor_locator(mdates.HourLocator(interval=6))

## Rotate date labels automatically
f.autofmt_xdate()
#plt.show()
plt.savefig('erp-daily.png', transparent=True, bbox_inches='tight', dpi=300)
#sys.exit(0)

###################################################################### EOPs
## plot
f, axarr = plt.subplots(2, sharex=True)
axarr[0].scatter(dts2, dxps2, s=0.5, c='red', alpha=0.8)
axarr[0].scatter(dts, dxps, s=50, c='black', marker='x', alpha=0.8)
axarr[0].set_ylabel('x-pole [arcseconds]')
axarr[0].set_title('IERS C04 EOPs')
axarr[1].scatter(dts2, dyps2, s=0.5, c='red', alpha=0.8)
axarr[1].scatter(dts, dyps, s=50, c='black', marker='x', alpha=0.8)
axarr[1].set_ylabel('y-pole [arcseconds]')

## date tick-marks
myFmt = mdates.DateFormatter("%m/%d:%H")
axarr[0].xaxis.set_major_formatter(myFmt)
axarr[0].xaxis.set_major_locator(mdates.HourLocator(interval=12))
axarr[0].xaxis.set_minor_locator(mdates.HourLocator(interval=6))

## Rotate date labels automatically
f.autofmt_xdate()
#plt.show()
plt.savefig('eop-daily.png', transparent=True, bbox_inches='tight', dpi=300)
