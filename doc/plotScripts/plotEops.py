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

tstart = datetime.strptime('2000-01-01', '%Y-%m-%d')
tend = datetime.strptime('2022-01-01', '%Y-%m-%d')

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
axarr[0].scatter(dts, dut1s, s=1e0, c='black', alpha=0.8)
axarr[0].scatter(dts2, dut1s2, s=1e0, c='red', alpha=0.8)
axarr[0].set_ylabel('UT1-UTC [seconds]')
axarr[0].set_title('IERS C04 EOPs Vs \'Regularized\'')
axarr[1].scatter(dts, dlods, s=1e0, c='black', alpha=0.8)
axarr[1].scatter(dts2, dlods2, s=1e0, c='red', alpha=0.8)
axarr[1].set_ylabel('LOD [seconds/day]')

## date tick-marks
myFmt = mdates.DateFormatter("%Y")
axarr[0].xaxis.set_major_formatter(myFmt)
fmt_half_year = mdates.MonthLocator(interval=24)
axarr[0].xaxis.set_major_locator(fmt_half_year)
fmt_month = mdates.MonthLocator(interval=6)
axarr[0].xaxis.set_minor_locator(fmt_month)

## Rotate date labels automatically
f.autofmt_xdate()
#plt.show()
plt.savefig('erp-00to22.png', transparent=True, bbox_inches='tight', dpi=300)

###################################################################### EOPs
## plot
f, axarr = plt.subplots(2, sharex=True)
axarr[0].scatter(dts, dxps, s=1e0, c='black', alpha=0.8)
axarr[0].scatter(dts2, dxps2, s=1e0, c='red', alpha=0.8)
axarr[0].set_ylabel('x-pole [arcseconds]')
axarr[0].set_title('IERS C04 EOPs')
axarr[1].scatter(dts, dyps, s=1e0, c='black', alpha=0.8)
axarr[1].scatter(dts2, dyps2, s=1e0, c='red', alpha=0.8)
axarr[1].set_ylabel('y-pole [arcseconds]')

## date tick-marks
myFmt = mdates.DateFormatter("%Y")
axarr[0].xaxis.set_major_formatter(myFmt)
fmt_half_year = mdates.MonthLocator(interval=24)
axarr[0].xaxis.set_major_locator(fmt_half_year)
fmt_month = mdates.MonthLocator(interval=6)
axarr[0].xaxis.set_minor_locator(fmt_month)

## Rotate date labels automatically
f.autofmt_xdate()
#plt.show()
plt.savefig('eop-00to22.png', transparent=True, bbox_inches='tight', dpi=300)
