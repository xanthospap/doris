#! /usr/bin/python

import sys
C = 299792458e0

def foo(vm, vt, Ndop, dt, fen, frt):
  coef = C * (Ndop + frt) / fen
  return 1e0 * (vm - vt) / coef

sats = sys.argv[1]

with open('foo', 'r') as fin:
  for line in fin.readlines():
      l = line.split()
      if l[0] == sats:
        l2 = float(l[1].split("=")[1])
        l1 = float(l[2].split("=")[1])
        dt = float(l[3].split("=")[1])
        ion= float(l[4].split("=")[1])
        tro= float(l[5].split("=")[1])
        rel= float(l[6].split("=")[1])
        dfe= float(l[7].split("=")[1])
        fen= float(l[8].split("=")[1])
        frt= float(l[9].split("=")[1])
        rt2= float(l[10].split("=")[1])
        rt1= float(l[11].split("=")[1])

        Vm = (C/fen) * (fen - frt - (l2-l1)/dt)
        Vt = (rt2 - rt1) / dt
        cor = ion + rel # + tro
        res1 = Vm + Vt
        res2 = Vm + Vt + cor
        x = foo(Vm, Vt, l2-l1, dt, fen, frt)

        print('{:+9.3f}  {:+9.3f}  {:+7.3f} {:+7.3f}'.format(Vm, Vt, res1, res2))
