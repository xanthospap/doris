#! /usr/bin/python

import numpy as np
import datetime
from numpy.linalg import inv
import matplotlib.pyplot as plt
import sys, os

def parse_xy(fn):
    xl = []
    yl = []
    with open(fn, 'r') as fin:
        for line in fin.readlines():
            l = line.split()
            x, y = [float(x) for x in line.split()[0:2]]
            xl.append(x)
            yl.append(y)
    return xl, yl


def make_model(tstart, tstop, tref, state, sigma_noise):
    model = []
    with open('model.dat', 'w') as fout:
        while tstart < tstop:
            dt = (tstart-tref).total_seconds() / \
                datetime.timedelta(days=1).total_seconds()
            val = state[0] + state[1] * dt + \
                np.random.normal(0e0, sigma_noise*sigma_noise)
            model.append((tstart, val[0]))
            #print(val, val.shape)
            print('{:.8f} {:.8f}'.format(dt, val[0]), file=fout)
            tstart += datetime.timedelta(0, 30)
    return model


class LinearKalman:
    ## time in MJD
    def __init__(self, x0, P0, tref, sigma_noise, process_noise=np.zeros((2,2))):
        self.x = x0
        self.P = P0
        self.tref = tref
        self.R = sigma_noise * sigma_noise
        self.Q = process_noise
        ## for debugging ...
        self.max_k0_dif = -1e-15
        self.max_k1_dif = -1e-15
        self.xdeb = x0
        self.Pdeb = P0
    
    def recreate(self, t1, t2, step, fout=sys.stdout):
        t = []
        v = []
        while t1<t2:
            val = (t1-self.tref)*self.x[1] + self.x[0]
            print('{:.8f} {:.8f}'.format(t1, val[0,0]), file=fout)
            t.append(t1)
            v.append(val[0,0])
            t1 += step
        return t, v

    def update(self, z, t, fout=None):
        # F = np.mat([[1e0, t-self.tref],[0e0, 1e0]])
        F = np.eye(2)
        # H = np.mat([1e0, 0e0])
        H = np.mat([1e0, t-self.tref])
        newx = F * self.x
        newP = F * self.P * F.transpose() + self.Q

        y = z - H*newx
        S = H*newP*H.transpose() + self.R
        K = newP * H.transpose() * inv(S)
        self.x = newx + K*y
        self.P = (np.eye(2)-K*H)*newP
    
        ## Debugging version ...    
        dt = t-self.tref
        p0p1dt = self.Pdeb[0,0] + self.Pdeb[0,1]*dt
        p1p2dt = self.Pdeb[0,1] + self.Pdeb[1,1]*dt
        s = p0p1dt + dt*p1p2dt + self.R
        k0 = p0p1dt / s
        k1 = p1p2dt / s
        if abs(k0-K[0,0]) > self.max_k0_dif: self.max_k0_dif = abs(k0-K[0,0])
        if abs(k1-K[1,0]) > self.max_k1_dif: self.max_k1_dif = abs(k1-K[1,0])
        self.xdeb[0,0] += k0 * y[0,0]
        self.xdeb[1,0] += k1 * y[0,0]
        #if abs(x0 - self.x[0,0]) > 1e12:
        #  print('>> Oops!! diff in x0 = {.12f}'.format(abs(x0 - self.x[0,0])))
        #if abs(x1 - self.x[1,0]) > 1e10:
        #  print('>> Oops!! diff in x1 = {.12f}'.format(abs(x1 - self.x[1,0])))
        p0 = (1e0-k0)*self.Pdeb[0,0] - self.Pdeb[0,1]*k0*dt
        p1 = ( ((1e0-k0)*self.Pdeb[0,1] - self.Pdeb[1,1]*k0*dt) +
          (-k1*self.Pdeb[0,0] + (1e0-k1*dt)*self.Pdeb[0,1]) ) / 2e0
        p2 = -k1*self.Pdeb[0,1] + (1e0-k1*dt) * self.Pdeb[1,1]
        self.Pdeb[0,0] = p0
        self.Pdeb[0,1] = self.Pdeb[1,0] = p1
        self.Pdeb[1,1] = p2
        #print(">> dt={:.10f} k0={:.8f} k1={:.8f} x0={:.8f} x1={:.8f} P0={:.8f} P1={:.8f} P2={:.8f}".format(dt, k0, k1, self.xdeb[0,0], self.xdeb[1,0], self.Pdeb[0,0], self.Pdeb[0,1], self.Pdeb[1,1]))

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print('Usage test_kalman [FILE]', file=sys.stderr)
        sys.exit(1)

    if not os.path.isfile(sys.argv[1]):
        print('ERROR File {:} does not exist!'.format(
            sys.argv[1]), file=sys.stderr)
        sys.exit(1)
    
    x, y = parse_xy(sys.argv[1])
    filter = LinearKalman(np.mat([[y[0]], [5e0]]), 
      np.mat([[2.5e0, 0e0],[0e0, 2e0]]), 
      x[0], 5e0
      #np.mat([[2e0, 0e0],[0e0, 1e0]])
      )

    for xy in zip(x,y):
        filter.update(xy[1], xy[0])
    
    with open('kmodel.dat', 'w') as fout:
        t, v = filter.recreate(x[0],x[-1],0.1e0,fout)

    plt.plot(x, y, "o", t, v)
    plt.show()

    sigmas = np.sqrt(filter.P)
    sigmasDeb = np.sqrt(filter.Pdeb)
    print('Estimates:\n\tx0 = {:.4f} +/- {:.5f}\n\tx1 = {:.4f} +/- {:.5f}'.format(
        filter.x[0,0], sigmas[0, 0], filter.x[1,0], sigmas[1, 1]))
    print('Estimates:\n\tx0 = {:.4f} +/- {:.5f}\n\tx1 = {:.4f} +/- {:.5f}'.format(
        filter.xdeb[0,0], sigmasDeb[0, 0], filter.xdeb[1,0], sigmasDeb[1, 1]))
    
    print('Implementation details:')    
    print('>> Max Diffs: K0-k0={:.12f} K1-k1={:.12f}'.format(
    filter.max_k0_dif, filter.max_k1_dif))
