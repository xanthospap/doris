#! /usr/bin/python

import math

class TzPoint:
  def __init__(self,x=0e0,y=0e0):
    self.x = x
    self.y = y
  def __add__(self, other):
    return TzPoint(self.x+other.x, self.y+other.y)
  def __sub__(self, other):
    return TzPoint(self.x-other.x, self.y-other.y)
  def __str__(self):
    return "({:.4f}, {:.4f})".format(self.x, self.y)

def outer_tangents(A,ra,B,rb):
  beta = math.asin((ra-rb) / ((A.x-B.x)*(A.x-B.x) + (A.y-B.y)*(A.y-B.y)))
  gamma = math.atan2(A.y-B.y, A.x-B.x)
  alpha = gamma - beta
  sa = math.sin(alpha)
  ca = math.cos(alpha)
  x31 = A.x + ra * sa
  y31 = A.y + ra * ca
  x32 = A.x - ra * sa
  y32 = A.y - ra * ca
  x41 = B.x + rb * sa
  y41 = B.y + rb * ca
  x42 = B.x - rb * sa
  y42 = B.y - rb * ca
  return TzPoint(x31,y31), TzPoint(x32,y32), TzPoint(x41,y41), TzPoint(x42,y42)

def points2line(A, B):
  slope = (B.y - A.y) / (B.x - A.x)
  coef = A.y - slope * A.x
  return slope, coef

def vertex(T1,T2):
  slope, coef = points2line(T1,T2)
  xv1 = -coef / slope
  if E1 is not None:
    slope, coef = points2line(E1,E2)
    xv2 = -coef / slope
    xv1 = (xv1+xv2) / 2e0
  return TzPoint(xv1,0e0)

def penumbra_up(S, rs, E, re, Vx):
  S1 = S + TzPoint(0e0, rs)
  S2 = E - TzPoint(0e0, re)
  slope, coef = points2line(S1, S2)
  Vy = slope*Vx + coef
  print('Penumbra Up at: {:.3f}, {:.3f}'.format(Vx, Vy))
  return TzPoint(Vx, Vy)

def penumbra_down(S, rs, E, re, Vx):
  S1 = S - TzPoint(0e0, rs)
  S2 = E + TzPoint(0e0, re)
  slope, coef = points2line(S1, S2)
  Vy = slope * Vx + coef
  print('Penumbra Down at: {:.3f}, {:.3f}'.format(Vx, Vy))
  return TzPoint(Vx, Vy)

S = TzPoint(-1e0, 0e0); rs = 2e0
print("\\coordinate (Sun) at {:}".format(S))
E = TzPoint(6e0, 0e0); re = 0.7e0;
print("\\coordinate (Earth) at {:}".format(E))
TS1, TS2, TE1, TE2 = outer_tangents(S,rs,E,re)
print("\\coordinate (TOEU) at {:}".format(TE1))
print("\\coordinate (TOED) at {:}".format(TE2))
print("\\coordinate (TOSU) at {:}".format(TS1))
print("\\coordinate (TOSD) at {:}".format(TS2))
V = vertex(TS1, TE1, TS2, TE2)
print("\\coordinate (V) at {:}".format(V))

PU = penumbra_up(S,rs,E,re,V.x)
PD = penumbra_down(S,rs,E,re,V.x)
