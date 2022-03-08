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
  def distance(self, other):
    return math.sqrt((self.x-other.x)*(self.x-other.x) + (self.y-other.y)*(self.y-other.y))

def inner_tangents(A,ra,B,rb,ud='u'):
  hyp = A.distance(B)
  short = ra+rb
  if ud == 'u':
    phi = math.atan2(B.y-A.y,B.x-A.x) + math.asin(short / hyp) - math.pi/2e0
  else:
    phi = math.atan2(B.y-A.y,B.x-A.x) - math.asin(short / hyp) + math.pi/2e0
  cp = math.cos(phi)
  sp = math.sin(phi)
  C1 = TzPoint(A.x+ra*cp, A.y+ra*sp)
  cp = math.cos(phi+math.pi)
  sp = math.sin(phi+math.pi)
  C2 = TzPoint(B.x+rb*cp, B.y+rb*sp)
  return C1,C2

def outer_tangents(A,ra,B,rb,ud='u'):
  hyp = A.distance(B)
  short = ra-rb
  if ud == 'u':
    phi = math.atan2(B.y-A.y,B.x-A.x) + math.acos(short/hyp)
  else:
    phi = math.atan2(B.y-A.y,B.x-A.x) - math.acos(short/hyp)
  sp = math.sin(phi)
  cp = math.cos(phi)
  C1 = TzPoint(A.x+ra*cp, A.y+ra*sp)
  C2 = TzPoint(B.x+rb*cp, B.y+rb*sp)
  return C1,C2

def points2line(A, B):
  slope = (B.y - A.y) / (B.x - A.x)
  coef = A.y - slope * A.x
  return slope, coef

def vertex(T1,E1,T2,E2):
  slope1, coef1 = points2line(T1,E1)
  slope2, coef2 = points2line(T2,E2)
  x = (coef2 - coef1) / (slope1-slope2)
  y1 = slope1 * x + coef1
  y2 = slope2 * x + coef2
  y = (y1+y2)/2e0
  return TzPoint(x,y)

def penumbra_proj(TS,TE,V):
  slope, coef = points2line(TS,TE)
  Vy = slope * V.x + coef
  return TzPoint(V.x, Vy)

def lineIntersectCircle(A1,B1,C,r):
  slope, coef = points2line(A1,B1)
  # ad = 1e0
  ad = 1e0 + slope*slope
  # bd = 2e0*(slope*coef-C.x-slope*C.y) / (1e0+slope*slope)
  bd = 2e0*(-C.x+slope*coef-slope*C.y)
  # gd = - (r*r - C.x*C.x-C.y*C.y-coef*coef+2e0*coef*C.y) / (1e0+slope*slope)
  gd = - (r*r-C.x*C.x-(coef-C.y)*(coef-C.y))
  x = (-bd + math.sqrt(gd)) / (2e0*ad)
  y = slope * x + coef
  sol1 = TzPoint(x,y)
  x = (-bd - math.sqrt(gd)) / (2e0*ad)
  y = slope * x + coef
  sol2 = TzPoint(x,y)
  return sol1, sol2

S = TzPoint(-1e0, 0e0); rs = 2e0
E = TzPoint(6e0, 0e0); re = 0.7e0;
OTS1, OTE1 = outer_tangents(S,rs,E,re)
OTS2, OTE2 = outer_tangents(S,rs,E,re,'d')
V = vertex(OTS1, OTE1, OTS2, OTE2)
ITS1, ITE1 = inner_tangents(S,rs,E,re)
ITS2, ITE2 = inner_tangents(S,rs,E,re, 'd')
V1 = penumbra_proj(ITS1, ITE1, V)
V2 = penumbra_proj(ITS2, ITE2, V)
SL = S - TzPoint(rs,0e0)
K = vertex(ITS1, ITE1, ITS2, ITE2)
print("\\coordinate (Sun) at {:};".format(S))
print("\\coordinate (Earth) at {:};".format(E))
print("\\coordinate (OTS1) at {:};".format(OTS1))
print("\\coordinate (OTE1) at {:};".format(OTE1))
print("\\coordinate (OTS2) at {:};".format(OTS2))
print("\\coordinate (OTE2) at {:};".format(OTE2))
print("\\coordinate (V) at {:};".format(V))
print("\\coordinate (ITS1) at {:};".format(ITS1))
print("\\coordinate (ITE1) at {:};".format(ITE1))
print("\\coordinate (ITS2) at {:};".format(ITS2))
print("\\coordinate (ITE2) at {:};".format(ITE2))
print("\\coordinate (V1) at {:};".format(V1))
print("\\coordinate (V2) at {:};".format(V2))
print("\\coordinate (K) at {:};".format(K))

print("\\draw[fill=black!5] (ITE1) -- (V1) -- (V2) -- (ITE2);")
print("\\draw[fill=black!20] (OTE1) -- (V) -- (OTE2);")

print("\\node[] at (Sun) {Sun};")
print("\\filldraw[color=yellow!60, fill=yellow!4, thick] (Sun) circle (2.0);")
print("\\filldraw[gray] (Sun) circle (1pt);")
print("\\node[] at (Earth) {Earth};")
print("\\filldraw[color=red!60, fill=red!5, very thick] (Earth) circle (0.7);")
print("\\filldraw[gray] (Earth) circle (1pt);")
print("%\\node[] at (OTS1) {OTS1};")
print("%\\filldraw[gray] (OTS1) circle (1pt);")
print("%\\node[] at (OTE1) {OTE1};")
print("%\\filldraw[gray] (OTE1) circle (1pt);")
print("%\\node[] at (OTS2) {OTS2};")
print("%\\filldraw[gray] (OTS2) circle (1pt);")
print("%\\node[] at (OTE2) {OTE2};")
print("%\\filldraw[gray] (OTE2) circle (1pt);")
print("\\node[] at (V) {V};")
print("\\filldraw[gray] (V) circle (1pt);")
print("%\\node[] at (ITS1) {ITS1};")
print("%\\filldraw[gray] (ITS1) circle (1pt);")
print("%\\node[] at (ITE1) {ITE1};")
print("%\\filldraw[gray] (ITE1) circle (1pt);")
print("%\\node[] at (ITS2) {ITS2};")
print("%\\filldraw[gray] (ITS2) circle (1pt);")
print("%\\node[] at (ITE2) {ITE2};")
print("%\\filldraw[gray] (ITE2) circle (1pt);")
print("\\draw[] (OTS1) -- (OTE1) -- (V);")
print("\\draw[] (OTS2) -- (OTE2) -- (V);")
print("%\\node[] at (V1) {V1};")
print("%\\filldraw[gray] (V1) circle (1pt);")
print("\\draw[] (ITS1) -- (ITE1) -- (V1);")
print("%\\node[] at (V2) {V2};")
print("%\\filldraw[gray] (V2) circle (1pt);")
print("\\draw[] (ITS2) -- (ITE2) -- (V2);")
print("\\draw[] {:} -- (V);".format(SL))
print("\\draw[] (V1) -- (V2);")
print("\\node[] at (K) {K};")
print("\\filldraw[gray] (K) circle (1pt);")
yputitle = (V1.y - V.y)*(2e0/3e0) + V.y
xputitle = V.x - (V.x - E.x)*(1e0/3e0)
Penumbra_tilte = TzPoint(xputitle, yputitle)
print("\\node[] at {:} {{PenUmbra}};".format(Penumbra_tilte))

print("\\pic[orange, \"${\\alpha}_{umb}$\", draw=orange, <->, angle eccentricity=1.2, angle radius=1cm]{angle = OTE1--V--Earth};")
print("\\pic[orange, \"${\\alpha}_{pen}$\", draw=orange, <->, angle eccentricity=1.2, angle radius=1cm]{angle = Earth--K--ITE1};")
Kproj1 = TzPoint(K.x, OTS2.y)
print("\\coordinate (KAumb) at {:};".format(Kproj1))
print("\\draw[thin,gray,dashed] (OTS2) -- (KAumb);")
print("\\pic[orange, \"${\\alpha}_{umb}$\", draw=orange, <->, angle eccentricity=1.2, angle radius=1cm]{angle = KAumb--OTS2--OTE2};")
Kproj2 = TzPoint(K.x, ITS1.y)
print("\\coordinate (KApen) at {:};".format(Kproj2))
print("\\draw[thin,gray,dashed] (ITS1) -- (KApen);")
print("\\pic[orange, \"${\\alpha}_{pen}$\", draw=orange, <->, angle eccentricity=1.2, angle radius=1cm]{angle = KApen--ITS1--ITE1};")
print("\\draw[thin,gray,dashed] (Sun) -- (OTS1);")
print("\\tkzMarkRightAngle[draw=gray,size=.2](OTE1,OTS1,Sun);")
print("\\draw[thin,gray,dashed] (Sun) -- (ITS1);")
print("\\tkzMarkRightAngle[draw=gray,size=.2](K,ITS1,Sun);")
print("\\draw[thin,gray,dashed] (Earth) -- (OTE1);")
print("\\tkzMarkRightAngle[draw=gray,size=.1](V,OTE1,Earth);")

RsRe1,RsRe2 = lineIntersectCircle(S,OTS1,OTS1,re)
#slope_, coef = points2line(OTS1,OTE1)
## must pass through Earth ...
#b_ = E.y - slope_*E.x
#slope, coef = points2line(S,OTS1)
#x_ = (coef - b_) / (slope_-slope)
#y_ = slope_ * x_ + b_
#RsRe1 = TzPoint(x_,y_)
print("\\draw[thin,gray,dashed] (Earth) -- {:};".format(RsRe1))
print("\\draw[thin,gray,dashed] (Earth) -- {:};".format(RsRe2))
print("\\filldraw[color=gray] {:} circle ({:});".format(S,re))
