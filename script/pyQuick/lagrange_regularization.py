#! /usr/bin/python

## for n = 1:N (inclusive)
##     for m = 1:M (inclusive)
## factor = sqrt[ (n-m)! (2n+1) (2-Kronecker(0,m) / (n+m)! ]
##

from math import factorial
from math import sqrt 

def Kronecker(n1,n2): return int(n1 == n2)

def raw(n,m):
  return sqrt(factorial(n-m) * (2*n+1) * (2-Kronecker(m,0)) / factorial(n+m))

def recursive(n,m):
  sqrt2_ = sqrt(2e0)
  for order in range(1, n+1):
    tnp1 = 2 * order + 1
    factor = sqrt(tnp1)
    print('F({:},{:}) = {:.12f}, '.format(order,0,factor),end='')
    factor *= sqrt2_
    for degree in range(1, order+1):
      #factor *= (1e0/(order+degree)) * (1e0/(order-(degree-1)))
      factor /= sqrt((order+degree))
      factor /= sqrt((order-(degree-1)))
      print('F({:},{:}) = {:.12f}, '.format(order,degree,factor),end='')
    print()


def results(func, n, m):
  for order in range(1, n+1):
    for degree in range(0, order+1):
      print('F({:},{:}) = {:.12f}, '.format(order,degree,func(order,degree)),end='')
    print()

if __name__ == "__main__" :
  results(raw,5,5)
  recursive(5,5)
