""" Module Fitting.py

by:      Rami A. Kishek
Created: July 17, 2001

Last Modified: 7/17/01

This module contains the following fitting functions of general use:

lsqfit ... Least Squares fit of a straight line to approximate a set of pts.
"""
Fitting_version = "$Id: Fitting.py,v 1.2 2002/08/14 21:06:46 ramiak Exp $"
import Numeric

def Fittingdoc():
  import Fitting
  print Fitting.__doc__

def lsqfit(x, y):
    """ (a, b, L) = lsqfit(x, y)
        Given arrays 'x' and 'y', fits a straight line of the form
        y = ax + b in order to minimize the least squares parameter
        L = (1/N)sum[(yn - a*xn - b)^2, n=1..N]

        Note that x and y are converted to arrays of float before operation.
    """
    N = len(x); x=Numeric.array(x, 'd'); y=Numeric.array(y, 'd')
    xm = Numeric.sum(x)/N; ym = Numeric.sum(y)/N
    a = ((Numeric.sum(x*y)/N)-xm*ym)/((Numeric.sum(x*x)/N) - xm**2)
    b = ym - a*xm
    L = Numeric.sum((y - a*x - b)**2)/N
    return (a, b, L)

