"""Pepper-pot style diagnostics
Module PepperGrid.py
by Rami A. Kishek
created: Mar. 7, 2002

Last Modified:  03/08/02

This module contains functions which can apply a variety of masks,
including a pepper-pot array, and grids with square holes.  The masks
are applied by removing particles in the 2-D code, but the boundary
conditions are not applied self-consistently.

Functions Available:
apply_pp    ... Pepper-pot
apply_grid  ... Rectangular Grid

"""
from warp import *

PepperGrid_version = "$Id: PepperGrid.py,v 1.4 2007/12/20 00:36:40 dave Exp $"
def PepperGriddoc():
  import PepperGrid
  print PepperGrid.__doc__

def apply_pp(h_size=0.0002, h_sep=0.002, peakx=None, peaky=None):
    """ apply_pp(h_size=0.0002, h_sep=0.002, peakx=None, peaky=None)
        Applies a rectangular pepper-pot mask with circular holes of
        radius 'h_size' and separation 'h_sep'.
        'peakx' and 'peaky' define the maximum extent of the mask from
        axis (default to xmmax and ymmax).
    """
    if peakx is None: peakx = w3d.xmmax
    if peaky is None: peaky = w3d.ymmax
    nix = nint(peakx/h_sep)             # Number of holes in x, y
    niy = nint(peaky/h_sep)
    tag = zeros([top.pgroup.npmax],'l')            # Array to tag particles
    for xh in range(0, nix+1):
        for yh in range(0, niy+1):
            # tag particles inside any pepper-pot hole
            tag = where( int(((abs(top.pgroup.xp)-xh*h_sep)**2 +
                     (abs(top.pgroup.yp)-yh*h_sep)**2)/(h_size**2)), tag, 1 )
    # kill all particles with tag unset (i.e., outside holes)
    top.pgroup.uzp = where( tag, top.pgroup.uzp, 0.0 )


def apply_grid(xwid=0.0001, ywid=0.0001, xsep=0.0005, ysep=0.0005, peakx=None, peaky=None):
    """ apply_grid(xwid=0.0001, ywid=0.0001, xsep=0.0005, ysep=0.0005, peakx=None, peaky=None)
        Applies a grid with rectangular holes of
        dimensions 'xsep' by 'ysep' and of thickness 'xwid' and 'ywid'.
        'peakx' and 'peaky' define the maximum extent of the mask from
        axis (default to xmmax and ymmax).
    """
    if peakx is None: peakx = w3d.xmmax
    if peaky is None: peaky = w3d.ymmax
    nix = nint(peakx/xsep)             # Number of lines in x, y
    niy = nint(peaky/ysep)
    # --- Remove particles hidding any gridlines
    for xh in range(0, nix+1):
        top.pgroup.uzp = where(int(abs(abs(top.pgroup.xp)-xh*xsep)/(0.5*xwid)),
                               top.pgroup.uzp, 0.0 )
    for yh in range(0, niy+1):
        top.pgroup.uzp = where(int(abs(abs(top.pgroup.yp)-yh*ysep)/(0.5*ywid)),
                               top.pgroup.uzp, 0.0 )

