"""Generates Thermal KV distribution
Module ParaKV.py

 by: Rami A. Kishek
 Date: June 5, 2001

    Last Modified: 3/14/2002

 This module contains functions to generate a "thermal" KV distribution,
 i.e., a semiGaussian with a Parabolic Gaussian distribution.  The function
 can also be used to adjust the initial temperature distribution via the
 parameter "delta_hp".  Note that in either case, the variable distrbtn must
 be set to "semigauss". The functions need be called after the "generate"
 but before stepping through the lattice.

 The formula used is :

            T(r) = To*C*[1 - delta_hp*(r/R0)^2]

 The variable delta_hp can vary from -inf to 1.0. Zero is a semi-Gaussian.
 0..1.0 gives a bulged profile.  Negative values give a hollow profile.
 delta_hp=1.0 corresponds to the thermal KV.

 The resultant distribution is normalized to produce the correct rms beam
 size and emittance.  The normalization constant C is the scaling at the
 beam center and is equal to 2/(2-delta_hp).

 10/03/01:  Modified the normalization to allow for a nonuniform beam
            in space.
 3/14/02:   Added functions 'dualgauss' and 'loadvels'

 Three functions are provided:

 para_temp ...  A general function to use for producing a parabolic temperature
                profile
 thermal_kv ... calls para_temp with the right parameter for a zero temperature
                at the edge.
 test_T_prof ... function to check temperature profile.
 dualgauss  ... generates a dual-Guassian profile
 loadvels   ... loads a given temperature distribution
"""

from warp import *

ParaKV_version = "$Id: ParaKV.py,v 1.3 2003/08/07 19:00:57 dave Exp $"
def ParaKVdoc():
  import ParaKV
  print ParaKV.__doc__

def rad(xxx, yyy):
    """ rad(xxx,yyy) -> radius r = sqrt(x^2+y^2) """
    return sqrt((xxx**2) + (yyy**2))

def para_temp(delta_hp=0.3, delta_h=None):
    """ para_temp(delta_hp=0.3, delta_h=w3d.hollow_h)

        Produces a beam with a parabolic temperature profile, according
        to the formula:

            T(r) = To*C*[1 - delta_hp*(r/R0)^2]

        The variable delta_hp can vary from -inf to 1.0. Zero is a
        semi-Gaussian.  0..1.0 gives a bulged profile.  Negative values
        give a hollow profile.  delta_hp=1.0 corresponds to the thermal KV.

        The resultant distribution is normalized to produce the correct
        rms beam  size and emittance.  The normalization constant
        C is the scaling at the beam center and is equal
        to 2/(2-delta_hp) for a uniform beam.  For nonuniform beam, parameter
        w3d.hollow_h is used to adjust normalization.

        Should be called after the generate, and the variable 'w3d.distrbtn'
        must be set to "semigauss".
    """
    if delta_h is None: delta_h = w3d.hollow_h
    def parabolic(radius, hhp, delta_h):
        """ parabolic(radius)
            returns scaling factor for a particular radius which results in
            the parabolic profile
        """
        hh = 1.0-(1.0/delta_h)
        c_scale = (1-0.5*hh)/(1.0-0.5*(hh+hhp)+(hh*hhp/3.0))     # Normalization factor
        return sqrt( c_scale*(1 - hhp*(radius**2)) )

    # --- Adjust particle velocities after first subtracting flow term
    #     i.e., scale only the thermal part of the velocities.

    top.uxp = top.uzp*(top.ap0/top.a0)*top.xp + (top.uxp -
                top.uzp*(top.ap0/top.a0)*top.xp)*parabolic(
                    rad(top.xp/top.a0, top.yp/top.b0), delta_hp, delta_h)
    top.uyp = top.uzp*(top.bp0/top.b0)*top.yp + (top.uyp -
                top.uzp*(top.bp0/top.b0)*top.yp)*parabolic(
                    rad(top.xp/top.a0, top.yp/top.b0), delta_hp, delta_h)


def thermal_kv():
    """ thermal_kv()

        Produces a "thermal KV" beam with a parabolic temperature profile
        that goes to zero at the beam edge.
        Should be called after the generate, and the variable 'w3d.distrbtn'
        must be set to "semigauss".
    """
    para_temp(delta_h=1.0)


def dualgauss(v, vm, vs):
    """ dualgauss(v, vm, vs)
        Profile of dual-Gausssian in velocity space.
    """
    return 2.0*exp(-0.5*(vm/vs)**2)*exp(-0.5*(v/vs)**2)*cosh(v*vm/(vs**2))


def interp(temp, F, v):
    """ interp(temp, F, v)
        finds values of v corresponding to value 'temp' from given F(v)
        using LINEAR interpolation
    """
    i1 = searchsorted(F, temp)
    return take(v,i1-1)+(temp-take(F,i1-1))*(take(v,i1)-take(v,i1-1))/(take(F,i1)-take(F,i1-1))


def loadvels(func=dualgauss, vm=0.02, vs=0.006, cut=0.038, npts=1000):
    """ loadvels(func=dualgauss, vm=0.02, vs=0.006, cut=0.038, npts=1000)
        Generates a 2-D random velocity distribution top.uxp, top.uyp that
        is uncorrelated in space and fits a profile in velocity space given
        by the function 'func'.
        The velocity is discretized by 'npts' points from 0 up to 'cut'.
    """
    dv = cut/npts
    v = arange(0, cut+dv, dv)
    f = array(map(func, v, vm*ones(v.shape), vs*ones(v.shape)))
    F = cumsum(2*pi*v*f)        # Calculate the c.d.f.
    F = F/F[-1];  f = f/F[-1]   # Normalize
    # -- Generate random arrays for creating velocities
    temp = ranf(top.uxp);  rvang = 2*pi*ranf(top.uxp)
    # -- Interpolate 'temp' in F to find total velocities
    rvt = interp(temp, F, v)
    # -- Set Particle Velocities
    top.uxp = top.uzp*rvt*cos(rvang)
    top.uyp = top.uzp*rvt*sin(rvang)


def test_T_prof(ncpb=2):
    """ test_T_prof(ncpb=2)

        Plots temperature profile Tx, Ty,  vs x and y, as well as
        particle plot of phase space.  Useful to verify correct operation
        and scaling of other functions in ParaKV.py

        ncpb = Number of Cells Per Bin
    """
    def findtemp(vels):
        npart = len(vels)
        if npart == 0:
            return 0
        else:
            return sqrt((sum(vels**2)/npart) - (sum(vels)/npart)**2)

    xcent = ycent = 0
    if w3d.l2symtry: xcent = w3d.nx/2
    elif w3d.l4symtry: pass
    else:
        xcent = w3d.nx/2
        ycent = w3d.ny/2

    temps_x = zeros([w3d.nx/ncpb,3],'d')
    temps_y = zeros([w3d.ny/ncpb,3],'d')
    for inx in range(0,w3d.nx/ncpb):
        temps_x[inx,0] = findtemp(getxp(iy=ycent, ix=inx*ncpb, wx=ncpb))**2
        temps_x[inx,1] = findtemp(getyp(iy=ycent, ix=inx*ncpb, wx=ncpb))**2
        temps_x[inx,2] = findtemp(getrp(iy=ycent, ix=inx*ncpb, wx=ncpb))**2
    for iny in range(0,w3d.ny/ncpb):
        temps_y[iny,0] = findtemp(getxp(ix=xcent, iy=iny*ncpb, wy=ncpb))**2
        temps_y[iny,1] = findtemp(getyp(ix=xcent, iy=iny*ncpb, wy=ncpb))**2
        temps_y[iny,2] = findtemp(getrp(ix=xcent, iy=iny*ncpb, wy=ncpb))**2
    plg(temps_x[:,0]); ptitles("Temperature Profile", "X", "Tx"); fma()
    plg(temps_x[:,1]); ptitles("Temperature Profile", "X", "Ty"); fma()
    plg(temps_x[:,2]); ptitles("Temperature Profile", "X", "Tr"); fma()
    plg(temps_y[:,0]); ptitles("Temperature Profile", "Y", "Tx"); fma()
    plg(temps_y[:,1]); ptitles("Temperature Profile", "Y", "Ty"); fma()
    plg(temps_y[:,2]); ptitles("Temperature Profile", "Y", "Tr"); fma()
    palette("gray.gp")
    pptrace(filled=1, particles=0, contours=7); fma()

