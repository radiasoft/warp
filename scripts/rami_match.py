""" Module rami_match.py
by: Rami A. Kishek  (based on Grote's wxy_match.py)
Created:  Sep. 27, 2000

    Last Revision: 10/17/01

CAUTION: MODULE NEEDS BE IMPORTED AT BEGINNING OF RUN BEFORE ANY STEPPING !!!
differs from wxy_match in that the starting point for matching is in between
quadrupoles when the beam is round, but the env slopes are different.

##############################################################################
# Iterate toward a matched beam using the slice code (instead of the env     #
# package).  It is currently setup to match a beam starting from the space   #
# b/w quadrupoles.  To match a beam under other conditions, change the lines #
# indicated below.                                                           #
##############################################################################
"""
from warp import *
import LinearAlgebra
import singleparticle

rami_match_version = "$Id: rami_match.py,v 1.1 2001/12/03 17:48:03 ramiak Exp $"
def rami_matchdoc():
  import rami_match
  print rami_match.__doc__

# --- Set diagnostic parameters to minimize output

# --- This is the critical routine which resets certain variables and
# --- redoes a generate.
zbeam_start = top.zbeam
time_start = top.time
it_start = top.it
def minit():
  """minit() -> Re-inits the slice package"""
  top.lprntpara = false
  top.verbosity = 0
  top.zbeam = zbeam_start
  top.zgrid = top.zbeam
  top.zgridprv = top.zbeam
  top.it = it_start
  top.time = time_start
  generate()


# --- This is the function which actually does the iteration
def match(imtch=1,s=None):
  """ match(imtch=1,s=None)
Matches the envelope using WARPxy, the beam parameters are set to the
average of the initial and final values
  - imtch=1 the number of iterations to perform
  - s=160 the number of time steps across the region to be matched
  """
  if s is None:
    try: s = 2*nstpp
    except: s = 320
  top.a0=0.5*(top.a0+2.*top.xrms[0])
  top.b0=0.5*(top.b0+2.*top.yrms[0])
  top.ap0=sign(0.5*(abs(top.ap0)+2.*top.vxrms[0]/top.vbeam),top.ap0)
  top.bp0=sign(0.5*(abs(top.bp0)+2.*top.vyrms[0]/top.vbeam),top.bp0)
  print ("a0 = %9.6f b0 = %9.6f ap0 = %9.6f bp0 = %9.6f" %
         (top.a0, top.b0, top.ap0, top.bp0))
  for i in xrange(imtch):
    minit()
    step(s)
    print "-"*70
    print ("a error = %18.13e a' error = %18.13e"%
           (2.*top.xrms[0]-top.a0,2.*top.vxrms[0]/top.vzbar[0]-abs(top.ap0)))
    print ("b error = %18.13e b' error = %18.13e"%
           (2.*top.yrms[0]-top.b0,2.*top.vyrms[0]/top.vzbar[0]-abs(top.bp0)))
    print "-"*70
    top.a0=0.5*(top.a0+2.*top.xrms[0])
    top.b0=0.5*(top.b0+2.*top.yrms[0])
    top.ap0=sign(0.5*(abs(top.ap0)+2.*top.vxrms[0]/top.vbeam),top.ap0)
    top.bp0=sign(0.5*(abs(top.bp0)+2.*top.vyrms[0]/top.vbeam),top.bp0)
    print ("a0 = %9.6f b0 = %9.6f ap0 = %9.6f bp0 = %9.6f" %
           (top.a0, top.b0, top.ap0, top.bp0))


def match1(imtch=1,s=128):
  """match1(imtch=1,s=128)
Matches the envelope using WARPxy, the X beam parameters are set to the
average of the initial and final values, then Y is set equal to X
  - imtch=1 the number of iterations to perform
  - s=128 the number of time steps across the region to be matched
  """
  for i in xrange(imtch):
    top.a0=0.5*(top.a0+2.*top.xrms[0]) #xrms[0]+yrms[0]
    top.b0=top.a0
    top.ap0=sign(0.5*(abs(top.ap0)+2.*top.vxrms[0]/top.vbeam),top.ap0)
    top.bp0=sign(0.5*(abs(top.bp0)+2.*top.vyrms[0]/top.vbeam),top.bp0)
    minit()
    step(s)
    print ("a error = %18.13e a' error = %18.13e"%
           (2.*top.xrms[0]-top.a0,2.*top.vxrms[0]/top.vzbar[0]-top.ap0))
    print ("b error = %18.13e b' error = %18.13e"%
           (2.*top.yrms[0]-top.b0,2.*top.vyrms[0]/top.vzbar[0]+top.bp0))


# --- This is an unrelated function which can be used to calculate sigma0
# --- from the particle in the full applied field.  It uses the calc_sig0
# --- script and resets several variables and turns the moments
# --- calculation off to reduce computation time.
#from calc_sig0 import *
#def go():
#  top.itmomnts = 0
#  top.zbeam = 0.;top.zgrid = 0.;top.it = 0;top.time = 0
#  setlatt()
#  calc_sig0()

########################################################################
def matchx(xf=0.,xpf=0.,yf=0.,ypf=0.,zs=None,ze=None,s=None,
           maxiter=100,tol=1.e-10,vary=0.01,
	   xi=0.,xpi=0.,yi=0.,ypi=0.,xs=1.,xps=1.,ys=1.,yps=1.):
  """matchx(xf=0.,xpf=0.,yf=0.,ypf=0.,zs=None,ze=None,s=None,
           maxiter=100,tol=1.e-10,vary=0.01,
	   xi=0.,xpi=0.,yi=0.,ypi=0.,xs=1.,xps=1.,ys=1.,yps=1.)
    Matches a single particle orbit to the specified values
     - xf=0. final value of x
     - xpf=0. final value of x'
     - yf=0. final value of y
     - ypf=0. final value of y'
     - zs=env.zl starting value of z
     - ze=env.zu final value of z
     - s=env.nenv number of time steps across region of interest
     - maxiter=100 maximum number of iterations
     - tol=1.e-10 tolerance position is matched to"""
  mat = zeros((4,4),'d')
  dsp = zeros(4,'d')

  # --- Get number of steps if not specified
  if not s:
    if not zs:
      zs = env.zl
    if not ze:
      ze = env.zu
    s = nint((ze - zs)/wxy.ds)

  notdone = 1
  iter = 0
  while notdone and iter < maxiter:
    iter = iter + 1

    xx = top.x0
    yy = top.y0
    vx = top.xp0*top.vbeam
    vy = top.yp0*top.vbeam
    singleparticle.spinit([xx,xx*(1. + vary),xx,xx,xx],
                          [yy,yy,yy,yy*(1. + vary),yy],
                          [zs,zs,zs,zs,zs],
                          [vx,vx,vx*(1. + vary),vx,vx],
                          [vy,vy,vy,vy,vy*(1. + vary)],
                          5*[top.vbeam],5*[1.],s)
    step(s)

    for ip in xrange(4):
      mat[0,ip]=(top.xp[ip+1]  - top.xp[0])/((xx-xi)/xs*vary)
      mat[1,ip]=(top.uxp[ip+1] - top.uxp[0])/((vx-xpi*top.vbeam)/xps*vary)
      mat[2,ip]=(top.yp[ip+1]  - top.yp[0])/((yy-yi)/ys*vary)
      mat[3,ip]=(top.uyp[ip+1] - top.uyp[0])/((vy-ypi*top.vbeam)/yps*vary)

    # --- Invert the matrix
    mati = LinearAlgebra.inverse(mat)

    # --- Get next starting point
    dsp[:] = [xf - top.xp[0],xpf - top.uxp[0]/top.uzp[0],
              yf - top.yp[0],ypf - top.uyp[0]/top.uzp[0]]
    top.x0 = top.x0 + sum(mati[0,:]*dsp)*xs
    top.xp0 = top.xp0 + sum(mati[1,:]*dsp)*xps
    top.y0 = top.y0 + sum(mati[2,:]*dsp)*ys
    top.yp0 = top.yp0 + sum(mati[3,:]*dsp)*yps

    # --- Check for convergence
    print "error => x = %10.3e x' = %10.3e"%tuple(dsp[:2])
    print "error => y = %10.3e y' = %10.3e"%tuple(dsp[2:])
    if max(abs(dsp)) < tol : notdone = 0

  if iter == maxiter:
    print 'Warning: Maximum number of iterations reached'
  else:
    print "top.x0 = %20.15e;top.xp0 = %20.15e"%(top.x0,top.xp0)
    print "top.y0 = %20.15e;top.yp0 = %20.15e"%(top.y0,top.yp0)

