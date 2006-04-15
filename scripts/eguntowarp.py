"""Inject particles into a simulation, based on egun results.
"""
from warp import *
import getzmom

eguntowarp_version = '$Id: eguntowarp.py,v 1.3 2006/04/15 00:13:37 dave Exp $'
def eguntowarpdoc():
  import eguntowarp
  print eguntowarp.__doc__

class EgunToWarp:
  """
EgunToWarp
  Inject particles into a simulation, basing the distribution off of an
  inputted list of EGUN trajectories. To use it, all that needs to be done is
  to create an instance of the class. There is an additional method, reset,
  which allows resetting the list of EGUN trajectories used.
 
  Just before a field solve is done, the particles are created and their
  charge density is deposited. Then, at the end of the step, the new
  particles are made live. At that time, their position and velocities are
  synchronized in time.

  Creator arguments:
  - r Array contains radius of egun trajectories
  - vr Array contains radial velocity of egun trajectories
  - vz Array contains axial velocity of egun trajectories
  - ibeam Array contains current of egun trajectories
  - slice=false When true, creates only a single slice of particles
  - npmax=1000000 Amount of space allocated for particles
  - seed1=1341392  Random number seed
  - seed2=9823748  Random number seed

  The class also relies on several WARP variables, including:
  top.zion, top.ibeam, top.dt, top.npinject, top.aion
  """

  def __init__(s,r,vr,vz,ibeam,slice=false,npmax=1000000,
                 seed1=1341392,seed2=9823748):
    # --- Allocate space for the particles
    if top.npmax == 0: top.npmax = npmax
    alotpart()
    top.ins[0] = top.npmax+1
    # --- Set the basic particle parameters
    top.sq[0] = top.zion*echarge
    top.sw[0] = (top.ibeam*top.dt/(echarge*top.zion))/top.npinject
    top.sm[0] = top.aion*amu
    # --- Now save the input egun data
    s.reset(r,vr,vz,ibeam,slice,seed1,seed2)
    # --- Install the routines in the proper places
    if not slice:
      installbeforefs(s.getparticles)
      installafterstep(s.makeparticleslive)
      # --- Force special steps
      top.allspecl = true
    # --- Make sure that normal injection is turned off
    top.inject = 0

  def reset(s,r,vr,vz,ibeam,slice=false,seed1=1341392,seed2=9823748):
    # --- Input quantities
    s.r = r
    s.vr = vr
    s.vz = vz
    s.ibeam = ibeam
    s.slice = slice
    # --- Initialize random number generator
    RandomArray.seed(seed1,seed2)
    # --- Derived quantities
    s.nn = len(s.r)
    s.ninject = s.ibeam/s.vz/top.sq[0]*s.vz*top.dt/top.sw[0]
    s.totalninject = sum(s.ninject)
    # --- If only loading a slice, do it now.
    # --- Also redo field solve and save new moments in the history arrays
    if slice:
      s.getparticles()
      s.makeparticleslive()
      getzmom.zmmnt()
      top.jhist = -1
      savehist(0.)

  # --- Routines which generate the particle coordinates and velocities.
  def gettime(s,i,ip):
    # --- Set starting time of the particles.
    # --- Times are chosen uniformly in the range 0 to top.dt.
    # --- When using slice only, times are all set to zero.
    if s.slice:
      return zeros(ip,'d')
    else:
      return ranf(zeros(ip,'d'))*top.dt
  def getr(s,i,ip,tt):
    # --- Get initial radius, including effect of radial velocity
    return s.r[i]*ones(ip,'d') + s.vr[i]*tt
  def getvr(s,i,ip):
    return s.vr[i]*ones(ip,'d')
  def gettheta(s,i,ip):
    return 2.*pi*ranf(zeros(ip,'d'))
  def getvtheta(s,i,ip):
    return zeros(ip,'d')
  def getz(s,i,ip,tt):
    return w3d.zmmin + s.vz[i]*tt
  def getvz(s,i,ip):
    return s.vz[i]*ones(ip,'d')

  # --- Generate the particles
  def getparticles(s):
    # --- Make sure there is enough space
    chckpart(top.pgroup,1,int(s.totalninject+1),0,true)
    print top.ins,top.npmax,int(s.totalninject+1)
    # --- Loop over the egun trajectories, load a ring of particles for each.
    injectedsofar = 0
    tempins = top.ins[0] - 1
    for i in range(s.nn):
      # --- Get number of particles to inject in this ring. Make sure that
      # --- on the last ring, enough particles are injected to get exactly
      # --- a total of totalninject of them.
      if i < s.nn-1:
        ip = nint(s.ninject[i])
        injectedsofar = injectedsofar + ip
      else:
        ip = nint(s.totalninject - injectedsofar)
      # --- Get the particle position and velocities
      tt = s.gettime(i,ip)
      rr = s.getr(i,ip,tt)
      vr = s.getvr(i,ip)
      theta = s.gettheta(i,ip)
      vtheta = s.getvtheta(i,ip)
      zz = s.getz(i,ip,tt)
      vz = s.getvz(i,ip)
      xx = rr*cos(theta)
      yy = rr*sin(theta)
      vx = vr*cos(theta) + vtheta*sin(theta)
      vy = vr*sin(theta) - vtheta*cos(theta)
      # --- Load the data into the WARP arrays
      print tempins,ip
      print shape(top.xp[tempins-ip:tempins]),shape(xx)
      print xx,zz
      top.xp[tempins-ip:tempins] = xx
      top.yp[tempins-ip:tempins] = yy
      top.zp[tempins-ip:tempins] = zz
      top.uxp[tempins-ip:tempins] = vx
      top.uyp[tempins-ip:tempins] = vy
      top.uzp[tempins-ip:tempins] = vz
      top.gaminv[tempins-ip:tempins] = 1.
      tempins = tempins - ip
    # --- Now add in the new particles contribution to the charge density
    loadrho(tempins+1,s.totalninject,1,0)

  def makeparticleslive(s):
    # --- Adjust WARP's particle indices to include new particles
    top.ins[0] = top.ins[0] - s.totalninject
    top.nps[0] = top.nps[0] + s.totalninject
    print "makeparticleslive",top.ins,top.nps

