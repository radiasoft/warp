from warp import *
from appendablearray import *
singleparticle_version = "$Id: singleparticle.py,v 1.3 2001/07/10 22:40:04 dave Exp $"

# --- Special code is needed here to make sure that top.ins and top.nps
# --- are set properly the first time an instance is created
_first_instance = 1

class SingleParticle:
  """
Class for running WARP in a single particle mode - no self-fields are
done and no diagnostic moments are calculated. Creator arguments...
 - x,y,z,vx,vy,vz: 6 coordinates. They can either be
     scalars (for one particle) or sequences (for multiple particles). The
     sequences must all be of the same length. If some are sequences and some
     scalars, the scalar values will be used for all particles. If vz
     is not specified, all the particles will use top.vbeam.
 - maxsteps=1000: an estimate of the number of steps taken. OK if value too
                  small
 - savedata=1: flag to turn on saving of particle trajectories
 - zerophi=1: when true, w3d.phi is zero out
 - resettime=0: when true, time and beam frame location reset to initial values

Available methods...
 - gett(i=0):  returns history of time for i'th particle
 - getx(i=0):  returns history of x for i'th particle
 - gety(i=0):  returns histort of y for i'th particle
 - getz(i=0):  returns histort of z for i'th particle
 - getvx(i=0): returns histort of vx for i'th particle
 - getvy(i=0): returns histort of vy for i'th particle
 - getvz(i=0): returns histort of vz for i'th particle
 - getgi(i=0): returns histort of gamma inverse for i'th particle

 - pxt(i=0), pyt(i=0), pzt(i=0), pvxt(i=0), pvyt(i=0), pvzt(i=0), pgit(i=0)
   pxy(i=0), pzx(i=0), pzy(i=0), pzvx(i=0), pzvy(i=0), pzvz(i=0)
   plots pairs of quantities for the i'th particle. All of these also
   take the same optional arguments as the plg command.

 - disable(): Clears out the particles. The history data is maintained, but
              the particles are no longer advanced.
 - enable(): Re-enables particles which have been disabled. They will
             continue from the place they were disabled.
  """

  #----------------------------------------------------------------------
  def __init__(self,x=0.,y=0.,z=0.,vx=0.,vy=0.,vz=None,
                    maxsteps=1000,savedata=1,zerophi=1,resettime=0):
    global _first_instance
    # --- Do some global initialization
    if _first_instance:
      top.ins[0] = 1
      top.nps[0] = 0
      top.sq[0] = top.zion*top.echarge
      top.sm[0] = top.aion*top.amu
      top.sw[0] = 0.
      _first_instance = 0
    self.savedata = savedata
    self.enabled = 0
    # --- Do some initialization
    self.spsetup(zerophi)
    # --- Setup particles
    self.spinit(x,y,z,vx,vy,vz)
    # --- Setup the lattice
    resetlat()
    setlatt()
    # --- Reset time if requested
    if resettime: self.resettime(resettime)
    # --- Setup history arrays
    self.setuphistory(maxsteps)

  #----------------------------------------------------------------------
  def spsetup(self,zerophi=1):
    """Sets up for a single particle run: turns diagnostics off, saves some
initial data, and redefines the step command
  - zerophi=1: when true, zeros out the w3d.phi array"""
    # --- Turn off all of the diagnostics
    wxy.ldiag = 0
    top.ifzmmnt = 0
    top.itmomnts = 0
    top.itplps = 0
    top.itplfreq = 0
    top.zzmomnts = 0
    top.zzplps = 0
    top.zzplfreq = 0
    top.nhist = top.nt
    top.iflabwn = 0
    w3d.lrhodia3d = false
    # --- Turn off the field solver
    top.fstype = -1
    # --- Zero out phi if requested.
    if zerophi: w3d.phi = 0.
    # --- Turn off charge deposition
    top.depos = 'none'
    # --- Save initial time-step data
    self.itsave = top.it
    self.timesave = top.time
    self.zbeamsave = top.zbeam
    self.zgridsave = top.zgrid
    self.zgridprvsave = top.zgridprv
    self.vbeamsave = top.vbeam
    self.vbeamfrmsave = top.vbeamfrm
    self.dtsave = top.dt
    self.dssave = wxy.ds
    # --- Add routine after step to save data
    installafterstep(self.spsavedata)

  #----------------------------------------------------------------------
  # --- Initialize the single particle.
  def spinit(self,x=0.,y=0.,z=0.,vx=0.,vy=0.,vz=None,
             maxsteps=1000):
    "Initializes one of more particles to run independently"
    # --- Set default value of vz and make sure it is not zero
    if not vz: vz = top.vbeam
    if vz == 0.: vz = top.smallpos
    # --- Make sure that the coordinates are not of list or tuple type
    if type(x)  in [ListType,TupleType]: x = array(x)
    if type(y)  in [ListType,TupleType]: y = array(y)
    if type(z)  in [ListType,TupleType]: z = array(z)
    if type(vx) in [ListType,TupleType]: vx = array(vx)
    if type(vy) in [ListType,TupleType]: vy = array(vy)
    if type(vz) in [ListType,TupleType]: vz = array(vz)
    # --- Find number of particles
    self.nn = 1
    if type(x)  is ArrayType: self.nn = len(x)
    if type(y)  is ArrayType: self.nn = len(y)
    if type(z)  is ArrayType: self.nn = len(z)
    if type(vx) is ArrayType: self.nn = len(vx)
    if type(vy) is ArrayType: self.nn = len(vy)
    if type(vz) is ArrayType: self.nn = len(vz)
    # --- Store the starting values
    self.x = x
    self.y = y
    self.z = z
    self.vx = vx
    self.vy = vy
    self.vz = vz
    # --- Set gamma inverse
    if top.lrelativ:
      self.gi = sqrt(1.- (self.vx**2+self.vy**2+self.vz**2)/clight**2)
      gamma = 1./self.gi
      self.vx = self.vx*gamma
      self.vy = self.vy*gamma
      self.vz = self.vz*gamma
    else:
      self.gi = self.nn*[1.]
    # --- Initialize particle coordinates
    self.enable()

  #----------------------------------------------------------------------
  def enable(self):
    """Load data into fortran arrays"""
    if self.enabled: return
    self.enabled = 1
    # --- make sure there is space
    chckpart(0,0,self.nn,false)
    # --- load the data
    ip1 = top.ins[0] - 1 + top.nps[0]
    top.nps[0] = top.nps[0] + self.nn
    ip2 = top.ins[0] + top.nps[0] - 1
    self.ip1 = ip1
    self.ip2 = ip2
    top.xp[ip1:ip2] = self.x
    top.yp[ip1:ip2] = self.y
    top.zp[ip1:ip2] = self.z
    top.uxp[ip1:ip2] = self.vx
    top.uyp[ip1:ip2] = self.vy
    top.uzp[ip1:ip2] = self.vz
    top.gaminv[ip1:ip2] = self.gi

  #----------------------------------------------------------------------
  def disable(self):
    """Disables the particles"""
    if not self.enabled: return
    self.enabled = 0
    # --- Save the last value in case the particles are re-enabled
    self.x = top.xp[self.ip1:self.ip2]
    self.y = top.yp[self.ip1:self.ip2]
    self.z = top.zp[self.ip1:self.ip2]
    self.vx = top.uxp[self.ip1:self.ip2]
    self.vy = top.uyp[self.ip1:self.ip2]
    self.vz = top.uzp[self.ip1:self.ip2]
    self.gi = top.gaminv[self.ip1:self.ip2]
    # --- Set uzp to zero, single of dead particles
    top.uzp[ip1:ip2] = 0.
    # --- Remove particles from looping if they are at either end of the
    # --- particle data in the arrays
    if self.ip1 == top.ins[0]:
      top.ins[0] = top.ins[0] + self.nn
      top.nps[0] = top.nps[0] - self.nn
    if self.ip2 == top.ins[0] + top.nps[0] - 1:
      top.nps[0] = top.nps[0] - self.nn
    # --- remove routine from after step
    uninstallafterstep(self.spsavedata)

  #----------------------------------------------------------------------
  def resettime(self):
    # --- Initialize time stepping
    if 0:
      top.it = self.itsave
      top.time = self.timesave
      top.zbeam = self.zbeamsave
      top.zgrid = self.zgridsave
      top.zgridprv = self.zgridprvsave
      top.vbeam = self.vbeamsave
      top.vbeamfrm = self.vbeamfrmsave
      top.dt = self.dtsave
      wxy.ds = self.dssave

  #----------------------------------------------------------------------
  def setuphistory(self,maxsteps=1000):
    # --- Create arrays for saving trajectory
    if self.savedata:
      self.spt = []
      self.spx = []
      self.spy = []
      self.spz = []
      self.spvx = []
      self.spvy = []
      self.spvz = []
      self.spgi = []
      for i in xrange(self.nn):
        self.spt.append(AppendableArray(maxsteps,'d'))
        self.spx.append(AppendableArray(maxsteps,'d'))
        self.spy.append(AppendableArray(maxsteps,'d'))
        self.spz.append(AppendableArray(maxsteps,'d'))
        self.spvx.append(AppendableArray(maxsteps,'d'))
        self.spvy.append(AppendableArray(maxsteps,'d'))
        self.spvz.append(AppendableArray(maxsteps,'d'))
        self.spgi.append(AppendableArray(maxsteps,'d'))
      self.spsavedata()

  #----------------------------------------------------------------------
  def spsavedata(self):
    """Saves data"""
    if not self.savedata: return
    for i in xrange(self.nn):
      self.spt[i].append(top.time)
      self.spx[i].append(top.xp[self.ip1 + i])
      self.spy[i].append(top.yp[self.ip1 + i])
      self.spz[i].append(top.zp[self.ip1 + i])
      self.spvx[i].append(top.uxp[self.ip1 + i])
      self.spvy[i].append(top.uyp[self.ip1 + i])
      self.spvz[i].append(top.uzp[self.ip1 + i])
      self.spgi[i].append(top.gaminv[self.ip1 + i])

  #----------------------------------------------------------------------
  def getsavedata(self):
    "Retrieves the save particle trajectories on a single array"
    if self.nn == 1:
      return array((self.spx[0].data(),self.spy[0].data(),self.spz[0].data(),
                    self.spvx[0].data(),self.spvy[0].data(),self.spvz[0].data(),
                    self.spgi[0].data()))
    else:
      r = zeros((self.nn,7,self.spx[0].len()),'d')
      for i in xrange (self.nn):
        r[i,0,:] = self.spx[i].data()
        r[i,1,:] = self.spy[i].data()
        r[i,2,:] = self.spz[i].data()
        r[i,3,:] = self.spvx[i].data()
        r[i,4,:] = self.spvy[i].data()
        r[i,5,:] = self.spvz[i].data()
        r[i,6,:] = self.spgi[i].data()
      return r

  #----------------------------------------------------------------------
  def gett(self,i=0):  return self.spt[i].data()
  def getx(self,i=0):  return self.spx[i].data()
  def gety(self,i=0):  return self.spy[i].data()
  def getz(self,i=0):  return self.spz[i].data()
  def getvx(self,i=0): return self.spvx[i].data()*self.spgi[i].data()
  def getvy(self,i=0): return self.spvy[i].data()*self.spgi[i].data()
  def getvz(self,i=0): return self.spvz[i].data()*self.spgi[i].data()
  def getgi(self,i=0): return self.spgi[i].data()

  #----------------------------------------------------------------------
  def pxt(self,i=0,**kw): apply(plg,(self.getx(i),self.gett(i)),kw)
  def pyt(self,i=0,**kw): apply(plg,(self.gety(i),self.gett(i)),kw)
  def pzt(self,i=0,**kw): apply(plg,(self.getz(i),self.gett(i)),kw)
  def pvxt(self,i=0,**kw): apply(plg,(self.getvx(i),self.gett(i)),kw)
  def pvyt(self,i=0,**kw): apply(plg,(self.getvy(i),self.gett(i)),kw)
  def pvzt(self,i=0,**kw): apply(plg,(self.getvz(i),self.gett(i)),kw)
  def pgit(self,i=0,**kw): apply(plg,(self.getgi(i),self.gett(i)),kw)
  def pxy(self,i=0,**kw): apply(plg,(self.gety(i),self.getx(i)),kw)
  def pzx(self,i=0,**kw): apply(plg,(self.getx(i),self.getz(i)),kw)
  def pzy(self,i=0,**kw): apply(plg,(self.gety(i),self.getz(i)),kw)
  def pzvx(self,i=0,**kw): apply(plg,(self.getvx(i),self.getz(i)),kw)
  def pzvy(self,i=0,**kw): apply(plg,(self.getvy(i),self.getz(i)),kw)
  def pzvz(self,i=0,**kw): apply(plg,(self.getvz(i),self.getz(i)),kw)

