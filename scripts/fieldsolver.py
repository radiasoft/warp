"""Generic class describing the interface needed for a field solver.
"""
from warp import *
import __main__

#=============================================================================
def loadrho(pgroup=None,ins_i=-1,nps_i=-1,is_i=-1,lzero=true):
  """
loadrho(pgroup=None,ins_i=-1,nps_i=-1,is_i=-1,lzero=1)
This routine provides a simple call from the interpreter to load the
rho array.  All of the arguments are optional.
If the species is not specified, all species are loaded, except
when ins or nps are specified, then only species 1 is loaded.
lzero is used to set whether or not rho is zeroed before the load.
The default is to zero out rho.
  """

  # --- Use top.pgroup as the default
  if pgroup is None: pgroup = top.pgroup

  # --- if particle location is specified but species is not, set so
  # --- only species number 1 is included
  if (ins_i != -1 and is_i == -1): is_i = 1

  # --- set number of particles
  if (ins_i != -1 and nps_i == -1):
    # --- if particle number is omitted but particle location is specified,
    # --- set nps to get rest of active particles of species
    nps_i = pgroup.nps[is_i] + pgroup.ins[is_i] - ins_i

  # --- if particle number is specified but species is not, set so
  # --- only species number 1 is included
  if (nps_i != -1 and is_i == -1): is_i = 1

  # --- Now call the appropriate compiled interface routine based on the
  # --- current package
  currpkg = package()[0]
  if (currpkg == "w3d"):
    loadrho3d(pgroup,ins_i,nps_i,is_i,lzero)
  elif (currpkg == "wxy"):
    loadrhoxy(pgroup,ins_i,nps_i,is_i,lzero)

#=============================================================================
def fieldsol(iwhich=0,lbeforefs=false,lafterfs=false):
  """
This routine provides a simple call from the interpreter to do the fieldsol.
It calls the appropriate compiled interface routine based on the current
package. Only w3d and wxy have field solves defined.
 - iwhich=0: specifies what action to take
 - lbeforefs=false: when true, call functions installed be installbeforefs
 - lafterfs=false:  when true, call functions installed be installafterfs
  """
  if lbeforefs: controllers.beforefs()
  currpkg = package()[0]
  if   (currpkg == "w3d"): fieldsol3d(iwhich)
  elif (currpkg == "wxy"): fieldsolxy(iwhich)
  if lafterfs: controllers.afterfs()
  # --- Now do extra work, updating arrays which depend directly on phi,
  # --- but only when a complete field solve was done.
  if iwhich == -1 or iwhich == 0:
    if (top.efetch == 3 and
        (w3d.solvergeom == w3d.XYZgeom or
         w3d.solvergeom == w3d.RZgeom or
         w3d.solvergeom == w3d.XZgeom or
         w3d.solvergeom == w3d.Rgeom  or
         w3d.solvergeom == w3d.Zgeom)):
      getselfe3d(w3d.phi,w3d.nx,w3d.ny,w3d.nz,w3d.selfe,
                 w3d.nx_selfe,w3d.ny_selfe,w3d.nz_selfe,
                 w3d.dx,w3d.dy,w3d.dz,
                 top.pboundxy,top.pboundxy,top.pboundxy,top.pboundxy)
    if top.inject > 0:
      try:
        # --- This routine is not defined in pywarp77
        getinj_phi()
      except NameError:
        pass

#=============================================================================
def loadj(pgroup=None,ins_i=-1,nps_i=-1,is_i=-1,lzero=true):
  """
loadj(ins_i=-1,nps_i=-1,is_i=-1,lzero=1)
This routine provides a simple call from the interpreter to load the
current density.  All of the arguments are optional.
If the species is not specified, all species are loaded, except
when ins or nps are specified, then only species 1 is loaded.
lzero is used to set whether or not rho is zeroed before the load.
The default is to zero out rho.
  """

  # --- if particle location is specified but species is not, set so
  # --- only species number 1 is included
  if (ins_i != -1 and is_i == -1): is_i = 1

  # --- set number of particles
  if (ins_i != -1 and nps_i == -1):
    # --- if particle number is omitted but particle location is specified,
    # --- set nps to get rest of active particles of species
    nps_i = pgroup.nps[is_i] + pgroup.ins[is_i] - ins_i

  # --- if particle number is specified but species is not, set so
  # --- only species number 1 is included
  if (nps_i != -1 and is_i == -1): is_i = 1

  # --- If pgroup is not given, then use the default one in top.
  if pgroup is None: pgroup = top.pgroup

  # --- Now call the appropriate compiled interface routine based on the
  # --- current package
  currpkg = package()[0]
  if (currpkg == "w3d"):
    loadj3d(pgroup,ins_i,nps_i,is_i,lzero)
  elif (currpkg == "wxy"):
    #loadrhoxy(ins_i,nps_i,is_i,lzero)
    print "loadj  not support in wxy yet"

#=============================================================================
# --- Setup routines which give access to fortran any field solver
_fieldsolvers = []
def registersolver(solver):
  """
Registers solvers to be used in the particle simulation.
 - solver: is the solver object. It must have the methods loadrho, solve, and
           fetche defined. fetchb and fetchphi will also be needed in some
           cases.

  """
  _fieldsolvers.append(solver)
  top.fstype = 12
def getregisteredsolver(i=0):
  if len(_fieldsolvers) == 0: return None
  return _fieldsolvers[i]
def unregistersolver(solver):
  assert solver in _fieldsolvers,"Specified solver was not registered"
  _fieldsolvers.remove(solver)
def loadrhoregistered():
  assert len(_fieldsolvers) > 0,"No solver has been registered"
  for f in _fieldsolvers:
    f.loadrho(lzero=not top.laccumulate_rho)
def loadjregistered():
  assert len(_fieldsolvers) > 0,"No solver has been registered"
  for f in _fieldsolvers:
    f.loadj(lzero=not top.laccumulate_rho)
def fieldsolregistered():
  assert len(_fieldsolvers) > 0,"No solver has been registered"
  for f in _fieldsolvers:
    f.solve()
def fetcheregistered():
  assert len(_fieldsolvers) > 0,"No solver has been registered"
  for f in _fieldsolvers:
    f.fetche()
def fetchbregistered():
  assert len(_fieldsolvers) > 0,"No solver has been registered"
  for f in _fieldsolvers:
    f.fetchb()
def fetchphiregistered():
  assert len(_fieldsolvers) > 0,"No solver has been registered"
  for f in _fieldsolvers:
    f.fetchphi()
def fetcharegistered():
  assert len(_fieldsolvers) > 0,"No solver has been registered"
  for f in _fieldsolvers:
    f.fetcha()
def rhodiaregistered():
  assert len(_fieldsolvers) > 0,"No solver has been registered"
  for f in _fieldsolvers:
    f.rhodia()
def gtlchgregistered():
  assert len(_fieldsolvers) > 0,"No solver has been registered"
  for f in _fieldsolvers:
    f.gtlchg()
def srhoaxregistered():
  assert len(_fieldsolvers) > 0,"No solver has been registered"
  for f in _fieldsolvers:
    f.srhoax()
def geteseregistered():
  assert len(_fieldsolvers) > 0,"No solver has been registered"
  for f in _fieldsolvers:
    f.getese()
def sphiaxregistered():
  assert len(_fieldsolvers) > 0,"No solver has been registered"
  for f in _fieldsolvers:
    f.sphiax()
def sezaxregistered():
  assert len(_fieldsolvers) > 0,"No solver has been registered"
  for f in _fieldsolvers:
    f.sezax()

def initfieldsolver():
    if w3d.AMRlevels>0:
      if 'AMRtree' not in __main__.__dict__:
        import AMR
        AMRtree=AMR.AMRTree()
        __main__.__dict__['AMRtree'] = AMRtree
        gchange('AMR')
__main__.__dict__['loadrhoregistered'] = loadrhoregistered
__main__.__dict__['loadjregistered'] = loadjregistered
__main__.__dict__['fieldsolregistered'] = fieldsolregistered
__main__.__dict__['fetcheregistered'] = fetcheregistered
__main__.__dict__['fetchbregistered'] = fetchbregistered
__main__.__dict__['fetchphiregistered'] = fetchphiregistered
__main__.__dict__['fetcharegistered'] = fetcharegistered
__main__.__dict__['rhodiaregistered'] = rhodiaregistered
__main__.__dict__['gtlchgregistered'] = gtlchgregistered
__main__.__dict__['srhoaxregistered'] = srhoaxregistered
__main__.__dict__['geteseregistered'] = geteseregistered
__main__.__dict__['sphiaxregistered'] = sphiaxregistered
__main__.__dict__['sezaxregistered'] = sezaxregistered
__main__.__dict__['initfieldsolver'] = initfieldsolver

#=============================================================================
# --- Setup routines which give access to fortran any B field solver
_bfieldsolver = [None]
def registerbsolver(bsolver):
  """
Registers the B field solver to be used in the particle simulation.
 - bsolver: is the solver object. It must have the methods loadj, solve, and
            fetchb defined. fetcha and fetchj will also be needed in some
            cases.

  """
  _bfieldsolver[0] = bsolver
  top.bfstype = 12
def getregisteredbsolver():
  return _bfieldsolver[0]
def bloadjregistered():
  assert _bfieldsolver[0] is not None,"No B solver has been registered"
  _bfieldsolver[0].loadj(lzero=not top.laccumulate_rho)
def bfieldsolregistered():
  assert _bfieldsolver[0] is not None,"No B solver has been registered"
  _bfieldsolver[0].solve()
def bfetchbregistered():
  assert _bfieldsolver[0] is not None,"No B solver has been registered"
  _bfieldsolver[0].fetchb()
def bfetcharegistered():
  assert _bfieldsolver[0] is not None,"No B solver has been registered"
  _bfieldsolver[0].fetcha()
def initbfieldsolver():
   pass
__main__.__dict__['bloadjregistered'] = bloadjregistered
__main__.__dict__['bfieldsolregistered'] = bfieldsolregistered
__main__.__dict__['bfetchbregistered'] = bfetchbregistered
__main__.__dict__['bfetcharegistered'] = bfetcharegistered
__main__.__dict__['initbfieldsolver'] = initbfieldsolver

#=============================================================================
#=============================================================================
#=============================================================================
class FieldSolver(object):
  """
Base class for a field solver that can be registered with Warp.  The
first block of routines must be defined since they are called by Warp
during a time step. Note that the load and fetch routines are written
out and can be used, but the methods called by them (which are shown at
the bottom of the class) must be redefined.  They are primarily written
out to show how to access the particle data and where to put the results
for the fetch routines.

The installconductor routine only needs to be redefined if the
conductors are used. The diagnostic routines only need to be defined if
the diagnostic is of interest and is meaningfull.
  """

  __flaginputs__ = {'lreducedpickle':0,'lnorestoreonpickle':0}

  def __init__(self,**kw):

    for name,defvalue in FieldSolver.__flaginputs__.iteritems():
      if name not in self.__dict__:
        #self.__dict__[name] = kw.pop(name,getattr(top,name)) # Python2.3
        self.__dict__[name] = kw.get(name,defvalue)
      if kw.has_key(name): del kw[name]

  def __getstate__(self):
    dict = self.__dict__.copy()
    return dict

  def __setstate__(self,dict):
    self.__dict__.update(dict)

  # ---------------------------------------------------------------------
  # --- These routines must at least be defined.
  def loadrho(self,lzero=true,**kw):
    'Charge deposition, uses particles from top directly'
    if lzero: self.zerorhop()

    for js,i,n,q,w in zip(arange(top.pgroup.ns),top.pgroup.ins-1,
                          top.pgroup.nps,top.pgroup.sq,top.pgroup.sw):
      if n > 0:
        self.setrhop(top.pgroup.xp[i:i+n],top.pgroup.yp[i:i+n],
                     top.pgroup.zp[i:i+n],top.pgroup.uzp[i:i+n],
                     q,w*top.pgroup.dtscale[js])

  def loadj(self,lzero=true,**kw):
    'Charge deposition, uses particles from top directly'
    if lzero: self.zeroj()
    for js in range(top.pgroup.ns):
      i = top.pgroup.ins-1
      n = top.pgroup.nps
      q = top.pgroup.sq
      w = top.pgroup.sw
      x = top.pgroup.xp[i:i+n]
      y = top.pgroup.yp[i:i+n]
      z = top.pgroup.zp[i:i+n]
      ux = top.pgroup.uxp[i:i+n]
      uy = top.pgroup.uyp[i:i+n]
      uz = top.pgroup.uzp[i:i+n]
      gaminv = top.pgroup.gaminv[i:i+n]
      if top.wpid > 0: wght = top.pgroup[i:i+n,top.wpid-1]
      else:            wght = None
      self.setj(x,y,z,ux,uy,uz,gaminv,wght,q,w)
    self.makejperiodic()
    self.getjforfieldsolve()

  def fetche(self,**kw):
    'Fetches the E field, uses arrays from w3d module FieldSolveAPI'
    if w3d.npfsapi == 0: return
    ipmin = w3d.ipminfsapi
    x = w3d.pgroupfsapi.xp[ipmin-1:ipmin-1+w3d.npfsapi]
    y = w3d.pgroupfsapi.yp[ipmin-1:ipmin-1+w3d.npfsapi]
    z = w3d.pgroupfsapi.zp[ipmin-1:ipmin-1+w3d.npfsapi]
    ex = w3d.pgroupfsapi.ex[ipmin-1:ipmin-1+w3d.npfsapi]
    ey = w3d.pgroupfsapi.ey[ipmin-1:ipmin-1+w3d.npfsapi]
    ez = w3d.pgroupfsapi.ez[ipmin-1:ipmin-1+w3d.npfsapi]
    self.fetchefrompositions(x,y,z,ex,ey,ez,w3d.pgroupfsapi)

  def fetchb(self,**kw):
    'Fetches the B field, uses arrays from w3d module FieldSolveAPI'
    if w3d.npfsapi == 0: return
    ipmin = w3d.ipminfsapi
    x = w3d.pgroupfsapi.xp[ipmin-1:ipmin-1+w3d.npfsapi]
    y = w3d.pgroupfsapi.yp[ipmin-1:ipmin-1+w3d.npfsapi]
    z = w3d.pgroupfsapi.zp[ipmin-1:ipmin-1+w3d.npfsapi]
    bx = w3d.pgroupfsapi.bx[ipmin-1:ipmin-1+w3d.npfsapi]
    by = w3d.pgroupfsapi.by[ipmin-1:ipmin-1+w3d.npfsapi]
    bz = w3d.pgroupfsapi.bz[ipmin-1:ipmin-1+w3d.npfsapi]
    self.fetchbfrompositions(x,y,z,ex,ey,ez,w3d.pgroupfsapi)

  def fetchphi(self):
    'Fetches the potential, uses arrays from w3d module FieldSolveAPI'
    if w3d.npfsapi == 0: return
    x = w3d.xfsapi
    y = w3d.yfsapi
    z = w3d.zfsapi
    self.fetchphifrompositions(x,y,z,w3d.phifsapi)

  def fetcha(self):
    'Fetches the magnetostatic potential, uses arrays from w3d module FieldSolveAPI'
    if w3d.npfsapi == 0: return
    ipmin = w3d.ipminfsapi
    x = w3d.pgroupfsapi.xp[ipmin-1:ipmin-1+w3d.npfsapi]
    y = w3d.pgroupfsapi.yp[ipmin-1:ipmin-1+w3d.npfsapi]
    z = w3d.pgroupfsapi.zp[ipmin-1:ipmin-1+w3d.npfsapi]
    self.fetchafrompositions(x,y,z,w3d.afsapi)

  def solve(self,iwhich=0):
    'do the field solve'
    pass

  def installconductor(self,conductor):
    pass

  def find_mgparam(self,lsavephi=false,resetpasses=0):
    find_mgparam(lsavephi=lsavephi,resetpasses=resetpasses,
                 solver=self,pkg3d=self)

  # --- Diagnostic routines
  def rhodia(self):
    pass
  def gtlchg(self):
    pass
  def srhoax(self):
    pass
  def getese(self):
    pass
  def sphiax(self):
    pass
  def sezax(self):
    pass

  def pfzx(self,**kw):
    'Plots potential in z-x plane. Does nothing if not defined.'
    pass
  def pfzy(self,**kw):
    'Plots potential in z-y plane. Does nothing if not defined.'
    pass
  def pfxy(self,**kw):
    'Plots potential in x-y plane. Does nothing if not defined.'
    pass

  # ---------------------------------------------------------------------
  # --- These can optionally be defined, making use of some of the above
  # --- routines.
  def zerorhop(self):
    pass

  def setrhop(self,x,y,z,vz,q,w,**kw):
    pass

  def zeroj(self):
    pass

  def setj(self,x,y,z,ux,uy,uz,gaminv,wght,q,w):
    pass

  def makejperiodic(self):
    pass

  def getjforfieldsolve(self):
    pass

  def fetchefrompositions(self,x,y,z,ex,ey,ez):
    'Fetches E at the given positions'
    pass

  def fetchbfrompositions(self,x,y,z,bx,by,bz):
    'Fetches B at the given positions'
    pass

  def fetchphifrompositions(self,x,y,z,phi):
    'Fetches the potential, given a list of positions'
    pass

  def fetchafrompositions(self,x,y,z,a):
    'Fetches the magnetostatic potential, given a list of positions'
    pass

#=============================================================================
#=============================================================================
#=============================================================================
class SubcycledPoissonSolver(FieldSolver):

  def __getstate__(self):
    dict = FieldSolver.__getstate__(self)
    if self.lreducedpickle:
      if 'rhoparray' in dict: del dict['rhoparray']
      if 'phiparray' in dict: del dict['phiparray']
    return dict

  def __setstate__(self,dict):
    FieldSolver.__setstate__(self,dict)
    if self.lreducedpickle and not self.lnorestoreonpickle:
      self.allocatedataarrays()

  def setrhopforparticles(self,irhopndtscopies,indts,iselfb):
    if irhopndtscopies is not None:
      rhopndts = self.getrhop(lndts=1)
      self.rhop = rhopndts[...,irhopndtscopies,indts]
    else:
      rhopselfb = self.getrhop(lselfb=1)
      self.rhop = rhopselfb[...,iselfb]

  def aftersetrhop(self,lzero):
    "Anything that needs to be done to rhop after the deposition"
    pass

  def loadrho(self,lzero=None,**kw):
    'Charge deposition, uses particles from top directly'
    if lzero is None: lzero = w3d.lzerorhofsapi
    self.allocatedataarrays()
    if lzero: self.zerorhop()

    for js,i,n,q,w in zip(arange(top.pgroup.ns),top.pgroup.ins-1,
                          top.pgroup.nps,top.pgroup.sq,top.pgroup.sw):
      if n > 0:
        args = [top.pgroup.xp[i:i+n],top.pgroup.yp[i:i+n],
                top.pgroup.zp[i:i+n],top.pgroup.uzp[i:i+n],
                q,w*top.pgroup.dtscale[js],top.zgrid]
        if top.pgroup.ldts[js]:
          indts = top.ndtstorho[top.pgroup.ndts[js]-1]
          self.setrhopforparticles(0,indts,None)
          zgrid = top.zgrid + (top.zgrid-top.zgridprv)*(top.pgroup.ndts[js]-1)
          args[-1] = zgrid
          self.setrhop(*args)
        if top.pgroup.iselfb[js] > -1:
          iselfb = top.pgroup.iselfb[js]
          self.setrhopforparticles(None,None,iselfb)
          self.setrhop(*args)

    if lzero:
      for indts in range(top.nsndts):
        if top.ldts[indts]:
          self.setrhopforparticles(0,indts,None)
          self.aftersetrhop(lzero)
      for iselfb in range(top.nsselfb):
        self.setrhopforparticles(None,None,iselfb)
        self.aftersetrhop(lzero)

      self.averagerhopwithsubcycling()

  def loadj(self,lzero=true,**kw):
    pass

  def getphip(self,lndts=None,lselfb=None):
    if lndts is None and lselfb is None:
      return self.phiparray
    elif lndts is not None:
      return self.phiparray[...,:top.nsndtsphi]
    elif lselfb is not None:
      return self.phiparray[...,top.nsndtsphi:]

  def getrhop(self,lndts=None,lselfb=None):
    if lndts is None and lselfb is None:
      return self.rhoparray
    elif lndts is not None:
      rhop = self.rhoparray[...,:top.nrhopndtscopies*top.nsndts]
      trhop = transpose(rhop)
      sss = list(shape(rhop)[:-1])+[top.nrhopndtscopies,top.nsndts]
      sss.reverse()
      trhop.shape = sss
      return transpose(trhop)
    elif lselfb is not None:
      return self.rhoparray[...,top.nrhopndtscopies*top.nsndts:]

  def setphipforparticles(self,indts,iselfb):
    "Sets self.phip to the currently active phi"
    if indts is not None:
      phipndts = self.getphip(lndts=1)
      self.phip = phipndts[...,min(indts,w3d.nsndtsphi3d-1)]
    else:
      phipselfb = self.getphip(lselfb=1)
      self.phip = phipselfb[...,iselfb]

  def fetche(self,**kw):
    'Fetches the E field, uses arrays from w3d module FieldSolveAPI'
    if w3d.npfsapi == 0: return
    ipmin = w3d.ipminfsapi
    x = w3d.pgroupfsapi.xp[ipmin-1:ipmin-1+w3d.npfsapi]
    y = w3d.pgroupfsapi.yp[ipmin-1:ipmin-1+w3d.npfsapi]
    z = w3d.pgroupfsapi.zp[ipmin-1:ipmin-1+w3d.npfsapi]
    ex = w3d.pgroupfsapi.ex[ipmin-1:ipmin-1+w3d.npfsapi]
    ey = w3d.pgroupfsapi.ey[ipmin-1:ipmin-1+w3d.npfsapi]
    ez = w3d.pgroupfsapi.ez[ipmin-1:ipmin-1+w3d.npfsapi]
    args = [x,y,z,ex,ey,ez,w3d.pgroupfsapi]

    js = w3d.jsfsapi
    if js < 0: ndtstorho = 0
    else:      ndtstorho = top.ndtstorho[top.pgroup.ndts[js]-1]

    tmpnsndts = getnsndtsforsubcycling()
    indts = min(tmpnsndts-1,ndtstorho)
    self.setphipforparticles(indts,None)
    self.fetchefrompositions(*args,**kw)

    # add self B correction as needed
    if js >= 0 and top.pgroup.iselfb[js] > -1:
      extemp = zeros(shape(x),'d')
      eytemp = zeros(shape(y),'d')
      eztemp = zeros(shape(z),'d')
      self.setphipforparticles(None,top.pgroup.iselfb[js])
      args = [x,y,z,extemp,eytemp,eztemp,w3d.pgroupfsapi]
      self.fetchefrompositions(*args,**kw)
      ex[...] += extemp
      ey[...] += eytemp

  def fetchphi(self):
    'Fetches the potential, uses arrays from w3d module FieldSolveAPI'
   #if w3d.npfsapi == 0: return
   #ipmin = w3d.ipminfsapi
   #x = w3d.pgroupfsapi.xp[ipmin-1:ipmin-1+w3d.npfsapi]
   #y = w3d.pgroupfsapi.yp[ipmin-1:ipmin-1+w3d.npfsapi]
   #z = w3d.pgroupfsapi.zp[ipmin-1:ipmin-1+w3d.npfsapi]
    if w3d.npfsapi == 0: return
    x = w3d.xfsapi
    y = w3d.yfsapi
    z = w3d.zfsapi

    js = w3d.jsfsapi
    if js < 0: ndtstorho = 0
    else:      ndtstorho = top.ndtstorho[top.pgroup.ndts[js]-1]
    tmpnsndts = getnsndtsforsubcycling()
    indts = min(tmpnsndts-1,ndtstorho)
    self.setphipforparticles(indts,None)
    self.fetchphifrompositions(x,y,z,w3d.phifsapi)

  def setrhoandphiforfieldsolve(self,irhopndtscopies,indts,iselfb):
    if irhopndtscopies is not None:
      rhopndts = self.getrhop(lndts=1)
      phipndts = self.getphip(lndts=1)
      self.rho = rhopndts[...,irhopndtscopies,indts]
      self.phi = phipndts[...,min(indts,w3d.nsndtsphi3d-1)]
    else:
      rhopselfb = self.getrhop(lselfb=1)
      phipselfb = self.getphip(lselfb=1)
      self.rho = rhopselfb[...,iselfb]
      self.phi = phipselfb[...,iselfb]

  def getphipforparticles(self,indts,iselfb):
    "Copies from phi to phip"
    if indts is not None:
      phipndts = self.getphip(lndts=1)
      phipndts[...,min(indts,w3d.nsndtsphi3d-1)] = self.phi[...]
    else:
      phipselfb = self.getphip(lselfb=1)
      phipselfb[...,iselfb] = self.phi[...]

  def selfbcorrection(self,js):
    "scale phispecies by -(1-1/gamma*2) store into top.fselfb"
    phipselfb = self.getphip(lselfb=1)
    phipselfb[...,js] *= top.fselfb[js]

  def dosolveonphi(self,iwhich,irhopndtscopies,indts,iselfb):
    "points rho and phi appropriately and call the solving routine"
    self.setrhoandphiforfieldsolve(irhopndtscopies,indts,iselfb)
    self.dosolve(iwhich,irhopndtscopies,indts,iselfb)
    if iselfb is not None: self.selfbcorrection(iselfb)
    self.getphipforparticles(indts,iselfb)

  def solve(self,iwhich=0):
    self.allocatedataarrays()
    # --- Loop over the subcyling groups and do any field solves that
    # --- are necessary.
    # --- Do loop in reverse order so that rho and phi end up with the arrays
    # --- for the speices with the smallest timestep.
    tmpnsndts = getnsndtsforsubcycling()
    for indts in range(tmpnsndts-1,-1,-1):
      if (not top.ldts[indts] and
          (top.ndtsaveraging == 0 and not sum(top.ldts))): continue
      self.dosolveonphi(iwhich,top.nrhopndtscopies-1,indts,None)

    # --- Solve for phi for groups which require the self B correction
    for js in range(top.nsselfb):
      self.dosolveonphi(iwhich,None,None,js)

  def getpdims(self):
    raise """getpdims must be supplied - it should return a list of the dimensions
of the rho and field arrays"""

  def allocatedataarrays(self):
    # --- Get base dimension of the arrays
    rhopdims,phipdims = self.getpdims()

    # --- Setup rho and phi arrays, including extra copies for subcycling
    # --- and self B corrections.
    setupSubcycling(top.pgroup)
    setupSelfB(top.pgroup)

    extrarhopdim = top.nrhopndtscopies*top.nsndts + top.nsselfb
    dims = list(rhopdims) + [extrarhopdim]
    if 'rhoparray' not in self.__dict__ or shape(self.rhoparray) != tuple(dims):
      self.rhoparray = fzeros(dims,'d')

    extraphidim = top.nsndtsphi + top.nsselfb
    dims = list(phipdims) + [extraphidim]
    if 'phiparray' not in self.__dict__ or shape(self.phiparray) != tuple(dims):
      self.phiparray = fzeros(dims,'d')

  def zerorhop(self):
    if top.ndtsaveraging == 0:
      self.zerorhopwithsampledsubcycling()
    elif top.ndtsaveraging == 1:
      self.zerorhopwithfullvsubcycling()
    elif top.ndtsaveraging == 2:
      self.zerorhopwithhalfvsubcycling()
    self.zerorhopselfb()

  def zerorhopwithsampledsubcycling(self):
    # --- Zero the rhop copy for species when the positions
    # --- are advanced.
    # --- rhopndts(...,2,in1) holds the old rhop which is still needed
    # --- for the faster advanced groups.
    rhopndts = self.getrhop(lndts=1)
    for in1 in range(top.nsndts):
      if top.ldts[in1]:
        if top.nrhopndtscopies == 2:
          rhopndts[...,1,in1] = rhopndts[...,0,in1]
        rhopndts[...,0,in1] = 0.

  def zerorhopwithfullvsubcycling(self):
    raise "fullv subcycling not yet implemented"

  def zerorhopwithhalfvsubcycling(self):
    raise "halfv subcycling not yet implemented"

  def averagerhopwithsubcycling(self):
    if top.ndtsaveraging == 0:
      self.averagerhopwithsampledsubcycling()
    elif top.ndtsaveraging == 1:
      self.averagerhopwithfullvsubcycling()
    elif top.ndtsaveraging == 2:
      self.averagerhopwithhalfvsubcycling()

  def averagerhopwithsampledsubcycling(self):
    if top.ndtsmax == 1: return

    rhopndts = self.getrhop(lndts=1)

    # --- During the generate, do the copy of the new rhop to the old rhop
    # --- for group 0, which is normally done during the zerorhop call.
    rhopndts[...,1,0] = rhopndts[...,0,0]

    # --- Save the rhop where the fastest particle's rhop is. For now,
    # --- assume that this is in1=0
    for in1 in range(1,top.nsndts):
      if top.it == 0:
        # --- At top.it==0, before the first step, always add the new rhop.
        rhopndts[...,1,0] = (rhopndts[...,1,0] + rhopndts[...,0,in1])
      elif top.rhotondts[in1]%2 == 1:
        # --- Use the rhop that is closest in time to the current time.
        if ((top.it-1)%top.rhotondts[in1] > top.rhotondts[in1]/2.-1.):
          rhopndts[...,1,0] = (rhopndts[...,1,0] + rhopndts[...,0,in1])
        else:
          rhopndts[...,1,0] = (rhopndts[...,1,0] + rhopndts[...,1,in1])
      else:
        # --- When ndts is even, at the mid point of the step, take the
        # --- average of the old and the new
        # --- Otherwise, use the rhop that is closest in time to the current
        # --- time.
        if (top.it-1)%top.rhotondts[in1] == top.rhotondts[in1]/2.-1.:
          rhopndts[...,1,0] = (rhopndts[...,1,0] +
               0.5*(rhopndts[...,0,in1] + rhopndts[...,1,in1]))
        elif ((top.it-1)%top.rhotondts[in1] > top.rhotondts[in1]/2.-1.):
          rhopndts[...,1,0] = (rhopndts[...,1,0] + rhopndts[...,0,in1])
        else:
          rhopndts[...,1,0] = (rhopndts[...,1,0] + rhopndts[...,1,in1])

  def zerorhopselfb(self):
    rhopselfb = self.getrhop(lselfb=1)
    rhopselfb[...] = 0.

