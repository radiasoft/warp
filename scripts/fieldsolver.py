"""Generic class describing the interface needed for a field solver.
"""
from warp import *
import __main__
import gc

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
      getselfe3d(w3d.phip,w3d.nxp,w3d.nyp,w3d.nzp,w3d.selfe,
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
    print "loadj not support in wxy yet"

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

  __w3dinputs__ = ['nx','ny','nz','dx','dy','dz','nzfull','nzpguard',
                   'xmmin','xmmax','ymmin','ymmax','zmmin','zmmax',
                   'zmminglobal','zmmaxglobal',
                   'bound0','boundnz','boundxy','l2symtry','l4symtry',
                   'solvergeom']
  __topinputs__ = ['pbound0','pboundnz','pboundxy','efetch',
                   'my_index','nslaves','lfsautodecomp','zslave','lautodecomp']
  __flaginputs__ = {'forcesymmetries':1,'lzerorhointerior':0,
                    'lreducedpickle':1,'lnorestoreonpickle':0}

  def __init__(self,**kw):
    try:
      kw['kwdict'].update(kw)
      kw = kw['kwdict']
      del kw['kwdict']
    except KeyError:
      pass

    # --- Save input parameters
    for name in FieldSolver.__w3dinputs__:
      if name not in self.__dict__:
        #self.__dict__[name] = kw.pop(name,getattr(w3d,name)) # Python2.3
        self.__dict__[name] = kw.get(name,getattr(w3d,name))
      if kw.has_key(name): del kw[name]
    for name in FieldSolver.__topinputs__:
      if name not in self.__dict__:
        #self.__dict__[name] = kw.pop(name,getattr(top,name)) # Python2.3
        self.__dict__[name] = kw.get(name,getattr(top,name))
      if kw.has_key(name): del kw[name]
    for name,defvalue in FieldSolver.__flaginputs__.iteritems():
      if name not in self.__dict__:
        #self.__dict__[name] = kw.pop(name,getattr(top,name)) # Python2.3
        self.__dict__[name] = kw.get(name,defvalue)
      if kw.has_key(name): del kw[name]

    # --- Make sure the top.nparpgrp is a large number. If it becomes too
    # --- small, fetche becomes inefficient since it is called many times,
    # --- once per each group. The not insignificant function call overhead
    # --- of python begins to use up a sizable chunk of time.
    top.nparpgrp = 100000

    # --- bounds is special since it will sometimes be set from the
    # --- variables bound0, boundnz, boundxy, l2symtry, and l4symtry
    if 'bounds' not in self.__dict__:
      if 'bounds' in kw:
        self.bounds = kw['bounds']
        del kw['bounds']
      else:
        self.bounds = zeros(6)
        self.bounds[0] = self.boundxy
        self.bounds[1] = self.boundxy
        self.bounds[2] = self.boundxy
        self.bounds[3] = self.boundxy
        self.bounds[4] = self.bound0
        self.bounds[5] = self.boundnz
        if self.l2symtry:
          self.bounds[2] = neumann
          if self.boundxy == periodic: self.bounds[3] = neumann
          if self.forcesymmetries: self.ymmin = 0.
        elif self.l4symtry:
          self.bounds[0] = neumann
          self.bounds[2] = neumann
          if self.boundxy == periodic: self.bounds[1] = neumann
          if self.boundxy == periodic: self.bounds[3] = neumann
          if self.forcesymmetries: self.xmmin = 0.
          if self.forcesymmetries: self.ymmin = 0.

    # --- pbounds is special since it will sometimes be set from the
    # --- variables pbound0, pboundnz, pboundxy, l2symtry, and l4symtry
    if 'pbounds' not in self.__dict__:
      if 'pbounds' in kw:
        self.pbounds = kw['pbounds']
        del kw['pbounds']
      else:
        self.pbounds = zeros(6)
        self.pbounds[0] = self.pboundxy
        self.pbounds[1] = self.pboundxy
        self.pbounds[2] = self.pboundxy
        self.pbounds[3] = self.pboundxy
        self.pbounds[4] = self.pbound0
        self.pbounds[5] = self.pboundnz
        if self.l2symtry:
          self.pbounds[2] = reflect
          if self.pboundxy == periodic: self.pbounds[3] = reflect
        elif self.l4symtry:
          self.pbounds[0] = reflect
          self.pbounds[2] = reflect
          if self.pboundxy == periodic: self.pbounds[1] = reflect
          if self.pboundxy == periodic: self.pbounds[3] = reflect

    # --- Set parallel related parameters and calculate mesh sizes
    if self.nslaves <= 1:
      self.my_index = 0
      if self.zmminglobal == self.zmmaxglobal:
        self.nzfull = self.nz
        self.zmminglobal = self.zmmin
        self.zmmaxglobal = self.zmmax
      self.izfsslave = zeros(1)
      self.nzfsslave = zeros(1) + self.nz
      self.nxp = self.nx
      self.nyp = self.ny
      self.nzp = self.nz
      self.xmminp = self.xmmin
      self.xmmaxp = self.xmmax
      self.ymminp = self.ymmin
      self.ymmaxp = self.ymmax
      self.zmminp = self.zmmin
      self.zmmaxp = self.zmmax
    else:
      self.my_index = me
      if self.zmminglobal == self.zmmaxglobal:
        self.nzfull = self.nz
        self.zmminglobal = self.zmmin
        self.zmmaxglobal = self.zmmax
      self.izfsslave = zeros(self.nslaves)
      self.nzfsslave = zeros(self.nslaves)
      self.grid_overlap = array([2])
      top.grid_overlap = 2
      domaindecomposefields(self.nzfull,self.nslaves,self.lfsautodecomp,
                            self.izfsslave,self.nzfsslave,self.grid_overlap)

      self.nz = self.nzfsslave[self.my_index]
      if self.dz == 0.: self.dz = (self.zmmaxglobal - self.zmminglobal)/self.nzfull
      self.zmmin = self.zmminglobal + self.izfsslave[self.my_index]*self.dz
      self.zmmax = self.zmminglobal + (self.izfsslave[self.my_index] + self.nzfsslave[self.my_index])*self.dz

      self.izpslave = zeros(self.nslaves)
      self.nzpslave = zeros(self.nslaves)
      self.zpslmin = zeros(self.nslaves,'d')
      self.zpslmax = zeros(self.nslaves,'d')
      domaindecomposeparticles(self.nzfull,self.nslaves,self.izfsslave,self.nzfsslave,
                               self.grid_overlap,self.nzpguard,
                               self.zmminglobal,self.zmmaxglobal,self.dz,self.zslave[:self.nslaves],
                               self.lautodecomp,self.izpslave,self.nzpslave,
                               self.zpslmin,self.zpslmax)

      self.nxp = self.nx
      self.nyp = self.ny
      self.nzp = self.nzpslave[self.my_index]
      self.xmminp = self.xmmin
      self.xmmaxp = self.xmmax
      self.ymminp = self.ymmin
      self.ymmaxp = self.ymmax
      self.zmminp = self.zmminglobal + self.izpslave[self.my_index]*self.dz
      self.zmmaxp = self.zmminglobal + (self.izpslave[self.my_index] + self.nzpslave[self.my_index])*self.dz

    if self.dx == 0.: self.dx = (self.xmmax - self.xmmin)/self.nx
    if self.dy == 0.: self.dy = (self.ymmax - self.ymmin)/self.ny
    if self.dz == 0.: self.dz = (self.zmmaxglobal - self.zmminglobal)/self.nzfull
    self.xsymmetryplane = 0.
    self.ysymmetryplane = 0.
    self.xmesh = self.xmmin + arange(0,self.nx+1)*self.dx
    self.ymesh = self.ymmin + arange(0,self.ny+1)*self.dy
    self.zmesh = self.zmmin + arange(0,self.nz+1)*self.dz
    self.zmeshglobal = self.zmminglobal + arange(0,self.nzfull+1)*self.dz

    self.ix_axis = nint(-self.xmmin/self.dx)
    self.iy_axis = nint(-self.ymmin/self.dy)
    self.iz_axis = nint(-self.zmminglobal/self.dz)

  def __getstate__(self):
    dict = self.__dict__.copy()
    # --- Flag whether this is the registered solver so it knows whether
    # --- to reregister itself upon the restore. The instance
    # --- is not registered if it is not going to be restored.
    if self is getregisteredsolver() and not self.lnorestoreonpickle:
      dict['iamtheregisteredsolver'] = 1
    else:
      dict['iamtheregisteredsolver'] = 0
    return dict

  def __setstate__(self,dict):
    self.__dict__.update(dict)
    if self.iamtheregisteredsolver and not self.lnorestoreonpickle:
      del self.iamtheregisteredsolver
      registersolver(self)

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
    return
    #if w3d.npfsapi == 0: return
    #ipmin = w3d.ipminfsapi
    #x = w3d.pgroupfsapi.xp[ipmin-1:ipmin-1+w3d.npfsapi]
    #y = w3d.pgroupfsapi.yp[ipmin-1:ipmin-1+w3d.npfsapi]
    #z = w3d.pgroupfsapi.zp[ipmin-1:ipmin-1+w3d.npfsapi]
    #ex = w3d.pgroupfsapi.ex[ipmin-1:ipmin-1+w3d.npfsapi]
    #ey = w3d.pgroupfsapi.ey[ipmin-1:ipmin-1+w3d.npfsapi]
    #ez = w3d.pgroupfsapi.ez[ipmin-1:ipmin-1+w3d.npfsapi]
    #self.fetchefrompositions(x,y,z,ex,ey,ez,w3d.pgroupfsapi)

  def fetchb(self,**kw):
    'Fetches the B field, uses arrays from w3d module FieldSolveAPI'
    return
    #if w3d.npfsapi == 0: return
    #ipmin = w3d.ipminfsapi
    #x = w3d.pgroupfsapi.xp[ipmin-1:ipmin-1+w3d.npfsapi]
    #y = w3d.pgroupfsapi.yp[ipmin-1:ipmin-1+w3d.npfsapi]
    #z = w3d.pgroupfsapi.zp[ipmin-1:ipmin-1+w3d.npfsapi]
    #bx = w3d.pgroupfsapi.bx[ipmin-1:ipmin-1+w3d.npfsapi]
    #by = w3d.pgroupfsapi.by[ipmin-1:ipmin-1+w3d.npfsapi]
    #bz = w3d.pgroupfsapi.bz[ipmin-1:ipmin-1+w3d.npfsapi]
    #self.fetchbfrompositions(x,y,z,ex,ey,ez,w3d.pgroupfsapi)

  def fetchphi(self):
    'Fetches the potential, uses arrays from w3d module FieldSolveAPI'
    return
    #if w3d.npfsapi == 0: return
    #x = w3d.xfsapi
    #y = w3d.yfsapi
    #z = w3d.zfsapi
    #self.fetchphifrompositions(x,y,z,w3d.phifsapi)

  def fetcha(self):
    'Fetches the magnetostatic potential, uses arrays from w3d module FieldSolveAPI'
    return
    #if w3d.npfsapi == 0: return
    #ipmin = w3d.ipminfsapi
    #x = w3d.pgroupfsapi.xp[ipmin-1:ipmin-1+w3d.npfsapi]
    #y = w3d.pgroupfsapi.yp[ipmin-1:ipmin-1+w3d.npfsapi]
    #z = w3d.pgroupfsapi.zp[ipmin-1:ipmin-1+w3d.npfsapi]
    #self.fetchafrompositions(x,y,z,w3d.afsapi)

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

  def fetchefrompositions(self,x,y,z,ex,ey,ez,pgroup=None):
    'Fetches E at the given positions'
    pass

  def fetchbfrompositions(self,x,y,z,bx,by,bz,pgroup=None):
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
      if 'sourceparray' in dict: del dict['sourceparray']
      if 'potentialparray' in dict: del dict['potentialparray']
      if 'fieldparray' in dict: del dict['fieldparray']
      if 'sourcep' in dict: del dict['sourcep']
      if 'potentialp' in dict: del dict['potentialp']
      if 'fieldp' in dict: del dict['fieldp']
      if 'source' in dict: del dict['source']
      if 'potential' in dict: del dict['potential']
      if 'field' in dict: del dict['field']
    return dict

  def __setstate__(self,dict):
    FieldSolver.__setstate__(self,dict)
    if self.lreducedpickle and not self.lnorestoreonpickle:
      installafterrestart(self.allocatedataarrays)

  def returnfieldp(self,lndts=None,lselfb=None):
    try:
      fieldparray = self.fieldparray
    except AttributeError:
      return None
    if lndts is None and lselfb is None:
      return fieldparray
    elif lndts is not None:
      return fieldparray[...,:top.nsndtsphi]
    elif lselfb is not None:
      return fieldparray[...,top.nsndtsphi:]

  def returnpotentialp(self,lndts=None,lselfb=None):
    if lndts is None and lselfb is None:
      return self.potentialparray
    elif lndts is not None:
      return self.potentialparray[...,:top.nsndtsphi]
    elif lselfb is not None:
      return self.potentialparray[...,top.nsndtsphi:]

  def returnsourcep(self,lndts=None,lselfb=None):
    if lndts is None and lselfb is None:
      return self.sourceparray
    elif lndts is not None:
      sourcep = self.sourceparray[...,:top.nrhopndtscopies*top.nsndts]
      tsourcep = transpose(sourcep)
      sss = list(shape(sourcep)[:-1])+[top.nrhopndtscopies,top.nsndts]
      sss.reverse()
      tsourcep.shape = sss
      return transpose(tsourcep)
    elif lselfb is not None:
      return self.sourceparray[...,top.nrhopndtscopies*top.nsndts:]

  def setsourcepforparticles(self,isourcepndtscopies,indts,iselfb):
    if isourcepndtscopies is not None:
      sourcepndts = self.returnsourcep(lndts=1)
      self.sourcep = sourcepndts[...,isourcepndtscopies,indts]
    else:
      sourcepselfb = self.returnsourcep(lselfb=1)
      self.sourcep = sourcepselfb[...,iselfb]

  def setpotentialpforparticles(self,isourcepndtscopies,indts,iselfb):
    "Sets potential reference to the currently active potential"
    if indts is not None:
      potentialpndts = self.returnpotentialp(lndts=1)
      self.potentialp = potentialpndts[...,min(indts,w3d.nsndtsphi3d-1)]
    else:
      potentialpselfb = self.returnpotentialp(lselfb=1)
      self.potentialp = potentialpselfb[...,iselfb]

  def setfieldpforparticles(self,isourcepndtscopies,indts,iselfb):
    "Sets field reference to the currently active field"
    if indts is not None:
      fieldpndts = self.returnfieldp(lndts=1)
      if fieldpndts is not None:
        self.fieldp = fieldpndts[...,min(indts,w3d.nsndtsphi3d-1)]
      else:
        self.fieldp = 0.
    else:
      fieldpselfb = self.returnfieldp(lselfb=1)
      if fieldpselfb is not None:
        self.fieldp = fieldpselfb[...,iselfb]
      else:
        self.fieldp = 0.

  def setarraysforfieldsolve(self,isourcepndtscopies,indts,iselfb):
    if isourcepndtscopies is not None:
      sourcepndts = self.returnsourcep(lndts=1)
      fieldpndts = self.returnfieldp(lndts=1)
      potentialpndts = self.returnpotentialp(lndts=1)
      self.source = sourcepndts[...,isourcepndtscopies,indts]
      self.potential = potentialpndts[...,min(indts,w3d.nsndtsphi3d-1)]
      if fieldpndts is not None:
        self.field = fieldpndts[...,min(indts,w3d.nsndtsphi3d-1)]
      else:
        self.field = 0.
    else:
      sourcepselfb = self.returnsourcep(lselfb=1)
      fieldpselfb = self.returnfieldp(lselfb=1)
      potentialpselfb = self.returnpotentialp(lselfb=1)
      self.source = sourcepselfb[...,iselfb]
      self.potential = potentialpselfb[...,iselfb]
      if fieldpselfb is not None:
        self.field = fieldpselfb[...,iselfb]
      else:
        self.field = 0.

  def getpotentialpforparticles(self,isourcepndtscopies,indts,iselfb):
    "Copies from potential to potentialp"
    if indts is not None:
      potentialpndts = self.returnpotentialp(lndts=1)
      potentialpndts[...,min(indts,w3d.nsndtsphi3d-1)] = self.potential[...]
    else:
      potentialpselfb = self.returnpotentialp(lselfb=1)
      potentialpselfb[...,iselfb] = self.potential[...]

  def getfieldpforparticles(self,isourcepndtscopies,indts,iselfb):
    "Copies from field to fieldp"
    if indts is not None:
      fieldpndts = self.returnfieldp(lndts=1)
      fieldpndts[...,min(indts,w3d.nsndtsphi3d-1)] = self.field[...]
    else:
      fieldpselfb = self.returnfieldp(lselfb=1)
      fieldpselfb[...,iselfb] = self.field[...]

  def loadsource(self,lzero=None,**kw):
    'Charge deposition, uses particles from top directly'
    if lzero is None: lzero = w3d.lzerorhofsapi
    self.allocatedataarrays()
    if lzero: self.zerosourcep()

    for js,n in zip(arange(top.pgroup.ns),top.pgroup.nps):
      if n == 0: continue
      if top.pgroup.ldts[js]:
        indts = top.ndtstorho[top.pgroup.ndts[js]-1]
        self.setsourcepforparticles(0,indts,None)
        #zgrid = top.zgrid + (top.zgrid-top.zgridprv)*(top.pgroup.ndts[js]-1)
        self.setsourcep(js,top.pgroup,top.zgridndts[indts])
      if top.pgroup.iselfb[js] > -1:
        iselfb = top.pgroup.iselfb[js]
        self.setsourcepforparticles(None,None,iselfb)
        self.setsourcep(js,top.pgroup,top.zgrid)

    if lzero:
      for indts in range(top.nsndts):
        if top.ldts[indts]:
          self.setsourcepforparticles(0,indts,None)
          self.aftersetsourcep(lzero)
          self.makesourceperiodic()
      for iselfb in range(top.nsselfb):
        self.setsourcepforparticles(None,None,iselfb)
        self.aftersetsourcep(lzero)
        self.makesourceperiodic()

      self.averagesourcepwithsubcycling()

  def aftersetsourcep(self,lzero):
    "Anything that needs to be done to sourcep after the deposition"
    pass

  def loadrho(self,lzero=true,**kw):
    pass

  def loadj(self,lzero=true,**kw):
    pass

  def fetchfield(self,**kw):
    'Fetches the field, uses arrays from w3d module FieldSolveAPI'
    if w3d.npfsapi == 0: return
    # --- First, check how the data is being passed.
    if w3d.getpyobject('xfsapi') is not None:
      # --- If xfsapi is being used, pass the api arrays in directly and
      # --- don't pass in a pgroup. Note the try is used since
      # --- in many cases, either the E's or B's will not be associated,
      # --- in which case an exception is raised.
      try:    ex = w3d.exfsapi
      except: ex = zeros((0,), 'd')
      try:    ey = w3d.eyfsapi
      except: ey = zeros((0,), 'd')
      try:    ez = w3d.ezfsapi
      except: ez = zeros((0,), 'd')
      try:    bx = w3d.bxfsapi
      except: bx = zeros((0,), 'd')
      try:    by = w3d.byfsapi
      except: by = zeros((0,), 'd')
      try:    bz = w3d.bzfsapi
      except: bz = zeros((0,), 'd')
      args = [w3d.xfsapi,w3d.yfsapi,w3d.zfsapi,
              ex,ey,ez,bx,by,bz,w3d.getpyobject('pgroupfsapi')]
    else:
      # --- Otherwise, use data from w3d.pgroupfsapi.
      ipmin = w3d.ipminfsapi
      x = w3d.pgroupfsapi.xp[ipmin-1:ipmin-1+w3d.npfsapi]
      y = w3d.pgroupfsapi.yp[ipmin-1:ipmin-1+w3d.npfsapi]
      z = w3d.pgroupfsapi.zp[ipmin-1:ipmin-1+w3d.npfsapi]
      ex = w3d.pgroupfsapi.ex[ipmin-1:ipmin-1+w3d.npfsapi]
      ey = w3d.pgroupfsapi.ey[ipmin-1:ipmin-1+w3d.npfsapi]
      ez = w3d.pgroupfsapi.ez[ipmin-1:ipmin-1+w3d.npfsapi]
      bx = w3d.pgroupfsapi.bx[ipmin-1:ipmin-1+w3d.npfsapi]
      by = w3d.pgroupfsapi.by[ipmin-1:ipmin-1+w3d.npfsapi]
      bz = w3d.pgroupfsapi.bz[ipmin-1:ipmin-1+w3d.npfsapi]
      args = [x,y,z,ex,ey,ez,bx,by,bz,w3d.pgroupfsapi]

    js = w3d.jsfsapi
    if js < 0: ndtstosource = 0
    else:      ndtstosource = top.ndtstorho[top.pgroup.ndts[js]-1]

    tmpnsndts = getnsndtsforsubcycling()
    indts = min(tmpnsndts-1,ndtstosource)
    self.setpotentialpforparticles(None,indts,None)
    self.setfieldpforparticles(None,indts,None)
    self.fetchfieldfrompositions(*args,**kw)

    # add self B correction as needed
    if js >= 0 and top.pgroup.iselfb[js] > -1:
      extemp = zeros(shape(x),'d')
      eytemp = zeros(shape(y),'d')
      eztemp = zeros(shape(z),'d')
      self.setpotentialpforparticles(None,None,top.pgroup.iselfb[js])
      self.setfieldpforparticles(None,None,top.pgroup.iselfb[js])
      args = [x,y,z,extemp,eytemp,eztemp,None,None,None,w3d.pgroupfsapi]
      self.fetchfieldfrompositions(*args,**kw)
      ex[...] += extemp
      ey[...] += eytemp

  def fetchpotential(self):
    'Fetches the potential, uses arrays from w3d module FieldSolveAPI'
    if w3d.npfsapi == 0: return
    x = w3d.xfsapi
    y = w3d.yfsapi
    z = w3d.zfsapi
    # --- One of w3d.phifsapi or w3d.afsapi must be associated.
    try:    potential = w3d.phifsapi
    except: potential = w3d.afsapi

    js = w3d.jsfsapi
    if js < 0: ndtstosource = 0
    else:      ndtstosource = top.ndtstorho[top.pgroup.ndts[js]-1]
    tmpnsndts = getnsndtsforsubcycling()
    indts = min(tmpnsndts-1,ndtstosource)
    self.setpotentialpforparticles(None,indts,None)
    self.fetchpotentialfrompositions(x,y,z,potential)

  def selfbcorrection(self,js):
    "scale phispecies by -(1-1/gamma*2) store into top.fselfb"
    phipselfb = self.getpotentialp(lselfb=1)
    phipselfb[...,js] *= top.fselfb[js]

  def dosolveonpotential(self,iwhich,isourcepndtscopies,indts,iselfb):
    "points source and potential appropriately and call the solving routine"
    self.setarraysforfieldsolve(isourcepndtscopies,indts,iselfb)
    self.dosolve(iwhich,isourcepndtscopies,indts,iselfb)
    if iselfb is not None: self.selfbcorrection(iselfb)
    self.getpotentialpforparticles(isourcepndtscopies,indts,iselfb)

  def solve(self,iwhich=0):
    self.allocatedataarrays()
    # --- Loop over the subcyling groups and do any field solves that
    # --- are necessary.
    # --- Do loop in reverse order so that source and potential end up with the arrays
    # --- for the speices with the smallest timestep.
    tmpnsndts = getnsndtsforsubcycling()
    for indts in range(tmpnsndts-1,-1,-1):
      if (not top.ldts[indts] and
          (top.ndtsaveraging == 0 and not sum(top.ldts))): continue
      self.dosolveonpotential(iwhich,top.nrhopndtscopies-1,indts,None)

    # --- Solve for phi for groups which require the self B correction
    for js in range(top.nsselfb):
      self.dosolveonpotential(iwhich,None,None,js)

    gc.collect()

  def getpdims(self):
    raise """getpdims must be supplied - it should return a list of the dimensions
of the arrays used by the particles"""

  def allocatedataarrays(self):
    # --- Get base dimension of the arrays
    pdims = self.getpdims()

    # --- Setup arrays, including extra copies for subcycling
    # --- and self B corrections.
    setupSubcycling(top.pgroup)
    setupSelfB(top.pgroup)

    extrasourcepdim = top.nrhopndtscopies*top.nsndts + top.nsselfb
    dims = list(pdims[0]) + [extrasourcepdim]
    if 'sourceparray' not in self.__dict__ or shape(self.sourceparray) != tuple(dims):
      self.sourceparray = fzeros(dims,'d')

    extrapotentialdim = top.nsndtsphi + top.nsselfb
    dims = list(pdims[-1]) + [extrapotentialdim]
    if 'potentialparray' not in self.__dict__ or shape(self.potentialparray) != tuple(dims):
      self.potentialparray = fzeros(dims,'d')

    if len(pdims) == 3:
      # --- Also, create fieldparray
      dims = list(pdims[1]) + [extrapotentialdim]
      if 'fieldparray' not in self.__dict__ or shape(self.fieldparray) != tuple(dims):
        self.fieldparray = fzeros(dims,'d')

  def zerosourcep(self):
    if top.ndtsaveraging == 0:
      self.zerosourcepwithsampledsubcycling()
    elif top.ndtsaveraging == 1:
      self.zerosourcepwithfullvsubcycling()
    elif top.ndtsaveraging == 2:
      self.zerosourcepwithhalfvsubcycling()
    self.zerosourcepselfb()

  def zerosourcepwithsampledsubcycling(self):
    # --- Zero the sourcep copy for species when the positions
    # --- are advanced.
    # --- sourcepndts(...,2,in1) holds the old sourcep which is still needed
    # --- for the faster advanced groups.
    # --- Note the operation is faster on the transposed arrays (since the
    # --- looping is relative to the C ordering).
    sourcepndts = self.returnsourcep(lndts=1)
    for in1 in range(top.nsndts):
      if top.ldts[in1]:
        tsourcepndts = transpose(sourcepndts)
        if top.nrhopndtscopies == 2:
          #sourcepndts[...,1,in1] = sourcepndts[...,0,in1]
          tsourcepndts[in1,1,...] = tsourcepndts[in1,0,...]
        #sourcepndts[...,0,in1] = 0.
        tsourcepndts[in1,0,...] = 0.

  def zerosourcepwithfullvsubcycling(self):
    raise "fullv subcycling not yet implemented"

  def zerosourcepwithhalfvsubcycling(self):
    raise "halfv subcycling not yet implemented"

  def averagesourcepwithsubcycling(self):
    if top.ndtsaveraging == 0:
      self.averagesourcepwithsampledsubcycling()
    elif top.ndtsaveraging == 1:
      self.averagesourcepwithfullvsubcycling()
    elif top.ndtsaveraging == 2:
      self.averagesourcepwithhalfvsubcycling()

  def averagesourcepwithsampledsubcycling(self):
    if top.ndtsmax == 1: return

    sourcepndts = self.returnsourcep(lndts=1)

    # --- Note the operation is faster on the transposed arrays (since the
    # --- looping is relative to the C ordering).
    tsourcepndts = transpose(sourcepndts)

    # --- During the generate, do the copy of the new sourcep to the old sourcep
    # --- for group 0, which is normally done during the zerosourcep call.
    #sourcepndts[...,1,0] = sourcepndts[...,0,0]
    tsourcepndts[0,1,...] = tsourcepndts[0,0,...]

    # --- Save the sourcep where the fastest particle's sourcep is. For now,
    # --- assume that this is in1=0
    for in1 in range(1,top.nsndts):
      if top.it == 0:
        # --- At top.it==0, before the first step, always add the new sourcep.
        #sourcepndts[...,1,0] = (sourcepndts[...,1,0] + sourcepndts[...,0,in1])
        tsourcepndts[0,1,...] = (tsourcepndts[0,1,...] + tsourcepndts[in1,0,...])
      elif top.rhotondts[in1]%2 == 1:
        # --- Use the sourcep that is closest in time to the current time.
        if ((top.it-1)%top.rhotondts[in1] > top.rhotondts[in1]/2.-1.):
          #sourcepndts[...,1,0] = (sourcepndts[...,1,0] + sourcepndts[...,0,in1])
          tsourcepndts[0,1,...] = (tsourcepndts[0,1,...] + tsourcepndts[in1,0,...])
        else:
          #sourcepndts[...,1,0] = (sourcepndts[...,1,0] + sourcepndts[...,1,in1])
          tsourcepndts[0,1,...] = (tsourcepndts[0,1,...] + tsourcepndts[in1,1,...])
      else:
        # --- When ndts is even, at the mid point of the step, take the
        # --- average of the old and the new
        # --- Otherwise, use the sourcep that is closest in time to the current
        # --- time.
        if (top.it-1)%top.rhotondts[in1] == top.rhotondts[in1]/2-1:
          #sourcepndts[...,1,0] = (sourcepndts[...,1,0] +
          #     0.5*(sourcepndts[...,0,in1] + sourcepndts[...,1,in1]))
          tsourcepndts[0,1,...] = (tsourcepndts[0,1,...] +
               0.5*(tsourcepndts[in1,0,...] + tsourcepndts[in1,1,...]))
        elif ((top.it-1)%top.rhotondts[in1] > top.rhotondts[in1]/2.-1.):
          #sourcepndts[...,1,0] = (sourcepndts[...,1,0] + sourcepndts[...,0,in1])
          tsourcepndts[0,1,...] = (tsourcepndts[0,1,...] + tsourcepndts[in1,0,...])
        else:
          #sourcepndts[...,1,0] = (sourcepndts[...,1,0] + sourcepndts[...,1,in1])
          tsourcepndts[0,1,...] = (tsourcepndts[0,1,...] + tsourcepndts[in1,1,...])

  def zerosourcepselfb(self):
    sourcepselfb = self.returnsourcep(lselfb=1)
    tsourcepselfb = transpose(sourcepselfb)
    #sourcepselfb[...] = 0.
    tsourcepselfb[...] = 0.

