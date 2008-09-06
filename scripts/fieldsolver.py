"""Generic class describing the interface needed for a field solver.
"""
from __future__ import generators
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
  if (currpkg == "wxy"):
    loadrhoxy(pgroup,ins_i,nps_i,is_i,lzero)
  else:
    # --- Note that this works for all other packages, not just 3d.
    loadrho3d(pgroup,ins_i,nps_i,is_i,lzero)

#=============================================================================
def fieldsolve(iwhich=0,lbeforefs=false,lafterfs=false):
  """
This routine provides a simple call from the interpreter to do the fieldsol.
It calls the appropriate compiled interface routine based on the current
package. Only w3d and wxy have field solves defined.
 - iwhich=0: specifies what action to take
 - lbeforefs=false: when true, call functions installed be installbeforefs
 - lafterfs=false:  when true, call functions installed be installafterfs
  """
  if lbeforefs: controllers.beforefs()

  if top.fstype == 12:
    if iwhich > 0: return
    starttime = wtime()
    fieldsolregistered()
    endtime = wtime()
    top.fstime += endtime - starttime
  else:
    currpkg = package()[0]
    if (currpkg == "wxy"): fieldsolxy(iwhich)
    else:                  fieldsol3d(iwhich)

  if lafterfs: controllers.afterfs()

  # --- Now do extra work, updating arrays which depend directly on phi,
  # --- but only when a complete field solve was done.
  if iwhich == -1 or iwhich == 0:
    if (sometrue(top.efetch == 3) and top.fstype != 12 and
        (w3d.solvergeom == w3d.XYZgeom or
         w3d.solvergeom == w3d.RZgeom or
         w3d.solvergeom == w3d.XZgeom or
         w3d.solvergeom == w3d.Rgeom  or
         w3d.solvergeom == w3d.Zgeom)):
      getselfe3d(w3d.phip,w3d.nxp,w3d.nyp,w3d.nzp,w3d.selfe,
                 w3d.nx_selfe,w3d.ny_selfe,w3d.nz_selfe,
                 w3d.dx,w3d.dy,w3d.dz,true,1,1,1)
    # --- Get the phi needed for injection
    if top.inject > 0: getinj_phi()

# --- Define the old name.
fieldsol = fieldsolve

#=============================================================================
def fetche(pgroup=None,ipmin=None,ip=None,js=None):
  """Fetches the E field for particles in the given pgroup"""
  if pgroup is None: pgroup = top.pgroup
  if js is None: js = 0
  if ipmin is None: ipmin = pgroup.ins[js]
  if ip is None: ip = pgroup.nps[js]
  fetche3d(pgroup,ipmin,ip,js+1)

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
  if (currpkg == "wxy"):
    #loadrhoxy(ins_i,nps_i,is_i,lzero)
    print "loadj not support in wxy yet"
  else:
    loadj3d(pgroup,ins_i,nps_i,is_i,lzero)

#=============================================================================
# --- These routines are used to handle registered field solvers.
_registeredfieldsolvers = []
class RegisteredSolversContainer(object):
  """This class is needed so that the list of registered field solvers
can be properly restored after a restart. An instance of this object will be
saved in a restart dump. Upon restoration, it re-registers the solvers.  This
is needed since _registeredfieldsolvers will not be included in __main__ and
therefore would not be otherwise saved. Also, in warp.py, the instance of
this class is explicitly removed from the list of python variables that are
not saved in dump files.
  """
  def __init__(self,slist):
    # --- Note that self.slist should be a reference to _registeredfieldsolvers
    self.slist = slist
  def __setstate__(self,dict):
    self.__dict__.update(dict)
    # --- Clear out any solvers that are already registered
    for i in range(len(_registeredfieldsolvers)):
      del _registeredfieldsolvers[0]
    # --- Re-register the solvers.
    for solver in self.slist:
      registersolver(solver)
    # --- Get the original slist from this module, overwriting the one that
    # --- was saved in the dump file.
    self.slist = _registeredfieldsolvers
registeredsolverscontainer = RegisteredSolversContainer(_registeredfieldsolvers)

def registersolver(solver):
  """
Registers solvers to be used in the particle simulation.
 - solver: is the solver object. It must have the methods loadrho, solve, and
           fetche defined. fetchb and fetchphi will also be needed in some
           cases.

  """
  _registeredfieldsolvers.append(solver)
  top.fstype = 12
def getregisteredsolver(i=0):
  if len(_registeredfieldsolvers) == 0: return None
  return _registeredfieldsolvers[i]
def getregisteredsolvers():
  "Return the list of all registered field solver"
  # --- A copy is returned to prevent the list from being mucked up.
  return copy.copy(_registeredfieldsolvers)
def registeredsolvers():
  for solver in _registeredfieldsolvers:
    yield solver
def unregistersolver(solver=None,i=None):
  if i is not None:
    del _registeredfieldsolvers[i]
  else:
    if solver is None: solver = _registeredfieldsolvers[0]
    _registeredfieldsolvers.remove(solver)
def loadrhoregistered():
  assert len(_registeredfieldsolvers) > 0,"No solver has been registered"
  for f in _registeredfieldsolvers:
    f.loadrho()
def loadjregistered():
  assert len(_registeredfieldsolvers) > 0,"No solver has been registered"
  for f in _registeredfieldsolvers:
    f.loadj()
def fieldsolregistered():
  assert len(_registeredfieldsolvers) > 0,"No solver has been registered"
  for f in _registeredfieldsolvers:
    f.solve()
def fetcheregistered():
  assert len(_registeredfieldsolvers) > 0,"No solver has been registered"
  for f in _registeredfieldsolvers:
    f.fetche()
def fetchbregistered():
  assert len(_registeredfieldsolvers) > 0,"No solver has been registered"
  for f in _registeredfieldsolvers:
    f.fetchb()
def fetchphiregistered():
  assert len(_registeredfieldsolvers) > 0,"No solver has been registered"
  for f in _registeredfieldsolvers:
    f.fetchphi()
def fetcharegistered():
  assert len(_registeredfieldsolvers) > 0,"No solver has been registered"
  for f in _registeredfieldsolvers:
    f.fetcha()
def rhodiaregistered():
  assert len(_registeredfieldsolvers) > 0,"No solver has been registered"
  for f in _registeredfieldsolvers:
    f.rhodia()
def gtlchgregistered():
  assert len(_registeredfieldsolvers) > 0,"No solver has been registered"
  for f in _registeredfieldsolvers:
    f.gtlchg()
def srhoaxregistered():
  assert len(_registeredfieldsolvers) > 0,"No solver has been registered"
  for f in _registeredfieldsolvers:
    f.srhoax()
def geteseregistered():
  assert len(_registeredfieldsolvers) > 0,"No solver has been registered"
  for f in _registeredfieldsolvers:
    f.getese()
def sphiaxregistered():
  assert len(_registeredfieldsolvers) > 0,"No solver has been registered"
  for f in _registeredfieldsolvers:
    f.sphiax()
def sezaxregistered():
  assert len(_registeredfieldsolvers) > 0,"No solver has been registered"
  for f in _registeredfieldsolvers:
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
  _bfieldsolver[0].loadj()
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

  __w3dinputs__ = ['nx','ny','nz','dx','dy','dz','nzlocal','nzpguard',
                   'xmmin','xmmax','ymmin','ymmax','zmminlocal','zmmaxlocal',
                   'zmmin','zmmax',
                   'bound0','boundnz','boundxy','l2symtry','l4symtry',
                   'solvergeom']
  __topinputs__ = ['pbound0','pboundnz','pboundxy',
                   'my_index','nslaves','lfsautodecomp','zslave','lautodecomp',
                   'debug']
  __flaginputs__ = {'forcesymmetries':1,
                    'lreducedpickle':1,'lnorestoreonpickle':0,
                    'ldosolve':1,'l_internal_dosolve':1,
                    'gridvz':None,
                    }

  def __init__(self,**kw):
    try:
      kw['kwdict'].update(kw)
      kw = kw['kwdict']
      del kw['kwdict']
    except KeyError:
      pass

    # --- Save input parameters
    self.processdefaultsfrompackage(FieldSolver.__w3dinputs__,w3d,kw)
    self.processdefaultsfrompackage(FieldSolver.__topinputs__,top,kw)
    self.processdefaultsfromdict(FieldSolver.__flaginputs__,kw)

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
        self.bounds = zeros(6,'l')
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
        if self.solvergeom == w3d.RZgeom:
          self.bounds[0] = neumann
          self.bounds[2] = neumann
          self.bounds[3] = neumann
          if self.xmmin < 0.: self.xmmin = 0.
        elif self.solvergeom == w3d.XZgeom:
          self.bounds[2] = neumann
          self.bounds[3] = neumann

    # --- pbounds is special since it will sometimes be set from the
    # --- variables pbound0, pboundnz, pboundxy, l2symtry, and l4symtry
    if 'pbounds' not in self.__dict__:
      if 'pbounds' in kw:
        self.pbounds = kw['pbounds']
        del kw['pbounds']
      else:
        self.pbounds = zeros(6,'l')
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
        if self.solvergeom == w3d.RZgeom:
          self.pbounds[0] = reflect
          self.pbounds[2] = reflect
          self.pbounds[3] = reflect
        elif self.solvergeom == w3d.XZgeom:
          self.pbounds[2] = reflect
          self.pbounds[3] = reflect

    # --- Check for zero length dimensions
    if self.nx == 0:
      self.xmmin = 0.
      self.xmmax = 0.
    elif self.xmmin == self.xmmax:
      self.nx = 0
    if self.ny == 0:
      self.ymmin = 0.
      self.ymmax = 0.
    elif self.ymmin == self.ymmax:
      self.ny = 0
    if self.nz == 0:
      self.zmmin = 0.
      self.zmmax = 0.
    elif self.zmmin == self.zmmax:
      self.nz = 0

    # --- Set parallel related parameters and calculate mesh sizes
    self.lparallel = (self.nslaves > 1)
    if not self.lparallel:
      self.my_index = 0
      self.nzlocal = self.nz
      self.zmminlocal = self.zmmin
      self.zmmaxlocal = self.zmmax
      self.izfsslave = zeros(1,'l')
      self.nzfsslave = zeros(1,'l') + self.nz
      self.izpslave = self.izfsslave
      self.nzpslave = self.nzfsslave
      self.nxp = self.nx
      self.nyp = self.ny
      self.nzp = self.nz
      self.xmminp = self.xmmin
      self.xmmaxp = self.xmmax
      self.ymminp = self.ymmin
      self.ymmaxp = self.ymmax
      self.zmminp = self.zmminlocal
      self.zmmaxp = self.zmmaxlocal
    else:
      self.my_index = me
      self.izfsslave = zeros(self.nslaves,'l')
      self.nzfsslave = zeros(self.nslaves,'l')
      # --- Note that self.grid_overlap must be set by the inheriting class.
      top.grid_overlap = self.grid_overlap
      domaindecomposefields(self.nz,self.nslaves,self.lfsautodecomp,
                            self.izfsslave,self.nzfsslave,self.grid_overlap)

      self.nzlocal = self.nzfsslave[self.my_index]
      if self.dz == 0.: self.dz = (self.zmmax - self.zmmin)/self.nz
      self.zmminlocal = self.zmmin + self.izfsslave[self.my_index]*self.dz
      self.zmmaxlocal = self.zmmin + (self.izfsslave[self.my_index] + self.nzfsslave[self.my_index])*self.dz

      self.izpslave = zeros(self.nslaves,'l')
      self.nzpslave = zeros(self.nslaves,'l')
      # --- This should only be called after the particle decomposition
      # --- has been done.
      #self.setparticledomains()

    if self.dx == 0.: self.dx = (self.xmmax - self.xmmin)/self.nx
    if self.dy == 0.:
      if self.ny > 0: self.dy = (self.ymmax - self.ymmin)/self.ny
      else:           self.dy = self.dx
    if self.dz == 0.: self.dz = (self.zmmax - self.zmmin)/self.nz

    # --- Check the mesh consistency
    self.checkmeshconsistency(self.xmmin,self.xmmax,self.nx,self.dx,'x')
    self.checkmeshconsistency(self.ymmin,self.ymmax,self.ny,self.dy,'y')
    self.checkmeshconsistency(self.zmmin,self.zmmax,self.nz,self.dz,'z')
    self.checkmeshconsistency(self.zmminlocal,self.zmmaxlocal,self.nzlocal,self.dz,'z')

    self.xsymmetryplane = 0.
    self.ysymmetryplane = 0.
    self.xmesh = self.xmmin + arange(0,self.nx+1)*self.dx
    self.ymesh = self.ymmin + arange(0,self.ny+1)*self.dy
    self.zmesh = self.zmmin + arange(0,self.nz+1)*self.dz
    self.zmeshlocal = self.zmminlocal + arange(0,self.nzlocal+1)*self.dz

    self.ix_axis = nint(-self.xmmin/self.dx)
    self.iy_axis = nint(-self.ymmin/self.dy)
    self.iz_axis = nint(-self.zmmin/self.dz)

    # --- Some flags
    self.sourcepfinalized = 1

  def processdefaultsfrompackage(self,defaults,package,kw):
    for name in defaults:
      if name not in self.__dict__:
        #self.__dict__[name] = kw.pop(name,getattr(w3d,name)) # Python2.3
        self.__dict__[name] = kw.get(name,getattr(package,name))
      if kw.has_key(name): del kw[name]

  def processdefaultsfromdict(self,dict,kw):
    for name,defvalue in dict.iteritems():
      if name not in self.__dict__:
        #self.__dict__[name] = kw.pop(name,getattr(top,name)) # Python2.3
        self.__dict__[name] = kw.get(name,defvalue)
      if kw.has_key(name): del kw[name]

  def checkmeshconsistency(self,min,max,nn,dd,axis):
    'Checks if the mesh quantities are consistent'
    # --- Note that the factor of 1.e-5 is somewhat arbitrary
    assert abs((max-min) - nn*dd) < dd*1.e-5,\
      'The grid quantities along the '+axis+' axis are inconsistent'

  def __getstate__(self):
    dict = self.__dict__.copy()

# --- This is no longer needed since the list of registered solvers is
# --- now directly saved in a restart dump.
#   # --- Flag whether this is the registered solver so it knows whether
#   # --- to reregister itself upon the restore. The instance
#   # --- is not registered if it is not going to be restored.
#   if self in getregisteredsolvers() and not self.lnorestoreonpickle:
#     dict['iamtheregisteredsolver'] = 1
#   else:
#     dict['iamtheregisteredsolver'] = 0

    return dict

  def __setstate__(self,dict):
    self.__dict__.update(dict)

    # --- Set now z quantities is reading in an old dump file
    if 'nzlocal' not in self.__dict__:
      self.nzlocal = self.nz
      self.zmminlocal = self.zmmin
      self.zmmaxlocal = self.zmmax
      self.nz = self.nzfull
      self.zmmin = self.zmminglobal
      self.zmmax = self.zmmaxglobal
      del self.nzfull
      del self.zmminglobal
      del self.zmmaxglobal

    # --- Make sure that the new attribute l_internal_dosolve is defined.
    if 'l_internal_dosolve' not in self.__dict__:
      self.l_internal_dosolve = 1

    # --- This is only needed when an old dump file is being restored.
    # --- Now, the list of registered solvers is saved directly in the
    # --- dump file.
    if 'iamtheregisteredsolver' in self.__dict__:
      if self.iamtheregisteredsolver and not self.lnorestoreonpickle:
        del self.iamtheregisteredsolver
        registersolver(self)

  # ---------------------------------------------------------------------
  def advancezgrid(self):
    if self.gridvz is None: return
    # --- Advance the grid with its own velocity. This is only called
    # --- when gridvz is not None.
    try:                   self.itprevious
    except AttributeError: self.itprevious = top.it
    try:
      self._zgrid
    except AttributeError:
      self.setzgrid(top.zgrid)
    if self.itprevious < top.it:
      self.itprevious = top.it
      # --- This a new step, so advance zgrid
      self._zgrid += top.dt*self.gridvz
      self._zgridprv = self._zgrid

  def setzgrid(self,zgrid):
    if self.gridvz is None:
      self.gridvz = top.vbeamfrm
    self._zgrid = zgrid
    self._zgridprv = zgrid

  def getzgrid(self):
    if self.gridvz is None: return top.zgrid
    else:                   return self._zgrid

  def getzgridprv(self):
    if self.gridvz is None: return top.zgridprv
    else:                   return self._zgridprv

  def getzgridndts(self):
    if self.gridvz is None: return top.zgridndts
    else:                   return self._zgridndts

  def setgridvz(self,gridvz):
    self.gridvz = gridvz
    self._zgrid = top.zgrid
    self._zgridprv = top.zgrid

  # ---------------------------------------------------------------------
  # --- These routines must at least be defined.
  def loadrho(self,pgroup=None,lzero=true,**kw):
    'Charge deposition, uses particles from top directly'
    if pgroup is None: pgroup = top.pgroup
    self.advancezgrid()
    if lzero: self.zerorhop()

    for js,i,n,q,w in zip(arange(pgroup.ns),pgroup.ins-1,
                          pgroup.nps,pgroup.sq,pgroup.sw):
      if n > 0:
        self.setrhop(pgroup.xp[i:i+n],pgroup.yp[i:i+n],
                     pgroup.zp[i:i+n],pgroup.uzp[i:i+n],
                     q,w*pgroup.dtscale[js])

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

  def setparticledomains(self):
    if not self.lparallel: return

    # --- Set iz and nz. This is done so that zmesh[izpslave] < zpslmin, and
    # --- zmesh[izpslave+nzpslave] > zpslmax.
    # --- NOTE: There may be an issue in some cases with round-off since
    # --- sometimes (top.zpslmin - self.zmmin)/self.dz will be an integer.
    self.izpslave[:] = int((top.zpslmin - self.zmmin)/self.dz) - self.nzpguard
    self.nzpslave[:] = (int((top.zpslmax - self.zmmin)/self.dz) -
                       self.izpslave + 1 + 2*self.nzpguard)

    # --- Make sure that the processors don't have grid cells
    # --- sticking out the end.
    self.nzpslave[:] = where(self.izpslave<0,self.nzpslave+self.izpslave,
                                             self.nzpslave)
    self.izpslave[:] = where(self.izpslave<0,0,self.izpslave)
    self.nzpslave[:] = where(self.izpslave+self.nzpslave > self.nz,
                             self.nz - self.izpslave,
                             self.nzpslave)

    self.nxp = self.nx
    self.nyp = self.ny
    self.nzp = self.nzpslave[self.my_index]
    self.xmminp = self.xmmin
    self.xmmaxp = self.xmmax
    self.ymminp = self.ymmin
    self.ymmaxp = self.ymmax
    self.zmminp = self.zmmin + self.izpslave[self.my_index]*self.dz
    self.zmmaxp = self.zmminp + self.nzp*self.dz

    self.checkmeshconsistency(self.xmminp,self.xmmaxp,self.nxp,self.dx,'x')
    self.checkmeshconsistency(self.ymminp,self.ymmaxp,self.nyp,self.dy,'y')
    self.checkmeshconsistency(self.zmminp,self.zmmaxp,self.nzp,self.dz,'z')

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
      if self.lnorestoreonpickle or top.nrhopndtscopies == 1:
        # --- Note that the sourcep must be preserved if there
        # --- is subcycling since it will contain old data
        # --- which cannot be restored. If lnorestoreonpickle
        # --- is true, then sourcep can be deleted anyway since
        # --- it will not be used.
        if 'sourceparray' in dict: del dict['sourceparray']
      if 'potentialparray' in dict: del dict['potentialparray']
      if 'fieldparray' in dict: del dict['fieldparray']
      if 'sourcearray' in dict: del dict['sourcearray']
      if 'potentialarray' in dict: del dict['potentialarray']
      if 'fieldarray' in dict: del dict['fieldarray']
      if 'sourcep' in dict: del dict['sourcep']
      if 'potentialp' in dict: del dict['potentialp']
      if 'fieldp' in dict: del dict['fieldp']
      if 'source' in dict: del dict['source']
      if 'potential' in dict: del dict['potential']
      if 'field' in dict: del dict['field']
    return dict

  def __setstate__(self,dict):
    FieldSolver.__setstate__(self,dict)
    if 'sourceparray' in dict:
      self.sourceparray = makefortranordered(self.sourceparray)
    if self.lreducedpickle and not self.lnorestoreonpickle:
      installafterrestart(self.allocatedataarrays)

  def returnfieldp(self,indts,iselfb):
    indts = min(indts,top.nsndtsphi-1)
    try:
      return self.fieldparray[...,indts]
    except AttributeError:
      return 0.

  def returnpotentialp(self,indts,iselfb):
    indts = min(indts,top.nsndtsphi-1)
    try:
      return self.potentialparray[...,indts,iselfb]
    except AttributeError:
      return 0.

  def returnsourcep(self,isourcepndtscopies,indts,iselfb):
    try:
      return self.sourceparray[...,isourcepndtscopies,indts,iselfb]
    except AttributeError:
      return 0.

  def returnfield(self,indts,iselfb):
    indts = min(indts,top.nsndtsphi-1)
    try:
      return self.fieldarray[...,indts]
    except AttributeError:
      return 0.

  def returnpotential(self,indts,iselfb):
    indts = min(indts,top.nsndtsphi-1)
    try:
      return self.potentialarray[...,indts,iselfb]
    except AttributeError:
      return 0.

  def returnsource(self,indts,iselfb):
    indts = min(indts,top.nsndtsphi-1)
    try:
      return self.sourcearray[...,indts,iselfb]
    except AttributeError:
      return 0.

  def setsourcepforparticles(self,isourcepndtscopies,indts,iselfb):
    self.sourcep = self.returnsourcep(isourcepndtscopies,indts,iselfb)

  def setpotentialpforparticles(self,isourcepndtscopies,indts,iselfb):
    "Sets potential reference to the currently active potential"
    self.potentialp = self.returnpotentialp(indts,iselfb)

  def setfieldpforparticles(self,isourcepndtscopies,indts,iselfb):
    "Sets field reference to the currently active field"
    self.fieldp = self.returnfieldp(indts,iselfb)

  def setsourceforfieldsolve(self,isourcepndtscopies,indts,iselfb):
    # --- This is called at the end of loadrho just before the b.c.'s are set
    self.source = self.returnsource(indts,iselfb)

  def setarraysforfieldsolve(self,isourcepndtscopies,indts,iselfb):
    # --- This is called at the beginning of the field solve
    self.source    = self.returnsource(indts,iselfb)
    self.field     = self.returnfield(indts,iselfb)
    self.potential = self.returnpotential(indts,iselfb)

  def getpotentialpforparticles(self,isourcepndtscopies,indts,iselfb):
    "Copies from potential to potentialp"
    # --- In the serial case, potentialp and self.potential point to the
    # --- same memory.
    if lparallel:
      potentialp = self.returnpotentialp(indts,iselfb)
      potentialp[...] = self.potential

  def getfieldpforparticles(self,isourcepndtscopies,indts,iselfb):
    "Copies from field to fieldp"
    # --- In the serial case, fieldp and self.field point to the
    # --- same memory.
    if lparallel:
      fieldp = self.returnfieldp(indts,iselfb)
      fieldp[...] = self.field

  # ---------------------------------------------------------------------
  def setupzgridndts(self):
    # --- Check to len if zgridndts and update it appropriately if it
    # --- has changed.
    if len(self._zgridndts) < top.nsndts:
      # --- Update ndtstozgrid
      ndtstozgridnew = zeros(top.ndtsmax,'d')
      nn = min(top.ndtsmax,len(self._ndtstozgrid))
      ndtstozgridnew[:nn] = self._ndtstozgrid[:nn]
      ndtstozgridnew[nn:] = self._zgrid
      self._ndtstozgrid = ndtstozgridnew
      # --- Now get zgridndts
      self._zgridndts = take(self._ndtstozgrid,top.ndts-1)

  def advancezgrid(self):
    if self.gridvz is None: return
    # --- This routine is called at the start of the loadsource routine
    # --- Get initial values from top
    try:
      self.itprevious
    except AttributeError:
      self.itprevious = top.it
    try:
      self._zgrid
    except AttributeError:
      self.setzgrid(top.zgrid)
    self.setupzgridndts()
    if self.itprevious < top.it:
      # --- This a new step, so advance zgrid.
      self.itprevious = top.it
      # --- Note that zgridprv is set
      # --- to the advanced value of zgrid. This is done since this happens
      # --- at the start of a loadsource. This loadsource only uses zgrid
      # --- (actually zgridndts), but sete3d needs zgridprv. In a normal step,
      # --- zgridprv is set to zgrid at the end of the particle advance.
      # --- Setting it here is equivalent, since zgridprv is not used anyway
      # --- until the next call to sete3d. A more precise version would set
      # --- zgridprv after the field solve when the fields are then aligned
      # --- with the rho (which is at zgrid).
      if top.nsndts > 0:
        for ndts in top.ndts:
          if (top.it-1)%ndts == 0:
            self._ndtstozgrid[ndts-1] = (self._zgrid + top.dt*self.gridvz*ndts)
        self._zgridndts = take(self._ndtstozgrid,top.ndts-1)
      self._zgrid += top.dt*self.gridvz
      self._zgridprv = self._zgrid

  def setzgrid(self,zgrid):
    if self.gridvz is None:
      self.gridvz = top.vbeamfrm
    self._zgrid = zgrid
    self._zgridprv = zgrid
    self._zgridndts = []
    self._ndtstozgrid = []
    self.setupzgridndts()

  def getzgrid(self):
    if self.gridvz is None: return top.zgrid
    else:                   return self._zgrid

  def getzgridprv(self):
    if self.gridvz is None: return top.zgridprv
    else:                   return self._zgridprv

  def getzgridndts(self):
    if self.gridvz is None: return top.zgridndts
    else:                   return self._zgridndts

  def setgridvz(self,gridvz):
    self.gridvz = gridvz
    self._zgrid = top.zgrid
    self._zgridprv = top.zgrid
    self._zgridndts = []
    self._ndtstozgrid = []
    self.setupzgridndts()

  # ---------------------------------------------------------------------
  def loadsource(self,lzero=None,pgroups=None,**kw):
    'Charge deposition, uses particles from top directly'
    # --- Note that the grid location is advanced even if no field solve
    # --- is being done.
    self.advancezgrid()
    # --- If ldosolve is false, then skip the gather of rho, unless
    # --- lzero is also false, in which case the solver is assumed to
    # --- be gathering the source (for example during an EGUN iteration).
    if not self.ldosolve and lzero: return
    if lzero is None: lzero = w3d.lzerorhofsapi

    self.setparticledomains()
    self.allocatedataarrays()
    if lzero: self.zerosourcep()

    if pgroups is None: pgroups = [top.pgroup]
    for pgroup in pgroups:

      if w3d.js1fsapi >= 0: js1 = w3d.js1fsapi
      else:                 js1 = 0
      if w3d.js2fsapi >= 0: js2 = w3d.js2fsapi+1
      else:                 js2 = pgroup.ns

      for js in range(js1,js2):
        n = pgroup.nps[js]
        if n == 0: continue
        if pgroup.ldts[js]:
          indts = top.ndtstorho[pgroup.ndts[js]-1]
          iselfb = pgroup.iselfb[js]
          self.setsourcepforparticles(0,indts,iselfb)

          if self.debug:
            i1 = pgroup.ins[js]-1
            i2 = pgroup.ins[js]+pgroup.nps[js]-1
            if self.nx > 0:
              x = pgroup.xp[i1:i2]
              assert min(abs(x-self.xmmin)) >= 0.,\
                     "Particles in species %d have x below the grid when depositing the source, min x = %e"%(js,min(x))
              assert max(x) < self.xmmax,\
                     "Particles in species %d have x above the grid when depositing the source, max x = %e"%(js,max(x))
            if self.ny > 0:
              y = pgroup.yp[i1:i2]
              assert min(abs(y-self.ymmin)) >= 0.,\
                     "Particles in species %d have y below the grid when depositing the source, min y = %e"%(js,min(y))
              assert max(y) < self.ymmax,\
                     "Particles in species %d have y above the grid when depositing the source, max y = %e"%(js,max(y))
            if self.nzlocal > 0:
              z = pgroup.zp[i1:i2]
              assert min(z) >= self.zmminp+self.getzgridndts()[indts],\
                     "Particles in species %d have z below the grid when depositing the source, min z = %e"%(js,min(z))
              assert max(z) < self.zmmaxp+self.getzgridndts()[indts],\
                     "Particles in species %d have z above the grid when depositing the source, max z = %e"%(js,max(z))

          self.setsourcep(js,pgroup,self.getzgridndts()[indts])

    # --- Only finalize the source if lzero is true, which means the this
    # --- call to loadsource should be a complete operation.
    self.sourcepfinalized = 0
    if lzero: self.finalizesourcep()

  def finalizesourcep(self):
    if self.sourcepfinalized: return
    self.sourcepfinalized = 1
    for indts in range(top.nsndts):
      if top.ldts[indts]:
        for iselfb in range(top.nsselfb):
          self.setsourcepforparticles(0,indts,iselfb)
          self.aftersetsourcep()

    self.averagesourcepwithsubcycling()

    tmpnsndts = getnsndtsforsubcycling()
    for indts in range(tmpnsndts-1,-1,-1):
      if (not top.ldts[indts] and
          ((top.ndtsaveraging == 0 or top.ndtsaveraging == 1)
           and not sum(top.ldts))): cycle
      for iselfb in range(top.nsselfb):
        isndts = min(indts,top.nsndtsphi)
        self.setsourceforfieldsolve(top.nrhopndtscopies-1,isndts,iselfb)
        self.makesourceperiodic()

  def aftersetsourcep(self):
    "Anything that needs to be done to sourcep after the deposition"
    pass

  def loadrho(self,lzero=true,**kw):
    pass

  def loadj(self,lzero=true,**kw):
    pass

  def fetchfield(self,*args,**kw):
    'Fetches the field, uses arrays from w3d module FieldSolveAPI'
    if w3d.npfsapi == 0: return
    # --- First, check how the data is being passed.
    if w3d.getpyobject('xfsapi') is not None:
      # --- If xfsapi is being used, pass the api arrays in directly and
      # --- don't pass in a pgroup. Note the try is used since
      # --- in many cases, either the E's or B's will not be associated,
      # --- in which case an exception is raised.
      x = w3d.xfsapi
      y = w3d.yfsapi
      z = w3d.zfsapi
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
      pgroup = w3d.getpyobject('pgroupfsapi')
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
      pgroup = w3d.pgroupfsapi

    jsid = w3d.jsfsapi
    if jsid < 0: js = 0
    else:        js = jsid

    if self.debug and top.efetch[js] != 5:
      if self.nx > 0:
        assert min(abs(x-self.xmmin)) >= 0.,\
               "Particles in species %d have x below the grid when fetching the field"%jsid
        assert max(x) < self.xmmax,\
               "Particles in species %d have x above the grid when fetching the field"%jsid
      if self.ny > 0:
        assert min(abs(y-self.ymmin)) >= 0.,\
               "Particles in species %d have y below the grid when fetching the field"%jsid
        assert max(y) < self.ymmax,\
               "Particles in species %d have y above the grid when fetching the field"%jsid
      if self.nzlocal > 0:
        assert min(z) >= self.zmminp+self.getzgridprv(),\
               "Particles in species %d have z below the grid when fetching the field"%jsid
        assert max(z) < self.zmmaxp+self.getzgridprv(),\
               "Particles in species %d have z above the grid when fetching the field"%jsid

    args = [x,y,z,ex,ey,ez,bx,by,bz,jsid,pgroup]

    if jsid < 0: indts = 0
    else:        indts = top.ndtstorho[w3d.ndtsfsapi-1]

    tmpnsndts = getnsndtsforsubcycling()
    indts = min(tmpnsndts-1,indts)
    iselfb = top.iselfb[jsid]
    self.setpotentialpforparticles(None,indts,iselfb)
    self.setfieldpforparticles(None,indts,iselfb)
    self.fetchfieldfrompositions(*args)

  def fetchpotential(self,*args,**kw):
    'Fetches the potential, uses arrays from w3d module FieldSolveAPI'
    if w3d.npfsapi == 0: return
    x = w3d.xfsapi
    y = w3d.yfsapi
    z = w3d.zfsapi
    # --- One of w3d.phifsapi or w3d.afsapi must be associated.
    try:    potential = w3d.phifsapi
    except: potential = w3d.afsapi

    if self.debug:
      if self.nx > 0:
        assert min(abs(x-self.xmmin)) >= 0.,\
               "Particles have x below the grid when fetching the potential"
        assert max(x) <= self.xmmax,\
               "Particles have x above the grid when fetching the potential"
      if self.ny > 0:
        assert min(abs(y-self.ymmin)) >= 0.,\
               "Particles have y below the grid when fetching the potential"
        assert max(y) <= self.ymmax,\
               "Particles have y above the grid when fetching the potential"
      if self.nzlocal > 0:
        assert min(z) >= self.zmminlocal,\
               "Particles have z below the grid when fetching the potential"
        assert max(z) <= self.zmmaxlocal,\
               "Particles have z above the grid when fetching the potential"

    jsid = w3d.jsfsapi
    if jsid < 0: indts = 0
    else:        indts = top.ndtstorho[w3d.ndtsfsapi-1]

    tmpnsndts = getnsndtsforsubcycling()
    indts = min(tmpnsndts-1,indts)
    iselfb = top.iselfb[jsid]
    self.setpotentialpforparticles(None,indts,iselfb)
    self.fetchpotentialfrompositions(x,y,z,potential)

  def dosolveonpotential(self,iwhich,isourcepndtscopies,indts,iselfb):
    "points source and potential appropriately and call the solving routine"
    self.setarraysforfieldsolve(isourcepndtscopies,indts,iselfb)
    self.dosolve(iwhich,isourcepndtscopies,indts,iselfb)
    self.getpotentialpforparticles(isourcepndtscopies,indts,iselfb)

  def solve(self,iwhich=0):
    if not self.ldosolve: return
    self.allocatedataarrays()
    # --- This is only needed in cases when the source is accumulated over
    # --- multiple steps, and can only be finalized (e.g. made periodic)
    # --- at this point.
    self.finalizesourcep()
    # --- Loop over the subcyling groups and do any field solves that
    # --- are necessary.
    # --- Do loop in reverse order so that source and potential end up
    # --- with the arrays
    # --- for the species with the smallest timestep.
    tmpnsndts = getnsndtsforsubcycling()
    for indts in range(tmpnsndts-1,-1,-1):
      if (not top.ldts[indts] and
          (top.ndtsaveraging == 0 and not sum(top.ldts))): continue
      # --- Note that the field solve is done even if there are no species
      # --- (i.e. when top.nsselfb==0)
      for iselfb in range(max(1,top.nsselfb)-1,-1,-1):
        self.dosolveonpotential(iwhich,top.nrhopndtscopies-1,indts,iselfb)

    # --- Is this still needed? It seems to slow things down alot.
    #gc.collect()

  def getallpotentialpforparticles(self,iwhich=0,lforce=0):
    "This transfers data from the potential array to the potentialp array"
    if not self.ldosolve and not lforce: return
    self.allocatedataarrays()
    # --- Loop over the subcyling groups and get any potentialp that
    # --- are necessary.
    # --- Do loop in reverse order so that source and potential end up
    # --- with the arrays
    # --- for the species with the smallest timestep.
    tmpnsndts = getnsndtsforsubcycling()
    for indts in range(tmpnsndts-1,-1,-1):
      if (not top.ldts[indts] and
          (top.ndtsaveraging == 0 and not sum(top.ldts))): continue
      for iselfb in range(top.nsselfb-1,-1,-1):
        self.setarraysforfieldsolve(top.nrhopndtscopies-1,indts,iselfb)
        self.getpotentialpforparticles(top.nrhopndtscopies-1,indts,iselfb)

  def getpdims(self):
    raise """getpdims must be supplied - it should return a list of the dimensions
of the arrays used by the particles"""

  def getdims(self):
    raise """getdims must be supplied - it should return a list of the dimensions
of the arrays used by the field solve"""

  def allocatedataarrays(self):
    # --- Setup arrays, including extra copies for subcycling
    # --- and self B corrections.
    setupSubcycling(top.pgroup)
    setupSelfB(top.pgroup)

    # --- Get base dimension of the arrays for the particles
    pdims = self.getpdims()

    # --- This ensures that the arrays would be set up for a field solve even
    # --- if there were no species (and top.nsselfb==0).
    nsselfb = max(1,top.nsselfb)

    sourcedims = list(pdims[0]) + [top.nrhopndtscopies,top.nsndts,nsselfb]
    if 'sourceparray' not in self.__dict__ or shape(self.sourceparray) != tuple(sourcedims):
      self.sourceparray = fzeros(sourcedims,'d')

    potentialdims = list(pdims[-1]) + [top.nsndtsphi,nsselfb]
    if 'potentialparray' not in self.__dict__ or shape(self.potentialparray) != tuple(potentialdims):
      self.potentialparray = fzeros(potentialdims,'d')

    if len(pdims) == 3:
      # --- Also, create fieldparray
      fielddims = list(pdims[1]) + [top.nsndtsphi]
      if 'fieldparray' not in self.__dict__ or shape(self.fieldparray) != tuple(fielddims):
        self.fieldparray = fzeros(fielddims,'d')

    if not self.lparallel:
      # --- For the serial case, the array for the field solve is the same as
      # --- the array for the particles.
      self.sourcearray = self.sourceparray[...,top.nrhopndtscopies-1,:,:]
      self.potentialarray = self.potentialparray[...,:,:]
      if len(pdims) == 3:
        self.fieldarray = self.fieldparray[...,:,:]
    else:
      # --- In parallel, the arrays for the field solver are separate arrays

      # --- Get base dimension of the arrays for the field solver
      dims = self.getdims()

      sourcedims = list(dims[0]) + [top.nsndtsphi,nsselfb]
      if 'sourcearray' not in self.__dict__ or shape(self.sourcearray) != tuple(sourcedims):
        self.sourcearray = fzeros(sourcedims,'d')

      potentialdims = list(dims[-1]) + [top.nsndtsphi,nsselfb]
      if 'potentialarray' not in self.__dict__ or shape(self.potentialarray) != tuple(potentialdims):
        self.potentialarray = fzeros(potentialdims,'d')

      if len(dims) == 3:
        # --- Also, create fieldarray
        fielddims = list(dims[1]) + [top.nsndtsphi]
        if 'fieldarray' not in self.__dict__ or shape(self.fieldarray) != tuple(fielddims):
          self.fieldarray = fzeros(fielddims,'d')

  def resetparticledomains(self):
    self.setparticledomains()
    # --- Note that this will call allocatedataarrays.
    self.getallpotentialpforparticles(lforce=1)
    # --- Make sure the sourcep gets set to the updated sourceparray.
    self.setsourcepforparticles(0,0,0)

  def zerosourcep(self):
    if top.ndtsaveraging == 0:
      self.zerosourcepwithsampledsubcycling()
    elif top.ndtsaveraging == 1:
      self.zerosourcepwithfullvsubcycling()
    elif top.ndtsaveraging == 2:
      self.zerosourcepwithhalfvsubcycling()

  def zerosourcepwithsampledsubcycling(self):
    # --- Zero the sourcep copy for species when the positions
    # --- are advanced.
    # --- sourcepndts(...,2,indts,:) holds the old sourcep which is still needed
    # --- for the faster advanced groups.
    # --- Note the operation is faster on the transposed arrays (since the
    # --- looping is relative to the C ordering).
    tsourcep = transpose(self.sourceparray)
    for indts in range(top.nsndts):
      if top.ldts[indts]:
        if top.nrhopndtscopies == 2:
          tsourcep[:,indts,1,...] = tsourcep[:,indts,0,...]
        tsourcep[:,indts,0,...] = 0.

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

    # --- Note the operation is faster on the transposed arrays (since the
    # --- looping is relative to the C ordering).
    tsourcep = transpose(self.sourceparray)

    # --- Do the copy of the new sourcep to the old sourcep for group 0,
    # --- the fastest group. Note that the old rho for this group is never
    # --- used so that space in the array is used during the field solve.
    tsourcep[:,0,1,...] = tsourcep[:,0,0,...]

    # --- Save the sourcep where the fastest particle's sourcep is. For now,
    # --- assume that this is in1=0
    for in1 in range(1,top.nsndts):
      if top.it == 0:
        # --- At top.it==0, before the first step, always add the new sourcep.
        tsourcep[:,0,1,...] = tsourcep[:,0,1,...] + tsourcep[:,in1,0,...]
      elif top.ndts[in1]%2 == 1:
        # --- Use the sourcep that is closest in time to the current time.
        if ((top.it-1)%top.ndts[in1] > top.ndts[in1]/2.-1.):
          tsourcep[:,0,1,...] = tsourcep[:,0,1,...] + tsourcep[:,in1,0,...]
        else:
          tsourcep[:,0,1,...] = tsourcep[:,0,1,...] + tsourcep[:,in1,1,...]
      else:
        # --- When ndts is even, at the mid point of the step, take the
        # --- average of the old and the new
        # --- Otherwise, use the sourcep that is closest in time to the current
        # --- time.
        if (top.it-1)%top.ndts[in1] == top.ndts[in1]/2-1:
          tsourcep[:,0,1,...] = (tsourcep[:,0,1,...] +
               0.5*(tsourcep[:,in1,0,...] + tsourcep[:,in1,1,...]))
        elif ((top.it-1)%top.ndts[in1] > top.ndts[in1]/2.-1.):
          tsourcep[:,0,1,...] = tsourcep[:,0,1,...] + tsourcep[:,in1,0,...]
        else:
          tsourcep[:,0,1,...] = tsourcep[:,0,1,...] + tsourcep[:,in1,1,...]

  def saveprevioussource(self,lfinalize=1):
    # --- This is needed by the EGUN method, which needs the previous rho. Note
    # --- that the subycling and selfb are ignored here since those models
    # --- don't make sense with the EGUN mode.
    # --- Do any finalization calculation on the source if requested.
    if lfinalize: self.finalizesourcep()
    self.sourceprevious = self.returnsource(0,0).copy()

  def averagewithprevioussource(self,param,lfinalize=1):
    # --- This is used by the EGUN method, to average the source over multiple
    # --- iterations.
    # --- Do any finalization calculation on the source if requested.
    if lfinalize: self.finalizesourcep()
    source = self.returnsource(0,0)
    source[...] = (1.-param)*source + param*self.sourceprevious

  def savepreviouspotential(self):
    # --- This is needed by the implicit algorithm.
    self.potentialprevious = self.potentialarray.copy()

  def addinpreviouspotential(self):
    # --- This is used by the implicit algorithm, to average the potential
    # --- over multiple iterations.
    self.potentialarray += self.potentialprevious

  def debugparticlebounds(self,text=''):
    pgroups = [top.pgroup]
    for pgroup in pgroups:

      for js in range(pgroup.ns):
        n = pgroup.nps[js]
        if n == 0: continue
        indts = top.ndtstorho[pgroup.ndts[js]-1]

        i1 = pgroup.ins[js]-1
        i2 = pgroup.ins[js]+pgroup.nps[js]-1
        if self.nx > 0:
          x = pgroup.xp[i1:i2]
          assert min(abs(x-self.xmmin)) >= 0.,\
                 text+"Particles in species %d have x below the grid when depositing the source, min x = %e"%(js,min(x))
          assert max(x) < self.xmmax,\
                 text+"Particles in species %d have x above the grid when depositing the source, min x = %e"%(js,max(x))
        if self.ny > 0:
          y = pgroup.yp[i1:i2]
          assert min(abs(y-self.ymmin)) >= 0.,\
                 text+"Particles in species %d have y below the grid when depositing the source, min x = %e"%(js,min(y))
          assert max(y) < self.ymmax,\
                 text+"Particles in species %d have y above the grid when depositing the source, min x = %e"%(js,max(y))
        if self.nzlocal > 0:
          z = pgroup.zp[i1:i2]
          assert min(z) >= self.zmminp+self.getzgridndts()[indts],\
                 text+"Particles in species %d have z below the grid when depositing the source, min x = %e"%(js,min(z))
          assert max(z) < self.zmmaxp+self.getzgridndts()[indts],\
                 text+"Particles in species %d have z above the grid when depositing the source, min x = %e"%(js,max(z))

