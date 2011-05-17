"""Generic class describing the interface needed for a field solver.
"""
from __future__ import generators
from warp import *
import __main__
import gc

fieldsolver_version = "$Id: fieldsolver.py,v 1.96 2011/05/17 19:19:12 grote Exp $"

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
      allocateselfepforparticles(true)
      getselfe3d(w3d.phip,w3d.nxp,w3d.nyp,w3d.nzp,w3d.selfep,
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
  currpkg = package()[0]
  if (currpkg == "wxy"):
    fetchexy(pgroup,ipmin,ip,js+1,top.pgroup.ex,top.pgroup.ey,top.pgroup.ez)
  else:
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
def _checkfstypeforinvalidvalue():
  if len(_registeredfieldsolvers) > 0 and top.fstype != 12:
    raise RuntimeError('top.fstype should not be changed when a field solver has been registered')
def getregisteredsolver(i=0):
  if len(_registeredfieldsolvers) == 0: return None
  _checkfstypeforinvalidvalue()
  return _registeredfieldsolvers[i]
def getregisteredsolvers():
  "Return the list of all registered field solver"
  _checkfstypeforinvalidvalue()
  # --- A copy is returned to prevent the list from being mucked up.
  return copy.copy(_registeredfieldsolvers)
def registeredsolvers():
  _checkfstypeforinvalidvalue()
  for solver in _registeredfieldsolvers:
    yield solver
def unregistersolver(solver=None,i=None):
  _checkfstypeforinvalidvalue()
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

#=============================================================================
#=============================================================================
class W3DFieldsolver(object):
  """
This is a wrapper class around the built in field solver in the w3d package
(with some references to top as well). This class stores no information
itself, but always gets data from the Fortran packages.
  """
  def __getattr__(self,name):
    "Get any names from w3d or top"
    # --- Maybe this should be limited to those related to the field solver?
    # --- First check w3d
    try:
      return getattr(w3d,name)
    except AttributeError:
      pass
    # --- Then check top, raising AttributeError is not found.
    return getattr(top,name)

  def __setattr__(self,name,val):
    "Set any names from w3d or top"
    # --- Maybe this should be limited to those related to the field solver?
    # --- First check w3d
    try:
      return setattr(w3d,name,val)
    except AttributeError:
      pass
    # --- Then check top, raising AttributeError is not found.
    return getattr(top,name,val)

  def __getstate__(self):
    # --- Nothing needs to be included in the pickle
    return {}

  def getrho(self):
    if solver.solvergeom in [w3d.RZgeom,w3d.XZgeom,w3d.Zgeom]:
      return frz.basegrid.rho
    else:
      return w3d.rho

  def getphi(self):
    if solver.solvergeom in [w3d.RZgeom,w3d.XZgeom,w3d.Zgeom]:
      return frz.basegrid.phi[1:-1,1:-1]
    else:
      return w3d.phi[1:-1,1:-1,1:-1]

  def getrhop(self):
    if solver.solvergeom in [w3d.RZgeom,w3d.XZgeom,w3d.Zgeom]:
      return frz.basegrid.rhop
    else:
      return w3d.rhop

  def getphip(self):
    if solver.solvergeom in [w3d.RZgeom,w3d.XZgeom,w3d.Zgeom]:
      return frz.basegrid.phi[1:-1,1:-1]
    else:
      return w3d.phip[1:-1,1:-1,1:-1]

  def resetparticledomains(self):
    if(w3d.solvergeom == w3d.XYZgeom):
      w3d.nxp = ppdecomp.nx[top.ixproc]
      w3d.nyp = ppdecomp.ny[top.iyproc]
      w3d.nzp = ppdecomp.nz[top.izproc]
      gchange("Fields3dParticles")
    else:
      gchange_rhop_phip_rz()
   ## --- Redistribute phi to the particle arrays if a field solve is
   ## --- not done.
   #if not dofs:
   #  if getregisteredsolver() is None:
   #    for i in range(getnsndtsforsubcycling()):
   #      getphipforparticles(i)

  def zerosource(self):
    if solver.solvergeom in [w3d.RZgeom,w3d.XZgeom,w3d.Zgeom]:
      frz.basegrid.rho = 0.
    else:
      w3d.rho = 0.

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

  __w3dinputs__ = ['nx','ny','nz','dx','dy','dz','nxlocal','nylocal','nzlocal',
                   'nxpextra','nypextra','nzpextra',
                   'xmmin','xmmax','ymmin','ymmax','zmmin','zmmax',
                   'xmminlocal','xmmaxlocal',
                   'ymminlocal','ymmaxlocal',
                   'zmminlocal','zmmaxlocal',
                   'bound0','boundnz','boundxy','l2symtry','l4symtry',
                   'solvergeom']
  __topinputs__ = ['pbound0','pboundnz','pboundxy',
                   'nprocs','nxprocs','nyprocs','nzprocs',
                   'userdecompx','userdecompy','userdecompz','lautodecomp',
                   'lfsautodecomp','zslave','debug']
  __flaginputs__ = {'forcesymmetries':1,
                    'lreducedpickle':1,'lnorestoreonpickle':0,
                    'ldosolve':1,'l_internal_dosolve':1,
                    'gridvz':None,'lchild':False,
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

    if self.solvergeom == w3d.RZgeom:
      # --- Turn off these flags just in case...
      self.l2symtry = false
      self.l4symtry = false

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

    if self.dx == 0.: self.dx = (self.xmmax - self.xmmin)/self.nx
    if self.dy == 0.:
      if self.ny > 0: self.dy = (self.ymmax - self.ymmin)/self.ny
      else:           self.dy = self.dx
    if self.dz == 0.: self.dz = (self.zmmax - self.zmmin)/self.nz
    # --- Set parallel related parameters and calculate mesh sizes
    self.lparallel = (self.nprocs>1)
    if not self.lparallel or self.lchild:
      self.my_index = 0
      self.nprocs = 1
      self.nxprocs = 1
      self.nyprocs = 1
      self.nzprocs = 1
      self.ixproc = self.ixproc = self.ixproc = 0
      self.nxlocal = self.nx
      self.nylocal = self.ny
      self.nzlocal = self.nz
      self.xmminlocal = self.xmmin
      self.xmmaxlocal = self.xmmax
      self.ymminlocal = self.ymmin
      self.ymmaxlocal = self.ymmax
      self.zmminlocal = self.zmmin
      self.zmmaxlocal = self.zmmax
      self.fsdecomp = Decomposition()
      self.initializeDecomposition(self.fsdecomp)
      self.fsdecomp.ix = 0
      self.fsdecomp.nx = self.nx
      self.fsdecomp.xmin = self.xmmin
      self.fsdecomp.xmax = self.xmmax
      self.fsdecomp.iy = 0
      self.fsdecomp.ny = self.ny
      self.fsdecomp.ymin = self.ymmin
      self.fsdecomp.ymax = self.ymmax
      self.fsdecomp.iz = 0
      self.fsdecomp.nz = self.nz
      self.fsdecomp.zmin = self.zmmin
      self.fsdecomp.zmax = self.zmmax
      self.nxp = self.nx
      self.nyp = self.ny
      self.nzp = self.nz
      self.xmminp = self.xmmin
      self.xmmaxp = self.xmmax
      self.ymminp = self.ymmin
      self.ymmaxp = self.ymmax
      self.zmminp = self.zmmin
      self.zmmaxp = self.zmmax
      self.ppdecomp = Decomposition()
      self.initializeDecomposition(self.ppdecomp)
      self.ppdecomp.ix = 0
      self.ppdecomp.nx = self.nx
      self.ppdecomp.xmin = self.xmmin
      self.ppdecomp.xmax = self.xmmax
      self.ppdecomp.iy = 0
      self.ppdecomp.ny = self.ny
      self.ppdecomp.ymin = self.ymmin
      self.ppdecomp.ymax = self.ymmax
      self.ppdecomp.iz = 0
      self.ppdecomp.nz = self.nz
      self.ppdecomp.zmin = self.zmmin
      self.ppdecomp.zmax = self.zmmax
    else:
      self.my_index = me
      self.nprocs = npes
      self.fsdecomp = Decomposition()
      fsdecomp = self.fsdecomp
      self.initializeDecomposition(fsdecomp)
      # --- Note that self.grid_overlap must be set by the inheriting class.
      top.grid_overlap = self.grid_overlap

      domaindecomposefields(self.nx,self.nxprocs,self.lfsautodecomp,
                            fsdecomp.ix,fsdecomp.nx,self.grid_overlap)
      domaindecomposefields(self.ny,self.nyprocs,self.lfsautodecomp,
                            fsdecomp.iy,fsdecomp.ny,self.grid_overlap)
      domaindecomposefields(self.nz,self.nzprocs,self.lfsautodecomp,
                            fsdecomp.iz,fsdecomp.nz,self.grid_overlap)

      fsdecomp.xmin[:] = self.xmmin + fsdecomp.ix*self.dx
      fsdecomp.ymin[:] = self.ymmin + fsdecomp.iy*self.dy
      fsdecomp.zmin[:] = self.zmmin + fsdecomp.iz*self.dz
      fsdecomp.xmax[:] = self.xmmin + (fsdecomp.ix + fsdecomp.nx)*self.dx
      fsdecomp.ymax[:] = self.ymmin + (fsdecomp.iy + fsdecomp.ny)*self.dy
      fsdecomp.zmax[:] = self.zmmin + (fsdecomp.iz + fsdecomp.nz)*self.dz

      self.nxlocal = fsdecomp.nx[self.ixproc]
      self.nylocal = fsdecomp.ny[self.iyproc]
      self.nzlocal = fsdecomp.nz[self.izproc]

      self.xmminlocal = fsdecomp.xmin[self.ixproc]
      self.xmmaxlocal = fsdecomp.xmax[self.ixproc]
      self.ymminlocal = fsdecomp.ymin[self.iyproc]
      self.ymmaxlocal = fsdecomp.ymax[self.iyproc]
      self.zmminlocal = fsdecomp.zmin[self.izproc]
      self.zmmaxlocal = fsdecomp.zmax[self.izproc]

      self.ppdecomp = Decomposition()
      self.initializeDecomposition(self.ppdecomp)
      # --- This should only be called after the particle decomposition
      # --- has been done.
      #self.setparticledomains()

    # --- Check the mesh consistency
    self.checkmeshconsistency(self.xmmin,self.xmmax,self.nx,self.dx,'x')
    self.checkmeshconsistency(self.ymmin,self.ymmax,self.ny,self.dy,'y')
    self.checkmeshconsistency(self.zmmin,self.zmmax,self.nz,self.dz,'z')
    self.checkmeshconsistency(self.xmminlocal,self.xmmaxlocal,self.nxlocal,self.dx,'x')
    self.checkmeshconsistency(self.ymminlocal,self.ymmaxlocal,self.nylocal,self.dy,'y')
    self.checkmeshconsistency(self.zmminlocal,self.zmmaxlocal,self.nzlocal,self.dz,'z')

    self.xsymmetryplane = 0.
    self.ysymmetryplane = 0.
    self.xmesh = self.xmmin + arange(0,self.nx+1)*self.dx
    self.ymesh = self.ymmin + arange(0,self.ny+1)*self.dy
    self.zmesh = self.zmmin + arange(0,self.nz+1)*self.dz
    self.xmeshlocal = self.xmminlocal + arange(0,self.nxlocal+1)*self.dx
    self.ymeshlocal = self.ymminlocal + arange(0,self.nylocal+1)*self.dy
    self.zmeshlocal = self.zmminlocal + arange(0,self.nzlocal+1)*self.dz

    self.ix_axis = nint(-self.xmmin/self.dx)
    self.iy_axis = nint(-self.ymmin/self.dy)
    if self.nz == 0:
      self.iz_axis = 0
    else:
      self.iz_axis = nint(-self.zmmin/self.dz)

    # --- Generate a list of the neighboring processors.
    self.neighborpes = [self.convertindextoproc(ix=self.ixproc-1),
                        self.convertindextoproc(ix=self.ixproc+1),
                        self.convertindextoproc(iy=self.iyproc-1),
                        self.convertindextoproc(iy=self.iyproc+1),
                        self.convertindextoproc(iz=self.izproc-1),
                        self.convertindextoproc(iz=self.izproc+1)]

    # --- Generate a list of unique neighbors, removing -1 and duplicates.
    self.neighborpeslist = []
    for pe in self.neighborpes:
      if pe >= 0 and pe not in self.neighborpeslist:
        self.neighborpeslist.append(pe)

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

  def initializeDecomposition(self,decomp):

    # --- First, get the number of processors along each decomposition
    # --- direction.
    self.nzprocs = self.nprocs/(self.nxprocs*self.nyprocs)
    assert (self.nxprocs*self.nyprocs*self.nzprocs == self.nprocs),\
           'nxprocs*nyprocs*nzprocs must be equal to nprocs'

    # --- Find the location of each processor along each direction
    self.izproc = int(self.my_index/(self.nxprocs*self.nyprocs))
    self.iyproc = int((self.my_index - self.izproc*(self.nxprocs*self.nyprocs))/self.nxprocs)
    self.ixproc = self.my_index - self.izproc*(self.nxprocs*self.nyprocs) - self.iyproc*self.nxprocs

    decomp.my_index = self.my_index
    decomp.nxglobal = self.nx
    decomp.nyglobal = self.ny
    decomp.nzglobal = self.nz
    decomp.iprocgrid[:] = [self.ixproc,self.iyproc,self.izproc]
    decomp.nprocgrid[:] = [self.nxprocs,self.nyprocs,self.nzprocs]
    decomp.ixproc = self.ixproc
    decomp.iyproc = self.iyproc
    decomp.izproc = self.izproc
    decomp.nxprocs = self.nxprocs
    decomp.nyprocs = self.nyprocs
    decomp.nzprocs = self.nzprocs
    decomp.gchange()
    initializedecomp(decomp)

  def checkmeshconsistency(self,min,max,nn,dd,axis):
    'Checks if the mesh quantities are consistent'
    # --- Note that the factor of 1.e-5 is a somewhat arbitrary fudge factor
    assert abs((max-min) - nn*dd) < dd*1.e-5,\
      'The grid quantities along the '+axis+' axis are inconsistent:\nmin = %e\nmax = %e\nnn = %d\ndd = %e\n'%(min,max,nn,dd)

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

   # --- Need to add conversion for reading in old dump files - yuck!

   ## --- Set new z quantities if reading in an old dump file
   #if 'nzlocal' not in self.__dict__:
   #  self.nzlocal = self.nz
   #  self.zmminlocal = self.zmmin
   #  self.zmmaxlocal = self.zmmax
   #  self.nz = self.nzfull
   #  self.zmmin = self.zmminglobal
   #  self.zmmax = self.zmmaxglobal
   #  del self.nzfull
   #  del self.zmminglobal
   #  del self.zmmaxglobal

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

    # --- Setup the MPI communicators if the arrays are to be restored
    if not self.lnorestoreonpickle:
      initializedecomp(self.fsdecomp)
      initializedecomp(self.ppdecomp)

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

  def convertindextoproc(self,ix=None,iy=None,iz=None,
                         bounds=None,isperiodic=None):
    if ix is None: ix = self.ixproc
    if iy is None: iy = self.iyproc
    if iz is None: iz = self.izproc
    nx = self.nxprocs
    ny = self.nyprocs
    nz = self.nzprocs

    if bounds is None: bounds = self.bounds
    if isperiodic is None:
      isperiodic = [periodic in bounds[0:2] and nx > 1,
                    periodic in bounds[2:4] and ny > 1,
                    periodic in bounds[4:]  and nz > 1]

    if not isperiodic[0] and (ix < 0 or ix > nx-1): return -1
    if not isperiodic[1] and (iy < 0 or iy > ny-1): return -1
    if not isperiodic[2] and (iz < 0 or iz > nz-1): return -1

    if ix < 0:    ix = nx - 1
    if ix > nx-1: ix = 0
    if iy < 0:    iy = ny - 1
    if iy > ny-1: iy = 0
    if iz < 0:    iz = nz - 1
    if iz > nz-1: iz = 0

    return ix + iy*nx + iz*nx*ny

  def setparticledomains(self):
    if not self.lparallel: return

    ppdecomp = self.ppdecomp

    ppdecomp.xmin[:] = top.ppdecomp.xmin
    ppdecomp.xmax[:] = top.ppdecomp.xmax
    ppdecomp.ymin[:] = top.ppdecomp.ymin
    ppdecomp.ymax[:] = top.ppdecomp.ymax
    ppdecomp.zmin[:] = top.ppdecomp.zmin
    ppdecomp.zmax[:] = top.ppdecomp.zmax

    domaindecomposeparticles(self.nx,self.nxprocs,self.nxpextra,
                             self.xmmin,self.xmmax,self.dx,
                             zeros(self.nxprocs,'d'),true,
                             ppdecomp.ix,ppdecomp.nx,
                             ppdecomp.xmin,ppdecomp.xmax)

    domaindecomposeparticles(self.ny,self.nyprocs,self.nypextra,
                             self.ymmin,self.ymmax,self.dy,
                             zeros(self.nyprocs,'d'),true,
                             ppdecomp.iy,ppdecomp.ny,
                             ppdecomp.ymin,ppdecomp.ymax)

    domaindecomposeparticles(self.nz,self.nzprocs,self.nzpextra,
                             self.zmmin,self.zmmax,self.dz,
                             zeros(self.nzprocs,'d'),true,
                             ppdecomp.iz,ppdecomp.nz,
                             ppdecomp.zmin,ppdecomp.zmax)

    self.nxp = ppdecomp.nx[self.ixproc]
    self.nyp = ppdecomp.ny[self.iyproc]
    self.nzp = ppdecomp.nz[self.izproc]
    self.xmminp = self.xmmin + ppdecomp.ix[self.ixproc]*self.dx
    self.xmmaxp = self.xmmin + (ppdecomp.ix[self.ixproc]+self.nxp)*self.dx
    self.ymminp = self.ymmin + ppdecomp.iy[self.iyproc]*self.dy
    self.ymmaxp = self.ymmin + (ppdecomp.iy[self.iyproc]+self.nyp)*self.dy
    self.zmminp = self.zmmin + ppdecomp.iz[self.izproc]*self.dz
    self.zmmaxp = self.zmmin + (ppdecomp.iz[self.izproc]+self.nzp)*self.dz

    self.xpminlocal = top.xpminlocal
    self.ypminlocal = top.ypminlocal
    self.zpminlocal = top.zpminlocal
    self.xpmaxlocal = top.xpmaxlocal
    self.ypmaxlocal = top.ypmaxlocal
    self.zpmaxlocal = top.zpmaxlocal

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
    """Returns the slice of fieldparray determined by the inputs. If
    fieldparray was not created, returns 0."""
    indts = min(indts,top.nsndtsphi-1)
    try:
      return self.fieldparray[...,indts,iselfb]
    except AttributeError:
      return 0.

  def returnpotentialp(self,indts,iselfb):
    """Returns the slice of potentialparray determined by the inputs. If
    potentialparray was not created, returns 0."""
    indts = min(indts,top.nsndtsphi-1)
    try:
      return self.potentialparray[...,indts,iselfb]
    except AttributeError:
      return 0.

  def returnsourcep(self,isourcepndtscopies,indts,iselfb):
    """Returns the slice of sourceparray determined by the inputs. If
    sourceparray was not created, returns 0."""
    try:
      return self.sourceparray[...,isourcepndtscopies,indts,iselfb]
    except AttributeError:
      return 0.

  def returnfield(self,indts,iselfb):
    """Returns the slice of fieldarray determined by the inputs. If
    fieldarray was not created, returns 0."""
    indts = min(indts,top.nsndtsphi-1)
    try:
      return self.fieldarray[...,indts,iselfb]
    except AttributeError:
      return 0.

  def returnpotential(self,indts,iselfb):
    """Returns the slice of potentialarray determined by the inputs. If
    potentialarray was not created, returns 0."""
    indts = min(indts,top.nsndtsphi-1)
    try:
      return self.potentialarray[...,indts,iselfb]
    except AttributeError:
      return 0.

  def returnsource(self,indts,iselfb):
    """Returns the slice of sourcearray determined by the inputs. If
    sourcearray was not created, returns 0."""
    indts = min(indts,top.nsndtsphi-1)
    try:
      return self.sourcearray[...,indts,iselfb]
    except AttributeError:
      return 0.

  def setsourcepforparticles(self,isourcepndtscopies,indts,iselfb):
    "Sets sourcep attribute to the currently active source slice"
    self.sourcep = self.returnsourcep(isourcepndtscopies,indts,iselfb)

  def setpotentialpforparticles(self,isourcepndtscopies,indts,iselfb):
    "Sets potentialp attribute to the currently active potential slice"
    self.potentialp = self.returnpotentialp(indts,iselfb)

  def setfieldpforparticles(self,isourcepndtscopies,indts,iselfb):
    "Sets fieldp attribute to the currently active field slice"
    self.fieldp = self.returnfieldp(indts,iselfb)

  def setsourceforfieldsolve(self,isourcepndtscopies,indts,iselfb):
    "Sets source attribute to the currently active source slice"
    # --- This is called at the end of loadrho just before the b.c.'s are set
    self.source = self.returnsource(indts,iselfb)

  def setpotentialforfieldsolve(self,isourcepndtscopies,indts,iselfb):
    "Sets potential attribute to the currently active potential slice"
    self.potential = self.returnpotential(indts,iselfb)

  def setfieldforfieldsolve(self,isourcepndtscopies,indts,iselfb):
    "Sets field attribute to the currently active field slice"
    self.field = self.returnfield(indts,iselfb)

  def setarraysforfieldsolve(self,isourcepndtscopies,indts,iselfb):
    "Sets source, field and potential attribute to the currently active slices"
    # --- This is called at the beginning of the field solve
    self.source    = self.returnsource(indts,iselfb)
    self.field     = self.returnfield(indts,iselfb)
    self.potential = self.returnpotential(indts,iselfb)

  def getpotentialpforparticles(self,isourcepndtscopies,indts,iselfb):
    "Copies from potential to potentialp"
    # --- In the serial case, potentialp and self.potential point to the
    # --- same memory so no copy is needed.
    if lparallel:
      # --- This is dummy code and it never executed! This should probably
      # --- be deleted to avoid confusion. getpotentialpforparticles is
      # --- redefined in the inheriting classes.
      if type(self.potential)<>type(0.):
        potentialp = self.returnpotentialp(indts,iselfb)
        potentialp[...] = self.potential
        
  def getfieldpforparticles(self,isourcepndtscopies,indts,iselfb):
    "Copies from field to fieldp"
    # --- In the serial case, fieldp and self.field point to the
    # --- same memory so no copy is needed.
    # --- Note that, so far, this method is never used.
    if lparallel:
      # --- This is dummy code and it never executed! This should probably
      # --- be deleted to avoid confusion.
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
    '''Charge deposition, uses particles from top directly
      - jslist: option list of species to load'''
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

      jslist = kw.get('jslist',None)
      if jslist is None: jslist = range(js1,js2)

      for js in jslist:
        n = pgroup.nps[js]
        if n == 0: continue
        if pgroup.ldts[js]:
          indts = top.ndtstorho[pgroup.ndts[js]-1]
          iselfb = pgroup.iselfb[js]
          self.setsourcepforparticles(0,indts,iselfb)

          if self.debug:
            i1 = pgroup.ins[js]-1
            i2 = pgroup.ins[js]+pgroup.nps[js]-1
            if self.nxlocal > 0:
              x = pgroup.xp[i1:i2]
              if self.l4symtry: x = abs(x)
              if self.solvergeom == w3d.RZgeom:
                y = pgroup.yp[i1:i2]
                x = sqrt(x**2 + y**2)
              assert x.min() >= self.xmminp,\
                     "Particles in species %d have x below the grid when depositing the source, min x = %e"%(js,x.min())
              assert x.max() < self.xmmaxp,\
                     "Particles in species %d have x above the grid when depositing the source, max x = %e"%(js,x.max())
            if self.nylocal > 0:
              y = pgroup.yp[i1:i2]
              if self.l4symtry or self.l2symtry: y = abs(y)
              assert y.min() >= self.ymminp,\
                     "Particles in species %d have y below the grid when depositing the source, min y = %e"%(js,y.min())
              assert y.max() < self.ymmaxp,\
                     "Particles in species %d have y above the grid when depositing the source, max y = %e"%(js,y.max())
            if self.nzlocal > 0:
              z = pgroup.zp[i1:i2]
              assert z.min() >= self.zmminp+self.getzgridndts()[indts],\
                     "Particles in species %d have z below the grid when depositing the source, min z = %e"%(js,z.min())
              assert z.max() < self.zmmaxp+self.getzgridndts()[indts],\
                     "Particles in species %d have z above the grid when depositing the source, max z = %e"%(js,z.max())

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
           and not sum(top.ldts))): continue
      for iselfb in range(top.nsselfb):
        isndts = min(indts,top.nsndtsphi)
        self.setsourceforfieldsolve(top.nrhopndtscopies-1,isndts,iselfb)
        self.applysourceboundaryconditions()

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
      pgroup = w3d.pgroupfsapi
      x = pgroup.xp[ipmin-1:ipmin-1+w3d.npfsapi]
      y = pgroup.yp[ipmin-1:ipmin-1+w3d.npfsapi]
      z = pgroup.zp[ipmin-1:ipmin-1+w3d.npfsapi]
      ex = pgroup.ex[ipmin-1:ipmin-1+w3d.npfsapi]
      ey = pgroup.ey[ipmin-1:ipmin-1+w3d.npfsapi]
      ez = pgroup.ez[ipmin-1:ipmin-1+w3d.npfsapi]
      bx = pgroup.bx[ipmin-1:ipmin-1+w3d.npfsapi]
      by = pgroup.by[ipmin-1:ipmin-1+w3d.npfsapi]
      bz = pgroup.bz[ipmin-1:ipmin-1+w3d.npfsapi]

    jsid = w3d.jsfsapi
    if jsid < 0: js = 0
    else:        js = jsid

    if self.debug and top.efetch[js] != 5:
      if self.nxlocal > 0:
        xdebug = x
        if self.l4symtry: xdebug = abs(x)
        if self.solvergeom == w3d.RZgeom:
          ydebug = y
          xdebug = sqrt(xdebug**2 + ydebug**2)
        assert xdebug.min() >= self.xmminp,\
               "Particles in species %d have x below the grid when fetching the field, min x = %e < %e"%(jsid,xdebug.min(),self.xmminp)
        assert xdebug.max() < self.xmmaxp,\
               "Particles in species %d have x above the grid when fetching the field, max x = %e > %e"%(jsid,xdebug.max(),self.xmmaxp)
      if self.nylocal > 0:
        ydebug = y
        if self.l4symtry or self.l2symtry: ydebug = abs(y)
        assert ydebug.min() >= self.ymminp,\
               "Particles in species %d have y below the grid when fetching the field, min y = %e < %e"%(jsid,ydebug.min(),self.ymminp)
        assert ydebug.max() < self.ymmaxp,\
               "Particles in species %d have y above the grid when fetching the field, max y = %e > %e"%(jsid,ydebug.max(),self.ymmaxp)
      if self.nzlocal > 0:
        assert z.min() >= self.zmminp+self.getzgridprv(),\
               "Particles in species %d have z below the grid when fetching the field, min z = %e < %e+%e"%(jsid,z.min(),self.zmminp,self.getzgridprv())
        assert z.max() < self.zmmaxp+self.getzgridprv(),\
               "Particles in species %d have z above the grid when fetching the field, max z = %e > %e+%e"%(jsid,z.max(),self.zmmaxp,self.getzgridprv())

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
      if self.nxlocal > 0:
        xdebug = x
        if self.l4symtry: xdebug = abs(x)
        if self.solvergeom == w3d.RZgeom:
          ydebug = y
          xdebug = sqrt(xdebug**2 + ydebug**2)
        assert xdebug.min() >= self.xpminlocal,\
               "Particles have x below the grid when fetching the potential"
        assert xdebug.max() < self.xpmaxlocal,\
               "Particles have x above the grid when fetching the potential"
      if self.nylocal > 0:
        ydebug = y
        if self.l4symtry or self.l2symtry: ydebug = abs(y)
        assert ydebug.min() >= self.ypminlocal,\
               "Particles have y below the grid when fetching the potential"
        assert ydebug.max() < self.ypmaxlocal,\
               "Particles have y above the grid when fetching the potential"
      if self.nzlocal > 0:
        assert z.min() >= self.zpminlocal,\
               "Particles have z below the grid when fetching the potential"
        assert z.max() <= self.zpmaxlocal,\
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
    raise NotImplementedError,"getpdims must be supplied - it should return a list of the dimensions of the arrays used by the particles"

  def getdims(self):
    raise NotImplementedError,"getdims must be supplied - it should return a list of the dimensions of the arrays used by the field solve"

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

    sourcepdims = list(pdims[0]) + [top.nrhopndtscopies,top.nsndts,nsselfb]
    if ('sourceparray' not in self.__dict__ or
        shape(self.sourceparray) != tuple(sourcepdims)):
      self.sourceparray = fzeros(sourcepdims,'d')

    potentialpdims = list(pdims[-1]) + [top.nsndtsphi,nsselfb]
    if ('potentialparray' not in self.__dict__ or
        shape(self.potentialparray) != tuple(potentialpdims)):
      self.potentialparray = fzeros(potentialpdims,'d')

    if len(pdims) == 3:
      # --- Also, create fieldparray
      fieldpdims = list(pdims[1]) + [top.nsndtsphi,nsselfb]
      if ('fieldparray' not in self.__dict__ or
          shape(self.fieldparray) != tuple(fieldpdims)):
        self.fieldparray = fzeros(fieldpdims,'d')

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
      if ('sourcearray' not in self.__dict__ or
          shape(self.sourcearray) != tuple(sourcedims)):
        self.sourcearray = fzeros(sourcedims,'d')

      potentialdims = list(dims[-1]) + [top.nsndtsphi,nsselfb]
      if ('potentialarray' not in self.__dict__ or
          shape(self.potentialarray) != tuple(potentialdims)):
        self.potentialarray = fzeros(potentialdims,'d')

      if len(dims) == 3:
        # --- Also, create fieldarray
        fielddims = list(dims[1]) + [top.nsndtsphi,nsselfb]
        if ('fieldarray' not in self.__dict__ or
            shape(self.fieldarray) != tuple(fielddims)):
          self.fieldarray = fzeros(fielddims,'d')

  def resetparticledomains(self):
    self.setparticledomains()
    # --- Note that this will call allocatedataarrays.
    #self.getallpotentialpforparticles(lforce=1)
    self.allocatedataarrays()
    # --- Make sure the sourcep gets set to the updated sourceparray.
    self.setsourcepforparticles(0,0,0)

  def zerosource(self):
    try:
      self.sourcearray[...] = 0.
    except AttributeError:
      pass

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
    raise NotImplementedError,"fullv subcycling not yet implemented"

  def zerosourcepwithhalfvsubcycling(self):
    raise NotImplementedError,"halfv subcycling not yet implemented"

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
          if self.solvergeom == w3d.RZgeom:
            y = pgroup.yp[i1:i2]
            x = sqrt(x**2 + y**2)
          assert abs(x-self.xmmin).min() >= 0.,\
                 text+"Particles in species %d have x below the grid when depositing the source, min x = %e"%(js,x.min())
          assert x.max() < self.xmmax,\
                 text+"Particles in species %d have x above the grid when depositing the source, min x = %e"%(js,x.max())
        if self.ny > 0:
          y = pgroup.yp[i1:i2]
          assert abs(y-self.ymmin).min() >= 0.,\
                 text+"Particles in species %d have y below the grid when depositing the source, min x = %e"%(js,y.min())
          assert y.max() < self.ymmax,\
                 text+"Particles in species %d have y above the grid when depositing the source, min x = %e"%(js,y.max())
        if self.nzlocal > 0:
          z = pgroup.zp[i1:i2]
          assert z.min() >= self.zmminp+self.getzgridndts()[indts],\
                 text+"Particles in species %d have z below the grid when depositing the source, min x = %e"%(js,z.min())
          assert z.max() < self.zmmaxp+self.getzgridndts()[indts],\
                 text+"Particles in species %d have z above the grid when depositing the source, min x = %e"%(js,z.max())


##########################################################################
##########################################################################
# --- These functions returns or sets slices of any decomposed array whose
# --- shape is the same as rho.
##########################################################################
def getdecomposedarray(arr,ix=None,iy=None,iz=None,bcast=1,local=0,
                       fullplane=0,xyantisymmetric=0,solver=None):
  """Returns slices of a decomposed array, The shape of
the object returned depends on the number of ix, iy and iz specified, which
can be from none to all three. Note that the values of ix, iy and iz are
relative to the fortran indexing, meaning that 0 is the lower boundary
of the domain.
  - ix = None:
  - iy = None: Defaults to 0 except when using 3-D geometry.
  - iz = None:
  - bcast=1: When 1, the result is broadcast to all of the processors
             (otherwise returns None to all but PE0
  - local=0: When 1, in the parallel version, each process will get its local
             value of array - no communication is done. Has no effect for
             serial version.
  - fullplane=0: When 1 and with transverse symmetries, the data is
                 replicated to fill the symmetric regions of the plane.
  - xyantisymmetric=0: When 1 and with fullplane=1, the data is treated as
                       anti-symmetric in the transverse plane, resulting
                       in a sign change during the replication.
  """
  if solver is None: solver = (getregisteredsolver() or w3d)
  if solver == w3d: decomp = top.fsdecomp
  else:             decomp = solver.fsdecomp
  if len(arr.shape) == 2: iy = None

  nx = decomp.nxglobal
  ny = decomp.nyglobal
  nz = decomp.nzglobal
  nxlocal = decomp.nx[decomp.ixproc]
  nylocal = decomp.ny[decomp.iyproc]
  nzlocal = decomp.nz[decomp.izproc]

  if (decomp.nxprocs <= 1 and decomp.nyprocs <=1 and decomp.nzprocs <= 1):
    local = 1

  if local or not lparallel:
    if ix is None     and iy is None     and iz is None    :
      result = arr[...]
    if ix is not None and iy is None     and iz is None    :
      result = arr[ix,...]
    if ix is None     and iy is not None and iz is None    :
      result = arr[:,iy,:]
    if ix is None     and iy is None     and iz is not None:
      result = arr[...,iz]
    if ix is not None and iy is not None and iz is None    :
      result = arr[ix,iy,:]
    if ix is not None and iy is None     and iz is not None:
      result = arr[ix,...,iz]
    if ix is None     and iy is not None and iz is not None:
      result = arr[:,iy,iz]
    if ix is not None and iy is not None and iz is not None:
      result = arr[ix,iy,iz]
  else:

    # --- Get the local extent of each processor.
    my_ixpp = decomp.ix[decomp.ixproc]
    my_nxpp = decomp.nx[decomp.ixproc]
    my_iypp = decomp.iy[decomp.iyproc]
    my_nypp = decomp.ny[decomp.iyproc]
    my_izpp = decomp.iz[decomp.izproc]
    my_nzpp = decomp.nz[decomp.izproc]

    # --- If ix,iy or iz was given, check if it is within the local domain.
    if ((ix is None or my_ixpp <= ix and ix <= my_ixpp+my_nxpp) and
        (iy is None or my_iypp <= iy and iy <= my_iypp+my_nypp) and
        (iz is None or my_izpp <= iz and iz <= my_izpp+my_nzpp)):
      # --- If so, grab the appropriate slice of array.
      sss = [slice(1+nxlocal),
             slice(1+nylocal),
             slice(1+nzlocal)]
      if ix is not None: sss[0] = slice(ix-my_ixpp,ix-my_ixpp+1)
      if iy is not None: sss[1] = slice(iy-my_iypp,iy-my_iypp+1)
      if iz is not None: sss[2] = slice(iz-my_izpp,iz-my_izpp+1)
      if nx == 0: sss[0] = Ellipsis
      if ny == 0: sss[1] = Ellipsis
      if nz == 0: sss[2] = Ellipsis
      result = arr[sss[0],sss[1],sss[2]]
    else:
      # --- Otherwise, use None
      result = None

    # --- Get the data (or None) from all of the processors.
    resultlist = gather(result)

    if me == 0:
      # --- Setup the size of the array to be returned and create it.
      sss = [1+nx,1+ny,1+nz]
      if ix is not None: sss[0] = 1
      if iy is not None: sss[1] = 1
      if iz is not None: sss[2] = 1
      if nz == 0: del sss[2]
      if ny == 0: del sss[1]
      if nx == 0: del sss[0]
      resultglobal = fzeros(sss,'d')

      # --- Loop over all processors and grab the data sent, putting it into
      # --- the appropriate place in the array.
      iproc = 0
      ix1,ix2 = 0,1
      iy1,iy2 = 0,1
      iz1,iz2 = 0,1
      sss = [1,1,1]
      for izproc in range(decomp.nzprocs):
        for iyproc in range(decomp.nyprocs):
          for ixproc in range(decomp.nxprocs):
            if resultlist[iproc] is not None:
              if ix is None:
                ix1 = decomp.ix[ixproc]
                ix2 = decomp.ix[ixproc] + decomp.nx[ixproc] + 1
              if iy is None:
                iy1 = decomp.iy[iyproc]
                iy2 = decomp.iy[iyproc] + decomp.ny[iyproc] + 1
              if iz is None:
                iz1 = decomp.iz[izproc]
                iz2 = decomp.iz[izproc] + decomp.nz[izproc] + 1
              sss[0] = slice(ix1,ix2)
              sss[1] = slice(iy1,iy2)
              sss[2] = slice(iz1,iz2)
              if nx == 0: sss[0] = Ellipsis
              if ny == 0: sss[1] = Ellipsis
              if nz == 0: sss[2] = Ellipsis
              resultglobal[sss[0],sss[1],sss[2]] = resultlist[iproc]
            iproc += 1

      # --- Now remove any of the reduced dimensions.
      if ix is None: ix = slice(None)
      else:          ix = 0
      if iy is None: iy = slice(None)
      else:          iy = 0
      if iz is None: iz = slice(None)
      else:          iz = 0
      if nx == 0: ix = Ellipsis
      if ny == 0: iy = Ellipsis
      if nz == 0: iz = Ellipsis
      result = resultglobal[ix,iy,iz]

    if bcast:
      result = parallel.broadcast(result)
    else:
      if me > 0: return None

  if not fullplane:
    return result
  else:
    ii = 0
    if ix is None and (solver.l4symtry or solver.solvergeom == w3d.RZgeom):
      ss = array(shape(result))
      nn = ss[ii] - 1
      ss[ii] = 2*nn + 1
      result1 = empty(tuple(ss),'d')
      if xyantisymmetric: fsign = -1
      else:               fsign = +1
      result1[nn:,...] = result
      result1[nn::-1,...] = fsign*result
      result = result1
    if ix is None: ii = ii + 1
    if iy is None and (solver.l2symtry or solver.l4symtry):
      ss = array(shape(result))
      nn = ss[ii] - 1
      ss[ii] = 2*nn + 1
      result1 = empty(tuple(ss),'d')
      if xyantisymmetric: fsign = -1
      else:               fsign = +1
      if ii == 0:
        result1[nn:,...] = result
        result1[nn::-1,...] = fsign*result
      else:
        result1[:,nn:,...] = result
        result1[:,nn::-1,...] = fsign*result
      result = result1
    return result

# --------------------------------------------------------------------------
def setdecomposedarray(arr,val,ix=None,iy=None,iz=None,local=0,solver=None):
  """Sets slices of a decomposed array. The shape of
the input object depends on the number of arguments specified, which can
be from none to all three.
  - val: input array (must be supplied)
  - ix = None:
  - iy = None: Defaults to 0 except when using 3-D geometry.
  - iz = None:
  """
  if solver is None: solver = (getregisteredsolver() or w3d)
  if solver == w3d: decomp = top.fsdecomp
  else:             decomp = solver.fsdecomp
  if len(arr.shape) == 2: iy = None

  nx = decomp.nxglobal
  ny = decomp.nyglobal
  nz = decomp.nzglobal
  nxlocal = decomp.nx[decomp.ixproc]
  nylocal = decomp.ny[decomp.iyproc]
  nzlocal = decomp.nz[decomp.izproc]

  if (decomp.nxprocs <= 1 and decomp.nyprocs <=1 and decomp.nzprocs <= 1):
    local = 1

  if local or not lparallel:
    if ix is None     and iy is None     and iz is None    :
      arr[...] = val
    if ix is not None and iy is None     and iz is None    :
      arr[ix,...] = val
    if ix is None     and iy is not None and iz is None    :
      arr[:,iy,:] = val
    if ix is None     and iy is None     and iz is not None:
      arr[...,iz] = val
    if ix is not None and iy is not None and iz is None    :
      arr[ix,iy,:] = val
    if ix is not None and iy is None     and iz is not None:
      arr[ix,:,iz] = val
    if ix is None     and iy is not None and iz is not None:
      arr[:,iy,iz] = val
    if ix is not None and iy is not None and iz is not None:
      arr[ix,iy,iz] = val
  else:

    ppplist = []
    if me == 0:

      # --- Add extra dimensions so that the input has the same number of
      # --- dimensions as array.
      ppp = array(val)
      sss = list(ppp.shape)
      if ix is not None and nx > 0: sss[0:0] = [1]
      if iy is not None and ny > 0: sss[1:1] = [1]
      if iz is not None and nz > 0: sss[2:2] = [1]
      ppp.shape = sss

      # --- Loop over all processors and grab the chunk of the input that
      # --- overlaps each of the domains.
      ix1,ix2 = 0,1
      iy1,iy2 = 0,1
      iz1,iz2 = 0,1
      sss = [1,1,1]
      for izproc in range(decomp.nzprocs):
        for iyproc in range(decomp.nyprocs):
          for ixproc in range(decomp.nxprocs):
            if ix is None:
              ix1 = decomp.ix[ixproc]
              ix2 = decomp.ix[ixproc] + decomp.nx[ixproc] + 1
            if iy is None:
              iy1 = decomp.iy[iyproc]
              iy2 = decomp.iy[iyproc] + decomp.ny[iyproc] + 1
            if iz is None:
              iz1 = decomp.iz[izproc]
              iz2 = decomp.iz[izproc] + decomp.nz[izproc] + 1
            sss[0] = slice(ix1,ix2)
            sss[1] = slice(iy1,iy2)
            sss[2] = slice(iz1,iz2)
            if nx == 0: sss[0] = Ellipsis
            if ny == 0: sss[1] = Ellipsis
            if nz == 0: sss[2] = Ellipsis
            ppplist.append(ppp[sss[0],sss[1],sss[2]])

    # --- Send the data to each of the processors
    ppp = comm_world.scatter(ppplist)[0]

    # --- Get the local extent of each processor.
    my_ixpp = decomp.ix[decomp.ixproc]
    my_nxpp = decomp.nx[decomp.ixproc]
    my_iypp = decomp.iy[decomp.iyproc]
    my_nypp = decomp.ny[decomp.iyproc]
    my_izpp = decomp.iz[decomp.izproc]
    my_nzpp = decomp.nz[decomp.izproc]

    # --- If ix,iy or iz was given, check if it is within the local domain.
    if ((ix is None or my_ixpp <= ix and ix <= my_ixpp+my_nxpp) and
        (iy is None or my_iypp <= iy and iy <= my_iypp+my_nypp) and
        (iz is None or my_izpp <= iz and iz <= my_izpp+my_nzpp)):
      # --- If so, set the appropriate slice of array.
      sss = [slice(1+nxlocal),
             slice(1+nylocal),
             slice(1+nzlocal)]
      if ix is not None: sss[0] = slice(ix-my_ixpp,ix-my_ixpp+1)
      if iy is not None: sss[1] = slice(iy-my_iypp,iy-my_iypp+1)
      if iz is not None: sss[2] = slice(iz-my_izpp,iz-my_izpp+1)
      if nx == 0: sss[0] = Ellipsis
      if ny == 0: sss[1] = Ellipsis
      if nz == 0: sss[2] = Ellipsis
      arr[sss[0],sss[1],sss[2]] = ppp

# --------------------------------------------------------------------------
def getrho(ix=None,iy=None,iz=None,bcast=1,local=0,fullplane=0,solver=None):
  """Returns slices of rho, the charge density array. The shape of the object
returned depends on the number of ix, iy and iz specified, which can be from
none to all three. If no components are given, then a 3-D array is returned. 
With one, a 2-D array, with two, 1-D array and with 3 a scalar is returned.
Note that 0 is the lower edge of the domain and nx, ny or nz is the upper edge.
  - ix = None:
  - iy = None: Defaults to 0 except when using 3-D geometry.
  - iz = None:
  - bcast=1: When 1, the result is broadcast to all of the processors
             (otherwise returns None to all but PE0
  - local=0: When 1, in the parallel version, each process will get its local
             value of rho - no communication is done. Has no effect for serial
             version.
  - fullplane=0: When 1 and with transverse symmetries, the data is
                 replicated to fill the symmetric regions of the plane.
  """
  if solver is None: solver = (getregisteredsolver() or w3d)
  if solver is w3d:
    if solver.solvergeom in [w3d.RZgeom,w3d.XZgeom,w3d.Zgeom]:
      rho = frz.basegrid.rho
    else:
      rho = w3d.rho
  else:
    rho = solver.getrho()

  return getdecomposedarray(rho,ix=ix,iy=iy,iz=iz,bcast=bcast,local=local,
                            fullplane=fullplane,solver=solver)

# --------------------------------------------------------------------------
def setrho(val,ix=None,iy=None,iz=None,local=0,solver=None):
  """Sets slices of rho, the charge density array. The shape of the object
returned depends on the number of ix, iy and iz specified, which can be from
none to all three. If no components are given, then a 3-D array is returned. 
With one, a 2-D array, with two, 1-D array and with 3 a scalar is returned.
Note that 0 is the lower edge of the domain and nx, ny or nz is the upper edge.
  - val: input array (must be supplied)
  - ix = None:
  - iy = None: Defaults to 0 except when using 3-D geometry.
  - iz = None:
  """
  if solver is None: solver = (getregisteredsolver() or w3d)
  if solver is w3d:
    if solver.solvergeom in [w3d.RZgeom,w3d.XZgeom,w3d.Zgeom]:
      rho = frz.basegrid.rho
    else:
      rho = w3d.rho
  else:
    rho = solver.getrho()

  setdecomposedarray(rho,val,ix=ix,iy=iy,iz=iz,local=local,solver=solver)

# --------------------------------------------------------------------------
def getphi(ix=None,iy=None,iz=None,bcast=1,local=0,fullplane=0,solver=None):
  """Returns slices of phi, the electrostatic potential array. The shape of
the object returned depends on the number of ix, iy and iz specified, which
can be from none to all three. If no components are given, then a 3-D array
is returned.  With one, a 2-D array, with two, 1-D array and with 3 a scalar
is returned.  Note that 0 is the lower edge of the domain and nx, ny or nz is
the upper edge.  
  - ix = None:
  - iy = None: Defaults to 0 except when using 3-D geometry.
  - iz = None:
               from -1 to nz+1
  - bcast=1: When 1, the result is broadcast to all of the processors
             (otherwise returns None to all but PE0
  - local=0: When 1, in the parallel version, each process will get its local
             value of phi - no communication is done. Has no effect for serial
             version.
  - fullplane=0: When 1 and with transverse symmetries, the data is
                 replicated to fill the symmetric regions of the plane.
  """
  if solver is None: solver = (getregisteredsolver() or w3d)
  if solver is w3d:
    if solver.solvergeom in [w3d.RZgeom,w3d.XZgeom,w3d.Zgeom]:
      phi = frz.basegrid.phi[1:-1,1:-1]
      iy = None
    else:
      phi = w3d.phi[1:-1,1:-1,1:-1]
  else:
    phi = solver.getphi()

  return getdecomposedarray(phi,ix=ix,iy=iy,iz=iz,bcast=bcast,local=local,
                            fullplane=fullplane,solver=solver)

# --------------------------------------------------------------------------
def setphi(val,ix=None,iy=None,iz=None,local=0,solver=None):
  """Sets slices of phi, the electrostatic potential array. The shape of
the input object depends on the number of arguments specified, which can
be from none to all three.
  The shape of the object
returned depends on the number of ix, iy and iz specified, which can be from
none to all three. If no components are given, then a 3-D array is returned. 
With one, a 2-D array, with two, 1-D array and with 3 a scalar is returned.
Note that 0 is the lower edge of the domain and nx, ny or nz is the upper edge.
  - val: input array (must be supplied)
  - ix = None:
  - iy = None: Defaults to 0 except when using 3-D geometry.
  - iz = None: Value is relative to the fortran indexing, so iz ranges
               from -1 to nz+1
  """
  if solver is None: solver = (getregisteredsolver() or w3d)
  if solver is w3d:
    if solver.solvergeom in [w3d.RZgeom,w3d.XZgeom,w3d.Zgeom]:
      phi = frz.basegrid.phi[1:-1,1:-1]
      iy = None
    else:
      phi = w3d.phi[1:-1,1:-1,1:-1]
  else:
    phi = solver.getphi()

  setdecomposedarray(phi,val,ix=ix,iy=iy,iz=iz,local=local,solver=solver)

# --------------------------------------------------------------------------
def getselfe(comp=None,ix=None,iy=None,iz=None,bcast=1,local=0,fullplane=0,
             solver=None):
  """Returns slices of selfe, the electrostatic field array. The shape of the
object returned depends on the number of ix, iy and iz specified, which can
be from none to all three. If no components are given, then a 3-D array is
returned.  With one, a 2-D array, with two, 1-D array and with 3 a scalar is
returned.  Note that 0 is the lower edge of the domain and nx, ny or nz is
the upper edge.  
  - comp: field component to get, either 'x', 'y', 'z', or 'E', must be given.
          Use 'E' to get the field magnitude.
  - ix = None:
  - iy = None: Defaults to 0 except when using 3-D geometry.
  - iz = None:
  - bcast=1: When 1, the result is broadcast to all of the processors
             (otherwise returns None to all but PE0
  - local=0: When 1, in the parallel version, each process will get its local
             value of E - no communication is done. Has no effect for serial
             version.
  - fullplane=0: When 1 and with transverse symmetries, the data is
                 replicated to fill the symmetric regions of the plane.
  """
  assert comp in ['x','y','z','E'],"comp must be one of 'x', 'y', 'z' or 'E'"
  if solver is None: solver = (getregisteredsolver() or w3d)
  if iy is None and solver.solvergeom in [w3d.RZgeom,w3d.XZgeom,w3d.Zgeom]: iy=0
  if alltrue(top.efetch != 3) or not w3d.allocated('selfe'):
    # --- If not already using selfe, then allocate it and set it.
    # --- Note that this could be an unexpected expense for a user.
    if solver is w3d:
      allocateselfeforfieldsolve()
      nx,ny,nz = array(w3d.phi.shape) - 1
      getselfe3d(w3d.phi,w3d.nxlocal,w3d.nylocal,w3d.nzlocal,
                 w3d.selfe,w3d.dx,w3d.dy,w3d.dz,
                 true,(nx-w3d.nxlocal)/2,(ny-w3d.nylocal)/2,(nz-w3d.nzlocal)/2)
      selfe = w3d.selfe
    else:
      selfe = solver.getselfe()
  if type(comp) == IntType: ic = comp
  else:                     ic = ['x','y','z','E'].index(comp)

  if comp == 'E':
    Ex = getdecomposedarray(selfe[0,...],ix=ix,iy=iy,iz=iz,
                            bcast=bcast,local=local,fullplane=fullplane,
                            xyantisymmetric=(ic in [0,1]),
                            solver=solver)
    Ey = getdecomposedarray(selfe[1,...],ix=ix,iy=iy,iz=iz,
                            bcast=bcast,local=local,fullplane=fullplane,
                            xyantisymmetric=(ic in [0,1]),
                            solver=solver)
    Ez = getdecomposedarray(selfe[2,...],ix=ix,iy=iy,iz=iz,
                            bcast=bcast,local=local,fullplane=fullplane,
                            xyantisymmetric=(ic in [0,1]),
                            solver=solver)
    return sqrt(Ex**2 + Ey**2 + Ez**2)
  else:
    return getdecomposedarray(selfe[ic,...],ix=ix,iy=iy,iz=iz,
                              bcast=bcast,local=local,fullplane=fullplane,
                              xyantisymmetric=(ic in [0,1]),
                              solver=solver)

# --------------------------------------------------------------------------
def getj(comp=None,ix=None,iy=None,iz=None,bcast=1,local=0,fullplane=0,
         solver=None):
  """Returns slices of J, the current density array. The shape of the object
returned depends on the number of ix, iy and iz specified, which can be from
none to all three. If no components are given, then a 3-D array is returned. 
With one, a 2-D array, with two, 1-D array and with 3 a scalar is returned.
Note that 0 is the lower edge of the domain and nx, ny or nz is the upper edge.
  - comp: field component to get, either 'x', 'y', 'z' or. 'J', must be given.
          Use 'J' to get the magnitude.
  - ix = None:
  - iy = None: Defaults to 0 except when using 3-D geometry.
  - iz = None:
  - bcast=1: When 1, the result is broadcast to all of the processors
             (otherwise returns None to all but PE0
  - local=0: When 1, in the parallel version, each process will get its local
             value of j - no communication is done. Has no effect for serial
             version.
  - fullplane=0: When 1 and with transverse symmetries, the data is
                 replicated to fill the symmetric regions of the plane.
  """
  assert comp in ['x','y','z','J'],"comp must be one of 'x', 'y', 'z' or 'J'"
  if type(comp) == IntType: ic = comp
  else:                     ic = ['x','y','z','J'].index(comp)
  if solver is None: solver = (getregisteredsolver() or w3d)
  if solver == w3d: bfield = f3d.bfield
  else:             bfield = solver

  if comp == 'J':
    Jx = getdecomposedarray(bfield.j[0,...],ix=ix,iy=iy,iz=iz,
                            bcast=bcast,local=local,fullplane=fullplane,
                            xyantisymmetric=(ic in [0,1]),
                            solver=solver)
    Jy = getdecomposedarray(bfield.j[1,...],ix=ix,iy=iy,iz=iz,
                            bcast=bcast,local=local,fullplane=fullplane,
                            xyantisymmetric=(ic in [0,1]),
                            solver=solver)
    Jz = getdecomposedarray(bfield.j[2,...],ix=ix,iy=iy,iz=iz,
                            bcast=bcast,local=local,fullplane=fullplane,
                            xyantisymmetric=(ic in [0,1]),
                            solver=solver)
    return sqrt(Jx**2 + Jy**2 + Jz**2)
  else:
    return getdecomposedarray(bfield.j[ic,...],ix=ix,iy=iy,iz=iz,
                              bcast=bcast,local=local,fullplane=fullplane,
                              xyantisymmetric=(ic in [0,1]),
                              solver=solver)

# --------------------------------------------------------------------------
def getb(comp=None,ix=None,iy=None,iz=None,bcast=1,local=0,fullplane=0,
         solver=None):
  """Returns slices of B, the magnetic field array. The shape of the object
returned depends on the number of ix, iy and iz specified, which can be from
none to all three. If no components are given, then a 3-D array is returned. 
With one, a 2-D array, with two, 1-D array and with 3 a scalar is returned.
Note that 0 is the lower edge of the domain and nx, ny or nz is the upper edge.
  - comp: field component to get, either 'x', 'y', 'z' or 'B', must be given.
          Use 'B' to get the magnitude.
  - ix = None:
  - iy = None: Defaults to 0 except when using 3-D geometry.
  - iz = None:
  - bcast=1: When 1, the result is broadcast to all of the processors
             (otherwise returns None to all but PE0
  - local=0: When 1, in the parallel version, each process will get its local
             value of j - no communication is done. Has no effect for serial
             version.
  - fullplane=0: When 1 and with transverse symmetries, the data is
                 replicated to fill the symmetric regions of the plane.
  """
  assert comp in ['x','y','z','B'],"comp must be one of 'x', 'y', 'z' or 'B'"
  if type(comp) == IntType: ic = comp
  else:                     ic = ['x','y','z','B'].index(comp)
  if solver is None: solver = (getregisteredsolver() or w3d)
  if solver == w3d: bfield = f3d.bfield
  else:             bfield = solver

  if comp == 'B':
    Bx = getdecomposedarray(bfield.b[0,...],ix=ix,iy=iy,iz=iz,
                            bcast=bcast,local=local,fullplane=fullplane,
                            xyantisymmetric=(ic in [0,1]),
                            solver=solver)
    By = getdecomposedarray(bfield.b[1,...],ix=ix,iy=iy,iz=iz,
                            bcast=bcast,local=local,fullplane=fullplane,
                            xyantisymmetric=(ic in [0,1]),
                            solver=solver)
    Bz = getdecomposedarray(bfield.b[2,...],ix=ix,iy=iy,iz=iz,
                            bcast=bcast,local=local,fullplane=fullplane,
                            xyantisymmetric=(ic in [0,1]),
                            solver=solver)
    return sqrt(Bx**2 + By**2 + Bz**2)
  else:
    return getdecomposedarray(bfield.b[ic,...],ix=ix,iy=iy,iz=iz,
                              bcast=bcast,local=local,fullplane=fullplane,
                              xyantisymmetric=(ic in [0,1]),
                              solver=solver)

# --------------------------------------------------------------------------
def geta(comp=None,ix=None,iy=None,iz=None,bcast=1,local=0,fullplane=0,
         solver=None):
  """Returns slices of B, the magnetic vector potential array. The shape of
the object returned depends on the number of ix, iy and iz specified, which
can be from none to all three. If no components are given, then a 3-D array
is returned.  With one, a 2-D array, with two, 1-D array and with 3 a scalar
is returned.  Note that 0 is the lower edge of the domain and nx, ny or nz is
the upper edge.  
  - comp: field component to get, either 'x', 'y', 'z' or 'A', must be given.
          Use 'A' to get the magnitude.
  - ix = None:
  - iy = None: Defaults to 0 except when using 3-D geometry.
  - iz = None:
  - bcast=1: When 1, the result is broadcast to all of the processors
             (otherwise returns None to all but PE0
  - local=0: When 1, in the parallel version, each process will get its local
             value of a - no communication is done. Has no effect for serial
             version.
  - fullplane=0: When 1 and with transverse symmetries, the data is
                 replicated to fill the symmetric regions of the plane.
  """
  assert comp in ['x','y','z','A'],"comp must be one of 'x', 'y', 'z' or 'A'"
  if type(comp) == IntType: ic = comp
  else:                     ic = ['x','y','z','A'].index(comp)
  if solver is None: solver = (getregisteredsolver() or w3d)
  if solver == w3d: bfield = f3d.bfield
  else:             bfield = solver

  if comp == 'A':
    Ax = getdecomposedarray(bfield.a[0,1:-1,1:-1,1:-1],ix=ix,iy=iy,iz=iz,
                            bcast=bcast,local=local,fullplane=fullplane,
                            xyantisymmetric=(ic in [0,1]),
                            solver=solver)
    Ay = getdecomposedarray(bfield.a[1,1:-1,1:-1,1:-1],ix=ix,iy=iy,iz=iz,
                            bcast=bcast,local=local,fullplane=fullplane,
                            xyantisymmetric=(ic in [0,1]),
                            solver=solver)
    Az = getdecomposedarray(bfield.a[2,1:-1,1:-1,1:-1],ix=ix,iy=iy,iz=iz,
                            bcast=bcast,local=local,fullplane=fullplane,
                            xyantisymmetric=(ic in [0,1]),
                            solver=solver)
    return sqrt(Ax**2 + Ay**2 + Az**2)
  else:
    return getdecomposedarray(bfield.a[ic,1:-1,1:-1,1:-1],ix=ix,iy=iy,iz=iz,
                              bcast=bcast,local=local,fullplane=fullplane,
                              xyantisymmetric=(ic in [0,1]),
                              solver=solver)

# --------------------------------------------------------------------------
def seta(val,comp=None,ix=None,iy=None,iz=None,local=0,solver=None):
  """Sets slices of a, the electrostatic potential array. The shape of the
object returned depends on the number of ix, iy and iz specified, which can
be from none to all three. If no components are given, then a 3-D array is
returned.  With one, a 2-D array, with two, 1-D array and with 3 a scalar is
returned.  Note that 0 is the lower edge of the domain and nx, ny or nz is
the upper edge.  
  - val: input array (must be supplied)
  - comp: field component to get, either 'x', 'y', or 'z', must be given
  - ix = None:
  - iy = None: Defaults to 0 except when using 3-D geometry.
  - iz = None:
  """
  assert comp in ['x','y','z'],"comp must be one of 'x', 'y', or 'z'"
  if type(comp) == IntType: ic = comp
  else:                     ic = ['x','y','z'].index(comp)
  if solver is None: solver = (getregisteredsolver() or w3d)
  if solver == w3d: bfield = f3d.bfield
  else:             bfield = solver

  setdecomposedarray(bfield.a[ic,1:-1,1:-1,1:-1],val,ix=ix,iy=iy,iz=iz,
                     local=local,solver=solver)

