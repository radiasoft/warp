"""Class for doing complete magnetostatic multigrid field solve"""
# ToDo:
#  - modify setj to check if particles are within grid
from warp import *
from fieldsolver import SubcycledPoissonSolver
from generateconductors import installconductors
from find_mgparam import find_mgparam
import MA

try:
  import psyco
except ImportError:
  pass

##############################################################################
class MagnetostaticMG(SubcycledPoissonSolver):
  
  __bfieldinputs__ = ['mgparam','downpasses','uppasses',
                      'mgmaxiters','mgtol','mgmaxlevels','mgform',
                      'lcndbndy','icndbndy','laddconductor',
                      'lcylindrical','lanalyticbtheta'] 
  __f3dinputs__ = ['gridmode','mgparam','downpasses','uppasses',
                   'mgmaxiters','mgtol','mgmaxlevels','mgform','mgverbose',
                   'lcndbndy','icndbndy','laddconductor'] 

  def __init__(self,**kw):
    SubcycledPoissonSolver.__init__(self,kwdict=kw)
    self.ncomponents = 3
    self.nxguard = 1
    self.nyguard = 1
    self.nzguard = 1
    self.lusevectorpotential = true

    # --- Kludge - make sure that the multigrid3df routines never sets up
    # --- any conductors.
    f3d.gridmode = 1

    # --- Save input parameters
    for name in MagnetostaticMG.__f3dinputs__:
      if name not in self.__dict__:
        #self.__dict__[name] = kw.pop(name,getattr(f3d,name)) # Python2.3
        self.__dict__[name] = kw.get(name,getattr(f3d,name))
      if kw.has_key(name): del kw[name]
    for name in MagnetostaticMG.__bfieldinputs__:
      if name not in self.__dict__:
        #self.__dict__[name] = kw.pop(name,getattr(f3d.bfield,name))#Python2.3
        self.__dict__[name] = kw.get(name,getattr(f3d.bfield,name))
      if kw.has_key(name): del kw[name]

    # --- If there are any remaning keyword arguments, raise an error.
    assert len(kw.keys()) == 0,"Bad keyword arguemnts %s"%kw.keys()

    # --- Check for cylindrical geometry
    self.lcylindrical = (self.solvergeom == w3d.RZgeom) or self.lcylindrical

    # --- Fix the values of xmmin and ymmin as necessary
    if self.lcylindrical:
      self.ny = 0
      if self.xmmin  < 0.: self.xmmin  = 0.
      if self.ymmin  < 0.: self.ymmin  = 0.
      if self.xmminp < 0.: self.xmminp = 0.
      if self.ymminp < 0.: self.ymminp = 0.
      self.l2symtry = false
      self.l4symtry = false
      self.bounds[0] = neumann
      self.bounds[2] = dirichlet
      self.bounds[3] = dirichlet

    # --- Create a conductor object, which by default is empty.
    self.conductors = ConductorType()
    self.conductorlist = []

    # --- Give these variables dummy initial values.
    self.mgiters = zeros(3)
    self.mgerror = zeros(3,'d')

    # --- Make sure that these are arrays
    self.mgmaxiters = ones(3)*self.mgmaxiters
    self.mgmaxlevels = ones(3)*self.mgmaxlevels
    self.mgparam = ones(3)*self.mgparam
    self.mgform = ones(3)*self.mgform
    self.mgtol = ones(3)*self.mgtol
    self.mgverbose = ones(3)*self.mgverbose
    self.downpasses = ones(3)*self.downpasses
    self.uppasses = ones(3)*self.uppasses

    # --- At the start, assume that there are no bends. This is corrected
    # --- in the solve method when there are bends.
    self.linbend = false

    # --- Turn of build quads option
    self.lbuildquads = false

  def __getstate__(self):
    dict = SubcycledPoissonSolver.__getstate__(self)
    if self.lreducedpickle:
      del dict['conductors']
    return dict

  def __setstate__(self,dict):
    SubcycledPoissonSolver.__setstate__(self,dict)
    if self.lreducedpickle and not self.lnorestoreonpickle:
      # --- Regenerate the conductor data
      self.conductors = ConductorType()
      conductorlist = self.conductorlist
      self.conductorlist = []
      for conductor in conductorlist:
        self.installconductor(conductor)

  def getpdims(self):
    # --- Returns the dimensions of the jp, bp, and ap arrays
    return ((3,1+self.nxp,1+self.nyp,1+self.nzp),
            (3,1+self.nxp,1+self.nyp,1+self.nzp),
            (3,3+self.nxp,3+self.nyp,3+self.nzp))

  def getdims(self):
    # --- Returns the dimensions of the j, b, and a arrays
    return ((3,1+self.nx,1+self.ny,1+self.nz),
            (3,1+self.nx,1+self.ny,1+self.nz),
            (3,3+self.nx,3+self.ny,3+self.nz))

  def loadj(self,lzero=None,**kw):
    SubcycledPoissonSolver.loadsource(self,lzero,**kw)

  def fetchb(self,*args):
    SubcycledPoissonSolver.fetchfield(self,*args)

  def setsourcep(self,js,pgroup,zgrid):
    n  = pgroup.nps[js]
    if n == 0: return
    i  = pgroup.ins[js] - 1
    x  = pgroup.xp[i:i+n]
    y  = pgroup.yp[i:i+n]
    z  = pgroup.zp[i:i+n]
    ux = pgroup.uxp[i:i+n]
    uy = pgroup.uyp[i:i+n]
    uz = pgroup.uzp[i:i+n]
    gaminv = pgroup.gaminv[i:i+n]
    q  = pgroup.sq[js]
    w  = pgroup.sw[js]*top.pgroup.dtscale[js]
    if top.wpid > 0: wght = top.pgroup[i:i+n,top.wpid-1]
    else:            wght = zeros((0,),'d')
    self.setsourcepatposition(x,y,z,ux,uy,uz,gaminv,wght,q,w,zgrid)

  def setsourcepatposition(self,x,y,z,ux,uy,uz,gaminv,wght,q,w,zgrid):
    n = len(x)
    if n == 0: return
    if len(wght) > 0:
      nw = len(wght)
    else:
      nw = 0.
      wght = zeros(1,'d')
    setj3d(self.sourcep,self.sourcep,n,x,y,z,zgrid,ux,uy,uz,gaminv,
           q,w,nw,wght,top.depos,
           self.nxp,self.nyp,self.nzp,self.dx,self.dy,self.dz,
           self.xmminp,self.ymminp,self.zmminp,
           self.l2symtry,self.l4symtry,self.lcylindrical)

  def fetchfieldfrompositions(self,x,y,z,ex,ey,ez,bx,by,bz,pgroup=None):
    n = len(x)
    if n == 0: return
    setb3d(self.fieldp,n,x,y,z,top.zgridprv,bx,by,bz,
           self.nxp,self.nyp,self.nzp,self.dx,self.dy,self.dz,
           self.xmminp,self.ymminp,self.zmminp,
           self.l2symtry,self.l4symtry,self.lcylindrical)

  def fetchpotentialfrompositions(self,x,y,z,a):
    n = len(x)
    if n == 0: return
    fetchafrompositions3d(self.potentialp,n,x,y,z,a,top.zgrid,
                          self.nxp,self.nyp,self.nzp,self.dx,self.dy,self.dz,
                          self.xmminp,self.ymminp,self.zmminp,
                          self.l2symtry,self.l4symtry,self.lcylindrical)

  def setarraysforfieldsolve(self,*args):
    if self.nslaves <= 1:
      SubcycledPoissonSolver.setarraysforfieldsolve(self,*args)
    else:
      # --- This needs checking.
      sourcedims,fielddims,potentialdims = self.getdims()
      if 'source' not in self.__dict__ or shape(self.source) != tuple(sourcedims):
        self.source = fzeros(sourcedims,'d')
      if 'field' not in self.__dict__ or shape(self.field) != tuple(fielddims):
        self.field = fzeros(fielddims,'d')
      if 'potential' not in self.__dict__ or shape(self.potential) != tuple(potentialdims):
        self.potential = fzeros(potentialdims,'d')
      SubcycledPoissonSolver.setsourcepforparticles(self,*args)
     #setrhoforfieldsolve3d(self.nx,self.ny,self.nz,self.source,
     #                      self.nxp,self.nyp,self.nzp,self.sourcep,self.nzpguard,
     #                      self.my_index,self.nslaves,self.izpslave,self.nzpslave,
     #                      self.izfsslave,self.nzfsslave)
      setjforfieldsolve3d(self.nx,self.ny,self.nz,self.source,
                          self.nxp,self.nyp,self.nzp,self.sourcep,self.nzpguard,
                          self.my_index,self.nslaves,self.izpslave,self.nzpslave,
                          self.izfsslave,self.nzfsslave)

  def getpotentialpforparticles(self,*args):
    """Despite the name, this actually gets the field instead, since that is
       always used in the magnetostatic solver"""
    if self.nslaves <= 1:
      SubcycledPoissonSolver.getfieldpforparticles(self,*args)
    else:
      self.setfieldpforparticles(*args)
      getphipforparticles3d(3,self.nx,self.ny,self.nz,self.field,
                            self.nxp,self.nyp,self.nzp,self.fieldp,0)

  def makesourceperiodic(self):
    if self.pbounds[0] == 2 or self.pbounds[1] == 2:
      self.source[:,0,:,:] = self.source[:,0,:,:] + self.source[:,-1,:,:]
      self.source[:,-1,:,:] = self.source[:,0,:,:]
    if self.pbounds[2] == 2 or self.pbounds[3] == 2:
      self.source[:,:,0,:] = self.source[:,:,0,:] + self.source[:,:,-1,:]
      self.source[:,:,-1,:] = self.source[:,:,0,:]
    if self.pbounds[0] == 1 and not self.l4symtry:
       self.source[:,0,:,:] = 2.*self.source[:,0,:,:]
    if self.pbounds[1] == 1: self.source[:,-1,:,:] = 2.*self.source[:,-1,:,:]
    if self.pbounds[2] == 1 and not (self.l2symtry or self.l4symtry):
       self.source[:,:,0,:] = 2.*self.source[:,:,0,:]
    if self.pbounds[3] == 1: self.source[:,:,-1,:] = 2.*self.source[:,:,-1,:]
    if self.pbounds[4] == 2 or self.pbounds[5] == 2:
      if self.nslaves > 1:
        self.makesourceperiodic_parallel()
      else:
        self.source[:,:,:,0] = self.source[:,:,:,0] + self.source[:,:,:,-1]
        self.source[:,:,:,-1] = self.source[:,:,:,0]
    if self.pbounds[4] == 1: self.source[:,:,:,0] = 2.*self.source[:,:,:,0]
    if self.pbounds[5] == 1: self.source[:,:,:,-1] = 2.*self.source[:,:,:,-1]

  def makesourceperiodic_parallel(self):
    tag = 70
    if me == self.nslaves-1:
      request = mpi.isend(self.source[:,:,:,self.nz],0,tag)
      self.source[:,:,:,self.nz],status = mpi.recv(0,tag)
    elif me == 0:
      sourcetemp,status = mpi.recv(self.nslaves-1,tag)
      self.source[:,:,:,0] = self.source[:,:,:,0] + sourcetemp
      request = mpi.isend(self.source[:,:,:,0],self.nslaves-1,tag)
    if me == 0 or me == self.nslaves-1:
      status = request.wait()

  def installconductor(self,conductor,
                            xmin=None,xmax=None,
                            ymin=None,ymax=None,
                            zmin=None,zmax=None,
                            dfill=top.largepos):
    if conductor in self.conductorlist: return
    self.conductorlist.append(conductor)
    installconductors(conductor,xmin,xmax,ymin,ymax,zmin,zmax,dfill,
                      top.zbeam,
                      self.nx,self.ny,self.nz,self.nzfull,
                      self.xmmin,self.xmmax,self.ymmin,self.ymmax,
                      self.zmmin,self.zmmax,self.l2symtry,self.l4symtry,
                      solvergeom=self.solvergeom,
                      conductors=self.conductors)

  def hasconductors(self):
    return (self.conductors.interior.n > 0 or
            self.conductors.evensubgrid.n > 0 or
            self.conductors.oddsubgrid.n > 0)

  def clearconductors(self):
    self.conductors.interior.n = 0
    self.conductors.evensubgrid.n = 0
    self.conductors.oddsubgrid.n = 0

  def find_mgparam(self,lsavephi=false,resetpasses=1):
    find_mgparam(lsavephi=lsavephi,resetpasses=resetpasses,
                 solver=self,pkg3d=self)

  def dosolve(self,iwhich=0,*args):
    # --- Setup data for bends.
    rstar = fzeros(3+self.nz,'d')
    if top.bends:
      setrstar(rstar,self.nz,self.dz,self.zmmin,top.zgrid)
      self.linbend = min(rstar) < largepos

    self.source[...] = self.source*mu0*eps0

    if self.lcylindrical:
      init_bworkgrid(self.nx,self.nz,self.dx,self.dz,
                     self.xmmin,self.zmmin,self.bounds,
                     self.potential[0,:,0,:],self.source[0,:,0,:])

    # --- Note that the arrays being passed in are not contiguous, which means
    # --- that copies are being done.
    # --- If only initialization is being done (iwhich==1) then the bvp3d_work
    # --- routine only needs to be called once. Proper arrays are still passed
    # --- though they should never be needed during initialization.
    idmax = 2
    if iwhich == 1: idmax = 0
    for id in range(idmax+1):
      if (self.lanalyticbtheta and
          ((self.lusevectorpotential and (id == 0 or id == 2)) or
          (not self.lusevectorpotential and id == 1))): continue

      if self.lcylindrical:
        multigridrzb(iwhich,id,self.potential[id,:,1,:],
                     self.source[id,:,0,:],
                     self.nx,self.nz,self.mgtol[id])
      else:
        multigrid3dsolve(iwhich,self.nx,self.ny,self.nz,self.nzfull,
                         self.dx,self.dy,self.dz,
                         self.potential[id,1:-1,1:-1,:],
                         self.source[id,:,:,:],
                         rstar,self.linbend,self.bounds,
                         self.xmmin,self.ymmin,self.zmmin,
                         self.zmminglobal,
                         top.zbeam,top.zgrid,
                         self.mgparam[id],self.mgform[id],
                         self.mgiters[id],self.mgmaxiters[id],
                         self.mgmaxlevels[id],self.mgerror[id],
                         self.mgtol[id],self.mgverbose[id],
                         self.downpasses[id],self.uppasses[id],
                         self.lcndbndy,self.laddconductor,
                         self.icndbndy,false,
                         self.gridmode,self.conductors,
                         self.my_index,self.nslaves,self.izfsslave,self.nzfsslave)

    # --- This is slightly inefficient in some cases, since for example, the
    # --- MG solver already takes care of the longitudinal BC's.
    setaboundaries3d(self.potential,self.nx,self.ny,self.nz,
                     self.zmmin,self.zmmax,self.zmminglobal,self.zmmaxglobal,
                     self.bounds,self.lcylindrical,false)

    # --- Now take the curl of A to get B.
    getbfroma3d(self.potential,self.field,self.nx,self.ny,self.nz,self.dx,self.dy,self.dz,self.xmmin,
                self.lcylindrical,self.lusevectorpotential)

    # --- If using the analytic form of Btheta, calculate it here.
    if self.lanalyticbtheta:
      getanalyticbtheta(self.field,self.source,
                        self.nx,self.ny,self.nz,self.dx,self.xmmin)

    # --- Unscale the current density
    self.source[...] = self.source/(mu0*eps0)


  ##########################################################################
  # Define the basic plot commands
  def genericpf(self,kw,pffunc):
    #kw['conductors'] = self.conductors
    kw['solver'] = self
    # --- This is a temporary kludge until the plot routines are updated to
    # --- use source and potential instead of rho and phi.
    self.j = self.source
    self.b = self.field
    self.a = self.potential
    pffunc(**kw)
    del self.j
    del self.b
    del self.a
  def pcjzy(self,**kw): self.genericpf(kw,pcjzy)
  def pcjzx(self,**kw): self.genericpf(kw,pcjzx)
  def pcjxy(self,**kw): self.genericpf(kw,pcjxy)
  def pcbzy(self,**kw): self.genericpf(kw,pcbzy)
  def pcbzx(self,**kw): self.genericpf(kw,pcbzx)
  def pcbxy(self,**kw): self.genericpf(kw,pcbxy)
  def pcazy(self,**kw): self.genericpf(kw,pcazy)
  def pcazx(self,**kw): self.genericpf(kw,pcazx)
  def pcaxy(self,**kw): self.genericpf(kw,pcaxy)

# --- This can only be done after MagnetostaticMG is defined.
try:
  psyco.bind(MagnetostaticMG)
except NameError:
  pass







