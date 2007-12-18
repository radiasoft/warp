"""Class for doing complete multigrid field solve"""
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
class MultiGridRZ(SubcycledPoissonSolver):
  
  __w3dinputs__ = ['iondensity','electrontemperature','plasmapotential',
                   'electrondensitymaxscale']
  __f3dinputs__ = ['gridmode','mgparam','downpasses','uppasses',
                   'mgmaxiters','mgtol','mgmaxlevels','mgform','mgverbose',
                   'lcndbndy','icndbndy','laddconductor'] 
  __frzinputs__ = ['mgridrz_ncycles','nguardx','nguardz']

  def __init__(self,lreducedpickle=1,**kw):
    kw['lreducedpickle'] = lreducedpickle

    # --- Force ny (which is not used here)
    self.ny = 0

    SubcycledPoissonSolver.__init__(self,kwdict=kw)
    if (self.solvergeom != w3d.RZgeom and self.solvergeom != w3d.XZgeom):
      self.solvergeom = w3d.RZgeom
    self.ncomponents = 1

    # --- Kludge - make sure that the multigrid3df routines never sets up
    # --- any conductors. This is not really needed here.
    f3d.gridmode = 1

    # --- Save input parameters
    self.processdefaultsfrompackage(MultiGridRZ.__w3dinputs__,w3d,kw)
    self.processdefaultsfrompackage(MultiGridRZ.__f3dinputs__,f3d,kw)
    self.processdefaultsfrompackage(MultiGridRZ.__frzinputs__,frz,kw)

    # --- If there are any remaning keyword arguments, raise an error.
    assert len(kw.keys()) == 0,"Bad keyword arguemnts %s"%kw.keys()

    self.nxguard = self.nguardx
    self.nyguard = 0
    self.nzguard = self.nguardz

    # --- Create a conductor object, which by default is empty.
    self.conductors = None
    self.conductorlist = []
    self.newconductorlist = []

    # --- Give these variables dummy initial values.
    self.mgiters = 0
    self.mgerror = 0.

    # --- Initialize the grid object
    self.initializegrid()

  def __getstate__(self):
    dict = SubcycledPoissonSolver.__getstate__(self)
    if self.lreducedpickle:
      del dict['grid']
    return dict

  def __setstate__(self,dict):
    SubcycledPoissonSolver.__setstate__(self,dict)
    if 'newconductorlist' not in self.__dict__:
      self.newconductorlist = self.conductorlist
      self.conductorlist = []
    if self.lreducedpickle and not self.lnorestoreonpickle:
      self.initializegrid()

  def getconductorobject(self):
    # --- Regenerate the conductor data
    for conductor in self.newconductorlist:
      self.installconductor(conductor)
    self.newconductorlist = []
    return None

  def initializegrid(self):
    self.grid = GRIDtype()
    bg = self.grid

    inveps0 = 1./eps0

    if w3d.solvergeom==w3d.Zgeom: self.nguardx = 0
    if w3d.solvergeom==w3d.Rgeom: self.nguardz = 0

    bg.nzp = self.nz
    bg.nrpar = 0
    bg.nzpar = 0
    bg.nr = self.nx
    bg.dr = self.dx
    bg.rmin = self.xmmin
    bg.rmax = self.xmmax
    bg.xmin = self.xmmin
    bg.xmax = self.xmmax
    bg.nz = self.nz
    bg.zminp = self.zmmin
    bg.dz = self.dz
    bg.zmin = self.zmmin
    bg.zmax = self.zmmax
    bg.jmin = 1
    bg.jmax = bg.nr+1
    bg.lmin = 1
    bg.lmax = bg.nr+1
    bg.nguardx = self.nguardx
    bg.nguardz = self.nguardz

    bg.mgparam = self.mgparam
    bg.npre = self.downpasses
    bg.npost = self.uppasses
    bg.ncycles = self.mgridrz_ncycles
    bg.ncmax = self.mgmaxiters
    bg.npmin = self.mgmaxlevels
    bg.transit_min_r = 0
    bg.transit_max_r = 0
    bg.transit_min_z = 0
    bg.transit_max_z = 0
    bg.invdr = 1./bg.dr
    bg.invdz = 1./bg.dz

    bg.gid = array([1])
    bg.loc_part = ones((bg.nr+1,bg.nzp+1))
    bg.loc_part_fd = ones((bg.nr+1,bg.nzp+1))
    bg.invvol = zeros((bg.nr+1),'d')

    frz.ngrids = 1
    frz.level_del_grid = 0
    frz.n_avail_ids = 0
    frz.avail_ids = -1

    if self.solvergeom == w3d.RZgeom or self.solvergeom == w3d.Rgeom:
#     --- computes divider by cell volumes to get density
      if bg.rmin == 0.:
#       --- the factor 0.75 corrects for overdeposition due to linear
#       --- weighting (for uniform distribution)
#       --- see Larson et al., Comp. Phys. Comm., 90:260-266, 1995
#       --- and Verboncoeur, J. of Comp. Phys.,
        bg.invvol[0] = 0.75/(pi*(0.5*0.5*bg.dr*bg.dr*bg.dz))
        for j in range(1,bg.nr+1):
          bg.invvol[j] = 1./(2.*pi*j*bg.dr*bg.dr*bg.dz)
      else:
        for j in range(0,bg.nr+1):
          bg.invvol[j] = 1./(2.*pi*(bg.rmin+j*bg.dr)*bg.dr*bg.dz)
      if self.solvergeom == w3d.Rgeom: bg.invvol = bg.invvol*bg.dz
    else:
#     --- solvergeom==XZgeom or solvergeom==XYgeom
      bg.invvol = 1./(bg.dr*bg.dz)

    if self.solvergeom == w3d.RZgeom or self.solvergeom == w3d.XZgeom:
      ixlbndi = self.boundxy
      ixrbndi = self.boundxy
      izlbndi = self.bound0
      izrbndi = self.boundnz
      ixlbnd = self.boundxy
      ixrbnd = self.boundxy
      izlbnd = self.bound0
      izrbnd = self.boundnz
      if bg.rmin == 0.:
        ixlbndi = neumann
        ixlbnd  = neumann
    else:
#     --- solvergeom==XYgeom
      ixlbndi = self.boundxy
      ixrbndi = self.boundxy
      izlbndi = self.boundxy
      izrbndi = self.boundxy
      ixlbnd = self.boundxy
      ixrbnd = self.boundxy
      izlbnd = self.boundxy
      izrbnd = self.boundxy
      if self.l4symtry:
        ixlbndi = neumann
        ixlbnd  = neumann
      if (w3d.l2symtry or w3d.l4symtry) and self.solvergeom == w3d.XYgeom:
        izlbndi = neumann
        izlbnd  = neumann

    bg.ixlbnd = ixlbndi
    bg.ixrbnd = ixrbndi
    bg.izlbnd = izlbndi
    bg.izrbnd = izrbndi

    if not (self.solvergeom == w3d.Zgeom or self.solvergeom == w3d.Rgeom):
      init_gridbnd(bg)
    bg.nlevels = self.mgmaxlevels

    for i in range(bg.nlevels):
      if i == 0:
        b = bg.bndfirst
      else:
        b = b.next
      if b.izlbnd == dirichlet:  b.v[:,0]    = 3 # v_dirichlet
      if b.izrbnd == dirichlet:  b.v[:,b.nz] = 3 # v_dirichlet
      if bg.ixlbnd == dirichlet: b.v[0,:]    = 3 # v_dirichlet
      if bg.ixrbnd == dirichlet: b.v[b.nr,:] = 3 # v_dirichlet
    setmglevels_rz(bg)
  
  def getpdims(self):
    # --- Returns the dimensions of the arrays used by the particles
    return ((1+self.nxp,1+self.nzp),
            (1+self.nxp+2*self.nguardx,1+self.nzp+2*self.nguardz))

  def getdims(self):
    # --- Returns the dimensions of the arrays used by the field solver
    return ((1+self.nx,1+self.nz),
            (1+self.nx+2*self.nguardx,1+self.nz+2*self.nguardz))

  def getrho(self):
    return self.source

  def getphi(self):
    'Returns the phi array without the guard cells'
    ix1 = self.nxguard
    if ix1 == 0: ix1 = None
    ix2 = -self.nxguard
    if ix2 == 0: ix2 = None
    ix = slice(ix1,ix2)
    iz1 = self.nzguard
    if iz1 == 0: iz1 = None
    iz2 = -self.nzguard
    if iz2 == 0: iz2 = None
    iz = slice(iz1,iz2)
    return self.potential[ix,iz]

  def getfield(self):
    return self.field

  def loadrho(self,lzero=None,**kw):
    SubcycledPoissonSolver.loadsource(self,lzero,**kw)

  def fetche(self,*args,**kw):
    SubcycledPoissonSolver.fetchfield(self,*args,**kw)

  def setsourcep(self,js,pgroup,zgrid):
    n = pgroup.nps[js]
    if n == 0: return
    i  = pgroup.ins[js] - 1
    x  = pgroup.xp[i:i+n]
    y  = pgroup.yp[i:i+n]
    z  = pgroup.zp[i:i+n]
    ux = zeros((0,), 'd')
    uy = zeros((0,), 'd')
    uz = pgroup.uzp[i:i+n]
    gaminv = zeros((0,), 'd')
    q  = pgroup.sq[js]
    w  = pgroup.sw[js]*top.pgroup.dtscale[js]
    if top.wpid == 0: wght = zeros((0,), 'd')
    else:             wght = pgroup.pid[i:i+n,top.wpid-1]
    self.setsourcepatposition(x,y,z,ux,uy,uz,gaminv,wght,q,w,zgrid)

  def setsourcepatposition(self,x,y,z,ux,uy,uz,gaminv,wght,q,w,zgrid):
    n = len(x)
    if n == 0: return
    sourcep = transpose(self.sourcep)
    sourcep.shape = (1+self.nzp,1,1+self.nxp)
    sourcep = transpose(sourcep)
    if top.wpid == 0:
      setrho3d(sourcep,n,x,y,z,zgrid,uz,q,w,top.depos,
               self.nxp,self.nyp,self.nzp,self.dx,1.,self.dz,
               self.xmminp,self.ymminp,self.zmminp,self.l2symtry,self.l4symtry,
               self.solvergeom==w3d.RZgeom)
    else:
      # --- Need top.pid(:,top.wpid)
      setrho3dw(sourcep,n,x,y,z,zgrid,uz,wght,q,w,top.depos,
                self.nxp,self.nyp,self.nzp,self.dx,1.,self.dz,
                self.xmminp,self.ymminp,self.zmminp,self.l2symtry,self.l4symtry,
                self.solvergeom==w3d.RZgeom)

  def fetchfieldfrompositions(self,x,y,z,ex,ey,ez,bx,by,bz,js=0,pgroup=None):
    # --- Only sets the E field from the potential
    n = len(x)
    if n == 0: return
    sete3d(self.potentialp,self.fieldp,n,x,y,z,self.getzgridprv(),
           self.xmminp,self.ymminp,self.zmminp,
           self.dx,self.dy,self.dz,
           self.nx,self.ny,self.nz,top.efetch[js],
           ex,ey,ez,self.l2symtry,self.l4symtry,self.solvergeom==w3d.RZgeom,
           self.nxguard,self.nyguard,self.nzguard)
    #ey[...] = 0.

  def fetchpotentialfrompositions(self,x,y,z,potential):
    n = len(x)
    if n == 0: return
    if self.solvergeom==w3d.RZgeom: r = sqrt(x**2 + y**2)
    else:                           r = x
    nx = self.nx + 2*self.nguardx
    nz = self.nz + 2*self.nguardz
    xmmin = self.xmmin - self.nguardx*self.dx
    xmmax = self.xmmax + self.nguardx*self.dx
    zmminlocal = self.zmminlocal - self.nguardz*self.dz
    zmmaxlocal = self.zmmaxlocal + self.nguardz*self.dz
    getgrid2d(n,r,z,potential,nx,nz,self.potential,
              xmmin,xmmax,zmminlocal,zmmaxlocal)

  def setarraysforfieldsolve(self,*args):
    SubcycledPoissonSolver.setarraysforfieldsolve(self,*args)
    if self.lparallel:
      raise "MultiGridRZ not parallelized"

  def getpotentialpforparticles(self,*args):
    SubcycledPoissonSolver.getpotentialpforparticles(self,*args)
    if self.lparallel:
      raise "MultiGridRZ not parallelized"

  def makesourceperiodic(self):
    if self.pbounds[0] == 2 or self.pbounds[1] == 2:
      self.source[0,:] = self.source[0,:] + self.source[-1,:]
      self.source[-1,:] = self.source[0,:]
    if self.pbounds[0] == 1 and not self.l4symtry:
       self.source[0,:] = 2.*self.source[0,:]
    if self.pbounds[1] == 1: self.source[-1,:] = 2.*self.source[-1,:]
    if self.pbounds[4] == 2 or self.pbounds[5] == 2:
      if self.lparallel:
        raise "MultiGridRZ not parallelized"
        #self.makesourceperiodic_parallel()
      else:
        self.source[:,0] = self.source[:,0] + self.source[:,-1]
        self.source[:,-1] = self.source[:,0]
    if self.pbounds[4] == 1: self.source[:,0] = 2.*self.source[:,0]
    if self.pbounds[5] == 1: self.source[:,-1] = 2.*self.source[:,-1]

  def installconductor(self,conductor,
                            xmin=None,xmax=None,
                            ymin=None,ymax=None,
                            zmin=None,zmax=None,
                            dfill=top.largepos):
    if conductor in self.conductorlist: return
    self.conductorlist.append(conductor)
    installconductors(conductor,xmin,xmax,ymin,ymax,zmin,zmax,dfill,
                      self.getzgrid(),
                      self.nx,self.ny,self.nz,self.nz,
                      self.xmmin,self.xmmax,self.ymmin,self.ymmax,
                      self.zmmin,self.zmmax,1.,self.l2symtry,self.l4symtry,
                      solvergeom=self.solvergeom,gridrz=self.grid)

  def hasconductors(self):
    "This is not used anywhere"
    return 1

  def clearconductors(self):
    "This is only used by realboundaries"
    pass

  def find_mgparam(self,lsavephi=false,resetpasses=1):
    "This needs to be thought through"
    pass
    #find_mgparam(lsavephi=lsavephi,resetpasses=resetpasses,
    #             solver=self,pkg3d=self)

  def dosolve(self,iwhich=0,*args):
    self.getconductorobject()
    self.grid.rho = self.source
    self.grid.phi = self.potential
    solve_mgridrz(self.grid,self.mgtol,false)
    self.mgiters = nb_iters
    #self.mgerror = mgerror[0] # not saved anywhere

  ##########################################################################
  # Define the basic plot commands
  def genericpf(self,kw,pffunc):
    kw['conductors'] = self.getconductorobject()
    kw['solver'] = self
    # --- This is a temporary kludge until the plot routines are updated to
    # --- use source and potential instead of rho and phi.
    self.rho = self.source
    self.phi = self.potential
    pffunc(**kw)
  def pfzx(self,**kw): self.genericpf(kw,pfzx)
  def pfzxg(self,**kw): self.genericpf(kw,pfzxg)
  def pfzr(self,**kw): self.genericpf(kw,pfzx)
  def pfzrg(self,**kw): self.genericpf(kw,pfzxg)

##############################################################################
##############################################################################
##############################################################################
##############################################################################
class MultiGrid2D(MultiGrid):
  
  def __init__(self,lreducedpickle=1,**kw):
    kw['lreducedpickle'] = lreducedpickle

    # --- Force ny (which is not used here)
    self.ny = 0

    SubcycledPoissonSolver.__init__(self,kwdict=kw)
    if (self.solvergeom != w3d.RZgeom and self.solvergeom != w3d.XZgeom):
      self.solvergeom = w3d.RZgeom
    self.ncomponents = 1
    self.nxguard = 1
    self.nyguard = 0
    self.nzguard = 1

    # --- Kludge - make sure that the multigrid3df routines never sets up
    # --- any conductors. This is not really needed here.
    f3d.gridmode = 1

    # --- Save input parameters
    self.processdefaultsfrompackage(MultiGrid2D.__w3dinputs__,w3d,kw)
    self.processdefaultsfrompackage(MultiGrid2D.__f3dinputs__,f3d,kw)

    # --- If there are any remaning keyword arguments, raise an error.
    assert len(kw.keys()) == 0,"Bad keyword arguemnts %s"%kw.keys()

    # --- Create conductor objects
    self.initializeconductors()

    # --- Give these variables dummy initial values.
    self.mgiters = 0
    self.mgerror = 0.

  def getrho(self):
    return self.source[:,0,:]

  def getphi(self):
    'Returns the phi array without the guard cells'
    ix1 = self.nxguard
    if ix1 == 0: ix1 = None
    ix2 = -self.nxguard
    if ix2 == 0: ix2 = None
    ix = slice(ix1,ix2)
    iz1 = self.nzguard
    if iz1 == 0: iz1 = None
    iz2 = -self.nzguard
    if iz2 == 0: iz2 = None
    iz = slice(iz1,iz2)
    return self.potential[ix,0,iz]

  def fetchpotentialfrompositions(self,x,y,z,potential):
    n = len(x)
    if n == 0: return
    if self.solvergeom==w3d.RZgeom: r = sqrt(x**2 + y**2)
    else:                           r = x
    nx = self.nx + 2*self.nxguard
    nzlocal = self.nzlocal + 2*self.nzguard
    xmmin = self.xmmin - self.nxguard*self.dx
    xmmax = self.xmmax + self.nxguard*self.dx
    zmminlocal = self.zmminlocal - self.nzguard*self.dz
    zmmaxlocal = self.zmmaxlocal + self.nzguard*self.dz
    getgrid2d(n,r,z,potential,nx,nzlocal,self.potential[:,0,:],
              xmmin,xmmax,zmminlocal,zmmaxlocal)

  def dosolve(self,iwhich=0,*args):
    if not self.l_internal_dosolve: return
    # --- set for longitudinal relativistic contraction
    iselfb = args[2]
    beta = top.pgroup.fselfb[iselfb]/clight
    zfact = 1./sqrt((1.-beta)*(1.+beta))

    # --- This is only done for convenience.
    self.phi = self.potential
    self.rho = self.source
    if isinstance(self.potential,FloatType): return

    if self.izfsslave is None: self.izfsslave = top.izfsslave
    if self.nzfsslave is None: self.nzfsslave = top.nzfsslave
    mgiters = zeros(1,'l')
    mgerror = zeros(1,'d')
    conductorobject = self.getconductorobject(top.pgroup.fselfb[iselfb])
    self.lbuildquads = false
    multigrid2dsolve(iwhich,self.nx,self.nzlocal,self.nz,self.dx,self.dz*zfact,
                     self.phi[:,0,:],self.rho[:,0,:],self.bounds,
                     self.xmmin,self.zmminlocal*zfact,self.zmmin*zfact,
                     self.getzgrid()*zfact,self.getzgrid()*zfact,
                     self.mgparam,mgiters,self.mgmaxiters,
                     self.mgmaxlevels,mgerror,self.mgtol,self.mgverbose,
                     self.downpasses,self.uppasses,
                     self.lcndbndy,self.laddconductor,self.icndbndy,self.lbuildquads,
                     self.gridmode,conductorobject,self.solvergeom==w3d.RZgeom,
                     self.my_index,self.nslaves,self.izfsslave,self.nzfsslave)

    self.mgiters = mgiters[0]
    self.mgerror = mgerror[0]

  ##########################################################################
  # Define the basic plot commands
  def pfzr(self,**kw): self.genericpf(kw,pfzx)
  def pfzrg(self,**kw): self.genericpf(kw,pfzxg)

##############################################################################
##############################################################################
##############################################################################
##############################################################################
class MultiGridImplicit2D(MultiGrid):
  """
This solves the modified Poisson equation which includes the suseptibility
tensor that appears from the direct implicit scheme.
It currently uses the generic sparse matrix solver SuperLU. The various
multigrid input parameters are maintained for future use, but are ignored now.

Initially, conductors are not implemented.
  """
  
  def __init__(self,lreducedpickle=1,**kw):
    kw['lreducedpickle'] = lreducedpickle

    # --- Force ny (which is not used here)
    self.ny = 0

    SubcycledPoissonSolver.__init__(self,kwdict=kw)
    if (self.solvergeom != w3d.RZgeom and self.solvergeom != w3d.XZgeom):
      self.solvergeom = w3d.RZgeom
    self.ncomponents = 1
    self.nxguard = 1
    self.nyguard = 0
    self.nzguard = 1

    # --- Kludge - make sure that the multigrid3df routines never sets up
    # --- any conductors. This is not really needed here.
    f3d.gridmode = 1

    # --- Save input parameters
    self.processdefaultsfrompackage(MultiGrid2D.__w3dinputs__,w3d,kw)
    self.processdefaultsfrompackage(MultiGrid2D.__f3dinputs__,f3d,kw)

    # --- If there are any remaning keyword arguments, raise an error.
    assert len(kw.keys()) == 0,"Bad keyword arguemnts %s"%kw.keys()

    # --- Create conductor objects
    self.initializeconductors()

    # --- Give these variables dummy initial values.
    self.mgiters = 0
    self.mgerror = 0.

    # --- Turn on the chi kludge, where chi is set to be an average value
    # --- of chi for grid cells where is it zero.
    self.chikludge = 1

  def __getstate__(self):
    dict = MultiGrid.__getstate__(self)
    if self.lreducedpickle:
      if 'chi0' in dict: del dict['chi0']
    return dict

  def getpdims(self):
    # --- This is needed to set the top.nsimplicit variable.
    setupImplicit(top.pgroup)
    dims = MultiGrid.getpdims(self)
    # --- The extra dimension is to hold the charge density and the chi's
    # --- for the implicit groups.
    dims = (tuple(list(dims[0])+[1+top.nsimplicit]),)+dims[1:]
    return dims

  def getdims(self):
    # --- This is needed to set the top.nsimplicit variable.
    setupImplicit(top.pgroup)
    dims = MultiGrid.getdims(self)
    # --- The extra dimension is to hold the charge density and the chi's
    # --- for the implicit groups.
    dims = (tuple(list(dims[0])+[1+top.nsimplicit]),)+dims[1:]
    return dims

  def getrho(self):
    return self.source[:,0,:,0]

  def getphi(self):
    'Returns the phi array without the guard cells'
    return MultiGrid.getphi(self)[:,0,:]

  def loadrho(self,lzero=None,**kw):
    # --- top.laccumulate_rho is used as a flag by the implicit stepper.
    # --- When true, the load rho is skipped - it is not needed at some
    # --- points during a step.
    if top.laccumulate_rho: return
    MultiGrid.loadsource(self,lzero,**kw)

  def fetche(self,*args,**kw):
    # --- lresetparticlee is used as a flag in the implicit stepper.
    # --- When false, skip the fetche since the field is calculated
    # --- from existing data.
    if not top.lresetparticlee: return
    MultiGrid.fetchfield(self,*args,**kw)

  def setsourcep(self,js,pgroup,zgrid):
    n  = pgroup.nps[js]
    if n == 0: return
    i  = pgroup.ins[js] - 1
    x  = pgroup.xp[i:i+n]
    y  = pgroup.yp[i:i+n]
    z  = pgroup.zp[i:i+n]
    ux = zeros((0,), 'd')
    uy = zeros((0,), 'd')
    uz = pgroup.uzp[i:i+n]
    gaminv = zeros((0,), 'd')
    q  = pgroup.sq[js]
    m  = pgroup.sm[js]
    w  = pgroup.sw[js]*top.pgroup.dtscale[js]
    iimp = pgroup.iimplicit[js]
    if top.wpid == 0: wght = zeros((0,), 'd')
    else:             wght = pgroup.pid[i:i+n,top.wpid-1]
    self.setsourcepatposition(x,y,z,ux,uy,uz,gaminv,wght,q,m,w,iimp,zgrid)

  def setsourcepatposition(self,x,y,z,ux,uy,uz,gaminv,wght,q,m,w,iimp,zgrid):
    n  = len(x)
    if n == 0: return
    # --- Create a temporary array to pass into setrho3d. This contributes
    # --- differently to the charge density and to chi. Also, make it a
    # --- 3-D array so it is accepted by setrho3d.
    sourcep = fzeros(self.sourcep.shape[:-1],'d')
    if top.wpid == 0:
      setrho3d(sourcep,n,x,y,z,zgrid,uz,q,w,top.depos,
               self.nxp,self.nyp,self.nzp,self.dx,1.,self.dz,
               self.xmminp,self.ymminp,self.zmminp,self.l2symtry,self.l4symtry,
               self.solvergeom==w3d.RZgeom)
    else:
      # --- Need top.pid(:,top.wpid)
      setrho3dw(sourcep,n,x,y,z,zgrid,uz,wght,q,w,top.depos,
                self.nxp,self.nyp,self.nzp,self.dx,1.,self.dz,
                self.xmminp,self.ymminp,self.zmminp,self.l2symtry,self.l4symtry,
                self.solvergeom==w3d.RZgeom)
    self.sourcep[...,0] += sourcep
    if iimp >= 0:
      # --- The extra terms convert rho to chi
      self.sourcep[...,iimp+1] += 0.5*sourcep*q/m*top.dt**2/eps0

  def setsourceforfieldsolve(self,*args):
    # --- A separate copy is needed since self.source has an extra dimension
    # --- which must be looped over.
    SubcycledPoissonSolver.setsourceforfieldsolve(self,*args)
    if self.lparallel:
      SubcycledPoissonSolver.setsourcepforparticles(self,*args)
      if isinstance(self.source,FloatType): return
      if isinstance(self.sourcep,FloatType): return
      for iimp in range(top.nsimplicit):
        setrhoforfieldsolve3d(self.nx,self.ny,self.nzlocal,self.source[...,iimp],
                              self.nxp,self.nyp,self.nzp,self.sourcep[...,iimp],
                              self.nzpguard,
                              self.my_index,self.nslaves,self.izpslave,self.nzpslave,
                              self.izfsslave,self.nzfsslave)

  def fetchfieldfrompositions(self,x,y,z,ex,ey,ez,bx,by,bz,js=0,pgroup=None):
    MultiGrid.fetchfieldfrompositions(self,x,y,z,ex,ey,ez,bx,by,bz,js,pgroup)
    # --- Force ey to zero (is this really needed?)
    #ey[...] = 0.

  def fetchpotentialfrompositions(self,x,y,z,potential):
    n = len(x)
    if n == 0: return
    if self.solvergeom==w3d.RZgeom: r = sqrt(x**2 + y**2)
    else:                           r = x
    nx = self.nx + 2*self.nxguard
    nzlocal = self.nzlocal + 2*self.nzguard
    xmmin = self.xmmin - self.nxguard*self.dx
    xmmax = self.xmmax + self.nxguard*self.dx
    zmminlocal = self.zmminlocal - self.nzguard*self.dz
    zmmaxlocal = self.zmmaxlocal + self.nzguard*self.dz
    getgrid2d(n,r,z,potential,nx,nzlocal,self.potential[:,0,:],
              xmmin,xmmax,zmminlocal,zmmaxlocal)

  def dosolve(self,iwhich=0,*args):
    if not self.l_internal_dosolve: return
    # --- Do the solve, including chi
    #self.dosolvesuperlu(iwhich,*args)
    self.dosolvemg(iwhich,*args)

  def dosolvemg(self,iwhich=0,*args):
    # --- set for longitudinal relativistic contraction
    iselfb = args[2]
    beta = top.pgroup.fselfb[iselfb]/clight
    zfact = 1./sqrt((1.-beta)*(1.+beta))

    # --- This is only done for convenience.
    self.phi = self.potential
    self.rho = self.source[...,0]
    if isinstance(self.potential,FloatType): return

    if self.izfsslave is None: self.izfsslave = top.izfsslave
    if self.nzfsslave is None: self.nzfsslave = top.nzfsslave
    mgiters = zeros(1,'l')
    mgerror = zeros(1,'d')
    conductorobject = self.getconductorobject(top.pgroup.fselfb[iselfb])
    self.lbuildquads = false

    # --- Setup implicit chi
    qomdt = top.implicitfactor*top.dt # implicitfactor = q/m
    #--- chi0 = 0.5*rho*q/m*top.dt**2/eps0
    self.chi0 = self.source[...,1:]
    # --- Kludge alart!!!
    if self.chikludge:
      for js in range(self.source.shape[-1]-1):
        if maxnd(abs(self.chi0[...,js])) == 0.: continue
        avechi = sumnd(self.chi0[...,js])/sumnd(where(self.chi0[...,js] == 0.,0.,1.))
        self.chi0[...,js] = where(self.chi0[...,js]==0.,avechi,self.chi0[...,js])
    """
    # --- Test a linearly varying chi and parabolic phi
    c1 = 10.
    c2 = 2.
    alpha = 10.
    for iz in range(self.nzlocal+1):
      self.chi0[...,iz] = (c1 + c2*self.zmesh[iz])
      self.source[...,iz] = -(2.*alpha + 2.*c1*alpha + 4.*c2*alpha*w3d.zmesh[iz])*eps0
    """

    mgsolveimplicites2d(iwhich,self.nx,self.nzlocal,self.nz,self.dx,self.dz*zfact,
                        self.potential,self.source,
                        top.nsimplicit,qomdt,self.chi0,
                        self.bounds,self.xmmin,self.zmminlocal*zfact,self.zmmin*zfact,
                        self.getzgrid()*zfact,self.getzgrid()*zfact,
                        self.mgparam,mgiters,self.mgmaxiters,
                        self.mgmaxlevels,mgerror,self.mgtol,self.mgverbose,
                        self.downpasses,self.uppasses,
                        self.lcndbndy,self.laddconductor,self.icndbndy,self.lbuildquads,
                        self.gridmode,conductorobject,self.solvergeom==w3d.RZgeom,
                        self.my_index,self.nslaves,self.izfsslave,self.nzfsslave)

    self.mgiters = mgiters[0]
    self.mgerror = mgerror[0]

  def dosolvesuperlu(self,iwhich=0,*args):
    "Note that this does not actually include the implicit susecptibility"
    #self.grid.rho = self.source
    #self.grid.phi = self.potential
    #solve_mgridrz(self.grid,self.mgtol,false)
    #self.mgiters = nb_iters
    ##self.mgerror = mgerror[0] # not saved anywhere
    self.getconductorobject()

    # --- Use direct matrix solver
    t0 = wtime()
    n = self.nx*self.nzlocal
    nrhs = 1
    b = -self.source[:-1,:-1]/eps0
    phi = self.potential[1:-1,1:-1]
    info = zeros(1,'l')

    values = fzeros((5,n),'d')
    rowind = fzeros((5,n),'l')
    colptr = arange(n+1)*5 + 1
    rowcnt = zeros(n,'l')
    rmmin = self.xmmin
    dr = self.dx
    dz = self.dz
    drsqi = 1./dr**2
    dzsqi = 1./dz**2
    nr = self.nx
    nzlocal = self.nzlocal
    coeffikm1 = dzsqi
    coeffikp1 = dzsqi
    for iz in range(0,nzlocal):
      for ix in range(0,nr):
        icol = iz*nr + ix
        r = rmmin + ix*dr
        if r == 0.:
          coeffik = - 4.*drsqi - 2.*dzsqi
          coeffim1k = 0
          coeffip1k = 4.*drsqi
        else:
          coeffik = - 2.*drsqi - 2.*dzsqi
          coeffim1k = (r-0.5*dr)/r*drsqi
          coeffip1k = (r+0.5*dr)/r*drsqi
          if ix == nr-1:
            b[ix,iz] += -coeffip1k*phi[ix+1,iz]
            coeffip1k = 0.

        vtemp = [coeffikm1,coeffim1k,coeffik,coeffip1k,coeffikp1]
        rtemp = [-nr+icol,-1+icol,0+icol,+1+icol,+nr+icol]
        if rtemp[0] < 0:
          # --- Periodic Z boundary condition
          rtemp = rtemp[1:] + [rtemp[0] + nzlocal*nr]
          vtemp = vtemp[1:] + [vtemp[0]]
        if rtemp[0] < 0:
          # --- Throw away point "below" r=0 axis at iz=0.
          del rtemp[0]
          del vtemp[0]
        if rtemp[-1] >= n:
          # --- Periodic Z boundary condition
          rtemp = [rtemp[-1] - nzlocal*nr] + rtemp[:-1]
          vtemp = [vtemp[-1]]         + vtemp[:-1]
        if rtemp[-1] >= n:
          # --- Throw away point beyond r=nr, iz=nzlocal
          del rtemp[-1]
          del vtemp[-1]
  
        for i in range(len(rtemp)):
          irow = rowcnt[rtemp[i]]
          values[irow,rtemp[i]] = vtemp[i]
          rowind[irow,rtemp[i]] = icol + 1
          rowcnt[rtemp[i]] += 1

    # --- There are two values of rowind that are unset, the (i-1) term
    # --- for (ix,iz)=(0,0) and the (i+1) term for (ix,iz)=(nx,nzlocal).
    # --- Give the first one a fake value (since the coefficient is zero
    # --- anway.
    rowind[-1,0] = rowind[-2,0] + 1
    # --- The other is ignored by decrementing the last value of colptr.
    colptr[-1] -= 1

    self.values = values
    self.rowind = rowind
    self.colptr = colptr
    nnz = colptr[-1] - 1

    t1 = wtime()
    superlu_dgssv(n,nnz,nrhs,values,rowind,colptr,b,info)
    t2 = wtime()

    self.potential[1:-2,1:-2] = b
    self.potential[0,1:-2] = self.potential[2,1:-2]
    self.potential[-1,1:-2] = 2*self.potential[-2,1:-2]-self.potential[-3,1:-2]
    self.potential[:,-2:] = self.potential[:,1:3]
    self.potential[:,0] = self.potential[:,-3]
    t3 = wtime()

    print "Solve time = ",t2 - t1
    print "Total time = ",t3 - t0
    self.fstime = t2 - t1
    self.tottime = t3 - t0

  ##########################################################################
  # Define the basic plot commands
  def pfzr(self,**kw): self.genericpf(kw,pfzx)
  def pfzrg(self,**kw): self.genericpf(kw,pfzxg)


# --- This can only be done after MultiGridRZ and MultiGridImplicit2D are defined.
try:
  psyco.bind(MultiGridRZ)
  psyco.bind(MultiGridImplicit2D)
except NameError:
  pass

