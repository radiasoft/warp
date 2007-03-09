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
    self.solvergeom = w3d.RZgeom
    kw['lreducedpickle'] = lreducedpickle

    # --- Force xmmin to be zero if using RZgeom
    if self.solvergeom == w3d.RZgeom:
      self.xmmin = 0.

    # --- Force ny (which is not used here)
    self.ny = 0

    SubcycledPoissonSolver.__init__(self,kwdict=kw)
    self.solvergeom = w3d.solvergeom
    if (self.solvergeom != w3d.RZgeom and self.solvergeom != w3d.XZgeom):
      self.solvergeom = w3d.RZgeom
    self.ncomponents = 1

    # --- Kludge - make sure that the multigrid3df routines never sets up
    # --- any conductors. This is not really needed here.
    f3d.gridmode = 1

    # --- Save input parameters
    for name in MultiGridRZ.__w3dinputs__:
      if name not in self.__dict__:
        #self.__dict__[name] = kw.pop(name,getattr(w3d,name)) # Python2.3
        self.__dict__[name] = kw.get(name,getattr(w3d,name))
      if kw.has_key(name): del kw[name]
    for name in MultiGridRZ.__f3dinputs__:
      if name not in self.__dict__:
        #self.__dict__[name] = kw.pop(name,getattr(f3d,name)) # Python2.3
        self.__dict__[name] = kw.get(name,getattr(f3d,name))
      if kw.has_key(name): del kw[name]
    for name in MultiGridRZ.__frzinputs__:
      if name not in self.__dict__:
        #self.__dict__[name] = kw.pop(name,getattr(frz,name)) # Python2.3
        self.__dict__[name] = kw.get(name,getattr(frz,name))
      if kw.has_key(name): del kw[name]

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

  def loadrho(self,lzero=None,**kw):
    SubcycledPoissonSolver.loadsource(self,lzero,**kw)

  def fetche(self,*args,**kw):
    SubcycledPoissonSolver.fetchfield(self,*args,**kw)

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
    w  = pgroup.sw[js]*top.pgroup.dtscale[js]
    if top.wpid == 0: wght = zeros((0,), 'd')
    else:             wght = pgroup.pid[i:i+n,top.wpid-1]
    self.setsourcepatposition(x,y,z,ux,uy,uz,gaminv,wght,q,w,zgrid)

  def setsourcepatposition(self,x,y,z,ux,uy,uz,gaminv,wght,q,w,zgrid):
    n  = len(x)
    if n == 0: return
    sourcep = transpose(self.sourcep)
    sourcep.shape = (1+self.nzp,1,1+self.nxp)
    sourcep = transpose(sourcep)
    if top.wpid == 0:
      setrho3d(sourcep,n,x,y,z,zgrid,uz,q,w,top.depos,
               self.nxp,self.nyp,self.nzp,self.dx,self.dy,self.dz,
               self.xmminp,self.ymminp,self.zmminp,self.l2symtry,self.l4symtry,
               self.solvergeom==w3d.RZgeom)
    else:
      # --- Need top.pid(:,top.wpid)
      setrho3dw(sourcep,n,x,y,z,zgrid,uz,wght,q,w,top.depos,
                self.nxp,self.nyp,self.nzp,self.dx,self.dy,self.dz,
                self.xmminp,self.ymminp,self.zmminp,self.l2symtry,self.l4symtry,
                self.solvergeom==w3d.RZgeom)

  def fetchfieldfrompositions(self,x,y,z,ex,ey,ez,bx,by,bz,pgroup=None):
    # --- Only sets the E field from the potential
    n = len(x)
    if n == 0: return
    sete3d(self.potential,self.field,n,x,y,z,top.zgridprv,
           self.xmmin-self.dx*self.nguardx,self.ymmin,self.zmmin,
           self.dx,self.dy,self.dz,
           self.nx+2*self.nguardx,self.ny,self.nz,self.efetch,
           ex,ey,ez,self.l2symtry,self.l4symtry,self.solvergeom==w3d.RZgeom)
    ey[...] = 0.

  def fetchpotentialfrompositions(self,x,y,z,potential):
    n = len(x)
    if n == 0: return
    if self.solvergeom==w3d.RZgeom: r = sqrt(x**2 + y**2)
    else:                           r = x
    nx = self.nx + 2*self.nguardx
    nz = self.nz + 2*self.nguardz
    xmmin = self.xmmin - self.nguardx*self.dx
    xmmax = self.xmmax + self.nguardx*self.dx
    zmmin = self.zmmin - self.nguardz*self.dz
    zmmax = self.zmmax + self.nguardz*self.dz
    getgrid2d(n,r,z,potential,nx,nz,self.potential,
              xmmin,xmmax,zmmin,zmmax,
              self.l2symtry,self.l4symtry)

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
                      top.zbeam,
                      self.nx,self.ny,self.nz,self.nzfull,
                      self.xmmin,self.xmmax,self.ymmin,self.ymmax,
                      self.zmmin,self.zmmax,self.l2symtry,self.l4symtry,
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
    del self.rho,self.phi
  def pfzx(self,**kw): self.genericpf(kw,pfzx)
  def pfzxg(self,**kw): self.genericpf(kw,pfzxg)
  def pfzr(self,**kw): self.genericpf(kw,pfzx)
  def pfzrg(self,**kw): self.genericpf(kw,pfzxg)

##############################################################################
class MultiGridImplicit2D(SubcycledPoissonSolver):
  """
This solves the modified Poisson equation which includes the suseptibility
tensor that appears from the direct implicit scheme.
It currently uses the generic sparse matrix solver SuperLU. The various
multigrid input parameters are maintained for future use, but are ignored now.

Initially, conductors are not implemented.
  """
  
  __w3dinputs__ = ['iondensity','electrontemperature','plasmapotential',
                   'electrondensitymaxscale']
  __f3dinputs__ = ['gridmode','mgparam','downpasses','uppasses',
                   'mgmaxiters','mgtol','mgmaxlevels','mgform','mgverbose',
                   'lcndbndy','icndbndy','laddconductor'] 
  __frzinputs__ = ['nguardx','nguardz']

  def __init__(self,lreducedpickle=1,**kw):
    self.solvergeom = w3d.solvergeom
    kw['lreducedpickle'] = lreducedpickle
    self.nguardx = 0
    self.nguardz = 1

    # --- Force xmmin to be zero if using RZgeom
    if self.solvergeom == w3d.RZgeom:
      self.xmmin = 0.

    # --- Force ny (which is not used here)
    self.ny = 0

    SubcycledPoissonSolver.__init__(self,kwdict=kw)
    self.solvergeom = w3d.solvergeom
    if (self.solvergeom != w3d.RZgeom and self.solvergeom != w3d.XZgeom):
      self.solvergeom = w3d.RZgeom
    self.ncomponents = 1

    # --- Kludge - make sure that the multigrid3df routines never sets up
    # --- any conductors. This is not really needed here.
    f3d.gridmode = 1

    # --- Save input parameters
    for name in MultiGridImplicit2D.__w3dinputs__:
      if name not in self.__dict__:
        #self.__dict__[name] = kw.pop(name,getattr(w3d,name)) # Python2.3
        self.__dict__[name] = kw.get(name,getattr(w3d,name))
      if kw.has_key(name): del kw[name]
    for name in MultiGridImplicit2D.__f3dinputs__:
      if name not in self.__dict__:
        #self.__dict__[name] = kw.pop(name,getattr(f3d,name)) # Python2.3
        self.__dict__[name] = kw.get(name,getattr(f3d,name))
      if kw.has_key(name): del kw[name]
    for name in MultiGridImplicit2D.__frzinputs__:
      if name not in self.__dict__:
        #self.__dict__[name] = kw.pop(name,getattr(frz,name)) # Python2.3
        self.__dict__[name] = kw.get(name,getattr(frz,name))
      if kw.has_key(name): del kw[name]

    # --- If there are any remaning keyword arguments, raise an error.
    assert len(kw.keys()) == 0,"Bad keyword arguemnts %s"%kw.keys()

    self.nxguard = self.nguardx
    self.nyguard = 0
    self.nzguard = self.nguardz

    # --- Create a conductor object, which by default is empty.
    self.conductors = ConductorType()
    self.conductorlist = []
    self.newconductorlist = []

    # --- Give these variables dummy initial values.
    self.mgiters = 0
    self.mgerror = 0.

    # --- Turn of build quads option
    self.lbuildquads = false

  def __getstate__(self):
    dict = SubcycledPoissonSolver.__getstate__(self)
    if self.lreducedpickle:
      # --- Delete the conductorobject since it can be big
      del dict['conductors']
      # --- Put all of the conductors in the newconductorlist so that they
      # --- will be reinstalled after the restore.
      dict['newconductorlist'] += self.conductorlist
      dict['conductorlist'] = []
      if 'rho' in dict: del dict['rho']
      if 'phi' in dict: del dict['phi']
      if 'chi0' in dict: del dict['chi0']
    return dict

  def __setstate__(self,dict):
    SubcycledPoissonSolver.__setstate__(self,dict)
    if 'newconductorlist' not in self.__dict__:
      # --- For backwards compatibility
      self.newconductorlist = self.conductorlist
      self.conductorlist = []
    if self.lreducedpickle and not self.lnorestoreonpickle:
      # --- Create a new (and now empty) conductor object.
      # --- Any conductors will be installed when it is referenced.
      self.conductors = ConductorType()

  def getconductorobject(self):
    # --- Regenerate the conductor data
    for conductor in self.newconductorlist:
      self.installconductor(conductor)
    self.newconductorlist = []
    return self.conductors

  def getpdims(self):
    # --- This is needed to set the top.nsimplicit variable.
    setupImplicit(top.pgroup)
    # --- Returns the dimensions of the arrays used by the particles
    # --- The extra dimension is to hold the charge density and the chi's
    # --- for the implicit groups.
    return ((1+self.nxp,1+self.nzp,1+top.nsimplicit),
            (1+self.nxp+2*self.nguardx,1+self.nzp+2*self.nguardz))

  def getdims(self):
    # --- This is needed to set the top.nsimplicit variable.
    setupImplicit(top.pgroup)
    # --- Returns the dimensions of the arrays used by the field solver
    return ((1+self.nx,1+self.nz,1+top.nsimplicit),
            (1+self.nx+2*self.nguardx,1+self.nz+2*self.nguardz))

  def loadrho(self,lzero=None,**kw):
    SubcycledPoissonSolver.loadsource(self,lzero,**kw)

  def fetche(self,*args,**kw):
    if not top.lresetparticlee: return
    SubcycledPoissonSolver.fetchfield(self,*args,**kw)

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
    ss = self.sourcep.shape
    sourcep = fzeros((ss[0],1,ss[1]),'d')
    if top.wpid == 0:
      setrho3d(sourcep,n,x,y,z,zgrid,uz,q,w,top.depos,
               self.nxp,self.nyp,self.nzp,self.dx,self.dy,self.dz,
               self.xmminp,self.ymminp,self.zmminp,self.l2symtry,self.l4symtry,
               self.solvergeom==w3d.RZgeom)
    else:
      # --- Need top.pid(:,top.wpid)
      setrho3dw(sourcep,n,x,y,z,zgrid,uz,wght,q,w,top.depos,
                self.nxp,self.nyp,self.nzp,self.dx,self.dy,self.dz,
                self.xmminp,self.ymminp,self.zmminp,self.l2symtry,self.l4symtry,
                self.solvergeom==w3d.RZgeom)
    self.sourcep[...,0] += sourcep[:,0,:]
    self.sourcep[...,iimp+1] += sourcep[:,0,:]*q/m

  def fetchfieldfrompositions(self,x,y,z,ex,ey,ez,bx,by,bz,pgroup=None):
    # --- Only sets the E field from the potential
    n = len(x)
    if n == 0: return
    sete3d(self.potential,self.field,n,x,y,z,top.zgridprv,
           self.xmmin-self.dx*self.nguardx,self.ymmin,self.zmmin,
           self.dx,self.dy,self.dz,
           self.nx+2*self.nguardx,self.ny,self.nz,self.efetch,
           ex,ey,ez,self.l2symtry,self.l4symtry,self.solvergeom==w3d.RZgeom)
    ey[...] = 0.

  def fetchpotentialfrompositions(self,x,y,z,potential):
    n = len(x)
    if n == 0: return
    if self.solvergeom==w3d.RZgeom: r = sqrt(x**2 + y**2)
    else:                           r = x
    nx = self.nx + 2*self.nguardx
    nz = self.nz + 2*self.nguardz
    xmmin = self.xmmin - self.nguardx*self.dx
    xmmax = self.xmmax + self.nguardx*self.dx
    zmmin = self.zmmin - self.nguardz*self.dz
    zmmax = self.zmmax + self.nguardz*self.dz
    getgrid2d(n,r,z,potential,nx,nz,self.potential,
              xmmin,xmmax,zmmin,zmmax,
              self.l2symtry,self.l4symtry)

  def setarraysforfieldsolve(self,*args):
    SubcycledPoissonSolver.setarraysforfieldsolve(self,*args)
    if self.lparallel:
      raise "MultiGridImplicit2D not parallelized"

  def getpotentialpforparticles(self,*args):
    SubcycledPoissonSolver.getpotentialpforparticles(self,*args)
    if self.lparallel:
      raise "MultiGridImplicit2D not parallelized"

  def makesourceperiodic(self):
    if self.pbounds[0] == 2 or self.pbounds[1] == 2:
      self.source[0,:,:] = self.source[0,:,:] + self.source[-1,:,:]
      self.source[-1,:,:] = self.source[0,:,:]
    if self.pbounds[0] == 1 and not self.l4symtry:
       self.source[0,:,:] = 2.*self.source[0,:,:]
    if self.pbounds[1] == 1: self.source[-1,:,:] = 2.*self.source[-1,:,:]
    if self.pbounds[4] == 2 or self.pbounds[5] == 2:
      if self.lparallel:
        raise "MultiGridImplicit2D not parallelized"
        #self.makesourceperiodic_parallel()
      else:
        self.source[:,0,:] = self.source[:,0,:] + self.source[:,-1,:]
        self.source[:,-1,:] = self.source[:,0,:]
    if self.pbounds[4] == 1: self.source[:,0,:] = 2.*self.source[:,0,:]
    if self.pbounds[5] == 1: self.source[:,-1,:] = 2.*self.source[:,-1,:]

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
                      solvergeom=self.solvergeom,conductors=self.conductors,
                      my_index=self.my_index,nslaves=self.nslaves,
                      izfsslave=self.izfsslave,nzfsslave=self.nzfsslave)

  def hasconductors(self):
    "This is not used anywhere"
    return 1

  def clearconductors(self):
    "This is only used by realboundaries"
    pass

  def find_mgparam(self,lsavephi=false,resetpasses=1):
    "This needs to be thought through"
    self.phi = self.potential
    # --- This is a kludge to get around the fieldsolve routine in
    # --- find_mgparam having to know about 2d versus 3d arrays.
    lsavephi = true
    self.phi[...] = 0.
    find_mgparam(lsavephi=lsavephi,resetpasses=resetpasses,
                 solver=self,pkg3d=self)

  def dosolve(self,iwhich=0,*args):
    # --- Do the solve, including chi
    #self.dosolvesuperlu(iwhich,*args)
    self.dosolvemg(iwhich,*args)

  def dosolvemg(self,iwhich=0,*args):

    qomdt = top.implicitfactor*top.dt # implicitfactor = q/m
    chi0 = 0.5*self.source[...,1:]*top.dt**2/eps0
    """
    # --- Test a linearly varying chi and parabolic phi
    c1 = 10.
    c2 = 2.
    alpha = 10.
    for iz in range(self.nz+1):
      chi0[...,iz] = (c1 + c2*self.zmesh[iz])
      self.source[...,iz] = -(2.*alpha + 2.*c1*alpha + 4.*c2*alpha*w3d.zmesh[iz])*eps0
    """

    # --- This is only done for convenience.
    self.phi = self.potential
    self.rho = self.source[...,0]
    self.chi0 = chi0
    if isinstance(self.potential,FloatType): return

    mgiters = zeros(1)
    mgerror = zeros(1,'d')
    conductorobject = self.getconductorobject()

    mgsolveimplicites2d(iwhich,self.nx,self.nz,self.nzfull,self.dx,self.dz,
                        self.potential,self.source,
                        top.nsimplicit,qomdt,chi0,
                        self.bounds,self.xmmin,self.zmmin,self.zmminglobal,top.zbeam,top.zgrid,
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
    n = self.nx*self.nz
    nrhs = 1
    b = -self.source[:-1,:-1]/eps0
    phi = self.potential[1:-1,1:-1]
    info = zeros(1)

    values = fzeros((5,n),'d')
    rowind = fzeros((5,n))
    colptr = arange(n+1)*5 + 1
    rowcnt = zeros(n)
    rmmin = self.xmmin
    dr = self.dx
    dz = self.dz
    drsqi = 1./dr**2
    dzsqi = 1./dz**2
    nr = self.nx
    nz = self.nz
    coeffikm1 = dzsqi
    coeffikp1 = dzsqi
    for iz in range(0,nz):
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
          rtemp = rtemp[1:] + [rtemp[0] + nz*nr]
          vtemp = vtemp[1:] + [vtemp[0]]
        if rtemp[0] < 0:
          # --- Throw away point "below" r=0 axis at iz=0.
          del rtemp[0]
          del vtemp[0]
        if rtemp[-1] >= n:
          # --- Periodic Z boundary condition
          rtemp = [rtemp[-1] - nz*nr] + rtemp[:-1]
          vtemp = [vtemp[-1]]         + vtemp[:-1]
        if rtemp[-1] >= n:
          # --- Throw away point beyond r=nr, iz=nz
          del rtemp[-1]
          del vtemp[-1]
  
        for i in range(len(rtemp)):
          irow = rowcnt[rtemp[i]]
          values[irow,rtemp[i]] = vtemp[i]
          rowind[irow,rtemp[i]] = icol + 1
          rowcnt[rtemp[i]] += 1

    # --- There are two values of rowind that are unset, the (i-1) term
    # --- for (ix,iz)=(0,0) and the (i+1) term for (ix,iz)=(nx,nz).
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
  def genericpf(self,kw,pffunc):
    kw['conductors'] = self.getconductorobject()
    kw['solver'] = self
    # --- This is a temporary kludge until the plot routines are updated to
    # --- use source and potential instead of rho and phi.
    sdims,pdims = self.getdims()
    sdims = [sdims[1],1,sdims[0]]
    pdims = [pdims[1],1,pdims[0]]
    source = transpose(self.source[...,0])
    source.shape = sdims
    source = transpose(source)
    potential = transpose(self.potential)
    potential.shape = pdims
    potential = transpose(potential)
    self.rho = source
    self.phi = potential
    pffunc(**kw)
  def pfzx(self,**kw): self.genericpf(kw,pfzx)
  def pfzxg(self,**kw): self.genericpf(kw,pfzxg)
  def pfzr(self,**kw): self.genericpf(kw,pfzx)
  def pfzrg(self,**kw): self.genericpf(kw,pfzxg)


# --- This can only be done after MultiGridRZ and MultiGridImplicit2D are defined.
try:
  psyco.bind(MultiGridRZ)
  psyco.bind(MultiGridImplicit2D)
except NameError:
  pass

