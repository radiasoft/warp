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
    self.conductorlist = []

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
    if self.lreducedpickle and not self.lnorestoreonpickle:
      self.initializegrid()
      # --- Regenerate the conductor data
      conductorlist = self.conductorlist
      self.conductorlist = []
      for conductor in conductorlist:
        self.installconductor(conductor)

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
            (1+self.nx+2*self.nguardx,3+self.nz+2*self.nguardz))

  def loadrho(self,lzero=None,**kw):
    SubcycledPoissonSolver.loadsource(self,lzero,**kw)

  def fetche(self,*args):
    SubcycledPoissonSolver.fetchfield(self,*args)

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
    if top.wpid == 0:
      rhoweightrzgrid(self.grid,x,y,z,n,q*w,self.nx,self.nz,self.dx,self.dz,
                      self.xmmin,top.zgrid)
    else:
      # --- Need top.pid(:,top.wpid)
      rhoweightrzgrid_weights(self.grid,x,y,z,wght,n,q*w,self.nx,self.nz,self.dx,self.dz,
                      self.xmmin,top.zgrid)
      

  def fetchfieldfrompositions(self,x,y,z,ex,ey,ez,bx,by,bz,pgroup=None):
    # --- Only sets the E field from the potential
    n = len(x)
    if n == 0: return
    sete3d(self.phi,self.selfe,n,x,y,z,top.zgridprv,
           self.xmmin,self.ymmin,self.zmmin,
           self.dx,self.dy,self.dz,self.nx,self.ny,self.nz,self.efetch,
           ex,ey,ez,self.l2symtry,self.l4symtry,self.solvergeom==w3d.RZgeom)

  def fetchpotentialfrompositions(self,x,y,z,phi):
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
    getgrid2d(n,r,z,phi,nx,nz,self.phi,
              xmmin,xmmax,zmmin,zmmax,
              self.l2symtry,self.l4symtry)

  def setarraysforfieldsolve(self,*args):
    if self.nslaves <= 1:
      SubcycledPoissonSolver.setarraysforfieldsolve(self,*args)
    else:
      raise "MultiGridRZ not parallelized"

  def getpotentialpforparticles(self,*args):
    if self.nslaves <= 1:
      SubcycledPoissonSolver.getpotentialpforparticles(self,*args)
    else:
      raise "MultiGridRZ not parallelized"

  def makesourceperiodic(self):
    if self.pbounds[0] == 2 or self.pbounds[1] == 2:
      self.source[0,:] = self.source[0,:] + self.source[-1,:]
      self.source[-1,:] = self.source[0,:]
    if self.pbounds[0] == 1 and not self.l4symtry:
       self.source[0,:] = 2.*self.source[0,:]
    if self.pbounds[1] == 1: self.source[-1,:] = 2.*self.source[-1,:]
    if self.pbounds[4] == 2 or self.pbounds[5] == 2:
      if self.nslaves > 1:
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
   #return (self.conductors.interior.n > 0 or
   #        self.conductors.evensubgrid.n > 0 or
   #        self.conductors.oddsubgrid.n > 0)

  def clearconductors(self):
    "This is only used by realboundaries"
    pass
    #self.conductors.interior.n = 0
    #self.conductors.evensubgrid.n = 0
    #self.conductors.oddsubgrid.n = 0

  def find_mgparam(self,lsavephi=false,resetpasses=1):
    "This needs to be thought through"
    pass
    #find_mgparam(lsavephi=lsavephi,resetpasses=resetpasses,
    #             solver=self,pkg3d=self)

  def dosolve(self,iwhich=0,*args):
    self.grid.rho = self.source
    self.grid.phi = self.potential
    solve_mgridrz(self.grid,self.mgtol,false)
    self.mgiters = nb_iters
    #self.mgerror = mgerror[0] # not saved anywhere

  ##########################################################################
  # Define the basic plot commands
  def genericpf(self,kw,pffunc):
    kw['conductors'] = self.conductors
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


# --- This can only be done after MultiGridRZ is defined.
try:
  psyco.bind(MultiGridRZ)
except NameError:
  pass

