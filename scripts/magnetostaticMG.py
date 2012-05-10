"""Class for doing complete magnetostatic multigrid field solve"""
# ToDo:
#  - modify setj to check if particles are within grid
from warp import *
from find_mgparam import find_mgparam

try:
  import psyco
except ImportError:
  pass

magnetostaticMG_version = "$Id: magnetostaticMG.py,v 1.46 2011/11/22 20:27:37 grote Exp $"

##############################################################################
class MagnetostaticMG(SubcycledPoissonSolver):
  
  __bfieldinputs__ = ['mgparam','downpasses','uppasses',
                      'mgmaxiters','mgtol','mgmaxlevels','mgform','mgverbose',
                      'lcndbndy','icndbndy','laddconductor',
                      'lcylindrical','lanalyticbtheta'] 
  __f3dinputs__ = ['gridmode','mgparam','downpasses','uppasses',
                   'mgmaxiters','mgtol','mgmaxlevels','mgform','mgverbose',
                   'lcndbndy','icndbndy','laddconductor','lprecalccoeffs'] 

  def __init__(self,luse2D=false,**kw):
    self.grid_overlap = 2
    self.luse2D = luse2D

    # --- Save input parameters
    self.processdefaultsfrompackage(MagnetostaticMG.__f3dinputs__,f3d,kw)
    self.processdefaultsfrompackage(MagnetostaticMG.__bfieldinputs__,
                                    f3d.bfield,kw)

    # --- Check for cylindrical geometry
    if self.lcylindrical: self.solvergeom = w3d.RZgeom

    SubcycledPoissonSolver.__init__(self,kwdict=kw)
    if self.ny == 0 and self.solvergeom not in [w3d.RZgeom,w3d.XZgeom]:
      self.solvergeom = w3d.RZgeom
    if self.solvergeom in [w3d.RZgeom,w3d.XZgeom]:
      self.nyguardphi = 0
      self.nyguardrho = 0
      self.nyguarde   = 0

    self.ncomponents = 3
    self.lusevectorpotential = true

    # --- Kludge - make sure that the multigrid3df routines never sets up
    # --- any conductors.
    f3d.gridmode = 1

    # --- If there are any remaning keyword arguments, raise an error.
    assert len(kw.keys()) == 0,"Bad keyword arguemnts %s"%kw.keys()

    # --- Create a conductor object, which by default is empty.
    self.conductors = ConductorType()
    self.conductorlist = []
    self.newconductorlist = []

    # --- Give these variables dummy initial values.
    self.mgiters = zeros(3,'l')
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

  def __getstate__(self):
    dict = SubcycledPoissonSolver.__getstate__(self)
    if self.lreducedpickle:
      del dict['conductors']
      dict['newconductorlist'] += self.conductorlist
      dict['conductorlist'] = []
    return dict

  def __setstate__(self,dict):
    SubcycledPoissonSolver.__setstate__(self,dict)
    if 'newconductorlist' not in self.__dict__:
      self.newconductorlist = self.conductorlist
      self.conductorlist = []
    if self.lreducedpickle and not self.lnorestoreonpickle:
      # --- Regenerate the conductor data
      self.conductors = ConductorType()

  def getconductorobject(self):
    for conductor in self.newconductorlist:
      self.installconductor(conductor)
    self.newconductorlist = []
    return self.conductors

  def getpdims(self):
    # --- Returns the dimensions of the jp, bp, and ap arrays
    return ((3,1+self.nxp+2*self.nxguardrho,
               1+self.nyp+2*self.nyguardrho,
               1+self.nzp+2*self.nzguardrho),
            (3,1+self.nxp+2*self.nxguarde,
               1+self.nyp+2*self.nyguarde,
               1+self.nzp+2*self.nzguarde),
            (3,1+self.nxp+2*self.nxguardphi,
               1+self.nyp+2*self.nyguardphi,
               1+self.nzp+2*self.nzguardphi))

  def getdims(self):
    # --- Returns the dimensions of the j, b, and a arrays
    return ((3,1+self.nxlocal+2*self.nxguardrho,
               1+self.nylocal+2*self.nyguardrho,
               1+self.nzlocal+2*self.nzguardrho),
            (3,1+self.nxlocal+2*self.nxguarde,
               1+self.nylocal+2*self.nyguarde,
               1+self.nzlocal+2*self.nzguarde),
            (3,1+self.nxlocal+2*self.nxguardphi,
               1+self.nylocal+2*self.nyguardphi,
               1+self.nzlocal+2*self.nzguardphi))

  def getj(self):
    'Returns the current density array'
    return self.source[:,self.nxguardrho:-self.nxguardrho or None,
                         self.nyguardrho:-self.nyguardrho or None,
                         self.nzguardrho:-self.nzguardrho or None]

  def getb(self):
    'Returns the B field array'
    return self.field[:,self.nxguarde:-self.nxguarde or None,
                        self.nyguarde:-self.nyguarde or None,
                        self.nzguarde:-self.nzguarde or None]
  
  def geta(self):
    'Returns the a array without the guard cells'
    return self.potential[:,self.nxguardphi:-self.nxguardphi or None,
                            self.nyguardphi:-self.nyguardphi or None,
                            self.nzguardphi:-self.nzguardphi or None]

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
    if top.wpid > 0: wght = top.pgroup.pid[i:i+n,top.wpid-1]
    else:            wght = zeros((0,),'d')
    self.setsourcepatposition(x,y,z,ux,uy,uz,gaminv,wght,zgrid,q,w)

  def setsourcepatposition(self,x,y,z,ux,uy,uz,gaminv,wght,zgrid,q,w):
    n = len(x)
    if n == 0: return
    if len(wght) > 0:
      nw = len(wght)
    else:
      nw = 0.
      wght = zeros(1,'d')
    setj3d(self.sourcep,self.sourcep,n,x,y,z,zgrid,ux,uy,uz,gaminv,
           q,w,nw,wght,top.depos,
           self.nxp,self.nyp,self.nzp,
           self.nxguardrho,self.nyguardrho,self.nzguardrho,
           self.dx,self.dy,self.dz,
           self.xmminp,self.ymminp,self.zmminp,
           self.l2symtry,self.l4symtry,self.solvergeom==w3d.RZgeom)

  def fetchfieldfrompositions(self,x,y,z,ex,ey,ez,bx,by,bz,js=0,pgroup=None):
    n = len(x)
    if n == 0: return
    setb3d(self.fieldp,n,x,y,z,self.getzgridprv(),bx,by,bz,
           self.nxp,self.nyp,self.nzp,
           self.nxguarde,self.nyguarde,self.nzguarde,
           self.dx,self.dy,self.dz,
           self.xmminp,self.ymminp,self.zmminp,
           self.l2symtry,self.l4symtry,self.solvergeom==w3d.RZgeom)

  def fetchpotentialfrompositions(self,x,y,z,a):
    n = len(x)
    if n == 0: return
    fetchafrompositions3d(self.potentialp,n,x,y,z,a,self.getzgrid(),
                          self.nxp,self.nyp,self.nzp,
                          self.nxguardphi,self.nyguardphi,self.nzguardphi,
                          self.dx,self.dy,self.dz,
                          self.xmminp,self.ymminp,self.zmminp,
                          self.l2symtry,self.l4symtry,
                          self.solvergeom==w3d.RZgeom)

  def setsourceforfieldsolve(self,*args):
    SubcycledPoissonSolver.setsourceforfieldsolve(self,*args)
    if self.lparallel:
      SubcycledPoissonSolver.setsourcepforparticles(self,*args)
      setjforfieldsolve3d(self.nxlocal,self.nylocal,self.nzlocal,self.source,
                          self.nxp,self.nyp,self.nzp,self.sourcep,
                          self.nxguardrho,self.nyguardrho,self.nzguardrho,
                          self.fsdecomp,self.ppdecomp)

  def getpotentialpforparticles(self,*args):
    """Despite the name, this actually gets the field instead, since that is
       always used in the magnetostatic solver"""
    if not self.lparallel:
      SubcycledPoissonSolver.getfieldpforparticles(self,*args)
    else:
      self.setfieldpforparticles(*args)
      getphipforparticles3d(3,self.nxlocal,self.nylocal,self.nzlocal,
                            self.nxguarde,self.nyguarde,self.nzguarde,
                            self.field,
                            self.nxp,self.nyp,self.nzp,self.fieldp,
                            self.fsdecomp,self.ppdecomp)

  def applysourceboundaryconditions(self):
    applyrhoboundaryconditions3d(self.ncomponents,
                             self.nxlocal,self.nylocal,self.nzlocal,
                             self.nxguardrho,self.nyguardrho,self.nzguardrho,
                             self.source,self.bounds,self.fsdecomp,
                             self.solvergeom==w3d.RZgeom)

  def installconductor(self,conductor,
                            xmin=None,xmax=None,
                            ymin=None,ymax=None,
                            zmin=None,zmax=None,
                            dfill=top.largepos):
    if conductor in self.conductorlist: return
    self.conductorlist.append(conductor)
    installconductors(conductor,xmin,xmax,ymin,ymax,zmin,zmax,dfill,
                      self.getzgrid(),
                      self.nx,self.ny,self.nz,
                      self.nxlocal,self.nylocal,self.nzlocal,
                      self.xmmin,self.xmmax,self.ymmin,self.ymmax,
                      self.zmmin,self.zmmax,1.,self.l2symtry,self.l4symtry,
                      solvergeom=self.solvergeom,
                      conductors=self.conductors,decomp=self.fsdecomp)

  def hasconductors(self):
    conductorobject = self.getconductorobject()
    return (conductorobject.interior.n > 0 or
            conductorobject.evensubgrid.n > 0 or
            conductorobject.oddsubgrid.n > 0)

  def clearconductors(self):
    self.conductors.interior.n = 0
    self.conductors.evensubgrid.n = 0
    self.conductors.oddsubgrid.n = 0

  def find_mgparam(self,lsavephi=false,resetpasses=1):
    find_mgparam(lsavephi=lsavephi,resetpasses=resetpasses,
                 solver=self,pkg3d=self)

  def dosolve(self,iwhich=0,*args):
    # --- Setup data for bends.
    rstar = fzeros(3+self.nzlocal,'d')
    if top.bends:
      setrstar(rstar,self.nzlocal,self.dz,self.zmminlocal,self.getzgrid())
      self.linbend = min(rstar) < largepos

    self.source[...] = self.source*mu0*eps0
    conductorobject = self.getconductorobject()

    if self.solvergeom == w3d.RZgeom and not self.luse2D:
      init_bworkgrid(self.nxlocal,self.nzlocal,self.dx,self.dz,
                     self.xmminlocal,self.zmminlocal,self.bounds,
                     self.lparallel)

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

      if ((self.luse2D and self.solvergeom == w3d.RZgeom) or
          self.solvergeom == w3d.XZgeom):
        bounds = self.bounds.copy()
        if self.solvergeom == w3d.RZgeom and id < 2 and self.xmminlocal == 0.:
          lmagnetostaticrz = (id < 2)
          bounds[0] = dirichlet
          self.potential[id,:self.nxguardphi+1,:,:] = 0.
        else:
          lmagnetostaticrz = false
        multigrid2dsolve(iwhich,self.nx,self.nz,self.nxlocal,self.nzlocal,
                         self.nxguardphi,self.nzguardphi,
                         self.nxguardrho,self.nzguardrho,
                         self.dx,self.dz,
                         self.potential[id,:,self.nyguardphi,:],
                         self.source[id,:,self.nyguardrho,:],
                         bounds,self.xmminlocal,
                         self.mgparam[id],self.mgform[id],
                         self.mgiters[id],self.mgmaxiters[id],
                         self.mgmaxlevels[id],self.mgerror[id],
                         self.mgtol[id],self.mgverbose[id],
                         self.downpasses[id],self.uppasses[id],
                         self.lcndbndy,self.laddconductor,self.icndbndy,
                         self.gridmode,conductorobject,self.solvergeom==w3d.RZgeom,
                         lmagnetostaticrz,self.fsdecomp)
      elif (not self.luse2D) and self.solvergeom == w3d.RZgeom:
          multigridrzb(iwhich,id,self.potential[id,
                                self.nxguardphi-1:(1-self.nxguardphi) or None,self.nyguardphi,
                                self.nzguardphi-1:(1-self.nzguardphi) or None],
                     self.source[id,self.nxguardrho:-self.nxguardrho or None,self.nyguardrho,
                                    self.nzguardrho:-self.nzguardrho or None],
                     self.nxlocal,self.nzlocal,self.mgtol[id])
      else:
        multigrid3dsolve(iwhich,self.nx,self.ny,self.nz,
                         self.nxlocal,self.nylocal,self.nzlocal,
                         self.nxguardphi,self.nyguardphi,self.nzguardphi,
                         self.nxguardrho,self.nyguardrho,self.nzguardrho,
                         self.dx,self.dy,self.dz,
                         self.potential[id,:,:,:],
                         self.source[id,:,:,:],
                         rstar,self.linbend,self.bounds,
                         self.xmmin,self.ymmin,self.zmmin,
                         self.mgparam[id],self.mgform[id],
                         self.mgiters[id],self.mgmaxiters[id],
                         self.mgmaxlevels[id],self.mgerror[id],
                         self.mgtol[id],self.mgverbose[id],
                         self.downpasses[id],self.uppasses[id],
                         self.lcndbndy,self.laddconductor,self.icndbndy,
                         self.gridmode,conductorobject,self.lprecalccoeffs,
                         self.fsdecomp)

  # # --- This is slightly inefficient in some cases, since for example, the
  # # --- MG solver already takes care of the longitudinal BC's.
  # setaboundaries3d(self.potential,self.nx,self.ny,self.nzlocal,
  #                  self.zmminlocal,self.zmmaxlocal,self.zmmin,self.zmmax,
  #                  self.bounds,self.solvergeom==w3d.RZgeom,false)

    # --- Now take the curl of A to get B.
    getbfroma3d(self.potential,self.field,
                self.nxlocal,self.nylocal,self.nzlocal,
                self.nxguardphi,self.nyguardphi,self.nzguardphi,
                self.nxguarde,self.nyguarde,self.nzguarde,
                self.dx,self.dy,self.dz,self.xmminlocal,
                self.solvergeom==w3d.RZgeom,self.lusevectorpotential)

    # --- If using the analytic form of Btheta, calculate it here.
    if self.lanalyticbtheta:
      getanalyticbtheta(self.field,self.source,
                        self.nxlocal,self.nylocal,self.nzlocal,
                        self.nxguarde,self.nyguarde,self.nzguarde,
                        self.nxguardrho,self.nyguardrho,self.nzguardrho,
                        self.dx,self.xmminlocal)

    # --- Unscale the current density
    self.source[...] = self.source/(mu0*eps0)


  ##########################################################################
  # Define the basic plot commands
  def genericpf(self,kw,pffunc):
    #kw['conductors'] = self.getconductorobject()
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

##############################################################################
class MagnetostaticFFT(MagnetostaticMG):
  
  def __init__(self,**kw):

    # --- Force periodic boundary conditions in z, since that is the only
    # --- boundary conditions currently supported.
    self.bound0 = periodic
    self.boundnz = periodic

    MagnetostaticMG.__init__(self,**kw)

    # --- Only cylindrical geometry is supported
    if self.solvergeom != w3d.RZgeom:
      raise Exception("MagnetostaticFFT only supports cylindrial geometry")

  def dosolve(self,iwhich=0,*args):

    self.source[...] = self.source*mu0*eps0

    lr = self.xmmax - self.xmmin
    lz = self.zmmax - self.zmmin
    kzsq = fzeros(1+self.nz)
    attz = fzeros(1+self.nz/2)
    filt = w3d.filt[:,2]
    rfsmat = fzeros((1+self.nx,3,1+self.nz))
    scrtch2 = fzeros((1+self.nz))

    # --- Initialize kzsq and attz.
    vpoisrzb(1,self.potential[0,...],kzsq,attz,filt,lr,lz,self.nx,self.nz,
             rfsmat,scrtch2,0)
    if iwhich == 1: return

    for axis in range(3):
      if (self.lanalyticbtheta and
          ((self.lusevectorpotential and (axis == 0 or axis == 2)) or
          (not self.lusevectorpotential and axis == 1))): continue
      # --- vpoisrzb converts J to A in place.
      a = self.source[axis,self.nxguardrho:-self.nxguardrho or None,
                           self.nyguardrho:-self.nyguardrho or None,
                           self.nzguardrho:-self.nzguardrho or None].copy()
      vpoisrzb(-1,a,kzsq,attz,filt,lr,lz,self.nx,self.nz,
               rfsmat,scrtch2,axis)
      self.potential[axis,self.nxguardphi:-self.nxguardphi or None,
                          self.nyguardrho:-self.nyguardrho or None,
                          self.nzguardphi:-self.nzguardphi or None] = a

      # --- This is needed since vpoisrzb doesn't have access to the guard cells.
      applyboundaryconditions3d(self.nx,self.ny,self.nz,
                                self.nxguardphi,self.nyguardphi,self.nzguardphi,
                                self.potential[axis,...],1,
                                array([1,0,1,1,2,2]),true,false)

    # --- Now take the curl of A to get B.
    getbfroma3d(self.potential,self.field,
                self.nxlocal,self.nylocal,self.nzlocal,
                self.nxguardphi,self.nyguardphi,self.nzguardphi,
                self.nxguarde,self.nyguarde,self.nzguarde,
                self.dx,self.dy,self.dz,self.xmminlocal,
                self.solvergeom==w3d.RZgeom,self.lusevectorpotential)

    # --- If using the analytic form of Btheta, calculate it here.
    if self.lanalyticbtheta:
      getanalyticbtheta(self.field,self.source,
                        self.nxlocal,self.nylocal,self.nzlocal,
                        self.nxguarde,self.nyguarde,self.nzguarde,
                        self.nxguardrho,self.nyguardrho,self.nzguardrho,
                        self.dx,self.xmminlocal)

    # --- Unscale the current density
    self.source[...] = self.source/(mu0*eps0)

# --- This can only be done after MagnetostaticMG is defined.
try:
  psyco.bind(MagnetostaticMG)
  psyco.bind(MagnetostaticFFT)
except NameError:
  pass


