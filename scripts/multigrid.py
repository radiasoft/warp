"""Class for doing complete multigrid field solve"""
# ToDo:
#  - modify setrhop to check if particles are within grid
#  - incorporate instances into the particle mover, so charge is deposited and
#    the E fields gather appropriately.
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
class MultiGrid(SubcycledPoissonSolver):
  
  __w3dinputs__ = ['iondensity','electrontemperature','plasmapotential',
                   'electrondensitymaxscale']
  __f3dinputs__ = ['gridmode','mgparam','downpasses','uppasses',
                   'mgmaxiters','mgtol','mgmaxlevels','mgform','mgverbose',
                   'lcndbndy','icndbndy','laddconductor'] 

  def __init__(self,lreducedpickle=1,**kw):
    kw['lreducedpickle'] = lreducedpickle
    SubcycledPoissonSolver.__init__(self,kwdict=kw)
    self.solvergeom = w3d.XYZgeom
    self.ncomponents = 1
    self.nxguard = 0
    self.nyguard = 0
    self.nzguard = 1

    # --- Kludge - make sure that the multigrid3df routines never sets up
    # --- any conductors.
    f3d.gridmode = 1

    # --- Save input parameters
    for name in MultiGrid.__w3dinputs__:
      if name not in self.__dict__:
        #self.__dict__[name] = kw.pop(name,getattr(w3d,name)) # Python2.3
        self.__dict__[name] = kw.get(name,getattr(w3d,name))
      if kw.has_key(name): del kw[name]
    for name in MultiGrid.__f3dinputs__:
      if name not in self.__dict__:
        #self.__dict__[name] = kw.pop(name,getattr(f3d,name)) # Python2.3
        self.__dict__[name] = kw.get(name,getattr(f3d,name))
      if kw.has_key(name): del kw[name]

    # --- If there are any remaning keyword arguments, raise an error.
    assert len(kw.keys()) == 0,"Bad keyword arguemnts %s"%kw.keys()

    # --- Create a conductor object, which by default is empty.
    self.conductors = ConductorType()
    self.conductorlist = []
    self.newconductorlist = []

    # --- Give these variables dummy initial values.
    self.mgiters = 0
    self.mgerror = 0.

    # --- At the start, assume that there are no bends. This is corrected
    # --- in the solve method when there are bends.
    self.linbend = false

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
    "Checks for and installs any new conductors before returning the object"
    # --- This method is needed since during a restore from a pickle, this
    # --- object may be restored before the conductors. This delays the
    # --- installation of the conductors until they are really needed,
    # --- which will only happen after the restoration is complete.
    for conductor in self.newconductorlist:
      self.installconductor(conductor)
    self.newconductorlist = []
    return self.conductors

  def getpdims(self):
    # --- Returns the dimensions of the arrays used by the particles

    # --- If there are any relativistic groups, them turn on the code
    # --- which uses the selfe array.
    if max(top.fselfb) > 0.:
      self.efetch = 3
      # --- Number of fields (E and B)
      nfields = 2
    else:
      # --- Number of fields (E only)
      nfields = 1

    if self.efetch == 3:
      return ((1+self.nxp,1+self.nyp,1+self.nzp),
              (3,1+self.nxp,1+self.nyp,1+self.nzp,nfields),
              (1+self.nxp,1+self.nyp,3+self.nzp))
    else:
      return ((1+self.nxp,1+self.nyp,1+self.nzp),
              (1+self.nxp,1+self.nyp,3+self.nzp))

  def getdims(self):
    # --- Returns the dimensions of the arrays used by the field solver
    return ((1+self.nx,1+self.ny,1+self.nz),
            (1+self.nx,1+self.ny,3+self.nz))

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
    wght = zeros((0,), 'd')
    self.setsourcepatposition(x,y,z,ux,uy,uz,gaminv,wght,q,w,zgrid)

  def setsourcepatposition(self,x,y,z,ux,uy,uz,gaminv,wght,q,w,zgrid):
    n  = len(x)
    if n == 0: return
    setrho3d(self.sourcep,n,x,y,z,zgrid,uz,q,w,top.depos,
             self.nxp,self.nyp,self.nzp,self.dx,self.dy,self.dz,
             self.xmminp,self.ymminp,self.zmminp,self.l2symtry,self.l4symtry,
             self.solvergeom==w3d.RZgeom)
    if self.lzerorhointerior:
      conductorobject = self.getconductorobject()
      cond_zerorhointerior(conductorobject.interior,self.nxp,self.nyp,self.nzp,
                           self.sourcep)

  def fetchfieldfrompositions(self,x,y,z,ex,ey,ez,bx,by,bz,pgroup=None):
    # --- Only sets the E field from the potential
    n = len(x)
    if n == 0: return
    sete3d(self.potentialp,self.fieldp,n,x,y,z,top.zgridprv,
           self.xmminp,self.ymminp,self.zmminp,
           self.dx,self.dy,self.dz,self.nxp,self.nyp,self.nzp,self.efetch,
           ex,ey,ez,self.l2symtry,self.l4symtry,self.solvergeom==w3d.RZgeom)
    if max(top.fselfb) > 0.:
      assert len(bx) == n,"The multigrid needs to be fixed so the B fields can be fetched with other than fetche3d"
      setb3d(self.fieldp[:,:,:,:,1],n,x,y,z,top.zgridprv,bx,by,bz,
             self.nxp,self.nyp,self.nzp,self.dx,self.dy,self.dz,
             self.xmminp,self.ymminp,self.zmminp,
             self.l2symtry,self.l4symtry,self.solvergeom==w3d.RZgeom)

  def fetchpotentialfrompositions(self,x,y,z,phi):
    n = len(x)
    getgrid3d(n,x,y,z,phi,self.nxp,self.nyp,self.nzp,self.potentialp[:,:,1:-1],
              self.xmminp,self.xmmaxp,self.ymminp,self.ymmaxp,self.zmminp,self.zmmaxp,
              self.l2symtry,self.l4symtry)

  def setsourceforfieldsolve(self,*args):
    SubcycledPoissonSolver.setsourceforfieldsolve(self,*args)
    if self.lparallel:
      SubcycledPoissonSolver.setsourcepforparticles(self,*args)
      setrhoforfieldsolve3d(self.nx,self.ny,self.nz,self.source,
                            self.nxp,self.nyp,self.nzp,self.sourcep,self.nzpguard,
                            self.my_index,self.nslaves,self.izpslave,self.nzpslave,
                            self.izfsslave,self.nzfsslave)

  def getpotentialpforparticles(self,*args):
    if not self.lparallel:
      SubcycledPoissonSolver.getpotentialpforparticles(self,*args)
    else:
      self.setpotentialpforparticles(*args)
      getphipforparticles3d(1,self.nx,self.ny,self.nz,self.potential,
                            self.nxp,self.nyp,self.nzp,self.potentialp,1)
    if self.efetch == 3:
      self.setpotentialpforparticles(*args)
      self.setfieldpforparticles(*args)
      indts = args[1]
      iselfb = args[2]
      # --- If this is the first group, set make sure that fieldp gets
      # --- zeroed out. Otherwise, the data in fieldp is accumulated.
      # --- This coding relies on the fact that fieldsolver does the
      # --- loops in descending order.
      tmpnsndts = getnsndtsforsubcycling()
      lzero = ((indts == tmpnsndts-1) and (iselfb == top.nsselfb-1))
      if lzero:
        tfieldp = transpose(self.fieldp)
        tfieldp[...] = 0.
      self.getselfe(recalculate=1)
      # --- If top.fslefb(iselfb) > 0, then calculate and include the
      # --- approximate correction terms A and dA/dt.
      self.getselfb(self.fieldp,top.fselfb[iselfb],self.potentialp)
      self.adddadttoe(self.fieldp,top.fselfb[iselfb],self.potentialp)

  def makesourceperiodic(self):
    if self.pbounds[0] == 2 or self.pbounds[1] == 2:
      self.source[0,:,:] = self.source[0,:,:] + self.source[-1,:,:]
      self.source[-1,:,:] = self.source[0,:,:]
    if self.pbounds[2] == 2 or self.pbounds[3] == 2:
      self.source[:,0,:] = self.source[:,0,:] + self.source[:,-1,:]
      self.source[:,-1,:] = self.source[:,0,:]
    if self.pbounds[0] == 1 and not self.l4symtry:
       self.source[0,:,:] = 2.*self.source[0,:,:]
    if self.pbounds[1] == 1: self.source[-1,:,:] = 2.*self.source[-1,:,:]
    if self.pbounds[2] == 1 and not (self.l2symtry or self.l4symtry):
       self.source[:,0,:] = 2.*self.source[:,0,:]
    if self.pbounds[3] == 1: self.source[:,-1,:] = 2.*self.source[:,-1,:]
    if self.pbounds[4] == 2 or self.pbounds[5] == 2:
      if self.lparallel:
        self.makesourceperiodic_parallel()
      else:
        self.source[:,:,0] = self.source[:,:,0] + self.source[:,:,-1]
        self.source[:,:,-1] = self.source[:,:,0]
    if self.pbounds[4] == 1: self.source[:,:,0] = 2.*self.source[:,:,0]
    if self.pbounds[5] == 1: self.source[:,:,-1] = 2.*self.source[:,:,-1]

  def makesourceperiodic_parallel(self):
    tag = 70
    if self.my_index == self.nslaves-1:
      request = mpi.isend(self.source[:,:,self.nz],0,tag)
      self.source[:,:,self.nz],status = mpi.recv(0,tag)
    elif self.my_index == 0:
      sourcetemp,status = mpi.recv(self.nslaves-1,tag)
      self.source[:,:,0] = self.source[:,:,0] + sourcetemp
      request = mpi.isend(self.source[:,:,0],self.nslaves-1,tag)
    if self.my_index == 0 or self.my_index == self.nslaves-1:
      status = request.wait()

  def getselfe(self,recalculate=0):
    if type(self.fieldp) != ArrayType:
      # --- This should only ever be done by an external routine, such as
      # --- a plotting function.
      self.fieldp = fzeros((3,1+self.nxp,1+self.nyp,1+self.nzp,1),'d')
    if recalculate:
      getselfe3d(self.potentialp,self.nxp,self.nyp,self.nzp,
                 self.fieldp[:,:,:,:,0],self.nxp,self.nyp,self.nzp,
                 self.dx,self.dy,self.dz,
                 self.bounds[0],self.bounds[1],self.bounds[2],self.bounds[3])
    return self.fieldp

  def getselfb(self,fieldp,fselfb,potentialp):
    Az = (fselfb/clight**2)*potentialp[:,:,1:-1]
    fieldp[0,:,1:-1,:,1] += (Az[:,2:,:] - Az[:,:-2,:])/(2.*self.dy)
    fieldp[1,1:-1,:,:,1] -= (Az[2:,:,:] - Az[:-2,:,:])/(2.*self.dx)
    if self.bounds[2] == 1 or self.l2symtry or self.l4symtry:
      pass
    elif self.bounds[2] == 0:
      fieldp[0,:,0,:,1] += (Az[:,1,:] - Az[:,0,:])/(self.dy)
    elif self.bounds[2] == 2:
      fieldp[0,:,0,:,1] += (Az[:,1,:] - Az[:,-2,:])/(2.*self.dy)
    if self.bounds[3] == 0:
      fieldp[0,:,-1,:,1] += (Az[:,-1,:] - Az[:,-2,:])/(self.dy)
    elif self.bounds[3] == 2:
      fieldp[0,:,-1,:,1] += (Az[:,1,:] - Az[:,-2,:])/(2.*self.dy)
    if self.bounds[0] == 1 or self.l4symtry:
      pass
    elif self.bounds[0] == 0:
      fieldp[1,0,:,:,1] -= (Az[1,:,:] - Az[0,:,:])/(self.dx)
    elif self.bounds[0] == 2:
      fieldp[1,0,:,:,1] -= (Az[1,:,:] - Az[-2,:,:])/(2.*self.dx)
    if self.bounds[1] == 0:
      fieldp[1,-1,:,:,1] -= (Az[-1,:,:] - Az[-2,:,:])/(self.dx)
    elif self.bounds[1] == 2:
      fieldp[1,-1,:,:,1] -= (Az[1,:,:] - Az[-2,:,:])/(2.*self.dx)
    
  def adddadttoe(self,fieldp,fselfb,potentialp):
    """Ez = -dA/dt = -beta**2 dphi/dz"""
    Ez = (fselfb/clight)**2*(potentialp[:,:,2:]-potentialp[:,:,:-2])/(2.*self.dz)
    fieldp[2,:,:,:,0] += Ez

  def installconductor(self,conductor,
                            xmin=None,xmax=None,
                            ymin=None,ymax=None,
                            zmin=None,zmax=None,
                            dfill=top.largepos):
    if conductor in self.conductorlist: return
    self.conductorlist.append(conductor)
    # --- Note that this uses the conductorobject directly since it is called
    # --- by getconductorobject.
    installconductors(conductor,xmin,xmax,ymin,ymax,zmin,zmax,dfill,
                      top.zbeam,
                      self.nx,self.ny,self.nz,self.nzfull,
                      self.xmmin,self.xmmax,self.ymmin,self.ymmax,
                      self.zmminglobal,self.zmmaxglobal,self.l2symtry,self.l4symtry,
                      solvergeom=self.solvergeom,conductors=self.conductors,
                      my_index=self.my_index,nslaves=self.nslaves,
                      izfsslave=self.izfsslave,nzfsslave=self.nzfsslave)

  def hasconductors(self):
    conductorobject = self.getconductorobject()
    return (conductorobject.interior.n > 0 or
            conductorobject.evensubgrid.n > 0 or
            conductorobject.oddsubgrid.n > 0)

  def clearconductors(self):
    "Should the conductorlist and newconductorlist be cleared also?"
    # --- This uses the conductor object directly since there is not point in
    # --- installing any new conductors (which would be done by
    # --- getconductorobject).
    self.conductors.interior.n = 0
    self.conductors.evensubgrid.n = 0
    self.conductors.oddsubgrid.n = 0

  def find_mgparam(self,lsavephi=false,resetpasses=1):
    # --- This is a temporary kludge, the same as is done in genericpf
    self.phi = self.potential
    find_mgparam(lsavephi=lsavephi,resetpasses=resetpasses,
                 solver=self,pkg3d=self)

  def dosolve(self,iwhich=0,*args):
    # --- set for longitudinal relativistic contraction 
    iselfb = args[2]
    beta = top.pgroup.fselfb[iselfb]/clight
    zfact = 1./sqrt((1.-beta)*(1.+beta))

    # --- This is only done for convenience.
    self.phi = self.potential
    self.rho = self.source

    # --- Setup data for bends.
    rstar = fzeros(3+self.nz,'d')
    if top.bends:

      # --- This commented out code does the same thing as the line below
      # --- setting linbend but is a bit more complicated. It is preserved
      # --- in case of some unforeseen problem with the code below.
      #ii = (top.cbendzs <= self.zmmax+top.zgrid and
      #                     self.zmmin+top.zgrid <= top.cbendze)
      #self.linbend = sometrue(ii)

      setrstar(rstar,self.nz,self.dz,self.zmmin,top.zgrid)
      self.linbend = min(rstar) < largepos

    if self.izfsslave is None: self.izfsslave = top.izfsslave
    if self.nzfsslave is None: self.nzfsslave = top.nzfsslave
    mgiters = zeros(1)
    mgerror = zeros(1,'d')
    conductorobject = self.getconductorobject()
    if self.electrontemperature == 0:
      multigrid3dsolve(iwhich,self.nx,self.ny,self.nz,self.nzfull,
                       self.dx,self.dy,self.dz*zfact,self.potential,self.source,
                       rstar,self.linbend,self.bounds,
                       self.xmmin,self.ymmin,self.zmmin*zfact,
                       self.zmminglobal*zfact,top.zbeam*zfact,top.zgrid*zfact,
                       self.mgparam,self.mgform,mgiters,self.mgmaxiters,
                       self.mgmaxlevels,mgerror,self.mgtol,self.mgverbose,
                       self.downpasses,self.uppasses,
                       self.lcndbndy,self.laddconductor,self.icndbndy,
                       self.lbuildquads,self.gridmode,conductorobject,
                       self.my_index,self.nslaves,self.izfsslave,self.nzfsslave)
    else:
      multigridbe3dsolve(iwhich,self.nx,self.ny,self.nz,self.nzfull,
                         self.dx,self.dy,self.dz*zfact,self.potential,self.source,
                         star,self.linbend,self.bounds,
                         self.xmmin,self.ymmin,self.zmmin*zfact,
                         self.zmminglobal*zfact,top.zbeam*zfact,top.zgrid*zfact,
                         self.mgparam,mgiters,self.mgmaxiters,
                         self.mgmaxlevels,mgerror,self.mgtol,self.mgverbose,
                         self.downpasses,self.uppasses,
                         self.lcndbndy,self.laddconductor,self.icndbndy,
                         self.lbuildquads,self.gridmode,conductorobject,
                         self.my_index,self.nslaves,self.izfsslave,
                         self.nzfsslave,
                         self.iondensity,self.electrontemperature,
                         self.plasmapotential,self.electrondensitymaxscale)
    self.mgiters = mgiters[0]
    self.mgerror = mgerror[0]

  ##########################################################################
  # Define the basic plot commands
  def genericpf(self,kw,pffunc):
    kw['conductors'] = self.getconductorobject()
    kw['solver'] = self
    pffunc(**kw)
  def pfxy(self,**kw): self.genericpf(kw,pfxy)
  def pfzx(self,**kw): self.genericpf(kw,pfzx)
  def pfzy(self,**kw): self.genericpf(kw,pfzy)
  def pfxyg(self,**kw): self.genericpf(kw,pfxyg)
  def pfzxg(self,**kw): self.genericpf(kw,pfzxg)
  def pfzyg(self,**kw): self.genericpf(kw,pfzyg)








  #===========================================================================
  def solve1(self,iwhich=0):
    # --- No initialization needed
    if iwhich == 1: return

    # --- Create temp arrays
    phisave = fzeros(shape(self.phi),'d')
    bendx = fzeros(((self.nx+1)*(self.ny+1)),'d')

    # --- Initialize temporaries
    nxy    = (self.nx+1)*(self.ny+1)
    nxyz   = (self.nx+1)*(self.ny+1)*(self.nz+1)
    dxsqi  = 1./self.dx**2
    dysqi  = 1./self.dy**2
    dzsqi  = 1./self.dz**2
    reps0c = self.mgparam/(eps0*2.*(dxsqi+dysqi+dzsqi))
    rdel   = dzsqi/(dxsqi + dysqi + dzsqi)

    conductorobject = self.getconductorobject()
    checkconductors(self.nx,self.ny,self.nz,self.nzfull,
                    self.dx,self.dy,self.dz,conductorobject,
                    top.my_index,top.nslaves,top.izfsslave,top.nzfsslave)

    # --- Preset rho to increase performance (reducing the number of
    # --- multiplies in the main SOR sweep loop).
    if not self.linbend:
      # --- Do the operation in place (to avoid temp arrays)
      multiply(self.rho,reps0c,self.rho)
    else:
      raise "Bends not yet supported"

    # --- If using residual correction form, need to save the original rho.
    # --- Also setup parallel arrays.
    if self.mgform == 2:
      rhosave = self.rho + 0.
      res = fzeros(shape(self.phi),'d')
      localbounds = bounds.copy()

    #ifdef MPIPARALLEL
    #  lparity = 0
    #  rparity = 0
    #  mggetexchangepes(nslaves,izfsslave,nzfsslave,my_index,
    #                   bounds,nzfull,
    #                   lparity,rparity,
    #                   whosendingleft,izsendingleft,
    #                   whosendingright,izsendingright)
    #  if (izfsslave(my_index) > 0) localbounds[4] = -1
    #  if (izfsslave(my_index)+nz < nzfull) localbounds[5] = -1
    #endif

    #   --- Main multigrid v-cycle loop. Calculate error each iteration since
    #   --- very few iterations are done.
    self.mgiters = 0
    self.mgerror = 2.*self.mgtol + 1.
    while (self.mgerror > self.mgtol and self.mgiters < self.mgmaxiters):
      self.mgiters = self.mgiters + 1

      # --- Save current value of phi
      phisave[:,:,:] = self.phi + 0.

      # --- If using residual correction form, calculate the residual and
      # --- copy it into rhosave, zero phisave (the initial error).
      # --- In the calls to cond_potmg and residual, the last argument
      # --- is true, telling the routines to use the actual value of
      # --- voltages rather than zero as is done otherwise for residual
      # --- correction form since it is operating on the error.
      if self.mgform == 2:
        cond_potmg(conductorobject.interior,
                   self.nx,self.ny,self.nz,phisave,0,self.mgform,true)
        residual(self.nx,self.ny,self.nz,self.nzfull,dxsqi,dysqi,dzsqi,
                 phisave,rhosave,res,0,localbounds,
                 self.mgparam,self.mgform,true,
                 self.lcndbndy,self.icndbndy,conductorobject)
    #ifdef MPIPARALLEL
    #  mgexchange_phi(nx,ny,nz,nzfull,res,localbounds,-1,
    #                 my_index,nslaves,izfsslave,nzfsslave,
    #                 whosendingleft,izsendingleft,
    #                 whosendingright,izsendingright)
    #  mgexchange_phiperiodic(nx,ny,nz,nzfull,res,localbounds,0,
    #                         my_index,nslaves,izfsslave,
    #                         whosendingleft,whosendingright)
    #endif
        self.rho[:,:,:] = res[:,:,1:-1]
        self.phi[:,:,:] = 0.

      # --- Do one vcycle.
      self.vcycle(0,self.nx,self.ny,self.nz,self.nzfull,
                  self.dx,self.dy,self.dz,self.phi,self.rho,
                  self.rstar,self.linbend,bendx,self.bounds,
                  self.mgparam,self.mgform,self.mgmaxlevels,
                  self.downpasses,self.uppasses,self.lcndbndy,
                  self.icndbndy,conductorobject)

      # --- If using residual correction form, add the resulting error to phi.
      if self.mgform == 2: add(self.phi,phisave,self.phi)

      # --- When using residual correction form, the other planes do need
      # --- to be set when using other than Dirichlet boundaries since
      # --- those planes are only set with the error of phi.
      if self.mgform == 2:
        if localbounds[4] == 1: self.phi[:,:,0] = self.phi[:,:,2]
        if localbounds[5] == 1: self.phi[:,:,-1] = self.phi[:,:,-3]
      #ifndef MPIPARALLEL
        if localbounds[4] == 2: self.phi[:,:,0] = self.phi[:,:,-3]
        if localbounds[5] == 2: self.phi[:,:,-1] = self.phi[:,:,2]
      #else
      # mgexchange_phi(nx,ny,nz,nzfull,phi,localbounds,0,
      #                my_index,nslaves,izfsslave,nzfsslave,
      #                whosendingleft,izsendingleft,
      #                whosendingright,izsendingright)
      # mgexchange_phi(nx,ny,nz,nzfull,phi,localbounds,-1,
      #                my_index,nslaves,izfsslave,nzfsslave,
      #                whosendingleft,izsendingleft,
      #                whosendingright,izsendingright)
      #endif

      # --- Calculate the change in phi.
      subtract(phisave,self.phi,phisave)
      absolute(phisave,phisave)
      self.mgerror = MA.maximum(phisave)

    #ifdef MPIPARALLEL
    #     --- calculate global sorerror
    #     call parallelmaxrealarray(self.mgerror,1)

    # --- For Dirichlet boundary conditions, copy data into guard planes
    # --- For other boundary conditions, the guard planes are used during
    # --- the solve are so are already set.
    if (self.bounds[4] == 0): self.phi[:,:,0] = self.phi[:,:,1]
    if (self.bounds[5] == 0): self.phi[:,:,-1] = self.phi[:,:,-2]

    # --- Make a print out.
    if (self.mgerror > self.mgtol):
      print "Multigrid: Maximum number of iterations reached"
    print ("Multigrid: Error converged to %11.3e in %4d v-cycles"%
           (self.mgerror,self.mgiters))

    # --- If using residual correction form, restore saved rho
    if self.mgform == 2: self.rho[:,:,:] = rhosave

    # --- Restore rho
    if (not self.linbend):
      multiply(self.rho,1./reps0c,self.rho)

  #===========================================================================
  def vcycle(self,mglevel,nx,ny,nz,nzfull,dx,dy,dz,
                  phi,rho,rstar,linbend,bendx,bounds,mgparam,mgform,
                  mgmaxlevels,downpasses,uppasses,lcndbndy,icndbndy,conductors):
   
    res = fzeros(shape(phi),'d')

    dxsqi = 1./dx**2
    dysqi = 1./dy**2
    dzsqi = 1./dz**2

    localbounds = bounds.copy()

    #ifdef MPIPARALLEL
    #lparityall = 0
    #rparityall = 0
    #mggetexchangepes(nslaves,izfsslave,nzfsslave,my_index,
    #                 bounds,nzfull,
    #                 lparityall,rparityall,
    #                 whosendingleft,izsendingleft,
    #                 whosendingright,izsendingright)
    #if (izfsslave(my_index) > 0) localbounds[4] = -1
    #if (izfsslave(my_index)+nz < nzfull) localbounds[5] = -1
    #endif

    for i in range(downpasses):
      self.sorpass3d(mglevel,nx,ny,nz,nzfull,phi,rho,rstar,
                     dxsqi,dysqi,dzsqi,linbend,bendx,
                     localbounds,mgparam,mgform,lcndbndy,icndbndy,conductors)

    # --- Check if this is the finest level. If so, then don't do any further
    # --- coarsening. This is the same check that is done in getmglevels.
    # --- If grid is not at its coarsest level in any of the axis or and
    # --- all dimensions are even, continue the coarsening.
    if ((nx%4) == 0 and (ny%4) == 0 and (nzfull%4) == 0 and
        mglevel < mgmaxlevels):

      # --- Get the residual on the current grid.
      residual(nx,ny,nz,nzfull,dxsqi,dysqi,dzsqi,phi,rho,res,
               mglevel,localbounds,mgparam,mgform,false,
               lcndbndy,icndbndy,conductors)
      #ifdef MPIPARALLEL
      #mgexchange_phi(nx,ny,nz,nzfull,res,localbounds,-1,
      #                  my_index,nslaves,izfsslave,nzfsslave,
      #                  whosendingleft,izsendingleft,
      #                  whosendingright,izsendingright)
      #mgexchange_phiperiodic(nx,ny,nz,nzfull,res,localbounds,0,
      #                       my_index,nslaves,izfsslave,
      #                       whosendingleft,whosendingright)
      #endif

      # --- If dz > 4/3 dx then only coarsen transversely, otherwise coarsen
      # --- all axis.  This is the same check that is done in getmglevels.
      # --- dz > 4/3 dx <=> (9/16) / dx^2 < 1 / dz^2
      partialcoarsening = (dz > 4./3.*dx)
      #ifdef MPIPARALLEL
      # --- This must be a global operation since, due to roundoff, each
      # --- processor can get a different value if dz == 4./3.*dx.
      #parallellor(partialcoarsening)
      #endif
      if partialcoarsening:

        # --- Allocate new work space
        phi2 = fzeros((1+nx/2,1+ny/2,2+nz+1),'d')
        rho2 = fzeros((1+nx/2,1+ny/2,1+nz),'d')

        # --- Ratio of old to new constant needed to scale the residual for
        # --- the restriction.
        ff = (dxsqi+dysqi+dzsqi)/(dxsqi*0.25 + dysqi*0.25 + dzsqi)
        restrict2d(nx,ny,nz,nzfull,res,rho2,ff,localbounds)

        # --- Continue at the next coarsest level.
        self.vcycle(mglevel+1,nx/2,ny/2,nz,nzfull,
                    dx*2,dy*2,dz,phi2,rho2,rstar,linbend,bendx,bounds,
                    mgparam,mgform,mgmaxlevels,downpasses,uppasses,
                    lcndbndy,icndbndy,conductors)

        # --- Add in resulting error.
        expand2d(nx/2,ny/2,nz,nzfull,phi2,phi,localbounds)

      else:

        localbounds2 = bounds.copy()

        #ifdef MPIPARALLEL
        # --- Find domains in coarser grid
        # call mgdividenz(nslaves,izfsslave,nzfsslave,izfsslave2,nzfsslave2,
        #                 nzfull)
        # --- Set new value of nz
        # nznew = nzfsslave2(my_index)
        # --- Difference between starts and ends of coarse and fine grids.
        # --- Should only be in the range 0-2.
        # lparityall = izfsslave - 2*izfsslave2
        # rparityall = 2*(izfsslave2 + nzfsslave2) - (izfsslave + nzfsslave)
        # --- Note that the lparityall and rparityall can only be used in
        # --- MPIPARALLEL sections since they will be unallocated in the
        # --- serial code. So, separate scalars are used in code which is
        # --- used in the serial version.
        # lparity = lparityall(my_index)
        # rparity = rparityall(my_index)
        # --- Get processor with which to exchange data on coarse grid
        # call mggetexchangepes(nslaves,izfsslave2,nzfsslave2,my_index,
        #                       globalb0,globalbnz,nzfull/2,
        #                       lparityall,rparityall,
        #                       whosendingleft2,izsendingleft2,
        #                       whosendingright2,izsendingright2)
        # if (izfsslave2(my_index) > 0) localbounds2[4] = -1
        # if (izfsslave2(my_index) + nznew < nzfull/2) localbounds2[5] = -1
        #else
        nznew = nz/2
        lparity = 0
        rparity = 0
        #endif

        # --- Alloate new work space
        phi2 = fzeros((1+nx/2,1+ny/2,2+nznew+1),'d')
        rho2 = fzeros((1+nx/2,1+ny/2,1+nznew),'d')

        # --- Restriction - note that scaling factor for residual is always
        # --- 4 for full-coarsening and is compiled into the restriction
        # --- routine.
        restrict3d(nx,ny,nz,nznew,nzfull,res,rho2,localbounds2,localbounds,
                   lparity,rparity)
        #ifdef MPIPARALLEL
        # mgexchange_phi(nx/2,ny/2,nznew,nzfull/2,rho2plusmorespace,
        #                localbounds2,0,
        #                my_index,nslaves,izfsslave2,nzfsslave2,
        #                whosendingleft2,izsendingleft2,
        #                whosendingright2,izsendingright2)
        #endif

        # --- Continue at the next coarsest level.
        self.vcycle(mglevel+1,nx/2,ny/2,nznew,nzfull/2,
                    dx*2,dy*2,dz*2,phi2,rho2,rstar,linbend,bendx,bounds,
                    mgparam,mgform,mgmaxlevels,downpasses,uppasses,
                    lcndbndy,icndbndy,conductors)

        # --- Add in resulting error.
        expand3d(nx/2,ny/2,nznew,nz,nzfull/2,phi2,phi,localbounds,
                 lparity,rparity)
        #ifdef MPIPARALLEL
        # mgexchange_phiperiodic(nx,ny,nz,nzfull,phi,
        #                        localbounds,1,
        #                        my_index,nslaves,izfsslave,
        #                        whosendingleft,whosendingright)
        #endif

    # --- Do final SOR passes.
    for i in range(uppasses):
      self.sorpass3d(mglevel,nx,ny,nz,nzfull,phi,rho,rstar,
                     dxsqi,dysqi,dzsqi,linbend,bendx,localbounds,
                     mgparam,mgform,lcndbndy,icndbndy,conductors)

  #===========================================================================
  def sorpass3d(self,mglevel,nx,ny,nz,nzfull,phi,rho,rstar,
                     rdx2,rdy2,rdz2,linbend,bendx,bounds,mgparam,mgform,
                     lcndbndy,icndbndy,conductors):

    # --- Put desired potential onto conductors in phi array.
    cond_potmg(conductors.interior,nx,ny,nz,phi,mglevel,mgform,false)

    # --- Set starting and ending parity.
    #ifdef MPIPARALLEL
    # s_parity = mod(izfsslave(my_index),2)
    # e_parity = mod(s_parity+1,2)
    #else
    s_parity = 0
    e_parity = 1
    #endif

    # --- do loop to cover even and odd points
    for parity in [s_parity,e_parity]:

      sorhalfpass3d(parity,mglevel,nx,ny,nz,nzfull,phi,rho,rstar,
                    rdx2,rdy2,rdz2,linbend,bendx,bounds,mgparam,mgform,
                    lcndbndy,icndbndy,conductors)

    #ifndef MPIPARALLEL
      if (bounds[4] == 2): phi[:,:,0] = phi[:,:,-3]
      if (bounds[5] == 2): phi[:,:,-2:] = phi[:,:,1:3]
    #else
    # mgexchange_phi(nx,ny,nz,nzfull,phi,bounds,0,
    #                my_index,nslaves,izfsslave,nzfsslave,
    #                whosendingleft,izsendingleft,
    #                whosendingright,izsendingright)
    # mgexchange_phiperiodic(nx,ny,nz,nzfull,phi,bounds,1,
    #                        my_index,nslaves,izfsslave,
    #                        whosendingleft,whosendingright)
    #endif

    #ifdef MPIPARALLEL
    # --- Exchange phi in the z guard planes
    #mgexchange_phi(nx,ny,nz,nzfull,phi,bounds,-1,
    #               my_index,nslaves,izfsslave,nzfsslave,
    #               whosendingleft,izsendingleft,
    #               whosendingright,izsendingright)
    #endif

  def getresidual(self):
    res = zeros(shape(self.phi),'d')
    dxsqi  = 1./self.dx**2
    dysqi  = 1./self.dy**2
    dzsqi  = 1./self.dz**2
    reps0c = self.mgparam/(eps0*2.*(dxsqi+dysqi+dzsqi))
    rho = self.rho*reps0c
    residual(self.nx,self.ny,self.nz,self.nzfull,dxsqi,dysqi,dzsqi,
             self.phi,rho,res,0,self.bounds,self.mgparam,self.mgform,false,
             self.lcndbndy,self.icndbndy,self.conductors)
    return res

# --- This can only be done after MultiGrid is defined.
try:
  psyco.bind(MultiGrid)
except NameError:
  pass

