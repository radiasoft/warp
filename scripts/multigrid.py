"""Class for doing complete multigrid field solve"""
# ToDo:
#  - modify setrhop to check if particles are within grid
#  - incorporate instances into the particle mover, so charge is deposited and
#    the E fields gather appropriately.
from warp import *
from fieldsolver import SubcycledPoissonSolver
from generateconductors import installconductors
from find_mgparam import find_mgparam

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
    self.nxguard = 1
    self.nyguard = 1
    self.nzguard = 1

    # --- Kludge - make sure that the multigrid3df routines never sets up
    # --- any conductors.
    f3d.gridmode = 1

    # --- Save input parameters
    self.processdefaultsfrompackage(MultiGrid.__w3dinputs__,w3d,kw)
    self.processdefaultsfrompackage(MultiGrid.__f3dinputs__,f3d,kw)

    # --- If there are any remaning keyword arguments, raise an error.
    assert len(kw.keys()) == 0,"Bad keyword arguemnts %s"%kw.keys()

    self.initializeconductors()

    # --- Give these variables dummy initial values.
    self.mgiters = 0
    self.mgerror = 0.

    # --- At the start, assume that there are no bends. This is corrected
    # --- in the solve method when there are bends.
    self.linbend = false

    # --- Turn of build quads option
    self.lbuildquads = false

  def initializeconductors(self):
    # --- Create the attributes for holding information about conductors
    # --- and conductor objects.
    # --- Note that a conductor object will be created for each value of
    # --- fselfb. This is needed since fselfb effects how the coarsening
    # --- is done, and different conductor data sets are needed for
    # --- different coarsenings.

    # --- This stores the ConductorType objects. Note that the objects are
    # --- not actually created until getconductorobject is called.
    self.conductorobjects = {}

    # --- This stores the conductors that have been installed in each
    # --- of the conductor objects.
    self.installedconductorlists = {}

    # --- This is a list of conductors that have been added.
    # --- New conductors are not actually installed until the data is needed,
    # --- when getconductorobject is called.
    # --- Each element of this list contains all of the input to the
    # --- installconductor method.
    self.conductordatalist = []

  def __getstate__(self):
    dict = SubcycledPoissonSolver.__getstate__(self)
    if self.lreducedpickle:

      # --- Write out an empy conductorobjects since it can be big. Also,
      # --- write out an empty list of conductors so they will all be
      # --- reinstalled upon restoration.
      dict['conductorobjects'] = {}
      dict['installedconductorlists'] = {}

      if 'rho' in dict: del dict['rho']
      if 'phi' in dict: del dict['phi']

    return dict

  def __setstate__(self,dict):
    SubcycledPoissonSolver.__setstate__(self,dict)

    # --- Check if an old file is being restored
    if 'conductorobjects' not in self.__dict__:

      # --- Create the appropriate attributes that are now needed.
      # --- This is not the best thing, since is replicates code in
      # --- the __init__
      self.conductorobjects = {}
      self.installedconductorlists = {}
      self.conductordatalist = []

      # --- Get the list of conductors from old formats
      if 'newconductorlist' in self.__dict__:
        conductorlist = self.newconductorlist
        del self.newconductorlist
      elif 'conductorlist' in self.__dict__:
        conductorlist = self.conductorlist
        del self.conductorlist
      else:
        conductorlist = []

      for conductor in conductorlist:
        self.installconductor(conductor)

  def getconductorobject(self,fselfb=0.):
    "Checks for and installs any conductors not yet installed before returning the object"
    # --- This is the routine that does the creation of the ConductorType
    # --- objects if needed and ensures that all conductors are installed
    # --- into it.

    # --- This method is needed during a restore from a pickle, since this
    # --- object may be restored before the conductors. This delays the
    # --- installation of the conductors until they are really needed.

    # --- There is a special case, fselfb='p', which refers to the conductor
    # --- object that has the data generated relative to the particle domain,
    # --- which can be different from the field domain, especially in parallel.
    if fselfb == 'p':
      # --- In serial, just use a reference to the conductor object for the
      # --- first iselfb group.
      if not lparallel and 'p' not in self.conductorobjects:
        self.conductorobjects['p'] = self.conductorobjects[top.fselfb[0]]
        self.installedconductorlists['p'] = self.installedconductorlists[top.fselfb[0]]
      # --- In parallel, a whole new instance is created (using the
      # --- setdefaults below).
      # --- Check to make sure that the grid the conductor uses is consistent
      # --- with the particle grid. This is needed so that the conductor
      # --- data is updated when particle load balancing is done. If the
      # --- data is not consistent, delete the conductor object so that
      # --- everything is reinstalled.
      try:
        conductorobject = self.conductorobjects['p']
        if (conductorobject.leveliz[0] != self.izpslave[self.my_index] or
            conductorobject.levelnz[0] != self.nzpslave[self.my_index]):
          del self.conductorobjects['p']
          del self.installedconductorlists['p']
      except KeyError:
        # --- 'p' object has not yet been created anyway, so do nothing.
        pass

    conductorobject = self.conductorobjects.setdefault(fselfb,ConductorType())
    installedconductorlist = self.installedconductorlists.setdefault(fselfb,[])

    # --- Now, make sure that the conductors are installed into the object.
    # --- This may be somewhat inefficient, since it loops over all of the
    # --- conductors everytime. This makes the code more robust, though, since
    # --- it ensures that all conductors will be properly installed into
    # --- the conductor object.
    for conductordata in self.conductordatalist:
      self._installconductor(conductorobject,installedconductorlist,
                             conductordata,fselfb)
 
    # --- Return the desired conductor object
    return conductorobject

  def setconductorvoltage(self,voltage,condid=0,discrete=false,
                          setvinject=false):
    'calls setconductorvoltage'
    # --- Loop over all of the selfb groups to that all conductor objects
    # --- are handled.
    for iselfb in range(top.nsselfb):
      conductorobject = self.getconductorobject(top.fselfb[iselfb])
      setconductorvoltage(voltage,condid,discrete,setvinject,
                          conductors=conductorobject)

  def getpdims(self):
    # --- Returns the dimensions of the arrays used by the particles

    # --- If there are any relativistic groups, then turn on the code
    # --- which uses the selfe array.
    if max(abs(top.fselfb)) > 0.:
      # --- This is probably redundant, but it shouldn't hurt.
      # --- This forces all species to use the precalculated E field
      # --- if any have the B correction.
      top.efetch = 3
      # --- Number of fields (E and B)
      nfields = 2
    else:
      # --- Number of fields (E only)
      nfields = 1

    if sometrue(top.efetch == 3):
      return ((1+self.nxp,1+self.nyp,1+self.nzp),
              (3,1+self.nxp,1+self.nyp,1+self.nzp,nfields),
              (1+self.nxp+2*self.nxguard,1+self.nyp+2*self.nyguard,1+self.nzp+2*self.nzguard))
    else:
      return ((1+self.nxp,1+self.nyp,1+self.nzp),
              (1+self.nxp+2*self.nxguard,1+self.nyp+2*self.nyguard,1+self.nzp+2*self.nzguard))

  def getdims(self):
    # --- Returns the dimensions of the arrays used by the field solver
    return ((1+self.nx,1+self.ny,1+self.nzlocal),
            (1+self.nx+2*self.nxguard,1+self.ny+2*self.nyguard,1+self.nzlocal+2*self.nzguard))

  def getrho(self):
    return self.source

  def getphi(self):
    'Returns the phi array without the guard cells'
    ix1 = self.nxguard
    if ix1 == 0: ix1 = None
    ix2 = -self.nxguard
    if ix2 == 0: ix2 = None
    ix = slice(ix1,ix2)
    iy1 = self.nyguard
    if iy1 == 0: iy1 = None
    iy2 = -self.nyguard
    if iy2 == 0: iy2 = None
    iy = slice(iy1,iy2)
    iz1 = self.nzguard
    if iz1 == 0: iz1 = None
    iz2 = -self.nzguard
    if iz2 == 0: iz2 = None
    iz = slice(iz1,iz2)
    return self.potential[ix,iy,iz]

  def getfield(self):
    return self.field

  def loadrho(self,lzero=None,pgroups=None,**kw):
    SubcycledPoissonSolver.loadsource(self,lzero,pgroups,**kw)

  def fetche(self,*args,**kw):
    SubcycledPoissonSolver.fetchfield(self,*args,**kw)

  def fetchphi(self,*args,**kw):
    SubcycledPoissonSolver.fetchpotential(self,*args,**kw)

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
    w  = pgroup.sw[js]*pgroup.dtscale[js]
    wght = zeros((0,), 'd')
    if top.wpid==0:
      wfact = zeros((0,), 'd')
    else:
      wfact = pgroup.pid[i:i+n,top.wpid-1]
    self.setsourcepatposition(x,y,z,ux,uy,uz,gaminv,wfact,wght,q,w,zgrid)

  def setsourcepatposition(self,x,y,z,ux,uy,uz,gaminv,wfact,wght,q,w,zgrid):
    n = len(x)
    if n == 0: return
    if isinstance(self.sourcep,FloatType): return
    if top.wpid==0:
      setrho3d(self.sourcep,n,x,y,z,zgrid,uz,q,w,top.depos,
               self.nxp,self.nyp,self.nzp,self.dx,self.dy,self.dz,
               self.xmminp,self.ymminp,self.zmminp,self.l2symtry,self.l4symtry,
               self.solvergeom==w3d.RZgeom)
    else:
      setrho3dw(self.sourcep,n,x,y,z,zgrid,uz,wfact,q,w,top.depos,
                self.nxp,self.nyp,self.nzp,self.dx,self.dy,self.dz,
                self.xmminp,self.ymminp,self.zmminp,self.l2symtry,self.l4symtry,
                self.solvergeom==w3d.RZgeom)

  def fetchfieldfrompositions(self,x,y,z,ex,ey,ez,bx,by,bz,js=0,pgroup=None):
    # --- This is called by fetchfield from fieldsolver.py
    # --- Only sets the E field from the potential
    n = len(x)
    if n == 0: return
    if top.efetch[js] == 3 and isinstance(self.fieldp,FloatType): return
    if top.efetch[js] != 3 and isinstance(self.potentialp,FloatType): return
    iselfb = top.iselfb[js]
    if (not (self.getconductorobject(top.fselfb[iselfb]).lcorrectede or
             f3d.lcorrectede)):
      sete3d(self.potentialp,self.fieldp,n,x,y,z,self.getzgridprv(),
             self.xmminp,self.ymminp,self.zmminp,
             self.dx,self.dy,self.dz,self.nxp,self.nyp,self.nzp,top.efetch[js],
             ex,ey,ez,self.l2symtry,self.l4symtry,self.solvergeom==w3d.RZgeom,
             self.nxguard,self.nyguard,self.nzguard)
    else:
      sete3dwithconductor(self.getconductorobject('p'),
             self.potentialp,self.fieldp,n,x,y,z,self.getzgridprv(),
             self.xmminp,self.ymminp,self.zmminp,
             self.dx,self.dy,self.dz,self.nxp,self.nyp,self.nzp,top.efetch[js],
             ex,ey,ez,self.l2symtry,self.l4symtry,self.solvergeom==w3d.RZgeom,
             self.nxguard,self.nyguard,self.nzguard)
    if max(abs(top.fselfb)) > 0.:
      #assert len(bx) == n,"The multigrid needs to be fixed so the B fields can be fetched with other than fetche3d"
      # --- For now, just skip the gather of the self B field if this was
      # --- called directly from fetche3dfrompositions (in which case
      # --- len(bx)==0).
      if len(bx) != n: return
      setb3d(self.fieldp[:,:,:,:,1],n,x,y,z,self.getzgridprv(),bx,by,bz,
             self.nxp,self.nyp,self.nzp,self.dx,self.dy,self.dz,
             self.xmminp,self.ymminp,self.zmminp,
             self.l2symtry,self.l4symtry,self.solvergeom==w3d.RZgeom)

  def fetchpotentialfrompositions(self,x,y,z,phi):
    n = len(x)
    if n == 0: return
    if isinstance(self.potentialp,FloatType): return
    nxp = self.nxp + 2*self.nxguard
    nyp = self.nyp + 2*self.nyguard,
    nzp = self.nzp + 2*self.nzguard
    xmminp = self.xmminp - self.dx*self.nxguard
    xmmaxp = self.xmmaxp + self.dx*self.nxguard
    ymminp = self.ymminp - self.dy*self.nyguard
    ymmaxp = self.ymmaxp + self.dy*self.nyguard
    zmminp = self.zmminp - self.dz*self.nzguard
    zmmaxp = self.zmmaxp + self.dz*self.nzguard
    getgrid3d(n,x,y,z,phi,nxp,nyp,nzp,self.potentialp,
              xmminp,xmmaxp,ymminp,ymmaxp,zmminp,zmmaxp,
              self.l2symtry,self.l4symtry)

  def setsourceforfieldsolve(self,*args):
    SubcycledPoissonSolver.setsourceforfieldsolve(self,*args)
    if self.lparallel:
      SubcycledPoissonSolver.setsourcepforparticles(self,*args)
      if isinstance(self.source,FloatType): return
      if isinstance(self.sourcep,FloatType): return
      setrhoforfieldsolve3d(self.nx,self.ny,self.nzlocal,self.source,
                            self.nxp,self.nyp,self.nzp,self.sourcep,self.nzpguard,
                            self.my_index,self.nslaves,self.izpslave,self.nzpslave,
                            self.izfsslave,self.nzfsslave)

  def getpotentialpforparticles(self,*args):
    self.setpotentialpforparticles(*args)
    if not self.lparallel:
      SubcycledPoissonSolver.getpotentialpforparticles(self,*args)
    else:
      if isinstance(self.potential,FloatType): return
      if isinstance(self.potentialp,FloatType): return
      getphipforparticles3d(1,self.nx,self.ny,self.nzlocal,self.potential,
                            self.nxp,self.nyp,self.nzp,self.potentialp,
                            self.nxguard,self.nyguard,self.nzguard,
                            self.my_index,self.nslaves,self.izpslave,self.nzpslave,
                            self.izfsslave,self.nzfsslave)

    iselfb = args[2]
    if (iselfb == 0 and
        (self.getconductorobject(top.fselfb[iselfb]).lcorrectede or
         f3d.lcorrectede)):
      # --- This only needs to be calculated once, so is only done
      # --- when iselfb == 0.
      conductorobject = self.getconductorobject('p')
      # --- This sets up the icgrid
      setupconductorfielddata(self.nxp,self.nyp,self.nzp,self.nz,
                              self.dx,self.dy,self.dz,conductorobject,
                              self.my_index,self.nslaves,self.izpslave,self.nzpslave)
      # --- This calculates the field
      getefieldatconductorsubgrid(conductorobject,
                                  self.potentialp,self.dx,self.dy,self.dz,
                                  self.nxp,self.nyp,self.nzp,
                                  self.nxguard,self.nyguard,self.nzguard,
                                  self.bounds)
    if sometrue(top.efetch == 3):
      self.setfieldpforparticles(*args)
      indts = args[1]
      # --- If this is the first group, set make sure that fieldp gets
      # --- zeroed out. Otherwise, the data in fieldp is accumulated.
      # --- This coding relies on the fact that fieldsolver does the
      # --- loops in descending order.
      tmpnsndts = getnsndtsforsubcycling()
      lzero = ((indts == tmpnsndts-1) and (iselfb == top.nsselfb-1))
      if lzero:
        tfieldp = transpose(self.fieldp)
        tfieldp[...] = 0.
      self.getselfe(recalculate=1,lzero=lzero)
      if abs(top.fselfb[iselfb]) > 0:
        # --- If the self-B correction is nonzero, then calculate and include
        # --- the approximate correction terms A and dA/dt.
        self.getselfb(self.fieldp,top.fselfb[iselfb],self.potentialp)
        self.adddadttoe(self.fieldp,top.fselfb[iselfb],self.potentialp)

      if (iselfb == 0 and
          (self.getconductorobject(top.fselfb[iselfb]).lcorrectede or
           f3d.lcorrectede)):
        # --- Now correct the E field at conductor points. This matters when
        # --- the edge of conductors are aligned with the mesh and there are no
        # --- subgrid points there.
        # --- This is done when iselfb == 0 since that will be the last
        # --- species - this routine modifies fieldp in place.
        fixefieldatconductorpoints(conductorobject,self.fieldp,
                                   self.dx,self.dy,self.dz,
                                   self.nx,self.ny,self.nz)

  def makesourceperiodic(self):
    if self.pbounds[0] == 2 or self.pbounds[1] == 2:
      self.source[0,:,:] = self.source[0,:,:] + self.source[-1,:,:]
      self.source[-1,:,:] = self.source[0,:,:]
    if self.pbounds[2] == 2 or self.pbounds[3] == 2:
      self.source[:,0,:] = self.source[:,0,:] + self.source[:,-1,:]
      self.source[:,-1,:] = self.source[:,0,:]
    if self.pbounds[0] == 1 and not self.l4symtry and self.nx > 0 and self.solvergeom != w3d.RZgeom:
      self.source[0,:,:] = 2.*self.source[0,:,:]
    if self.pbounds[1] == 1 and self.nx > 0:
      self.source[-1,:,:] = 2.*self.source[-1,:,:]
    if self.pbounds[2] == 1 and not (self.l2symtry or self.l4symtry) and self.ny > 0:
      self.source[:,0,:] = 2.*self.source[:,0,:]
    if self.pbounds[3] == 1 and self.ny > 0:
      self.source[:,-1,:] = 2.*self.source[:,-1,:]
    if self.pbounds[4] == 2 or self.pbounds[5] == 2:
      if self.lparallel:
        self.makesourceperiodic_parallel()
      else:
        self.source[:,:,0] = self.source[:,:,0] + self.source[:,:,-1]
        self.source[:,:,-1] = self.source[:,:,0]
    if self.pbounds[4] == 1 and self.izfsslave[self.my_index] == 0:
      self.source[:,:,0] = 2.*self.source[:,:,0]
    if (self.pbounds[5] == 1 and
        self.izfsslave[self.my_index]+self.nzlocal == self.nz):
      self.source[:,:,-1] = 2.*self.source[:,:,-1]

  def makesourceperiodic_parallel(self):
    tag = 70
    if self.my_index == self.nslaves-1:
      mpi.send(self.source[:,:,self.nzlocal],0,tag)
      self.source[:,:,self.nzlocal],status = mpi.recv(0,tag)
    elif self.my_index == 0:
      sourcetemp,status = mpi.recv(self.nslaves-1,tag)
      self.source[:,:,0] = self.source[:,:,0] + sourcetemp
      mpi.send(self.source[:,:,0],self.nslaves-1,tag)

  def getselfe(self,recalculate=0,lzero=true):
    # --- Make sure that fieldp is at least defined.
    try: self.fieldp
    except AttributeError: self.setfieldpforparticles(0,0,0)

    if type(self.fieldp) != ArrayType:
      # --- This should only ever be done by an external routine, such as
      # --- a plotting function.
      self.fieldp = fzeros((3,1+self.nxp,1+self.nyp,1+self.nzp,1),'d')
    if recalculate:
      if isinstance(self.potentialp,FloatType): return
      if isinstance(self.fieldp,FloatType): return
      getselfe3d(self.potentialp,self.nxp,self.nyp,self.nzp,
                 self.fieldp[:,:,:,:,0],self.nxp,self.nyp,self.nzp,
                 self.dx,self.dy,self.dz,
                 lzero,self.nxguard,self.nyguard,self.nzguard)
    return self.fieldp

  def getslicewithguard(self,i1,i2,guard):
    if i1 is not None: i1 = i1 + guard
    if i2 is not None:
      if i2 < 0: i2 = i2 - guard
      else:      i2 = i2 + guard
    if i1 is None and guard > 0: i1 = +guard
    if i2 is None and guard > 0: i2 = -guard
    return slice(i1,i2)

  def getselfb(self,fieldp,fselfb,potentialp):
    ix = self.getslicewithguard(None,None,self.nxguard)
    iy = self.getslicewithguard(None,None,self.nyguard)
    iz = self.getslicewithguard(None,None,self.nzguard)
    Az = (fselfb/clight**2)*potentialp[ix,iy,iz]
    if self.ny > 0:
      fieldp[0,:,1:-1,:,1] += (Az[:,2:,:] - Az[:,:-2,:])/(2.*self.dy)
    if self.nx > 0:
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
    ix = self.getslicewithguard(None,None,self.nxguard)
    iy = self.getslicewithguard(None,None,self.nyguard)
    # --- This assumes that nzguard is always 1
    Ez = (fselfb/clight)**2*(potentialp[ix,iy,2:]-potentialp[ix,iy,:-2])/(2.*self.dz)
    fieldp[2,:,:,:,0] += Ez

  def installconductor(self,conductor,
                            xmin=None,xmax=None,
                            ymin=None,ymax=None,
                            zmin=None,zmax=None,
                            dfill=top.largepos):
    # --- This only adds the conductor to the list. The data is only actually
    # --- installed when it is needed, during a call to getconductorobject.
    self.conductordatalist.append((conductor,xmin,xmax,ymin,ymax,zmin,zmax,dfill))

  def _installconductor(self,conductorobject,installedlist,conductordata,fselfb):
    # --- This does that actual installation of the conductor into the
    # --- conductor object

    # --- Extract the data from conductordata (the arguments to installconductor)
    conductor,xmin,xmax,ymin,ymax,zmin,zmax,dfill = conductordata

    if conductor in installedlist: return
    installedlist.append(conductor)

    if fselfb == 'p':
      zscale = 1.
      nx,ny,nzlocal,nz = self.nxp,self.nyp,self.nzp,self.nz
      xmmin,xmmax = self.xmminp,self.xmmaxp
      ymmin,ymmax = self.ymminp,self.ymmaxp
      zmmin,zmmax = self.zmminp,self.zmmaxp
      mgmaxlevels = 1
    else:
      # --- Get relativistic longitudinal scaling factor
      # --- This is quite ready yet.
      beta = fselfb/clight
      zscale = 1./sqrt((1.-beta)*(1.+beta))
      nx,ny,nzlocal,nz = self.nx,self.ny,self.nzlocal,self.nz
      xmmin,xmmax = self.xmmin,self.xmmax
      ymmin,ymmax = self.ymmin,self.ymmax
      zmmin,zmmax = self.zmmin,self.zmmax
      mgmaxlevels = None

    installconductors(conductor,xmin,xmax,ymin,ymax,zmin,zmax,dfill,
                      self.getzgrid(),
                      nx,ny,nzlocal,nz,
                      xmmin,xmmax,ymmin,ymmax,zmmin,zmmax,
                      zscale,self.l2symtry,self.l4symtry,
                      installrz=0,
                      solvergeom=self.solvergeom,conductors=conductorobject,
                      mgmaxlevels=mgmaxlevels,
                      my_index=self.my_index,nslaves=self.nslaves,
                      izfsslave=self.izfsslave,nzfsslave=self.nzfsslave)

  def hasconductors(self):
    return len(self.conductordatalist) > 0

  def clearconductors(self):
    "Clear out the conductor data"
    for fselfb in top.fselfb:
      if fselfb in self.conductorobjects:
        conductorobject = self.conductorobjects[fselfb]
        conductorobject.interior.n = 0
        conductorobject.evensubgrid.n = 0
        conductorobject.oddsubgrid.n = 0
        self.installedconductorlists[fselfb] = []

  def find_mgparam(self,lsavephi=false,resetpasses=1):
    # --- This is a temporary kludge, the same as is done in genericpf
    self.phi = self.potential
    find_mgparam(lsavephi=lsavephi,resetpasses=resetpasses,
                 solver=self,pkg3d=self)

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

    # --- Setup data for bends.
    rstar = fzeros(3+self.nzlocal,'d')
    if top.bends:

      # --- This commented out code does the same thing as the line below
      # --- setting linbend but is a bit more complicated. It is preserved
      # --- in case of some unforeseen problem with the code below.
      #ii = (top.cbendzs <= self.zmmax+zgrid and
      #                     self.zmmin+zgrid <= top.cbendze)
      #self.linbend = sometrue(ii)

      setrstar(rstar,self.nzlocal,self.dz,self.zmminlocal,self.getzgrid())
      self.linbend = min(rstar) < largepos

    if self.izfsslave is None: self.izfsslave = top.izfsslave
    if self.nzfsslave is None: self.nzfsslave = top.nzfsslave
    mgiters = zeros(1,'l')
    mgerror = zeros(1,'d')
    conductorobject = self.getconductorobject(top.pgroup.fselfb[iselfb])
    if self.electrontemperature == 0:
      multigrid3dsolve(iwhich,self.nx,self.ny,self.nzlocal,self.nz,
                       self.dx,self.dy,self.dz*zfact,self.potential,self.source,
                       rstar,self.linbend,self.bounds,
                       self.xmmin,self.ymmin,self.zmminlocal*zfact,
                       self.zmmin*zfact,self.getzgrid()*zfact,self.getzgrid()*zfact,
                       self.mgparam,self.mgform,mgiters,self.mgmaxiters,
                       self.mgmaxlevels,mgerror,self.mgtol,self.mgverbose,
                       self.downpasses,self.uppasses,
                       self.lcndbndy,self.laddconductor,self.icndbndy,
                       self.lbuildquads,self.gridmode,conductorobject,
                       self.my_index,self.nslaves,self.izfsslave,self.nzfsslave)
    else:
      iondensitygrid3d = Grid3dtype()
      setupiondensitygrid3d(self.xmmin,self.ymmin,self.zmmin,
                            self.dx,self.dy,self.dz,
                            self.nx,self.ny,self.nzlocal,
                            self.rho,iondensitygrid3d)
      self.iondensitygrid3d = iondensitygrid3d
      multigridbe3dsolve(iwhich,self.nx,self.ny,self.nzlocal,self.nz,
                         self.dx,self.dy,self.dz*zfact,self.potential,self.source,
                         rstar,self.linbend,self.bounds,
                         self.xmmin,self.ymmin,self.zmminlocal*zfact,
                         self.zmmin*zfact,self.getzgrid()*zfact,self.getzgrid()*zfact,
                         self.mgparam,mgiters,self.mgmaxiters,
                         self.mgmaxlevels,mgerror,self.mgtol,self.mgverbose,
                         self.downpasses,self.uppasses,
                         self.lcndbndy,self.laddconductor,self.icndbndy,
                         self.lbuildquads,self.gridmode,conductorobject,
                         iondensitygrid3d,
                         self.my_index,self.nslaves,self.izfsslave,
                         self.nzfsslave)
    self.mgiters = mgiters[0]
    self.mgerror = mgerror[0]

  ##########################################################################
  # Define the basic plot commands
  def genericpf(self,kw,pffunc):
    fselfb = kw.get('fselfb',top.fselfb[0])
    if 'fselfb' in kw: del kw['fselfb']
    kw['conductors'] = self.getconductorobject(fselfb)
    kw['solver'] = self
    pffunc(**kw)
  def pfxy(self,**kw): self.genericpf(kw,pfxy)
  def pfzx(self,**kw): self.genericpf(kw,pfzx)
  def pfzy(self,**kw): self.genericpf(kw,pfzy)
  def pfxyg(self,**kw): self.genericpf(kw,pfxyg)
  def pfzxg(self,**kw): self.genericpf(kw,pfzxg)
  def pfzyg(self,**kw): self.genericpf(kw,pfzyg)

  def getresidual(self):
    res = zeros(shape(self.phi),'d')
    dxsqi  = 1./self.dx**2
    dysqi  = 1./self.dy**2
    dzsqi  = 1./self.dz**2
    reps0c = self.mgparam/(eps0*2.*(dxsqi+dysqi+dzsqi))
    rho = self.rho*reps0c
    residual(self.nx,self.ny,self.nzlocal,self.nz,dxsqi,dysqi,dzsqi,
             self.phi,rho,res,0,self.bounds,self.mgparam,self.mgform,false,
             self.lcndbndy,self.icndbndy,self.conductors)
    return res


##############################################################################
##############################################################################
##############################################################################
##############################################################################
class MultiGridImplicit3D(MultiGrid):
  """
This solves the modified Poisson equation which includes the suseptibility
tensor that appears from the direct implicit scheme.
  """
  
  def __init__(self,lreducedpickle=1,**kw):
    kw['lreducedpickle'] = lreducedpickle

    SubcycledPoissonSolver.__init__(self,kwdict=kw)
    self.solvergeom = w3d.XYZgeom
    self.ncomponents = 1
    self.nxguard = 1
    self.nyguard = 1
    self.nzguard = 1

    # --- Kludge - make sure that the multigrid3df routines never sets up
    # --- any conductors. This is not really needed here.
    f3d.gridmode = 1

    # --- Save input parameters
    self.processdefaultsfrompackage(MultiGrid.__w3dinputs__,w3d,kw)
    self.processdefaultsfrompackage(MultiGrid.__f3dinputs__,f3d,kw)

    # --- If there are any remaning keyword arguments, raise an error.
    assert len(kw.keys()) == 0,"Bad keyword arguemnts %s"%kw.keys()

    # --- Create conductor objects
    self.initializeconductors()

    # --- Give these variables dummy initial values.
    self.mgiters = 0
    self.mgerror = 0.

    # --- At the start, assume that there are no bends. This is corrected
    # --- in the solve method when there are bends.
    self.linbend = false

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
    return self.source[:,:,:,0]

  def getphi(self):
    'Returns the phi array without the guard cells'
    return MultiGrid.getphi(self)[:,:,:]

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
               self.nxp,self.nyp,self.nzp,self.dx,self.dy,self.dz,
               self.xmminp,self.ymminp,self.zmminp,self.l2symtry,self.l4symtry,
               self.solvergeom==w3d.RZgeom)
    else:
      # --- Need top.pid(:,top.wpid)
      setrho3dw(sourcep,n,x,y,z,zgrid,uz,wght,q,w,top.depos,
                self.nxp,self.nyp,self.nzp,self.dx,self.dy,self.dz,
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
        setrhoforfieldsolve3d(self.nx,self.ny,self.nzlocal,
                              self.source[...,iimp],
                              self.nxp,self.nyp,self.nzp,
                              self.sourcep[...,iimp],
                              self.nzpguard,
                              self.my_index,self.nslaves,
                              self.izpslave,self.nzpslave,
                              self.izfsslave,self.nzfsslave)

  def fetchpotentialfrompositions(self,x,y,z,potential):
    n = len(x)
    if n == 0: return
    nx = self.nx + 2*self.nxguard
    ny = self.ny + 2*self.nyguard
    nzlocal = self.nzlocal + 2*self.nzguard
    xmmin = self.xmmin - self.nxguard*self.dx
    xmmax = self.xmmax + self.nxguard*self.dx
    ymmin = self.ymmin - self.nyguard*self.dy
    ymmax = self.ymmax + self.nyguard*self.dy
    zmminlocal = self.zmminlocal - self.nzguard*self.dz
    zmmaxlocal = self.zmmaxlocal + self.nzguard*self.dz
    getgrid3d(n,x,y,z,potential,nx,ny,nzlocal,self.potential,
              xmmin,xmmax,ymmin,ymmax,zmminlocal,zmmaxlocal)

  def dosolve(self,iwhich=0,*args):
    if not self.l_internal_dosolve: return
    # --- set for longitudinal relativistic contraction
    iselfb = args[2]
    beta = top.pgroup.fselfb[iselfb]/clight
    zfact = 1./sqrt((1.-beta)*(1.+beta))

    # --- This is only done for convenience.
    self.phi = self.potential
    self.rho = self.source[...,0]
    if isinstance(self.potential,FloatType): return

    # --- Setup data for bends.
    rstar = fzeros(3+self.nzlocal,'d')
    if top.bends:

      # --- This commented out code does the same thing as the line below
      # --- setting linbend but is a bit more complicated. It is preserved
      # --- in case of some unforeseen problem with the code below.
      #ii = (top.cbendzs <= self.zmmax+zgrid and
      #                     self.zmmin+zgrid <= top.cbendze)
      #self.linbend = sometrue(ii)

      setrstar(rstar,self.nzlocal,self.dz,self.zmminlocal,self.getzgrid())
      self.linbend = min(rstar) < largepos

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

    mgsolveimplicites3d(iwhich,self.nx,self.ny,self.nzlocal,self.nz,
                        self.dx,self.dy,self.dz*zfact,
                        self.potential,self.rho,
                        top.nsimplicit,qomdt,self.chi0,
                        rstar,self.linbend,
                        self.bounds,self.xmmin,self.ymmin,
                        self.zmminlocal*zfact,self.zmmin*zfact,
                        self.getzgrid()*zfact,self.getzgrid()*zfact,
                        self.mgparam,mgiters,self.mgmaxiters,
                        self.mgmaxlevels,mgerror,self.mgtol,
                        self.mgverbose,
                        self.downpasses,self.uppasses,
                        self.lcndbndy,self.laddconductor,self.icndbndy,
                        self.lbuildquads,self.gridmode,conductorobject,
                        self.my_index,self.nslaves,
                        self.izfsslave,self.nzfsslave)

    self.mgiters = mgiters[0]
    self.mgerror = mgerror[0]

# --- This can only be done after MultiGrid is defined.
try:
  psyco.bind(MultiGrid)
  psyco.bind(MultiGridImplicit3D)
except NameError:
  pass


