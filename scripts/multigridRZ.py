"""Class for doing complete multigrid field solve"""
# ToDo:
#  - modify setrho to check if particles are within grid
#  - incorporate instances into the particle mover, so charge is deposited and
#    the E fields gather appropriately.
from warp import *
from generateconductors import installconductors
from find_mgparam import find_mgparam
import MA

try:
  import psyco
except ImportError:
  pass

##############################################################################
class MultiGridRZ(object):
  
  __w3dinputs__ = ['nx','nz','nzfull','nzpguard',
                   'xmmin','xmmax','zmmin','zmmax',
                   'zmminglobal','zmmaxglobal',
                   'bound0','boundnz','boundxy','l2symtry','l4symtry',
                   'solvergeom',
                   'iondensity','electrontemperature','plasmapotential',
                   'electrondensitymaxscale']
  __GRIDtypeinputs__ = ['nguardx','nguardz']
  __f3dinputs__ = ['gridmode','mgparam','downpasses','uppasses',
                   'mgmaxiters','mgtol','mgmaxlevels','mgform',
                   'lcndbndy','icndbndy','laddconductor'] 
  __frzinputs__ = ['mgridrz_ncycles']
  __topinputs__ = ['pbound0','pboundnz','pboundxy','efetch',
                   'my_index','nslaves','izfsslave','nzfsslave']
  __flaginputs__ = {'forcesymmetries':1,'lzerorhointerior':0,
                    'lreducedpickle':0}

  def __init__(self,**kw):
    self.solvergeom = w3d.RZgeom

    # --- Kludge - make sure that the multigrid3df routines never sets up
    # --- any conductors.
    f3d.gridmode = 1

    # --- Save input parameters
    for name in MultiGridRZ.__w3dinputs__:
      if name not in self.__dict__:
        #self.__dict__[name] = kw.pop(name,getattr(w3d,name)) # Python2.3
        self.__dict__[name] = kw.get(name,getattr(w3d,name))
      if kw.has_key(name): del kw[name]
    for name in MultiGridRZ.__GRIDtypeinputs__:
      if name not in self.__dict__:
        #self.__dict__[name] = kw.pop(name,getattr(w3d,name)) # Python2.3
        self.__dict__[name] = kw.get(name,getattr(w3d,name))
      if kw.has_key(name): del kw[name]
    for name in MultiGridRZ.__f3dinputs__:
      if name not in self.__dict__:
        #self.__dict__[name] = kw.pop(name,getattr(f3d,name)) # Python2.3
        self.__dict__[name] = kw.get(name,getattr(f3d,name))
      if kw.has_key(name): del kw[name]
    for name in MultiGridRZ.__topinputs__:
      if name not in self.__dict__:
        #self.__dict__[name] = kw.pop(name,getattr(top,name)) # Python2.3
        self.__dict__[name] = kw.get(name,getattr(top,name))
      if kw.has_key(name): del kw[name]
    for name,defvalue in MultiGridRZ.__flaginputs__.iteritems():
      if name not in self.__dict__:
        #self.__dict__[name] = kw.pop(name,getattr(top,name)) # Python2.3
        self.__dict__[name] = kw.get(name,defvalue)
      if kw.has_key(name): del kw[name]

    # --- bounds is special since it will sometimes be set from the
    # --- variables bound0, boundnz, boundxy, l2symtry, and l4symtry
    if 'bounds' not in self.__dict__:
      if 'bounds' in kw:
        self.bounds = kw['bounds']
      else:
        self.bounds = zeros(6)
        self.bounds[0] = self.boundxy
        self.bounds[1] = self.boundxy
        self.bounds[2] = neumann
        self.bounds[3] = neumann
        self.bounds[4] = self.bound0
        self.bounds[5] = self.boundnz
        if self.solvergeom != w3d.RZgeom:
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
      else:
        self.pbounds = zeros(6)
        self.pbounds[0] = self.pboundxy
        self.pbounds[1] = self.pboundxy
        self.pbounds[2] = neumann
        self.pbounds[3] = neumann
        self.pbounds[4] = self.pbound0
        self.pbounds[5] = self.pboundnz
        if self.solvergeom != w3d.RZgeom:
          if self.l2symtry:
            self.pbounds[2] = reflect
            if self.pboundxy == periodic: self.pbounds[3] = reflect
          elif self.l4symtry:
            self.pbounds[0] = reflect
            self.pbounds[2] = reflect
            if self.pboundxy == periodic: self.pbounds[1] = reflect
            if self.pboundxy == periodic: self.pbounds[3] = reflect

    # --- If there are any remaning keyword arguments, raise an error.
    assert len(kw.keys()) == 0,"Bad keyword arguemnts %s"%kw.keys()

    # --- Set set parallel related paramters and calculate mesh sizes
    if self.nslaves <= 1:
      self.nzfull = self.nz
      self.zmminglobal = self.zmmin
      self.zmmaxglobal = self.zmmax
      self.izfsslave = zeros(1)
      self.nzfsslave = zeros(1) + self.nz
    self.dx = (self.xmmax - self.xmmin)/self.nx
    self.dy = self.dx
    self.dz = (self.zmmaxglobal - self.zmminglobal)/self.nzfull
    self.ny = 0

    self.xmesh = self.xmmin + arange(0,self.nx+1)*self.dx
    self.ymesh = self.ymmin + arange(0,self.ny+1)*self.dy
    self.zmesh = self.zmminglobal + arange(0,self.nzfull+1)*self.dz
    self.zmeshlocal = self.zmmin + arange(0,self.nz+1)*self.dz

    # --- Create extra variables which are used in various places
    self.nxp = self.nx
    self.nyp = self.ny
    self.nzp = self.nz

    # --- Create phi and rho arrays and other arrays.
    self.allocatefieldarrays()

    # --- Create a conductor object, which by default is empty.
    self.conductorlist = []

    # --- These must be arrays since they are modified in the call to the
    # --- MG solver.
    self.mgiters = zeros(1)
    self.mgerror = zeros(1,'d')

    # --- At the start, assume that there are no bends. This is corrected
    # --- in the solve method when there are bends.
    self.linbend = false

  def __getstate__(self):
    dict = self.__dict__.copy()
    if self.lreducedpickle:
      del dict['rho']
      del dict['rhop']
      del dict['phi']
      del dict['phip']
    return dict

  def __setstate__(self,dict):
    self.__dict__.update(dict)
    if self.lreducedpickle:
      self.allocatefieldarrays()
      # --- Regenerate the conductor data
      conductorlist = self.conductorlist
      self.conductorlist = []
      for conductor in conductorlist:
        self.installconductor(conductor)

  def allocatefieldarrays(self):
    self.grid = GRIDtype()
    bg = self.grid

    inveps0 = 1./eps0

    if w3d.solvergeom==w3d.Zgeom:
      self.nguardx = 0
    END if
    if w3d.solvergeom==w3d.Rgeom:
      self.nguardz = 0
    END if

    bg.nzp = self.nz
    bg.nrpar = 0
    bg.nzpar = 0
    bg.nr = self.nx
    bg.dr = self.dx
    bg.rmin = self.xmmin
    bg.rmax = self.xmmin+bg.nr*bg.dr
    bg.xmin = self.xmmin
    bg.xmax = self.xmmin+bg.nr*bg.dr
    bg.nz = self.nz
    bg.zminp = self.zmmin
    bg.dz = self.dz
    bg.zmin = self.zmmin
    bg.zmax = self.zmmin+bg.nz*bg.dz
    bg.jmin = 1
    bg.jmax = bg.nr+1
    bg.lmin = 1
    bg.lmax = bg.nr+1
    bg.nguardx = self.nguardx
    bg.nguardz = self.nguardz
    gchange(bg)
    bg.gid = 1
    bg.loc_part = 1
    bg.loc_part_fd = 1
    bg.mgparam = self.mgparam
    bg.npre = self.downpasses
    bg.npost = self.uppasses
    bg.ncycles = self.mgridrz_ncycles
    bg.ncmax = self.mgmaxiters
    bg.npmin = self.mgmaxlevels
    bg.phi = 0.
    bg.rho = 0.
    bg.transit_min_r = 0
    bg.transit_max_r = 0
    bg.transit_min_z = 0
    bg.transit_max_z = 0
    frz.ngrids = 1
    frz.level_del_grid = 0
    frz.n_avail_ids = 0
    frz.avail_ids = -1
    bg.invdr = 1./bg.dr
    bg.invdz = 1./bg.dz
    if solvergeom==RZgeom or solvergeom==Rgeom:
#     --- computes divider by cell volumes to get density
      if bg.rmin==0.:
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
      if solvergeom==Rgeom: bg.invvol = bg.invvol*bg.dz
    else:
#     --- solvergeom==XZgeom or solvergeom==XYgeom
      bg.invvol = 1./(bg.dr*bg.dz)

    if solvergeom==RZgeom or solvergeom==XZgeom:
      ixlbndi = self.boundxy
      ixrbndi = self.boundxy
      izlbndi = self.bound0
      izrbndi = self.boundnz
      ixlbnd = self.boundxy
      ixrbnd = self.boundxy
      izlbnd = self.bound0
      izrbnd = self.boundnz
      if bg.rmin==0.:
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
      if (w3d.l2symtry or w3d.l4symtry) and solvergeom==w3d.XYgeom:
        izlbndi = neumann
        izlbnd  = neumann

    bg.ixlbnd = ixlbndi
    bg.ixrbnd = ixrbndi
    bg.izlbnd = izlbndi
    bg.izrbnd = izrbndi

    if not (self.solvergeom==w3d.Zgeom or solvergeom==w3d.Rgeom):
      init_bnd(bg,bg.nr,bg.nz,bg.dr,bg.dz,bg.zmin,bg.zmax)
    bg.nlevels = self.mgmaxlevels

    for i in range(bg.nlevels):
      if i==0:
        b = bg.bndfirst
      else:
        b = b.next
      if b.izlbnd==dirichlet:  b.v[:,0]    = 3 # v_dirichlet
      if b.izrbnd==dirichlet:  b.v[:,b.nz] = 3 # v_dirichlet
      if bg.ixlbnd==dirichlet: b.v[0,:]    = 3 # v_dirichlet
      if bg.ixrbnd==dirichlet: b.v[b.nr,:] = 3 # v_dirichlet
    setmglevels_rz(bg)
  
  def setrho(self,x,y,z,uz,q,w):
    n = len(x)
    if n == 0: return
    if top.wpid == 0:
      rhoweightrzgrid(self.grid,x,y,z,n,q*w,self.nx,self.nz,self.dx,self.dz,
                      self.xmmin,top.zgrid)
    else:
      # --- Need top.pid(:,top.wpid)
      rhoweightrzgrid_weights(self.grid,x,y,z,n,q*w,self.nx,self.nz,self.dx,self.dz,
                      self.xmmin,top.zgrid)
      

  def setrhoselect(self,x,y,z,uz,q,w):
    n = len(x)
    if n == 0: return
    setrho3dselect(self.rho,self.rho,n,x,y,z,top.zgrid,uz,q,w,top.depos,
             self.nx,self.ny,self.nz,self.dx,self.dy,self.dz,
             self.xmmin,self.ymmin,self.zmmin,self.l2symtry,self.l4symtry)

  def fetchefrompositions(self,x,y,z,ex,ey,ez):
    n = len(x)
    if n == 0: return
    sete3d(self.phi,self.selfe,n,x,y,z,top.zgridprv,
           self.xmmin,self.ymmin,self.zmmin,
           self.dx,self.dy,self.dz,self.nx,self.ny,self.nz,self.efetch,
           ex,ey,ez,self.l2symtry,self.l4symtry)

  def fetchphifrompositions(self,x,y,z,phi):
    n = len(x)
    getgrid3d(n,x,y,z,phi,self.nx,self.ny,self.nz,self.phi[:,:,1:-1],
              self.xmmin,self.xmmax,self.ymmin,self.ymmax,self.zmmin,self.zmmax,
              self.l2symtry,self.l4symtry)

  def loadrho(self,ins_i=-1,nps_i=-1,is_i=-1,lzero=true):
    if lzero: self.rho[...] = 0.
    for i,n,q,w in zip(top.pgroup.ins-1,top.pgroup.nps,
                       top.pgroup.sq,top.pgroup.sw):
      self.setrho(top.pgroup.xp[i:i+n],top.pgroup.yp[i:i+n],
                  top.pgroup.zp[i:i+n],top.pgroup.uzp[i:i+n],q,w)
    self.makerhoperiodic()
    self.getrhoforfieldsolve()
    if self.lzerorhointerior:
      cond_zerorhointerior(self.conductors.interior,self.nx,self.ny,self.nz,self.rho)

  def makerhoperiodic(self):
    if self.pbounds[0] == 2 or self.pbounds[1] == 2:
      self.rho[0,:,:] = self.rho[0,:,:] + self.rho[-1,:,:]
      self.rho[-1,:,:] = self.rho[0,:,:]
    if self.pbounds[2] == 2 or self.pbounds[3] == 2:
      self.rho[:,0,:] = self.rho[:,0,:] + self.rho[:,-1,:]
      self.rho[:,-1,:] = self.rho[:,0,:]
    if self.pbounds[0] == 1 and not self.l4symtry:
       self.rho[0,:,:] = 2.*self.rho[0,:,:]
    if self.pbounds[1] == 1: self.rho[-1,:,:] = 2.*self.rho[-1,:,:]
    if self.pbounds[2] == 1 and not (self.l2symtry or self.l4symtry):
       self.rho[:,0,:] = 2.*self.rho[:,0,:]
    if self.pbounds[3] == 1: self.rho[:,-1,:] = 2.*self.rho[:,-1,:]
    if self.pbounds[4] == 2 or self.pbounds[5] == 2:
      if self.nslaves > 1:
        self.makerhoperiodic_parallel()
      else:
        self.rho[:,:,0] = self.rho[:,:,0] + self.rho[:,:,-1]
        self.rho[:,:,-1] = self.rho[:,:,0]
    if self.pbounds[4] == 1: self.rho[:,:,0] = 2.*self.rho[:,:,0]
    if self.pbounds[5] == 1: self.rho[:,:,-1] = 2.*self.rho[:,:,-1]

  def getrhoforfieldsolve(self):
    if self.nslaves > 1:
      getrhoforfieldsolve3d(self.nx,self.ny,self.nz,self.rho,
                            self.nx,self.ny,self.nz,self.rho,self.nzpguard)

  def makerhoperiodic_parallel(self):
    tag = 70
    if me == self.nslaves-1:
      request = mpi.isend(self.rho[:,:,self.nz],0,tag)
      self.rho[:,:,self.nz],status = mpi.recv(0,tag)
    elif me == 0:
      rhotemp,status = mpi.recv(self.nslaves-1,tag)
      self.rho[:,:,0] = self.rho[:,:,0] + rhotemp
      request = mpi.isend(self.rho[:,:,0],self.nslaves-1,tag)
    if me == 0 or me == self.nslaves-1:
      status = request.wait()

  def fetche(self):
    if w3d.api_xlf2:
      w3d.xfsapi=top.pgroup.xp[w3d.ipminfsapi-1:w3d.ipminfsapi-1+w3d.ipfsapi]
      w3d.yfsapi=top.pgroup.yp[w3d.ipminfsapi-1:w3d.ipminfsapi-1+w3d.ipfsapi]
      w3d.zfsapi=top.pgroup.zp[w3d.ipminfsapi-1:w3d.ipminfsapi-1+w3d.ipfsapi]
    self.fetchefrompositions(w3d.xfsapi,w3d.yfsapi,w3d.zfsapi,
                             w3d.exfsapi,w3d.eyfsapi,w3d.ezfsapi)

  def fetchphi(self):
    self.fetchphifrompositions(w3d.xfsapi,w3d.yfsapi,w3d.zfsapi,w3d.phifsapi)

  def getselfe(self,recalculate=0):
    if type(self.selfe) != ArrayType:
      self.selfe = fzeros((3,1+self.nx,1+self.ny,1+self.nz),'d')
    if recalculate:
      getselfe3d(self.phi,self.nx,self.ny,self.nz,
                 self.selfe,self.nx,self.ny,self.nz,self.dx,self.dy,self.dz,
                 self.bounds[0],self.bounds[1],self.bounds[2],self.bounds[3])
    return self.selfe

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

  def clearconductors(self):
    self.conductors.interior.n = 0
    self.conductors.evensubgrid.n = 0
    self.conductors.oddsubgrid.n = 0

  def optimizeconvergence(self,resetpasses=1):
    find_mgparam(resetpasses=resetpasses,solver=self)

  def solve(self,iwhich=0):
    # --- Setup data for bends.
    if top.bends:

      # --- This commented out code does the same thing as the line below
      # --- setting linbend but is a bit more complicated. It is preserved
      # --- in case of some unforeseen problem with the code below.
      #ii = (top.cbendzs <= self.zmmax+top.zgrid and
      #                     self.zmmin+top.zgrid <= top.cbendze)
      #self.linbend = sometrue(ii)

      setrstar(self.rstar,self.nz,self.dz,self.zmmin,top.zgrid)
      self.linbend = min(self.rstar) < largepos

    if self.electrontemperature == 0:
      multigrid3dsolve(iwhich,self.nx,self.ny,self.nz,self.nzfull,
                       self.dx,self.dy,self.dz,self.phi,self.rho,
                       self.rstar,self.linbend,self.bounds,
                       self.xmmin,self.ymmin,self.zmmin,top.zbeam,top.zgrid,
                       self.mgparam,self.mgform,self.mgiters,self.mgmaxiters,
                       self.mgmaxlevels,self.mgerror,self.mgtol,
                       self.downpasses,self.uppasses,
                       self.lcndbndy,self.laddconductor,self.icndbndy,
                       self.lbuildquads,
                       self.gridmode,
                       self.conductors,
                       self.my_index,self.nslaves,self.izfsslave,self.nzfsslave)
    else:
      multigridbe3dsolve(iwhich,self.nx,self.ny,self.nz,self.nzfull,
                         self.dx,self.dy,self.dz,self.phi,self.rho,
                         star,self.linbend,self.bounds,
                         self.xmmin,self.ymmin,self.zmmin,top.zbeam,top.zgrid,
                         self.mgparam,self.mgiters,self.mgmaxiters,
                         self.mgmaxlevels,self.mgerror,self.mgtol,
                         self.downpasses,self.uppasses,
                         self.lcndbndy,self.laddconductor,self.icndbndy,
                         self.lbuildquads,self.gridmode,self.conductors,
                         self.my_index,self.nslaves,self.izfsslave,
                         self.nzfsslave,
                         self.iondensity,self.electrontemperature,
                         self.plasmapotential,self.electrondensitymaxscale)

    if self.efetch == 3:
      MultiGridRZ.getselfe(self,recalculate=1)

  ##########################################################################
  # Define the basic plot commands
  def genericpf(self,kw,pffunc):
    kw['conductors'] = self.conductors
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

    checkconductors(self.nx,self.ny,self.nz,self.nzfull,
                    self.dx,self.dy,self.dz,self.conductors,
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
    self.mgiters[0] = 0
    self.mgerror[0] = 2.*self.mgtol + 1.
    while (self.mgerror[0] > self.mgtol and self.mgiters[0] < self.mgmaxiters):
      self.mgiters[0] = self.mgiters[0] + 1

      # --- Save current value of phi
      phisave[:,:,:] = self.phi + 0.

      # --- If using residual correction form, calculate the residual and
      # --- copy it into rhosave, zero phisave (the initial error).
      # --- In the calls to cond_potmg and residual, the last argument
      # --- is true, telling the routines to use the actual value of
      # --- voltages rather than zero as is done otherwise for residual
      # --- correction form since it is operating on the error.
      if self.mgform == 2:
        cond_potmg(self.conductors.interior,
                   self.nx,self.ny,self.nz,phisave,0,self.mgform,true)
        residual(self.nx,self.ny,self.nz,self.nzfull,dxsqi,dysqi,dzsqi,
                 phisave,rhosave,res,0,localbounds,
                 self.mgparam,self.mgform,true,
                 self.lcndbndy,self.icndbndy,self.conductors)
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
                  self.icndbndy,self.conductors)

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
      self.mgerror[0] = MA.maximum(phisave)

    #ifdef MPIPARALLEL
    #     --- calculate global sorerror
    #     call parallelmaxrealarray(self.mgerror[0],1)

    # --- For Dirichlet boundary conditions, copy data into guard planes
    # --- For other boundary conditions, the guard planes are used during
    # --- the solve are so are already set.
    if (self.bounds[4] == 0): self.phi[:,:,0] = self.phi[:,:,1]
    if (self.bounds[5] == 0): self.phi[:,:,-1] = self.phi[:,:,-2]

    # --- Make a print out.
    if (self.mgerror[0] > self.mgtol):
      print "Multigrid: Maximum number of iterations reached"
    print ("Multigrid: Error converged to %11.3e in %4d v-cycles"%
           (self.mgerror[0],self.mgiters[0]))

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

# --- This can only be done after MultiGridRZ is defined.
try:
  psyco.bind(MultiGridRZ)
except NameError:
  pass

