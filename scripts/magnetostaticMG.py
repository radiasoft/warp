"""Class for doing complete magnetostatic multigrid field solve"""
# ToDo:
#  - modify setj to check if particles are within grid
from warp import *
from generateconductors import installconductors
from find_mgparam import find_mgparam
import MA

try:
  import psyco
except ImportError:
  pass

##############################################################################
class MagnetostaticMG(object):
  
  __w3dinputs__ = ['nx','ny','nz','nzfull','nzpguard',
                   'xmmin','xmmax','ymmin','ymmax','zmmin','zmmax',
                   'zmminglobal','zmmaxglobal',
                   'bound0','boundnz','boundxy','l2symtry','l4symtry',
                   'solvergeom']
  __bfieldinputs__ = ['mgparam','downpasses','uppasses',
                      'mgmaxiters','mgtol','mgmaxlevels','mgform',
                      'lcndbndy','icndbndy','laddconductor',
                      'lcylindrical','lanalyticbtheta'] 
  __f3dinputs__ = ['gridmode','mgparam','downpasses','uppasses',
                   'mgmaxiters','mgtol','mgmaxlevels','mgform',
                   'lcndbndy','icndbndy','laddconductor'] 
  __topinputs__ = ['pbound0','pboundnz','pboundxy','efetch',
                   'my_index','nslaves','izfsslave','nzfsslave']
  __flaginputs__ = {'forcesymmetries':1,'lzerorhointerior':0,
                    'lreducedpickle':0}

  def __init__(self,**kw):
    #self.solvergeom = w3d.XYZgeom

    # --- Kludge - make sure that the multigrid3df routines never sets up
    # --- any conductors.
    f3d.gridmode = 1

    # --- Save input parameters
    for name in MagnetostaticMG.__w3dinputs__:
      if name not in self.__dict__:
        #self.__dict__[name] = kw.pop(name,getattr(w3d,name)) # Python2.3
        self.__dict__[name] = kw.get(name,getattr(w3d,name))
      if kw.has_key(name): del kw[name]
    for name in MagnetostaticMG.__bfieldinputs__:
      if name not in self.__dict__:
        #self.__dict__[name] = kw.pop(name,getattr(f3d.bfield,name))#Python2.3
        self.__dict__[name] = kw.get(name,getattr(f3d.bfield,name))
      if kw.has_key(name): del kw[name]
    for name in MagnetostaticMG.__f3dinputs__:
      if name not in self.__dict__:
        #self.__dict__[name] = kw.pop(name,getattr(f3d,name)) # Python2.3
        self.__dict__[name] = kw.get(name,getattr(f3d,name))
      if kw.has_key(name): del kw[name]
    for name in MagnetostaticMG.__topinputs__:
      if name not in self.__dict__:
        #self.__dict__[name] = kw.pop(name,getattr(top,name)) # Python2.3
        self.__dict__[name] = kw.get(name,getattr(top,name))
      if kw.has_key(name): del kw[name]
    for name,defvalue in MagnetostaticMG.__flaginputs__.iteritems():
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

    # --- pbounds is special since it will sometimes be set from the
    # --- variables pbound0, pboundnz, pboundxy, l2symtry, and l4symtry
    if 'pbounds' not in self.__dict__:
      if 'pbounds' in kw:
        self.pbounds = kw['pbounds']
      else:
        self.pbounds = zeros(6)
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

    # --- If there are any remaning keyword arguments, raise an error.
    assert len(kw.keys()) == 0,"Bad keyword arguemnts %s"%kw.keys()

    # --- Check for cylindrical geometry
    self.lcylindrical = (self.solvergeom == w3d.RZgeom) or self.lcylindrical

    # --- Fix the values of xmmin and ymmin is necessary
    if self.lcylindrical:
      self.ny = 0
      if self.xmmin < 0.: self.xmmin = 0.
      if self.ymmin < 0.: self.ymmin = 0.
      self.l2symtry = false
      self.l4symtry = false

    # --- Set set parallel related paramters and calculate mesh sizes
    if self.nslaves <= 1:
      self.nzfull = self.nz
      self.zmminglobal = self.zmmin
      self.zmmaxglobal = self.zmmax
      self.izfsslave = zeros(1)
      self.nzfsslave = zeros(1) + self.nz
    if self.nx > 0: self.dx = (self.xmmax - self.xmmin)/self.nx
    if self.ny > 0: self.dy = (self.ymmax - self.ymmin)/self.ny
    if self.nzfull > 0: self.dz = (self.zmmaxglobal - self.zmminglobal)/self.nzfull
    self.xsymmetryplane = 0.
    self.ysymmetryplane = 0.

    if self.lcylindrical:
      self.dy = self.dx
      if w3d.l4symtry or (self.lcylindrical and self.xmmin==0.):
        self.bounds[0] = neumann
      self.bounds[2] = dirichlet
      self.bounds[3] = dirichlet

    self.xmesh = self.xmmin + arange(0,self.nx+1)*self.dx
    self.ymesh = self.ymmin + arange(0,self.ny+1)*self.dy
    self.zmesh = self.zmmin + arange(0,self.nz+1)*self.dz

    # --- Create a conductor object, which by default is empty.
    self.conductorlist = []

    # --- Create phi and rho arrays and other arrays.
    self.allocatefieldarrays()

    # --- At the start, assume that there are no bends. This is corrected
    # --- in the solve method when there are bends.
    self.linbend = false

    # --- Turn of build quads option
    self.lbuildquads = false

  def __getstate__(self):
    dict = self.__dict__.copy()
    if self.lreducedpickle:
      del dict['j']
      del dict['b']
      del dict['a']
      del dict['bfield']
      del dict['conductors']
    return dict

  def __setstate__(self,dict):
    self.__dict__.update(dict)
    if self.lreducedpickle:
      self.allocatefieldarrays()
      # --- Regenerate the conductor data
      self.conductors = ConductorType()
      conductorlist = self.conductorlist
      self.conductorlist = []
      for conductor in conductorlist:
        self.installconductor(conductor)

  def allocatefieldarrays(self):

    self.bfield = BFieldGridType()
    self.bfield.nx = self.nx
    self.bfield.ny = self.ny
    self.bfield.nz = self.nz
    self.bfield.nzfull = self.nzfull
    self.bfield.dx = self.dx
    self.bfield.dy = self.dy
    self.bfield.dz = self.dz
    self.bfield.xmmin = self.xmmin
    self.bfield.xmmax = self.xmmax
    self.bfield.ymmin = self.ymmin
    self.bfield.ymmax = self.ymmax
    self.bfield.zmmin = self.zmmin
    self.bfield.zmmax = self.zmmax
    self.bfield.zmminglobal = self.zmminglobal
    self.bfield.zmmaxglobal = self.zmmaxglobal
    self.bfield.lcndbndy = self.lcndbndy
    self.bfield.icndbndy = self.icndbndy
    self.bfield.laddconductor = self.laddconductor
    self.bfield.mgparam = self.mgparam
    self.bfield.mgmaxiters = self.mgmaxiters
    self.bfield.mgmaxlevels = self.mgmaxlevels
    self.bfield.mgiters = zeros(3)
    self.bfield.mgtol = self.mgtol
    self.bfield.mgerror = zeros(3,'d')
    self.bfield.mgform = self.mgform
    self.bfield.downpasses = self.downpasses
    self.bfield.uppasses = self.uppasses
    self.bfield.bounds = self.bounds
    self.bfield.lcylindrical = self.lcylindrical
    self.bfield.lusevectorpotential = true
    self.bfield.lanalyticbtheta = self.lanalyticbtheta
    self.bfield.j = fzeros((3,1+self.nx,1+self.ny,1+self.nz),'d')
    self.bfield.b = fzeros((3,1+self.nx,1+self.ny,1+self.nz),'d')
    self.bfield.a = fzeros((3,3+self.nx,3+self.ny,3+self.nz),'d')
    self.bfield.rstar = fzeros(3+self.nz,'d')
    self.bfield.conductors = ConductorType()

    # --- Delete the names from the dictionary so there is no confusion
#   del self.nx
#   del self.ny
#   del self.nz
#   del self.nzfull
#   del self.dx
#   del self.dy
#   del self.dz
#   del self.xmmin
#   del self.xmmax
#   del self.ymmin
#   del self.ymmax
#   del self.zmmin
#   del self.zmmax
#   del self.zmminglobal
#   del self.zmmaxglobal
#   del self.lcndbndy
#   del self.icndbndy
#   del self.laddconductor
#   del self.mgparam
#   del self.mgmaxiters
#   del self.mgmaxlevels
#   del self.mgtol
#   del self.mgform
#   del self.downpasses
#   del self.uppasses
#   del self.bounds
#   del self.lcylindrical

  def setj(self,x,y,z,ux,uy,uz,gaminv,wght,q,w):
    n = len(x)
    if n == 0: return
    if wght is not None:
      nw = len(wght)
    else:
      nw = 0.
      wght = 0.
    setj3d(self.bfield,self.bfield.j,n,x,y,z,top.zgrid,ux,uy,uz,gaminv,
           q*w,nw,wght,top.depos,self.l2symtry,self.l4symtry)

  def fetchbfrompositions(self,x,y,z,bx,by,bz):
    n = len(x)
    if n == 0: return
    setb3d(self.bfield,n,x,y,z,top.zgridprv,
           bx,by,bz,self.l2symtry,self.l4symtry)

  def fetchafrompositions(self,x,y,z,a):
    n = len(x)
    if n == 0: return
    fetchafrompositions3d(n,x,y,z,a,top.zgrid,self.bfield,
                          self.l2symtry,self.l4symtry)

  def loadj(self,ins_i=-1,nps_i=-1,is_i=-1,lzero=true):
    if lzero: self.bfield.j[...] = 0.
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

  def makejperiodic(self):
    j = self.bfield.j
    if self.pbounds[0] == 2 or self.pbounds[1] == 2:
      j[:,0,:,:] = j[:,0,:,:] + j[:,-1,:,:]
      j[:,-1,:,:] = j[:,0,:,:]
    if self.pbounds[2] == 2 or self.pbounds[3] == 2:
      j[:,:,0,:] = j[:,:,0,:] + j[:,:,-1,:]
      j[:,:,-1,:] = j[:,:,0,:]
    if self.pbounds[0] == 1 and not self.l4symtry:
       j[:,0,:,:] = 2.*j[:,0,:,:]
    if self.pbounds[1] == 1: j[:,-1,:,:] = 2.*j[:,-1,:,:]
    if self.pbounds[2] == 1 and not (self.l2symtry or self.l4symtry):
       j[:,:,0,:] = 2.*j[:,:,0,:]
    if self.pbounds[3] == 1: j[:,:,-1,:] = 2.*j[:,:,-1,:]
    if self.pbounds[4] == 2 or self.pbounds[5] == 2:
      if self.nslaves > 1:
        self.makejperiodic_parallel()
      else:
        j[:,:,:,0] = j[:,:,:,0] + j[:,:,:,-1]
        j[:,:,:,-1] = j[:,:,:,0]
    if self.pbounds[4] == 1: j[:,:,:,0] = 2.*j[:,:,:,0]
    if self.pbounds[5] == 1: j[:,:,:,-1] = 2.*j[:,:,:,-1]

  def getjforfieldsolve(self):
    if self.nslaves > 1:
      getjforfieldsolve3d(self.bfield,self.bfield)

  def makejperiodic_parallel(self):
    tag = 70
    j = self.bfield.j
    if me == self.nslaves-1:
      request = mpi.isend(j[:,:,:,self.nz],0,tag)
      j[:,:,:,self.nz],status = mpi.recv(0,tag)
    elif me == 0:
      jtemp,status = mpi.recv(self.nslaves-1,tag)
      j[:,:,:,0] = j[:,:,:,0] + jtemp
      request = mpi.isend(j[:,:,:,0],self.nslaves-1,tag)
    if me == 0 or me == self.nslaves-1:
      status = request.wait()

  def fetchb(self):
    if w3d.api_xlf2:
      w3d.xfsapi=top.pgroup.xp[w3d.ipminapi-1:w3d.ipminapi-1+w3d.ipapi]
      w3d.yfsapi=top.pgroup.yp[w3d.ipminapi-1:w3d.ipminapi-1+w3d.ipapi]
      w3d.zfsapi=top.pgroup.zp[w3d.ipminapi-1:w3d.ipminapi-1+w3d.ipapi]
    self.fetchbfrompositions(w3d.xfsapi,w3d.yfsapi,w3d.zfsapi,
                             w3d.bxfsapi,w3d.byfsapi,w3d.bzfsapi)

  def fetcha(self):
    self.fetchafrompositions(w3d.xfsapi,w3d.yfsapi,w3d.zfsapi,w3d.afsapi)

  def installconductor(self,conductor,
                            xmin=None,xmax=None,
                            ymin=None,ymax=None,
                            zmin=None,zmax=None,
                            dfill=top.largepos):
    if conductor in self.conductorlist: return
    self.conductorlist.append(conductor)
    bfield = self.bfield
    installconductors(conductor,xmin,xmax,ymin,ymax,zmin,zmax,dfill,
                      top.zbeam,
                      bfield.nx,bfield.ny,bfield.nz,bfield.nzfull,
                      bfield.xmmin,bfield.xmmax,bfield.ymmin,bfield.ymmax,
                      bfield.zmmin,bfield.zmmax,self.l2symtry,self.l4symtry,
                      solvergeom=self.solvergeom,
                      conductors=self.conductors)

  def clearconductors(self):
    self.conductors.interior.n = 0
    self.conductors.evensubgrid.n = 0
    self.conductors.oddsubgrid.n = 0

  def optimizeconvergence(self,resetpasses=1):
    pass
    #find_mgparam(resetpasses=resetpasses,solver=self)

  def solve(self,iwhich=0):
    bfield = self.bfield
    # --- Setup data for bends.
    if top.bends:
      setrstar(self.rstar,bfield.nz,bfield.dz,bfield.zmmin,
               top.zgrid)
      self.linbend = min(self.rstar) < largepos

    bfield.j[...] = bfield.j*mu0*eps0

    if self.lcylindrical:
      init_bworkgrid(bfield.nx,bfield.nz,bfield.dx,bfield.dz,
                     bfield.xmmin,bfield.zmmin,bfield.bounds,
                     bfield.a[0,:,0,:],bfield.j[0,:,0,:])

    # --- Note that the arrays being passed in are not contiguous, which means
    # --- that copies are being done.
    # --- If only initialization is being done (iwhich==1) then the bvp3d_work
    # --- routine only needs to be called once. Proper arrays are still passed
    # --- though they should never be needed during initialization.
    idmax = 2
    if iwhich == 1: idmax = 0
    for id in range(idmax+1):
      if (bfield.lanalyticbtheta and
          ((bfield.lusevectorpotential and (id == 0 or id == 2)) or
          (not bfield.lusevectorpotential and id == 1))): continue

      if self.lcylindrical:
        multigridrzb(iwhich,id,bfield.a[id,:,1,:],
                     bfield.j[id,:,0,:],
                     bfield.nx,bfield.nz,bfield.mgtol[id])
      else:
        multigrid3dsolve(iwhich,bfield.nx,bfield.ny,bfield.nz,bfield.nzfull,
                         bfield.dx,bfield.dy,bfield.dz,
                         bfield.a[id,1:-1,1:-1,:],
                         bfield.j[id,:,:,:],
                         bfield.rstar,self.linbend,bfield.bounds,
                         bfield.xmmin,bfield.ymmin,bfield.zmmin,
                         top.zbeam,top.zgrid,
                         bfield.mgparam[id],bfield.mgform[id],
                         bfield.mgiters[id],bfield.mgmaxiters[id],
                         bfield.mgmaxlevels[id],bfield.mgerror[id],
                         bfield.mgtol[id],
                         bfield.downpasses[id],bfield.uppasses[id],
                         bfield.lcndbndy,bfield.laddconductor,
                         bfield.icndbndy,false,
                         self.gridmode,bfield.conductors,
                         self.my_index,self.nslaves,self.izfsslave,self.nzfsslave)

    # --- This is slightly inefficient in some cases, since for example, the
    # --- MG solver already takes care of the longitudinal BC's.
    setaboundaries3d(bfield)

    # --- Now take the curl of A to get B.
    getbfroma3d(bfield)

    # --- If using the analytic form of Btheta, calculate it here.
    if bfield.lanalyticbtheta: getanalyticbtheta(bfield)

    # --- Unscale the current density
    bfield.j[...] = bfield.j/(mu0*eps0)


  ##########################################################################
  # Define the basic plot commands
  def genericpf(self,kw,pffunc):
    kw['solver'] = self
    pffunc(**kw)
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







