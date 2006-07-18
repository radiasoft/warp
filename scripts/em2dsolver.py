"""Class for doing 2 1/2 D electromagnetic solver using code adapted from the
emi code"""
from warp import *

try:
  import psyco
except ImportError:
  pass

##############################################################################
class EM2D(object):
  
  __w3dinputs__ = ['solvergeom','nx','ny','nz',
                   'xmmin','xmmax','ymmin','ymmax','zmmin','zmmax',
                   'bound0','boundnz','boundxy','l2symtry','l4symtry']
  __f3dinputs__ = ['gridmode','mgparam','downpasses','uppasses',
                   'mgmaxiters','mgtol','mgmaxlevels','mgform',
                   'lcndbndy','icndbndy','laddconductor'] 
  __topinputs__ = ['pbound0','pboundnz','pboundxy','efetch',
                   'my_index','nslaves','izfsslave','nzfsslave']
  __emiinputs__ = ['l_copyfields','teta',
                   'esirkepov','l_elaser_out_plane','dt','nmaxw','xgauss']
  __flaginputs__ = {'forcesymmetries':1,'nvect':1,'nbndx':10,'nbndy':10,'rap':1}

  def __init__(self,**kw):

    # --- Save input parameters
    for name in EM2D.__w3dinputs__:
      if name not in self.__dict__:
        #self.__dict__[name] = kw.pop(name,getattr(w3d,name)) # Python2.3
        self.__dict__[name] = kw.get(name,getattr(w3d,name))
      if kw.has_key(name): del kw[name]
    for name in EM2D.__f3dinputs__:
      if name not in self.__dict__:
        #self.__dict__[name] = kw.pop(name,getattr(f3d,name)) # Python2.3
        self.__dict__[name] = kw.get(name,getattr(f3d,name))
      if kw.has_key(name): del kw[name]
    for name in EM2D.__topinputs__:
      if name not in self.__dict__:
        #self.__dict__[name] = kw.pop(name,getattr(top,name)) # Python2.3
        self.__dict__[name] = kw.get(name,getattr(top,name))
      if kw.has_key(name): del kw[name]
    for name in EM2D.__emiinputs__:
      if name not in self.__dict__:
        #self.__dict__[name] = kw.pop(name,getattr(emi,name)) # Python2.3
        self.__dict__[name] = kw.get(name,getattr(emi,name))
      if kw.has_key(name): del kw[name]
    for name,defvalue in EM2D.__flaginputs__.iteritems():
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

    # --- Calculate mesh sizes
    # --- All cases use X
    self.dx = (self.xmmax - self.xmmin)/self.nx
    self.xsymmetryplane = 0.
    self.xmesh = self.xmmin + arange(0,self.nx+1)*self.dx
    self.xmin = self.xmmin
    if self.solvergeom in [w3d.XYgeom]:
      # --- When Y is used, nothing special is done
      self.dy = (self.ymmax - self.ymmin)/self.ny
      self.ymesh = self.ymmin + arange(0,self.ny+1)*self.dy
      self.ysymmetryplane = 0.
      self.ymin = self.ymmin
    if self.solvergeom in [w3d.RZgeom,w3d.XZgeom]:
      # --- When Z is used, set Z quantities, and set Y quantities to the same values.
      # --- Internally, the class always uses x and y. The user interface will depend
      # --- on solvergeom
      self.dz = (self.zmmaxglobal - self.zmminglobal)/self.nzfull
      self.zmesh = self.zmmin + arange(0,self.nz+1)*self.dz
      self.zmin = self.zmmin
      self.ny = self.nz
      self.dy = self.dz
      self.ymesh = self.zmesh
      self.ymin = self.zmin

    # --- Create phi and rho arrays and other arrays.
    self.allocatefieldarrays()

  def __getstate__(self):
    dict = self.__dict__.copy()
    return dict

  def __setstate__(self,dict):
    self.__dict__.update(dict)

  def allocatefieldarrays(self):
    # --- Code transcribed from init_fields
    self.field = NewEMIFIELDType()
    self.bnds = [[Newtype_bnd(),Newtype_bnd()]]
    f = self.field
    f.l_apply_pml=true
    f.nx = self.nx
    f.ny = self.ny
    f.xmin = self.xmin
    f.ymin = self.ymin
    f.rap = self.rap
    f.dx = self.dx
    f.dy = self.dy
    f.dxi = 1./self.dx
    f.dyi = 1./self.dy
    f.npulse = 300 # ???
    f.ntemp = 2*max(self.nx,self.ny)+4
    self.dtm = dt/nmaxw
    if self.l_elaser_out_plane:
      f.cst1 = (1.-dtm/self.dy)/(1.+dtm/self.dy)
      f.cst2 =  2.*dtm/self.dy /(1.+dtm/self.dy)
    else
      f.cst1 = (1.-dtm/self.dx)/(1.+dtm/self.dx)
      f.cst2 =  2.*dtm/self.dx /(1.+dtm/self.dx)

    f.nvect = self.nvect

    f.js = 1 # position of the source

    f.sinteta = sin(self.teta)

    f.n4x      = self.nx+4
    f.n4xp1    = f.n4x+1
    f.n4xp2    = f.n4x+2
    f.n4xp3    = f.n4x+3

    f.n4xf2    = f.n4x*2
    f.n4xf2p1  = f.n4x*2+1
    f.n4xf2p2  = f.n4x*2+2
    f.n4xf2p3  = f.n4x*2+3

    f.n4xf3    = f.n4x*3
    f.n4xf3p1  = f.n4x*3+1
    f.n4xf3p2  = f.n4x*3+2
    f.n4xf3p3  = f.n4x*3+3

    f.ngmax    = (self.nx+4)*(self.ny+3)
    f.ngmaxf2  = (self.nx+4)*(self.ny+3)*2
    f.ngmaxf3  = 3*f.ngmax

    if self.l_copyfields:
      f.nxcopy = f.nx
      f.nycopy = f.ny
    else:
      f.nxcopy = 0
      f.nycopy = 0

    f.allot()

    if self.esirkepov:

      for k in range(1,self.nvect*41+41+1,41):
        if k >= self.nvect*41: break
          f.id1(41+k-1)=(nx+3)+(k%ny)*f.n4x
          f.sd(41+k-1)=0.

      for k in range(1,nvect*19+19+1,19):
        if k >= self.nvect*19: break
          f.id4(19+k-1)=(nx+3)+(k%ny)*f.n4x

      for k in range(1,nvect*17+17+1,17):
        if k >= self.nvect*17: break
          f.id2(17+k-1)=(nx+3)+(k%ny)*f.n4x
          f.sd(17+k-1)=0.

      for k in range(1,nvect*49+49+1,49):
        if k >= self.nvect*49: break
          f.id3(49+k-1)=(nx+3)+(k%ny)*f.n4x

      f.Ex = 0.
      f.Ey = 0.
      f.Ez = 0.
      f.Bx = 0.
      f.By = 0.
      f.Bz = 0.

      f.Rhojxjy = 0.

      f.Bz_in = 0.
      f.Ey_in = 0.
      f.Ex_in = 0.
      f.Ez_in = 0.
      f.By_in = 0.
      f.Bx_in = 0.
      f.pulse = 0.
      f.tpulse = 0.

      if self.l_elaser_out_plane:
        for m in range(0,nx+3+1):
          f.profx[m]=exp(-((m*dx+f.xmin*self.field.dx-0.5*self.field.nx*self.field.dx)/self.xgauss)**2/2.)
      else:
        for m in range(0,ny+2+1):
          f.profy[m]=exp(-((m*dy+f.ymin*self.field.dy-0.5*self.field.ny*self.field.dy)/self.xgauss)**2/2.)

      create_bnd(self.bnds[0], f.nx, f.ny, self.nbndx, self.nbndy, self.dtm, f.dx, f.dy)
      create_bnd(self.bnds[1], f.nx, f.ny, self.nbndx, self.nbndy, self.dtm, f.dx, f.dy)



  def setrho(self,x,y,z,uz,q,w):
    n = len(x)
    if n == 0: return
    setrho3d(self.rho,self.rho,n,x,y,z,top.zgrid,uz,q,w,top.depos,
             self.nx,self.ny,self.nz,self.dx,self.dy,self.dz,
             self.xmmin,self.ymmin,self.zmmin,self.l2symtry,self.l4symtry)

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
      w3d.xfsapi=top.pgroup.xp[w3d.ipminapi-1:w3d.ipminapi-1+w3d.ipapi]
      w3d.yfsapi=top.pgroup.yp[w3d.ipminapi-1:w3d.ipminapi-1+w3d.ipapi]
      w3d.zfsapi=top.pgroup.zp[w3d.ipminapi-1:w3d.ipminapi-1+w3d.ipapi]
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
      EMI.getselfe(self,recalculate=1)

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






