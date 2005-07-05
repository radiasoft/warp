"""Class for doing complete FFT field solve"""
# ToDo:
#  - ???
from warp import *
from lattice import addnewbgrd,addnewbsqgrad
import MA

##############################################################################
class FieldSolver3dBase(object):
  
  __w3dinputs__ = ['nx','ny','nz','nzfull',
                   'xmmin','xmmax','ymmin','ymmax','zmmin','zmmax',
                   'zmminglobal','zmmaxglobal',
                   'bound0','boundnz','boundxy','l2symtry','l4symtry',
                   'solvergeom']
  __f3dinputs__ = []
  __topinputs__ = ['pbound0','pboundnz','pboundxy','efetch','nslaves']
  __flaginputs__ = {'forcesymmetries':1}

  def __init__(self,useselfb=0,**kw):
    self.solvergeom = w3d.XYZgeom
    self.useselfb = useselfb

    # --- Save input parameters
    for name in self.__class__.__w3dinputs__:
      if name not in self.__dict__:
        #self.__dict__[name] = kw.pop(name,getattr(w3d,name)) # Python2.3
        self.__dict__[name] = kw.get(name,getattr(w3d,name))
        if kw.has_key(name): del kw[name]
    for name in self.__class__.__f3dinputs__:
      if name not in self.__dict__:
        #self.__dict__[name] = kw.pop(name,getattr(f3d,name)) # Python2.3
        self.__dict__[name] = kw.get(name,getattr(f3d,name))
        if kw.has_key(name): del kw[name]
    for name in self.__class__.__topinputs__:
      if name not in self.__dict__:
        #self.__dict__[name] = kw.pop(name,getattr(top,name)) # Python2.3
        self.__dict__[name] = kw.get(name,getattr(top,name))
        if kw.has_key(name): del kw[name]
    for name,defvalue in self.__class__.__flaginputs__.iteritems():
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

    # --- Set set parallel related paramters and calculate mesh sizes
    if self.nslaves <= 1:
      self.nzfull = self.nz
      self.zmminglobal = self.zmmin
      self.zmmaxglobal = self.zmmax
    else:
      if self.nzfull == 0:
        self.nzfull = self.nz
        self.zmminglobal = self.zmmin
        self.zmmaxglobal = self.zmmax
        self.nz = self.nzfull/self.nslaves
        self.zmmin = self.zmminglobal + me*self.nz*w3d.dz
        self.zmmax = self.zmminglobal + (me+1)*self.nz*w3d.dz
        if me == self.nslaves-1: self.zmmax = self.zmmaxglobal

    self.dx = (self.xmmax - self.xmmin)/self.nx
    self.dy = (self.ymmax - self.ymmin)/self.ny
    self.dz = (self.zmmaxglobal - self.zmminglobal)/self.nzfull
    self.xsymmetryplane = 0.
    self.ysymmetryplane = 0.

    self.xmesh = self.xmmin + arange(0,self.nx+1)*self.dx
    self.ymesh = self.ymmin + arange(0,self.ny+1)*self.dy
    self.zmesh = self.zmmin + arange(0,self.nz+1)*self.dz

    # --- Create extra variables which are used in various places
    self.nxp = self.nx
    self.nyp = self.ny
    self.nzp = self.nz

    # --- Create phi and rho arrays and other arrays. These are created
    # --- with fortran ordering so no transpose and copy is needed when
    # --- they are passed to fortran.
    self.rho = fzeros((1+self.nx,1+self.ny,1+self.nz),'d')
    self.phi = fzeros((1+self.nx,1+self.ny,3+self.nz),'d')
    self.rstar = fzeros(3+self.nz,'d')
    if self.efetch == 3:
      self.selfe = fzeros((3,1+self.nx,1+self.ny,1+self.nz),'d')
      self.nx_selfe = self.nx
      self.ny_selfe = self.ny
      self.nz_selfe = self.nz
    else:
      self.selfe = 0.
      self.nx_selfe = 0
      self.ny_selfe = 0
      self.nz_selfe = 0
    self.rhop = self.rho
    self.phip = self.phi
    if self.useselfb:
      self.bx = fzeros((1+self.nx,1+self.ny,3+self.nz),'d')
      self.by = fzeros((1+self.nx,1+self.ny,3+self.nz),'d')
      self.bz = fzeros((1+self.nx,1+self.ny,3+self.nz),'d')

    # --- Create a conductor object, which by default is empty.
    self.conductors = None
    self.conductorlist = []

    # --- At the start, assume that there are no bends. This is corrected
    # --- in the solve method when there are bends.
    self.linbend = false

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

  def fetchbfrompositions(self,x,y,z,bx,by,bz):
    n = len(x)
    if n == 0: return
    getgrid3d(n,x,y,z,bx,self.nx,self.ny,self.nz,self.bx[:,:,1:-1],
              self.xmmin,self.xmmax,self.ymmin,self.ymmax,
              self.zmmin+top.zgridprv,self.zmmax+top.zgridprv,
              self.l2symtry,self.l4symtry)
    getgrid3d(n,x,y,z,by,self.nx,self.ny,self.nz,self.by[:,:,1:-1],
              self.xmmin,self.xmmax,self.ymmin,self.ymmax,
              self.zmmin+top.zgridprv,self.zmmax+top.zgridprv,
              self.l2symtry,self.l4symtry)
    getgrid3d(n,x,y,z,bz,self.nx,self.ny,self.nz,self.bz[:,:,1:-1],
              self.xmmin,self.xmmax,self.ymmin,self.ymmax,
              self.zmmin+top.zgridprv,self.zmmax+top.zgridprv,
              self.l2symtry,self.l4symtry)

  def fetchphifrompositions(self,x,y,z,phi):
    n = len(x)
    getgrid3d(n,x,y,z,phi,self.nx,self.ny,self.nz,self.phi[:,:,1:-1],
              self.xmmin,self.xmmax,self.ymmin,self.ymmax,self.zmmin,self.zmmax,
              self.l2symtry,self.l4symtry)

  def loadrho(self,ins_i=-1,nps_i=-1,is_i=-1,lzero=true):
    if lzero: self.rho[...] = 0.
    for i,n,q,w in zip(top.ins-1,top.nps,top.sq,top.sw):
      self.setrho(top.xp[i:i+n],top.yp[i:i+n],top.zp[i:i+n],top.uzp[i:i+n],q,w)
    self.makerhoperiodic()
    self.getrhoforfieldsolve()

  def makerhoperiodic(self):
    trho = transpose(self.rho)
    if self.pbounds[0] == 2 or self.pbounds[1] == 2:
      trho[:,:,0] = trho[:,:,0] + trho[:,:,-1]
      trho[:,:,-1] = trho[:,:,0]
    if self.pbounds[2] == 2 or self.pbounds[3] == 2:
      trho[:,0,:] = trho[:,0,:] + trho[:,-1,:]
      trho[:,-1,:] = trho[:,0,:]
    if self.pbounds[0] == 1 and not self.l4symtry:
       trho[:,:,0] = 2.*trho[:,:,0]
    if self.pbounds[1] == 1: trho[:,:,-1] = 2.*trho[:,:,-1]
    if self.pbounds[2] == 1 and not (self.l2symtry or self.l4symtry):
       trho[:,0,:] = 2.*trho[:,0,:]
    if self.pbounds[3] == 1: trho[:,-1,:] = 2.*trho[:,-1,:]
    if self.pbounds[4] == 2 or self.pbounds[5] == 2:
      if self.nslaves > 1:
        self.makerhoperiodic_parallel()
      else:
        trho[0,:,:] = trho[0,:,:] + trho[-1,:,:]
        trho[-1,:,:] = trho[0,:,:]
    if self.pbounds[4] == 1: trho[0,:,:] = 2.*trho[0,:,:]
    if self.pbounds[5] == 1: trho[-1,:,:] = 2.*trho[-1,:,:]

  def getrhoforfieldsolve(self):
    if self.nslaves > 1:
      getrhoforfieldsolve3d(self.nx,self.ny,self.nz,self.rho,
                            self.nx,self.ny,self.nz,self.rhop,
                            me,self.nslaves,
                          top.izfsslave,top.nzfsslave,top.izpslave,top.nzpslave)

  def makerhoperiodic_parallel(self):
    tag = 70
    if me == self.nslaves-1:
      request = mpi.isend(self.rho[:,:,nz],0,tag)
      self.rho[:,:,nz],status = mpi.recv(0,tag)
    elif me == 0:
      rhotemp,status = mpi.recv(self.nslaves-1,tag)
      self.rho[:,:,0] = self.rho[:,:,0] + rhotemp
      request = mpi.isend(self.rho[:,:,0],self.nslaves-1,tag)
    if me == 0 or me == self.nslaves-1:
      status = request.wait()

  def fetche(self):
    self.fetchefrompositions(w3d.xfsapi,w3d.yfsapi,w3d.zfsapi,
                             w3d.exfsapi,w3d.eyfsapi,w3d.ezfsapi)

  def fetchb(self):
    self.fetchbfrompositions(w3d.xfsapi,w3d.yfsapi,w3d.zfsapi,
                             w3d.bxfsapi,w3d.byfsapi,w3d.bzfsapi)

  def fetchphi(self):
    self.fetchphifrompositions(w3d.xfsapi,w3d.yfsapi,w3d.zfsapi,w3d.phifsapi)

  def getselfe(self,recalculate=0):
    if type(self.selfe) != ArrayType:
      self.selfe = fzeros((3,1+self.nx,1+self.ny,1+self.nz),'d')
      recalculate = 1
    if recalculate:
      getselfe3d(self.phi,self.nx,self.ny,self.nz,
                 self.selfe,self.nx,self.ny,self.nz,self.dx,self.dy,self.dz,
                 self.bounds[0],self.bounds[1],self.bounds[2],self.bounds[3])
    return self.selfe

  def solve(self,iwhich=0):
    raise "solve must be implemented"

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

  def __getstate__(self):
    """
Check whether this instance is the registered solver so that upon unpickling
it knows whether to re-register itself.
    """
    dict = self.__dict__.copy()
    if self is getregisteredsolver():
      dict['iamtheregisteredsolver'] = 1
    else:
      dict['iamtheregisteredsolver'] = 0
    return dict

  def __setstate__(self,dict):
    self.__dict__.update(dict)
    if self.iamtheregisteredsolver:
      registersolver(self)
      if self.nx == w3d.nx and self.ny == w3d.ny and self.nz == w3d.nz:
        w3d.rho = self.rho
        w3d.phi = self.phi
        w3d.nxp = self.nx
        w3d.nyp = self.ny
        w3d.nzp = self.nz
        w3d.rhop = self.rho
        w3d.phip = self.phi


##############################################################################
class FFTSolver2d(FieldSolver3dBase):
  __w3dinputs__ = FieldSolver3dBase.__w3dinputs__ + ['filt']
  def __init__(self,**kw):
    FieldSolver3dBase.__init__(self,**kw)

    # --- Allocate the ksq and att arrays
    self.kxsq = zeros(self.nx,'d')
    self.kysq = zeros(self.ny,'d')
    self.kzsq = zeros(self.nzfull+1,'d')
    self.attx = zeros(self.nx,'d')
    self.atty = zeros(self.ny,'d')
    self.attz = zeros(self.nzfull+1,'d')

    # --- Allocate work arrays
    nmxyz = max(self.nx,self.ny,self.nz)
    nmxy = max(self.nx,self.ny)
    self.scratch = fzeros((nmxyz,nmxy),'d')
    self.xywork = fzeros((1+self.nx,1+self.ny),'d')
    self.zwork = 0. # --- zwork isn't used during the solve

    # --- Initialize itself
    self.vp3d(1)

    # --- Clear out the z part so it doesn't contribute
    self.kzsq[:] = 0.
    self.attz[:] = 1.

  def vp3d(self,iwhich):
    lx = self.xmmax - self.xmmin
    ly = self.ymmax - self.ymmin
    lz = 1.
    vpois3d(iwhich,self.phi[:,:,1:-1],self.phi[:,:,1:-1],
            self.kxsq,self.kysq,self.kzsq,
            self.attx,self.atty,self.attz,self.filt,lx,ly,lz,
            self.nx,self.ny,self.nz,self.nz,
            self.scratch,self.xywork,self.zwork,0,self.l2symtry,self.l4symtry,
            self.bound0,self.boundnz,self.boundxy)

  def solve(self,iwhich=0):

    if iwhich == 1:
      # --- This has already been initialized so don't do anything.
      # --- If vp3d(1) were called, it would mess up kzsq.
      return

    ikxmin = 1
    if self.l4symtry: ikxmin = 0
    ikymin = 1
    if self.l2symtry or self.l4symtry: ikymin = 0

    # --- This is much faster when transposed to C ordering.
    #self.phi[:,:,1:-1] = self.rho
    tphi = transpose(self.phi[:,:,1:-1])
    trho = transpose(self.rho)
    tphi[:,:,:] = trho

    # --- The temporary setting of bound0 forces vp3d to loop from
    # --- iz=0,nz
    s0 = self.bound0
    self.bound0 = dirichlet
    self.vp3d(12)
    attenuate(self.nx,self.ny,self.nz,self.phi[:,:,1:-1],
              self.attx,self.atty,self.attz,ikxmin,ikymin,1,1,0)
    rhotophi(self.nx,self.ny,self.nz,self.phi[:,:,1:-1],
             self.kxsq,self.kysq,self.kzsq,ikxmin,ikymin,1,1,0)
    self.vp3d(13)
    self.bound0 = s0

    tphi = transpose(self.phi)
    if self.bound0 == neumann: tphi[0,:,:] = tphi[1,:,:]
    if self.bound0 == periodic: tphi[0,:,:] = tphi[-3,:,:]
    if self.boundnz == neumann: tphi[-1,:,:] = tphi[-3,:,:]
    if self.boundnz == periodic: tphi[-2:,:,:] = tphi[:2,:,:]

##############################################################################
class RelativisticFFTSolver2d(object):
  """
Transverse 2-D field solver, ignores self Ez and Bz.
  beamspecies=0:
  backgroundspecies=[]:
  ignorebeamez=1:
  ignorebackgroundez=1:
  useselfb=0: When true, use B fields directly rather than just E/gamma**2
  """
  def __init__(self,beamspecies=0,backgroundspecies=[],
                    ignorebeamez=1,ignorebackgroundez=1,
                    useselfb=0,usebeamb=0,useselfbsqgrad=0,**kw):

    self.beamspecies = beamspecies
    self.backgroundspecies = backgroundspecies
    self.ignorebeamez = ignorebeamez
    self.ignorebackgroundez = ignorebackgroundez
    self.useselfb = useselfb
    self.usebeamb = usebeamb
    self.useselfbsqgrad = useselfbsqgrad and (max(w3d.interpdk) > 0.)

    # --- Create separate solvers for the two particle types and initialize them
    self.beamsolver = FFTSolver2d(useselfb=self.useselfb or self.usebeamb,**kw)
    if len(self.backgroundspecies) > 0:
      self.backgroundsolver = FFTSolver2d(**kw)

    if self.useselfb:
      top.lcallfetchb = true

    if self.useselfbsqgrad:
      w3d.dx = (w3d.xmmax - w3d.xmmin)/w3d.nx
      w3d.dy = (w3d.ymmax - w3d.ymmin)/w3d.ny

      # --- Create bsqgrad element the same size as the field grid.
      self.bsqgradid = addnewbsqgrad(zs=w3d.zmmin,ze=w3d.zmmax,
                                     xs=w3d.xmmin,ys=w3d.ymmin,
                                     dx=w3d.dx,dy=w3d.dy,
                                     nx=w3d.nx,ny=w3d.ny,nz=w3d.nz)
      # --- Turn on the use of bsqgrad
      w3d.igradb = 1

      resetlat()
      setlatt()

    # --- Make the w3d arrays point to the ones in the beamsolver.
    # --- That was an arbitrary choice.
    if (self.beamsolver.nx == w3d.nx and
        self.beamsolver.ny == w3d.ny and
        self.beamsolver.nz == w3d.nz):
      w3d.rho = self.beamsolver.rho
      w3d.phi = self.beamsolver.phi
      w3d.nxp = self.beamsolver.nx
      w3d.nyp = self.beamsolver.ny
      w3d.nzp = self.beamsolver.nz
      w3d.rhop = self.beamsolver.rho
      w3d.phip = self.beamsolver.phi

  def loadrho(self,ins_i=-1,nps_i=-1,is_i=-1,lzero=true):
    if lzero:
      self.beamsolver.rho[...] = 0.
      if len(self.backgroundspecies) > 0:
        self.backgroundsolver.rho[...] = 0.

    js = self.beamspecies
    i,n,q,w = top.ins[js]-1,top.nps[js],top.sq[js],top.sw[js]
    if n > 0:
      self.beamsolver.setrho(top.xp[i:i+n],top.yp[i:i+n],top.zp[i:i+n],
                             top.uzp[i:i+n],q,w)
    self.beamsolver.makerhoperiodic()
    self.beamsolver.getrhoforfieldsolve()

    if len(self.backgroundspecies) > 0:
      for js in self.backgroundspecies:
        i,n,q,w = top.ins[js]-1,top.nps[js],top.sq[js],top.sw[js]
        if n > 0:
          self.backgroundsolver.setrho(top.xp[i:i+n],top.yp[i:i+n],
                                       top.zp[i:i+n],top.uzp[i:i+n],q,w)
      self.backgroundsolver.makerhoperiodic()
      self.backgroundsolver.getrhoforfieldsolve()


  def calcbeambfield(self):
    # --- Fill in the bgrd arrays
    # --- Note that array operations are done transposed for efficiency.
    phib = transpose(self.beamsolver.phi)

    vz = top.vbeam_s[self.beamspecies]
    Az = (+vz*eps0*mu0)*phib
    bx = transpose(self.beamsolver.bx)
    by = transpose(self.beamsolver.by)
    #bx[:,1:-1,:] = +(Az[:,2:,:] - Az[:,:-2,:])/(2.*w3d.dy)
    #by[:,:,1:-1] = -(Az[:,:,2:] - Az[:,:,:-2])/(2.*w3d.dx)

    bx[:,1:-1,:] = Az[:,2:,:]
    subtract(bx[:,1:-1,:],Az[:,:-2,:],bx[:,1:-1,:])
    multiply(bx[:,1:-1,:],1./(2.*w3d.dy),bx[:,1:-1,:])

    by[:,:,1:-1] = Az[:,:,:-2]
    subtract(by[:,:,1:-1],Az[:,:,2:],by[:,:,1:-1])
    multiply(by[:,:,1:-1],1./(2.*w3d.dx),by[:,:,1:-1])

  def calcbeambsqgrad(self):
    if self.useselfbsqgrad:
      bx = transpose(self.beamsolver.bx)
      by = transpose(self.beamsolver.by)

      bsq = bx**2 + by**2
      id = self.bsqgradid - 1
      bsqgrad = transpose(top.bsqgrad[:,:,:,:,id])
      bsqgrad[1:-1,:,:,0] = (bsq[2:,:,:] - bsq[:-2,:,:])/(2.*w3d.dz)
      bsqgrad[:,:,1:-1,1] = (bsq[:,:,2:] - bsq[:,:,:-2])/(2.*w3d.dx)
      bsqgrad[:,1:-1,:,2] = (bsq[:,2:,:] - bsq[:,:-2,:])/(2.*w3d.dy)

  def calcbfieldsumphi(self):
    # --- Calculate the B field from the phi from the beam and sum the phi's
    # --- of the beam and background
    # --- Also, calculate bsqgrad if needed
      
    self.calcbeambfield()
    self.calcbeambsqgrad()

    # --- Sum the phi's

    # --- The operations below are faster when the arrays are in C ordering
    if len(self.backgroundspecies) > 0:
      tbeackgroundphi = transpose(self.backgroundsolver.phi)
      tbeamphi = transpose(self.beamsolver.phi)

      add(tbeamphi,tbeackgroundphi,tbeamphi)
      tbeackgroundphi[...] = tbeamphi

  def gammacorrectsumphi(self):
    # --- The beam gets the background phi plus 1/gamma**2 * the beam phi
    # --- The background gets the sum of the beam and background phi's

    # --- Calculate the beam's beam field before applying the gamma correction
    if self.usebeamb:
      self.calcbeambfield()
      self.calcbeambsqgrad()

    # --- First make a copy of the beam's phi since it is overwritten below
    # --- The operations below are faster when the arrays are in C ordering
    if len(self.backgroundspecies) > 0:
      tbeamphitemp = transpose(self.beamsolver.phi.copy())
      tbeackgroundphi = transpose(self.backgroundsolver.phi)

    tbeamphi = transpose(self.beamsolver.phi)

    # --- Calculate gammabar, avoiding roundoff from subtraction of
    # --- similar numbers
    js = self.beamspecies
    ke = jperev*top.ekin_s[js]/dvnz(top.aion_s[js]*amu*clight**2)
    gammabar = 1. + ke

    # --- Now calculate the beam's phi in place
    multiply(tbeamphi,1./gammabar**2,tbeamphi)
    if len(self.backgroundspecies) > 0:
      add(tbeamphi,tbeackgroundphi,tbeamphi)

    if len(self.backgroundspecies) > 0:
      # --- Sum the phi's for the background
      add(tbeamphitemp,tbeackgroundphi,tbeackgroundphi)

  def solve(self):

    self.beamsolver.solve()
    if len(self.backgroundspecies) > 0:
      self.backgroundsolver.solve()

    if self.useselfb:
      self.calcbfieldsumphi()
    else:
      self.gammacorrectsumphi()

  def fetche(self):

    if w3d.isfsapi-1 == self.beamspecies:
      self.beamsolver.fetchefrompositions(w3d.xfsapi,w3d.yfsapi,w3d.zfsapi,
                                          w3d.exfsapi,w3d.eyfsapi,w3d.ezfsapi)
      if self.ignorebeamez: w3d.ezfsapi = 0.
    else:
      self.backgroundsolver.fetchefrompositions(w3d.xfsapi,w3d.yfsapi,
                                                w3d.zfsapi,
                                                w3d.exfsapi,w3d.eyfsapi,
                                                w3d.ezfsapi)
      if self.ignorebackgroundez: w3d.ezfsapi = 0.

  def fetchb(self):

    if w3d.isfsapi-1 == self.beamspecies:
      if self.useselfb:
        self.beamsolver.fetchbfrompositions(w3d.xfsapi,w3d.yfsapi,w3d.zfsapi,
                                            w3d.bxfsapi,w3d.byfsapi,w3d.bzfsapi)
      else:
        w3d.bxfsapi = 0.
        w3d.byfsapi = 0.
        w3d.bzfsapi = 0.
    else:
      if self.usebeamb or self.useselfb:
        self.beamsolver.fetchbfrompositions(w3d.xfsapi,w3d.yfsapi,w3d.zfsapi,
                                            w3d.bxfsapi,w3d.byfsapi,w3d.bzfsapi)
      else:
        w3d.bxfsapi = 0.
        w3d.byfsapi = 0.
        w3d.bzfsapi = 0.

      if self.ignorebackgroundez: w3d.bzfsapi = 0.

  def fetchphi(self):
    # --- For present purposes, the results from this routine are never
    # --- used. Also, presently, there is no way of telling which species
    # --- the data is for.
    w3d.phifsapi[:] = 0.

  def __getstate__(self):
    """
Check whether this instance is the registered solver so that upon unpickling
it knows whether to re-register itself.
    """
    dict = self.__dict__.copy()
    if self is getregisteredsolver():
      dict['iamtheregisteredsolver'] = 1
    else:
      dict['iamtheregisteredsolver'] = 0
    return dict

  def __setstate__(self,dict):
    self.__dict__.update(dict)
    if self.iamtheregisteredsolver:
      registersolver(self)
      if (self.beamsolver.nx == w3d.nx and
          self.beamsolver.ny == w3d.ny and
          self.beamsolver.nz == w3d.nz):
        w3d.rho = self.beamsolver.rho
        w3d.phi = self.beamsolver.phi
        w3d.nxp = self.beamsolver.nx
        w3d.nyp = self.beamsolver.ny
        w3d.nzp = self.beamsolver.nz
        w3d.rhop = self.beamsolver.rho
        w3d.phip = self.beamsolver.phi

  ##########################################################################
  # Define the basic plot commands
  def genericpf(self,kw,pffunc):
    kw['conductors'] = None
    kw['phicolor'] = blue
    kw['solver'] = self.beamsolver
    pffunc(**kw)
    if len(self.backgroundspecies) > 0:
      kw['phicolor'] = green
      kw['solver'] = self.backgroundsolver
      pffunc(**kw)
  def pfxy(self,**kw): self.genericpf(kw,pfxy)
  def pfzx(self,**kw): self.genericpf(kw,pfzx)
  def pfzy(self,**kw): self.genericpf(kw,pfzy)
  def pfxyg(self,**kw): self.genericpf(kw,pfxyg)
  def pfzxg(self,**kw): self.genericpf(kw,pfzxg)
  def pfzyg(self,**kw): self.genericpf(kw,pfzyg)

