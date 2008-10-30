"""Class for doing 2 1/2 D electromagnetic solver using code adapted from the
emi code"""
from warp import *
from getzmom import * 
import operator

try:
  import psyco
except ImportError:
  pass

##############################################################################
class EM2D(object):
  
  __w3dinputs__ = ['solvergeom','nx','ny','nzlocal',
                   'xmmin','xmmax','ymmin','ymmax','zmminlocal','zmmaxlocal',
                   'bound0','boundnz','boundxy','l2symtry','l4symtry',
                   'solvergeom']
  __em2dinputs__ = ['l_onegrid','l_copyfields','l_moving_window',
                    'tmin_moving_main_window','l_smoothdensity',
                    'l_elaser_out_plane','ndelta_t',
                   'ntamp_scatter','ntamp_gather']
  __flaginputs__ = {'l_apply_pml':true,'nbndx':10,'nbndy':10,
                    'l_particles_weight':false,'l_usecoeffs':false,
                    'l_verbose':1,
                    'laser_amplitude':1.,'laser_profile':None,
                    'laser_gauss_width':None,'laser_angle':0.,
                    'laser_wavelength':None,'laser_wavenumber':None,
                    'laser_frequency':None,'laser_source':2,
                    'laser_focus':None,'laser_focus_velocity':0.,
                    'density_1d':False,'nfield_subcycle':1,
                    'autoset_timestep':true,'dtcoef':0.99}

  def __init__(self,**kw):
#    top.allspecl = true
    top.lcallfetchb = true
    #top.bfstype = 12
    top.lcallfetchb = true
    top.lgridqnt = true
    self.zgridprv=top.zgrid
    # --- Make sure the refinement is turned off
#    em2d.l_onegrid = true

    # --- Save input parameters
    for name in EM2D.__w3dinputs__:
      if name not in self.__dict__:
        #self.__dict__[name] = kw.pop(name,getattr(w3d,name)) # Python2.3
        self.__dict__[name] = kw.get(name,getattr(w3d,name))
      if kw.has_key(name): del kw[name]
    for name in EM2D.__em2dinputs__:
      if name not in self.__dict__:
        #self.__dict__[name] = kw.pop(name,getattr(em2d,name)) # Python2.3
        self.__dict__[name] = kw.get(name,getattr(em2d,name))
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
        self.bounds = zeros(6,'l')
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

    # --- If there are any remaning keyword arguments, raise an error.
    assert len(kw.keys()) == 0,"Bad keyword arguemnts %s"%kw.keys()

    assert self.solvergeom == w3d.XZgeom or not self.l_moving_window,"The moving window can only be used with XZ geometry"

    # --- Calculate mesh sizes
    if self.solvergeom in [w3d.XYgeom]:
      # --- When Y is used, nothing special is done
      self.dx = (self.xmmax - self.xmmin)/self.nx
      self.xmesh = self.xmmin + arange(0,self.nx+1)*self.dx
      self.dy = (self.ymmax - self.ymmin)/self.ny
      self.ymesh = self.ymmin + arange(0,self.ny+1)*self.dy

    if self.solvergeom in [w3d.XZgeom]:
      # --- When Z is used, set Z quantities, and set Y quantities to the same
      # --- values. Internally, the class always uses x and y. The user
      # --- interface will depend on solvergeom
      self.ny = self.nx
      self.ymmin = self.xmmin
      self.ymmax = self.xmmax
      self.dy = (self.ymmax - self.ymmin)/self.ny
      self.ymesh = self.ymmin + arange(0,self.ny+1)*self.dy
      self.dz = (self.zmmaxlocal - self.zmminlocal)/self.nzlocal
      self.zmesh = self.zmminlocal + arange(0,self.nzlocal+1)*self.dz
      self.zmeshlocal = self.zmminlocal + arange(0,self.nzlocal+1)*self.dz
      self.nx = self.nzlocal
      self.dx = self.dz
      self.xmesh = self.zmeshlocal
      self.xmmin = self.zmminlocal
      self.xmmax = self.zmmaxlocal
      self.bounds[0] = self.bounds[4]
      self.bounds[1] = self.bounds[5]
#     self.ndelta_t = top.vbeam*top.dt/w3d.dz
#     em2d.ndelta_t = self.ndelta_t
#     assert (self.ndelta_t - top.vbeam*top.dt/w3d.dz) < 1.e-6,"The input ndelta_t is not commensurate with the values of top.vbeam, top.dt, and w3d.dz"

    # --- set time step as a fraction of Courant condition
    # --- also set self.nfield_subcycle if top.dt over Courant condition times dtcoef
    if self.autoset_timestep:
      dt=self.dtcoef/(clight*sqrt(1./self.dx**2+1./self.dy**2))
      if top.dt==0.:
        top.dt=dt
      else:
        if top.dt<>dt:
          self.nfield_subcycle=int(top.dt/dt)+1
    self.dtinit = top.dt

    # --- Create field and source arrays and other arrays.
    self.allocatefieldarrays()

    # --- Handle laser inputs
    self.setuplaser()
    
    # --- set xyz oldpid, if needed
    if top.xoldpid==0:top.xoldpid=nextpid()
    if top.yoldpid==0:top.yoldpid=nextpid()
    if top.zoldpid==0:top.zoldpid=nextpid()

    # --- register solver
    registersolver(self)

  def allocatefieldarrays(self):
    # --- Code transcribed from init_fields
    self.field = EM2D_FIELDtype()
    self.field.l_usecoeffs = self.l_usecoeffs
    self.field.nxf = 0
    self.field.nyf = 0
    init_fields(self.field,self.nx,self.ny,self.nbndx,self.nbndx,top.dt/self.nfield_subcycle,
                self.dx,self.dy,clight,mu0,
                self.xmmin,self.ymmin,1,
                self.bounds[0],self.bounds[1],self.bounds[2],self.bounds[3])
    self.field.js = self.laser_source
    self.fpatches  = []


  def addpatch(self,ixpatch,iypatch,nxpatch,nypatch,rap):
    # --- Initializes refinement patch
    em2d.l_onegrid = false
    self.l_onegrid = em2d.l_onegrid
    self.fpatches.append([EM2D_FIELDtype(),EM2D_FIELDtype()])
    self.fpatchcoarse = self.fpatches[-1][0]
    self.fpatchfine   = self.fpatches[-1][1]
    self.fpatchfine.rap = rap
    xlb = xrb = ylb = yrb = absorb
    if ixpatch==0:xlb = self.bounds[0]
    if iypatch==0:ylb = self.bounds[1]
    if ixpatch+nxpatch==self.nx:xrb = self.bounds[2]
    if iypatch+nypatch==self.ny:yrb = self.bounds[3]
    init_fields(self.fpatchcoarse,nxpatch,nypatch,self.nbndx,self.nbndx,top.dt/self.nfield_subcycle,
                self.dx,self.dy,clight,mu0,
                self.xmmin+ixpatch*self.dx,self.ymmin+iypatch*self.dy,1,
                xlb,ylb,xrb,yrb)
    init_fields(self.fpatchfine,nxpatch*rap,nypatch*rap,self.nbndx,self.nbndx,top.dt/self.nfield_subcycle,
                self.dx/rap,self.dy/rap,clight,mu0,
                self.xmmin+ixpatch*self.dx,self.ymmin+iypatch*self.dy,rap,
                xlb,ylb,xrb,yrb)
    self.fpatchcoarse.js = self.laser_source
    self.fpatchfine.js = self.laser_source*rap
    self.fpatchfine.xminpatch_scatter = self.fpatchfine.xmin+em2d.ntamp_scatter*rap*self.fpatchfine.dx
    self.fpatchfine.xmaxpatch_scatter = self.fpatchfine.xmax-em2d.ntamp_scatter*rap*self.fpatchfine.dx
    self.fpatchfine.yminpatch_scatter = self.fpatchfine.ymin+em2d.ntamp_scatter*rap*self.fpatchfine.dy
    self.fpatchfine.ymaxpatch_scatter = self.fpatchfine.ymax-em2d.ntamp_scatter*rap*self.fpatchfine.dy
    self.fpatchfine.xminpatch_gather = self.fpatchfine.xmin+em2d.ntamp_gather*rap*self.fpatchfine.dx
    self.fpatchfine.xmaxpatch_gather = self.fpatchfine.xmax-em2d.ntamp_gather*rap*self.fpatchfine.dx
    self.fpatchfine.yminpatch_gather = self.fpatchfine.ymin+em2d.ntamp_gather*rap*self.fpatchfine.dy
    self.fpatchfine.ymaxpatch_gather = self.fpatchfine.ymax-em2d.ntamp_gather*rap*self.fpatchfine.dy
    self.fpatchfine.nxfsum = self.fpatchfine.nx
    self.fpatchfine.nyfsum = self.fpatchfine.ny
    self.fpatchfine.gchange()
    self.setuplaser_profile(self.fpatchcoarse)
    self.setuplaser_profile(self.fpatchfine)

  def setuplaser(self):
    if self.laser_profile is not None:
      if self.laser_frequency is None:
        if self.laser_wavenumber is not None:
          self.laser_wavelength = clight*self.laser_wavenumber
        elif self.laser_wavelength is not None:
          self.laser_wavelength = 2.*pi*clight/self.laser_wavelength
      assert self.laser_frequency is not None,\
             "One of the frequency, wavenumber, or wavelength must be given"

    # --- Check if laser_amplitude is a function, table, or constant
    self.laser_amplitude_func = None
    self.laser_amplitude_table = None
    if operator.isSequenceType(self.laser_amplitude):
      assert len(self.laser_amplitude.shape) == 2 and \
             self.laser_amplitude.shape[1] == 2,\
             "The laser_amplitude table is not formatted properly"
      self.laser_amplitude_table = self.laser_amplitude
      self.laser_amplitude_table_i = -1
    elif callable(self.laser_amplitude):
      self.laser_amplitude_func = self.laser_amplitude
      
    self.setuplaser_profile(self.field)

  def setuplaser_profile(self,f):
    # --- Check if laser_profile has a type, is a function, or a table
    self.laser_profile_func = None
    # --- disable laser emission on processors id>0 
    if me>0:return
    if self.laser_profile == 'gaussian':
      assert self.laser_gauss_width is not None,\
             "For a gaussian laser, the width must be specified using laser_gauss_width"
#      xx = arange(self.nx+4)*f.dx+f.xmin*f.dx - 0.5*f.nx*f.dx
#      self.laser_profile = exp(-(xx/self.laser_gauss_width)**2/2.)
      yy = arange(f.ny+4)*f.dy+f.ymin*f.dy - 0.5*f.ny*f.dy
      f.laser_profile = exp(-0.5*(yy/self.laser_gauss_width)**2)
    elif operator.isSequenceType(self.laser_profile):
      assert len(f.laser_profile) == f.ny+4,"The specified profile must be of length ny+4"
    elif callable(self.laser_profile):
      self.laser_profile_func = self.laser_profile

  def transformparticles(self,x,y,z,ux=None,uy=None,uz=None,xold=None,yold=None,zold=None):
    if self.solvergeom == w3d.XYgeom:
      return x,y,ux,uy,uz,xold,yold
    elif self.solvergeom == w3d.XZgeom:
      return z,x,uz,ux,uy,zold,xold

  def transformfields(self,fx,fy,fz):
    if self.solvergeom == w3d.XYgeom:
      return fx,fy,fz
    elif self.solvergeom == w3d.XZgeom:
      if self.l_moving_window and not self.l_elaser_out_plane:
        return fz,fx,fy
      else:
#        return fx,fz,fy
        return fz,fx,fy

  def setj(self,x,y,xold,yold,uz,gaminv,q,w,wfact,dt):
    n = len(x)
    if n == 0: return
    if wfact is None:
      wfact = zeros(n,'d')
      l_particles_weight = false
    else:
      l_particles_weight = true
    if self.l_onegrid:
      depose_current_em2d(n,x,y,xold,yold,uz,gaminv,wfact,q*w,dt,l_particles_weight,self.field,self.field)
    else:
      depose_current_em2d(n,x,y,xold,yold,uz,gaminv,wfact,q*w,dt,l_particles_weight,self.field,self.fpatchfine)
    
  def setjpy(self,x,y,ux,uy,uz,gaminv,q,w):
    n = len(x)
    if n == 0: return
    wtmp = zeros(n,'d')
    if self.l_onegrid:
      em2d_depose_jxjy_esirkepov_linear_serial(self.field.J,n,x,y,ux,uy,uz,
             gaminv,wtmp,q*w,self.field.xmin,self.field.ymin,top.dt,
             self.field.dx,self.field.dy,self.field.nx,self.field.ny,
             self.l_particles_weight)
    else:
      inpatch = where((x>self.fpatchfine.xminpatch_scatter) & \
                      (x<self.fpatchfine.xmaxpatch_scatter) & \
                      (y>self.fpatchfine.yminpatch_scatter) & \
                      (y<self.fpatchfine.ymaxpatch_scatter), 1, 0)
      ii = arange(n)
      iin  = compress(inpatch,ii); iout = compress(1-inpatch,ii)
      nin = len(iin);              nout = len(iout)
      if nin>0:
        xin = take(x,iin)
        yin = take(y,iin)
        uxin = take(ux,iin)
        uyin = take(uy,iin)
        uzin = take(uz,iin)
        giin = take(gaminv,iin)
        field = self.fpatchfine
        wtmp = zeros(nin,'d')
        print min(xin),max(xin),field.xmin,field.xmax
        print min(yin),max(yin),field.ymin,field.ymax
        em2d_depose_jxjy_esirkepov_linear_serial(field.J,nin,xin,yin,uxin,uyin,uzin,
               giin,wtmp,q*w,field.xmin,field.ymin,top.dt,
               field.dx,field.dy,field.nx,field.ny,
               self.l_particles_weight)
      if nout>0:
        xout = take(x,iout)
        yout = take(y,iout)
        uxout = take(ux,iout)
        uyout = take(uy,iout)
        uzout = take(uz,iout)
        giout = take(gaminv,iout)
        field = self.field
        wtmp = zeros(nout,'d')
        print min(xout),max(xout),field.xmin,field.xmax
        print min(yout),max(yout),field.ymin,field.ymax
        em2d_depose_jxjy_esirkepov_linear_serial(field.J,nout,xout,yout,uxout,uyout,uzout,
               giout,wtmp,q*w,field.xmin,field.ymin,top.dt,
               field.dx,field.dy,field.nx,field.ny,
               self.l_particles_weight)
      
  def fetchefrompositions(self,x,y,ex,ey,ez):
    n = len(x)
    if n == 0: return
    self.fetchffrompositions(n,x,y,ex,ey,ez,'e')

  def fetchbfrompositions(self,x,y,bx,by,bz):
    n = len(x)
    if n == 0: return
    self.fetchffrompositions(n,x,y,bx,by,bz,'b')
    
  def fetchffrompositions(self,n,x,y,fx,fy,fz,which):  
    if which=='e':iwhich=1
    if which=='b':iwhich=2
    if self.l_onegrid:
      getf_em2d(n,x,y,fx,fy,fz,self.field,self.field,iwhich)
    else:
      getf_em2d(n,x,y,fx,fy,fz,self.field,self.fpatchfine,iwhich)

  def fetchffrompositionspy(self,n,x,y,fx,fy,fz,which='e'):  
    if which=='e':
      fxg = self.field.Ex
      fyg = self.field.Ey
      fzg = self.field.Ez
    if which=='b':
      fxg = self.field.Bx
      fyg = self.field.By
      fzg = self.field.Bz
    if self.l_onegrid:
      em2d_getf2d_linear_serial(n,x,y,fx,fy,fz,
                                self.field.xmin,self.field.ymin,
                                self.field.dx,self.field.dy,
                                self.field.nx,self.field.ny,
                                fxg,fyg,fzg)
    else:
      inpatch = where((x>self.fpatchfine.xminpatch_gather) & 
                      (x<self.fpatchfine.xmaxpatch_gather) &
                      (y>self.fpatchfine.yminpatch_gather) & 
                      (y<self.fpatchfine.ymaxpatch_gather), 1, 0)
      ii = arange(n)
      iin  = compress(inpatch,ii); iout = compress(1-inpatch,ii)
      nin = len(iin);              nout = len(iout)
      if nin>0:
        xin = take(x,iin)
        yin = take(y,iin)
        fxin = take(fx,iin)
        fyin = take(fy,iin)
        fzin = take(fz,iin)
        field = self.fpatchfine
        if which=='e':
          ffxg = field.Exfsum
          ffyg = field.Eyfsum
          ffzg = field.Ezfsum
        if which=='b':
          ffxg = field.Bxfsum
          ffyg = field.Byfsum
          ffzg = field.Bzfsum
        wtmp = zeros(nin,'d')
        em2d_getf2d_linear_serial(nin,xin,yin,fxin,fyin,fzin,
                                  field.xmin,field.ymin,
                                  field.dx,field.dy,
                                  field.nx,field.ny,
                                  ffxg,ffyg,ffzg)
        put(fx,iin,fxin)
        put(fy,iin,fyin)
        put(fz,iin,fzin)
      if nout>0:
        xout = take(x,iout)
        yout = take(y,iout)
        fxout = take(fx,iout)
        fyout = take(fy,iout)
        fzout = take(fz,iout)
        field = self.field
        wtmp = zeros(nout,'d')
        em2d_getf2d_linear_serial(nout,xout,yout,fxout,fyout,fzout,
                                  field.xmin,field.ymin,
                                  field.dx,field.dy,
                                  field.nx,field.ny,
                                  fxg,fyg,fzg)
        put(fx,iout,fxout)
        put(fy,iout,fyout)
        put(fz,iout,fzout)

  def fetchphifrompositions(self,x,z,phi):
    pass

  def loadrho(self,lzero=0):
    pass

  def loadj(self,ins_i=-1,nps_i=-1,is_i=-1,lzero=true):
    if self.l_onegrid:
      fields = [self.field]
    else:
      fields = [self.field,self.fpatchcoarse,self.fpatchfine]    

    # --- reallocate Jarray if needed
    if self.field.ntimes<>top.nsndts:
      self.field.ntimes=top.nsndts
      self.field.gchange()
      force_deposition=true
    else:
      force_deposition=false

    # --- zero proper portion of Jarray
    if lzero: 
      for i in range(top.nsndts-1,-1,-1):
        if force_deposition or (top.it-1)%(2**i)==0:
          for field in fields:
            field.Jarray[:,:,:,i] = 0.
        
    # --- loop over species
    for js,i,n,q,w in zip(arange(top.pgroup.ns),top.pgroup.ins-1,top.pgroup.nps,
                       top.pgroup.sq,top.pgroup.sw):
      if n == 0 or ((top.it-1)%top.pgroup.ndts[js]<>0 and not force_deposition): continue
      x,y,ux,uy,uz,xold,yold = self.transformparticles(
            top.pgroup.xp[i:i+n],top.pgroup.yp[i:i+n],top.pgroup.zp[i:i+n],
            top.pgroup.uxp[i:i+n],top.pgroup.uyp[i:i+n],top.pgroup.uzp[i:i+n],
            top.pgroup.pid[i:i+n,top.xoldpid-1],
            top.pgroup.pid[i:i+n,top.yoldpid-1],
            top.pgroup.pid[i:i+n,top.zoldpid-1]
            )
      if top.wpid==0:
        wfact = None
      else:
        wfact = top.pgroup.pid[i:i+n,top.wpid-1]
      # --- point J array to proper Jarray slice
      for field in fields:
        field.J = field.Jarray[:,:,:,top.ndtstorho[top.pgroup.ndts[js]-1]]
      # --- call routine performing current deposition
      self.setj(x,y,xold,yold,uz,top.pgroup.gaminv[i:i+n],q,w,wfact,top.dt*top.pgroup.ndts[js])

    # --- add slices
    if top.nsndts>1:
      for i in range(top.nsndts-2,-1,-1):
        if force_deposition or (top.it-1)%(2**i)==0:
          add_current_slice(self.field,i+1)
    
    # --- smooth current density 
    if self.l_smoothdensity:self.smoothdensity()

    # --- exchange slices of Jarray among processors
    self.apply_current_bc()

    # --- point J to first slice of Jarray
    for field in fields:
      field.J = field.Jarray[:,:,:,0]

    # --- get 1-D density
    if self.density_1d:
      for i in range(3):
       J = sum(self.field.J[:,:,i],1)
       for ii in range(shape(self.field.J[:,:,i])[1]):
         self.field.J[:,ii,i] = J


  def apply_current_bc(self):
    # --- apply periodic BC
    if self.field.xlbound==periodic:
      if npes<=1:
        self.field.Jarray[0:4,:,0:3,0]+=self.field.Jarray[-4:,:,0:3,0]
        self.field.Jarray[-4:,:,0:3,0] =self.field.Jarray[0:4,:,0:3,0]
    if self.field.ylbound==periodic:
      if npes<=1:
        self.field.Jarray[:,0:3,0:3,0]+=self.field.Jarray[:,-3:,0:3,0]
        self.field.Jarray[:,-3:,0:3,0] =self.field.Jarray[:,0:3,0:3,0]
    # --- exchange slices of Jarray among processors along z
    if npes>1:
      if me>0:
        mpi.send(self.field.Jarray[0:4,:,0:3,0],me-1)
      if me<npes-1:
        recv,status=mpi.recv(me+1)
        mpi.send(self.field.Jarray[-4:-2,:,0:3,0],me+1)
        self.field.Jarray[-4:,:,0:3,0]+=recv
      if me>0:
        recv,status=mpi.recv(me-1)
        self.field.Jarray[0:2,:,0:3,0]+=recv
    return
    if npes>1:
      if me>0:
        mpi.send(self.field.Jarray[0:2,:,2,0],me-1)
      if me<npes-1:
        recv,status=mpi.recv(me+1)
        mpi.send(self.field.Jarray[-4:-2,:,2,0],me+1)
        self.field.Jarray[-5:-3,:,2,0]+=recv
      if me>0:
        recv,status=mpi.recv(me-1)
        self.field.Jarray[0:2,:,2,0]+=recv

  def apply_current_bc_old(self):
    # --- apply periodic BC
    if self.field.xlbound==periodic:
      if npes<=1:
        self.field.Jarray[0:4,:,0:3,:]+=self.field.Jarray[-4:,:,0:3,:]
        self.field.Jarray[-4:,:,0:3,:] =self.field.Jarray[0:4,:,0:3,:]
    if self.field.ylbound==periodic:
      if npes<=1:
        self.field.Jarray[:,0:3,0:3,:]+=self.field.Jarray[:,-3:,0:3,:]
        self.field.Jarray[:,-3:,0:3,:] =self.field.Jarray[:,0:3,0:3,:]
    # --- exchange slices of Jarray among processors along z
    if npes>1:
      if me>0:
        mpi.send(self.field.Jarray[0:4,:,0:3,:],me-1)
      if me<npes-1:
        recv,status=mpi.recv(me+1)
        mpi.send(self.field.Jarray[-4:-2,:,0:3,:],me+1)
        self.field.Jarray[-4:,:,0:3,:]+=recv
      if me>0:
        recv,status=mpi.recv(me-1)
        self.field.Jarray[0:2,:,0:3,:]+=recv
    return
    if npes>1:
      if me>0:
        mpi.send(self.field.Jarray[0:2,:,2,:],me-1)
      if me<npes-1:
        recv,status=mpi.recv(me+1)
        mpi.send(self.field.Jarray[-4:-2,:,2,:],me+1)
        self.field.Jarray[-5:-3,:,2,:]+=recv
      if me>0:
        recv,status=mpi.recv(me-1)
        self.field.Jarray[0:2,:,2,:]+=recv

  def smoothdensity(self):
    mysmooth = smooth2d_lindman
#    mysmooth = smooth2d_121
    if self.l_onegrid:
      fields=[self.field]
    else:
      fields=[self.field,self.fpatchcoarse,self.fpatchfine]
    for field in fields:
      mysmooth(field.Jarray[:,:,0,0],field.nx,field.ny)
      mysmooth(field.Jarray[:,:,1,0],field.nx,field.ny)
      mysmooth(field.Jarray[:,:,2,0],field.nx,field.ny)
      
    
  def fetche(self):
#    if w3d.api_xlf2:
    w3d.xfsapi=top.pgroup.xp[w3d.ipminfsapi-1:w3d.ipminfsapi-1+w3d.npfsapi]
    w3d.yfsapi=top.pgroup.yp[w3d.ipminfsapi-1:w3d.ipminfsapi-1+w3d.npfsapi]
    w3d.zfsapi=top.pgroup.zp[w3d.ipminfsapi-1:w3d.ipminfsapi-1+w3d.npfsapi]
#    ex,ey,ez = self.transformfields(w3d.exfsapi,w3d.eyfsapi,w3d.ezfsapi)
    ex,ey,ez = self.transformfields(top.pgroup.ex[w3d.ipminfsapi-1:w3d.ipminfsapi-1+w3d.npfsapi],
                                    top.pgroup.ey[w3d.ipminfsapi-1:w3d.ipminfsapi-1+w3d.npfsapi],
                                    top.pgroup.ez[w3d.ipminfsapi-1:w3d.ipminfsapi-1+w3d.npfsapi])
    ex[:] = 0.
    ey[:] = 0.
    ez[:] = 0.
    x,y,ux,uy,uz,xold,yold = self.transformparticles(w3d.xfsapi,w3d.yfsapi,w3d.zfsapi)
    self.fetchefrompositions(x,y,ex,ey,ez)

  def fetchb(self):
#    if w3d.api_xlf2:
    w3d.xfsapi=top.pgroup.xp[w3d.ipminfsapi-1:w3d.ipminfsapi-1+w3d.npfsapi]
    w3d.yfsapi=top.pgroup.yp[w3d.ipminfsapi-1:w3d.ipminfsapi-1+w3d.npfsapi]
    w3d.zfsapi=top.pgroup.zp[w3d.ipminfsapi-1:w3d.ipminfsapi-1+w3d.npfsapi]
#    bx,by,bz = self.transformfields(w3d.bxfsapi,w3d.byfsapi,w3d.bzfsapi)
    bx,by,bz = self.transformfields(top.pgroup.bx[w3d.ipminfsapi-1:w3d.ipminfsapi-1+w3d.npfsapi],
                                    top.pgroup.by[w3d.ipminfsapi-1:w3d.ipminfsapi-1+w3d.npfsapi],
                                    top.pgroup.bz[w3d.ipminfsapi-1:w3d.ipminfsapi-1+w3d.npfsapi])
    bx[:] = 0.
    by[:] = 0.
    bz[:] = 0.
    x,y,ux,uy,uz,xold,yold = self.transformparticles(w3d.xfsapi,w3d.yfsapi,w3d.zfsapi)
    self.fetchbfrompositions(x,y,bx,by,bz)

  def fetchphi(self):
    pass

  def fetcha(self):
    pass

  def installconductor(self,conductor,
                            xmin=None,xmax=None,
                            ymin=None,ymax=None,
                            zmin=None,zmax=None,
                            dfill=top.largepos):
    pass

  def clearconductors(self):
    pass

  def optimizeconvergence(self,resetpasses=1):
    pass

  def move_window_fields(self):
    if (w3d.solvergeom not in [w3d.XZgeom]) or \
       (abs(top.zgrid-self.zgridprv)<0.5*self.dz):return 
#    if not self.l_moving_window or ((top.it%self.ndelta_t)!=0): return
#    if top.time < self.tmin_moving_main_window: return
    move_window_field(self.field)
    self.zgridprv=top.zgrid
    if not self.l_onegrid:
      self.fpatchfine.xminscatter+=w3d.dz
      self.fpatchfine.xmaxscatter+=w3d.dz
      self.fpatchfine.xmingather+=w3d.dz
      self.fpatchfine.xmaxgather+=w3d.dz

  def add_laser(self,field):
    if self.laser_profile is None: return

    if self.laser_amplitude_func is not None:
      self.laser_amplitude = self.laser_amplitude_func(top.time)
    elif self.laser_amplitude_table is not None:
      if top.time < self.laser_amplitude_table[0,1]:
        self.laser_amplitude = self.laser_amplitude_table[0,0]
      elif top.time >= self.laser_amplitude_table[-1,1]:
        self.laser_amplitude = self.laser_amplitude_table[-1,0]
      else:
        i = self.laser_amplitude_table_i
        while top.time > self.laser_amplitude_table[i+1,1]:
          i = i + 1
        self.laser_amplitude_table_i = i
        ww = ((top.time - self.laser_amplitude_table[i,1])/
           (self.laser_amplitude_table[i+1,1]-self.laser_amplitude_table[i,1]))
        self.laser_amplitude = ((1.-ww)*self.laser_amplitude_table[i,0] +
                                    ww *self.laser_amplitude_table[i+1,0])

    if self.laser_profile_func is not None:
      self.laser_profile = self.laser_profile_func(top.time)
      assert len(self.laser_profile) == field.ny+4,"The specified profile must be of length ny+4"

#    if (self.l_elaser_out_plane):
#      xx = (arange(self.nx+4) - 0.5)*self.field.dx+self.field.xmin
#    else:
    xx = (arange(field.ny+3) - 0.5)*field.dy + field.ymin

    if self.laser_frequency is not None:
      if self.laser_focus is not None:

        Z_R = (self.laser_gauss_width**2)/(clight/self.laser_frequency) #>> angular freq.!!
        z0 = self.laser_focus + self.laser_focus_velocity*top.time
        z0 = -z0/Z_R  ## now measured in Rayleigh ranges; negative value means we're upstream
 
        phi0_z=-top.time*self.laser_frequency - 0.5*arctan(z0)
        omgi_z_sqr=(1. + z0*z0) * (self.laser_gauss_width**2)
        phifac=0.5*z0/omgi_z_sqr   #>> factor of 0.5 due to slab model
        phase=phi0_z + phifac*(xx**2)

 ##       print 'Z_R %g ; z0  %g  omgi_z_sqr %g '%(Z_R, z0, omgi_z_sqr)
        laser_amp_factor2 = 1./sqrt(1. + z0*z0)   
        self.laser_amplitude=self.laser_amplitude*sqrt(laser_amp_factor2)

        qqww = exp(-0.5*(xx**2)/omgi_z_sqr)
        field.laser_profile[:-1]=qqww[:]

      else:
        phase = (xx*sin(self.laser_angle)/clight-top.time)*self.laser_frequency
    else:
      phase = 0.

    if (self.l_elaser_out_plane):
      field.Ez_in = self.laser_amplitude*field.laser_profile[:-1]*cos(phase)
    else:
      field.Bz_in = self.laser_amplitude*field.laser_profile[:-1]*cos(phase)

  def solve(self,iwhich=0):
    if any(top.fselfb<>0.):raise('Error:EM solver does not work if fselfb<>0.')
    if top.dt<>self.dtinit:raise('Time step has been changed since initialization of EM2D.')
    # --- Set nxl and nyl if using large stencil
    if(not self.l_onegrid):
      project_j(self.field,self.fpatchcoarse,self.fpatchfine)
    if self.l_onegrid:
      fields = [self.field]
    else:
      fields = [self.field,self.fpatchcoarse,self.fpatchfine]    
    for field in fields:
      if field.l_uselargestencil and field.nxl<>field.nx:
        field.nxl=field.nx
        field.nyl=field.ny
        field.gchange()
      self.add_laser(field)
      grimax(field)
      dt = top.dt/self.nfield_subcycle
      for i in range(self.nfield_subcycle):
        push_em_b(field,0.5*dt)
        push_em_e(field,dt)
        push_em_b(field,0.5*dt)
    self.move_window_fields()
    for field in fields:
      griuni(field)
    if not self.l_onegrid:set_substitute_fields(self.field,self.fpatchcoarse,self.fpatchfine)

  ##########################################################################
  # Define the basic plot commands
  def genericfp(self,data,title,l_transpose=true,**kw):
    f=self.field
    if self.solvergeom == w3d.XYgeom:
      settitles(title,'X','Y','t = %gs'%(top.time))
    elif self.solvergeom == w3d.XZgeom:
      settitles(title,'Z','X','t = %gs'%(top.time))
    if l_transpose:
     if npes>1:
      kw.setdefault('xmin',w3d.zmmin)
      kw.setdefault('xmax',w3d.zmmax)
     else:
      kw.setdefault('xmin',self.field.xmin)
      kw.setdefault('xmax',self.field.xmax)
     kw.setdefault('ymin',self.field.ymin)
     kw.setdefault('ymax',self.field.ymax)
    else:
     kw.setdefault('xmin',self.field.ymin)
     kw.setdefault('xmax',self.field.ymax)
     if npes>0:
      kw.setdefault('ymin',w3d.zmmin)
      kw.setdefault('ymax',w3d.zmmax)
     else:
      kw.setdefault('ymin',self.field.xmin)
      kw.setdefault('ymax',self.field.xmax)
    ppgeneric(grid=data,**kw)
      
  def fpex(self,**kw):
    if self.solvergeom == w3d.XZgeom:
      self.genericfp(gatherarray(self.field.Ey[1:w3d.nzp+1,...]),'E_x',**kw)
    elif self.solvergeom == w3d.XYgeom:
      self.genericfp(gatherarray(self.field.Ex),'E_x',**kw)

  def fpey(self,**kw):
    if self.solvergeom == w3d.XZgeom:
      self.genericfp(gatherarray(self.field.Ez[1:w3d.nzp+1,...]),'E_y',**kw)
    elif self.solvergeom == w3d.XYgeom:
      self.genericfp(gatherarray(self.field.Ey),'E_y',**kw)

  def fpez(self,**kw):
    if self.solvergeom == w3d.XZgeom:
      self.genericfp(gatherarray(self.field.Ex[1:w3d.nzp+1,...]),'E_z',**kw)
    elif self.solvergeom == w3d.XYgeom:
      self.genericfp(gatherarray(self.field.Ez),'E_z',**kw)

  def fpbx(self,**kw):
    if self.solvergeom == w3d.XZgeom:
      self.genericfp(gatherarray(self.field.By[1:w3d.nzp+1,...]),'B_x',**kw)
    elif self.solvergeom == w3d.XYgeom:
      self.genericfp(gatherarray(self.field.Bx),'B_x',**kw)

  def fpby(self,**kw):
    if self.solvergeom == w3d.XZgeom:
      self.genericfp(gatherarray(self.field.Bz[1:w3d.nzp+1,...]),'B_y',**kw)
    elif self.solvergeom == w3d.XYgeom:
      self.genericfp(gatherarray(self.field.By),'B_y',**kw)

  def fpbz(self,**kw):
    if self.solvergeom == w3d.XZgeom:
      self.genericfp(gatherarray(self.field.Bx[1:w3d.nzp+1,...]),'B_z',**kw)
    elif self.solvergeom == w3d.XYgeom:
      self.genericfp(gatherarray(self.field.Bz),'B_z',**kw)

  def fpjx(self,**kw):
    if self.solvergeom == w3d.XZgeom:
      self.genericfp(gatherarray(self.field.J[1:w3d.nzp+1,:,1]),'J_x',**kw)
    elif self.solvergeom == w3d.XYgeom:
      self.genericfp(gatherarray(self.field.J[:,:,0]),'J_x',**kw)

  def fpjy(self,**kw):
    if self.solvergeom == w3d.XZgeom:
      self.genericfp(gatherarray(self.field.J[1:w3d.nzp+1,:,2]),'J_y',**kw)
    elif self.solvergeom == w3d.XYgeom:
      self.genericfp(gatherarray(self.field.J[:,:,1]),'J_y',**kw)

  def fpjz(self,**kw):
    if self.solvergeom == w3d.XZgeom:
      self.genericfp(gatherarray(self.field.J[1:w3d.nzp+1,:,0]),'J_z',**kw)
    elif self.solvergeom == w3d.XYgeom:
      self.genericfp(gatherarray(self.field.J[:,:,2]),'J_z',**kw)

  def sezax(self):
    pass

  def sphiax(self):
    pass

  def rhodia(self):
    pass

  def gtlchg(self):
    pass

  def srhoax(self):
    pass

  def getese(self):
    pass

  def plbnd(self,h,fin=None):
    f = self.field
    if em2d.l_pml_cummer:
      b = f.bndexeybz_cummer
    else:
      b = f.bndexeybz
    g = fzeros([b.nx+1,b.ny+1],'d')

    for k in range(1, b.nbndy+1):
      jk1 = b.ntop1 + k * b.n1x
      for j in range(1, b.nx+1):
        jk = jk1 + j
        g[j-1,k-1] = h[jk-1]

    for k in range(  b.ny-b.nbndy+1, b.ny+1):
      jk1 = b.ntop2 + k * b.n1x
      for j in range(1, b.nx+1):
        jk = jk1 + j
        g[j-1,k-1] = h[jk-1]

    for k in range(b.nbndy+2, b.ny-b.nbndy-2+1):
      for j in range(1, b.nbndx+1):
        jk= b.nbot1 + j + k * b.nint
        g[j-1,k-1] = h[jk-1]

    k=b.nbndy+1
    for j in range(1, b.nbndx+1):
      jk=bndijk(f,j,k)
      g[j-1,k-1] = h[jk-1]

    for k in range(b.ny-b.nbndy-1,b.ny-b.nbndy+1):
      for j in range(1, b.nbndx+1):
        jk=bndijk(f,j,k)
        g[j-1,k-1] = h[jk-1]

    for k in range(b.nbndy+2, b.ny-b.nbndy-2+1):
      for j in range(b.nx-b.nbndx+1, b.nx+1):
        jk = b.nbot2+ j + k * b.nint
        g[j-1,k-1] = h[jk-1]

    k=b.nbndy+1
    for j in range(b.nx-b.nbndx+1, b.nx+1):
      jk=bndijk(f,j,k)
      g[j-1,k-1] = h[jk-1]

    for k in range(b.ny-b.nbndy-1,b.ny-b.nbndy+1):
      for j in range( b.nx-b.nbndx+1, b.nx+1):
        jk=bndijk(f,j,k)
        g[j-1,k-1] = h[jk-1]

    if fin is not None:
      grimax(f)
      g[b.nbndx-1:b.nbndx+f.nx+1,b.nbndy-1:b.nbndy+f.ny+1]=fin[:f.nx+2,:f.ny+2]
      griuni(f)
      
    ppgeneric(g)
    
  def fpezall(self,**kw):
    f = self.field
    if em2d.l_pml_cummer:
      h = f.bndbxbyez_cummer.Bz
    else:
      h = f.bndbxbyez.Bzx+f.bndbxbyez.Bzy
    self.plbnd(-h*clight,f.Ez)
    
  def fpbxall(self,**kw):
    f = self.field
    if em2d.l_pml_cummer:
      h = f.bndbxbyez_cummer.Ex
    else:
      h = f.bndbxbyez.Ex
    self.plbnd(h,f.Bx)
        
  def fpbyall(self,**kw):
    f = self.field
    if em2d.l_pml_cummer:
      h = f.bndbxbyez_cummer.Ey
    else:
      h = f.bndbxbyez.Ey
    self.plbnd(h,f.By)
        
  def step(self,n=1):
    for i in range(n):
       self.onestep()
       
  def onestep(self):
    # --- call beforestep functions
    callbeforestepfuncs.callfuncsinlist()
    
    top.zgrid+=top.vbeamfrm*top.dt
    top.zbeam=top.zgrid

    for js in range(top.pgroup.ns):
      self.fetcheb(js)
    for js in range(top.pgroup.ns):
      self.push_velocity_second_half(js)
      self.set_gamma(js)
      self.push_positions(js)

    particleboundaries3d(top.pgroup)

    # --- call beforeloadrho functions
    beforelr.callfuncsinlist()
    self.loadrho()
    self.loadj()
      
 #   self.solve2ndhalf()
    self.solve()
    
    for js in range(top.pgroup.ns):
      self.fetcheb(js)
    
    for js in range(top.pgroup.ns):
      self.push_velocity_first_half(js)
      self.set_gamma(js)

    # --- update time, time counter
    top.time+=top.dt
    if top.it%top.nhist==0:
       zmmnt()
       minidiag(top.it,top.time,top.lspecial)
    top.it+=1

    # --- call afterstep functions
    callafterstepfuncs.callfuncsinlist()

  def fetcheb(self,js):
    if self.l_verbose:print me,'enter fetcheb'
    pg = top.pgroup
    np = pg.nps[js]
    if np==0:return
    w3d.ipminfsapi=pg.ins[js]
    w3d.npfsapi=pg.nps[js]
    self.fetche()
    self.fetchb()

  def push_velocity_first_half(self,js):
    if self.l_verbose:print me,'enter push_ions_velocity_first_half'
    pg = top.pgroup
    np = pg.nps[js]
    if np==0:return
    il = pg.ins[js]-1
    iu = il+pg.nps[js]
    if pg.lebcancel_pusher:
      ebcancelpush3d(np,pg.uxp[il:iu],pg.uyp[il:iu],pg.uzp[il:iu],pg.gaminv[il:iu],
                        pg.ex[il:iu], pg.ey[il:iu], pg.ez[il:iu],
                        pg.bx[il:iu], pg.by[il:iu], pg.bz[il:iu],
                        pg.sq[js],pg.sm[js],top.dt,1)
    else:
      # --- push velocity from electric field (half step)
      epush3d(np,pg.uxp[il:iu],pg.uyp[il:iu],pg.uzp[il:iu],
                 pg.ex[il:iu], pg.ey[il:iu], pg.ez[il:iu], 
                 pg.sq[js],pg.sm[js],0.5*top.dt)
      # --- update gamma
      self.set_gamma(js)
      # --- push velocity from magnetic field
      bpush3d (np,pg.uxp[il:iu],pg.uyp[il:iu],pg.uzp[il:iu],pg.gaminv[il:iu],
                  pg.bx[il:iu], pg.by[il:iu], pg.bz[il:iu], 
                  pg.sq[js],pg.sm[js],0.5*top.dt, top.ibpush)

    if self.l_verbose:print me,'exit push_ions_velocity_first_half'
    
  def push_velocity_second_half(self,js):
    if self.l_verbose:print me,'enter push_ions_velocity_second_half'
    pg = top.pgroup
    np = pg.nps[js]
    if np==0:return
    il = pg.ins[js]-1
    iu = il+pg.nps[js]
    if pg.lebcancel_pusher:
      ebcancelpush3d(np,pg.uxp[il:iu],pg.uyp[il:iu],pg.uzp[il:iu],pg.gaminv[il:iu],
                        pg.ex[il:iu], pg.ey[il:iu], pg.ez[il:iu],
                        pg.bx[il:iu], pg.by[il:iu], pg.bz[il:iu],
                        pg.sq[js],pg.sm[js],top.dt,2)
    else:
      # --- push velocity from magnetic field
      bpush3d (np,pg.uxp[il:iu],pg.uyp[il:iu],pg.uzp[il:iu],pg.gaminv[il:iu],
                  pg.bx[il:iu], pg.by[il:iu], pg.bz[il:iu], 
                  pg.sq[js],pg.sm[js],0.5*top.dt, top.ibpush)
      # --- push velocity from electric field (half step)
      epush3d(np,pg.uxp[il:iu],pg.uyp[il:iu],pg.uzp[il:iu],
                 pg.ex[il:iu], pg.ey[il:iu], pg.ez[il:iu], 
                 pg.sq[js],pg.sm[js],0.5*top.dt)
      # --- update gamma
      self.set_gamma(js)

    if self.l_verbose:print me,'exit push_ions_velocity_second_half'
    
  def set_gamma(self,js):
    if self.l_verbose:print me,'enter set_gamma'
    pg = top.pgroup
    np = pg.nps[js]
    if np==0:return
    il = pg.ins[js]-1
    iu = il+pg.nps[js]
    # --- update gamma
    gammaadv(np,pg.gaminv[il:iu],pg.uxp[il:iu],pg.uyp[il:iu],pg.uzp[il:iu],
             top.gamadv,top.lrelativ)

    if self.l_verbose:print me,'exit push_ions_velocity_second_half'
    
  def push_positions(self,js):
    if self.l_verbose:print me,'enter push_ions_positions'
    pg = top.pgroup
    np = pg.nps[js]
    if np==0:return
    il = pg.ins[js]-1
    iu = il+pg.nps[js]
    xpush3d(np,pg.xp[il:iu],pg.yp[il:iu],pg.zp[il:iu],
                   pg.uxp[il:iu],pg.uyp[il:iu],pg.uzp[il:iu],
                   pg.gaminv[il:iu],top.dt)      

    if self.l_verbose:print me,'exit push_ions_positions'

  def apply_bndconditions(self,js):
    if self.l_verbose:print me,'enter apply_ions_bndconditions'
    # --- apply boundary conditions
    pg = top.pgroup
    if pg.nps[js]==0:return
    self.apply_bnd_conditions(js)
    if self.l_verbose:print me,'exit apply_ions_bndconditions'
    
  def apply_bnd_conditions(self,js):
    if self.l_verbose:print me,'enter apply_bnd_conditions'
    pg = top.pgroup
    if pg.nps[js]==0:return
    il = pg.ins[js]-1
    iu = il+pg.nps[js]
    stckxy3d(pg.nps[js],pg.xp[il:iu],w3d.xmmax,w3d.xmmin,w3d.dx,
                  pg.yp[il:iu],w3d.ymmax,w3d.ymmin,w3d.dy,
                  pg.zp[il:iu],w3d.zmminlocal,w3d.dz,
                  pg.uxp[il:iu],pg.uyp[il:iu],pg.uzp[il:iu],pg.gaminv[il:iu],
                  top.zgrid,top.zbeam,w3d.l2symtry,w3d.l4symtry,top.pboundxy,true)
    if js==0 or js==w3d.nzp-1:
      if js==0:top.pboundnz=-1
      if js==w3d.nzp-1:top.pbound0=-1
      zpartbndwithdata(pg.nps[js],pg.zp[il:iu],pg.uzp[il:iu],pg.gaminv[il:iu],
                       w3d.zmmaxlocal,w3d.zmminlocal,w3d.dz,top.zgrid)
      if js==0:top.pboundnz=0
      if js==w3d.nzp-1:top.pbound0=0
    if self.scraper is not None:self.scraper.scrape(js)
    processlostpart(pg,js+1,top.clearlostpart,top.time+top.dt*pg.ndts[js],top.zbeam)
    if self.l_verbose:print me,'enter apply_bnd_conditions'

  def initfrompoisson(self):
    tmpbound0  = w3d.bound0
    tmpboundnz = w3d.boundnz
    tmpboundxy = w3d.boundxy
    w3d.bound0  = dirichlet
    w3d.boundnz = dirichlet
    w3d.boundxy = dirichlet
    init_base(w3d.nx,w3d.nzlocal,w3d.dx,w3d.dz,w3d.xmmin,w3d.zmminlocal,lparallel)
    fstypecp = top.fstype
    top.fstype = 10
    loadrho()
    vp3d(0)
    top.fstype = fstypecp
    w3d.bound0 = tmpbound0
    w3d.boundnz = tmpboundnz
    w3d.boundxy = tmpboundxy
    grimax(self.field)
#    self.field.Ex[:w3d.nzlocal,:w3d.nx+1] = -transpose(frz.basegrid.phi[1:-1:,2:-1]-frz.basegrid.phi[1:-1,1:-2])/w3d.dz
#    self.field.Ey[:w3d.nzlocal+1,:w3d.nx] = -transpose(frz.basegrid.phi[2:-1,1:-1]-frz.basegrid.phi[1:-2,1:-1])/w3d.dx
    self.field.Ex[:w3d.nzlocal+2,:w3d.nx+3] = -transpose(frz.basegrid.phi[::,1:]-frz.basegrid.phi[:,:-1])/w3d.dz
    self.field.Ey[:w3d.nzlocal+3,:w3d.nx+2] = -transpose(frz.basegrid.phi[1:,:]-frz.basegrid.phi[:-1,:])/w3d.dx
    griuni(self.field)
  
##############################################################################
# --- This can only be done after the class is defined.
#try:
#  psyco.bind(EM1D)
#except NameError:
#  pass

