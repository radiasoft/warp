"""Class for doing 2 1/2 D electromagnetic solver using code adapted from the
emi code"""
from warp import *
import operator

try:
  import psyco
except ImportError:
  pass

##############################################################################
class EM2D(object):
  
  __w3dinputs__ = ['solvergeom','nx','ny','nz',
                   'xmmin','xmmax','ymmin','ymmax','zmmin','zmmax',
                   'bound0','boundnz','boundxy','l2symtry','l4symtry',
                   'solvergeom']
  __em2dinputs__ = ['l_onegrid','l_copyfields','l_moving_window',
                    'tmin_moving_main_window','l_smoothdensity',
                    'l_elaser_out_plane','ndelta_t',
                   'ntamp_scatter','ntamp_gather']
  __flaginputs__ = {'l_apply_pml':true,'nbndx':10,'nbndy':10,
                    'l_particles_weight':false,'l_usecoeffs':false,
                    'laser_amplitude':1.,'laser_profile':None,
                    'laser_gauss_width':None,'laser_angle':0.,
                    'laser_wavelength':None,'laser_wavenumber':None,
                    'laser_frequency':None,'laser_source':2,
                    'density_1d':False,'nfield_subcycle':1,
                    'autoset_timestep':true,'dtcoef':0.99}

  def __init__(self,**kw):
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
      self.dz = (self.zmmax - self.zmmin)/self.nz
      self.zmesh = self.zmmin + arange(0,self.nz+1)*self.dz
      self.zmminlocal = self.zmmin
      self.zmmaxlocal = self.zmmax
      self.zmeshlocal = self.zmminlocal + arange(0,self.nz+1)*self.dz
      self.nx = self.nz
      self.dx = self.dz
      self.xmesh = self.zmeshlocal
      self.xmmin = self.zmmin
      self.xmmax = self.zmmax
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
        self.nfield_subcycle=int(top.dt/dt)+1
    self.dtinit = top.dt
    # --- Create field and source arrays and other arrays.
    self.allocatefieldarrays()

    # --- Handle laser inputs
    self.setuplaser()

  def allocatefieldarrays(self):
    # --- Code transcribed from init_fields
    self.field = EM2D_FIELDtype()
    self.field.l_usecoeffs = self.l_usecoeffs
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
    if self.laser_profile == 'gaussian':
      assert self.laser_gauss_width is not None,\
             "For a gaussian laser, the width must be specified using laser_gauss_width"
#      xx = arange(self.nx+4)*f.dx+f.xmin*f.dx - 0.5*f.nx*f.dx
#      self.laser_profile = exp(-(xx/self.laser_gauss_width)**2/2.)
      yy = arange(f.ny+4)*f.dy+f.ymin*f.dy - 0.5*f.ny*f.dy
      f.laser_profile = exp(-(yy/self.laser_gauss_width)**2/2.)
    elif operator.isSequenceType(self.laser_profile):
      assert len(f.laser_profile) == f.ny+4,"The specified profile must be of length ny+4"
    elif callable(self.laser_profile):
      self.laser_profile_func = self.laser_profile

  def transformparticles(self,x,y,z,ux=None,uy=None,uz=None):
    if self.solvergeom == w3d.XYgeom:
      return x,y,ux,uy,uz
    elif self.solvergeom == w3d.XZgeom:
      if self.l_moving_window and not self.l_elaser_out_plane:
        return z,x,uz,ux,uy
      else:
#        return x,z,ux,uz,uy
        return z,x,uz,ux,uy

  def transformfields(self,fx,fy,fz):
    if self.solvergeom == w3d.XYgeom:
      return fx,fy,fz
    elif self.solvergeom == w3d.XZgeom:
      if self.l_moving_window and not self.l_elaser_out_plane:
        return fz,fx,fy
      else:
#        return fx,fz,fy
        return fz,fx,fy

  def setj(self,x,y,ux,uy,uz,gaminv,q,w,wfact,dt):
    n = len(x)
    if n == 0: return
    if wfact is None:
      wfact = zeros(n,'d')
      l_particles_weight = false
    else:
      l_particles_weight = true
    if self.l_onegrid:
      depose_current_em2d(n,x,y,ux,uy,uz,gaminv,wfact,q*w,dt,l_particles_weight,self.field,self.field)
    else:
      depose_current_em2d(n,x,y,ux,uy,uz,gaminv,wfact,q*w,dt,l_particles_weight,self.field,self.fpatchfine)
    
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
      print 'nin,nout',len(iin),len(iout)
      print min(x),max(x),min(y),max(y)
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
          self.field.Jarray[:,:,:,i] = 0.
          if not self.l_onegrid:
            self.fpatchcoarse.Jarray[:,:,:,i] = 0.
            self.fpatchfine.Jarray[:,:,:,i] = 0.
        
    # --- loop over species
    for js,i,n,q,w in zip(arange(top.pgroup.ns),top.pgroup.ins-1,top.pgroup.nps,
                       top.pgroup.sq,top.pgroup.sw):
      if n == 0 or ((top.it-1)%top.pgroup.ndts[js]<>0 and not force_deposition): continue
      x,y,ux,uy,uz = self.transformparticles(
            top.pgroup.xp[i:i+n],top.pgroup.yp[i:i+n],top.pgroup.zp[i:i+n],
            top.pgroup.uxp[i:i+n],top.pgroup.uyp[i:i+n],top.pgroup.uzp[i:i+n])
      if top.wpid==0:
        wfact = None
      else:
        wfact = top.pgroup.pid[i:i+n,top.wpid-1]
      # --- point J array to proper Jarray slice
      self.field.J = self.field.Jarray[:,:,:,top.ndtstorho[top.pgroup.ndts[js]-1]]
      # --- call routine performing current deposition
      self.setj(x,y,ux,uy,uz,top.pgroup.gaminv[i:i+n],q,w,wfact,top.dt*top.pgroup.ndts[js])

    # --- add slices
    if top.nsndts>1:
      for i in range(top.nsndts-2,-1,-1):
        if force_deposition or (top.it-1)%(2**i)==0:
          add_current_slice(self.field,i+1)

    # --- point J to first slice of Jarray
    self.field.J = self.field.Jarray[:,:,:,0]
    
    # --- smooth current density 
    if self.l_smoothdensity:self.smoothdensity()

    # --- get 1-D density
    if self.density_1d:
      for i in range(3):
       J = sum(self.field.J[:,:,i],1)
       for ii in range(shape(self.field.J[:,:,i])[1]):
         self.field.J[:,ii,i] = J

  def smoothdensity(self):
    smooth2d_lindman(self.field.J[:,:,0],self.field.nx,self.field.ny)
    smooth2d_lindman(self.field.J[:,:,1],self.field.nx,self.field.ny)
    smooth2d_lindman(self.field.J[:,:,2],self.field.nx,self.field.ny)

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
    x,y,ux,uy,uz = self.transformparticles(w3d.xfsapi,w3d.yfsapi,w3d.zfsapi)
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
    x,y,ux,uy,uz = self.transformparticles(w3d.xfsapi,w3d.yfsapi,w3d.zfsapi)
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
    xx = (arange(field.ny+3) - 0.5)*field.dy+field.ymin

    if self.laser_frequency is not None:
      phase = (xx*sin(self.laser_angle)/clight-top.time)*self.laser_frequency
    else:
      phase = 0.

    if (self.l_elaser_out_plane):
      field.Ez_in = self.laser_amplitude*field.laser_profile[:-1]*cos(phase)
    else:
      field.Bz_in = self.laser_amplitude*field.laser_profile[:-1]*cos(phase)

  def solve(self,iwhich=0):
    if top.dt<>self.dtinit:raise('Time step has been changed since initialization of EM2D.')
    if(not self.l_onegrid):
      project_j(self.field,self.fpatchcoarse,self.fpatchfine)
    if self.l_onegrid:
      fields = [self.field]
    else:
      fields = [self.field,self.fpatchcoarse,self.fpatchfine]    
    for field in fields:
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
      kw.setdefault('xmin',self.field.xmin)
      kw.setdefault('xmax',self.field.xmax)
      kw.setdefault('ymin',self.field.ymin)
      kw.setdefault('ymax',self.field.ymax)
    else:
      kw.setdefault('xmin',self.field.ymin)
      kw.setdefault('xmax',self.field.ymax)
      kw.setdefault('ymin',self.field.xmin)
      kw.setdefault('ymax',self.field.xmax)
    ppgeneric(grid=data,**kw)
      
  def fpex(self,**kw):
    if self.solvergeom == w3d.XZgeom:
      self.genericfp(self.field.Ey,'E_x',**kw)
    elif self.solvergeom == w3d.XYgeom:
      self.genericfp(self.field.Ex,'E_x',**kw)

  def fpey(self,**kw):
    if self.solvergeom == w3d.XZgeom:
      self.genericfp(self.field.Ez,'E_y',**kw)
    elif self.solvergeom == w3d.XYgeom:
      self.genericfp(self.field.Ey,'E_y',**kw)

  def fpez(self,**kw):
    if self.solvergeom == w3d.XZgeom:
      self.genericfp(self.field.Ex,'E_z',**kw)
    elif self.solvergeom == w3d.XYgeom:
      self.genericfp(self.field.Ez,'E_z',**kw)

  def fpbx(self,**kw):
    if self.solvergeom == w3d.XZgeom:
      self.genericfp(self.field.By,'B_x',**kw)
    elif self.solvergeom == w3d.XYgeom:
      self.genericfp(self.field.Bx,'B_x',**kw)

  def fpby(self,**kw):
    if self.solvergeom == w3d.XZgeom:
      self.genericfp(self.field.Bz,'B_y',**kw)
    elif self.solvergeom == w3d.XYgeom:
      self.genericfp(self.field.By,'B_y',**kw)

  def fpbz(self,**kw):
    if self.solvergeom == w3d.XZgeom:
      self.genericfp(self.field.Bx,'B_z',**kw)
    elif self.solvergeom == w3d.XYgeom:
      self.genericfp(self.field.Bz,'B_z',**kw)

  def fpjx(self,**kw):
    if self.solvergeom == w3d.XZgeom:
      self.genericfp(self.field.J[:,:,1],'J_x',**kw)
    elif self.solvergeom == w3d.XYgeom:
      self.genericfp(self.field.J[:,:,0],'J_x',**kw)

  def fpjy(self,**kw):
    if self.solvergeom == w3d.XZgeom:
      self.genericfp(self.field.J[:,:,2],'J_y',**kw)
    elif self.solvergeom == w3d.XYgeom:
      self.genericfp(self.field.J[:,:,1],'J_y',**kw)

  def fpjz(self,**kw):
    if self.solvergeom == w3d.XZgeom:
      self.genericfp(self.field.J[:,:,0],'J_z',**kw)
    elif self.solvergeom == w3d.XYgeom:
      self.genericfp(self.field.J[:,:,2],'J_z',**kw)

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
        
# --- This can only be done after the class is defined.
try:
  psyco.bind(EM2D)
except NameError:
  pass

