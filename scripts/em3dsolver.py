"""Class for doing 3 D electromagnetic solver """
from warp import *
import operator

try:
  import psyco
except ImportError:
  pass

class EM3D(object):
  
  __w3dinputs__ = ['solvergeom','nx','ny','nzlocal','nz',
                   'xmmin','xmmax','ymmin','ymmax','zmmin','zmmax',
                   'zmminlocal','zmmaxlocal',
                   'bound0','boundnz','boundxy','l2symtry','l4symtry',
                   'solvergeom']
  __em3dinputs__ = []
#  __em3dinputs__ = ['l_onegrid','l_copyfields','l_moving_window',
#                    'tmin_moving_main_window','l_smoothdensity',
#                    'l_elaser_out_plane','ndelta_t',
#                   'ntamp_scatter','ntamp_gather']
  __flaginputs__ = {'l_apply_pml':true,'nbndx':10,'nbndy':10,'nbndz':10,
                    'nxguard':1,'nyguard':1,'nzguard':1,
                    'l_particles_weight':false,'l_usecoeffs':false,
                    'l_pushf':false,'l_verbose':false,
                    'laser_amplitude':1.,'laser_profile':None,
                    'laser_gauss_width':None,'laser_angle':0.,
                    'laser_wavelength':None,'laser_wavenumber':None,
                    'laser_frequency':None,'laser_source':2,
                    'laser_focus':None,'laser_focus_velocity':0.,
                    'nfield_subcycle':1,
                    'autoset_timestep':true,'dtcoef':0.99}

  def __init__(self,**kw):
    top.allspecl = true
    top.lcallfetchb = true
    top.lgridqnt = true
    self.zgridprv=top.zgrid
    # --- Make sure the refinement is turned off

    # --- Save input parameters
    for name in EM3D.__w3dinputs__:
      if name not in self.__dict__:
        #self.__dict__[name] = kw.pop(name,getattr(w3d,name)) # Python2.3
        self.__dict__[name] = kw.get(name,getattr(w3d,name))
      if kw.has_key(name): del kw[name]
    for name in EM3D.__em3dinputs__:
      if name not in self.__dict__:
        #self.__dict__[name] = kw.pop(name,getattr(em3d,name)) # Python2.3
        self.__dict__[name] = kw.get(name,getattr(em3d,name))
      if kw.has_key(name): del kw[name]
    for name,defvalue in EM3D.__flaginputs__.iteritems():
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

    # --- Calculate mesh sizes
    self.dx = (self.xmmax - self.xmmin)/self.nx
    self.xmesh = self.xmmin + arange(0,self.nx+1)*self.dx
    self.dy = (self.ymmax - self.ymmin)/self.ny
    self.ymesh = self.ymmin + arange(0,self.ny+1)*self.dy
    self.dz = (self.zmmaxlocal - self.zmminlocal)/self.nzlocal
    self.zmesh = self.zmminlocal + arange(0,self.nzlocal+1)*self.dz
    self.zmeshlocal = self.zmminlocal + arange(0,self.nzlocal+1)*self.dz

    # --- set time step as a fraction of Courant condition
    # --- also set self.nfield_subcycle if top.dt over Courant condition times dtcoef
    if self.autoset_timestep:
      dt=self.dtcoef/(clight*sqrt(1./self.dx**2+1./self.dy**2+1./self.dz**2))
      if top.dt==0.:
        top.dt=dt
      else:
        self.nfield_subcycle=int(top.dt/dt)+1
    self.dtinit = top.dt
    
    # --- Create field and source arrays and other arrays.
    self.allocatefieldarrays()

    # --- Handle laser inputs
#    self.setuplaser()

    installbeforeloadrho(self.solve2ndhalf)
    registersolver(self)

  def allocatefieldarrays(self):
    self.block = pyinit_3dem_block(self.nx, 
                                   self.ny, 
                                   self.nzlocal, 
                                   self.nbndx, 
                                   self.nbndy, 
                                   self.nbndz, 
                                   self.nxguard, 
                                   self.nyguard, 
                                   self.nzguard, 
                                   top.dt, 
                                   self.dx, 
                                   self.dy, 
                                   self.dz, 
                                   clight, 
                                   mu0, 
                                   self.xmmin, 
                                   self.ymmin, 
                                   self.zmminlocal, 
                                   self.bounds[0], 
                                   self.bounds[2], 
                                   self.bounds[4], 
                                   self.bounds[1], 
                                   self.bounds[3], 
                                   self.bounds[5])
    self.field = self.block.core.yf
    
  def addpatch(self,ixpatch,iypatch,nxpatch,nypatch,rap):
    # --- Initializes refinement patch
    em3d.l_onegrid = false
    self.l_onegrid = em3d.l_onegrid
    self.fpatches.append([EM3D_FIELDtype(),EM3D_FIELDtype()])
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
    self.fpatchfine.xminpatch_scatter = self.fpatchfine.xmin+em3d.ntamp_scatter*rap*self.fpatchfine.dx
    self.fpatchfine.xmaxpatch_scatter = self.fpatchfine.xmax-em3d.ntamp_scatter*rap*self.fpatchfine.dx
    self.fpatchfine.yminpatch_scatter = self.fpatchfine.ymin+em3d.ntamp_scatter*rap*self.fpatchfine.dy
    self.fpatchfine.ymaxpatch_scatter = self.fpatchfine.ymax-em3d.ntamp_scatter*rap*self.fpatchfine.dy
    self.fpatchfine.xminpatch_gather = self.fpatchfine.xmin+em3d.ntamp_gather*rap*self.fpatchfine.dx
    self.fpatchfine.xmaxpatch_gather = self.fpatchfine.xmax-em3d.ntamp_gather*rap*self.fpatchfine.dx
    self.fpatchfine.yminpatch_gather = self.fpatchfine.ymin+em3d.ntamp_gather*rap*self.fpatchfine.dy
    self.fpatchfine.ymaxpatch_gather = self.fpatchfine.ymax-em3d.ntamp_gather*rap*self.fpatchfine.dy
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
      f.laser_profile = exp(-(yy/self.laser_gauss_width)**2/2.)
    elif operator.isSequenceType(self.laser_profile):
      assert len(f.laser_profile) == f.ny+4,"The specified profile must be of length ny+4"
    elif callable(self.laser_profile):
      self.laser_profile_func = self.laser_profile
    
  def fetchefrompositions(self,x,y,z,ex,ey,ez):
    n = len(x)
    if n == 0: return
    f = self.block.core.yf
    getf3d_linear(n,x,y,z,ex,ey,ez,
                  f.xmin,f.ymin,f.zmin,
                  f.dx,f.dy,f.dz,
                  f.nx,f.ny,f.nz,
                  f.Ex,f.Ey,f.Ez)

  def fetchbfrompositions(self,x,y,z,bx,by,bz):
    n = len(x)
    if n == 0: return
    f = self.block.core.yf
    getf3d_linear(n,x,y,z,bx,by,bz,
                  f.xmin,f.ymin,f.zmin,
                  f.dx,f.dy,f.dz,
                  f.nx,f.ny,f.nz,
                  f.Bx,f.By,f.Bz)
    
  def fetchphifrompositions(self,x,z,phi):
    pass

  def loadrho(self,lzero=true):
    if not self.l_pushf:return
    if self.l_verbose:print 'loadrho',lzero
    # --- reallocate Jarray if needed
    if self.field.ntimes<>top.nsndts:
      self.field.ntimes=top.nsndts
      self.field.gchange()
      force_deposition=true
    else:
      force_deposition=false

    # --- zero proper portion of Rhoarray
    if lzero: 
      for i in range(top.nsndts-1,-1,-1):
        if force_deposition or (top.it-1)%(2**i)==0:
          self.field.Rhoarray[:,:,:,i] = 0.

    # --- loop over species
    for js,i,n,q,w in zip(arange(top.pgroup.ns),top.pgroup.ins-1,top.pgroup.nps,
                       top.pgroup.sq,top.pgroup.sw):
      if n == 0 or ((top.it-1)%top.pgroup.ndts[js]<>0 and not force_deposition): continue
      if top.wpid==0:
        wfact = zeros(n,'d')
        l_particles_weight = false
      else:
        wfact = top.pgroup.pid[i:i+n,top.wpid-1]
        l_particles_weight = true
      # --- point rho array to proper Jarray slice
      self.field.Rho = self.field.Rhoarray[:,:,:,top.ndtstorho[top.pgroup.ndts[js]-1]]
      # --- call routine performing current deposition
      depose_rho_linear_serial(self.field.Rho,n,
                               top.pgroup.xp[i:i+n],
                               top.pgroup.yp[i:i+n],
                               top.pgroup.zp[i:i+n],
                               wfact,q*w,
                               self.block.xmin,self.block.ymin,self.block.zmin,
                               self.block.core.yf.dx,self.block.core.yf.dy,self.block.core.yf.dz,
                               self.block.core.yf.nx,self.block.core.yf.ny,self.block.core.yf.nz,
                               l_particles_weight)

    # --- add slices
    if top.nsndts>1:
      for i in range(top.nsndts-2,-1,-1):
        if force_deposition or (top.it-1)%(2**i)==0:
          add_rho_slice_3d(self.field,i+1)

    # --- apply boundary condition on current
    self.apply_rho_bc()

    # --- point Rho to first slice of Rhoarray
    self.field.Rho = self.field.Rhoarray[:,:,:,0]

  def apply_rho_bc(self):
    # --- point Rho to first slice of Rhoarray
    self.field.Rho = self.field.Rhoarray[:,:,:,0]
    em3d_exchange_rho(self.block)

  def loadj(self,ins_i=-1,nps_i=-1,is_i=-1,lzero=true):
    if self.l_verbose:print 'loadj'
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
          self.field.Jarray[:,:,:,:,i] = 0.
        
    # --- loop over species
    for js,i,n,q,w in zip(arange(top.pgroup.ns),top.pgroup.ins-1,top.pgroup.nps,
                       top.pgroup.sq,top.pgroup.sw):
      if n == 0 or ((top.it-1)%top.pgroup.ndts[js]<>0 and not force_deposition): continue
      if top.wpid==0:
        wfact = zeros(n,'d')
        l_particles_weight = false
      else:
        wfact = top.pgroup.pid[i:i+n,top.wpid-1]
        l_particles_weight = true
      # --- point J array to proper Jarray slice
      self.field.J = self.field.Jarray[:,:,:,:,top.ndtstorho[top.pgroup.ndts[js]-1]]
      # --- call routine performing current deposition
      depose_jxjyjz_esirkepov_linear_serial(self.field.J,n,
                                            top.pgroup.xp[i:i+n],
                                            top.pgroup.yp[i:i+n],
                                            top.pgroup.zp[i:i+n],
                                            top.pgroup.uxp[i:i+n],
                                            top.pgroup.uyp[i:i+n],
                                            top.pgroup.uzp[i:i+n],
                                            top.pgroup.gaminv[i:i+n],wfact,q*w,
                                            self.block.xmin,self.block.ymin,self.block.zmin,
                                            top.dt*top.pgroup.ndts[js],
                                            self.block.core.yf.dx,self.block.core.yf.dy,self.block.core.yf.dz,
                                            self.block.core.yf.nx,self.block.core.yf.ny,self.block.core.yf.nz,
                                            l_particles_weight)

    # --- add slices
    if top.nsndts>1:
      for i in range(top.nsndts-2,-1,-1):
        if force_deposition or (top.it-1)%(2**i)==0:
          add_current_slice_3d(self.field,i+1)

    # --- apply boundary condition on current
    self.apply_current_bc()

    # --- point J to first slice of Jarray
    self.field.J = self.field.Jarray[:,:,:,:,0]
    
    # --- smooth current density 
#    if self.l_smoothdensity:self.smoothdensity()


  def apply_current_bc(self):
    # --- point J to first slice of Jarray
    self.field.J = self.field.Jarray[:,:,:,:,0]
    em3d_exchange_j(self.block)

  def smoothdensity(self):
    smooth2d_lindman(self.field.J[:,:,0],self.field.nx,self.field.ny)
    smooth2d_lindman(self.field.J[:,:,1],self.field.nx,self.field.ny)
    smooth2d_lindman(self.field.J[:,:,2],self.field.nx,self.field.ny)

  def fetche(self):
    x=top.pgroup.xp[w3d.ipminfsapi-1:w3d.ipminfsapi-1+w3d.npfsapi]
    y=top.pgroup.yp[w3d.ipminfsapi-1:w3d.ipminfsapi-1+w3d.npfsapi]
    z=top.pgroup.zp[w3d.ipminfsapi-1:w3d.ipminfsapi-1+w3d.npfsapi]
    ex = top.pgroup.ex[w3d.ipminfsapi-1:w3d.ipminfsapi-1+w3d.npfsapi]
    ey = top.pgroup.ey[w3d.ipminfsapi-1:w3d.ipminfsapi-1+w3d.npfsapi]
    ez = top.pgroup.ez[w3d.ipminfsapi-1:w3d.ipminfsapi-1+w3d.npfsapi]
    ex[:] = 0.
    ey[:] = 0.
    ez[:] = 0.
    self.fetchefrompositions(x,y,z,ex,ey,ez)
#    print 'e',ex,ey,ez

  def fetchb(self):
    x=top.pgroup.xp[w3d.ipminfsapi-1:w3d.ipminfsapi-1+w3d.npfsapi]
    y=top.pgroup.yp[w3d.ipminfsapi-1:w3d.ipminfsapi-1+w3d.npfsapi]
    z=top.pgroup.zp[w3d.ipminfsapi-1:w3d.ipminfsapi-1+w3d.npfsapi]
    bx = top.pgroup.bx[w3d.ipminfsapi-1:w3d.ipminfsapi-1+w3d.npfsapi]
    by = top.pgroup.by[w3d.ipminfsapi-1:w3d.ipminfsapi-1+w3d.npfsapi]
    bz = top.pgroup.bz[w3d.ipminfsapi-1:w3d.ipminfsapi-1+w3d.npfsapi]
    bx[:] = 0.
    by[:] = 0.
    bz[:] = 0.
    self.fetchbfrompositions(x,y,z,bx,by,bz)
#    print 'b',bx,by,bz

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
      if self.laser_focus is not None:
        z0 = self.laser_focus+self.laser_focus_velocity*top.time
        if self.laser_focus>0.:
          phase = (-(sqrt(xx**2+z0**2)-z0)/clight-top.time)*self.laser_frequency
        else:
          phase = ((sqrt(xx**2+z0**2)-z0)/clight-top.time)*self.laser_frequency
      else:
        phase = (xx*sin(self.laser_angle)/clight-top.time)*self.laser_frequency
    else:
      phase = 0.

    if (self.l_elaser_out_plane):
      field.Ez_in = self.laser_amplitude*field.laser_profile[:-1]*cos(phase)
    else:
      field.Bz_in = self.laser_amplitude*field.laser_profile[:-1]*cos(phase)

  def solve2ndhalf(self):
    if self.l_verbose:print 'solve 2nd half'
    if top.dt<>self.dtinit:raise('Time step has been changed since initialization of EM3D.')
    # --- Set nxl and nyl if using large stencil
#    if(not self.l_onegrid):
#      project_j(self.field,self.fpatchcoarse,self.fpatchfine)
#    for field in fields:
#      if field.l_uselargestencil and field.nxl<>field.nx:
#        field.nxl=field.nx
#        field.nyl=field.ny
#        field.gchange()
#      self.add_laser(field)
#      grimax(field)
    node2yee3d(self.block.core.yf)
    dt = top.dt/self.nfield_subcycle
    push_em3d_bf(self.block,dt,2,self.l_pushf)
    if self.l_verbose:print 'solve 2nd half done'
#    yee2node3d(self.block.core.yf)
 #   self.move_window_fields()
 #   for field in fields:
 #     griuni(field)
 #   if not self.l_onegrid:set_substitute_fields(self.field,self.fpatchcoarse,self.fpatchfine)

  def solve(self):
    if self.l_verbose:print 'solve 1st half'
    if top.dt<>self.dtinit:raise('Time step has been changed since initialization of EM3D.')
    # --- Set nxl and nyl if using large stencil
#    if(not self.l_onegrid):
#      project_j(self.field,self.fpatchcoarse,self.fpatchfine)
#    for field in fields:
#      if field.l_uselargestencil and field.nxl<>field.nx:
#        field.nxl=field.nx
#        field.nyl=field.ny
#        field.gchange()
#      self.add_laser(field)
#      grimax(field)
#    node2yee3d(self.block.core.yf)
    dt = top.dt/self.nfield_subcycle
    push_em3d_eef(self.block,dt,0,self.l_pushf)
    push_em3d_bf(self.block,dt,1,self.l_pushf)
    for i in range(self.nfield_subcycle-1):
      push_em3d_bf(self.block,dt,2,self.l_pushf)
      push_em3d_eef(self.block,dt,0,self.l_pushf)
      push_em3d_bf(self.block,dt,1,self.l_pushf)
    yee2node3d(self.block.core.yf)
 #   self.move_window_fields()
 #   for field in fields:
 #     griuni(field)
 #   if not self.l_onegrid:set_substitute_fields(self.field,self.fpatchcoarse,self.fpatchfine)

  ##########################################################################
  # Define the basic plot commands
  def genericfp(self,data,title,l_transpose=true,direction=2,slice=w3d.nz/2,**kw):
    f=self.block.core.yf
    if direction==0:
      settitles(title,'Y','Z','t = %gs'%(top.time))
      if l_transpose:
        kw.setdefault('xmin',f.zmin)
        kw.setdefault('xmax',f.zmax)
        kw.setdefault('ymin',f.ymin)
        kw.setdefault('ymax',f.ymax)
      else:
        kw.setdefault('xmin',f.ymin)
        kw.setdefault('xmax',f.ymax)
        kw.setdefault('ymin',f.zmin)
        kw.setdefault('ymax',f.zmax)
#      data=data[slice+f.nxguard,f.nyguard:-f.nyguard,f.nzguard:-f.nzguard]
      data=data[slice,:,:]
    if direction==1:
      settitles(title,'Z','X','t = %gs'%(top.time))
      if l_transpose:
        kw.setdefault('xmin',f.zmin)
        kw.setdefault('xmax',f.zmax)
        kw.setdefault('ymin',f.xmin)
        kw.setdefault('ymax',f.xmax)
      else:
        kw.setdefault('xmin',f.xmin)
        kw.setdefault('xmax',f.xmax)
        kw.setdefault('ymin',f.zmin)
        kw.setdefault('ymax',f.zmax)
#      data=data[f.nxguard:-f.nxguard,slice+f.nyguard,f.nzguard:-f.nzguard]
      data=data[:,slice,:]
    if direction==2:
      settitles(title,'X','Y','t = %gs'%(top.time))
      if l_transpose:
        kw.setdefault('xmin',f.ymin)
        kw.setdefault('xmax',f.ymax)
        kw.setdefault('ymin',f.xmin)
        kw.setdefault('ymax',f.xmax)
      else:
        kw.setdefault('xmin',f.xmin)
        kw.setdefault('xmax',f.xmax)
        kw.setdefault('ymin',f.ymin)
        kw.setdefault('ymax',f.ymax)
#      data=data[f.nxguard:-f.nxguard,f.nyguard:-f.nyguard,slice+f.nzguard]
      data=data[:,:,slice]
    ppgeneric(grid=data,**kw)
      
  def gatherarray(self,g):
    f=self.field
    if me<npes-1:
      return transpose(gatherarray(transpose(g[f.nxguard:-f.nxguard,f.nyguard:-f.nyguard,f.nzguard:-f.nzguard-1],[2,1,0])),[2,1,0])
    else:
      return transpose(gatherarray(transpose(g[f.nxguard:-f.nxguard,f.nyguard:-f.nyguard,f.nzguard:-f.nzguard],[2,1,0])),[2,1,0])
      
  def fpex(self,**kw):
      self.genericfp(self.gatherarray(self.field.Ex),'E_x',**kw)

  def fpey(self,**kw):
      self.genericfp(self.gatherarray(self.field.Ey),'E_y',**kw)

  def fpez(self,**kw):
      self.genericfp(self.gatherarray(self.field.Ez),'E_z',**kw)

  def fpbx(self,**kw):
      self.genericfp(self.gatherarray(self.field.Bx),'B_x',**kw)

  def fpby(self,**kw):
      self.genericfp(self.gatherarray(self.field.By),'B_y',**kw)

  def fpbz(self,**kw):
      self.genericfp(self.gatherarray(self.field.Bz),'B_z',**kw)

  def fpjx(self,**kw):
      self.genericfp(self.gatherarray(self.field.J[:,:,:,0]),'J_x',**kw)

  def fpjy(self,**kw):
      self.genericfp(self.gatherarray(self.field.J[:,:,:,1]),'J_y',**kw)

  def fpjz(self,**kw):
      self.genericfp(self.gatherarray(self.field.J[:,:,:,2]),'J_z',**kw)

  def fprho(self,**kw):
      self.genericfp(self.gatherarray(self.field.Rho),'Rho',**kw)

  def fpf(self,**kw):
      self.genericfp(self.gatherarray(self.field.F),'F',**kw)

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

  def getjx(self):
      return self.gatherarray(self.field.J[:,:,:,0])

  def getjy(self):
      return self.gatherarray(self.field.J[:,:,:,1])

  def getjz(self):
      return self.gatherarray(self.field.J[:,:,:,2])

  def getex(self):
      return self.gatherarray(self.field.Ex)

  def getey(self):
      return self.gatherarray(self.field.Ey)
        
  def getez(self):
      return self.gatherarray(self.field.Ez)
        
  def getbx(self):
      return self.gatherarray(self.field.Bx)

  def getby(self):
      return self.gatherarray(self.field.By)
        
  def getbz(self):
      return self.gatherarray(self.field.Bz)
        
  def getrho(self):
      return self.gatherarray(self.field.Rho)

  def getf(self):
      return self.gatherarray(self.field.F)

def allocatesf(f):
    f.syf = EM3D_SPLITYEEFIELDtype()
    f.fieldtype = f.syf.fieldtype
    f.proc=me

def allocatef(f):
    f.yf = EM3D_YEEFIELDtype()
    f.fieldtype = f.yf.fieldtype
    f.proc=me

def pyinit_3dem_block(nx, ny, nz, nbndx, nbndy, nbndz, nxguard, nyguard, nzguard, dt, dx, dy, dz, clight, mu0, xmin, ymin, zmin, xlb, ylb, zlb, xrb, yrb, zrb):
  
  l_1d_decomposition=true
  
  b = EM3D_BLOCKtype()
  b.core = EM3D_FIELDtype()
  b.core.yf = EM3D_YEEFIELDtype()
  b.core.fieldtype = b.core.yf.fieldtype
  b.core.proc = me

  if npes>1:
    if l_1d_decomposition:
      if me==0:
        zrb = em3d.otherproc
        if zlb==periodic:
          zlb = em3d.otherproc
      elif me==npes-1:
        zlb = em3d.otherproc
        if zrb==periodic:
          zrb = em3d.otherproc
      else:
        zlb=zrb=em3d.otherproc

  b.nx = nx
  b.ny = ny
  b.nz = nz
  b.nbndx = nbndx
  b.nbndy = nbndy
  b.nbndz = nbndz
  b.xmin = xmin
  b.ymin = ymin
  b.zmin = zmin
  b.dx = dx
  b.dy = dy
  b.dz = dz
  b.xmax = xmin+dx*nx
  b.ymax = ymin+dy*ny
  b.zmax = zmin+dz*nz
  b.dxi = 1./dx
  b.dyi = 1./dy
  b.dzi = 1./dz
  b.xlbnd = xlb
  b.xrbnd = xrb
  b.ylbnd = ylb
  b.yrbnd = yrb
  b.zlbnd = zlb
  b.zrbnd = zrb

  b.core.yf.nx = nx
  b.core.yf.ny = ny
  b.core.yf.nz = nz
  b.core.yf.nxguard = nxguard
  b.core.yf.nyguard = nyguard
  b.core.yf.nzguard = nzguard
  b.core.yf.xmin = xmin
  b.core.yf.ymin = ymin
  b.core.yf.zmin = zmin
  b.core.yf.dx = dx
  b.core.yf.dy = dy
  b.core.yf.dz = dz
  b.core.yf.xmax = xmin+dx*nx
  b.core.yf.ymax = ymin+dy*ny
  b.core.yf.zmax = zmin+dz*nz
  b.core.yf.dxi = 1./dx
  b.core.yf.dyi = 1./dy
  b.core.yf.dzi = 1./dz
  b.core.yf.clight = clight
  b.core.yf.mu0    = mu0

  b.core.yf.gchange()
  
  nnx=em3d.nn
  nny=em3d.nn
  nnz=em3d.nn
  smaxx=em3d.s_max_init
  smaxy=em3d.s_max_init
  smaxz=em3d.s_max_init
  sdeltax=em3d.s_delta
  sdeltay=em3d.s_delta
  sdeltaz=em3d.s_delta
  
# *** sides
# x
  b.sidexl = EM3D_FIELDtype()
  if xlb==openbc:
    allocatesf(b.sidexl)
    init_splitfield(b.sidexl.syf,nbndx,ny,nz,nxguard,nyguard,nzguard, dt, dx, dy, dz, clight,-1, 0, 0, nnx, smaxx, sdeltax, nny, smaxy, sdeltay, nnz, smaxz, sdeltaz)
  if xlb==periodic:
    b.sidexl=b.core
  b.sidexr = EM3D_FIELDtype()
  if xrb==openbc:
    allocatesf(b.sidexr)
    init_splitfield(b.sidexr.syf,nbndx,ny,nz,nxguard,nyguard,nzguard, dt, dx, dy, dz, clight, 1, 0, 0, nnx, smaxx, sdeltax, nny, smaxy, sdeltay, nnz, smaxz, sdeltaz)
  if xrb==periodic:
    b.sidexr=b.core
# y
  b.sideyl = EM3D_FIELDtype()
  if ylb==openbc:
    allocatesf(b.sideyl)
    init_splitfield(b.sideyl.syf,nx,nbndy,nz,nxguard,nyguard,nzguard, dt, dx, dy, dz, clight, 0,-1, 0, nnx, smaxx, sdeltax, nny, smaxy, sdeltay, nnz, smaxz, sdeltaz)
  if ylb==periodic:
    b.sideyl=b.core
  b.sideyr = EM3D_FIELDtype()
  if yrb==openbc:
    allocatesf(b.sideyr)
    init_splitfield(b.sideyr.syf,nx,nbndy,nz,nxguard,nyguard,nzguard, dt, dx, dy, dz, clight, 0, 1, 0, nnx, smaxx, sdeltax, nny, smaxy, sdeltay, nnz, smaxz, sdeltaz)
  if yrb==periodic:
    b.sideyr=b.core
# z
  b.sidezl = EM3D_FIELDtype()
  if zlb==openbc:
    allocatesf(b.sidezl)
    init_splitfield(b.sidezl.syf,nx,ny,nbndz,nxguard,nyguard,nzguard, dt, dx, dy, dz, clight, 0, 0,-1, nnx, smaxx, sdeltax, nny, smaxy, sdeltay, nnz, smaxz, sdeltaz)
  if zlb==periodic:
    b.sidezl=b.core
  if zlb==em3d.otherproc:
    allocatef(b.sidezl)
    if l_1d_decomposition:
      b.sidezl.proc=(me-1)%npes      
  b.sidezr = EM3D_FIELDtype()
  if zrb==openbc:
    allocatesf(b.sidezr)
    init_splitfield(b.sidezr.syf,nx,ny,nbndz,nxguard,nyguard,nzguard, dt, dx, dy, dz, clight, 0, 0, 1, nnx, smaxx, sdeltax, nny, smaxy, sdeltay, nnz, smaxz, sdeltaz)
  if zrb==periodic:
    b.sidezr=b.core
  if zrb==em3d.otherproc:
    allocatef(b.sidezr)
    if l_1d_decomposition:
      b.sidezr.proc=(me+1)%npes      
    
# *** edges
# xy
  b.edgexlyl = EM3D_FIELDtype()
  if xlb==openbc and ylb==openbc:
    allocatesf(b.edgexlyl)
    init_splitfield(b.edgexlyl.syf,nbndx,nbndy,nz,nxguard,nyguard,nzguard, dt, dx, dy, dz, clight,-1,-1, 0, nnx, smaxx, sdeltax, nny, smaxy, sdeltay, nnz, smaxz, sdeltaz)
# partial periodic bc to be completed 
#  if xlb==openbc and ylb==periodic:
#    b.edgexlyl=b.sidexl
#  if ylb==openbc and xlb==periodic:
#    b.edgexlyl=b.sideyl
  b.edgexryl = EM3D_FIELDtype()
  if xrb==openbc and ylb==openbc:
    allocatesf(b.edgexryl)
    init_splitfield(b.edgexryl.syf,nbndx,nbndy,nz,nxguard,nyguard,nzguard, dt, dx, dy, dz, clight, 1,-1, 0, nnx, smaxx, sdeltax, nny, smaxy, sdeltay, nnz, smaxz, sdeltaz)
  b.edgexlyr = EM3D_FIELDtype()
  if xlb==openbc and yrb==openbc:
    allocatesf(b.edgexlyr)
    init_splitfield(b.edgexlyr.syf,nbndx,nbndy,nz,nxguard,nyguard,nzguard, dt, dx, dy, dz, clight,-1, 1, 0, nnx, smaxx, sdeltax, nny, smaxy, sdeltay, nnz, smaxz, sdeltaz)
  b.edgexryr = EM3D_FIELDtype()
  if xrb==openbc and yrb==openbc:
    allocatesf(b.edgexryr)
    init_splitfield(b.edgexryr.syf,nbndx,nbndy,nz,nxguard,nyguard,nzguard, dt, dx, dy, dz, clight, 1, 1, 0, nnx, smaxx, sdeltax, nny, smaxy, sdeltay, nnz, smaxz, sdeltaz)
# xz
  b.edgexlzl = EM3D_FIELDtype()
  if xlb==openbc and zlb==openbc:
    allocatesf(b.edgexlzl)
    init_splitfield(b.edgexlzl.syf,nbndx,ny,nbndz,nxguard,nyguard,nzguard, dt, dx, dy, dz, clight,-1, 0,-1, nnx, smaxx, sdeltax, nny, smaxy, sdeltay, nnz, smaxz, sdeltaz)
  if zlb==em3d.otherproc:
    allocatesf(b.edgexlzl)
    if l_1d_decomposition:
      b.edgexlzl.proc=(me-1)%npes      
  b.edgexrzl = EM3D_FIELDtype()
  if xrb==openbc and ylb==openbc:
    allocatesf(b.edgexrzl)
    init_splitfield(b.edgexrzl.syf,nbndx,ny,nbndz,nxguard,nyguard,nzguard, dt, dx, dy, dz, clight, 1, 0,-1, nnx, smaxx, sdeltax, nny, smaxy, sdeltay, nnz, smaxz, sdeltaz)
  if zlb==em3d.otherproc:
    allocatesf(b.edgexrzl)
    if l_1d_decomposition:
      b.edgexrzl.proc=(me-1)%npes      
  b.edgexlzr = EM3D_FIELDtype()
  if xlb==openbc and zrb==openbc:
    allocatesf(b.edgexlzr)
    init_splitfield(b.edgexlzr.syf,nbndx,ny,nbndz,nxguard,nyguard,nzguard, dt, dx, dy, dz, clight,-1, 0, 1, nnx, smaxx, sdeltax, nny, smaxy, sdeltay, nnz, smaxz, sdeltaz)
  if zrb==em3d.otherproc:
    allocatesf(b.edgexlzr)
    if l_1d_decomposition:
      b.edgexlzr.proc=(me+1)%npes      
  b.edgexrzr = EM3D_FIELDtype()
  if xrb==openbc and zrb==openbc:
    allocatesf(b.edgexrzr)
    init_splitfield(b.edgexrzr.syf,nbndx,ny,nbndz,nxguard,nyguard,nzguard, dt, dx, dy, dz, clight, 1, 0, 1, nnx, smaxx, sdeltax, nny, smaxy, sdeltay, nnz, smaxz, sdeltaz)
  if zrb==em3d.otherproc:
    allocatesf(b.edgexrzr)
    if l_1d_decomposition:
      b.edgexrzr.proc=(me+1)%npes      
# yz
  b.edgeylzl = EM3D_FIELDtype()
  if ylb==openbc and zlb==openbc:
    allocatesf(b.edgeylzl)
    init_splitfield(b.edgeylzl.syf,nx,nbndy,nbndz,nxguard,nyguard,nzguard, dt, dx, dy, dz, clight, 0,-1,-1, nnx, smaxx, sdeltax, nny, smaxy, sdeltay, nnz, smaxz, sdeltaz)
  if zlb==em3d.otherproc:
    allocatesf(b.edgeylzl)
    if l_1d_decomposition:
      b.edgeylzl.proc=(me-1)%npes      
  b.edgeyrzl = EM3D_FIELDtype()
  if yrb==openbc and zlb==openbc:
    allocatesf(b.edgeyrzl)
    init_splitfield(b.edgeyrzl.syf,nx,nbndy,nbndz,nxguard,nyguard,nzguard, dt, dx, dy, dz, clight, 0, 1,-1, nnx, smaxx, sdeltax, nny, smaxy, sdeltay, nnz, smaxz, sdeltaz)
  if zlb==em3d.otherproc:
    allocatesf(b.edgeyrzl)
    if l_1d_decomposition:
      b.edgeyrzl.proc=(me-1)%npes      
  b.edgeylzr = EM3D_FIELDtype()
  if ylb==openbc and zrb==openbc:
    allocatesf(b.edgeylzr)
    init_splitfield(b.edgeylzr.syf,nx,nbndy,nbndz,nxguard,nyguard,nzguard, dt, dx, dy, dz, clight, 0,-1, 1, nnx, smaxx, sdeltax, nny, smaxy, sdeltay, nnz, smaxz, sdeltaz)
  if zrb==em3d.otherproc:
    allocatesf(b.edgeylzr)
    if l_1d_decomposition:
      b.edgeylzr.proc=(me+1)%npes      
  b.edgeyrzr = EM3D_FIELDtype()
  if yrb==openbc and zrb==openbc:
    allocatesf(b.edgeyrzr)
    init_splitfield(b.edgeyrzr.syf,nx,nbndy,nbndz,nxguard,nyguard,nzguard, dt, dx, dy, dz, clight, 0, 1, 1, nnx, smaxx, sdeltax, nny, smaxy, sdeltay, nnz, smaxz, sdeltaz)
  if zrb==em3d.otherproc:
    allocatesf(b.edgeyrzr)
    if l_1d_decomposition:
      b.edgeyrzr.proc=(me+1)%npes      

# *** corners
  b.cornerxlylzl = EM3D_FIELDtype()
  if xlb==openbc and ylb==openbc and zlb==openbc:
    allocatesf(b.cornerxlylzl)
    init_splitfield(b.cornerxlylzl.syf,nbndx,nbndy,nbndz,nxguard,nyguard,nzguard, dt, dx, dy, dz, clight,-1,-1,-1, nnx, smaxx, sdeltax, nny, smaxy, sdeltay, nnz, smaxz, sdeltaz)
  if zlb==em3d.otherproc:
    allocatesf(b.cornerxlylzl)
    if l_1d_decomposition:
      b.cornerxlylzl.proc=(me-1)%npes      
  b.cornerxrylzl = EM3D_FIELDtype()
  if xrb==openbc and ylb==openbc and zlb==openbc:
    allocatesf(b.cornerxrylzl)
    init_splitfield(b.cornerxrylzl.syf,nbndx,nbndy,nbndz,nxguard,nyguard,nzguard, dt, dx, dy, dz, clight, 1,-1,-1, nnx, smaxx, sdeltax, nny, smaxy, sdeltay, nnz, smaxz, sdeltaz)
  if zlb==em3d.otherproc:
    allocatesf(b.cornerxrylzl)
    if l_1d_decomposition:
      b.cornerxrylzl.proc=(me-1)%npes      
  b.cornerxlyrzl = EM3D_FIELDtype()
  if xlb==openbc and yrb==openbc and zlb==openbc:
    allocatesf(b.cornerxlyrzl)
    init_splitfield(b.cornerxlyrzl.syf,nbndx,nbndy,nbndz,nxguard,nyguard,nzguard, dt, dx, dy, dz, clight,-1, 1,-1, nnx, smaxx, sdeltax, nny, smaxy, sdeltay, nnz, smaxz, sdeltaz)
  if zlb==em3d.otherproc:
    allocatesf(b.cornerxlyrzl)
    if l_1d_decomposition:
      b.cornerxlyrzl.proc=(me-1)%npes      
  b.cornerxryrzl = EM3D_FIELDtype()
  if xrb==openbc and yrb==openbc and zlb==openbc:
    allocatesf(b.cornerxryrzl)
    init_splitfield(b.cornerxryrzl.syf,nbndx,nbndy,nbndz,nxguard,nyguard,nzguard, dt, dx, dy, dz, clight, 1, 1,-1, nnx, smaxx, sdeltax, nny, smaxy, sdeltay, nnz, smaxz, sdeltaz)
  if zlb==em3d.otherproc:
    allocatesf(b.cornerxryrzl)
    if l_1d_decomposition:
      b.cornerxryrzl.proc=(me-1)%npes      
  b.cornerxlylzr = EM3D_FIELDtype()
  if xlb==openbc and ylb==openbc and zrb==openbc:
    allocatesf(b.cornerxlylzr)
    init_splitfield(b.cornerxlylzr.syf,nbndx,nbndy,nbndz,nxguard,nyguard,nzguard, dt, dx, dy, dz, clight,-1,-1, 1, nnx, smaxx, sdeltax, nny, smaxy, sdeltay, nnz, smaxz, sdeltaz)
  if zrb==em3d.otherproc:
    allocatesf(b.cornerxlylzr)
    if l_1d_decomposition:
      b.cornerxlylzr.proc=(me+1)%npes      
  b.cornerxrylzr = EM3D_FIELDtype()
  if xrb==openbc and ylb==openbc and zrb==openbc:
    allocatesf(b.cornerxrylzr)
    init_splitfield(b.cornerxrylzr.syf,nbndx,nbndy,nbndz,nxguard,nyguard,nzguard, dt, dx, dy, dz, clight, 1,-1, 1, nnx, smaxx, sdeltax, nny, smaxy, sdeltay, nnz, smaxz, sdeltaz)
  if zrb==em3d.otherproc:
    allocatesf(b.cornerxrylzr)
    if l_1d_decomposition:
      b.cornerxrylzr.proc=(me+1)%npes      
  b.cornerxlyrzr = EM3D_FIELDtype()
  if xlb==openbc and yrb==openbc and zrb==openbc:
    allocatesf(b.cornerxlyrzr)
    init_splitfield(b.cornerxlyrzr.syf,nbndx,nbndy,nbndz,nxguard,nyguard,nzguard, dt, dx, dy, dz, clight,-1, 1, 1, nnx, smaxx, sdeltax, nny, smaxy, sdeltay, nnz, smaxz, sdeltaz)
  if zrb==em3d.otherproc:
    allocatesf(b.cornerxlyrzr)
    if l_1d_decomposition:
      b.cornerxlyrzr.proc=(me+1)%npes      
  b.cornerxryrzr = EM3D_FIELDtype()
  if xrb==openbc and yrb==openbc and zrb==openbc:
    allocatesf(b.cornerxryrzr)
    init_splitfield(b.cornerxryrzr.syf,nbndx,nbndy,nbndz,nxguard,nyguard,nzguard, dt, dx, dy, dz, clight, 1, 1, 1, nnx, smaxx, sdeltax, nny, smaxy, sdeltay, nnz, smaxz, sdeltaz)
  if zrb==em3d.otherproc:
    allocatesf(b.cornerxryrzr)
    if l_1d_decomposition:
      b.cornerxryrzr.proc=(me+1)%npes      

  return b

##############################################################################
# --- This can only be done after the class is defined.
#try:
#  psyco.bind(EM3D)
#except NameError:
#  pass

