"""Class for doing 3 D electromagnetic solver """
from warp import *
import operator

try:
  from Opyndx import *
  l_opyndx = 1
except:
  l_opyndx = 0
  
try:
  import psyco
except ImportError:
  pass

class EM3D(SubcycledPoissonSolver):
  
  __em3dinputs__ = []
  __flaginputs__ = {'l_apply_pml':true,'nbndx':10,'nbndy':10,'nbndz':10,
                    'nxguard':1,'nyguard':1,'nzguard':1,
                    'l_particles_weight':false,'l_usecoeffs':false,
                    'l_pushf':false,'l_pushpot':false,'l_verbose':false,
                    'laser_amplitude':1.,'laser_profile':None,
                    'laser_gauss_widthx':None,
                    'laser_gauss_widthy':None,
                    'laser_anglex':0.,'laser_angley':0.,
                    'laser_polangle':0.,
                    'laser_wavelength':None,'laser_wavenumber':None,
                    'laser_frequency':None,'laser_source_z':None,
                    'laser_focus':None,'laser_focus_velocity':0.,
                    'nfield_subcycle':1,'ncyclesperstep':None,
                    'l_2dxz':0,'npass_smooth':array([[0],[0],[0]]),
                    'alpha_smooth':array([[0.5],[0.5],[0.5]]),
                    'l_smooth_particle_fields':True,
                    'autoset_timestep':true,'dtcoef':1.,#0.99,
                    'deposit_energy_density':false,'refinement':None,
                    'l_force_nzlocal2nz':false,'inactive':false,
                    'stencil':false}

  def __init__(self,**kw):
    self.solveroff=False # flag to turn off the solver, for testing purpose
    assert top.grid_overlap<>2,"The EM solver needs top.grid_overlap==1!"
    self.grid_overlap = 1
    FieldSolver.__init__(self,kwdict=kw)
#    top.allspecl = true
    top.lcallfetchb = true
    top.lgridqnt = true
    self.zgrid=top.zgrid
    self.nzshifts=0
    # --- Save input parameters
#    self.processdefaultsfrompackage(EM3D.__w3dinputs__,w3d,kw)
#    self.processdefaultsfrompackage(EM3D.__topinputs__,top,kw)
    self.processdefaultsfrompackage(EM3D.__em3dinputs__,em3d,kw)
    self.processdefaultsfromdict(EM3D.__flaginputs__,kw)
    
    if self.inactive:self.isactive=False
        
    # --- When initializing self.field_coarse, there is some inconsistency 
    # --- with the way that FieldSolver.__init__ resize nzlocal in parallel. 
    # --- This flag takes care of it for now.
    if self.l_force_nzlocal2nz:
      self.nxlocal=self.nx
      self.nylocal=self.ny
      self.nzlocal=self.nz

    # --- set number of guard cells to appropriate value depending on order of 
    # --- current deposition, smoothing and stencil.
    self.npass_smooth  = array(self.npass_smooth)
    self.alpha_smooth  = array(self.alpha_smooth)
    minguards = array([1+int(top.depos_order.max(1)/2),self.npass_smooth.sum(1)]).max(0)
    if self.nxguard==1:self.nxguard = minguards[0]
    if self.nyguard==1:self.nyguard = minguards[1]
    if self.nzguard==1:self.nzguard = minguards[2]
    if self.stencil>0:
      if self.nxguard<2:self.nxguard=2
      if self.nyguard<2:self.nyguard=2
      if self.nzguard<2:self.nzguard=2
      self.npass_smooth = where(self.npass_smooth==0,1,self.npass_smooth)
    
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
    self.bounds = self.bounds.copy()

    # --- removes bounds from self.kw to prevent conflict with bounds created 
    # --- by the MR class, when adding MR children.
    if self.kw.has_key('bounds'):self.kw.pop('bounds')
    
    # --- If there are any remaning keyword arguments, raise an error.
    assert len(kw.keys()) == 0,"Bad keyword arguemnts %s"%kw.keys()

    # --- Calculate mesh sizes
    self.dx = (self.xmmax - self.xmmin)/self.nx
    self.xmesh = self.xmmin + arange(0,self.nx+1)*self.dx
    if self.l_2dxz:
      self.dy = 1.e36
      self.ymesh = array([0.])
    else:
      self.dy = (self.ymmax - self.ymmin)/self.ny
      self.ymesh = self.ymmin + arange(0,self.ny+1)*self.dy
    self.dz = (self.zmmax - self.zmmin)/self.nz
    self.zmesh = self.zmmin + arange(0,self.nz+1)*self.dz
#    self.zmeshlocal = self.zmminlocal + arange(0,self.nzlocal+1)*self.dz

    # --- set time step as a fraction of Courant condition
    # --- also set self.nfield_subcycle if top.dt over Courant condition times dtcoef
    if self.autoset_timestep:
      if self.stencil==0:
        dt=self.dtcoef/(clight*sqrt(1./self.dx**2+1./self.dy**2+1./self.dz**2))
      else:
        dt=self.dtcoef*min(self.dx,self.dy,self.dz)/clight 
      if top.dt==0.:
        top.dt=dt
      else:
        if top.dt>dt:
          self.nfield_subcycle=int(top.dt/dt)+1
    if self.ncyclesperstep is None:self.ncyclesperstep=self.nfield_subcycle
    self.dtinit = top.dt
    
    self.setbcparallel(0) # x
    self.setbcparallel(1) # y
    self.setbcparallel(2) # z
 
    if self.refinement is not None:
      ref=self.refinement
      self.field_coarse = self.__class__(l_force_nzlocal2nz=True,
                               nx=self.nx/ref[0],dx=self.dx*ref[0],
                               ny=self.ny/ref[1],dy=self.dy*ref[1],
                               nz=self.nz/ref[2],dz=self.dz*ref[2],
                               nxlocal=self.nxlocal/ref[0],
                               nylocal=self.nylocal/ref[1],
                               nzlocal=self.nzlocal/ref[2],
                               xmmin=self.xmmin,xmmax=self.xmmax,
                               ymmin=self.ymmin,ymmax=self.ymmax,
                               zmmin=self.zmmin,zmmax=self.zmmax,
                               xmminlocal=self.xmminlocal,xmmaxlocal=self.xmmaxlocal,
                               ymminlocal=self.ymminlocal,ymmaxlocal=self.ymmaxlocal,
                               zmminlocal=self.zmminlocal,zmmaxlocal=self.zmmaxlocal,
                               bounds=self.bounds,inactive=not self.isactive,
                               **self.kw)

    if self.xmmin>w3d.xmmaxlocal or self.xmmax<w3d.xmminlocal:return
    if self.ymmin>w3d.ymmaxlocal or self.ymmax<w3d.ymminlocal:return
    if self.zmmin>w3d.zmmaxlocal or self.zmmax<w3d.zmminlocal:return
    if self.refinement is not None:
      self.block_coarse = self.field_coarse.block
    
    # --- Create field and source arrays and other arrays.
    self.allocatefieldarrays()
    self.initializeconductors()
    # --- Handle laser inputs
    self.setuplaser()

    installbeforeloadrho(self.solve2ndhalf)

  def processdefaultsfrompackage(self,defaults,package,kw):
    for name in defaults:
      if name not in self.__dict__:
        #self.__dict__[name] = kw.pop(name,getattr(w3d,name)) # Python2.3
        self.__dict__[name] = kw.get(name,getattr(package,name))
      if kw.has_key(name): del kw[name]

  def processdefaultsfromdict(self,dict,kw):
    for name,defvalue in dict.iteritems():
      if name not in self.__dict__:
        #self.__dict__[name] = kw.pop(name,getattr(top,name)) # Python2.3
        self.__dict__[name] = kw.get(name,defvalue)
      if kw.has_key(name): del kw[name]

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
    return ((1+self.nxlocal,1+self.nylocal,1+self.nzlocal),
            (1+self.nxlocal+2*self.nxguard,
             1+self.nylocal+2*self.nyguard,
             1+self.nzlocal+2*self.nzguard))

  def allocatefieldarrays(self):
    self.block = pyinit_3dem_block(self.nxlocal, 
                                   self.nylocal, 
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
                                   self.xmminlocal, 
                                   self.ymminlocal, 
                                   self.zmminlocal, 
                                   int(self.bounds[0]), 
                                   int(self.bounds[2]), 
                                   int(self.bounds[4]), 
                                   int(self.bounds[1]), 
                                   int(self.bounds[3]), 
                                   int(self.bounds[5]),
                                   self.deposit_energy_density,
                                   self.refinement,
                                   self.l_pushf,
                                   self.stencil,
                                   self.npass_smooth,
                                   self.l_smooth_particle_fields,
                                   self.l_2dxz)
    self.fields = self.block.core.yf
    
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
      
    self.setuplaser_profile(self.fields)

    if self.laser_source_z is None:
      self.laser_source_z = w3d.zmmin
    self.laser_source_z = max(min(self.laser_source_z,w3d.zmmax),w3d.zmmin)

  def setuplaser_profile(self,f):
    # --- Check if laser_profile has a type, is a function, or a table
    self.laser_profile_func = None
    if self.laser_profile == 'gaussian':
      assert self.laser_gauss_widthx is not None,\
             "For a gaussian laser, the width in X must be specified using laser_gauss_widthx"
      assert self.laser_gauss_widthy is not None,\
             "For a gaussian laser, the width in Y must be specified using laser_gauss_widthy"
      f = self.block.core.yf
      if not self.l_2dxz:
        xxey,yyey = getmesh2d(f.xmin,f.dx,f.nx,f.ymin+0.5*f.dy,f.dy,f.ny-1)
        xxex,yyex = getmesh2d(f.xmin+0.5*f.dx,f.dx,f.nx-1,f.ymin,f.dy,f.ny)
      else:
        xxex = (arange(f.nx) + 0.5)*f.dx + f.xmin
        xxey = arange(f.nx+1)*f.dx + f.xmin
        yyex = yyey = 0.
      self.laser_profile = [exp(-((xxex/self.laser_gauss_widthx)**2+(yyex/self.laser_gauss_widthy)**2)/2.),
                            exp(-((xxey/self.laser_gauss_widthx)**2+(yyey/self.laser_gauss_widthy)**2)/2.)]
    elif operator.isSequenceType(self.laser_profile):
      assert len(self.laser_profile[:,0]) == f.nx,"The specified profile must be of length nx"
      assert len(self.laser_profile[0,:]) == f.ny,"The specified profile must be of length ny"
    elif callable(self.laser_profile):
      self.laser_profile_func = self.laser_profile

  def add_laser(self,field):
    if self.laser_profile is None: 
      e_inz_pos=-1
      return

    if self.laser_source_z>=w3d.zmminlocal and self.laser_source_z<=w3d.zmmaxlocal:
      self.block.core.yf.E_inz_pos = nint((self.laser_source_z-w3d.zmminlocal)/w3d.dz)
    else:
      return

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
      assert len(self.laser_profile[:,0]) == field.nx,"The specified profile must be of length nx"
      assert len(self.laser_profile[0,:]) == field.ny,"The specified profile must be of length ny"

#    if (self.l_elaser_out_plane):
#      xx = (arange(self.nx+4) - 0.5)*self.fields.dx+self.fields.xmin
#    else:
    f = self.block.core.yf
    if not self.l_2dxz:
      xxey,yyey = getmesh2d(f.xmin,f.dx,f.nx,f.ymin+0.5*f.dy,f.dy,f.ny-1)
      xxex,yyex = getmesh2d(f.xmin+0.5*f.dx,f.dx,f.nx-1,f.ymin,f.dy,f.ny)
    else:
      xxex = (arange(f.nx) + 0.5)*f.dx + f.xmin
      xxey = arange(f.nx+1)*f.dx + f.xmin
      yyex = yyey = 0.
    if self.laser_frequency is not None:
      if self.laser_focus is not None:
        z0 = self.laser_focus+self.laser_focus_velocity*top.time
        if self.laser_focus>0.:
          phaseex = (-(sqrt(xxex**2+yyex**2+z0**2)-z0)/clight-top.time)*self.laser_frequency
          phaseey = (-(sqrt(xxey**2+yyey**2+z0**2)-z0)/clight-top.time)*self.laser_frequency
        else:
          phaseex = ((sqrt(xxex**2+yyex**2+z0**2)-z0)/clight-top.time)*self.laser_frequency
          phaseey = ((sqrt(xxey**2+yyey**2+z0**2)-z0)/clight-top.time)*self.laser_frequency
      else:
        phaseex = ((xxex*sin(self.laser_anglex)+yyex*sin(self.laser_angley))/clight-top.time)*self.laser_frequency
        phaseey = ((xxey*sin(self.laser_anglex)+yyey*sin(self.laser_angley))/clight-top.time)*self.laser_frequency
    else:
      phase = 0.
    f = self.block.core.yf
    if self.l_2dxz:
      f.Ex_inz[f.jxmin:f.jxmax  ,f.jymin] = self.laser_amplitude*self.laser_profile[0]*cos(phaseex)*cos(self.laser_polangle)
      f.Ey_inz[f.jxmin:f.jxmax+1,f.jymin] = self.laser_amplitude*self.laser_profile[1]*cos(phaseey)*sin(self.laser_polangle)
    else:      
      f.Ex_inz[f.jxmin:f.jxmax  ,f.jymin:f.jymax+1] = self.laser_amplitude*self.laser_profile[0]*cos(phaseex)*cos(self.laser_polangle)
      f.Ey_inz[f.jxmin:f.jxmax+1,f.jymin:f.jymax  ] = self.laser_amplitude*self.laser_profile[1]*cos(phaseey)*sin(self.laser_polangle)
    
  def fetchfieldfrompositions(self,x,y,z,ex,ey,ez,bx,by,bz,js=0,pgroup=None):
    # --- This is called by fetchfield from fieldsolver.py
    n = len(x)
    if n == 0: return
    nox = top.depos_order[0,w3d.jsfsapi]
    noy = top.depos_order[1,w3d.jsfsapi]
    noz = top.depos_order[2,w3d.jsfsapi]
#    if (not (self.getconductorobject(top.fselfb[iselfb]).lcorrectede or
#    else:
    f = self.block.core.yf
    # --- fetch e
    if top.efetch[w3d.jsfsapi] in [1,3,5]:
      if self.l_2dxz:
        getf2dxz_n(n,x,z,ex,ey,ez,
                 f.xmin,f.zmin+self.zgrid,
                 f.dx,f.dz,
                 f.nx,f.ny,f.nz,
                 f.nxguard,f.nyguard,f.nzguard,
                 nox,noz,
                 f.Exp,f.Eyp,f.Ezp,
                 w3d.l4symtry)
      else:
        getf3d_n(n,x,y,z,ex,ey,ez,
                 f.xmin,f.ymin,f.zmin+self.zgrid,
                 f.dx,f.dy,f.dz,
                 f.nx,f.ny,f.nz,
                 f.nxguard,f.nyguard,f.nzguard,
                 nox,noy,noz,
                 f.Exp,f.Eyp,f.Ezp,
                 w3d.l4symtry)
    elif top.efetch[w3d.jsfsapi]==4:
      gete3d_n_energy_conserving(n,x,y,z,ex,ey,ez,
                    f.xmin,f.ymin,f.zmin+self.zgrid,
                    f.dx,f.dy,f.dz,
                    f.nx,f.ny,f.nz,
                    f.nxguard,f.nyguard,f.nzguard,
                    nox,noy,noz,
                    f.Exp,f.Eyp,f.Ezp,
                    w3d.l4symtry)
    # --- fetch b
    if top.efetch[w3d.jsfsapi] in [1,3,5]:
      if self.l_2dxz:
        getf2dxz_n(n,x,z,bx,by,bz,
                 f.xmin,f.zmin+self.zgrid,
                 f.dx,f.dz,
                 f.nx,f.ny,f.nz,
                 f.nxguard,f.nyguard,f.nzguard,
                 nox,noz,
                 f.Bxp,f.Byp,f.Bzp,
                 w3d.l4symtry)
      else:
        getf3d_n(n,x,y,z,bx,by,bz,
                 f.xmin,f.ymin,f.zmin+self.zgrid,
                 f.dx,f.dy,f.dz,
                 f.nx,f.ny,f.nz,
                 f.nxguard,f.nyguard,f.nzguard,
                 nox,noy,noz,
                 f.Bxp,f.Byp,f.Bzp,
                 w3d.l4symtry)
    elif top.efetch[w3d.jsfsapi]==4:
      getb3d_n_energy_conserving(n,x,y,z,bx,by,bz,
                    f.xmin,f.ymin,f.zmin+self.zgrid,
                    f.dx,f.dy,f.dz,
                    f.nx,f.ny,f.nz,
                    f.nxguard,f.nyguard,f.nzguard,
                    nox,noy,noz,
                    f.Bxp,f.Byp,f.Bzp,
                    w3d.l4symtry)
      
  def fetchphifrompositions(self,x,z,phi):
    pass

  def setsourcep(self,js,pgroup,zgrid):
    if self.l_verbose:print 'setsourcep, species ',js
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
    w  = pgroup.sw[js]*pgroup.dtscale[js]
    if top.wpid==0:
      wfact = zeros((0,), 'd')
    else:
      wfact = top.pgroup.pid[i:i+n,top.wpid-1]
    w3d.jsfsapi=js
    self.setsourcepatposition(x,y,z,ux,uy,uz,gaminv,wfact,zgrid,q,w)

  def setsourcepatposition(self,x,y,z,ux,uy,uz,gaminv,wfact,zgrid,q,w):
    n = x.shape[0]
    if n == 0: return
    # --- call routine performing current deposition
    f = self.block.core.yf
    js = w3d.jsfsapi
    if top.wpid==0:
      wfact = ones((1,),'d')
      l_particles_weight = false
    else:
      l_particles_weight = true
    if not self.deposit_energy_density:
      if self.l_2dxz:
        j = self.fields.J[:,self.fields.nyguard,:,:].copy()
        depose_jxjyjz_esirkepov_n_2d(j,n,
                                            x,z,ux,uy,uz,
                                            gaminv,wfact,q*w,
                                            f.xmin,f.zmin+self.zgrid,
                                            top.dt*top.pgroup.ndts[js],
                                            f.dx,f.dz,
                                            f.nx,f.nz,
                                            f.nxguard,f.nzguard,
                                            top.depos_order[0,js],
                                            top.depos_order[2,js],
                                            l_particles_weight,w3d.l4symtry)
        self.fields.Jarray[:,self.fields.nyguard,:,:,0]+=j[:,:,:].copy()
      else:
          depose_jxjyjz_esirkepov_n(self.fields.J,n,
                                            x,y,z,ux,uy,uz,
                                            gaminv,wfact,q*w,
                                            f.xmin,f.ymin,f.zmin+self.zgrid,
                                            top.dt*top.pgroup.ndts[js],
                                            f.dx,f.dy,f.dz,
                                            f.nx,f.ny,f.nz,
                                            f.nxguard,f.nyguard,f.nzguard,
                                            top.depos_order[0,js],
                                            top.depos_order[1,js],
                                            top.depos_order[2,js],
                                            l_particles_weight,w3d.l4symtry)
    else:
          depose_jxjyjz_pxpypz_esirkepov_linear_serial(self.fields.J,self.fields.Mp,n,
                                            x,y,z,ux,uy,uz,
                                            gaminv,wfact,q*w,m,
                                            f.xmin,f.ymin,f.zmin+self.zgrid,
                                            top.dt*top.pgroup.ndts[js],
                                            f.dx,f.dy,f.dz,
                                            f.nx,f.ny,f.nz,
                                            f.nxguard,f.nyguard,f.nzguard,
                                            l_particles_weight,
                                            top.lrelativ)
    if self.l_pushf:
      depose_rho_n(self.fields.Rho,n,
                               x,y,z,
                               wfact,q*w,
                               f.xmin,f.ymin,f.zmin+self.zgrid,
                               f.dx,f.dy,f.dz,
                               f.nx,f.ny,f.nz,
                               f.nxguard,f.nyguard,f.nzguard,
                               top.depos_order[0,js],
                               top.depos_order[1,js],
                               top.depos_order[2,js],
                               l_particles_weight,w3d.l4symtry)

  def allocatedataarrays(self):
    if self.l_verbose:print 'allocatedataarrays'
    # --- reallocate Jarray if needed
    if self.fields.ntimes<>top.nsndts:
      self.fields.ntimes=top.nsndts
      self.fields.gchange()

  def zerosourcep(self):
    if self.l_verbose:print 'zerosourcep',self

    # --- copy rho to rhoold if needed
    if self.l_pushf and self.ncyclesperstep>1:
      self.fields.Rhoold = self.fields.Rho.copy()
      if self.refinement is not None:
        self.field_coarse.fields.Rhoold = self.field_coarse.fields.Rho.copy()

    # --- zero proper portion of Jarray
    for indts in range(top.nsndts):
      if top.ldts[indts]:
        self.fields.Jarray[...,indts] = 0.
        if self.refinement is not None:
          self.field_coarse.fields.Jarray[...,indts] = 0.
        if self.l_pushf:
          self.fields.Rhoarray[...,indts] = 0.
          if self.refinement is not None:
            self.field_coarse.fields.Rho[...,indts] = 0.
        if self.deposit_energy_density:
          self.fields.Mp[...] = 0.
  
  def setsourcepforparticles(self,isourcepndtscopies,indts,iselfb):
    if self.l_verbose:print 'setsourcepforparticles'
    # --- point J array to proper Jarray slice
    self.fields.J = self.fields.Jarray[:,:,:,:,indts]
    if self.l_pushf: self.fields.Rho = self.fields.Rhoarray[:,:,:,indts]

  def finalizesourcep(self):
    if self.sourcepfinalized: return
    self.sourcepfinalized = 1
    if self.l_verbose:print 'finalizesourcep'
    # --- add slices
    if top.nsndts>1:
      if self.refinement is not None:raise('Error in finalizesourcep:nsndts>1 not fully implemented yet with MR')
      for indts in range(top.nsndts-2,-1,-1):
        if top.ldts[indts]:
          add_current_slice_3d(self.fields,indts+1)
          if self.l_pushf:add_rho_slice_3d(self.fields,indts+1)
    # --- smooth current density 
    if any(self.npass_smooth)>0:self.smoothdensity()
    self.aftersetsourcep()
    self.applysourceboundaryconditions()

    if self.l_verbose:print 'finalizesourcep done'

  def apply_rho_bc(self,block):
    # --- point Rho to first slice of Rhoarray
    block.core.yf.Rho = block.core.yf.Rhoarray[:,:,:,0]
    em3d_applybc_rho(block.core.yf, 
                     block.xlbnd,
                     block.xrbnd,
                     block.ylbnd,
                     block.yrbnd,
                     block.zlbnd,
                     block.zrbnd)
    em3d_exchange_rho(block)

  def apply_current_bc(self,block):
    # --- point J to first slice of Jarray
    block.core.yf.J = block.core.yf.Jarray[:,:,:,:,0]
    em3d_applybc_j(block.core.yf, 
                   block.xlbnd,
                   block.xrbnd,
                   block.ylbnd,
                   block.yrbnd,
                   block.zlbnd,
                   block.zrbnd)
    em3d_exchange_j(block)

  def applysourceboundaryconditions(self):
    # --- apply boundary condition on current
    self.apply_current_bc(self.block)
    if self.refinement is not None:
      self.apply_current_bc(self.block_coarse)
    if self.l_pushf:
      self.apply_rho_bc(self.block)
      if self.refinement is not None:
        self.apply_rho_bc(self.block_coarse)

    # --- point J to first slice of Jarray
    self.fields.J = self.fields.Jarray[:,:,:,:,0]
    if self.l_pushf:self.fields.Rho = self.fields.Rhoarray[:,:,:,0]
    
  def smoothdensity(self):
    nx,ny,nz = shape(self.fields.J[...,0])
    nsm = shape(self.npass_smooth)[1]
    for js in range(nsm):
      smooth3d_121(self.fields.J[...,0],nx-1,ny-1,nz-1,self.npass_smooth[:,js],self.alpha_smooth[:,js])
      smooth3d_121(self.fields.J[...,1],nx-1,ny-1,nz-1,self.npass_smooth[:,js],self.alpha_smooth[:,js])
      smooth3d_121(self.fields.J[...,2],nx-1,ny-1,nz-1,self.npass_smooth[:,js],self.alpha_smooth[:,js])
      if self.l_pushf:
        smooth3d_121(self.fields.Rho[...],nx-1,ny-1,nz-1,self.npass_smooth[:,js],self.alpha_smooth[:,js])

  def smoothfields(self):
    nx,ny,nz = shape(self.fields.J[...,0])
    nsm = shape(self.npass_smooth)[1]
    for js in range(nsm):
      smooth3d_121(self.fields.Exp,nx-1,ny-1,nz-1,self.npass_smooth[:,js],self.alpha_smooth[:,js])
      smooth3d_121(self.fields.Eyp,nx-1,ny-1,nz-1,self.npass_smooth[:,js],self.alpha_smooth[:,js])
      smooth3d_121(self.fields.Ezp,nx-1,ny-1,nz-1,self.npass_smooth[:,js],self.alpha_smooth[:,js])
      smooth3d_121(self.fields.Bxp,nx-1,ny-1,nz-1,self.npass_smooth[:,js],self.alpha_smooth[:,js])
      smooth3d_121(self.fields.Byp,nx-1,ny-1,nz-1,self.npass_smooth[:,js],self.alpha_smooth[:,js])
      smooth3d_121(self.fields.Bzp,nx-1,ny-1,nz-1,self.npass_smooth[:,js],self.alpha_smooth[:,js])

  def fetche(self,*args,**kw):
    SubcycledPoissonSolver.fetchfield(self,*args,**kw)

  def loadrho(self,lzero=None,pgroups=None,**kw):
    if self.l_verbose:print 'loadrho',self
    SubcycledPoissonSolver.loadsource(self,lzero,pgroups,**kw)
    
  def fetchphi(self):
    pass

  def fetcha(self):
    pass

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
#    self.mgmaxlevels=0
#    self.mglevels=0

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
      nxlocal,nylocal,nzlocal,nx,ny,nz = self.nxp,self.nyp,self.nzp,self.nx,self.ny,self.nz
      xmmin,xmmax = self.xmminp,self.xmmaxp
      ymmin,ymmax = self.ymminp,self.ymmaxp
      zmmin,zmmax = self.zmminp,self.zmmaxp
      mgmaxlevels = 1
    else:
      # --- Get relativistic longitudinal scaling factor
      # --- This is quite ready yet.
      beta = fselfb/clight
      zscale = 1./sqrt((1.-beta)*(1.+beta))
      nxlocal,nylocal,nzlocal,nx,ny,nz = self.nxlocal,self.nylocal,self.nzlocal,self.nx,self.ny,self.nz
      xmmin,xmmax = self.xmmin,self.xmmax
      ymmin,ymmax = self.ymmin,self.ymmax
      zmmin,zmmax = self.zmmin,self.zmmax
      mgmaxlevels = None

    mgmaxlevels=1
    # to be updated
    installconductors(conductor,xmin,xmax,ymin,ymax,zmin,zmax,dfill,
                      self.getzgrid(),
                      nx,ny,nz,
                      nxlocal,nylocal,nzlocal,
                      xmmin,xmmax,ymmin,ymmax,zmmin,zmmax,
                      zscale,self.l2symtry,self.l4symtry,
                      installrz=0,
                      solvergeom=self.solvergeom,conductors=conductorobject,
                      mgmaxlevels=mgmaxlevels,
                      my_index=self.my_index,nslaves=self.nslaves,
                      izfsslave=self.izfsslave,nzfsslave=self.nzfsslave)

#installconductors(a,xmin=None,xmax=None,ymin=None,ymax=None,
#                        zmin=None,zmax=None,dfill=2.,
#                        zbeam=None,
#                        nx=None,ny=None,nz=None,
#                        nxlocal=None,nylocal=None,nzlocal=None,
#                        xmmin=None,xmmax=None,ymmin=None,ymmax=None,
#                        zmmin=None,zmmax=None,zscale=1.,l2symtry=None,l4symtry=None,
#                        installrz=None,gridmode=1,solvergeom=None,
#                        conductors=None,gridrz=None,mgmaxlevels=None,
#                        decomp=None):
                        
    self.fields.nconds = self.getconductorobject().interior.n
    self.fields.nxcond = self.fields.nx
    self.fields.nycond = self.fields.ny
    self.fields.nzcond = self.fields.nz
    self.fields.gchange()
    self.fields.incond=False
    if self.fields.nconds>0:
      set_incond(self.fields,self.fields.nconds,int(self.getconductorobject().interior.indx))
    if self.block.xlbnd==openbc:self.fields.incond[:3,:,:]=False
    if self.block.xrbnd==openbc:self.fields.incond[-3:,:,:]=False
    if self.block.ylbnd==openbc:self.fields.incond[:,:3,:]=False
    if self.block.yrbnd==openbc:self.fields.incond[:,-3:,:]=False
    if self.block.zlbnd==openbc:self.fields.incond[:,:,:3]=False
    if self.block.zrbnd==openbc:self.fields.incond[:,:,-3:]=False
    
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

  def optimizeconvergence(self,resetpasses=1):
    pass

  def move_window_fields(self):
    if (abs(top.zgrid-self.zgrid)<0.5*self.dz):return 
#    if not self.l_moving_window or ((top.it%self.ndelta_t)!=0): return
#    if top.time < self.tmin_moving_main_window: return
    shift_em3dblock_ncells_z(self.block,1)
    self.zgrid+=self.dz
    self.nzshifts+=1
   
  def solve2ndhalf(self):
    if self.solveroff:return
    if self.l_verbose:print 'solve 2nd half',self
    if top.dt<>self.dtinit:raise('Time step has been changed since initialization of EM3D.')
    if top.efetch[0]<>4:node2yee3d(self.block.core.yf)
    dt = top.dt/self.nfield_subcycle
    push_em3d_bf(self.block,dt,2,self.l_pushf,self.l_pushpot)
    if self.l_verbose:print 'solve 2nd half done'
    self.move_window_fields()

  def dosolve(self,iwhich=0,*args):
    if self.solveroff:return
    if any(top.fselfb<>0.):raise('Error:EM solver does not work if fselfb<>0.')
    if self.l_verbose:print 'solve 1st half'
    if top.dt<>self.dtinit:raise('Time step has been changed since initialization of EM3D.')
    dt = top.dt/self.nfield_subcycle
    if self.l_pushf and self.ncyclesperstep>1:
      w = 1./self.ncyclesperstep
      self.fields.Rho = (1.-w)*self.fields.Rhoold + w*self.fields.Rhoarray[...,0] 
    push_em3d_eef(self.block,dt,0,self.l_pushf,self.l_pushpot)
    self.add_laser(self.fields)
    for i in range(self.ncyclesperstep-1):
      push_em3d_bf(self.block,dt,0,self.l_pushf,self.l_pushpot)
      if self.l_pushf:
        w = float(i+2)/self.ncyclesperstep
        self.fields.Rho = (1.-w)*self.fields.Rhoold + w*self.fields.Rhoarray[...,0] 
      push_em3d_eef(self.block,dt,0,self.l_pushf,self.l_pushpot)
    push_em3d_bf(self.block,dt,1,self.l_pushf,self.l_pushpot)
    if not all(top.efetch==top.efetch[0]):raise('Error:top.efetch must have same value for every species when using EM solver.')
    if top.efetch[0]<>4:yee2node3d(self.block.core.yf)
    if self.l_smooth_particle_fields:
      if self.refinement is None and any(self.npass_smooth>0):
        self.fields.Exp[...] = self.fields.Ex[...]
        self.fields.Eyp[...] = self.fields.Ey[...]
        self.fields.Ezp[...] = self.fields.Ez[...]
        self.fields.Bxp[...] = self.fields.Bx[...]
        self.fields.Byp[...] = self.fields.By[...]
        self.fields.Bzp[...] = self.fields.Bz[...]
      if any(self.npass_smooth>0):self.smoothfields()
    
    if self.l_verbose:print 'solve 1st half done'
    if self.refinement is not None:
      self.__class__.__bases__[1].dosolve(self.field_coarse)

  ##########################################################################
  # Define the basic plot commands
  def genericpfem3d(self,data,title,titles=True,l_transpose=false,direction=None,slice=None,l_abs=False,
                    origins=None,deltas=None,display=1,scale=None,camera=None,l_box=1,
                    interactive=0,labels=['X','Y','Z'],
                    adjust=0,labelscale=0.5,color='auto',ncolor=None,cmin=None,cmax=None,
                    procs=None,
                    **kw):
    if direction is None and not l_opyndx:direction=2
    if kw.has_key('view'):
      view=kw['view']
    else:
      view=1
    if kw.has_key('xscale'):
      xscale=kw['xscale']
    else:
      xscale=1
    if kw.has_key('yscale'):
      yscale=kw['yscale']
    else:
      yscale=1
    if kw.has_key('zscale'):
      zscale=kw['zscale']
    else:
      zscale=1
    f=self.block.core.yf
    nxd,nyd,nzd=shape(data)
    if direction==0:
      if slice is None:
        if self.l4symtry:
          slice=0
        else:
          slice=self.nx/2
      xslice = w3d.xmmin+slice*w3d.dx
      slice = nint((xslice-self.block.xmin)/self.block.dx)
      if slice<0 or slice>nxd-1:
        data=None
        xmin=xmax=ymin=ymax=0.
      else:
        if l_transpose:
          xmin=self.block.zmin+self.zgrid
          xmax=self.block.zmax+self.zgrid
          ymin=self.block.ymin
          ymax=self.block.ymax
        else:
          xmin=self.block.ymin
          xmax=self.block.ymax
          ymin=self.block.zmin+self.zgrid
          ymax=self.block.zmax+self.zgrid
        if l_transpose:
          data=transpose(data[slice,:,:])
        else:
          data=data[slice,:,:]
    if direction==1:
      if slice is None:
        if self.l4symtry:
          slice=0
        else:
          slice=self.ny/2
      yslice = w3d.ymmin+slice*w3d.dy
      slice = nint((yslice-self.block.ymin)/self.block.dy)
      if slice<0 or slice>nyd-1:
        data=None
        xmin=xmax=ymin=ymax=0.
      else:
        if l_transpose:
          xmin=self.block.zmin+self.zgrid
          xmax=self.block.zmax+self.zgrid
          ymin=self.block.xmin
          ymax=self.block.xmax
        else:
          xmin=self.block.xmin
          xmax=self.block.xmax
          ymin=self.block.zmin+self.zgrid
          ymax=self.block.zmax+self.zgrid
        if l_transpose:
          data=transpose(data[:,slice,:])
        else:
          data=data[:,slice,:]
    if direction==2:
      if slice is None:slice=self.nz/2
      zslice = w3d.zmmin+slice*w3d.dz
      slice = nint((zslice-self.block.zmin)/self.block.dz)
      print zslice,slice,self.block.nz
      if slice<0 or slice>nzd-1:
        data=None
        xmin=xmax=ymin=ymax=0.
      else:
        if l_transpose:
          xmin=self.block.ymin
          xmax=self.block.ymax
          ymin=self.block.xmin
          ymax=self.block.xmax
        else:
          xmin=self.block.xmin
          xmax=self.block.xmax
          ymin=self.block.ymin
          ymax=self.block.ymax
        if l_transpose:
          data=transpose(data[:,:,slice])
        else:
          data=data[:,:,slice]
    if l_abs and data is not None:data=abs(data)
    if cmin is None:
      if data is None:
        mindata = 1.e36
      else:
        mindata = minnd(data)
      cmin=parallelmin(mindata)
    if cmax is None:
      if data is None:
        maxdata = -1.e36
      else:
        maxdata = maxnd(data)
      cmax=parallelmax(maxdata)
    if procs is None:procs=arange(npes)
    if direction in [0,1,2]:
      kw.setdefault('cmin',cmin)
      kw.setdefault('cmax',cmax)
      if me==0 and titles:
        if direction==0:
          if l_transpose:
            xtitle='Z';ytitle='Y'
          else:
            xtitle='Y';ytitle='Z'
        if direction==1:
          if l_transpose:
            xtitle='Z';ytitle='X'
          else:
            xtitle='X';ytitle='Z'
        if direction==2:
          if l_transpose:
            xtitle='Y';ytitle='X'
          else:
            xtitle='X';ytitle='Y'
        plsys(view)
        ptitles(title,xtitle,ytitle,'t = %gs'%(top.time))
      if type(procs) is type(1):procs=[procs]
      if me>0 and me in procs:
        mpi.send((xmin,xmax,ymin,ymax,data),0,3)
      else:
        if me in procs and data is not None:
          kw.setdefault('xmin',xmin)
          kw.setdefault('xmax',xmax)
          kw.setdefault('ymin',ymin)
          kw.setdefault('ymax',ymax)
          ppgeneric(grid=data,**kw)
        for i in range(1,npes):
          xminp,xmaxp,yminp,ymaxp,data=mpirecv(i,3)
          if data is not None:
            kw['xmin']=xminp
            kw['xmax']=xmaxp
            kw['ymin']=yminp
            kw['ymax']=ymaxp
            ppgeneric(grid=data,lcolorbar=0,**kw)
    else:
      xmin=self.block.xmin
      xmax=self.block.xmax
      ymin=self.block.ymin
      ymax=self.block.ymax
      zmin=self.block.zmin
      zmax=self.block.zmax
      if me>0 and me in procs:
        mpi.send((xmin,xmax,ymin,ymax,zmin,zmax,data),0,3)
      else:
        if ncolor is None:ncolor = 2
        dec = (cmax-cmin)/ncolor
        isos = cmin+arange(ncolor)*dec+dec/2
        origins = [xmin*xscale,ymin*yscale,zmin*zscale]
        deltas = [w3d.dx*xscale,w3d.dy*yscale,w3d.dz*zscale]
        e3d,colorbar = viewisosurface1(data,isos,color=color,display=0,
                        origins=origins,
                        deltas=deltas)
        dxob = [e3d,colorbar]
        for i in range(1,npes):
          xminp,xmaxp,yminp,ymaxp,zminp,zmaxp,data=mpirecv(i,3)
          origins = [xminp*xscale,yminp*yscale,zminp*zscale]
          e3d,colorbar = viewisosurface1(data,isos,color=color,display=0,
                          origins=origins,
                          deltas=deltas)
          dxob.append(e3d)
        if l_box:
          box = viewboundingbox(w3d.xmmin*xscale,
                                w3d.xmmax*xscale,
                                w3d.ymmin*yscale,
                                w3d.ymmax*yscale,
                                w3d.zmmin*zscale,
                                w3d.zmmax*zscale,color='yellow')
          dxob.append(box)
        dxob = DXCollect(dxob)
        if camera is None:
          camera = DXAutocamera(dxob,direction=[-1.,1.,-1.],width=100.,resolution=640,
                                aspect=2.,up=[0,1,0],perspective=1,angle=60.,
                                background='black')
        DXImage(dxob,camera=camera,
                     l_interactive=interactive,
                     labels=labels,
                     adjust=adjust,
                     scale=scale,
                     labelscale=labelscale)
    
  def getarray(self,g,guards=0):
    if guards:
      return g
    else:
      f=self.fields
      ox=oy=oz=0
      if self.block.xrbnd==em3d.otherproc:ox=1
      if self.block.yrbnd==em3d.otherproc:oy=1
      if self.block.zrbnd==em3d.otherproc:oz=1
      return g[f.nxguard:-f.nxguard-ox,f.nyguard:-f.nyguard-oy,f.nzguard:-f.nzguard-oz]
      
  def pfex(self,guards=0,**kw):
      self.genericpfem3d(self.getarray(self.fields.Ex,guards),'E_x',**kw)

  def pfey(self,**kw):
      self.genericpfem3d(self.getarray(self.fields.Ey,guards),'E_y',**kw)

  def pfez(self,**kw):
      self.genericpfem3d(self.getarray(self.fields.Ez,guards),'E_z',**kw)

  def pfbx(self,**kw):
      self.genericpfem3d(self.getarray(self.fields.Bx,guards),'B_x',**kw)

  def pfby(self,**kw):
      self.genericpfem3d(self.getarray(self.fields.By,guards),'B_y',**kw)

  def pfbz(self,**kw):
      self.genericpfem3d(self.getarray(self.fields.Bz,guards),'B_z',**kw)

  def pfexp(self,**kw):
      self.genericpfem3d(self.getarray(self.fields.Exp,guards),'Ep_x',**kw)

  def pfeyp(self,**kw):
      self.genericpfem3d(self.getarray(self.fields.Eyp,guards),'Ep_y',**kw)

  def pfezp(self,**kw):
      self.genericpfem3d(self.getarray(self.fields.Ezp,guards),'Ep_z',**kw)

  def pfbxp(self,**kw):
      self.genericpfem3d(self.getarray(self.fields.Bxp,guards),'Bp_x',**kw)

  def pfbyp(self,**kw):
      self.genericpfem3d(self.getarray(self.fields.Byp,guards),'Bp_y',**kw)

  def pfbzp(self,**kw):
      self.genericpfem3d(self.getarray(self.fields.Bzp,guards),'Bp_z',**kw)

  def pfjx(self,**kw):
      self.genericpfem3d(self.getarray(self.fields.J[:,:,:,0],guards),'J_x',**kw)

  def pfjy(self,**kw):
      self.genericpfem3d(self.getarray(self.fields.J[:,:,:,1],guards),'J_y',**kw)

  def pfjz(self,**kw):
      self.genericpfem3d(self.getarray(self.fields.J[:,:,:,2],guards),'J_z',**kw)

  def pfrho(self,**kw):
      self.genericpfem3d(self.getarray(self.fields.Rho,guards),'Rho',**kw)

  def pff(self,**kw):
      self.genericpfem3d(self.getarray(self.fields.F,guards),'F',**kw)

  def pfdive(self,**kw):
      self.genericpfem3d(self.getarray(self.getdive(),guards),'div(E)',**kw)

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

  def getjx(self,guards=0):
      return self.getarray(self.fields.J[:,:,:,0],guards)

  def getjy(self,guards=0):
      return self.getarray(self.fields.J[:,:,:,1],guards)

  def getjz(self,guards=0):
      return self.getarray(self.fields.J[:,:,:,2],guards)

  def getex(self,guards=0):
      return self.getarray(self.fields.Ex,guards)

  def getey(self,guards=0):
      return self.getarray(self.fields.Ey,guards)
        
  def getez(self,guards=0):
      return self.getarray(self.fields.Ez,guards)
        
  def getbx(self,guards=0):
      return self.getarray(self.fields.Bx,guards)

  def getby(self,guards=0):
      return self.getarray(self.fields.By,guards)
        
  def getbz(self,guards=0):
      return self.getarray(self.fields.Bz,guards)
        
  def getrho(self,guards=0):
      return self.getarray(self.fields.Rho,guards)

  def getf(self,guards=0):
      return self.getarray(self.fields.F,guards)

  def getdive(self,guards=0):
      dive = zeros(shape(self.fields.Ex),'d')
      f = self.fields
      if top.efetch[0]<>4:node2yee3d(f)
      dive[1:-1,1:-1,1:-1] = (f.Ex[1:-1,1:-1,1:-1]-f.Ex[:-2,1:-1,1:-1])/f.dx \
                           + (f.Ey[1:-1,1:-1,1:-1]-f.Ey[1:-1,:-2,1:-1])/f.dy \
                           + (f.Ez[1:-1,1:-1,1:-1]-f.Ez[1:-1,1:-1,:-2])/f.dz 
      if top.efetch[0]<>4:yee2node3d(f)
      return dive
      
  def sumdive(self):
    flist = [self.block.core.yf]
    if self.refinement is not None:
      flist += [self.field_coarse.block.core.yf]
    q = []
    for f in flist:
      if top.efetch[0]<>4:node2yee3d(f)
      Ex = f.Ex#self.getarray(f.Ex)
      Ey = f.Ey#self.getarray(f.Ey)
      Ez = f.Ez#self.getarray(f.Ez)
      q.append( (sum(Ex[1:-1,1:-1,1:-1]-Ex[:-2,1:-1,1:-1])/f.dx \
              +  sum(Ey[1:-1,1:-1,1:-1]-Ey[1:-1,:-2,1:-1])/f.dy \
              +  sum(Ez[1:-1,1:-1,1:-1]-Ez[1:-1,1:-1,:-2])/f.dz) \
              * eps0*(f.dx*f.dy*f.dz))
      if top.efetch[0]<>4:yee2node3d(f)
    return q

  def sumq(self):
    flist = [self.block.core.yf]
    if self.refinement is not None:
      flist += [self.field_coarse.block.core.yf]
    q = []
    for f in flist:
      q.append(sum(f.Rho)*(f.dx*f.dy*f.dz))
    return q

  def setbcparallelold(self):
    if top.nzprocs>1:
      down = top.procneighbors[0,2]
      up   = top.procneighbors[1,2]
      # --- check lower bound in z
      if self.bounds[5] == periodic:
        condu = up<>me
      else:
        condu = up>me
      if self.bounds[4] == periodic:
        condd = down<>me
      else:
        condd = down<me
      if condu:
        mpi.send(self.isactive,up)
      if condd:
        isactive,status = mpi.recv(down)
        if isactive and self.isactive:self.bounds[4]=em3d.otherproc    
      # --- check upper bound in z
      if condd:
        mpi.send(self.isactive,down)
      if condu:
        isactive,status = mpi.recv(up)
        if isactive and self.isactive:self.bounds[5]=em3d.otherproc    

  def setbcparallel(self,dir):
    if dir==0:nprocs = top.nxprocs
    if dir==1:nprocs = top.nyprocs
    if dir==2:nprocs = top.nzprocs
    ib = 2*dir
    if nprocs>1:
      down = top.procneighbors[0,dir]
      up   = top.procneighbors[1,dir]
      # --- check lower bound in z
      if self.bounds[ib+1] == periodic:
        condu = up<>me
      else:
        condu = up>me
      if self.bounds[ib] == periodic:
        condd = down<>me
      else:
        condd = down<me
      if condu:
        mpi.send(self.isactive,up)
      if condd:
        isactive,status = mpi.recv(down)
        if isactive and self.isactive:self.bounds[ib]=em3d.otherproc    
      # --- check upper bound in z
      if condd:
        mpi.send(self.isactive,down)
      if condu:
        isactive,status = mpi.recv(up)
        if isactive and self.isactive:self.bounds[ib+1]=em3d.otherproc    

def allocatesf(f,stencil):
    f.syf = EM3D_SPLITYEEFIELDtype()
    f.fieldtype = f.syf.fieldtype
    f.proc=me
    f.syf.stencil=stencil
    
def allocatef(f,stencil):
    f.yf = EM3D_YEEFIELDtype()
    f.fieldtype = f.yf.fieldtype
    f.proc=me
    f.yf.stencil=stencil
    
def pyinit_3dem_block(nx, ny, nz, 
                      nbndx, nbndy, nbndz, 
                      nxguard, nyguard, nzguard, 
                      dt, dx, dy, dz, 
                      clight, mu0, 
                      xmin, ymin, zmin, 
                      xlb, ylb, zlb, xrb, yrb, zrb, 
                      deposit_energy_density,
                      refinement,l_pushf,
                      stencil,
                      npass_smooth,
                      l_smooth_particle_fields,
                      l_2dxz):
  
  procxl = top.procneighbors[0,0]
  procxr = top.procneighbors[1,0]
  procyl = top.procneighbors[0,1]
  procyr = top.procneighbors[1,1]
  proczl = top.procneighbors[0,2]
  proczr = top.procneighbors[1,2]
  
  b = EM3D_BLOCKtype()
  b.core = EM3D_FIELDtype()
  b.core.yf = EM3D_YEEFIELDtype()
  b.core.fieldtype = b.core.yf.fieldtype
  b.core.proc = me
  b.core.yf.stencil=stencil
  
  b.nx = nx
  b.ny = ny
  b.nz = nz
  b.nbndx = nbndx
  b.nbndy = nbndy
  b.nbndz = nbndz
  b.nxguard = nxguard
  b.nyguard = nyguard
  b.nzguard = nzguard
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

  f=b.core.yf
  f.nx = nx
  f.ny = ny
  f.nz = nz
  if refinement is None and (all(npass_smooth==0) or not l_smooth_particle_fields):
    f.nxp = 0
    f.nyp = 0
    f.nzp = 0
  else:
    f.nxp = f.nx
    f.nyp = f.ny
    f.nzp = f.nz
  if not l_pushf:
    f.nxf = 0
    f.nyf = 0
    f.nzf = 0
  else:
    f.nxf = f.nx
    f.nyf = f.ny
    f.nzf = f.nz
  f.nxguard = nxguard
  f.nyguard = nyguard
  f.nzguard = nzguard
  # set min/max of cells positions with FORTRAN indexing
  f.ixmin = 0
  f.iymin = 0
  f.izmin = 0
  f.ixmax = f.nx
  f.iymax = f.ny
  f.izmax = f.nz
  f.ixming = -f.nxguard
  f.iyming = -f.nyguard
  f.izming = -f.nzguard
  f.ixmaxg = f.ixmax+f.nxguard
  f.iymaxg = f.iymax+f.nyguard
  f.izmaxg = f.izmax+f.nzguard
  # set min/max of cells positions with Python indexing
  f.jxmin = f.ixmin-f.ixming
  f.jymin = f.iymin-f.iyming
  f.jzmin = f.izmin-f.izming
  f.jxmax = f.ixmax-f.ixming
  f.jymax = f.iymax-f.iyming
  f.jzmax = f.izmax-f.izming
  f.jxming = 0
  f.jyming = 0
  f.jzming = 0
  f.jxmaxg = f.ixmaxg-f.ixming
  f.jymaxg = f.iymaxg-f.iyming
  f.jzmaxg = f.izmaxg-f.izming
  f.xmin = xmin
  f.ymin = ymin
  f.zmin = zmin
  f.dx = dx
  f.dy = dy
  f.dz = dz
  f.xmax = xmin+dx*nx
  f.ymax = ymin+dy*ny
  f.zmax = zmin+dz*nz
  f.dxi = 1./dx
  f.dyi = 1./dy
  f.dzi = 1./dz
  f.clight = clight
  f.mu0    = mu0
  f.l_2dxz = l_2dxz

  if deposit_energy_density:
    f.nxmp=f.nx
    f.nymp=f.ny
    f.nzmp=f.nz

  f.gchange()
  
  if refinement is None and  (all(npass_smooth==0) or not l_smooth_particle_fields):
    f.nxp = f.nx
    f.nyp = f.ny
    f.nzp = f.nz
    f.Exp = f.Ex
    f.Eyp = f.Ey
    f.Ezp = f.Ez
    f.Bxp = f.Bx
    f.Byp = f.By
    f.Bzp = f.Bz

  nnx=em3d.nn
  nny=em3d.nn
  nnz=em3d.nn
  smaxx=em3d.s_max_init
  smaxy=em3d.s_max_init
  smaxz=em3d.s_max_init
  sdeltax=em3d.s_delta
  sdeltay=em3d.s_delta
  sdeltaz=em3d.s_delta
  
# --- sides
# x
  b.sidexl = EM3D_FIELDtype()
  if xlb==openbc:
    allocatesf(b.sidexl,stencil)
    init_splitfield(b.sidexl.syf,nbndx,ny,nz,nxguard,nyguard,nzguard, dt, dx, dy, dz, clight,-1, 0, 0, nnx, smaxx, sdeltax, nny, smaxy, sdeltay, nnz, smaxz, sdeltaz, l_2dxz)
  if xlb==periodic:
    b.sidexl=b.core
  if xlb==em3d.otherproc:
    allocatef(b.sidexl,stencil)
    b.sidexl.proc=procxl
  b.sidexr = EM3D_FIELDtype()
  if xrb==openbc:
    allocatesf(b.sidexr,stencil)
    init_splitfield(b.sidexr.syf,nbndx,ny,nz,nxguard,nyguard,nzguard, dt, dx, dy, dz, clight, 1, 0, 0, nnx, smaxx, sdeltax, nny, smaxy, sdeltay, nnz, smaxz, sdeltaz, l_2dxz)
  if xrb==periodic:
    b.sidexr=b.core
  if xrb==em3d.otherproc:
    allocatef(b.sidexr,stencil)
    b.sidexr.proc=procxr
# y
  b.sideyl = EM3D_FIELDtype()
  if ylb==openbc:
    allocatesf(b.sideyl,stencil)
    init_splitfield(b.sideyl.syf,nx,nbndy,nz,nxguard,nyguard,nzguard, dt, dx, dy, dz, clight, 0,-1, 0, nnx, smaxx, sdeltax, nny, smaxy, sdeltay, nnz, smaxz, sdeltaz, l_2dxz)
  if ylb==periodic:
    b.sideyl=b.core
  if ylb==em3d.otherproc:
    allocatef(b.sideyl,stencil)
    b.sideyl.proc=procyl
  b.sideyr = EM3D_FIELDtype()
  if yrb==openbc:
    allocatesf(b.sideyr,stencil)
    init_splitfield(b.sideyr.syf,nx,nbndy,nz,nxguard,nyguard,nzguard, dt, dx, dy, dz, clight, 0, 1, 0, nnx, smaxx, sdeltax, nny, smaxy, sdeltay, nnz, smaxz, sdeltaz, l_2dxz)
  if yrb==periodic:
    b.sideyr=b.core
  if yrb==em3d.otherproc:
    allocatef(b.sideyr,stencil)
    b.sideyr.proc=procyr
# z
  b.sidezl = EM3D_FIELDtype()
  if zlb==openbc:
    allocatesf(b.sidezl,stencil)
    init_splitfield(b.sidezl.syf,nx,ny,nbndz,nxguard,nyguard,nzguard, dt, dx, dy, dz, clight, 0, 0,-1, nnx, smaxx, sdeltax, nny, smaxy, sdeltay, nnz, smaxz, sdeltaz, l_2dxz)
  if zlb==periodic:
    b.sidezl=b.core
  if zlb==em3d.otherproc:
    allocatef(b.sidezl,stencil)
    b.sidezl.proc=proczl
  b.sidezr = EM3D_FIELDtype()
  if zrb==openbc:
    allocatesf(b.sidezr,stencil)
    init_splitfield(b.sidezr.syf,nx,ny,nbndz,nxguard,nyguard,nzguard, dt, dx, dy, dz, clight, 0, 0, 1, nnx, smaxx, sdeltax, nny, smaxy, sdeltay, nnz, smaxz, sdeltaz, l_2dxz)
  if zrb==periodic:
    b.sidezr=b.core
  if zrb==em3d.otherproc:
    allocatef(b.sidezr,stencil)
    b.sidezr.proc=proczr
    
# --- edges
# xy
  b.edgexlyl = EM3D_FIELDtype()
  if xlb==openbc and ylb==openbc:
    allocatesf(b.edgexlyl,stencil)
    init_splitfield(b.edgexlyl.syf,nbndx,nbndy,nz,nxguard,nyguard,nzguard, dt, dx, dy, dz, clight,-1,-1, 0, nnx, smaxx, sdeltax, nny, smaxy, sdeltay, nnz, smaxz, sdeltaz, l_2dxz)
  else:
    if xlb==openbc:
      if ylb==periodic:
        b.edgexlyl=b.sidexl
      if ylb==em3d.otherproc:
        allocatesf(b.edgexlyl,stencil)
        b.edgexlyl.proc=procyl
    if ylb==openbc:
      if xlb==periodic:
        b.edgexlyl=b.sideyl
      if xlb==em3d.otherproc:
        allocatesf(b.edgexlyl,stencil)
        b.edgexlyl.proc=procxl
  b.edgexryl = EM3D_FIELDtype()
  if xrb==openbc and ylb==openbc:
    allocatesf(b.edgexryl,stencil)
    init_splitfield(b.edgexryl.syf,nbndx,nbndy,nz,nxguard,nyguard,nzguard, dt, dx, dy, dz, clight, 1,-1, 0, nnx, smaxx, sdeltax, nny, smaxy, sdeltay, nnz, smaxz, sdeltaz, l_2dxz)
  else:
    if xrb==openbc:
      if ylb==periodic:
        b.edgexryl=b.sidexr
      if ylb==em3d.otherproc:
        allocatesf(b.edgexryl,stencil)
        b.edgexryl.proc=procyl
    if ylb==openbc:
      if xrb==periodic:
        b.edgexryl=b.sideyl
      if xrb==em3d.otherproc:
        allocatesf(b.edgexryl,stencil)
        b.edgexryl.proc=procxr
  b.edgexlyr = EM3D_FIELDtype()
  if xlb==openbc and yrb==openbc:
    allocatesf(b.edgexlyr,stencil)
    init_splitfield(b.edgexlyr.syf,nbndx,nbndy,nz,nxguard,nyguard,nzguard, dt, dx, dy, dz, clight,-1, 1, 0, nnx, smaxx, sdeltax, nny, smaxy, sdeltay, nnz, smaxz, sdeltaz, l_2dxz)
  else:
    if xlb==openbc:
      if yrb==periodic:
        b.edgexlyr=b.sidexl
      if yrb==em3d.otherproc:
        allocatesf(b.edgexlyr,stencil)
        b.edgexlyr.proc=procyr
    if yrb==openbc:
      if xlb==periodic:
        b.edgexlyr=b.sideyr
      if xlb==em3d.otherproc:
        allocatesf(b.edgexlyr,stencil)
        b.edgexlyr.proc=procxl
  b.edgexryr = EM3D_FIELDtype()
  if xrb==openbc and yrb==openbc:
    allocatesf(b.edgexryr,stencil)
    init_splitfield(b.edgexryr.syf,nbndx,nbndy,nz,nxguard,nyguard,nzguard, dt, dx, dy, dz, clight, 1, 1, 0, nnx, smaxx, sdeltax, nny, smaxy, sdeltay, nnz, smaxz, sdeltaz, l_2dxz)
  else:
    if xrb==openbc:
      if yrb==periodic:
        b.edgexryr=b.sidexr
      if yrb==em3d.otherproc:
        allocatesf(b.edgexryr,stencil)
        b.edgexryr.proc=procyr
    if yrb==openbc:
      if xrb==periodic:
        b.edgexryr=b.sideyr
      if xrb==em3d.otherproc:
        allocatesf(b.edgexryr,stencil)
        b.edgexryr.proc=procxr
# xz
  b.edgexlzl = EM3D_FIELDtype()
  if xlb==openbc and zlb==openbc:
    allocatesf(b.edgexlzl,stencil)
    init_splitfield(b.edgexlzl.syf,nbndx,ny,nbndz,nxguard,nyguard,nzguard, dt, dx, dy, dz, clight,-1, 0,-1, nnx, smaxx, sdeltax, nny, smaxy, sdeltay, nnz, smaxz, sdeltaz, l_2dxz)
  else:
    if xlb==openbc:
      if zlb==periodic:
        b.edgexlzl=b.sidexl
      if zlb==em3d.otherproc:
        allocatesf(b.edgexlzl,stencil)
        b.edgexlzl.proc=proczl
    if zlb==openbc:
      if xlb==periodic:
        b.edgexlzl=b.sidezl
      if xlb==em3d.otherproc:
        allocatesf(b.edgexlzl,stencil)
        b.edgexlzl.proc=procxl
  b.edgexrzl = EM3D_FIELDtype()
  if xrb==openbc and zlb==openbc:
    allocatesf(b.edgexrzl,stencil)
    init_splitfield(b.edgexrzl.syf,nbndx,ny,nbndz,nxguard,nyguard,nzguard, dt, dx, dy, dz, clight, 1, 0,-1, nnx, smaxx, sdeltax, nny, smaxy, sdeltay, nnz, smaxz, sdeltaz, l_2dxz)
  else:
    if xrb==openbc:
      if zlb==periodic:
        b.edgexrzl=b.sidexr
      if zlb==em3d.otherproc:
        allocatesf(b.edgexrzl,stencil)
        b.edgexrzl.proc=proczl
    if zlb==openbc:
      if xrb==periodic:
        b.edgexrzl=b.sidezl
      if xrb==em3d.otherproc:
        allocatesf(b.edgexrzl,stencil)
        b.edgexrzl.proc=procxr
  b.edgexlzr = EM3D_FIELDtype()
  if xlb==openbc and zrb==openbc:
    allocatesf(b.edgexlzr,stencil)
    init_splitfield(b.edgexlzr.syf,nbndx,ny,nbndz,nxguard,nyguard,nzguard, dt, dx, dy, dz, clight,-1, 0, 1, nnx, smaxx, sdeltax, nny, smaxy, sdeltay, nnz, smaxz, sdeltaz, l_2dxz)
  else:
    if xlb==openbc:
      if zrb==periodic:
        b.edgexlzr=b.sidexl
      if zrb==em3d.otherproc:
        allocatesf(b.edgexlzr,stencil)
        b.edgexlzr.proc=proczr
    if zrb==openbc:
      if xlb==periodic:
        b.edgexlzr=b.sidezr
      if xlb==em3d.otherproc:
        allocatesf(b.edgexlzr,stencil)
        b.edgexlzr.proc=procxl
  b.edgexrzr = EM3D_FIELDtype()
  if xrb==openbc and zrb==openbc:
    allocatesf(b.edgexrzr,stencil)
    init_splitfield(b.edgexrzr.syf,nbndx,ny,nbndz,nxguard,nyguard,nzguard, dt, dx, dy, dz, clight, 1, 0, 1, nnx, smaxx, sdeltax, nny, smaxy, sdeltay, nnz, smaxz, sdeltaz, l_2dxz)
  else:
    if xrb==openbc:
      if zrb==periodic:
        b.edgexrzr=b.sidexr
      if zrb==em3d.otherproc:
        allocatesf(b.edgexrzr,stencil)
        b.edgexrzr.proc=proczr
    if zrb==openbc:
      if xrb==periodic:
        b.edgexrzr=b.sidezr
      if xrb==em3d.otherproc:
        allocatesf(b.edgexrzr,stencil)
        b.edgexrzr.proc=procxr

# yz
  b.edgeylzl = EM3D_FIELDtype()
  if ylb==openbc and zlb==openbc:
    allocatesf(b.edgeylzl,stencil)
    init_splitfield(b.edgeylzl.syf,nx,nbndy,nbndz,nxguard,nyguard,nzguard, dt, dx, dy, dz, clight, 0,-1,-1, nnx, smaxx, sdeltax, nny, smaxy, sdeltay, nnz, smaxz, sdeltaz, l_2dxz)
  else:
    if ylb==openbc:
      if zlb==periodic:
        b.edgeylzl=b.sideyl
      if zlb==em3d.otherproc:
        allocatesf(b.edgeylzl,stencil)
        b.edgeylzl.proc=proczl
    if zlb==openbc:
      if ylb==periodic:
        b.edgeylzl=b.sidezl
      if ylb==em3d.otherproc:
        allocatesf(b.edgeylzl,stencil)
        b.edgeylzl.proc=procyl
  b.edgeyrzl = EM3D_FIELDtype()
  if yrb==openbc and zlb==openbc:
    allocatesf(b.edgeyrzl,stencil)
    init_splitfield(b.edgeyrzl.syf,nx,nbndy,nbndz,nxguard,nyguard,nzguard, dt, dx, dy, dz, clight, 0, 1,-1, nnx, smaxx, sdeltax, nny, smaxy, sdeltay, nnz, smaxz, sdeltaz, l_2dxz)
  else:
    if yrb==openbc:
      if zlb==periodic:
        b.edgeyrzl=b.sideyr
      if zlb==em3d.otherproc:
        allocatesf(b.edgeyrzl,stencil)
        b.edgeyrzl.proc=proczl
    if zlb==openbc:
      if yrb==periodic:
        b.edgeyrzl=b.sidezl
      if yrb==em3d.otherproc:
        allocatesf(b.edgeyrzl,stencil)
        b.edgeyrzl.proc=procyr
  b.edgeylzr = EM3D_FIELDtype()
  if ylb==openbc and zrb==openbc:
    allocatesf(b.edgeylzr,stencil)
    init_splitfield(b.edgeylzr.syf,nx,nbndy,nbndz,nxguard,nyguard,nzguard, dt, dx, dy, dz, clight, 0,-1, 1, nnx, smaxx, sdeltax, nny, smaxy, sdeltay, nnz, smaxz, sdeltaz, l_2dxz)
  else:
    if ylb==openbc:
      if zrb==periodic:
        b.edgeylzr=b.sideyl
      if zrb==em3d.otherproc:
        allocatesf(b.edgeylzr,stencil)
        b.edgeylzr.proc=proczr
    if zrb==openbc:
      if ylb==periodic:
        b.edgeylzr=b.sidezr
      if ylb==em3d.otherproc:
        allocatesf(b.edgeylzr,stencil)
        b.edgeylzr.proc=procyl
  b.edgeyrzr = EM3D_FIELDtype()
  if yrb==openbc and zrb==openbc:
    allocatesf(b.edgeyrzr,stencil)
    init_splitfield(b.edgeyrzr.syf,nx,nbndy,nbndz,nxguard,nyguard,nzguard, dt, dx, dy, dz, clight, 0, 1, 1, nnx, smaxx, sdeltax, nny, smaxy, sdeltay, nnz, smaxz, sdeltaz, l_2dxz)
  else:
    if yrb==openbc:
      if zrb==periodic:
        b.edgeyrzr=b.sideyr
      if zrb==em3d.otherproc:
        allocatesf(b.edgeyrzr,stencil)
        b.edgeyrzr.proc=proczr
    if zrb==openbc:
      if yrb==periodic:
        b.edgeyrzr=b.sidezr
      if yrb==em3d.otherproc:
        allocatesf(b.edgeyrzr,stencil)
        b.edgeyrzr.proc=procyr

# *** corners
  b.cornerxlylzl = EM3D_FIELDtype()
  if xlb==openbc and ylb==openbc and zlb==openbc:
    allocatesf(b.cornerxlylzl,stencil)
    init_splitfield(b.cornerxlylzl.syf,nbndx,nbndy,nbndz,nxguard,nyguard,nzguard, dt, dx, dy, dz, clight,-1,-1,-1, nnx, smaxx, sdeltax, nny, smaxy, sdeltay, nnz, smaxz, sdeltaz, l_2dxz)
  else:
    if xlb==openbc and ylb==openbc:
      if zlb==periodic:
        b.cornerxlylzl=b.edgexlyl
      if zlb==em3d.otherproc:
        allocatesf(b.cornerxlylzl,stencil)
        b.cornerxlylzl.proc=proczl      
    if xlb==openbc and zlb==openbc:
      if ylb==periodic:
        b.cornerxlylzl=b.edgexlzl
      if ylb==em3d.otherproc:
        allocatesf(b.cornerxlylzl,stencil)
        b.cornerxlylzl.proc=procyl      
    if ylb==openbc and zlb==openbc:
      if xlb==periodic:
        b.cornerxlylzl=b.edgeylzl
      if xlb==em3d.otherproc:
        allocatesf(b.cornerxlylzl,stencil)
        b.cornerxlylzl.proc=procxl      
  b.cornerxrylzl = EM3D_FIELDtype()
  if xrb==openbc and ylb==openbc and zlb==openbc:
    allocatesf(b.cornerxrylzl,stencil)
    init_splitfield(b.cornerxrylzl.syf,nbndx,nbndy,nbndz,nxguard,nyguard,nzguard, dt, dx, dy, dz, clight, 1,-1,-1, nnx, smaxx, sdeltax, nny, smaxy, sdeltay, nnz, smaxz, sdeltaz, l_2dxz)
  else:
    if xrb==openbc and ylb==openbc:
      if zlb==periodic:
        b.cornerxrylzl=b.edgexryl
      if zlb==em3d.otherproc:
        allocatesf(b.cornerxrylzl,stencil)
        b.cornerxrylzl.proc=proczl      
    if xrb==openbc and zlb==openbc:
      if ylb==periodic:
        b.cornerxrylzl=b.edgexrzl
      if ylb==em3d.otherproc:
        allocatesf(b.cornerxrylzl,stencil)
        b.cornerxrylzl.proc=procyl      
    if ylb==openbc and zlb==openbc:
      if xrb==periodic:
        b.cornerxrylzl=b.edgeylzl
      if xrb==em3d.otherproc:
        allocatesf(b.cornerxrylzl,stencil)
        b.cornerxrylzl.proc=procxr      
  b.cornerxlyrzl = EM3D_FIELDtype()
  if xlb==openbc and yrb==openbc and zlb==openbc:
    allocatesf(b.cornerxlyrzl,stencil)
    init_splitfield(b.cornerxlyrzl.syf,nbndx,nbndy,nbndz,nxguard,nyguard,nzguard, dt, dx, dy, dz, clight,-1, 1,-1, nnx, smaxx, sdeltax, nny, smaxy, sdeltay, nnz, smaxz, sdeltaz, l_2dxz)
  else:
    if xlb==openbc and yrb==openbc:
      if zlb==periodic:
        b.cornerxlyrzl=b.edgexlyr
      if zlb==em3d.otherproc:
        allocatesf(b.cornerxlyrzl,stencil)
        b.cornerxlyrzl.proc=proczl      
    if xlb==openbc and zlb==openbc:
      if yrb==periodic:
        b.cornerxlyrzl=b.edgexlzl
      if yrb==em3d.otherproc:
        allocatesf(b.cornerxlyrzl,stencil)
        b.cornerxlyrzl.proc=procyr      
    if yrb==openbc and zlb==openbc:
      if xlb==periodic:
        b.cornerxlyrzl=b.edgeyrzl
      if xlb==em3d.otherproc:
        allocatesf(b.cornerxlyrzl,stencil)
        b.cornerxlyrzl.proc=procxl      
  b.cornerxryrzl = EM3D_FIELDtype()
  if xrb==openbc and yrb==openbc and zlb==openbc:
    allocatesf(b.cornerxryrzl,stencil)
    init_splitfield(b.cornerxryrzl.syf,nbndx,nbndy,nbndz,nxguard,nyguard,nzguard, dt, dx, dy, dz, clight, 1, 1,-1, nnx, smaxx, sdeltax, nny, smaxy, sdeltay, nnz, smaxz, sdeltaz, l_2dxz)
  else:
    if xrb==openbc and yrb==openbc:
      if zlb==periodic:
        b.cornerxryrzl=b.edgexryr
      if zlb==em3d.otherproc:
        allocatesf(b.cornerxryrzl,stencil)
        b.cornerxryrzl.proc=proczl      
    if xrb==openbc and zlb==openbc:
      if yrb==periodic:
        b.cornerxryrzl=b.edgexrzl
      if yrb==em3d.otherproc:
        allocatesf(b.cornerxryrzl,stencil)
        b.cornerxryrzl.proc=procyr      
    if yrb==openbc and zlb==openbc:
      if xrb==periodic:
        b.cornerxryrzl=b.edgeyrzl
      if xrb==em3d.otherproc:
        allocatesf(b.cornerxryrzl,stencil)
        b.cornerxryrzl.proc=procxr      
  b.cornerxlylzr = EM3D_FIELDtype()
  if xlb==openbc and ylb==openbc and zrb==openbc:
    allocatesf(b.cornerxlylzr,stencil)
    init_splitfield(b.cornerxlylzr.syf,nbndx,nbndy,nbndz,nxguard,nyguard,nzguard, dt, dx, dy, dz, clight,-1,-1, 1, nnx, smaxx, sdeltax, nny, smaxy, sdeltay, nnz, smaxz, sdeltaz, l_2dxz)
  else:
    if xlb==openbc and ylb==openbc:
      if zrb==periodic:
        b.cornerxlylzr=b.edgexlyl
      if zrb==em3d.otherproc:
        allocatesf(b.cornerxlylzr,stencil)
        b.cornerxlylzr.proc=proczr      
    if xlb==openbc and zrb==openbc:
      if ylb==periodic:
        b.cornerxlylzr=b.edgexlzr
      if ylb==em3d.otherproc:
        allocatesf(b.cornerxlylzr,stencil)
        b.cornerxlylzr.proc=procyl      
    if ylb==openbc and zrb==openbc:
      if xlb==periodic:
        b.cornerxlylzr=b.edgeylzr
      if xlb==em3d.otherproc:
        allocatesf(b.cornerxlylzr,stencil)
        b.cornerxlylzr.proc=procxl      
  b.cornerxrylzr = EM3D_FIELDtype()
  if xrb==openbc and ylb==openbc and zrb==openbc:
    allocatesf(b.cornerxrylzr,stencil)
    init_splitfield(b.cornerxrylzr.syf,nbndx,nbndy,nbndz,nxguard,nyguard,nzguard, dt, dx, dy, dz, clight, 1,-1, 1, nnx, smaxx, sdeltax, nny, smaxy, sdeltay, nnz, smaxz, sdeltaz, l_2dxz)
  else:
    if xrb==openbc and ylb==openbc:
      if zrb==periodic:
        b.cornerxrylzr=b.edgexryl
      if zrb==em3d.otherproc:
        allocatesf(b.cornerxrylzr,stencil)
        b.cornerxrylzr.proc=proczr      
    if xrb==openbc and zrb==openbc:
      if ylb==periodic:
        b.cornerxrylzr=b.edgexrzr
      if ylb==em3d.otherproc:
        allocatesf(b.cornerxrylzr,stencil)
        b.cornerxrylzr.proc=procyl      
    if ylb==openbc and zrb==openbc:
      if xrb==periodic:
        b.cornerxrylzr=b.edgeylzr
      if xrb==em3d.otherproc:
        allocatesf(b.cornerxrylzr,stencil)
        b.cornerxrylzr.proc=procxr      
  b.cornerxlyrzr = EM3D_FIELDtype()
  if xlb==openbc and yrb==openbc and zrb==openbc:
    allocatesf(b.cornerxlyrzr,stencil)
    init_splitfield(b.cornerxlyrzr.syf,nbndx,nbndy,nbndz,nxguard,nyguard,nzguard, dt, dx, dy, dz, clight,-1, 1, 1, nnx, smaxx, sdeltax, nny, smaxy, sdeltay, nnz, smaxz, sdeltaz, l_2dxz)
  else:
    if xlb==openbc and yrb==openbc:
      if zrb==periodic:
        b.cornerxlyrzr=b.edgexlyr
      if zrb==em3d.otherproc:
        allocatesf(b.cornerxlyrzr,stencil)
        b.cornerxlyrzr.proc=proczr      
    if xlb==openbc and zrb==openbc:
      if yrb==periodic:
        b.cornerxlyrzr=b.edgexlzr
      if yrb==em3d.otherproc:
        allocatesf(b.cornerxlyrzr,stencil)
        b.cornerxlyrzr.proc=procyr      
    if yrb==openbc and zrb==openbc:
      if xlb==periodic:
        b.cornerxlyrzr=b.edgeyrzr
      if xlb==em3d.otherproc:
        allocatesf(b.cornerxlyrzr,stencil)
        b.cornerxlyrzr.proc=procxl      
  b.cornerxryrzr = EM3D_FIELDtype()
  if xrb==openbc and yrb==openbc and zrb==openbc:
    allocatesf(b.cornerxryrzr,stencil)
    init_splitfield(b.cornerxryrzr.syf,nbndx,nbndy,nbndz,nxguard,nyguard,nzguard, dt, dx, dy, dz, clight, 1, 1, 1, nnx, smaxx, sdeltax, nny, smaxy, sdeltay, nnz, smaxz, sdeltaz, l_2dxz)
  else:
    if xrb==openbc and yrb==openbc:
      if zrb==periodic:
        b.cornerxryrzr=b.edgexryr
      if zrb==em3d.otherproc:
        allocatesf(b.cornerxryrzr,stencil)
        b.cornerxryrzr.proc=proczr      
    if xrb==openbc and zrb==openbc:
      if yrb==periodic:
        b.cornerxryrzr=b.edgexrzr
      if yrb==em3d.otherproc:
        allocatesf(b.cornerxryrzr,stencil)
        b.cornerxryrzr.proc=procyr      
    if yrb==openbc and zrb==openbc:
      if xrb==periodic:
        b.cornerxryrzr=b.edgeyrzr
      if xrb==em3d.otherproc:
        allocatesf(b.cornerxryrzr,stencil)
        b.cornerxryrzr.proc=procxr      

  return b

##############################################################################
# --- This can only be done after the class is defined.
#try:
#  psyco.bind(EM3D)
#except NameError:
#  pass

