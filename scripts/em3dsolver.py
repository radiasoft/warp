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
  __flaginputs__ = {'mode':1,'l_apply_pml':true,
                    'nbndx':10,'nbndy':10,'nbndz':10,
                    'nxguard':1,'nyguard':1,'nzguard':1,
                    'l_particles_weight':false,'l_usecoeffs':false,
                    'l_pushf':false,'l_pushpot':false,'l_verbose':false,
                    'laser_amplitude':1.,'laser_profile':None,
                    'laser_gauss_widthx':None,'laser_gauss_centerx':0.,
                    'laser_gauss_widthy':None,'laser_gauss_centery':0.,
                    'laser_anglex':0.,'laser_angley':0.,
                    'laser_polangle':0.,
                    'laser_wavelength':None,'laser_wavenumber':None,
                    'laser_frequency':None,
                    'laser_source_z':None,'laser_source_v':0.,
                    'laser_focus_z':None,'laser_focus_v':0.,
                    'ncyclesperstep':1,'ncyclesperstep':None,
                    'l_2dxz':0,'l_1dz':0,'l_sumjx':0,
                    'npass_smooth':array([[0],[0],[0]]),
                    'alpha_smooth':array([[0.5],[0.5],[0.5]]),
                    'stride_smooth':array([[1],[1],[1]]),
                    'l_smooth_particle_fields':True,
                    'n_smooth_fields':None,
                    'autoset_timestep':true,'dtcoef':1.,#0.99,
                    'deposit_energy_density':false,'refinement':None,
                    'l_force_nzlocal2nz':false,'isactiveem':None,
                    'l_coarse_patch':false,'stencil':false,
                    'l_esirkepov':true,'theta_damp':0.,
                    'sigmae':0.,'sigmab':0.,
                    'colecoefs':None}

  def __init__(self,**kw):
    self.solveroff=False # flag to turn off the solver, for testing purpose
    assert top.grid_overlap<2,"The EM solver needs top.grid_overlap<2!"
    self.grid_overlap = 1
    FieldSolver.__init__(self,kwdict=kw)
#    top.allspecl = true
    top.lcallfetchb = true
    top.lgridqnt = true
    self.xgrid=0.
    self.xgridcont=0.
    self.nxshifts=0
    self.vxgrid=0.
    self.zgrid=top.zgrid
    self.nzshifts=0
    self.odd=0
    # --- Save input parameters
#    self.processdefaultsfrompackage(EM3D.__w3dinputs__,w3d,kw)
#    self.processdefaultsfrompackage(EM3D.__topinputs__,top,kw)
    self.processdefaultsfrompackage(EM3D.__em3dinputs__,em3d,kw)
    self.processdefaultsfromdict(EM3D.__flaginputs__,kw)
    if self.isactiveem is not None:
      self.isactive=self.isactiveem
    else:
      try:
        self.isactive=self.isactive
      except:
        self.isactive=True
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
    self.stride_smooth  = array(self.stride_smooth)
#    minguards = array([1+int(top.depos_order.max(1)/2),self.npass_smooth.sum(1)]).max(0)
    minguards = 2+int(top.depos_order.max(1)/2)+(self.npass_smooth*self.stride_smooth).sum(1)
    if self.nxguard==1:self.nxguard = minguards[0]
    if self.nyguard==1:self.nyguard = minguards[1]
    if self.nzguard==1:self.nzguard = minguards[2]
    if self.stencil>0:
      if self.nxguard<2:self.nxguard=2
      if self.nyguard<2:self.nyguard=2
      if self.nzguard<2:self.nzguard=2
#      self.npass_smooth = where(self.npass_smooth==0,1,self.npass_smooth)
    if self.l_2dxz:
      self.nyguard=self.nylocal=self.ny=0

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
    try:
      if self.kw.has_key('bounds'):self.kw.pop('bounds')
    except:
      pass
      
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

    # --- sets coefficients of Cole solver
    if self.colecoefs is None:
      if self.l_2dxz:
        em3d.betax = 1./8.
        em3d.alphax = 1.-2*em3d.betax
      else:
        em3d.alphax = 7./12.
        em3d.betax  = 1./12.
        em3d.gammax = 1./48.
    else:
      em3d.alphax = self.colecoefs[0]
      em3d.betax  = self.colecoefs[1]
      if not self.l_2dxz:
        em3d.gammax = self.colecoefs[2]
    em3d.alphay = em3d.alphaz = em3d.alphax
    em3d.betay  = em3d.betaz  = em3d.betax
    em3d.gammay = em3d.gammaz = em3d.gammax

    # --- set time step as a fraction of Courant condition
    # --- also set self.ncyclesperstep if top.dt over Courant condition times dtcoef
    try:
      parentid = self.parents[0]
      parent = self.root.listofblocks[parentid]
      try:
        sibling = parent.children[0]
      except:
        sibling = None
    except:
      parent = None
      sibling = None
    if self.ncyclesperstep is None:
      self.ncyclesperstep = 1
      if sibling is None:
        if self.autoset_timestep:
          if self.mode==2 and self.dtcoef>0.5:self.dtcoef/=2
          if self.l_2dxz:
            if self.stencil==0:
              dtcourant=1./(clight*sqrt(1./self.dx**2+1./self.dz**2))
            else:
              dtcourant=min(self.dx,self.dz)/clight 
              Cx = em3d.alphax -2.*em3d.betax
              Cz = em3d.alphaz -2.*em3d.betaz
              dtcourant=1./(clight*sqrt(Cx/self.dx**2+Cz/self.dz**2))
          else:
            if self.stencil==0:
              dtcourant=1./(clight*sqrt(1./self.dx**2+1./self.dy**2+1./self.dz**2))
            else:
              dtcourant=min(self.dx,self.dy,self.dz)/clight 
              Cx = em3d.alphax -4.*em3d.betax+4.*em3d.gammax
              Cy = em3d.alphay -4.*em3d.betay+4.*em3d.gammay
              Cz = em3d.alphaz -4.*em3d.betaz+4.*em3d.gammaz
              dtcourant=1./(clight*sqrt(Cx/self.dx**2+Cy/self.dy**2+Cz/self.dz**2))
          if self.theta_damp>0.:
            dtcourant*=sqrt((2.+self.theta_damp)/(2.+3.*self.theta_damp))
          if top.dt==0.:
            top.dt=dtcourant*self.dtcoef
          else:
            if top.dt>(self.dtcoef*dtcourant):
#              self.ncyclesperstep = (nint(top.dt/(self.dtcoef*dtcourant))+0)
              self.ncyclesperstep = int(top.dt/(self.dtcoef*dtcourant)+1.)
              print '#1', self.ncyclesperstep,top.dt,self.dtcoef,dtcourant
            else:
              self.ncyclesperstep = 1./(nint((self.dtcoef*dtcourant)/top.dt)+0)
              print '#2', self.ncyclesperstep,top.dt,self.dtcoef,dtcourant
      else:
        self.ncyclesperstep=sibling.ncyclesperstep
    self.dtinit = top.dt
    
    if  top.vbeamfrm<>0.:self.bounds[-2:]=-1
    
    self.setbcparallel(0) # x
    self.setbcparallel(1) # y
    self.setbcparallel(2) # z
    
    if self.refinement is not None:
      ref=self.refinement
      self.field_coarse = self.__class__(l_force_nzlocal2nz=True,
                               l_coarse_patch=True,
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
                               bounds=self.bounds,
                               isactiveem=self.isactive,
                               ncyclesperstep=self.root.listofblocks[self.parents[0]].ncyclesperstep,
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
    if self.refinement is not None: # --- disable laser on MR patches
      self.laser_profile=None
      self.field_coarse.laser_profile=None
      
    # ---- install 2nd part of field solve (only for root)
    if self.refinement is None and not self.l_coarse_patch:
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
                                   self.ncyclesperstep,
                                   self.l_2dxz,
                                   self.theta_damp,
                                   self.sigmae,
                                   self.sigmab)
    self.fields = self.block.core.yf
    
    
  def setuplaser(self):
    if self.laser_profile is not None:
      if self.laser_frequency is None:
        if self.laser_wavenumber is not None:
          self.laser_frequency = clight*self.laser_wavenumber
        elif self.laser_wavelength is not None:
          self.laser_frequency = 2.*pi*clight/self.laser_wavelength
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
      xxex -= self.laser_gauss_centerx; xxex /= self.laser_gauss_widthx
      xxey -= self.laser_gauss_centerx; xxey /= self.laser_gauss_widthx
      yyex -= self.laser_gauss_centery; yyex /= self.laser_gauss_widthy
      yyey -= self.laser_gauss_centery; yyey /= self.laser_gauss_widthy
      self.laser_profile = [exp(-(xxex**2+yyex**2)/2.),
                            exp(-(xxey**2+yyey**2)/2.)]
    elif operator.isSequenceType(self.laser_profile):
      print f.nx,f.ny,shape(self.laser_profile)
      assert len(self.laser_profile[:,0]) == f.nx+1,"The specified profile must be of length nx+1"
      assert len(self.laser_profile[0,:]) == f.ny+1,"The specified profile must be of length ny+1"
      self.laser_profile_init = self.laser_profile.copy()
      self.laser_profile = [0.5*(self.laser_profile_init[1:,:]+self.laser_profile_init[:-1,:]),
                            0.5*(self.laser_profile_init[:,1:]+self.laser_profile_init[:,:-1])]
    elif callable(self.laser_profile):
      self.laser_profile_func = self.laser_profile

  def add_laser(self,field):
    if self.laser_profile is None: 
      self.block.core.yf.E_inz_pos=w3d.zmmin-(self.nzguard*2.)*self.dz
      return

    if 1:#self.laser_source_z>self.zmmin+self.zgrid and self.laser_source_z<=self.zmmax+self.zgrid:
      self.block.core.yf.E_inz_pos = self.laser_source_z-self.zgrid
      if self.laser_focus_z is not None:self.laser_focus_z+=self.laser_focus_v*top.dt/self.ncyclesperstep
      self.laser_source_z+=self.laser_source_v*top.dt/self.ncyclesperstep
    else:
      return

    if self.laser_amplitude_func is not None:
      self.laser_amplitude = self.laser_amplitude_func(top.time*(1.-self.laser_source_v/clight))
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
      if self.laser_focus_z is not None:
        z0 = self.laser_focus_z
        if self.laser_focus_z>0.:
          phaseex = (-(sqrt(xxex**2+yyex**2+z0**2)-z0)/clight-top.time*(1.-self.laser_source_v/clight))*self.laser_frequency
          phaseey = (-(sqrt(xxey**2+yyey**2+z0**2)-z0)/clight-top.time*(1.-self.laser_source_v/clight))*self.laser_frequency
        else:
          phaseex = ((sqrt(xxex**2+yyex**2+z0**2)-z0)/clight-top.time*(1.-self.laser_source_v/clight))*self.laser_frequency
          phaseey = ((sqrt(xxey**2+yyey**2+z0**2)-z0)/clight-top.time*(1.-self.laser_source_v/clight))*self.laser_frequency
      else:
        phaseex = ((xxex*sin(self.laser_anglex)+yyex*sin(self.laser_angley))/clight-top.time*(1.-self.laser_source_v/clight))*self.laser_frequency
        phaseey = ((xxey*sin(self.laser_anglex)+yyey*sin(self.laser_angley))/clight-top.time*(1.-self.laser_source_v/clight))*self.laser_frequency
    else:
      phase = 0.
    f = self.block.core.yf
    if self.l_2dxz:
      laser_amplitude=self.laser_amplitude*top.dt*clight/w3d.dz
      f.Ex_inz[f.jxmin:f.jxmax  ,f.jymin] = laser_amplitude*self.laser_profile[0]*sin(phaseex)*cos(self.laser_polangle)*(1.-self.laser_source_v/clight)
      f.Ey_inz[f.jxmin:f.jxmax+1,f.jymin] = laser_amplitude*self.laser_profile[1]*sin(phaseey)*sin(self.laser_polangle)*(1.-self.laser_source_v/clight)
    else:      
      laser_amplitude=self.laser_amplitude*top.dt*clight/w3d.dz
      f.Ex_inz[f.jxmin:f.jxmax  ,f.jymin:f.jymax+1] = laser_amplitude*self.laser_profile[0]*sin(phaseex)*cos(self.laser_polangle)*(1.-self.laser_source_v/clight)
      f.Ey_inz[f.jxmin:f.jxmax+1,f.jymin:f.jymax  ] = laser_amplitude*self.laser_profile[1]*sin(phaseey)*sin(self.laser_polangle)*(1.-self.laser_source_v/clight)
    
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
       if nox==1 and noy==1 and noz==1 and not w3d.l4symtry:
        getf3d_linear(n,x,y,z,ex,ey,ez,
                 f.xmin,f.ymin,f.zmin+self.zgrid,
                 f.dx,f.dy,f.dz,
                 f.nx,f.ny,f.nz,
                 f.nxguard,f.nyguard,f.nzguard,
                 f.Exp,f.Eyp,f.Ezp)
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
      if self.l_2dxz:
        gete2dxz_n_energy_conserving(n,x,z,ex,ey,ez,
                      f.xmin,f.zmin+self.zgrid,
                      f.dx,f.dz,
                      f.nx,f.nz,
                      f.nxguard,f.nzguard,
                      nox,noz,
                      f.Exp,f.Eyp,f.Ezp,
                      w3d.l4symtry)
      else:
       if nox==1 and noy==1 and noz==1 and not w3d.l4symtry:
        gete3d_linear_energy_conserving(n,x,y,z,ex,ey,ez,
                      f.xmin,f.ymin,f.zmin+self.zgrid,
                      f.dx,f.dy,f.dz,
                      f.nx,f.ny,f.nz,
                      f.nxguard,f.nyguard,f.nzguard,
                      f.Exp,f.Eyp,f.Ezp)
       else:
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
      if self.l_2dxz:
        getb2dxz_n_energy_conserving(n,x,z,bx,by,bz,
                    f.xmin,f.zmin+self.zgrid,
                    f.dx,f.dz,
                    f.nx,f.nz,
                    f.nxguard,f.nzguard,
                    nox,noz,
                    f.Bxp,f.Byp,f.Bzp,
                    w3d.l4symtry)
      else:
       if nox==1 and noy==1 and noz==1 and not w3d.l4symtry:
        getb3d_linear_energy_conserving(n,x,y,z,bx,by,bz,
                    f.xmin,f.ymin,f.zmin+self.zgrid,
                    f.dx,f.dy,f.dz,
                    f.nx,f.ny,f.nz,
                    f.nxguard,f.nyguard,f.nzguard,
                    f.Bxp,f.Byp,f.Bzp)
       else:
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

  def getfieldsfrompositions(self,x,y,z):
    # --- This returns e and b at positions x,y,z
    n = len(x)
    if n == 0: return None
    f = self.block.core.yf
    if self.l_2dxz:
      ilist = compress((x>=f.xmin) & (x<f.xmax) &
                       (z>=(f.zmin+self.zgrid)) & (z<(f.zmax+self.zgrid)),arange(n))
    else:
      ilist = compress((x>=f.xmin) & (x<f.xmax) &
                       (y>=f.ymin) & (y<f.ymax) &
                       (z>=(f.zmin+self.zgrid)) & (z<(f.zmax+self.zgrid)),arange(n))
    nlocal = len(ilist)
    if nlocal>0:
      x = take(x,ilist)
      y = take(y,ilist)
      z = take(z,ilist)
      ex=zeros(nlocal,'d')
      ey=zeros(nlocal,'d')
      ez=zeros(nlocal,'d')
      bx=zeros(nlocal,'d')
      by=zeros(nlocal,'d')
      bz=zeros(nlocal,'d')
      self.fetchfieldfrompositions(x,y,z,ex,ey,ez,bx,by,bz)

    exp=zeros(n,'d')
    eyp=zeros(n,'d')
    ezp=zeros(n,'d')
    bxp=zeros(n,'d')
    byp=zeros(n,'d')
    bzp=zeros(n,'d')
    if me>0:
      mpi.send(nlocal,0,3)
      if nlocal>0:
        mpi.send((ilist,ex,ey,ez,bx,by,bz),0,3)
    else:
      if nlocal>0:
        for j,i in enumerate(ilist):
          exp[i] = ex[j]
          eyp[i] = ey[j]
          ezp[i] = ez[j]
          bxp[i] = bx[j]
          byp[i] = by[j]
          bzp[i] = bz[j]
      for ip in range(1,npes):
        nlocal = mpirecv(ip,3)
        if nlocal>0:
          ilist,ex,ey,ez,bx,by,bz=mpirecv(ip,3)
          for j,i in enumerate(ilist):
            exp[i] = ex[j]
            eyp[i] = ey[j]
            ezp[i] = ez[j]
            bxp[i] = bx[j]
            byp[i] = by[j]
            bzp[i] = bz[j]

    return exp,eyp,ezp,bxp,byp,bzp
  
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
      wfact = pgroup.pid[i:i+n,top.wpid-1]
    w3d.jsfsapi=js
    self.setsourcepatposition(x,y,z,ux,uy,uz,gaminv,wfact,zgrid,q,w)
    
  def setsourcepatposition(self,x,y,z,ux,uy,uz,gaminv,wfact,zgrid,q,w):
    n = x.shape[0]
    if n == 0: return
    # --- call routine performing current deposition
    f = self.block.core.yf
    js = w3d.jsfsapi
    nox = top.depos_order[0,js]
    noy = top.depos_order[1,js]
    noz = top.depos_order[2,js]
    if top.wpid==0:
      wfact = ones((1,),'d')
      l_particles_weight = false
    else:
      l_particles_weight = true
    if not self.deposit_energy_density:
      if self.l_1dz:
          jx = self.fields.J[self.fields.nxguard,self.fields.nyguard,:,0]*0.
          jy = self.fields.J[self.fields.nxguard,self.fields.nyguard,:,0]*0.
          jz = self.fields.J[self.fields.nxguard,self.fields.nyguard,:,0]*0.
          depose_jxjyjz_esirkepov_linear_serial1d(jx,jy,jz,n,x*0.,y,z,ux,uy,uz,gaminv,
                                            wfact*w,q/w3d.dx,
                                            -0.5,-0.5,f.zmin+self.zgrid,
                                            top.dt*top.pgroup.ndts[js],
                                            1.,1.,self.fields.dz,
                                            self.fields.nz,
                                            l_particles_weight)
#          for j in range(1,shape(self.fields.Jarray)[0]-1):
          for j in range(shape(self.fields.Jarray)[0]):
            self.fields.Jarray[j,self.fields.nyguard,:,0,0]+=jx
            self.fields.Jarray[j,self.fields.nyguard,:,1,0]+=jy
            self.fields.Jarray[j,self.fields.nyguard,:,2,0]+=jz
      elif self.l_2dxz:
        j = self.fields.J[:,self.fields.nyguard,:,:]*0.
        if 0:
#          depose_jxjy_esirkepov_linear_serial_2d(j,n,z,x,z-gaminv*uz*top.dt,x-gaminv*ux*top.dt,
#                                                 uy,gaminv,
#                                            wfact,q*w,
#                                            f.zmin+self.zgrid,f.xmin,
#                                            top.dt*top.pgroup.ndts[js],
#                                            f.dz,f.dx,
#                                            f.nz,f.nx,
#                                           l_particles_weight)
          depose_jxjy_esirkepov_linear_serial_2d(j,n,x,z,x-gaminv*ux*top.dt,z-gaminv*uz*top.dt,
                                                 uy,gaminv,
                                            wfact,q*w,
                                            f.xmin,f.zmin+self.zgrid,
                                            top.dt*top.pgroup.ndts[js],
                                            f.dx,f.dz,
                                            f.nx,f.nz,
                                            l_particles_weight)
        else:
         if self.l_esirkepov:
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
         else:
          depose_j_n_2dxz(j,n,
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
        self.fields.Jarray[:,self.fields.nyguard,:,:,0]+=j[:,:,:]
      else:
       if 0:#nox==1 and noy==1 and noz==1 and not w3d.l4symtry:
          depose_jxjyjz_esirkepov_linear_serial(self.fields.J,n,
                                            x,y,z,ux,uy,uz,
                                            gaminv,wfact,q*w,
                                            f.xmin,f.ymin,f.zmin+self.zgrid,
                                            top.dt*top.pgroup.ndts[js],
                                            f.dx,f.dy,f.dz,
                                            f.nx,f.ny,f.nz,
                                            f.nxguard,f.nyguard,f.nzguard,
                                            l_particles_weight)
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
      if self.l_2dxz:
        depose_rho_n_2dxz(self.fields.Rho,n,
                               x,z,
                               wfact,q*w,
                               f.xmin,f.zmin+self.zgrid,
                               f.dx,f.dz,
                               f.nx,f.nz,
                               f.nxguard,f.nzguard,
                               top.depos_order[0,js],
                               top.depos_order[2,js],
                               l_particles_weight,w3d.l4symtry)
      else:
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
#        if self.refinement is not None:
#          self.field_coarse.fields.Jarray[...,indts] = 0.
        if self.l_pushf:
          self.fields.Rhoarray[...,indts] = 0.
#          if self.refinement is not None:
#            self.field_coarse.fields.Rho[...,indts] = 0.
        if self.deposit_energy_density:
          self.fields.Mp[...] = 0.
  
  def setsourcepforparticles(self,isourcepndtscopies,indts,iselfb):
    if self.l_verbose:print 'setsourcepforparticles'
    # --- point J array to proper Jarray slice
    self.fields.J = self.fields.Jarray[:,:,:,:,indts]
    if self.l_pushf: self.fields.Rho = self.fields.Rhoarray[:,:,:,indts]

  def add_source_ndts_slices(self):  
    # --- add slices
    if top.nsndts>1:
      if self.refinement is not None:raise('Error in finalizesourcep:nsndts>1 not fully implemented yet with MR')
      for indts in range(top.nsndts-2,-1,-1):
        if top.ldts[indts]:
          add_current_slice_3d(self.fields,indts+1)
          if self.l_pushf:add_rho_slice_3d(self.fields,indts+1)

  def finalizesourcep(self):
    if self.sourcepfinalized: return
    self.sourcepfinalized = 1
    if self.l_verbose:print 'finalizesourcep'
    # --- add slices
    self.add_source_ndts_slices()
    self.aftersetsourcep()
    # --- smooth current density 
    if any(self.npass_smooth>0):self.smoothdensity()
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

    if self.l_sumjx:
      j = self.fields.Jarray[0,:,:,:,0]*0.
      for i in range(shape(self.fields.Jarray)[0]):
        j+=self.fields.Jarray[i,:,:,:,0]
      for i in range(shape(self.fields.Jarray)[0]):
        self.fields.Jarray[i,:,:,:,0]=j.copy()
        
    # --- point J to first slice of Jarray
    self.fields.J = self.fields.Jarray[:,:,:,:,0]
    if self.l_pushf:self.fields.Rho = self.fields.Rhoarray[:,:,:,0]
    
  def smoothdensity(self):
    if all(self.npass_smooth==0):return
    nx,ny,nz = shape(self.fields.J[...,0])
    nsm = shape(self.npass_smooth)[1]
    for js in range(nsm):
#      smooth3d_121(self.fields.J[...,0],nx-1,ny-1,nz-1,self.npass_smooth[:,js],self.alpha_smooth[:,js])
#      smooth3d_121(self.fields.J[...,1],nx-1,ny-1,nz-1,self.npass_smooth[:,js],self.alpha_smooth[:,js])
#      smooth3d_121(self.fields.J[...,2],nx-1,ny-1,nz-1,self.npass_smooth[:,js],self.alpha_smooth[:,js])
      smooth3d_121_stride(self.fields.J[...,0],nx-1,ny-1,nz-1,self.npass_smooth[:,js],self.alpha_smooth[:,js],self.stride_smooth[:,js])
      smooth3d_121_stride(self.fields.J[...,1],nx-1,ny-1,nz-1,self.npass_smooth[:,js],self.alpha_smooth[:,js],self.stride_smooth[:,js])
      smooth3d_121_stride(self.fields.J[...,2],nx-1,ny-1,nz-1,self.npass_smooth[:,js],self.alpha_smooth[:,js],self.stride_smooth[:,js])
      if self.l_pushf:
        smooth3d_121_stride(self.fields.Rho[...],nx-1,ny-1,nz-1,self.npass_smooth[:,js],self.alpha_smooth[:,js],self.stride_smooth[:,js])

  def smoothfields(self):
    nx,ny,nz = shape(self.fields.J[...,0])
    nsm = shape(self.npass_smooth)[1]
    for js in range(nsm):
      smooth3d_121_stride(self.fields.Exp,nx-1,ny-1,nz-1,self.npass_smooth[:,js],self.alpha_smooth[:,js],self.stride_smooth[:,js])
      smooth3d_121_stride(self.fields.Eyp,nx-1,ny-1,nz-1,self.npass_smooth[:,js],self.alpha_smooth[:,js],self.stride_smooth[:,js])
      smooth3d_121_stride(self.fields.Ezp,nx-1,ny-1,nz-1,self.npass_smooth[:,js],self.alpha_smooth[:,js],self.stride_smooth[:,js])
      smooth3d_121_stride(self.fields.Bxp,nx-1,ny-1,nz-1,self.npass_smooth[:,js],self.alpha_smooth[:,js],self.stride_smooth[:,js])
      smooth3d_121_stride(self.fields.Byp,nx-1,ny-1,nz-1,self.npass_smooth[:,js],self.alpha_smooth[:,js],self.stride_smooth[:,js])
      smooth3d_121_stride(self.fields.Bzp,nx-1,ny-1,nz-1,self.npass_smooth[:,js],self.alpha_smooth[:,js],self.stride_smooth[:,js])
#      smooth3d_121(self.fields.Exp,nx-1,ny-1,nz-1,self.npass_smooth[:,js],self.alpha_smooth[:,js])
#      smooth3d_121(self.fields.Eyp,nx-1,ny-1,nz-1,self.npass_smooth[:,js],self.alpha_smooth[:,js])
#      smooth3d_121(self.fields.Ezp,nx-1,ny-1,nz-1,self.npass_smooth[:,js],self.alpha_smooth[:,js])
#      smooth3d_121(self.fields.Bxp,nx-1,ny-1,nz-1,self.npass_smooth[:,js],self.alpha_smooth[:,js])
#      smooth3d_121(self.fields.Byp,nx-1,ny-1,nz-1,self.npass_smooth[:,js],self.alpha_smooth[:,js])
#      smooth3d_121(self.fields.Bzp,nx-1,ny-1,nz-1,self.npass_smooth[:,js],self.alpha_smooth[:,js])
    if self.n_smooth_fields is not None:
     if top.it%self.n_smooth_fields==0:
#      m = array([0,0,1]) # smooth only in z
#      f = 1.#array([0,0,1./self.n_smooth_fields])
#      npass_smooth = array([[ 0 , 0 ],[ 0 , 0 ],[ 1 , 1 ]])
#      alpha_smooth = array([[ 1., 1.],[ 1., 1.],[0.5, 3./2.]])
#      npass_smooth = array([[ 0 , 0 ],[ 0 , 0 ],[ 4 , 1 ]])
#      alpha_smooth = array([[ 1., 1.],[ 1., 1.],[0.5, 3.]])
#      for i in range(1):
#       for js in range(2):
      nsm = shape(self.npass_smooth)[1]
      for js in range(nsm):
        smooth3d_121_stride(self.fields.Ex,nx-1,ny-1,nz-1,self.npass_smooth[:,js],self.alpha_smooth[:,js],self.stride_smooth[:,js])
        smooth3d_121_stride(self.fields.Ey,nx-1,ny-1,nz-1,self.npass_smooth[:,js],self.alpha_smooth[:,js],self.stride_smooth[:,js])
        smooth3d_121_stride(self.fields.Ez,nx-1,ny-1,nz-1,self.npass_smooth[:,js],self.alpha_smooth[:,js],self.stride_smooth[:,js])
        smooth3d_121_stride(self.fields.Bx,nx-1,ny-1,nz-1,self.npass_smooth[:,js],self.alpha_smooth[:,js],self.stride_smooth[:,js])
        smooth3d_121_stride(self.fields.By,nx-1,ny-1,nz-1,self.npass_smooth[:,js],self.alpha_smooth[:,js],self.stride_smooth[:,js])
        smooth3d_121_stride(self.fields.Bz,nx-1,ny-1,nz-1,self.npass_smooth[:,js],self.alpha_smooth[:,js],self.stride_smooth[:,js])
#        smooth3d_121(self.fields.Ex,nx-1,ny-1,nz-1,npass_smooth[:,js]*m,alpha_smooth[:,js]*f)
#        smooth3d_121(self.fields.Ey,nx-1,ny-1,nz-1,npass_smooth[:,js]*m,alpha_smooth[:,js]*f)
#        smooth3d_121(self.fields.Ez,nx-1,ny-1,nz-1,npass_smooth[:,js]*m,alpha_smooth[:,js]*f)
#        smooth3d_121(self.fields.Bx,nx-1,ny-1,nz-1,npass_smooth[:,js]*m,alpha_smooth[:,js]*f)
#        smooth3d_121(self.fields.By,nx-1,ny-1,nz-1,npass_smooth[:,js]*m,alpha_smooth[:,js]*f)
#        smooth3d_121(self.fields.Bz,nx-1,ny-1,nz-1,npass_smooth[:,js]*m,alpha_smooth[:,js]*f)
        if self.l_pushf:
         smooth3d_121(self.fields.F,nx-1,ny-1,nz-1,npass_smooth[:,js]*m,alpha_smooth[:,js])
       
  def getsmoothx(self):
    nx = w3d.nx
    ny = 1
    nz = 1
    a = zeros([nx,ny,nz],'d')
    a[nx/2,:,:]=1.
    nsm = shape(self.npass_smooth)[1]
    for js in range(nsm):
      smooth3d_121_stride(a,nx-1,ny-1,nz-1,self.npass_smooth[:,js],self.alpha_smooth[:,js],self.stride_smooth[:,js])
    return a[:,0,0]
    
  def getsmoothz(self):
    nz = w3d.nz
    ny = 1
    nx = 1
    a = zeros([nx,ny,nz],'d')
    a[:,:,nz/2]=1.
    nsm = shape(self.npass_smooth)[1]
    for js in range(nsm):
      smooth3d_121_stride(a,nx-1,ny-1,nz-1,self.npass_smooth[:,js],self.alpha_smooth[:,js],self.stride_smooth[:,js])
    return a[0,0,:]
    
  def getsmoothxfftth(self):
    nx = w3d.nx
    theta=(1.+arange(nx))*pi/nx
    nsm = shape(self.npass_smooth)[1]
    f = ones(nx,'d')
    print nx,shape(theta),shape(f)
    for js in range(nsm):
      w = 0.5*(1./self.alpha_smooth[0,js]-1.)
      f *= ((1.+2.*w*cos(theta*self.stride_smooth[0,js]))/(1.+2.*w))**self.npass_smooth[0,js]
    return f,theta
    
  def pltsmx(self,color=black,width=1.):
    pla(self.getsmoothx(),color=color,width=width)
    
  def pltsmxfft(self,color=black,width=1.):
    f=abs(fft.fft(self.getsmoothx()))
    theta=(1.+arange(shape(f)[0]))*2.*pi/shape(f)[0]
    pla(f,2.*pi/theta,color=color,width=width)
    
  def pltsmxfftth(self,color=black,width=1.):
    f,theta = self.getsmoothxfftth()
    pla(f,2.*pi/theta,color=color,width=width)
    
  def pltsmz(self,color=black,width=1.):
    pla(self.getsmoothz(),color=color,width=width)
    
  def pltsmzfft(self,color=black,width=1.):
    f=abs(fft.fft(self.getsmoothz()))
    theta=arange(shape(f)[0])*2.*pi/shape(f)[0]
    pla(f,theta,color=color,width=width)
    
  def binomial_expansion(self,x,y,n):
    # --- we assume that x is real and y is an array of coefficients
    # initializes f
    f = 0.*y
    # initializes Pascal triangle
    pt = zeros(n)
    pt[0] = 1.
    for i in range(n):
      f += pt[i]*x**(n-i)*self.getpower(y,i)
      if i<n-1:pt[1:] = pt[1:]+pt[:-1]
    return f
    
  def getpower(self,y,n):
    # computes coefficients of y**n where y is a list of coefficients
    o = shape(y)[0]
    f = y.copy()
    for k in range(n):
      x = f.copy()
      for i in range(o):
        for j in range(o):
          if (i+j)<o:
            f[i+j] += y[i]*x[j]
    return f

  def get_binomial_filter_factors(self,n,o):
    # set costheta as order n expansion of cos(theta)
    costheta = zeros(o,'d')
    costheta[0] = 1.
    s = -1.
    f = 2.
    for i in range(2,o,2):
      costheta[i] = s / f
      f = f * (f+1.) * (f+2.)
      s*=-1.
      
    # get ((1 + cos(theta))/2)^n
    cf = self.binomial_expansion(1.,costheta,n)/(2.**n)

    # get coeffs for 1 pass compensation
    a = cf[2]
    wcomp1 = -a/(2*a-1.)
    coefsm1 = 1./(1.+2.*wcomp1)

    # get coeffs for 2 pass compensation
    a = cf[0]*2**n
    b = cf[2]*2**n
    c = cf[4]*2**n
    A = a**2-4.*a*c-7.*a*b/3.+4.*b**2
    B = a*b/12.+a*c-b**2
    alpha = -a+2*b
    beta = b-4.*(B/A)*(b-a)
    gamma = (a-2.*b)*B/A
    delta = sqrt(beta**2-4.*alpha*gamma)
    print a,b,c,A,B,alpha,beta,gamma,delta
    
    print '&&&', beta**2-4.*alpha*gamma,beta**2,4.*alpha*gamma
    w1 = (-beta+delta)/(2.*alpha)
    w2 = -B/(w1*A)
    
    coefsm2_1 = 1./(1.+2.*w1)
    coefsm2_2 = 1./(1.+2.*w2)

    return cf,wcomp1,coefsm1,w1,w2,coefsm2_1,coefsm2_2
    
    
    
    
    
    
  def get_binomial_filter_factorsold(self,n,o):
    # set costheta as order n expansion of cos(theta)
    costheta = zeros(o,'d')
    costheta[0] = 1.
    s = -1.
    f = 2.
    for i in range(2,o,2):
      costheta[i] = s / f
      f = f * (f+1.) * (f+2.)
      s*=-1.
      
    # get ((1 + cos(theta))/2)^n
    cf = self.binomial_expansion(1.,costheta,n)/(2.**n)

    # get coeffs for 1 pass compensation
    a = cf[2]
    wcomp1 = -a/(2*a-1.)
    coefsm1 = 1./(1.+2.*wcomp1)

    # get coeffs for 2 pass compensation
    alpha = cf[2]
    beta  = cf[4]
    a1 = b1 = -1.-2.*alpha
    c1 = -4.*(1.+alpha)
    d1 = -alpha
    a2 = b2 = 1./12.-2.*beta
    c2 = 4.*(1./3.-beta)
    d2 = -beta
    A = (c1*b2-c2*b1)/(c2*a1-c1*a2)
    B = (c1*d2-c2*d1)/(c2*a1-c1*a2)
    a = c1*A
    b = a1*A+b1+c1*B
    c = a1*B+d1
    delta = sqrt(b**2-4.*a*c)
    print a1,b1,c1,a2,b2,c2
    print '***',delta,b**2-4.*a*c,b,a,c,A,B
    w2 = (-b+delta)/(2.*a)
    w1 = A*w2+B
    
    coefsm2_1 = 1./(1.+2.*w1)
    coefsm2_2 = 1./(1.+2.*w2)

    return cf,wcomp1,coefsm1,w1,w2,coefsm2_1,coefsm2_2
    
  def fetche(self,*args,**kw):
#    import traceback as tb
#    tb.print_stack()
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
                      solvergeom=self.solvergeom,conductors=conductorobject)

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

  def move_window_fieldsold(self):
    if self.refinement is None:
      # --- move window in x
      self.xgridcont+=self.vxgrid*top.dt
      while (abs(self.xgrid-self.xgridcont)>=0.5*self.dx):
        self.shift_cells_x(1)
        dx = self.dx*sign(self.vxgrid)
        w3d.xmmin+=dx
        w3d.xmmax+=dx
        w3d.xmminp+=dx
        w3d.xmmaxp+=dx
        w3d.xmminlocal+=dx
        w3d.xmmaxlocal+=dx
        w3d.xmminglobal+=dx
        w3d.xmmaxglobal+=dx
        top.xpmin+=dx
        top.xpmax+=dx
        top.xpminlocal+=dx
        top.xpmaxlocal+=dx
      # --- move window in z
      if (abs(top.zgrid-self.zgrid)>=0.5*self.dz):
        shift_em3dblock_ncells_z(self.block,1)
        self.zgrid+=self.dz
        self.nzshifts+=1
    else:
      # --- move window in x
      fc = self.field_coarse
      fc.vxgrid=self.vxgrid
      fc.xgridcont+=fc.vxgrid*top.dt
      while (abs(fc.xgrid-fc.xgridcont)>=0.5*fc.dx):
        fc.shift_cells_x(1)
        self.shift_cells_x(self.refinement[0])
      # --- move window in z
      if (abs(top.zgrid-fc.zgrid)>=0.5*fc.dz):
        shift_em3dblock_ncells_z(fc.block,1)
        shift_em3dblock_ncells_z(self.block,self.refinement[0])
        fc.zgrid+=fc.dz
        fc.nzshifts+=1
        self.zgrid+=self.dz*self.refinement[2]
        self.nzshifts+=self.refinement[2]

  def move_window_fields(self):
      # --- move window in x
      self.xgridcont+=self.vxgrid*top.dt
      while (abs(self.xgrid-self.xgridcont)>=0.5*self.dx):
        self.shift_cells_x(int(sign(self.vxgrid)))
        dx = self.dx*sign(self.vxgrid)
        w3d.xmmin+=dx
        w3d.xmmax+=dx
        w3d.xmminp+=dx
        w3d.xmmaxp+=dx
        w3d.xmminlocal+=dx
        w3d.xmmaxlocal+=dx
        w3d.xmminglobal+=dx
        w3d.xmmaxglobal+=dx
        top.xpmin+=dx
        top.xpmax+=dx
        top.xpminlocal+=dx
        top.xpmaxlocal+=dx
      # --- move window in z
#      while (abs(top.zgrid-self.zgrid)>=0.5*self.dz):
      while ((top.zgrid-self.zgrid)>=0.5*self.dz):
        self.shift_cells_z(1)

  def shift_cells_x(self,n):
        shift_em3dblock_ncells_x(self.block,n)
        dx = self.dx*n
        self.xgrid+=dx
        self.xmmin+=dx
        self.xmmax+=dx
        self.xmminlocal+=dx
        self.xmmaxlocal+=dx
        self.fields.xmin+=dx
        self.fields.xmax+=dx
        self.block.xmin+=dx
        self.block.xmax+=dx
        self.nxshifts+=n

  def shift_cells_z(self,n):
        shift_em3dblock_ncells_z(self.block,n)
        self.zgrid+=self.dz*n
        self.nzshifts+=1

  def solve2ndhalf(self):
    if self.solveroff:return
    if self.mode==2:
      self.solve2ndhalfmode2()
      return
    if self.l_verbose:print 'solve 2nd half',self
    if top.dt<>self.dtinit:raise('Time step has been changed since initialization of EM3D.')
    self.push_b_part_2()
    if self.l_pushf:self.exchange_f()
    em3d_exchange_b(self.block)
    self.move_window_fields()
    if self.l_verbose:print 'solve 2nd half done'

  def dosolve(self,iwhich=0,*args):
    if self.solveroff:return
    if self.mode==2:
      self.dosolvemode2()
      return
    if any(top.fselfb<>0.):raise('Error:EM solver does not work if fselfb<>0.')
    if self.l_verbose:print 'solve 1st half'
    if top.dt<>self.dtinit:raise('Time step has been changed since initialization of EM3D.')
    self.push_e()
    self.exchange_e()
    for i in range(int(self.ncyclesperstep)-1):
      self.push_b_full()
      if self.l_pushf:self.exchange_f()
      self.exchange_b()
      self.push_e_full(i)
      self.exchange_e()
    self.push_b_part_1()
    if self.l_pushf:self.exchange_f()
    self.exchange_b()
    self.setebp()
    if top.efetch[0]<>4:self.yee2node3d()
    if self.l_smooth_particle_fields and any(self.npass_smooth>0):
       self.smoothfields()
    # --- for fields that are overcycled, they need to be pushed backward every ncyclesperstep
    self.push_e(dir=-1)
    self.exchange_e(dir=-1)
    if self.l_verbose:print 'solve 1st half done'

  def push_e(self,dir=1.):  
    dt = dir*top.dt/self.ncyclesperstep
    if self.novercycle==1:
      if dir>0.:
        doit=True
      else:
        doit=False
    else:
      if self.icycle==0 or (self.icycle==self.novercycle-1 and dir>0.):
        doit=True
      else:
        doit=False
    if doit:
      self.add_laser(self.fields)
      if dir<0.:
        self.fields.Ex_inz*=-1.
        self.fields.Ey_inz*=-1.
      if self.l_verbose:print 'push_e',self,dt,top.it,self.icycle
      push_em3d_eef(self.block,dt,0,self.l_pushf,self.l_pushpot)
      if dir<0.:
        self.fields.Ex_inz*=-1.
        self.fields.Ey_inz*=-1.
    if self.refinement is not None:
      self.__class__.__bases__[1].push_e(self.field_coarse,dir)

  def setebp(self):  
#    if (self.refinement is not None) or \
#       (self.novercycle>1) or \
#       (self.l_smooth_particle_fields and any(self.npass_smooth>0)):
    setebp(self.fields,self.icycle,self.novercycle)
    if self.refinement is not None:
      self.__class__.__bases__[1].setebp(self.field_coarse)

  def exchange_e(self,dir=1.):  
    if self.novercycle==1:
      if dir>0.:
        doit=True
      else:
        doit=False
    else:
      if self.icycle==0 or (self.icycle==self.novercycle-1 and dir>0.):
        doit=True
      else:
        doit=False
    if doit:
      em3d_exchange_e(self.block)
    if self.refinement is not None:
      self.__class__.__bases__[1].exchange_e(self.field_coarse)

  def push_b_part_1(self,dir=1.):  
    dt = dir*top.dt/self.ncyclesperstep
    if self.novercycle==1:
      if dir>0.:
        doit=True
      else:
        doit=False
    else:
      if self.icycle==0 or (self.icycle==self.novercycle-1 and dir>0.):
        doit=True
      else:
        doit=False
    if doit:
      if self.l_verbose:print 'push_b part 1',self,dt,top.it,self.icycle,dir
      push_em3d_bf(self.block,dt,1,self.l_pushf,self.l_pushpot)
    if self.refinement is not None:
      self.__class__.__bases__[1].push_b_part_1(self.field_coarse,dir)

  def push_b_part_2(self):
#    if top.efetch[0]<>4 and (self.refinement is None) and not \
#       (self.l_smooth_particle_fields and any(self.npass_smooth>0)):node2yee3d(self.block.core.yf)
    dt = top.dt/self.ncyclesperstep
    if self.ncyclesperstep<1.:
      self.novercycle = nint(1./self.ncyclesperstep)
      self.icycle = (top.it-1)%self.novercycle
    else:
      self.novercycle = 1
      self.icycle = 0
    if self.icycle==0:
      if self.l_verbose:print 'push_b part 2',self,dt,top.it,self.icycle
      push_em3d_bf(self.block,dt,2,self.l_pushf,self.l_pushpot)
    if self.refinement is not None:
      self.__class__.__bases__[1].push_b_part_2(self.field_coarse)

  def yee2node3d(self):  
    yee2node3d(self.block.core.yf)
    if self.refinement is not None:
      self.__class__.__bases__[1].yee2node3d(self.field_coarse)

  def node2yee3d(self):  
    node2yee3d(self.block.core.yf)
    if self.refinement is not None:
      self.__class__.__bases__[1].node2yee3d(self.field_coarse)

  def exchange_b(self,dir=1.):  
    if self.novercycle==1:
      if dir>0.:
        doit=True
      else:
        doit=False
    else:
      if self.icycle==0 or (self.icycle==self.novercycle-1 and dir>0.):
        doit=True
      else:
        doit=False
    if doit:
      if self.l_verbose:print 'exchange_b',self,top.it,self.icycle
      em3d_exchange_b(self.block)
    if self.refinement is not None:
      self.__class__.__bases__[1].exchange_b(self.field_coarse,dir)

  def exchange_f(self,dir=1.):  
    if self.novercycle==1:
      if dir>0.:
        doit=True
      else:
        doit=False
    else:
      if self.icycle==0 or (self.icycle==self.novercycle-1 and dir>0.):
        doit=True
      else:
        doit=False
    if doit:
      em3d_exchange_f(self.block)
    if self.refinement is not None:
      self.__class__.__bases__[1].exchange_f(self.field_coarse,dir)

  def push_b_full(self):  
    dt = top.dt/self.ncyclesperstep
    if self.l_verbose:print 'push_b full',self,dt,top.it,self.icycle
    push_em3d_bf(self.block,dt,0,self.l_pushf,self.l_pushpot)

  def push_e_full(self,i):  
    dt = top.dt/self.ncyclesperstep
    if self.l_pushf:
      w = float(i+2)/self.ncyclesperstep
      self.fields.Rho = (1.-w)*self.fields.Rhoold + w*self.fields.Rhoarray[...,0] 
    if self.l_verbose:print 'push_e full',self,dt,top.it,self.icycle
    self.add_laser(self.fields)
    push_em3d_eef(self.block,dt,0,self.l_pushf,self.l_pushpot)

  def dosolvemode2(self,iwhich=0,*args):
    if self.solveroff:return
    if any(top.fselfb<>0.):raise('Error:EM solver does not work if fselfb<>0.')
    if self.l_verbose:print 'solve 1st half'
    if top.dt<>self.dtinit:raise('Time step has been changed since initialization of EM3D.')
    dt=top.dt*2
    if self.odd:
      push_em3d_bf(self.block,dt,1,self.l_pushf,self.l_pushpot)
      if self.l_pushf:self.exchange_f()
      self.exchange_b()
      self.add_laser(self.fields)
      push_em3d_eef(self.block,dt,2,self.l_pushf,self.l_pushpot)
      self.exchange_e()
    else:
      self.add_laser(self.fields)
      push_em3d_eef(self.block,dt,1,self.l_pushf,self.l_pushpot)
      self.exchange_e()
      push_em3d_bf(self.block,dt,2,self.l_pushf,self.l_pushpot)
      if self.l_pushf:self.exchange_f()
      self.exchange_b()
    self.odd = 1-self.odd
    if not all(top.efetch==top.efetch[0]):raise('Error:top.efetch must have same value for every species when using EM solver.')
    self.setebp()
    if top.efetch[0]<>4:self.yee2node3d()
    if self.l_smooth_particle_fields and any(self.npass_smooth>0):
       self.smoothfields()
    
  def solve2ndhalfmode2(self):
    if self.solveroff:return
    if self.l_verbose:print 'solve 2nd half',self
    if top.dt<>self.dtinit:raise('Time step has been changed since initialization of EM3D.')
#    if top.efetch[0]<>4:node2yee3d(self.block.core.yf)
    self.move_window_fields()
    if self.ncyclesperstep<1.:
      self.novercycle = nint(1./self.ncyclesperstep)
      self.icycle = (top.it-1)%self.novercycle
    else:
      self.novercycle = 1
      self.icycle = 0
    return
    dt = top.dt
    if not self.odd:
      push_em3d_eef(self.block,dt,2,self.l_pushf,self.l_pushpot)
      push_em3d_bf(self.block,dt,1,self.l_pushf,self.l_pushpot)
    else:
      push_em3d_bf(self.block,dt,2,self.l_pushf,self.l_pushpot)
      self.add_laser(self.fields)
      push_em3d_eef(self.block,dt,1,self.l_pushf,self.l_pushpot)
    self.odd = 1-self.odd
    if self.l_verbose:print 'solve 2nd half done'

  def dosolvemode2old(self,iwhich=0,*args):
    if self.solveroff:return
    if any(top.fselfb<>0.):raise('Error:EM solver does not work if fselfb<>0.')
    if self.l_verbose:print 'solve 1st half'
    if top.dt<>self.dtinit:raise('Time step has been changed since initialization of EM3D.')
    dt = top.dt*2
    if self.odd:
      push_em3d_bf(self.block,dt,1,self.l_pushf,self.l_pushpot)
      push_em3d_eef(self.block,dt,2,self.l_pushf,self.l_pushpot)
    else:
      self.add_laser(self.fields)
      push_em3d_eef(self.block,dt,1,self.l_pushf,self.l_pushpot)
      push_em3d_bf(self.block,dt,2,self.l_pushf,self.l_pushpot)
    self.odd = 1-self.odd
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
    
  def dosolvemode2vold(self,iwhich=0,*args):
    if self.solveroff:return
    if any(top.fselfb<>0.):raise('Error:EM solver does not work if fselfb<>0.')
    if self.l_verbose:print 'solve 1st half'
    if top.dt<>self.dtinit:raise('Time step has been changed since initialization of EM3D.')
    dt = top.dt*2
    if self.odd:
      push_em3d_eef(self.block,dt,2,self.l_pushf,self.l_pushpot)
      push_em3d_bf(self.block,dt,1,self.l_pushf,self.l_pushpot)
    else:
      push_em3d_bf(self.block,dt,2,self.l_pushf,self.l_pushpot)
      self.add_laser(self.fields)
      push_em3d_eef(self.block,dt,1,self.l_pushf,self.l_pushpot)
    self.odd = 1-self.odd
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
                    adjust=0,labelscale=0.5,color='auto',ncolor=None,cmin=None,cmax=None,cscale=1.,l_csym=0,
                    procs=None,
                    **kw):
    if direction is None and not l_opyndx:direction=2
    if self.l_2dxz:direction=1
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
    if kw.has_key('gridscale'):
      gridscale=kw['gridscale']
    else:
      gridscale=None
    if self.isactive:
     f=self.block.core.yf
     nxd,nyd,nzd=shape(data)
     if direction==0:
      if slice is None:
        if self.l4symtry:
          slice=0
        else:
          slice=self.nx/2
      xslice = w3d.xmmin+slice*w3d.dx
      selfslice = nint((xslice-self.block.xmin)/self.block.dx)
      if selfslice<0 or selfslice>nxd-1:
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
          data=transpose(data[selfslice,:,:])
        else:
          data=data[selfslice,:,:]
     if direction==1:
      if slice is None:
        if self.l4symtry:
          slice=0
        else:
          slice=self.ny/2
      yslice = w3d.ymmin+slice*w3d.dy
      selfslice = nint((yslice-self.block.ymin)/self.block.dy)
      if self.l_2dxz:slice=0
      if selfslice<0 or selfslice>nyd-1:
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
          data=transpose(data[:,selfslice,:])
        else:
          data=data[:,selfslice,:]
     if direction==2:
      if slice is None:slice=self.nz/2
      zslice = w3d.zmmin+slice*w3d.dz
      selfslice = nint((zslice-self.block.zmin)/self.block.dz)
      if selfslice<0 or selfslice>nzd-1:
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
          data=transpose(data[:,:,selfslice])
        else:
          data=data[:,:,selfslice]
    if gridscale is not None:
      data=data.copy()*gridscale
      kw['gridscale']=None
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
    if cscale is not None:
      cmin *= cscale
      cmax *= cscale
    if l_csym:
      cmax = max(abs(cmin),abs(cmax))
      cmin = -cmax
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
        mpi.send(self.isactive,0,3)
        if self.isactive:
          mpi.send((xmin,xmax,ymin,ymax,data),0,3)
      else:
        if me in procs and data is not None:
          kw.setdefault('xmin',xmin)
          kw.setdefault('xmax',xmax)
          kw.setdefault('ymin',ymin)
          kw.setdefault('ymax',ymax)
          ppgeneric(grid=data,**kw)
        for i in range(1,npes):
          isactive = mpirecv(i,3)
          if isactive:
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
      dx=self.block.dx
      dy=self.block.dy
      dz=self.block.dz
      if me>0 and me in procs:
        mpi.send((xmin,xmax,dx,ymin,ymax,dy,zmin,zmax,dz,data),0,3)
      else:
        if ncolor is None:ncolor = 2
        dec = (cmax-cmin)/ncolor
        isos = cmin+arange(ncolor)*dec+dec/2
        origins = [xmin*xscale,ymin*yscale,zmin*zscale]
        deltas = [dx*xscale,dy*yscale,dz*zscale]
        if color=='auto':
          color = []
          for i in range(ncolor):
            r = array([1.,0.,0.])
            g = array([0.,1.,0.])
            b = array([0.,0.,1.])
            if i<=ncolor/4:
              color.append(b+g*(i*4./ncolor))
            elif i<=ncolor/2:
              color.append(b*(2.-(i*4./ncolor))+g)
            elif i<=3*ncolor/4:
              color.append(g+r*((i*4./ncolor)-2.))
            else:
              color.append(g*(4.-(i*4./ncolor))+r)
        colormap,opacity  = DXColormap(data=isos,
                                       ncolors=ncolor,
                                       colors=color,
                                       opacitystart=None,opacityend=None,opacities=None)
        DXReference(colormap)
        e3d,colorbar = viewisosurface1(data,isos,color=color,display=0,
                        origins=origins,
                        deltas=deltas,
                        colormap=colormap)
        dxob = [e3d,colorbar]
        for i in range(1,npes):
          xminp,xmaxp,dxp,yminp,ymaxp,dyp,zminp,zmaxp,dzp,data=mpirecv(i,3)
          origins = [xminp*xscale,yminp*yscale,zminp*zscale]
          deltas = [dxp*xscale,dyp*yscale,dzp*zscale]
          DXReference(colormap)
          e3d,colorbar = viewisosurface1(data,isos,color=color,display=0,
                          origins=origins,
                          deltas=deltas,
                          colormap=colormap)
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
    return slice
    
  ##########################################################################
  # Gather requested array on processor 0
  def gatherarray(self,data,direction=None,slice=None,procs=None,guards=0,**kw):
    if self.l_2dxz:direction=1
    if self.isactive:
     f=self.block.core.yf
     nxd,nyd,nzd=shape(data)
     if direction==0:
      if slice is None:
        if self.l4symtry:
          slice=0
        else:
          slice=self.nx/2
      xslice = w3d.xmmin+slice*w3d.dx
      selfslice = nint((xslice-self.block.xmin)/self.block.dx)
      if selfslice<0 or selfslice>nxd-1:
        data=None
        xmin=xmax=ymin=ymax=0.
      else:
        xmin=self.block.ymin
        xmax=self.block.ymax
        ymin=self.block.zmin
        ymax=self.block.zmax
        data=data[selfslice,:,:]
     if direction==1:
      if slice is None:
        if self.l4symtry:
          slice=0
        else:
          slice=self.ny/2
      yslice = w3d.ymmin+slice*w3d.dy
      selfslice = nint((yslice-self.block.ymin)/self.block.dy)
      if self.l_2dxz:slice=0
      if selfslice<0 or selfslice>nyd-1:
        data=None
        xmin=xmax=ymin=ymax=0.
      else:
        xmin=self.block.xmin
        xmax=self.block.xmax
        ymin=self.block.zmin
        ymax=self.block.zmax
        data=data[:,selfslice,:]
     if direction==2:
      if slice is None:slice=self.nz/2
      zslice = w3d.zmmin+slice*w3d.dz
      selfslice = nint((zslice-self.block.zmin)/self.block.dz)
      if selfslice<0 or selfslice>nzd-1:
        data=None
        xmin=xmax=ymin=ymax=0.
      else:
        xmin=self.block.xmin
        xmax=self.block.xmax
        ymin=self.block.ymin
        ymax=self.block.ymax
        data=data[:,:,selfslice]
    if procs is None:procs=arange(npes)
    if me==0:
      if guards:
        nxg=self.nxguard
        nyg=self.nyguard
        nzg=self.nzguard
      else:
        nxg=0
        nyg=0
        nzg=0
    if direction in [0,1,2]:
      if me==0:
        if direction==0: datag = zeros([self.ny+1+nyg*2,self.nz+1+nzg*2],'d')
        if direction==1: datag = zeros([self.nx+1+nxg*2,self.nz+1+nzg*2],'d')
        if direction==2: datag = zeros([self.nx+1+nxg*2,self.ny+1+nyg*2],'d')
      else:
        datag=None
      if type(procs) is type(1):procs=[procs]
      validdata = gatherlist(self.isactive and data is not None,bcast=1)
      validprocs = compress(validdata,arange(len(validdata)))
      alldata = gatherlist([xmin,ymin,data],dest=0,procs=validprocs)
      barrier() # this ensures that processor 0 will not bet overflowed with messages
      if me==0:
        for i in range(len(alldata)):
          xminp = alldata[i][0]
          yminp = alldata[i][1]
          data  = alldata[i][-1]
          nx,ny = shape(data)
          if direction==0:
            ixmin = nint((xminp-self.ymmin)/self.dy)+nyg
            iymin = nint((yminp-self.zmmin)/self.dz)+nzg
          if direction==1:
            ixmin = nint((xminp-self.xmmin)/self.dx)+nxg
            iymin = nint((yminp-self.zmmin)/self.dz)+nzg
          if direction==2:
            ixmin = nint((xminp-self.xmmin)/self.dx)+nxg
            iymin = nint((yminp-self.ymmin)/self.dy)+nyg
          datag[ixmin:ixmin+nx,iymin:iymin+ny] = data[...]
    else:
      if me==0:
        datag = zeros([self.nx+1+nxg*2,self.ny+1+nyg*2,self.nz+1+nzg*2],'d')
      else:
        datag = None
      xmin=self.block.xmin
      xmax=self.block.xmax
      ymin=self.block.ymin
      ymax=self.block.ymax
      zmin=self.block.zmin
      zmax=self.block.zmax
      dx=self.block.dx
      dy=self.block.dy
      dz=self.block.dz
      if me>0 and me in procs:
        mpi.send((xmin,xmax,dx,ymin,ymax,dy,zmin,zmax,dz,data),0,3)
      else:
        for i in range(0,npes):
          if i<>me:
            xminp,xmaxp,dxp,yminp,ymaxp,dyp,zminp,zmaxp,dzp,data=mpirecv(i,3)
          else:
            xminp = xmin
            yminp = ymin
            zminp = zmin
          if data is not None:
            nx,ny,nz = shape(data)
            ixmin = nint((xminp-self.xmmin)/self.dx)+nxg
            iymin = nint((yminp-self.ymmin)/self.dy)+nyg
            izmin = nint((zminp-self.zmmin)/self.dz)+nzg
            datag[ixmin:ixmin+nx,iymin:iymin+ny,izmin:izmin+nz] = data[...]
      barrier() # this ensures that processor 0 will not bet overflowed with messages

    return datag
    
  def getarray(self,g,guards=0,overlap=0):
    if guards:
      return g
    else:
      f=self.fields
      ox=oy=oz=0
      if not overlap:
        if self.block.xrbnd==em3d.otherproc:ox=1
        if self.block.yrbnd==em3d.otherproc:oy=1
        if self.block.zrbnd==em3d.otherproc:oz=1
      if self.l_2dxz:
        return g[f.nxguard:-f.nxguard-ox,0:1,f.nzguard:-f.nzguard-oz]
      else:
        return g[f.nxguard:-f.nxguard-ox,f.nyguard:-f.nyguard-oy,f.nzguard:-f.nzguard-oz]
      
  def pfex(self,l_children=1,guards=0,direction=None,**kw):
      self.genericpfem3d(self.getarray(self.fields.Exp,guards,overlap=direction is None),'E_x',
      direction=direction,**kw)

  def pfey(self,l_children=1,guards=0,direction=None,**kw):
      self.genericpfem3d(self.getarray(self.fields.Eyp,guards,overlap=direction is None),'E_y',
      direction=direction,**kw)

  def pfez(self,l_children=1,guards=0,direction=None,**kw):
      self.genericpfem3d(self.getarray(self.fields.Ezp,guards,overlap=direction is None),'E_z',
      direction=direction,**kw)

  def pfbx(self,l_children=1,guards=0,direction=None,**kw):
      self.genericpfem3d(self.getarray(self.fields.Bxp,guards,overlap=direction is None),'B_x',
      direction=direction,**kw)

  def pfby(self,l_children=1,guards=0,direction=None,**kw):
      self.genericpfem3d(self.getarray(self.fields.Byp,guards,overlap=direction is None),'B_y',
      direction=direction,**kw)

  def pfbz(self,l_children=1,guards=0,direction=None,**kw):
      self.genericpfem3d(self.getarray(self.fields.Bzp,guards,overlap=direction is None),'B_z',
      direction=direction,**kw)

  def pfexp(self,l_children=1,guards=0,direction=None,**kw):
      self.genericpfem3d(self.getarray(self.fields.Ex,guards,overlap=direction is None),'Eg_x',
      direction=direction,**kw)

  def pfeyg(self,l_children=1,guards=0,direction=None,**kw):
      self.genericpfem3d(self.getarray(self.fields.Ey,guards,overlap=direction is None),'Eg_y',
      direction=direction,**kw)

  def pfezg(self,l_children=1,guards=0,direction=None,**kw):
      self.genericpfem3d(self.getarray(self.fields.Ez,guards,overlap=direction is None),'Eg_z',
      direction=direction,**kw)

  def pfbxg(self,l_children=1,guards=0,direction=None,**kw):
      self.genericpfem3d(self.getarray(self.fields.Bx,guards,overlap=direction is None),'Bg_x',
      direction=direction,**kw)

  def pfbyg(self,l_children=1,guards=0,direction=None,**kw):
      self.genericpfem3d(self.getarray(self.fields.By,guards,overlap=direction is None),'Bg_y',
      direction=direction,**kw)

  def pfbzg(self,l_children=1,guards=0,direction=None,**kw):
      self.genericpfem3d(self.getarray(self.fields.Bz,guards,overlap=direction is None),'Bg_z',
      direction=direction,**kw)

  def pfjx(self,l_children=1,guards=0,direction=None,**kw):
      self.genericpfem3d(self.getarray(self.fields.J[:,:,:,0],guards,overlap=direction is None),'J_x',
      direction=direction,**kw)

  def pfjy(self,l_children=1,guards=0,direction=None,**kw):
      self.genericpfem3d(self.getarray(self.fields.J[:,:,:,1],guards,overlap=direction is None),'J_y',
      direction=direction,**kw)

  def pfjz(self,l_children=1,guards=0,direction=None,**kw):
      self.genericpfem3d(self.getarray(self.fields.J[:,:,:,2],guards,overlap=direction is None),'J_z',
      direction=direction,**kw)

  def pfrho(self,l_children=1,guards=0,direction=None,**kw):
      self.genericpfem3d(self.getarray(self.fields.Rho,guards,overlap=direction is None),'Rho',
      direction=direction,**kw)

  def pff(self,l_children=1,guards=0,direction=None,**kw):
      self.genericpfem3d(self.getarray(self.fields.F,guards,overlap=direction is None),'F',
      direction=direction,**kw)

  def pfdive(self,l_children=1,guards=0,direction=None,**kw):
      self.genericpfem3d(self.getarray(self.getdive(),guards,overlap=direction is None),'div(E)',
      direction=direction,**kw)

  def pfe(self,l_children=1,guards=0,direction=None,**kw):
      e = self.getarray(self.fields.Exp**2+self.fields.Eyp**2+self.fields.Ezp**2,guards,overlap=direction is None)
      self.genericpfem3d(sqrt(e),'E',
      direction=direction,**kw)

  def pfb(self,l_children=1,guards=0,direction=None,**kw):
      b = self.getarray(self.fields.Bxp**2+self.fields.Byp**2+self.fields.Bzp**2,guards,overlap=direction is None)
      self.genericpfem3d(sqrt(b),'B',
      direction=direction,**kw)

  def pfw(self,l_children=1,guards=0,direction=None,**kw):
      e = self.getarray(self.fields.Exp**2+self.fields.Eyp**2+self.fields.Ezp**2,guards,overlap=direction is None)
      b = self.getarray(self.fields.Bxp**2+self.fields.Byp**2+self.fields.Bzp**2,guards,overlap=direction is None)
      self.genericpfem3d(sqrt(e+b),'W',
      direction=direction,**kw)

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

  def getjx(self,guards=0,overlap=0):
      return self.getarray(self.fields.J[:,:,:,0],guards,overlap)

  def getjy(self,guards=0,overlap=0):
      return self.getarray(self.fields.J[:,:,:,1],guards,overlap)

  def getjz(self,guards=0,overlap=0):
      return self.getarray(self.fields.J[:,:,:,2],guards,overlap)

  def getex(self,guards=0,overlap=0):
      return self.getarray(self.fields.Exp,guards,overlap)

  def getey(self,guards=0,overlap=0):
      return self.getarray(self.fields.Eyp,guards,overlap)
        
  def getez(self,guards=0,overlap=0):
      return self.getarray(self.fields.Ezp,guards,overlap)
        
  def getbx(self,guards=0,overlap=0):
      return self.getarray(self.fields.Bxp,guards,overlap)

  def getby(self,guards=0,overlap=0):
      return self.getarray(self.fields.Byp,guards,overlap)
        
  def getbz(self,guards=0,overlap=0):
      return self.getarray(self.fields.Bzp,guards,overlap)
        
  def getexg(self,guards=0,overlap=0):
      return self.getarray(self.fields.Ex,guards,overlap)

  def geteyg(self,guards=0,overlap=0):
      return self.getarray(self.fields.Ey,guards,overlap)
        
  def getezg(self,guards=0,overlap=0):
      return self.getarray(self.fields.Ez,guards,overlap)
        
  def getbxg(self,guards=0,overlap=0):
      return self.getarray(self.fields.Bx,guards,overlap)

  def getbyg(self,guards=0,overlap=0):
      return self.getarray(self.fields.By,guards,overlap)
        
  def getbzg(self,guards=0,overlap=0):
      return self.getarray(self.fields.Bz,guards,overlap)
        
  def getrho(self,guards=0,overlap=0):
      return self.getarray(self.fields.Rho,guards,overlap)

  def getf(self,guards=0,overlap=0):
      return self.getarray(self.fields.F,guards,overlap)

  def getdive(self,guards=0,overlap=0):
      dive = zeros(shape(self.fields.Ex),'d')
      f = self.fields
      if top.efetch[0]<>4:node2yee3d(f)
      if self.l_2dxz:
        dive[1:-1,1:-1,1:-1] = (f.Ex[1:-1,0,1:-1]-f.Ex[:-2,0,1:-1])/f.dx \
                             + (f.Ez[1:-1,0,1:-1]-f.Ez[1:-1,0,:-2])/f.dz 
      else:
        dive[1:-1,1:-1,1:-1] = (f.Ex[1:-1,1:-1,1:-1]-f.Ex[:-2,1:-1,1:-1])/f.dx \
                             + (f.Ey[1:-1,1:-1,1:-1]-f.Ey[1:-1,:-2,1:-1])/f.dy \
                             + (f.Ez[1:-1,1:-1,1:-1]-f.Ez[1:-1,1:-1,:-2])/f.dz 
      if top.efetch[0]<>4:yee2node3d(f)
      return dive
      
  def gete(self,guards=0,overlap=0):
      return self.getarray(self.fields.Exp**2+self.fields.Eyp**2+self.fields.Ezp**2,guards,overlap)

  def getb(self,guards=0,overlap=0):
      return self.getarray(self.fields.Bxp**2+self.fields.Byp**2+self.fields.Bzp**2,guards,overlap)

  def getw(self,guards=0,overlap=0):
      e2 = self.getarray(self.fields.Exp**2+self.fields.Eyp**2+self.fields.Ezp**2,guards,overlap)
      b2 = self.getarray(self.fields.Bxp**2+self.fields.Byp**2+self.fields.Bzp**2,guards,overlap)
      return sqrt(e2+clight**2*b2)

  def getw2(self,guards=0,overlap=0):
      e2 = self.getarray(self.fields.Exp**2+self.fields.Eyp**2+self.fields.Ezp**2,guards,overlap)
      b2 = self.getarray(self.fields.Bxp**2+self.fields.Byp**2+self.fields.Bzp**2,guards,overlap)
      return e2+clight**2*b2

  def geteg(self,guards=0,overlap=0):
      return self.getarray(self.fields.Ex**2+self.fields.Ey**2+self.fields.Ez**2,guards,overlap)

  def getbg(self,guards=0,overlap=0):
      return self.getarray(self.fields.Bx**2+self.fields.By**2+self.fields.Bz**2,guards,overlap)

  def getwg(self,guards=0,overlap=0):
      e2 = self.getarray(self.fields.Ex**2+self.fields.Ey**2+self.fields.Ez**2,guards,overlap)
      b2 = self.getarray(self.fields.Bx**2+self.fields.By**2+self.fields.Bz**2,guards,overlap)
      return sqrt(e2+clight**2*b2)

  def getwg2(self,guards=0,overlap=0):
      e2 = self.getarray(self.fields.Ex**2+self.fields.Ey**2+self.fields.Ez**2,guards,overlap)
      b2 = self.getarray(self.fields.Bx**2+self.fields.By**2+self.fields.Bz**2,guards,overlap)
      return e2+clight**2*b2

  def gatherex(self,guards=0,direction=None,**kw):
      return self.gatherarray(self.getex(guards,overlap=direction is None),direction=direction,**kw)

  def gatherey(self,guards=0,direction=None,**kw):
      return self.gatherarray(self.getey(guards,overlap=direction is None),direction=direction,**kw)

  def gatherez(self,guards=0,direction=None,**kw):
      return self.gatherarray(self.getez(guards,overlap=direction is None),direction=direction,**kw)

  def gatherbx(self,guards=0,direction=None,**kw):
      return self.gatherarray(self.getbx(guards,overlap=direction is None),direction=direction,**kw)

  def gatherby(self,guards=0,direction=None,**kw):
      return self.gatherarray(self.getby(guards,overlap=direction is None),direction=direction,**kw)

  def gatherbz(self,guards=0,direction=None,**kw):
      return self.gatherarray(self.getbz(guards,overlap=direction is None),direction=direction,**kw)

  def gatherjx(self,guards=0,direction=None,**kw):
      return self.gatherarray(self.getjx(guards,overlap=direction is None),direction=direction,**kw)

  def gatherjy(self,guards=0,direction=None,**kw):
      return self.gatherarray(self.getjy(guards,overlap=direction is None),direction=direction,**kw)

  def gatherjz(self,guards=0,direction=None,**kw):
      return self.gatherarray(self.getjz(guards,overlap=direction is None),direction=direction,**kw)

  def gatherrho(self,guards=0,direction=None,**kw):
      return self.gatherarray(self.getrho(guards,overlap=direction is None),direction=direction,**kw)

  def gatherf(self,guards=0,direction=None,**kw):
      return self.gatherarray(self.getf(guards,overlap=direction is None),direction=direction,**kw)

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
      if self.l_2dxz:
        q.append( (sum(Ex[1:-1,0,1:-1]-Ex[:-2,0,1:-1])/f.dx \
                +  sum(Ez[1:-1,0,1:-1]-Ez[1:-1,0,:-2])/f.dz) \
                * eps0*(f.dx*f.dz))
      else:
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
      if self.l_2dxz:
        q.append(sum(f.Rho)*(f.dx*f.dz))
      else:
        q.append(sum(f.Rho)*(f.dx*f.dy*f.dz))
    return q

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

  def connectblocks(self,blo,bup,dir):
    if dir==0:
      # --- X
      loxropen = blo.xrbnd == openbc
      loylopen = blo.ylbnd == openbc
      loyropen = blo.yrbnd == openbc
      lozlopen = blo.zlbnd == openbc
      lozropen = blo.zrbnd == openbc

      upxlopen = bup.xlbnd == openbc
      upylopen = bup.ylbnd == openbc
      upyropen = bup.yrbnd == openbc
      upzlopen = bup.zlbnd == openbc
      upzropen = bup.zrbnd == openbc

      # --- sides
      if loxropen: del blo.sidexr.syf # --- deallocate upper x PML block
      if upxlopen: del bup.sidexl.syf # --- deallocate lower x PML block
      blo.sidexr = bup.core       
#      bup.sidexl = blo.core
      bup.sidexl = EM3D_FIELDtype()
      # --- edges
      if loxropen and loylopen:del blo.edgexryl.syf
      if loxropen and loyropen:del blo.edgexryr.syf
      if loxropen and lozlopen:del blo.edgexrzl.syf
      if loxropen and lozropen:del blo.edgexrzr.syf
      if upxlopen and upylopen:del bup.edgexlyl.syf
      if upxlopen and upyropen:del bup.edgexlyr.syf
      if upxlopen and upzlopen:del bup.edgexlzl.syf
      if upxlopen and upzropen:del bup.edgexlzr.syf
      if loylopen:
        blo.edgexryl = bup.sideyl
      else:
        blo.edgexryl = EM3D_FIELDtype()
      if loyropen:
        blo.edgexryr = bup.sideyr
      else:
        blo.edgexryr = EM3D_FIELDtype()
      if lozlopen:
        blo.edgexrzl = bup.sidezl
      else:
        blo.edgexrzl = EM3D_FIELDtype()
      if lozropen:
        blo.edgexrzr = bup.sidezr
      else:
        blo.edgexrzr = EM3D_FIELDtype()
#      bup.edgexlyl = blo.sideyl
#      bup.edgexlyr = blo.sideyr
#      bup.edgexlzl = blo.sidezl
#      bup.edgexlzr = blo.sidezr
      bup.edgexlyl = EM3D_FIELDtype()
      bup.edgexlyr = EM3D_FIELDtype()
      bup.edgexlzl = EM3D_FIELDtype()
      bup.edgexlzr = EM3D_FIELDtype()
      # --- corners
      if loxropen and loylopen and lozlopen:del blo.cornerxrylzl.syf
      if loxropen and loyropen and lozlopen:del blo.cornerxryrzl.syf
      if loxropen and loylopen and lozropen:del blo.cornerxrylzr.syf
      if loxropen and loyropen and lozropen:del blo.cornerxryrzr.syf
      if upxlopen and upylopen and upzlopen:del bup.cornerxlylzl.syf
      if upxlopen and upyropen and upzlopen:del bup.cornerxlyrzl.syf
      if upxlopen and upylopen and upzropen:del bup.cornerxlylzr.syf
      if upxlopen and upyropen and upzropen:del bup.cornerxlyrzr.syf
      if loylopen and lozlopen:
        blo.cornerxrylzl = bup.edgeylzl
      else:
        blo.cornerxrylzl = EM3D_FIELDtype()
      if loyropen and lozlopen:
        blo.cornerxryrzl = bup.edgeyrzl
      else:
        blo.cornerxryrzl = EM3D_FIELDtype()
      if loylopen and lozropen:
        blo.cornerxrylzr = bup.edgeylzr
      else:
        blo.cornerxrylzr = EM3D_FIELDtype()
      if loyropen and lozropen:
        blo.cornerxryrzr = bup.edgeyrzr
      else:
        blo.cornerxryrzr = EM3D_FIELDtype()
#      bup.cornerxlylzl = blo.edgeylzl
#      bup.cornerxlyrzl = blo.edgeyrzl
#      bup.cornerxlylzr = blo.edgeylzr
#      bup.cornerxlyrzr = blo.edgeyrzr
      bup.cornerxlylzl = EM3D_FIELDtype()
      bup.cornerxlyrzl = EM3D_FIELDtype()
      bup.cornerxlylzr = EM3D_FIELDtype()
      bup.cornerxlyrzr = EM3D_FIELDtype()

      blo.xrbnd = em3d.otherblock
      bup.xlbnd = em3d.otherblock

    if dir==1:
      # --- Y
      loyropen = blo.yrbnd == openbc
      lozlopen = blo.zlbnd == openbc
      lozropen = blo.zrbnd == openbc
      loxlopen = blo.xlbnd == openbc
      loxropen = blo.xrbnd == openbc

      upylopen = bup.ylbnd == openbc
      upzlopen = bup.zlbnd == openbc
      upzropen = bup.zrbnd == openbc
      upxlopen = bup.xlbnd == openbc
      upxropen = bup.xrbnd == openbc

      # --- sides
      if loyropen: del blo.sideyr.syf # --- deallocate upper y PML block
      if upylopen: del bup.sideyl.syf # --- deallocate lower y PML block
      blo.sideyr = bup.core       
#      bup.sideyl = blo.core
      bup.sideyl = EM3D_FIELDtype()
      # --- edges
      if loyropen and lozlopen:del blo.edgeyrzl.syf
      if loyropen and lozropen:del blo.edgeyrzr.syf
      if loyropen and loxlopen:del blo.edgexlyr.syf
      if loyropen and loxropen:del blo.edgexryr.syf
      if upylopen and upzlopen:del bup.edgeylzl.syf
      if upylopen and upzropen:del bup.edgeylzr.syf
      if upylopen and upxlopen:del bup.edgexlyl.syf
      if upylopen and upxropen:del bup.edgexryr.syf
      if loxlopen:
        blo.edgexlyr = bup.sidexl
      else:
        blo.edgexlyr = EM3D_FIELDtype()
      if loxropen:
        blo.edgexryr = bup.sidexr
      else:
        blo.edgexryr = EM3D_FIELDtype()
      if lozlopen:
        blo.edgeyrzl = bup.sidezl
      else:
        blo.edgeyrzl = EM3D_FIELDtype()
      if lozropen:
        blo.edgeyrzr = bup.sidezr
      else:
        blo.edgeyrzr = EM3D_FIELDtype()
#      bup.edgexlyl = blo.sidexl
#      bup.edgexryl = blo.sidexr
#      bup.edgeylzl = blo.sidezl
#      bup.edgeylzr = blo.sidezr
      bup.edgexlyl = EM3D_FIELDtype()
      bup.edgexryl = EM3D_FIELDtype()
      bup.edgeylzl = EM3D_FIELDtype()
      bup.edgeylzr = EM3D_FIELDtype()
      # --- corners
      if loyropen and lozlopen and loxlopen:del blo.cornerxlyrzl.syf
      if loyropen and lozropen and loxlopen:del blo.cornerxlyrzr.syf
      if loyropen and lozlopen and loxropen:del blo.cornerxryrzl.syf
      if loyropen and lozropen and loxropen:del blo.cornerxryrzr.syf
      if upylopen and upzlopen and upxlopen:del bup.cornerxlylzl.syf
      if upylopen and upzropen and upxlopen:del bup.cornerxlylzr.syf
      if upylopen and upzlopen and upxropen:del bup.cornerxrylzl.syf
      if upylopen and upzropen and upxropen:del bup.cornerxrylzr.syf
      if loxlopen and lozlopen:
        blo.cornerxlyrzl = bup.edgexlzl
      else:
        blo.cornerxlyrzl = EM3D_FIELDtype()
      if loxropen and lozlopen:
        blo.cornerxryrzl = bup.edgexrzl
      else:
        blo.cornerxryrzl = EM3D_FIELDtype()
      if loxlopen and lozropen:
        blo.cornerxlyrzr = bup.edgexlzr
      else:
        blo.cornerxlyrzr = EM3D_FIELDtype()
      if loxropen and lozropen:
        blo.cornerxryrzr = bup.edgexrzr
      else:
        blo.cornerxryrzr = EM3D_FIELDtype()
#      bup.cornerxlylzl = blo.edgexlzl
#      bup.cornerxrylzl = blo.edgexrzl
#      bup.cornerxlylzr = blo.edgexlzr
#      bup.cornerxrylzr = blo.edgexrzr
      bup.cornerxlylzl = EM3D_FIELDtype()
      bup.cornerxrylzl = EM3D_FIELDtype()
      bup.cornerxlylzr = EM3D_FIELDtype()
      bup.cornerxrylzr = EM3D_FIELDtype()

      blo.yrbnd = em3d.otherblock
      bup.ylbnd = em3d.otherblock

    if dir==2:
      # --- Z
      lozropen = blo.zrbnd == openbc
      loxlopen = blo.xlbnd == openbc
      loxropen = blo.xrbnd == openbc
      loylopen = blo.ylbnd == openbc
      loyropen = blo.yrbnd == openbc

      upzlopen = bup.zlbnd == openbc
      upxlopen = bup.xlbnd == openbc
      upxropen = bup.xrbnd == openbc
      upylopen = bup.ylbnd == openbc
      upyropen = bup.yrbnd == openbc

      # --- sides
      if lozropen: del blo.sidezr.syf # --- deallocate upper y PML block
      if upzlopen: del bup.sidezl.syf # --- deallocate lower y PML block
      blo.sidezr = bup.core       
#      bup.sidezl = blo.core
      bup.sidezl = EM3D_FIELDtype()
      # --- edges
      if lozropen and loxlopen:del blo.edgexlzr.syf
      if lozropen and loxropen:del blo.edgexrzr.syf
      if lozropen and loylopen:del blo.edgeylzr.syf
      if lozropen and loyropen:del blo.edgeyrzr.syf
      if upzlopen and upxlopen:del bup.edgexlzl.syf
      if upzlopen and upxropen:del bup.edgexrzl.syf
      if upzlopen and upylopen:del bup.edgeylzl.syf
      if upzlopen and upyropen:del bup.edgeyrzl.syf
      if loxlopen:
        blo.edgexlzr = bup.sidexl
      else:
        blo.edgexlzr = EM3D_FIELDtype()
      if loxropen:
        blo.edgexrzr = bup.sidexr
      else:
        blo.edgexrzr = EM3D_FIELDtype()
      if loylopen:
        blo.edgeylzr = bup.sideyl
      else:
        blo.edgeylzr = EM3D_FIELDtype()
      if loyropen:
        blo.edgeyrzr = bup.sideyr
      else:
        blo.edgeyrzr = EM3D_FIELDtype()
 #     bup.edgexlzl = blo.sidexl
 #     bup.edgexrzl = blo.sidexr
 #     bup.edgeylzl = blo.sideyl
 #     bup.edgeyrzl = blo.sideyr
      bup.edgexlzl = EM3D_FIELDtype()
      bup.edgexrzl = EM3D_FIELDtype()
      bup.edgeylzl = EM3D_FIELDtype()
      bup.edgeyrzl = EM3D_FIELDtype()
      # --- corners
      if lozropen and loxlopen and loylopen:del blo.cornerxlylzr.syf
      if lozropen and loxropen and loylopen:del blo.cornerxrylzr.syf
      if lozropen and loxlopen and loyropen:del blo.cornerxlyrzr.syf
      if lozropen and loxropen and loyropen:del blo.cornerxryrzr.syf
      if upzlopen and upxlopen and upylopen:del bup.cornerxlylzl.syf
      if upzlopen and upxropen and upylopen:del bup.cornerxrylzl.syf
      if upzlopen and upxlopen and upyropen:del bup.cornerxlyrzl.syf
      if upzlopen and upxropen and upyropen:del bup.cornerxryrzl.syf
      if loxlopen and loylopen:
        blo.cornerxlylzr = bup.edgexlyl
      else:
        blo.cornerxlylzr = EM3D_FIELDtype()
      if loxropen and loylopen:
        blo.cornerxrylzr = bup.edgexryl
      else:
        blo.cornerxrylzr = EM3D_FIELDtype()
      if loxlopen and loyropen:
        blo.cornerxlyrzr = bup.edgexlyr
      else:
        blo.cornerxlyrzr = EM3D_FIELDtype()
      if loxropen and loyropen:
        blo.cornerxryrzr = bup.edgexryr
      else:
        blo.cornerxryrzr = EM3D_FIELDtype()
   #   bup.cornerxlylzl = blo.edgexlyl
   #   bup.cornerxrylzl = blo.edgexryl
   #   bup.cornerxlyrzl = blo.edgexlyr
   #   bup.cornerxryrzl = blo.edgexryr
      bup.cornerxlylzl = EM3D_FIELDtype()
      bup.cornerxrylzl = EM3D_FIELDtype()
      bup.cornerxlyrzl = EM3D_FIELDtype()
      bup.cornerxryrzl = EM3D_FIELDtype()

      blo.zrbnd = em3d.otherblock
      bup.zlbnd = em3d.otherblock

  def setdtinit(self):
    self.dtinit=top.dt

  def setbcoverlaps(self):
    xlbndmatch = False
    xrbndmatch = False
    ylbndmatch = False
    yrbndmatch = False
    zlbndmatch = False
    zrbndmatch = False
    for n in self.overlapshigher.keys():
      block = self.root.listofblocks[n].block
      blockcoarse = self.root.listofblocks[n].field_coarse.block
      lower = self.overlapshigher[n][0]
      upper = self.overlapshigher[n][1]
      # --- check grid matching in upper x
      if self.fullupper[0] == lower[0] and \
         all(self.fulllower[1:] == lower[1:]) and \
         all(self.fullupper[1:] == upper[1:]):
         if not (self.block.sidexr is block):
           self.connectblocks(self.block,block,0)
           self.connectblocks(self.field_coarse.block,blockcoarse,0)
         xrbndmatch = True
      # --- check grid matching in upper y
      if self.fullupper[1] == lower[1] and \
         all(self.fulllower[0::2] == lower[0::2]) and \
         all(self.fullupper[0::2] == upper[0::2]):
         if not (self.block.sideyr is block):
           self.connectblocks(self.block,block,1)
           self.connectblocks(self.field_coarse.block,blockcoarse,1)
         yrbndmatch = True
      # --- check grid matching in upper z
      if self.fullupper[2] == lower[2] and \
         all(self.fulllower[:2] == lower[:2]) and \
         all(self.fullupper[:2] == upper[:2]):
         if not (self.block.sidezr is block):
           self.connectblocks(self.block,block,2)
           self.connectblocks(self.field_coarse.block,blockcoarse,2)
         zrbndmatch = True
    for n in self.overlapslower.keys():
      block = self.root.listofblocks[n].block
      blockcoarse = self.root.listofblocks[n].field_coarse.block
      lower = self.overlapslower[n][0]
      upper = self.overlapslower[n][1]
      # --- check grid matching in lower x
      if self.fulllower[0] == upper[0] and \
         all(self.fulllower[1:] == lower[1:]) and \
         all(self.fullupper[1:] == upper[1:]):
         if not (self.block.sidexl is block):
           self.connectblocks(block,self.block,0)
           self.connectblocks(blockcoarse,self.field_coarse.block,0)
         xlbndmatch = True
      # --- check grid matching in lower y
      if self.fulllower[1] == upper[1] and \
         all(self.fulllower[0::2] == lower[0::2]) and \
         all(self.fullupper[0::2] == upper[0::2]):
         if not (self.block.sideyl is block):
           self.connectblocks(block,self.block,1)
           self.connectblocks(blockcoarse,self.field_coarse.block,1)
         ylbndmatch = True
      # --- check grid matching in lower z
      if self.fulllower[2] == upper[2] and \
         all(self.fulllower[:2] == lower[:2]) and \
         all(self.fullupper[:2] == upper[:2]):
         if not (self.block.sidezl is block):
           self.connectblocks(block,self.block,2)
           self.connectblocks(blockcoarse,self.field_coarse.block,2)
         zlbndmatch = True

    if not xlbndmatch and self.block.xlbnd == em3d.otherblock:
           self.block.xlbnd = openbc
           self.field_coarse.block.xlbnd = openbc
           # --- TODO: allocate openbc
    if not xrbndmatch and self.block.xrbnd == em3d.otherblock:
           self.block.xrbnd = openbc
           self.field_coarse.block.xrbnd = openbc
           # --- TODO: allocate openbc
    if not ylbndmatch and self.block.ylbnd == em3d.otherblock:
           self.block.ylbnd = openbc
           self.field_coarse.block.ylbnd = openbc
           # --- TODO: allocate openbc
    if not yrbndmatch and self.block.yrbnd == em3d.otherblock:
           self.block.yrbnd = openbc
           self.field_coarse.block.yrbnd = openbc
           # --- TODO: allocate openbc
    if not zlbndmatch and self.block.zlbnd == em3d.otherblock:
           self.block.zlbnd = openbc
           self.field_coarse.block.zlbnd = openbc
           # --- TODO: allocate openbc
    if not zrbndmatch and self.block.zrbnd == em3d.otherblock:
           self.block.zrbnd = openbc
           self.field_coarse.block.zrbnd = openbc
           # --- TODO: allocate openbc

    self.bounds = array([self.block.xlbnd,
                         self.block.xrbnd,
                         self.block.ylbnd,
                         self.block.yrbnd,
                         self.block.zlbnd,
                         self.block.zrbnd])
    if self.refinement is not None:
      self.field_coarse.bounds = array([self.block.xlbnd,
                           self.block.xrbnd,
                           self.block.ylbnd,
                           self.block.yrbnd,
                           self.block.zlbnd,
                           self.block.zrbnd])
      
  def checkconnections(self):
    yeefieldtype = -1
    if self.block.xlbnd==openbc:
      if self.block.edgexlyl.fieldtype==yeefieldtype:
         self.block.edgexlyl=EM3D_FIELDtype()
         self.block.cornerxlylzl=EM3D_FIELDtype()
         self.block.cornerxlylzr=EM3D_FIELDtype()
      if self.block.edgexlyr.fieldtype==yeefieldtype:
         self.block.edgexlyr=EM3D_FIELDtype()
         self.block.cornerxlyrzl=EM3D_FIELDtype()
         self.block.cornerxlyrzr=EM3D_FIELDtype()
      if self.block.edgexlzl.fieldtype==yeefieldtype:
         self.block.edgexlzl=EM3D_FIELDtype()
         self.block.cornerxlylzl=EM3D_FIELDtype()
         self.block.cornerxlyrzl=EM3D_FIELDtype()
      if self.block.edgexlzr.fieldtype==yeefieldtype:
         self.block.edgexlzr=EM3D_FIELDtype()
         self.block.cornerxlylzr=EM3D_FIELDtype()
         self.block.cornerxlyrzr=EM3D_FIELDtype()
    if self.block.xrbnd==openbc:
      if self.block.edgexryl.fieldtype==yeefieldtype:
         self.block.edgexryl=EM3D_FIELDtype()
         self.block.cornerxrylzl=EM3D_FIELDtype()
         self.block.cornerxrylzr=EM3D_FIELDtype()
      if self.block.edgexryr.fieldtype==yeefieldtype:
         self.block.edgexryr=EM3D_FIELDtype()
         self.block.cornerxryrzl=EM3D_FIELDtype()
         self.block.cornerxryrzr=EM3D_FIELDtype()
      if self.block.edgexrzl.fieldtype==yeefieldtype:
         self.block.edgexrzl=EM3D_FIELDtype()
         self.block.cornerxrylzl=EM3D_FIELDtype()
         self.block.cornerxryrzl=EM3D_FIELDtype()
      if self.block.edgexrzr.fieldtype==yeefieldtype:
         self.block.edgexrzr=EM3D_FIELDtype()
         self.block.cornerxrylzr=EM3D_FIELDtype()
         self.block.cornerxryrzr=EM3D_FIELDtype()
                 
    if self.block.ylbnd==openbc:
      if self.block.edgexlyl.fieldtype==yeefieldtype:
         self.block.edgexlyl=EM3D_FIELDtype()
         self.block.cornerxlylzl=EM3D_FIELDtype()
         self.block.cornerxlylzr=EM3D_FIELDtype()
      if self.block.edgexryl.fieldtype==yeefieldtype:
         self.block.edgexryl=EM3D_FIELDtype()
         self.block.cornerxrylzl=EM3D_FIELDtype()
         self.block.cornerxrylzr=EM3D_FIELDtype()
      if self.block.edgeylzl.fieldtype==yeefieldtype:
         self.block.edgeylzl=EM3D_FIELDtype()
         self.block.cornerxlylzl=EM3D_FIELDtype()
         self.block.cornerxrylzl=EM3D_FIELDtype()
      if self.block.edgeylzr.fieldtype==yeefieldtype:
         self.block.edgeylzr=EM3D_FIELDtype()
         self.block.cornerxlylzr=EM3D_FIELDtype()
         self.block.cornerxrylzr=EM3D_FIELDtype()
    if self.block.yrbnd==openbc:
      if self.block.edgexlyr.fieldtype==yeefieldtype:
         self.block.edgexlyr=EM3D_FIELDtype()
         self.block.cornerxlyrzl=EM3D_FIELDtype()
         self.block.cornerxlyrzr=EM3D_FIELDtype()
      if self.block.edgexryr.fieldtype==yeefieldtype:
         self.block.edgexryr=EM3D_FIELDtype()
         self.block.cornerxryrzl=EM3D_FIELDtype()
         self.block.cornerxryrzr=EM3D_FIELDtype()
      if self.block.edgeyrzl.fieldtype==yeefieldtype:
         self.block.edgeyrzl=EM3D_FIELDtype()
         self.block.cornerxlyrzl=EM3D_FIELDtype()
         self.block.cornerxryrzl=EM3D_FIELDtype()
      if self.block.edgeyrzr.fieldtype==yeefieldtype:
         self.block.edgeyrzr=EM3D_FIELDtype()
         self.block.cornerxlyrzr=EM3D_FIELDtype()
         self.block.cornerxryrzr=EM3D_FIELDtype()

    if self.block.zlbnd==openbc:
      if self.block.edgexlzl.fieldtype==yeefieldtype:
         self.block.edgexlzl=EM3D_FIELDtype()
         self.block.cornerxlylzl=EM3D_FIELDtype()
         self.block.cornerxlyrzl=EM3D_FIELDtype()
      if self.block.edgexrzl.fieldtype==yeefieldtype:
         self.block.edgexrzl=EM3D_FIELDtype()
         self.block.cornerxrylzl=EM3D_FIELDtype()
         self.block.cornerxryrzl=EM3D_FIELDtype()
      if self.block.edgeylzl.fieldtype==yeefieldtype:
         self.block.edgeylzl=EM3D_FIELDtype()
         self.block.cornerxlylzl=EM3D_FIELDtype()
         self.block.cornerxrylzl=EM3D_FIELDtype()
      if self.block.edgeyrzl.fieldtype==yeefieldtype:
         self.block.edgeyrzl=EM3D_FIELDtype()
         self.block.cornerxlyrzl=EM3D_FIELDtype()
         self.block.cornerxryrzl=EM3D_FIELDtype()
    if self.block.zrbnd==openbc:
      if self.block.edgexlzr.fieldtype==yeefieldtype:
         self.block.edgexlzr=EM3D_FIELDtype()
         self.block.cornerxlylzr=EM3D_FIELDtype()
         self.block.cornerxlyrzr=EM3D_FIELDtype()
      if self.block.edgexrzr.fieldtype==yeefieldtype:
         self.block.edgexrzr=EM3D_FIELDtype()
         self.block.cornerxrylzr=EM3D_FIELDtype()
         self.block.cornerxryrzr=EM3D_FIELDtype()
      if self.block.edgeylzr.fieldtype==yeefieldtype:
         self.block.edgeylzr=EM3D_FIELDtype()
         self.block.cornerxlylzr=EM3D_FIELDtype()
         self.block.cornerxrylzr=EM3D_FIELDtype()
      if self.block.edgeyrzr.fieldtype==yeefieldtype:
         self.block.edgeyrzr=EM3D_FIELDtype()
         self.block.cornerxlyrzr=EM3D_FIELDtype()
         self.block.cornerxryrzr=EM3D_FIELDtype()

    if self.block.xlbnd==openbc and  self.block.ylbnd==openbc:
      if self.block.cornerxlylzl.fieldtype==yeefieldtype:self.block.cornerxlylzl=EM3D_FIELDtype()
      if self.block.cornerxlylzr.fieldtype==yeefieldtype:self.block.cornerxlylzr=EM3D_FIELDtype()
    if self.block.xrbnd==openbc and  self.block.ylbnd==openbc:
      if self.block.cornerxrylzl.fieldtype==yeefieldtype:self.block.cornerxrylzl=EM3D_FIELDtype()
      if self.block.cornerxrylzr.fieldtype==yeefieldtype:self.block.cornerxrylzr=EM3D_FIELDtype()
    if self.block.xlbnd==openbc and  self.block.yrbnd==openbc:
      if self.block.cornerxlyrzl.fieldtype==yeefieldtype:self.block.cornerxlyrzl=EM3D_FIELDtype()
      if self.block.cornerxlyrzr.fieldtype==yeefieldtype:self.block.cornerxlyrzr=EM3D_FIELDtype()
    if self.block.xrbnd==openbc and  self.block.yrbnd==openbc:
      if self.block.cornerxryrzl.fieldtype==yeefieldtype:self.block.cornerxryrzl=EM3D_FIELDtype()
      if self.block.cornerxryrzr.fieldtype==yeefieldtype:self.block.cornerxryrzr=EM3D_FIELDtype()

    if self.block.xlbnd==openbc and  self.block.zlbnd==openbc:
      if self.block.cornerxlylzl.fieldtype==yeefieldtype:self.block.cornerxlylzl=EM3D_FIELDtype()
      if self.block.cornerxlyrzl.fieldtype==yeefieldtype:self.block.cornerxlyrzl=EM3D_FIELDtype()
    if self.block.xrbnd==openbc and  self.block.zlbnd==openbc:
      if self.block.cornerxrylzl.fieldtype==yeefieldtype:self.block.cornerxlylzl=EM3D_FIELDtype()
      if self.block.cornerxryrzl.fieldtype==yeefieldtype:self.block.cornerxlyrzl=EM3D_FIELDtype()
    if self.block.xlbnd==openbc and  self.block.zrbnd==openbc:
      if self.block.cornerxlylzr.fieldtype==yeefieldtype:self.block.cornerxlylzr=EM3D_FIELDtype()
      if self.block.cornerxlyrzr.fieldtype==yeefieldtype:self.block.cornerxlyrzr=EM3D_FIELDtype()
    if self.block.xrbnd==openbc and  self.block.zrbnd==openbc:
      if self.block.cornerxrylzr.fieldtype==yeefieldtype:self.block.cornerxlylzr=EM3D_FIELDtype()
      if self.block.cornerxryrzr.fieldtype==yeefieldtype:self.block.cornerxlyrzr=EM3D_FIELDtype()

    if self.block.ylbnd==openbc and  self.block.zlbnd==openbc:
      if self.block.cornerxlylzl.fieldtype==yeefieldtype:self.block.cornerxlylzl=EM3D_FIELDtype()
      if self.block.cornerxrylzl.fieldtype==yeefieldtype:self.block.cornerxrylzl=EM3D_FIELDtype()
    if self.block.yrbnd==openbc and  self.block.zlbnd==openbc:
      if self.block.cornerxlyrzl.fieldtype==yeefieldtype:self.block.cornerxlylzl=EM3D_FIELDtype()
      if self.block.cornerxryrzl.fieldtype==yeefieldtype:self.block.cornerxrylzl=EM3D_FIELDtype()
    if self.block.ylbnd==openbc and  self.block.zrbnd==openbc:
      if self.block.cornerxlylzr.fieldtype==yeefieldtype:self.block.cornerxlylzr=EM3D_FIELDtype()
      if self.block.cornerxrylzr.fieldtype==yeefieldtype:self.block.cornerxrylzr=EM3D_FIELDtype()
    if self.block.yrbnd==openbc and  self.block.zrbnd==openbc:
      if self.block.cornerxlyrzr.fieldtype==yeefieldtype:self.block.cornerxlylzr=EM3D_FIELDtype()
      if self.block.cornerxryrzr.fieldtype==yeefieldtype:self.block.cornerxrylzr=EM3D_FIELDtype()

    if self.refinement is not None:
      self.field_coarse.checkconnections()

  def fillchilddomains(self):
    parent = self.root.listofblocks[self.parents[0]]
    xl,yl,zl = self.fullloweroverrefinement-parent.fulllower
#    xu,yu,zu = self.fullupperoverrefinement-parent.fulllower+[1,1,1]
    xu,yu,zu = self.fullupperoverrefinement-parent.fulllower+[0,0,0]

    if self.bounds[0]==openbc:
      xlguard = self.nguard[0]
      xlguarddepos = self.nguarddepos[0]
    else:
      xlguard = 0
      xlguarddepos = 0
    if self.bounds[1]==openbc:
      xrguard = self.nguard[0]
      xrguarddepos = self.nguarddepos[0]
    else:
      xrguard = 0
      xrguarddepos = 0

    if self.bounds[2]==openbc:
      ylguard = self.nguard[1]
      ylguarddepos = self.nguarddepos[1]
    else:
      ylguard = 0
      ylguarddepos = 0
    if self.bounds[3]==openbc:
      yrguard = self.nguard[1]
      yrguarddepos = self.nguarddepos[1]
    else:
      yrguard = 0
      yrguarddepos = 0

    if self.bounds[4]==openbc:
      zlguard = self.nguard[2]
      zlguarddepos = self.nguarddepos[2]
    else:
      zlguard = 0
      zlguarddepos = 0
    if self.bounds[5]==openbc:
      zrguard = self.nguard[2]
      zrguarddepos = self.nguarddepos[2]
    else:
      zrguard = 0
      zrguarddepos = 0

    cd = parent.childdomains
    cd[xl:xu,yl:yu,zl:zu]=parent.blocknumber
    cd[xl+xlguarddepos:xu-xrguarddepos,
       yl+ylguarddepos:yu-yrguarddepos,
       zl+zlguarddepos:zu-zrguarddepos]=-self.blocknumber
    cd[xl+xlguard:xu-xrguard,
       yl+ylguard:yu-yrguard,
       zl+zlguard:zu-zrguard]=self.blocknumber
    
  def fillchilddomainsold(self):
    xl,yl,zl = self.fullloweroverrefinement
    xu,yu,zu = self.fullupperoverrefinement+[1,1,1]

    if self.bounds[0]==openbc:
      xlguard = self.nguard[0]
      xlguarddepos = self.nguarddepos[0]
    else:
      xlguard = 0
      xlguarddepos = 0
    if self.bounds[1]==openbc:
      xrguard = self.nguard[0]
      xrguarddepos = self.nguarddepos[0]
    else:
      xrguard = 0
      xrguarddepos = 0

    if self.bounds[2]==openbc:
      ylguard = self.nguard[1]
      ylguarddepos = self.nguarddepos[1]
    else:
      ylguard = 0
      ylguarddepos = 0
    if self.bounds[3]==openbc:
      yrguard = self.nguard[1]
      yrguarddepos = self.nguarddepos[1]
    else:
      yrguard = 0
      yrguarddepos = 0

    if self.bounds[4]==openbc:
      zlguard = self.nguard[2]
      zlguarddepos = self.nguarddepos[2]
    else:
      zlguard = 0
      zlguarddepos = 0
    if self.bounds[5]==openbc:
      zrguard = self.nguard[2]
      zrguarddepos = self.nguarddepos[2]
    else:
      zrguard = 0
      zrguarddepos = 0

    cd = self.root.listofblocks[self.parents[0]].childdomains
    cd[xl:xu,yl:yu,zl:zu]=0
    cd[xl+xlguarddepos:xu-xrguarddepos,
       yl+ylguarddepos:yu-yrguarddepos,
       zl+zlguarddepos:zu-zrguarddepos]=-self.blocknumber
    cd[xl+xlguard:xu-xrguard,
       yl+ylguard:yu-yrguard,
       zl+zlguard:zu-zrguard]=self.blocknumber
    

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
                      ncyclesperstep,
                      l_2dxz,
                      theta_damp,
                      sigmae,
                      sigmab):
                      
  if xlb == em3d.otherproc:
    procxl = top.procneighbors[0,0]
  else:
    procxl = me
  if xrb == em3d.otherproc:
    procxr = top.procneighbors[1,0]
  else:
    procxr = me
  if ylb == em3d.otherproc:
    procyl = top.procneighbors[0,1]
  else:
    procyl = me
  if yrb == em3d.otherproc:
    procyr = top.procneighbors[1,1]
  else:
    procyr = me
  if zlb == em3d.otherproc:
    proczl = top.procneighbors[0,2]
  else:
    proczl = me
  if zrb == em3d.otherproc:
    proczr = top.procneighbors[1,2]
  else:
    proczr = me
  
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
  f.theta_damp=theta_damp
  f.sigmae=sigmae
  f.sigmab=sigmab
  if 0:#refinement is None and (all(npass_smooth==0) or not l_smooth_particle_fields):
    f.nxp = 0
    f.nyp = 0
    f.nzp = 0
  else:
    f.nxp = f.nx
    f.nyp = f.ny
    f.nzp = f.nz
  if ncyclesperstep<0.75:
    f.nxpnext = f.nx
    f.nypnext = f.ny
    f.nzpnext = f.nz
  if not l_pushf:
    f.nxf = 0
    f.nyf = 0
    f.nzf = 0
  else:
    f.nxf = f.nx
    f.nyf = f.ny
    f.nzf = f.nz
  if f.theta_damp<>0.:
    f.nxdamp=f.nx
    f.nydamp=f.ny
    f.nzdamp=f.nz
  if f.sigmae<>0. or f.sigmab<>0.:
    f.nxext=f.nx
    f.nyext=f.ny
    f.nzext=f.nz
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
  
  if 0:#refinement is None and  (all(npass_smooth==0) or not l_smooth_particle_fields):
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
  b.sidexl.proc=me
  if xlb==periodic:
    b.sidexl=b.core
  if xlb==em3d.otherproc:
    allocatef(b.sidexl,stencil)
    b.sidexl.proc=procxl
  b.sidexr = EM3D_FIELDtype()
  if xrb==openbc:
    allocatesf(b.sidexr,stencil)
    init_splitfield(b.sidexr.syf,nbndx,ny,nz,nxguard,nyguard,nzguard, dt, dx, dy, dz, clight, 1, 0, 0, nnx, smaxx, sdeltax, nny, smaxy, sdeltay, nnz, smaxz, sdeltaz, l_2dxz)
  b.sidexr.proc=me
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
  b.sideyl.proc=me
  if ylb==periodic:
    b.sideyl=b.core
  if ylb==em3d.otherproc:
    allocatef(b.sideyl,stencil)
    b.sideyl.proc=procyl
  b.sideyr = EM3D_FIELDtype()
  if yrb==openbc:
    allocatesf(b.sideyr,stencil)
    init_splitfield(b.sideyr.syf,nx,nbndy,nz,nxguard,nyguard,nzguard, dt, dx, dy, dz, clight, 0, 1, 0, nnx, smaxx, sdeltax, nny, smaxy, sdeltay, nnz, smaxz, sdeltaz, l_2dxz)
  b.sideyr.proc=me
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
  b.sidezl.proc=me
  if zlb==periodic:
    b.sidezl=b.core
  if zlb==em3d.otherproc:
    allocatef(b.sidezl,stencil)
    b.sidezl.proc=proczl
  b.sidezr = EM3D_FIELDtype()
  if zrb==openbc:
    allocatesf(b.sidezr,stencil)
    init_splitfield(b.sidezr.syf,nx,ny,nbndz,nxguard,nyguard,nzguard, dt, dx, dy, dz, clight, 0, 0, 1, nnx, smaxx, sdeltax, nny, smaxy, sdeltay, nnz, smaxz, sdeltaz, l_2dxz)
  b.sidezr.proc=me
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

