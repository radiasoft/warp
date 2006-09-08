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
                    'tmin_moving_main_window',
                    'l_elaser_out_plane','ndelta_t',
                   'ntamp_scatter','ntamp_gather']
  __flaginputs__ = {'l_apply_pml':true,'nbndx':10,'nbndy':10,
                    'l_particles_weight':false,
                    'laser_amplitude':1.,'laser_profile':None,
                    'laser_gauss_width':None,'laser_angle':0.,
                    'laser_wavelength':None,'laser_wavenumber':None,
                    'laser_frequency':None,'laser_source':2}

  def __init__(self,**kw):
    top.lcallfetchb = true
    #top.bfstype = 12
    top.lcallfetchb = true
    # --- Make sure the refinement is turned off
    em2d.l_onegrid = true

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

    # --- If there are any remaning keyword arguments, raise an error.
    assert len(kw.keys()) == 0,"Bad keyword arguemnts %s"%kw.keys()

    assert self.solvergeom == w3d.XZgeom or not self.l_moving_window,"The moving window can only be used with XZ geometry"

    # --- Calculate mesh sizes
    # --- All cases use X
    self.dx = (self.xmmax - self.xmmin)/self.nx
    self.xsymmetryplane = 0.
    self.xmesh = self.xmmin + arange(0,self.nx+1)*self.dx
    if self.solvergeom in [w3d.XYgeom]:
      # --- When Y is used, nothing special is done
      self.dy = (self.ymmax - self.ymmin)/self.ny
      self.ymesh = self.ymmin + arange(0,self.ny+1)*self.dy
      self.ysymmetryplane = 0.
      self.ymin = self.ymmin
    if self.solvergeom in [w3d.XZgeom]:
      # --- When Z is used, set Z quantities, and set Y quantities to the same
      # --- values. Internally, the class always uses x and y. The user
      # --- interface will depend on solvergeom
      self.dz = (self.zmmax - self.zmmin)/self.nzfull
      self.zmesh = self.zmmin + arange(0,self.nzfull+1)*self.dz
      self.zmeshlocal = self.zmmin + arange(0,self.nz+1)*self.dz
      self.zmin = self.zmmin
      self.ny = self.nz
      self.dy = self.dz
      self.ymesh = self.zmeshlocal
      self.ymin = self.zmin
      self.bounds[2] = self.bounds[4]
      self.bounds[3] = self.bounds[5]
#     self.ndelta_t = top.vbeam*top.dt/w3d.dz
#     em2d.ndelta_t = self.ndelta_t
#     assert (self.ndelta_t - top.vbeam*top.dt/w3d.dz) < 1.e-6,"The input ndelta_t is not commensurate with the values of top.vbeam, top.dt, and w3d.dz"

    # --- Create field and source arrays and other arrays.
    self.allocatefieldarrays()

    # --- Handle laser inputs
    self.setuplaser()

  def __getstate__(self):
    dict = self.__dict__.copy()
    return dict

  def __setstate__(self,dict):
    self.__dict__.update(dict)

  def allocatefieldarrays(self):
    # --- Code transcribed from init_fields
    self.field = EM2D_FIELDtype()
    init_fields(self.field,self.nx,self.ny,self.nbndx,self.nbndx,top.dt,
                self.dx,self.dy,clight,mu0,
                self.xmmin,self.ymmin,1,
                self.bounds[0],self.bounds[1],self.bounds[2],self.bounds[3])
    self.field.js = self.laser_source

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

    # --- Check if laser_profile has a type, is a function, or a table
    self.laser_profile_func = None
    if self.laser_profile == 'gaussian':
      assert self.laser_gauss_width is not None,\
             "For a gaussian laser, the width must be specified using laser_gauss_width"
      f = self.field
      xx = arange(self.nx+4)*f.dx+f.xmin*f.dx - 0.5*f.nx*f.dx
      self.laser_profile = exp(-(xx/self.laser_gauss_width)**2/2.)
    elif operator.isSequenceType(self.laser_profile):
      assert len(self.laser_profile) == self.nx+4,"The specified profile must be of length nx+4"
    elif callable(self.laser_profile):
      self.laser_profile_func = self.laser_profile

  def transformparticles(self,x,y,z,ux=None,uy=None,uz=None):
    if self.solvergeom == w3d.XYgeom:
      return x,y,ux,uy,uz
    elif self.solvergeom == w3d.XZgeom:
      if self.l_moving_window and not self.l_elaser_out_plane:
        return z,x,uz,ux,uy
      else:
        return x,z,ux,uz,uy

  def transformfields(self,fx,fy,fz):
    if self.solvergeom == w3d.XYgeom:
      return fx,fy,fz
    elif self.solvergeom == w3d.XZgeom:
      if self.l_moving_window and not self.l_elaser_out_plane:
        return fz,fx,fy
      else:
        return fx,fz,fy

  def setj(self,x,y,ux,uy,uz,gaminv,q,w):
    n = len(x)
    if n == 0: return
    wtmp = zeros(n,'d')
    em2d_depose_jxjy_esirkepov_linear_serial(self.field.J,n,x,y,ux,uy,uz,
           gaminv,wtmp,q*w,self.field.xmin,self.field.ymin,top.dt,
           self.field.dx,self.field.dy,self.field.nx,self.field.ny,
           self.l_particles_weight)

  def fetchefrompositions(self,x,y,ex,ey,ez):
    n = len(x)
    if n == 0: return
    em2d_gete2d_linear_serial(n,x,y,ex,ey,ez,
                              self.field.xmin,self.field.ymin,
                              self.field.dx,self.field.dy,
                              self.field.nx,self.field.ny,
                              self.field.Ex,self.field.Ey,self.field.Ez)

  def fetchbfrompositions(self,x,y,bx,by,bz):
    # --- This assumes that fetchefrompositions was already called
    n = len(x)
    if n == 0: return
    em2d_getb2d_linear_serial(n,x,y,bx,by,bz,
                              self.field.xmin,self.field.ymin,
                              self.field.dx,self.field.dy,
                              self.field.nx,self.field.ny,
                              self.field.Bx,self.field.By,self.field.Bz)

  def fetchphifrompositions(self,x,z,phi):
    pass

  def loadrho(self,lzero=0):
    pass

  def loadj(self,ins_i=-1,nps_i=-1,is_i=-1,lzero=true):
    if lzero: self.field.J[...] = 0.
    for i,n,q,w in zip(top.pgroup.ins-1,top.pgroup.nps,
                       top.pgroup.sq,top.pgroup.sw):
      if n == 0: continue
      x,y,ux,uy,uz = self.transformparticles(
            top.pgroup.xp[i:i+n],top.pgroup.yp[i:i+n],top.pgroup.zp[i:i+n],
            top.pgroup.uxp[i:i+n],top.pgroup.uyp[i:i+n],top.pgroup.uzp[i:i+n])
      self.setj(x,y,ux,uy,uz,top.pgroup.gaminv[i:i+n],q,w)
    self.smoothdensity()

  def smoothdensity(self):
    smooth2d_lindman(self.field.J[:,:,0],self.field.nx,self.field.ny)
    smooth2d_lindman(self.field.J[:,:,1],self.field.nx,self.field.ny)
    smooth2d_lindman(self.field.J[:,:,2],self.field.nx,self.field.ny)

  def fetche(self):
    if w3d.api_xlf2:
      w3d.xfsapi=top.pgroup.xp[w3d.ipminapi-1:w3d.ipminapi-1+w3d.ipapi]
      w3d.yfsapi=top.pgroup.yp[w3d.ipminapi-1:w3d.ipminapi-1+w3d.ipapi]
      w3d.zfsapi=top.pgroup.zp[w3d.ipminapi-1:w3d.ipminapi-1+w3d.ipapi]
    ex,ey,ez = self.transformfields(w3d.exfsapi,w3d.eyfsapi,w3d.ezfsapi)
    ex[:] = 0.
    ey[:] = 0.
    ez[:] = 0.
    x,y,ux,uy,uz = self.transformparticles(w3d.xfsapi,w3d.yfsapi,w3d.zfsapi)
    self.fetchefrompositions(x,y,ex,ey,ez)

  def fetchb(self):
    if w3d.api_xlf2:
      w3d.xfsapi=top.pgroup.xp[w3d.ipminapi-1:w3d.ipminapi-1+w3d.ipapi]
      w3d.yfsapi=top.pgroup.yp[w3d.ipminapi-1:w3d.ipminapi-1+w3d.ipapi]
      w3d.zfsapi=top.pgroup.zp[w3d.ipminapi-1:w3d.ipminapi-1+w3d.ipapi]
    bx,by,bz = self.transformfields(w3d.bxfsapi,w3d.byfsapi,w3d.bzfsapi)
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
    if not self.l_moving_window or ((top.it%self.ndelta_t)!=0): return
    if top.time < self.tmin_moving_main_window: return
    move_window_field(self.field)

  def add_laser(self):
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
      assert len(self.laser_profile) == self.nx+4,"The specified profile must be of length nx+4"

    if (self.l_elaser_out_plane):
      xx = (arange(self.nx+4) - 0.5)*self.field.dx+self.field.xmin
    else:
      xx = (arange(self.nx+3) - 0.5)*self.field.dx+self.field.xmin

    if self.laser_frequency is not None:
      phase = (xx*sin(self.laser_angle)/clight-top.time)*self.laser_frequency
    else:
      phase = 0.

    if (self.l_elaser_out_plane):
      self.field.Ez_in = self.laser_amplitude*self.laser_profile*cos(phase)
    else:
      self.field.Bz_in = self.laser_amplitude*self.laser_profile[:-1]*cos(phase)

  def solve(self,iwhich=0):
    self.add_laser()
    grimax(self.field)
    push_em_b(self.field,0.5*top.dt)
    push_em_e(self.field,top.dt)
    push_em_b(self.field,0.5*top.dt)
    self.move_window_fields()
    griuni(self.field)


  ##########################################################################
  # Define the basic plot commands
  def genericfp(self,data,title,l_transpose=true,**kw):
      f=self.field
      settitles(title,'X','Y','t = %gs'%(top.time))
      if l_transpose:
        kw.setdefault('xmin',self.xmmin)
        kw.setdefault('xmax',self.xmmax)
        kw.setdefault('ymin',self.ymmin)
        kw.setdefault('ymax',self.ymmax)
        ppgeneric(gridt=data,**kw)
      else:
        kw.setdefault('xmin',self.ymmin)
        kw.setdefault('xmax',self.ymmax)
        kw.setdefault('ymin',self.xmmin)
        kw.setdefault('ymax',self.xmmax)
        ppgeneric(grid=data,**kw)
      
  def fpex(self,**kw):
      self.genericfp(self.field.Ex,'E_x',**kw)

  def fpey(self,**kw):
      self.genericfp(self.field.Ey,'E_y',**kw)

  def fpez(self,**kw):
      self.genericfp(self.field.Ez,'E_z',**kw)

  def fpbx(self,**kw):
      self.genericfp(self.field.Bx,'B_x',**kw)

  def fpby(self,**kw):
      self.genericfp(self.field.By,'B_y',**kw)

  def fpbz(self,**kw):
      self.genericfp(self.field.Bz,'B_z',**kw)

  def fpjx(self,**kw):
      self.genericfp(self.field.J[:,:,0],'J_x',**kw)

  def fpjy(self,**kw):
      self.genericfp(self.field.J[:,:,1],'J_y',**kw)

  def fpjz(self,**kw):
      self.genericfp(self.field.J[:,:,2],'J_z',**kw)

# --- This can only be done after MultiGrid is defined.
try:
  psyco.bind(EM2D)
except NameError:
  pass

