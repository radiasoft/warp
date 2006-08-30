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
                   'bound0','boundnz','boundxy','l2symtry','l4symtry',
                   'solvergeom']
  __em2dinputs__ = ['l_onegrid','l_copyfields','l_moving_window',
                    'l_elaser_out_plane','ndelta_t',
                   'ntamp_scatter','ntamp_gather']
  __flaginputs__ = {'l_apply_pml':true,'nbndx':10,'nbndy':10,
                    'l_particles_weight':false}

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
      self.dz = (self.zmmaxglobal - self.zmminglobal)/self.nzfull
      self.zmesh = self.zmminglobal + arange(0,self.nzfull+1)*self.dz
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
    wtmp = zerps(n,'d')
    depose_current_em2d(n,x,y,ux,uy,uz,gaminv,wtmp,q*w,top.dt,
                        self.l_particles_weight)

  def fetchefrompositions(self,x,y,ex,ey,ez):
    n = len(x)
    if n == 0: return
    self.bx = zeros(n,'d')
    self.by = zeros(n,'d')
    self.bz = zeros(n,'d')
    geteb_em2d(n,x,y,ex,ey,ez,self.bx,self.by,self.bz)

  def fetchbfrompositions(self,x,z,bx,by,bz):
    # --- This assumes that fetchefrompositions was already called
    n = len(x)
    if n == 0: return
    bx[:] = self.bx
    by[:] = self.by
    bz[:] = self.bz

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
    x,y,ux,uy,uz = self.transformparticles(w3d.xfsapi,w3d.yfsapi,w3d.zfsapi)
    self.fetchefrompositions(x,y,ex,ey,ez)

  def fetchb(self):
    if w3d.api_xlf2:
      w3d.xfsapi=top.pgroup.xp[w3d.ipminapi-1:w3d.ipminapi-1+w3d.ipapi]
      w3d.yfsapi=top.pgroup.yp[w3d.ipminapi-1:w3d.ipminapi-1+w3d.ipapi]
      w3d.zfsapi=top.pgroup.zp[w3d.ipminapi-1:w3d.ipminapi-1+w3d.ipapi]
    bx,by,bz = self.transformfields(w3d.bxfsapi,w3d.byfsapi,w3d.bzfsapi)
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

  def solve(self,iwhich=0):
    grimax(self.field)
    push_em_b(self.field,0.5*top.dt)
    push_em_e(self.field,top.dt)
    push_em_b(self.field,0.5*top.dt)
    move_window_field(self.field)
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

