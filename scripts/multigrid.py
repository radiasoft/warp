"""Class for doing complete multigrid field solve"""
# ToDo:
#  - modify setrho to check if particles are within grid
#  - allow exchanging data with other instances, giving Dirichlet boundaries
#    on common boundaries
#  - incorporate instances into the particle mover, so charge is deposited and
#    the E fields gather appropriately.
from warp import *
from generateconductors import installconductors


##############################################################################
class MultiGrid:
  
  __w3dinputs__ = ['nx','ny','nz',
                   'xmmin','xmmax','ymmin','ymmax','zmmin','zmmax',
                   'bound0','boundnz','boundxy','l2symtry','l4symtry',
                   'solvergeom']
  __f3dinputs__ = ['gridmode','mgparam','downpasses','uppasses',
                   'mgmaxiters','mgtol','mgmaxlevels','mgform',
                   'lcndbndy','icndbndy','laddconductor'] 

  def __init__(self,**kw):

    # --- Kludge - make sure that the multigrid3df routines never sets up
    # --- any conductors.
    f3d.gridmode = 1

    # --- Save input parameters
    for name in MultiGrid.__w3dinputs__:
      #self.__dict__[name] = kw.pop(name,getattr(w3d,name)) # Python2.3
      self.__dict__[name] = kw.get(name,getattr(w3d,name))
      if kw.has_key(name): del kw[name]
    for name in MultiGrid.__f3dinputs__:
      #self.__dict__[name] = kw.pop(name,getattr(f3d,name)) # Python2.3
      self.__dict__[name] = kw.get(name,getattr(f3d,name))
      if kw.has_key(name): del kw[name]

    # --- If there are any remaning keyword arguments, raise an error.
    assert len(kw.keys()) == 0,"Bad keyword arguemnts %s"%kw.keys()

    # --- Set set parallel related paramters and calculate mesh sizes
    self.nzfull = self.nz
    self.zmminglobal = self.zmmin
    self.zmmaxglobal = self.zmmax
    self.dx = (self.xmmax - self.xmmin)/self.nx
    self.dy = (self.ymmax - self.ymmin)/self.ny
    self.dz = (self.zmmax - self.zmmin)/self.nz

    # --- Create phi and rho arrays and other arrays. These are created
    # --- with fortran ordering so no transpose and copy is needed when
    # --- they are passed to fortran.
    self.rho = fzeros((1+self.nx,1+self.ny,1+self.nz),'d')
    self.phi = fzeros((1+self.nx,1+self.ny,3+self.nz),'d')
    self.rstar = fzeros(1+self.nz,'d')

    # --- Create a conductor object, which by default is empty.
    self.conductors = ConductorType()

    self.mgiters = 0
    self.mgerror = 0.

  def copypkgtodict(self,pkg,varlist,dict):
    for name in varlist:
      dict[name] = getattr(pkg,name)

  def copydicttopkg(self,pkg,varlist,dict):
    for name in varlist:
      if type(dict[name]) is ArrayType:
        pkg.forceassign(name,dict[name])
      else:
        setattr(pkg,name,dict[name])

  def loadrho(self,ins_i=-1,nps_i=-1,is_i=-1,lzero=true):
    # --- First, save reference to w3d.rho and other fortran variables.
    w3dvars = ['rho','nx','ny','nz','nzfull','dx','dy','dz',
               'xmmin','ymmin','zmmin','l4symtry','l2symtry']
    w3ddict = {}
    self.copypkgtodict(w3d,w3dvars,w3ddict)
    self.copydicttopkg(w3d,w3dvars,self.__dict__)
    # --- Now load rho
    loadrho(ins_i,nps_i,is_i,lzero)
    # --- Restore w3d variables
    self.copydicttopkg(w3d,w3dvars,w3ddict)

  def fetche(self,ipmin,ip,js,ex,ey,ez):
    # --- First, save reference to w3d.phi and other fortran variables.
    w3ddict = ['phi','nx','ny','nz','nzfull','dx','dy','dz',
               'xmmin','ymmin','zmmin','l4symtry','l2symtry']
    w3ddict = {}
    self.copypkgtodict(w3d,w3dvars,w3ddict)
    self.copydicttopkg(w3d,w3dvars,self.__dict__)
    # --- Now fetch the E field
    fetche3d(ipmin,ip,js,ex,ey,ez)
    # --- Restore w3d variables
    self.copydicttopkg(w3d,w3dvars,w3ddict)

  def installconductor(self,conductor,
                            xmin=None,xmax=None,
                            ymin=None,ymax=None,
                            zmin=None,zmax=None,
                            dfill=top.largepos):
    installconductors(conductor,xmin,xmax,ymin,ymax,zmin,zmax,dfill,
                      top.zbeam,
                      self.nx,self.ny,self.nz,self.nzfull,
                      self.xmmin,self.xmmax,self.ymmin,self.ymmax,
                      self.zmmin,self.zmmax,self.l2symtry,self.l4symtry,
                      solvergeom=self.solvergeom,
                      conductors=self.conductors)

  def solve(self,iwhich=0):
    # --- Note that false is passed in for linbend. At this time, bends
    # --- are not supported.
    multigrid3dsolve(iwhich,self.nx,self.ny,self.nz,self.nzfull,
                     self.dx,self.dy,self.dz,self.phi,self.rho,
                     self.rstar,false,
                     self.bound0,self.boundnz,self.boundxy,
                     self.l2symtry,self.l4symtry,
                     self.xmmin,self.ymmin,self.zmmin,self.zbeam,self.zgrid,
                     self.mgparam,self.mgform,self.mgiters,self.mgmaxiters,
                     self.mgmaxlevels,self.mgerror,self.mgtol,
                     self.downpasses,self.uppasses,
                     self.lcndbndy,self.laddconductor,self.icndbndy,
                     self.gridmode,
                     self.conductors)

  ##########################################################################
  # Define the basic plot commands
  def genericpf(self,kw,pffunc):
    kw['conductors'] = self.conductors
    kw['w3dgrid'] = self
    pffunc(**kw)
  def pfxy(self,**kw): self.genericpf(kw,pfxy)
  def pfzx(self,**kw): self.genericpf(kw,pfzx)
  def pfzy(self,**kw): self.genericpf(kw,pfzy)
  def pfxyg(self,**kw): self.genericpf(kw,pfxyg)
  def pfzxg(self,**kw): self.genericpf(kw,pfzxg)
  def pfzyg(self,**kw): self.genericpf(kw,pfzyg)

