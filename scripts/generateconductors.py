"""
This module contains classes for generating the conductor data from a
combination of simple geometrical elements.
The following elements are defined:

Box(xsize,ysize,zsize,voltage=0.,xcent=0.,ycent=0.,zcent=0.)
Cylinder(radius,length,theta=0.,phi=0.,voltage=0.,xcent=0.,ycent=0.,zcent=0.)
ZCylinder(radius,length,voltage=0.,xcent=0.,ycent=0.,zcent=0.)
YCylinder(radius,length,voltage=0.,xcent=0.,ycent=0.,zcent=0.)
XCylinder(radius,length,voltage=0.,xcent=0.,ycent=0.,zcent=0.)
Sphere(radius,voltage=0.,xcent=0.,ycent=0.,zcent=0.)
ZCone(r_zmin,r_zmax,length,voltage=0.,xcent=0.,ycent=0.,zcent=0.)
ZTorus(r1,r2,voltage=0.,xcent=0.,ycent=0.,zcent=0.)

installconductors(a): generates the data needed for the fieldsolve
"""

# The following classes are also defined but should not be directly used.
# Grid
# Assembly
# AssemblyNot
# AssemblyAnd
# AssemblyPlus
# AssemblyMinus
# Delta

from warp import *

generateconductorsversion = "$Id: generateconductors.py,v 1.2 2002/06/20 20:56:59 dave Exp $"
def generateconductors_doc():
  import generateconductors
  print generateconductors.__doc__

##############################################################################
def installconductors(a):
  # First, create a grid object
  g = Grid()
  # Generate the conductor data
  g.getdata(a)
  # Then install it
  g.installdata()

##############################################################################
##############################################################################
##############################################################################
class Assembly:
  """
Class to hold assemblies of conductors.  Base class of all conductors.
  """

  voltage = 0.
  xcent = 0.
  ycent = 0.
  zcent = 0.

  def __init__(self,v=0.,x=0.,y=0.,z=0.):
    self.voltage = v
    self.xcent = x
    self.ycent = y
    self.zcent = z

  def distance(self,ix,iy,iz,xx,yy,zz):
    "Return distances to conductor from point."
    d = Delta()
    return d

  # Operations which return an Assembly expression.
  def __mul__(self,right):
    return AssemblyAnd(self,right)
  def __add__(self,right):
    return AssemblyPlus(self,right)
  def __sub__(self,right):
    return AssemblyMinus(self,right)
  def __neg__(self):
    return AssemblyNot(self)


class AssemblyNot(Assembly):
  """
AssemblyNot class.  Represents 'not' of assemblies.
  """
  def __init__(self,l):
    Assembly.__init__(self,0.,l.xcent,l.ycent,l.zcent)
    self.left = l
  def distance(self,ix,iy,iz,xx,yy,zz):
    return (-(self.left.distance(ix,iy,iz,xx,yy,zz)))


class AssemblyAnd(Assembly):
  """
AssemblyAnd class.  Represents 'and' of assemblies.
  """
  def __init__(self,l,r):
    Assembly.__init__(self,0.,l.xcent,l.ycent,l.zcent)
    self.left = l
    self.right = r
  def distance(self,ix,iy,iz,xx,yy,zz):
    return (self.left.distance(ix,iy,iz,xx,yy,zz) *
            self.right.distance(ix,iy,iz,xx,yy,zz))


class AssemblyPlus(Assembly):
  """
AssemblyPlus class.  Represents 'or' of assemblies.
  """
  def __init__(self,l,r):
    Assembly.__init__(self,0.,l.xcent,l.ycent,l.zcent)
    self.left = l
    self.right = r
  def distance(self,ix,iy,iz,xx,yy,zz):
    return (self.left.distance(ix,iy,iz,xx,yy,zz) +
            self.right.distance(ix,iy,iz,xx,yy,zz))


class AssemblyMinus(Assembly):
  """
AssemblyMinus class.
  """
  def __init__(self,l,r):
    Assembly.__init__(self,0.,l.xcent,l.ycent,l.zcent)
    self.left = l
    self.right = r
  def distance(self,ix,iy,iz,xx,yy,zz):
    return (self.left.distance(ix,iy,iz,xx,yy,zz) -
            self.right.distance(ix,iy,iz,xx,yy,zz))

##############################################################################

class Delta:
  """
Class to hold the set of distances in each of the six directions.
Distances have the sign of the outward normal surface vector, i.e.
distances to outside the surface are positive, inside negative.
  """

  def __init__(self,ix=None,iy=None,iz=None,xx=None,yy=None,zz=None,
                    dels=None,vs=None,ns=None,
                    parity=None,voltage=0.,condid=0,generator=None,kwlist=[]):
    if ix is None:
      self.ndata = 0
      nn = 10000
      self.ix = zeros(nn)
      self.iy = zeros(nn)
      self.iz = zeros(nn)
      self.xx = zeros(nn,'d')
      self.yy = zeros(nn,'d')
      self.zz = zeros(nn,'d')
      self.dels = zeros((6,nn),'d')
      self.vs = zeros((6,nn),'d')
      self.ns = zeros((6,nn))
      self.parity = zeros(nn)
      self.levelxy = zeros(nn)
      self.levelz  = zeros(nn)
    elif generator is not None:
      self.ndata = len(ix)
      self.ix = ix
      self.iy = iy
      self.iz = iz
      self.xx = xx
      self.yy = yy
      self.zz = zz
      self.dels = zeros((6,self.ndata),'d')
      apply(generator,kwlist + [self.ndata,self.xx,self.yy,self.zz,
                                self.dels[0,:],self.dels[1,:],
                                self.dels[2,:],self.dels[3,:],
                                self.dels[4,:],self.dels[5,:]])
      self.setvoltages(voltage)
      self.setcondids(condid)
      self.setlevels(1,1)
    else:
      self.ndata = len(ix)
      self.ix = ix
      self.iy = iy
      self.iz = iz
      self.xx = xx
      self.yy = yy
      self.zz = zz
      self.dels = dels
      self.vs = vs
      self.ns = ns
      self.parity = parity
      self.setlevels(1,1)
   
  def setvoltages(self,voltage):
    "Routine to set appropriate voltages."
    self.vs = voltage + zeros((6,self.ndata),'d')
   
  def setcondids(self,condid):
    "Routine to setcondid condids."
    self.ns = condid + zeros((6,self.ndata))
   
  def setlevels(self,lxy,lz):
    self.levelxy = lxy + zeros(self.ndata)
    self.levelz  = lz  + zeros(self.ndata)

  def normalize(self,dx,dy,dz):
    """
Normalizes the data with respect to the grid cell sizes.
dx,dy,dz: the grid cell sizes
    """
    self.dels[:,:] = self.dels/array([dx,dx,dy,dy,dz,dz])[:,NewAxis]

  def setparity(self):
    """
Set parity. For points inside, this is set to -1. For points near the surface,
this is set to the parity of ix+iy+iz. Otherwise defaults to large integer.
This assumes that the data has already been normalized with respect to the
grid cell sizes.
    """
    self.parity = zeros(self.ndata) + 999
    fuzz = 1.e-9
    # --- A compiled routine is called for optimization
    setconductorparity(self.ndata,self.ix,self.iy,self.iz,
                       self.dels[0,:],self.dels[1,:],
                       self.dels[2,:],self.dels[3,:],
                       self.dels[4,:],self.dels[5,:],self.parity,fuzz)

  def clean(self):
    """
Removes the data which is far from any conductors. Assumes that setparity
has already been called.
    """
    ii = compress(self.parity < 2,arange(self.ndata))
    self.ix    = take(self.ix,ii)
    self.iy    = take(self.iy,ii)
    self.iz    = take(self.iz,ii)
    self.xx    = take(self.xx,ii)
    self.yy    = take(self.yy,ii)
    self.zz    = take(self.zz,ii)
    self.dels  = take(self.dels,ii,1)
    self.vs  = take(self.vs,ii,1)
    self.ns  = take(self.ns,ii,1)
    self.parity= take(self.parity,ii)
    self.ndata = len(self.ix)

  def append(self,d):
    n1 = self.ndata
    n2 = d.ndata
    if n1 + n2 > len(self.ix):
      ix = self.ix[:n1]
      iy = self.iy[:n1]
      iz = self.iz[:n1]
      xx = self.xx[:n1]
      yy = self.yy[:n1]
      zz = self.zz[:n1]
      dels = self.dels[:,:n1]
      vs = self.vs[:,:n1]
      ns = self.ns[:,:n1]
      parity = self.parity[:n1]
      levelxy = self.levelxy[:n1]
      levelz = self.levelz[:n1]

      newn = max(int(2*len(self.ix)),n1+n2)
      self.ix = zeros(newn)
      self.iy = zeros(newn)
      self.iz = zeros(newn)
      self.xx = zeros(newn,'d')
      self.yy = zeros(newn,'d')
      self.zz = zeros(newn,'d')
      self.dels = zeros((6,newn),'d')
      self.vs = zeros((6,newn),'d')
      self.ns = zeros((6,newn))
      self.parity = zeros(newn)
      self.levelxy = zeros(newn)
      self.levelz  = zeros(newn)

      self.ix[:n1] = ix
      self.iy[:n1] = iy
      self.iz[:n1] = iz
      self.xx[:n1] = xx
      self.yy[:n1] = yy
      self.zz[:n1] = zz
      self.dels[:,:n1] = dels
      self.vs[:,:n1] = vs
      self.ns[:,:n1] = ns
      self.parity[:n1] = parity
      self.levelxy[:n1] = levelxy
      self.levelz[:n1] = levelz

    self.ix[n1:n1+n2] = d.ix[:n2]
    self.iy[n1:n1+n2] = d.iy[:n2]
    self.iz[n1:n1+n2] = d.iz[:n2]
    self.xx[n1:n1+n2] = d.xx[:n2]
    self.yy[n1:n1+n2] = d.yy[:n2]
    self.zz[n1:n1+n2] = d.zz[:n2]
    self.dels[:,n1:n1+n2] = d.dels[:,:n2]
    self.vs[:,n1:n1+n2] = d.vs[:,:n2]
    self.ns[:,n1:n1+n2] = d.ns[:,:n2]
    self.parity[n1:n1+n2] = d.parity[:n2]
    self.levelxy[n1:n1+n2] = d.levelxy[:n2]
    self.levelz[n1:n1+n2] = d.levelz[:n2]
    self.ndata = n1 + n2

  def install(self):
    """
Installs the data into the WARP database
    """
    nc = f3d.ncond
    nn = sum(where(self.parity[:self.ndata] == -1,1,0))
    if nn > 0:
      if nc + nn > f3d.ncondmax:
        f3d.ncondmax = nn + nc
        gchange("Conductor3d")
      f3d.ncond = f3d.ncond + nn
      ii = compress(self.parity[:self.ndata] == -1,arange(self.ndata))
      f3d.ixcond[nc:nc+nn] = take(self.ix,ii)
      f3d.iycond[nc:nc+nn] = take(self.iy,ii)
      f3d.izcond[nc:nc+nn] = take(self.iz,ii)
      f3d.condvolt[nc:nc+nn] = take(self.vs[0,:],ii)
      f3d.condnumb[nc:nc+nn] = take(self.ns[0,:],ii)
      f3d.icondlxy[nc:nc+nn] = take(self.levelxy,ii)
      f3d.icondlz[nc:nc+nn] = take(self.levelz,ii)

    ne = f3d.necndbdy
    nn = sum(where(self.parity[:self.ndata] == 0,1,0))
    if nn > 0:
      if ne + nn > f3d.ncndmax:
        f3d.ncndmax = nn + ne
        gchange("Conductor3d")
      f3d.necndbdy = f3d.necndbdy + nn
      ii = compress(self.parity[:self.ndata] == 0,arange(self.ndata))
      f3d.iecndx[ne:ne+nn] = take(self.ix,ii)
      f3d.iecndy[ne:ne+nn] = take(self.iy,ii)
      f3d.iecndz[ne:ne+nn] = take(self.iz,ii)
      f3d.ecdelmx[ne:ne+nn] = take(self.dels[0,:],ii)
      f3d.ecdelpx[ne:ne+nn] = take(self.dels[1,:],ii)
      f3d.ecdelmy[ne:ne+nn] = take(self.dels[2,:],ii)
      f3d.ecdelpy[ne:ne+nn] = take(self.dels[3,:],ii)
      f3d.ecdelmz[ne:ne+nn] = take(self.dels[4,:],ii)
      f3d.ecdelpz[ne:ne+nn] = take(self.dels[5,:],ii)
      f3d.ecvoltmx[ne:ne+nn] = take(self.vs[0,:],ii)
      f3d.ecvoltpx[ne:ne+nn] = take(self.vs[1,:],ii)
      f3d.ecvoltmy[ne:ne+nn] = take(self.vs[2,:],ii)
      f3d.ecvoltpy[ne:ne+nn] = take(self.vs[3,:],ii)
      f3d.ecvoltmz[ne:ne+nn] = take(self.vs[4,:],ii)
      f3d.ecvoltpz[ne:ne+nn] = take(self.vs[5,:],ii)
      f3d.ecnumbmx[ne:ne+nn] = take(self.ns[0,:],ii)
      f3d.ecnumbpx[ne:ne+nn] = take(self.ns[1,:],ii)
      f3d.ecnumbmy[ne:ne+nn] = take(self.ns[2,:],ii)
      f3d.ecnumbpy[ne:ne+nn] = take(self.ns[3,:],ii)
      f3d.ecnumbmz[ne:ne+nn] = take(self.ns[4,:],ii)
      f3d.ecnumbpz[ne:ne+nn] = take(self.ns[5,:],ii)
      f3d.ecvolt[ne:ne+nn] = take(self.vs[0,:],ii)
      f3d.ecnumb[ne:ne+nn] = take(self.ns[0,:],ii)
      f3d.iecndlxy[ne:ne+nn] = take(self.levelxy,ii)
      f3d.iecndlz[ne:ne+nn] = take(self.levelz,ii)

    no = f3d.nocndbdy
    nn = sum(where(self.parity[:self.ndata] == 1,1,0))
    if nn > 0:
      if no + nn > f3d.ncndmax:
        f3d.ncndmax = nn + no
        gchange("Conductor3d")
      f3d.nocndbdy = f3d.nocndbdy + nn
      ii = compress(self.parity[:self.ndata] == 1,arange(self.ndata))
      f3d.iocndx[no:no+nn] = take(self.ix,ii)
      f3d.iocndy[no:no+nn] = take(self.iy,ii)
      f3d.iocndz[no:no+nn] = take(self.iz,ii)
      f3d.ocdelmx[no:no+nn] = take(self.dels[0,:],ii)
      f3d.ocdelpx[no:no+nn] = take(self.dels[1,:],ii)
      f3d.ocdelmy[no:no+nn] = take(self.dels[2,:],ii)
      f3d.ocdelpy[no:no+nn] = take(self.dels[3,:],ii)
      f3d.ocdelmz[no:no+nn] = take(self.dels[4,:],ii)
      f3d.ocdelpz[no:no+nn] = take(self.dels[5,:],ii)
      f3d.ocvoltmx[no:no+nn] = take(self.vs[0,:],ii)
      f3d.ocvoltpx[no:no+nn] = take(self.vs[1,:],ii)
      f3d.ocvoltmy[no:no+nn] = take(self.vs[2,:],ii)
      f3d.ocvoltpy[no:no+nn] = take(self.vs[3,:],ii)
      f3d.ocvoltmz[no:no+nn] = take(self.vs[4,:],ii)
      f3d.ocvoltpz[no:no+nn] = take(self.vs[5,:],ii)
      f3d.ocnumbmx[no:no+nn] = take(self.ns[0,:],ii)
      f3d.ocnumbpx[no:no+nn] = take(self.ns[1,:],ii)
      f3d.ocnumbmy[no:no+nn] = take(self.ns[2,:],ii)
      f3d.ocnumbpy[no:no+nn] = take(self.ns[3,:],ii)
      f3d.ocnumbmz[no:no+nn] = take(self.ns[4,:],ii)
      f3d.ocnumbpz[no:no+nn] = take(self.ns[5,:],ii)
      f3d.ocvolt[no:no+nn] = take(self.vs[0,:],ii)
      f3d.ocnumb[no:no+nn] = take(self.ns[0,:],ii)
      f3d.iocndlxy[no:no+nn] = take(self.levelxy,ii)
      f3d.iocndlz[no:no+nn] = take(self.levelz,ii)

  def __neg__(self):
    "Delta not operator."
    return Delta(self.ix,self.iy,self.iz,self.xx,self.yy,self.zz,
                 -self.dels,self.vs,self.ns)

  def __mul__(self,right):
    "'and' operator, returns maximum of distances to surfaces."
    c = less(self.dels,right.dels)
    return Delta(self.ix,self.iy,self.iz,self.xx,self.yy,self.zz,
                 choose(c,(self.dels,right.dels)),
                 choose(c,(self.vs,right.vs)),
                 choose(c,(self.ns,right.ns)))

  def __add__(self,right):
    "'or' operator, returns minimum of distances to surfaces."
    c = greater(self.dels,right.dels)
    return Delta(self.ix,self.iy,self.iz,self.xx,self.yy,self.zz,
                 choose(c,(self.dels,right.dels)),
                 choose(c,(self.vs,right.vs)),
                 choose(c,(self.ns,right.ns)))

  def __sub__(self,right):
    "'or' operator, returns minimum of distances to surfaces."
    rdels = -right.dels
    c = less(self.dels,rdels)
    return Delta(self.ix,self.iy,self.iz,self.xx,self.yy,self.zz,
                 choose(c,(self.dels,rdels)),
                 choose(c,(self.vs,right.vs)),
                 choose(c,(self.ns,right.ns)))

  def __str__(self):
    "Prints out delta"
    return repr(self.dels)+" "+repr(self.vs)+" "+repr(self.ns)

##############################################################################
##############################################################################
##############################################################################

class Grid:
  """
Class holding the grid info.
  """

  def __init__(self,xmin=None,xmax=None,ymin=None,ymax=None,
                    zmin=None,zmax=None):
    """
Creates a grid object which can generate conductor data.
Constructor arguments:
  - xmin,xmax,ymin,ymax,zmin,zmax: extent of conductors. Defaults to the
    mesh size. These only need to be set for optimization, to avoid looking
    for conductors where there are none.
Call getdata(a) to generate the conductor data. 'a' is a geometry object.
Call installdata() to install the data into the WARP database.
    """
    self.xmin = where(xmin==None,w3d.xmmin,xmin)[0]
    self.xmax = where(xmax==None,w3d.xmmax,xmax)[0]
    self.ymin = where(ymin==None,w3d.ymmin,ymin)[0]
    self.ymax = where(ymax==None,w3d.ymmax,ymax)[0]
    self.zmin = where(zmin==None,w3d.zmmin,zmin)[0]
    self.zmax = where(zmax==None,w3d.zmmax,zmax)[0]

    self.nx = w3d.nx
    self.ny = w3d.ny
    self.nz = w3d.nz
    self.nzfull = w3d.nzfull
    self.xmmin = w3d.xmmin
    self.ymmin = w3d.ymmin
    self.zmmin = w3d.zmmin
    if lparallel: self.zmmin = top.zmslmin[0]

    # --- Calculate dx, dy, and dz in case this is called before
    # --- the generate.
    self.dx = (w3d.xmmax - w3d.xmmin)/w3d.nx
    self.dy = (w3d.ymmax - w3d.ymmin)/w3d.ny
    self.dz = (w3d.zmmax - w3d.zmmin)/w3d.nz

    self.mglevels = zeros(1)
    self.mglevelsnx = zeros(20)
    self.mglevelsny = zeros(20)
    self.mglevelsnzfull = zeros(20)
    self.mglevelsiz = zeros(20)
    self.mglevelsnz = zeros(20)
    if top.fstype in [7,11]:
      getmglevels(self.nx,self.ny,self.nz,self.nzfull,self.dx,self.dy,self.dz,
                  f3d.mgminlevelxy,f3d.mgminlevelz,
                  self.mglevels,self.mglevelsnx,self.mglevelsny,
                  self.mglevelsnzfull,self.mglevelsiz,self.mglevelsnz)
      self.mglevels = self.mglevels[0]
    else:
      self.mglevels = 1
      self.mglevelsnx[0] = self.nx
      self.mglevelsny[0] = self.ny
      self.mglevelsnzfull[0] = self.nzfull
      self.mglevelsiz[0] = top.izfsslave[me]
      self.mglevelsnz[0] = self.nz
    self.mglevelsnx = self.mglevelsnx[:self.mglevels]
    self.mglevelsny = self.mglevelsny[:self.mglevels]
    self.mglevelsnzfull = self.mglevelsnzfull[:self.mglevels]
    self.mglevelsiz = self.mglevelsiz[:self.mglevels]
    self.mglevelsnz = self.mglevelsnz[:self.mglevels]

  def getdata(self,a):
    """
Given an Assembly, accumulate the appropriate data to represent that
Assembly on this grid.
    """
    starttime = wtime()
    tt2 = zeros(8,'d')
    self.dall = Delta()
    for i in range(self.mglevels):
      tt1 = wtime()
      nx = self.mglevelsnx[i]
      ny = self.mglevelsny[i]
      nzfull = self.mglevelsnzfull[i]
      iz = self.mglevelsiz[i]
      nz = self.mglevelsnz[i]

      dx = self.dx*self.nx/nx
      dy = self.dy*self.ny/ny
      dz = self.dz*self.nzfull/nzfull

      zmmin = self.zmmin + iz*dz

      xmesh = self.xmmin + dx*arange(nx+1)
      ymesh = self.ymmin + dy*arange(ny+1)
      zmesh =      zmmin + dz*arange(nz+1)
      xmesh = compress(logical_and(self.xmin-dx <= xmesh,
                                                   xmesh <= self.xmax+dx),xmesh)
      ymesh = compress(logical_and(self.ymin-dy <= ymesh,
                                                   ymesh <= self.ymax+dy),ymesh)
      zmesh = compress(logical_and(self.zmin-dz <= zmesh,
                                                   zmesh <= self.zmax+dz),zmesh)
      x = ravel(xmesh[:,NewAxis]*ones(len(ymesh)))
      y = ravel(ymesh*ones(len(xmesh))[:,NewAxis])
      z = zeros(len(xmesh)*len(ymesh),'d')
      ix = nint((x - self.xmmin)/dx)
      iy = nint((y - self.ymmin)/dy)
      iz = zeros(len(xmesh)*len(ymesh))
      tt2[0] = tt2[0] + wtime() - tt1
      if len(x) == 0: continue
      for zz in zmesh:
        tt1 = wtime()
        z[:] = zz
        iz[:] = nint((zz - zmmin)/dz)
        tt2[1] = tt2[1] + wtime() - tt1
        tt1 = wtime()
        d = a.distance(ix,iy,iz,x,y,z)
        tt2[2] = tt2[2] + wtime() - tt1
        tt1 = wtime()
        d.normalize(dx,dy,dz)
        tt2[3] = tt2[3] + wtime() - tt1
        tt1 = wtime()
        d.setparity()
        tt2[4] = tt2[4] + wtime() - tt1
        tt1 = wtime()
        d.clean()
        tt2[5] = tt2[5] + wtime() - tt1
        tt1 = wtime()
        d.setlevels(self.nx/nx,self.nz/nz)
        tt2[6] = tt2[6] + wtime() - tt1
        tt1 = wtime()
        self.dall.append(d)
        tt2[7] = tt2[7] + wtime() - tt1
    endtime = wtime()
    self.generatetime = endtime - starttime
    print tt2

  def installdata(self):
    """
Installs the conductor data into the fortran database
    """
    self.dall.install()

##############################################################################
##############################################################################
##############################################################################

#============================================================================
class Box(Assembly):
  """
Box class
  """

  xsize = 0.
  ysize = 0.
  zsize = 0.

  def __init__(self,xsize,ysize,zsize,voltage=0.,xcent=0.,ycent=0.,zcent=0.,
                    condid=0):
    Assembly.__init__(self,voltage,xcent,ycent,zcent)
    self.xsize = xsize
    self.ysize = ysize
    self.zsize = zsize
    self.condid = condid

  def distance(self,ix,iy,iz,xx,yy,zz):
    result = Delta(ix,iy,iz,xx,yy,zz,voltage=self.voltage,condid=self.condid,
                   generator=f3d.boxconductorf,
                   kwlist=[self.xsize,self.ysize,self.zsize,
                           self.xcent,self.ycent,self.zcent])
    return result

#============================================================================
class Cylinder(Assembly):
  """
Cylinder class
  """

  def __init__(self,radius,length,theta=0.,phi=0.,
                    voltage=0.,xcent=0.,ycent=0.,zcent=0.,condid=0):
    Assembly.__init__(self,voltage,xcent,ycent,zcent)
    self.radius = radius
    self.length = length
    self.theta = theta
    self.phi = phi
    self.condid = condid

  def distance(self,ix,iy,iz,xx,yy,zz):
    result = Delta(ix,iy,iz,xx,yy,zz,voltage=self.voltage,condid=self.condid,
                   generator=f3d.cylinderconductorf,
                   kwlist=[self.radius,self.length,self.theta,self.phi,
                           self.xcent,self.ycent,self.zcent])
    return result

class ZCylinder(Assembly):
  """
ZCylinder
  """

  def __init__(self,radius,length,voltage=0.,xcent=0.,ycent=0.,zcent=0.,
                    condid=0):
    Assembly.__init__(self,voltage,xcent,ycent,zcent)
    self.radius = radius
    self.length = length
    self.condid = condid

  def distance(self,ix,iy,iz,xx,yy,zz):
    result = Delta(ix,iy,iz,xx,yy,zz,voltage=self.voltage,condid=self.condid,
                   generator=f3d.zcylinderconductorf,
                   kwlist=[self.radius,self.length,
                           self.xcent,self.ycent,self.zcent])
    return result

class YCylinder(Assembly):
  """
YCylinder
  """

  def __init__(self,radius,length,voltage=0.,xcent=0.,ycent=0.,zcent=0.,
                    condid=0):
    Assembly.__init__(self,voltage,xcent,ycent,zcent)
    self.radius = radius
    self.length = length
    self.condid = condid

  def distance(self,ix,iy,iz,xx,yy,zz):
    result = Delta(ix,iy,iz,xx,yy,zz,voltage=self.voltage,condid=self.condid,
                   generator=f3d.ycylinderconductorf,
                   kwlist=[self.radius,self.length,
                           self.xcent,self.ycent,self.zcent])
    return result

class XCylinder(Assembly):
  """
XCylinder
  """

  def __init__(self,radius,length,voltage=0.,xcent=0.,ycent=0.,zcent=0.,
                    condid=0):
    Assembly.__init__(self,voltage,xcent,ycent,zcent)
    self.radius = radius
    self.length = length
    self.condid = condid

  def distance(self,ix,iy,iz,xx,yy,zz):
    result = Delta(ix,iy,iz,xx,yy,zz,voltage=self.voltage,condid=self.condid,
                   generator=f3d.xcylinderconductorf,
                   kwlist=[self.radius,self.length,
                           self.xcent,self.ycent,self.zcent])
    return result

#============================================================================
class Spheroid(Assembly):
  """
Spheroid class
  """

  def __init__(self,xradius,yradius,zradius,
                    voltage=0.,xcent=0.,ycent=0.,zcent=0.,condid=0):
    Assembly.__init__(self,voltage,xcent,ycent,zcent)
    self.xradius = xradius
    self.yradius = yradius
    self.zradius = zradius
    self.condid = condid

class Sphere(Spheroid):
  """
Sphere
  """

  def __init__(self,radius,voltage=0.,xcent=0.,ycent=0.,zcent=0.,
                    condid=0):
    Spheroid.__init__(self,radius,radius,radius,voltage,xcent,ycent,zcent)
    self.radius = radius

  def distance(self,ix,iy,iz,xx,yy,zz):
    result = Delta(ix,iy,iz,xx,yy,zz,voltage=self.voltage,condid=self.condid,
                   generator=f3d.sphereconductorf,
                   kwlist=[self.radius,self.xcent,self.ycent,self.zcent])
    return result

#============================================================================
class ZCone(Assembly):
  """
Cone
  """

  def __init__(self,r_zmin,r_zmax,length,voltage=0.,xcent=0.,ycent=0.,zcent=0.,
                    condid=0):
    Assembly.__init__(self,voltage,xcent,ycent,zcent)
    self.r_zmin = r_zmin
    self.r_zmax = r_zmax
    self.length = length
    self.condid = condid

  def distance(self,ix,iy,iz,xx,yy,zz):
    result = Delta(ix,iy,iz,xx,yy,zz,voltage=self.voltage,condid=self.condid,
                   generator=f3d.zconeconductorf,
                   kwlist=[self.r_zmin,self.r_zmax,self.length,
                           self.xcent,self.ycent,self.zcent])
    return result

#============================================================================
class ZTorus(Assembly):
  """
Cone
  """

  def __init__(self,r1,r2,voltage=0.,xcent=0.,ycent=0.,zcent=0.,condid=0):
    Assembly.__init__(self,voltage,xcent,ycent,zcent)
    self.r1 = r1
    self.r2 = r2
    self.condid = condid

  def distance(self,ix,iy,iz,xx,yy,zz):
    result = Delta(ix,iy,iz,xx,yy,zz,voltage=self.voltage,condid=self.condid,
                   generator=f3d.ztorusconductorf,
                   kwlist=[self.r1,self.r2,
                           self.xcent,self.ycent,self.zcent])
    return result

#============================================================================
#############################################################################

