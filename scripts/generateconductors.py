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

generateconductorsversion = "$Id: generateconductors.py,v 1.1 2002/06/04 23:30:35 dave Exp $"
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

_n = 10000
class Delta:
  """
Class to hold the set of distances in each of the six directions.
Distances have the sign of the outward normal surface vector, i.e.
distances to outside the surface are positive, inside negative.
  """

  def __init__(self,ix=zeros(_n),iy=zeros(_n),iz=zeros(_n),
                    xx=zeros(_n,'d'),yy=zeros(_n,'d'),zz=zeros(_n,'d'),
                    delmx=zeros(_n,'d'),delpx=zeros(_n,'d'),
                    delmy=zeros(_n,'d'),delpy=zeros(_n,'d'),
                    delmz=zeros(_n,'d'),delpz=zeros(_n,'d'),
                    vmx=zeros(_n,'d'),vpx=zeros(_n,'d'),
                    vmy=zeros(_n,'d'),vpy=zeros(_n,'d'),
                    vmz=zeros(_n,'d'),vpz=zeros(_n,'d'),
                    parity=zeros(_n),
                    voltage=0.,generator=None,kwlist=[]):
    if generator is None:
      self.ndata = len(ix)
      self.ix = ix
      self.iy = iy
      self.iz = iz
      self.xx = xx
      self.yy = yy
      self.zz = zz
      self.delmx = delmx
      self.delpx = delpx
      self.delmy = delmy
      self.delpy = delpy
      self.delmz = delmz
      self.delpz = delpz
      self.vmx = vmx
      self.vpx = vmx
      self.vmy = vmy
      self.vpy = vmy
      self.vmz = vmz
      self.vpz = vmz
      self.parity = parity
    else:
      self.ndata = len(ix)
      self.ix = ix
      self.iy = iy
      self.iz = iz
      self.xx = xx
      self.yy = yy
      self.zz = zz
      self.delmx = zeros(self.ndata,'d')
      self.delpx = zeros(self.ndata,'d')
      self.delmy = zeros(self.ndata,'d')
      self.delpy = zeros(self.ndata,'d')
      self.delmz = zeros(self.ndata,'d')
      self.delpz = zeros(self.ndata,'d')
      apply(generator,kwlist + [self.ndata,self.xx,self.yy,self.zz,
                                self.delmx,self.delpx,
                                self.delmy,self.delpy,
                                self.delmz,self.delpz])
      self.setvoltages(voltage)
   
  def setvoltages(self,voltage):
    "Routine to set appropriate voltages."
    self.vmx = voltage + zeros(self.ndata,'d')
    self.vpx = voltage + zeros(self.ndata,'d')
    self.vmy = voltage + zeros(self.ndata,'d')
    self.vpy = voltage + zeros(self.ndata,'d')
    self.vmz = voltage + zeros(self.ndata,'d')
    self.vpz = voltage + zeros(self.ndata,'d')
   
  def normalize(self,dx,dy,dz):
    """
Normalizes the data with respect to the grid cell sizes.
dx,dy,dz: the grid cell sizes
    """
    self.delmx[:] = self.delmx/dx
    self.delpx[:] = self.delpx/dx
    self.delmy[:] = self.delmy/dy
    self.delpy[:] = self.delpy/dy
    self.delmz[:] = self.delmz/dz
    self.delpz[:] = self.delpz/dz

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
                       self.delmx,self.delpx,self.delmy,self.delpy,
                       self.delmz,self.delpz,self.parity,fuzz)

  # iparity = (self.ix+self.iy+self.iz)%2
  # self.parity[:] = where(self.delmx < 1.-fuzz,iparity,self.parity)
  # self.parity[:] = where(self.delpx < 1.-fuzz,iparity,self.parity)
  # self.parity[:] = where(self.delmy < 1.-fuzz,iparity,self.parity)
  # self.parity[:] = where(self.delpy < 1.-fuzz,iparity,self.parity)
  # self.parity[:] = where(self.delmz < 1.-fuzz,iparity,self.parity)
  # self.parity[:] = where(self.delpz < 1.-fuzz,iparity,self.parity)

  # self.parity[:] = where(self.delmx <= 0.+fuzz,-1,self.parity)
  # self.parity[:] = where(self.delpx <= 0.+fuzz,-1,self.parity)
  # self.parity[:] = where(self.delmy <= 0.+fuzz,-1,self.parity)
  # self.parity[:] = where(self.delpy <= 0.+fuzz,-1,self.parity)
  # self.parity[:] = where(self.delmz <= 0.+fuzz,-1,self.parity)
  # self.parity[:] = where(self.delpz <= 0.+fuzz,-1,self.parity)

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
    self.delmx = take(self.delmx,ii)
    self.delpx = take(self.delpx,ii)
    self.delmy = take(self.delmy,ii)
    self.delpy = take(self.delpy,ii)
    self.delmz = take(self.delmz,ii)
    self.delpz = take(self.delpz,ii)
    self.vmx   = take(self.vmx,ii)
    self.vpx   = take(self.vpx,ii)
    self.vmy   = take(self.vmy,ii)
    self.vpy   = take(self.vpy,ii)
    self.vmz   = take(self.vmz,ii)
    self.vpz   = take(self.vpz,ii)
    self.parity= take(self.parity,ii)
    self.ndata = len(self.ix)

  def append(self,d):
    n1 = self.ndata
    n2 = len(d.ix)
    if n1 + n2 > len(self.ix):
      ix = self.ix[:n1]
      iy = self.iy[:n1]
      iz = self.iz[:n1]
      xx = self.xx[:n1]
      yy = self.yy[:n1]
      zz = self.zz[:n1]
      delmx = self.delmx[:n1]
      delpx = self.delpx[:n1]
      delmy = self.delmy[:n1]
      delpy = self.delpy[:n1]
      delmz = self.delmz[:n1]
      delpz = self.delpz[:n1]
      vmx = self.vmx[:n1]
      vpx = self.vpx[:n1]
      vmy = self.vmy[:n1]
      vpy = self.vpy[:n1]
      vmz = self.vmz[:n1]
      vpz = self.vpz[:n1]
      parity = self.parity[:n1]

      newn = max(int(2*len(self.ix)),n1+n2)
      self.ix = zeros(newn)
      self.iy = zeros(newn)
      self.iz = zeros(newn)
      self.xx = zeros(newn,'d')
      self.yy = zeros(newn,'d')
      self.zz = zeros(newn,'d')
      self.delmx = zeros(newn,'d')
      self.delpx = zeros(newn,'d')
      self.delmy = zeros(newn,'d')
      self.delpy = zeros(newn,'d')
      self.delmz = zeros(newn,'d')
      self.delpz = zeros(newn,'d')
      self.vmx = zeros(newn,'d')
      self.vpx = zeros(newn,'d')
      self.vmy = zeros(newn,'d')
      self.vpy = zeros(newn,'d')
      self.vmz = zeros(newn,'d')
      self.vpz = zeros(newn,'d')
      self.parity = zeros(newn)

      self.ix[:n1] = ix
      self.iy[:n1] = iy
      self.iz[:n1] = iz
      self.xx[:n1] = xx
      self.yy[:n1] = yy
      self.zz[:n1] = zz
      self.delmx[:n1] = delmx
      self.delpx[:n1] = delpx
      self.delmy[:n1] = delmy
      self.delpy[:n1] = delpy
      self.delmz[:n1] = delmz
      self.delpz[:n1] = delpz
      self.vmx[:n1] = vmx
      self.vpx[:n1] = vpx
      self.vmy[:n1] = vmy
      self.vpy[:n1] = vpy
      self.vmz[:n1] = vmz
      self.vpz[:n1] = vpz
      self.parity[:n1] = parity

    self.ix[n1:n1+n2] = d.ix
    self.iy[n1:n1+n2] = d.iy
    self.iz[n1:n1+n2] = d.iz
    self.xx[n1:n1+n2] = d.xx
    self.yy[n1:n1+n2] = d.yy
    self.zz[n1:n1+n2] = d.zz
    self.delmx[n1:n1+n2] = d.delmx
    self.delpx[n1:n1+n2] = d.delpx
    self.delmy[n1:n1+n2] = d.delmy
    self.delpy[n1:n1+n2] = d.delpy
    self.delmz[n1:n1+n2] = d.delmz
    self.delpz[n1:n1+n2] = d.delpz
    self.vmx[n1:n1+n2] = d.vmx
    self.vpx[n1:n1+n2] = d.vpx
    self.vmy[n1:n1+n2] = d.vmy
    self.vpy[n1:n1+n2] = d.vpy
    self.vmz[n1:n1+n2] = d.vmz
    self.vpz[n1:n1+n2] = d.vpz
    self.parity[n1:n1+n2] = d.parity
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
        gchange("PSOR3d")
        gchange("MultigridConductor3d")
      f3d.ncond = f3d.ncond + nn
      ii = compress(self.parity[:self.ndata] == -1,arange(self.ndata))
      f3d.ixcond[nc:nc+nn] = take(self.ix,ii)
      f3d.iycond[nc:nc+nn] = take(self.iy,ii)
      f3d.izcond[nc:nc+nn] = take(self.iz,ii)
      f3d.condvolt[nc:nc+nn] = take(self.vmx,ii)
      f3d.icondlxy[nc:nc+nn] = 1
      f3d.icondlz[nc:nc+nn] = 1

    ne = f3d.necndbdy
    nn = sum(where(self.parity[:self.ndata] == 0,1,0))
    if nn > 0:
      if ne + nn > f3d.ncndmax:
        f3d.ncndmax = nn + ne
        gchange("PSOR3d")
        gchange("MultigridConductor3d")
      f3d.necndbdy = f3d.necndbdy + nn
      ii = compress(self.parity[:self.ndata] == 0,arange(self.ndata))
      f3d.iecndx[ne:ne+nn] = take(self.ix,ii)
      f3d.iecndy[ne:ne+nn] = take(self.iy,ii)
      f3d.iecndz[ne:ne+nn] = take(self.iz,ii)
      f3d.ecdelmx[ne:ne+nn] = take(self.delmx,ii)
      f3d.ecdelpx[ne:ne+nn] = take(self.delpx,ii)
      f3d.ecdelmy[ne:ne+nn] = take(self.delmy,ii)
      f3d.ecdelpy[ne:ne+nn] = take(self.delpy,ii)
      f3d.ecdelmz[ne:ne+nn] = take(self.delmz,ii)
      f3d.ecdelpz[ne:ne+nn] = take(self.delpz,ii)
      f3d.ecvoltmx[ne:ne+nn] = take(self.vmx,ii)
      f3d.ecvoltpx[ne:ne+nn] = take(self.vpx,ii)
      f3d.ecvoltmy[ne:ne+nn] = take(self.vmy,ii)
      f3d.ecvoltpy[ne:ne+nn] = take(self.vpy,ii)
      f3d.ecvoltmz[ne:ne+nn] = take(self.vmz,ii)
      f3d.ecvoltpz[ne:ne+nn] = take(self.vpz,ii)
      f3d.ecvolt[ne:ne+nn] = take(self.vmx,ii)
      f3d.iecndlxy[ne:ne+nn] = 1
      f3d.iecndlz[ne:ne+nn] = 1

    no = f3d.nocndbdy
    nn = sum(where(self.parity[:self.ndata] == 1,1,0))
    if nn > 0:
      if no + nn > f3d.ncndmax:
        f3d.ncndmax = nn + no
        gchange("PSOR3d")
        gchange("MultigridConductor3d")
      f3d.nocndbdy = f3d.nocndbdy + nn
      ii = compress(self.parity[:self.ndata] == 1,arange(self.ndata))
      f3d.iocndx[no:no+nn] = take(self.ix,ii)
      f3d.iocndy[no:no+nn] = take(self.iy,ii)
      f3d.iocndz[no:no+nn] = take(self.iz,ii)
      f3d.ocdelmx[no:no+nn] = take(self.delmx,ii)
      f3d.ocdelpx[no:no+nn] = take(self.delpx,ii)
      f3d.ocdelmy[no:no+nn] = take(self.delmy,ii)
      f3d.ocdelpy[no:no+nn] = take(self.delpy,ii)
      f3d.ocdelmz[no:no+nn] = take(self.delmz,ii)
      f3d.ocdelpz[no:no+nn] = take(self.delpz,ii)
      f3d.ocvoltmx[no:no+nn] = take(self.vmx,ii)
      f3d.ocvoltpx[no:no+nn] = take(self.vpx,ii)
      f3d.ocvoltmy[no:no+nn] = take(self.vmy,ii)
      f3d.ocvoltpy[no:no+nn] = take(self.vpy,ii)
      f3d.ocvoltmz[no:no+nn] = take(self.vmz,ii)
      f3d.ocvoltpz[no:no+nn] = take(self.vpz,ii)
      f3d.ocvolt[no:no+nn] = take(self.vmx,ii)
      f3d.iocndlxy[no:no+nn] = 1
      f3d.iocndlz[no:no+nn] = 1

  def __neg__(self):
    "Delta not operator."
    return Delta(self.ix,self.iy,self.iz,self.xx,self.yy,self.zz,
                 -self.delmx,-self.delpx,
                 -self.delmy,-self.delpy,
                 -self.delmz,-self.delpz,
                 self.vmx,self.vpx,
                 self.vmy,self.vpy,
                 self.vmz,self.vpz)

  def __mul__(self,right):
    "'and' operator, returns maximum of distances to surfaces."
    return Delta(self.ix,self.iy,self.iz,self.xx,self.yy,self.zz,
                 maximum(self.delmx,right.delmx),
                 maximum(self.delpx,right.delpx),
                 maximum(self.delmy,right.delmy),
                 maximum(self.delpy,right.delpy),
                 maximum(self.delmz,right.delmz),
                 maximum(self.delpz,right.delpz),
                 where(self.delmx >= right.delmx,self.vmx,right.vmx),
                 where(self.delpx >= right.delpx,self.vpx,right.vpx),
                 where(self.delmy >= right.delmy,self.vmy,right.vmy),
                 where(self.delpy >= right.delpy,self.vpy,right.vpy),
                 where(self.delmz >= right.delmz,self.vmz,right.vmz),
                 where(self.delpz >= right.delpz,self.vpz,right.vpz))

  def __add__(self,right):
    "'or' operator, returns minimum of distances to surfaces."
    return Delta(self.ix,self.iy,self.iz,self.xx,self.yy,self.zz,
                 minimum(self.delmx,right.delmx),
                 minimum(self.delpx,right.delpx),
                 minimum(self.delmy,right.delmy),
                 minimum(self.delpy,right.delpy),
                 minimum(self.delmz,right.delmz),
                 minimum(self.delpz,right.delpz),
                 where(self.delmx <= right.delmx,self.vmx,right.vmx),
                 where(self.delpx <= right.delpx,self.vpx,right.vpx),
                 where(self.delmy <= right.delmy,self.vmy,right.vmy),
                 where(self.delpy <= right.delpy,self.vpy,right.vpy),
                 where(self.delmz <= right.delmz,self.vmz,right.vmz),
                 where(self.delpz <= right.delpz,self.vpz,right.vpz))

  def __sub__(self,right):
    "'or' operator, returns minimum of distances to surfaces."
    return Delta(self.ix,self.iy,self.iz,self.xx,self.yy,self.zz,
                 maximum(self.delmx,-right.delmx),
                 maximum(self.delpx,-right.delpx),
                 maximum(self.delmy,-right.delmy),
                 maximum(self.delpy,-right.delpy),
                 maximum(self.delmz,-right.delmz),
                 maximum(self.delpz,-right.delpz),
                 where(self.delmx >= -right.delmx,self.vmx,right.vmx),
                 where(self.delpx >= -right.delpx,self.vpx,right.vpx),
                 where(self.delmy >= -right.delmy,self.vmy,right.vmy),
                 where(self.delpy >= -right.delpy,self.vpy,right.vpy),
                 where(self.delmz >= -right.delmz,self.vmz,right.vmz),
                 where(self.delpz >= -right.delpz,self.vpz,right.vpz))

  def __str__(self):
    "Prints out delta"
    return repr(self.delmx)+" "+repr(self.delpx)+" "+ \
           repr(self.delmy)+" "+repr(self.delpy)+" "+ \
           repr(self.delmz)+" "+repr(self.delpz)+"\n"+ \
           repr(self.vmx)+" "+repr(self.vpx)+" "+ \
           repr(self.vmy)+" "+repr(self.vpy)+" "+ \
           repr(self.vmz)+" "+repr(self.vpz)

##############################################################################
##############################################################################
##############################################################################

class Grid:
  """
Class holding the grid info.
  """

  nx = 1
  ny = 1
  nz = 1
  xmin = 0.
  xmax = 0.
  ymin = 0.
  ymax = 0.
  zmin = 0.
  zmax = 0.
  dx = 0.
  dy = 0.
  dz = 0.

  def __init__(self,nx=None,ny=None,nz=None,
               xmin=None,xmax=None,ymin=None,ymax=None,zmin=None,zmax=None):
    self.nx = where(nx==None,w3d.nx,nx)[0]
    self.ny = where(ny==None,w3d.ny,ny)[0]
    self.nz = where(nz==None,w3d.nz,nz)[0]
    self.xmin = where(xmin==None,w3d.xmmin,xmin)[0]
    self.xmax = where(xmax==None,w3d.xmmax,xmax)[0]
    self.ymin = where(ymin==None,w3d.ymmin,ymin)[0]
    self.ymax = where(ymax==None,w3d.ymmax,ymax)[0]
    self.zmin = where(zmin==None,w3d.zmmin,zmin)[0]
    self.zmax = where(zmax==None,w3d.zmmax,zmax)[0]
    self.dx = (self.xmax - self.xmin)/float(self.nx)
    self.dy = (self.ymax - self.ymin)/float(self.ny)
    self.dz = (self.zmax - self.zmin)/float(self.nz)
    self.xmesh = self.xmin + self.dx*arange(self.nx+1)
    self.ymesh = self.ymin + self.dy*arange(self.ny+1)
    self.zmesh = self.zmin + self.dz*arange(self.nz+1)

  def getdata(self,a):
    """
Given an Assembly, accumulate the appropriate data to represent that
Assembly on this grid.
    """
    self.dall = Delta()
    ix = ravel(arange(1+self.nx)[:,NewAxis]*ones(1+self.ny))
    iy = ravel(arange(1+self.ny)*ones(1+self.nx)[:,NewAxis])
    iz = zeros((1+self.nx)*(1+self.ny))
    x = self.xmin + ix*self.dx
    y = self.ymin + iy*self.dy
    z = zeros((1+self.nx)*(1+self.ny),'d')
    ttt = zeros(7,'d')
    stt = zeros(6,'d')
    for jz in range(0,self.nz+1):
      iz[:] = jz
      ttt[0] = wtime()
      z[:] = self.zmin + jz*self.dz
      ttt[1] = wtime()
      d = a.distance(ix,iy,iz,x,y,z)
      ttt[2] = wtime()
      d.normalize(self.dx,self.dy,self.dz)
      ttt[3] = wtime()
      d.setparity()
      ttt[4] = wtime()
      d.clean()
      ttt[5] = wtime()
      self.dall.append(d)
      ttt[6] = wtime()
      stt[:] = stt + ttt[1:] - ttt[:-1]

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

  def __init__(self,xsize,ysize,zsize,voltage=0.,xcent=0.,ycent=0.,zcent=0.):
    Assembly.__init__(self,voltage,xcent,ycent,zcent)
    self.xsize = xsize
    self.ysize = ysize
    self.zsize = zsize

  def distance(self,ix,iy,iz,xx,yy,zz):
    result = Delta(ix,iy,iz,xx,yy,zz,self.voltage,
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
                    voltage=0.,xcent=0.,ycent=0.,zcent=0.):
    Assembly.__init__(self,voltage,xcent,ycent,zcent)
    self.radius = radius
    self.length = length
    self.theta = theta
    self.phi = phi

  def distance(self,ix,iy,iz,xx,yy,zz):
    result = Delta(ix,iy,iz,xx,yy,zz,self.voltage,
                   generator=f3d.cylinderconductorf,
                   kwlist=[self.radius,self.length,self.theta,self.phi,
                           self.xcent,self.ycent,self.zcent])
    return result

class ZCylinder(Assembly):
  """
ZCylinder
  """

  def __init__(self,radius,length,voltage=0.,xcent=0.,ycent=0.,zcent=0.):
    Assembly.__init__(self,voltage,xcent,ycent,zcent)
    self.radius = radius
    self.length = length

  def distance(self,ix,iy,iz,xx,yy,zz):
    result = Delta(ix,iy,iz,xx,yy,zz,self.voltage,
                   generator=f3d.zcylinderconductorf,
                   kwlist=[self.radius,self.length,
                           self.xcent,self.ycent,self.zcent])
    return result

class YCylinder(Assembly):
  """
YCylinder
  """

  def __init__(self,radius,length,voltage=0.,xcent=0.,ycent=0.,zcent=0.):
    Assembly.__init__(self,voltage,xcent,ycent,zcent)
    self.radius = radius
    self.length = length

  def distance(self,ix,iy,iz,xx,yy,zz):
    result = Delta(ix,iy,iz,xx,yy,zz,self.voltage,
                   generator=f3d.ycylinderconductorf,
                   kwlist=[self.radius,self.length,
                           self.xcent,self.ycent,self.zcent])
    return result

class XCylinder(Assembly):
  """
XCylinder
  """

  def __init__(self,radius,length,voltage=0.,xcent=0.,ycent=0.,zcent=0.):
    Assembly.__init__(self,voltage,xcent,ycent,zcent)
    self.radius = radius
    self.length = length

  def distance(self,ix,iy,iz,xx,yy,zz):
    result = Delta(ix,iy,iz,xx,yy,zz,self.voltage,
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
                    voltage=0.,xcent=0.,ycent=0.,zcent=0.):
    Assembly.__init__(self,voltage,xcent,ycent,zcent)
    self.xradius = xradius
    self.yradius = yradius
    self.zradius = zradius

class Sphere(Spheroid):
  """
Sphere
  """

  def __init__(self,radius,voltage=0.,xcent=0.,ycent=0.,zcent=0.):
    Spheroid.__init__(self,radius,radius,radius,voltage,xcent,ycent,zcent)
    self.radius = radius

  def distance(self,ix,iy,iz,xx,yy,zz):
    result = Delta(ix,iy,iz,xx,yy,zz,self.voltage,
                   generator=f3d.sphereconductorf,
                   kwlist=[self.radius,self.xcent,self.ycent,self.zcent])
    return result

#============================================================================
class ZCone(Assembly):
  """
Cone
  """

  def __init__(self,r_zmin,r_zmax,length,voltage=0.,xcent=0.,ycent=0.,zcent=0.):
    Assembly.__init__(self,voltage,xcent,ycent,zcent)
    self.r_zmin = r_zmin
    self.r_zmax = r_zmax
    self.length = length

  def distance(self,ix,iy,iz,xx,yy,zz):
    result = Delta(ix,iy,iz,xx,yy,zz,self.voltage,
                   generator=f3d.zconeconductorf,
                   kwlist=[self.r_zmin,self.r_zmax,self.length,
                           self.xcent,self.ycent,self.zcent])
    return result

#============================================================================
class ZTorus(Assembly):
  """
Cone
  """

  def __init__(self,r1,r2,voltage=0.,xcent=0.,ycent=0.,zcent=0.):
    Assembly.__init__(self,voltage,xcent,ycent,zcent)
    self.r1 = r1
    self.r2 = r2

  def distance(self,ix,iy,iz,xx,yy,zz):
    result = Delta(ix,iy,iz,xx,yy,zz,self.voltage,
                   generator=f3d.ztorusconductorf,
                   kwlist=[self.r1,self.r2,
                           self.xcent,self.ycent,self.zcent])
    return result

#============================================================================
#############################################################################

