"""
This module contains classes for generating the conductor data from a
combination of simple geometrical elements.
The following elements are defined:

Plane(z0=0.,zsign=1,theta=0.,phi=0.,...)
Box(xsize,ysize,zsize,...)
Cylinder(radius,length,theta=0.,phi=0.,...)
ZCylinder(radius,length,...)
ZCylinderOut(radius,length,...)
ZRoundedCylinder(radius,length,radius2,...)
ZRoundedCylinderOut(radius,length,radius2,...)
YCylinder(radius,length,...)
XCylinder(radius,length,...)
Sphere(radius,...)
Cone(r_zmin,r_zmax,length,theta=0.,phi=0.,...)
ZCone(r_zmin,r_zmax,length,...)
ZConeOut(r_zmin,r_zmax,length,...)
ConeSlope(slope,intercept,length,theta=0.,phi=0.,...)
ZConeSlope(slope,intercept,length,...)
ZConeOutSlope(slope,intercept,length,...)
ZTorus(r1,r2,...)
Beamletplate(za,zb,z0,thickness,...)
ZSrfrvOut(rofzfunc,zmin,zmax,rmax,...)
ZSrfrvIn(rofzfunc,zmin,zmax,rmin,...)
ZSrfrvInOut(rminofz,rmaxofz,zmin,zmax,...)

Note that all take the following additional arguments:
voltage=0.,xcent=0.,ycent=0.,zcent=0.,condid=1

installconductors(a): generates the data needed for the fieldsolve

A set of classes for generating and manipulating surfaces of revolution
using lines and arcs primitives (SRFRVLA) are available:
SRFRVLAfromfile(filename,voltages,condids,zshifts=None,rshifts=None,install=1)
SRFRVLAsystem(SRFRVLAconds,install=1)
SRFRVLAcond(name,parts,voltage,condid,install=1)
SRFRVLApart(name,data)
SRFRVLA_circle(name,c,r)
SRFRVLA_rectangle(name,c,l,h)
SRFRVLA_rnd_rectangle(name,c,l,h,r)
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
import operator
if not lparallel: import VPythonobjects
from string import *

generateconductorsversion = "$Id: generateconductors.py,v 1.40 2003/12/03 01:02:01 dave Exp $"
def generateconductors_doc():
  import generateconductors
  print generateconductors.__doc__


# --- Temporary hack to allow use with older versions of source
try:
  ycylinderconductord
except:
  def ycylinderconductord(rad,length,xcent,ycent,zcent,
                          n,x,y,z,delmx,delpx,delmy,delpy,delmz,delpz,fuzz):
    pass
  def xcylinderconductord(rad,length,xcent,ycent,zcent,
                          n,x,y,z,delmx,delpx,delmy,delpy,delmz,delpz,fuzz):
    pass

##############################################################################
def installconductors(a,xmin=None,xmax=None,ymin=None,ymax=None,
                        zmin=None,zmax=None,dfill=top.largepos,
                        installrz=1,gridmode=1):
  """
Installs the given conductors.
  - a: the assembly of conductors
  - xmin,xmax,ymin,ymax,zmin,zmax: extent of conductors. Defaults to the
    mesh size. These can be set for optimization, to avoid looking
    for conductors where there are none. Also, they can be used crop a
    conductor
  - dfill=largepos: points at a depth in the conductor greater than dfill
                    are skipped.
  """
  # First, create a grid object
  g = Grid(xmin,xmax,ymin,ymax,zmin,zmax)
  # Generate the conductor data
  g.getdata(a,dfill)
  # Then install it
  g.installdata(installrz,gridmode)
  
##############################################################################
##############################################################################
##############################################################################
class Assembly:
  """
Class to hold assemblies of conductors.  Base class of all conductors.
Should never be directly created by the user.
 - v=0.: voltage on conductor
 - x,y,z=0.,0.,0: center of conductor
 - condid=1: conductor identification number
 - kwlist=[]: list of string names of variable describing conductor
 - generatorf=None: function which generates the distances between the points
                    and the conductors along the axis.
 - generatord=None: function which generates the smallest distance between the
                    points and the conductor surface.
  """

  voltage = 0.
  xcent = 0.
  ycent = 0.
  zcent = 0.

  def __init__(self,v=0.,x=0.,y=0.,z=0.,condid=1,kwlist=[],
                    generatorf=None,generatord=None,generatori=None):
    self.voltage = v
    self.xcent = x
    self.ycent = y
    self.zcent = z
    self.condid = condid
    self.kwlist = kwlist
    self.generatorf = generatorf
    self.generatord = generatord
    self.generatori = generatori

  def getkwlist(self):
    kwlist = []
    for k in self.kwlist:
      kwlist.append(self.__dict__[k])
    kwlist.append(self.__dict__['xcent'])
    kwlist.append(self.__dict__['ycent'])
    kwlist.append(self.__dict__['zcent'])
    return kwlist

  def griddistance(self,ix,iy,iz,xx,yy,zz):
    result = Delta(ix,iy,iz,xx,yy,zz,voltage=self.voltage,condid=self.condid,
                   generator=self.generatorf,kwlist=self.getkwlist())
    return result

  def distance(self,xx,yy,zz):
    result = Distance(xx,yy,zz,generator=self.generatord,
                      kwlist=self.getkwlist())
    return result

  def isinside(self,xx,yy,zz):
    result = IsInside(xx,yy,zz,generator=self.generatord,
                      condid=self.condid,kwlist=self.getkwlist())
    return result

  def intercept(self,xx,yy,zz,vx,vy,vz):
    result = Intercept(xx,yy,zz,vx,vy,vz,generator=self.generatori,
                       condid=self.condid,kwlist=self.getkwlist())
    return result

  def visualize(self,xmin=None,xmax=None,ymin=None,ymax=None):
    pass  # virtual function

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
    Assembly.__init__(self,0.,l.xcent,l.ycent,l.zcent,l.condid)
    self.left = l
  def griddistance(self,ix,iy,iz,xx,yy,zz):
    return (-(self.left.griddistance(ix,iy,iz,xx,yy,zz)))
  def distance(self,xx,yy,zz):
    return (-(self.left.distance(xx,yy,zz)))
  def isinside(self,xx,yy,zz):
    return (-(self.left.isinside(xx,yy,zz)))
  def visualize(self,xmin=None,xmax=None,ymin=None,ymax=None):
    self.left.visualize(xmin,xmax,ymin,ymax)


class AssemblyAnd(Assembly):
  """
AssemblyAnd class.  Represents 'and' of assemblies.
  """
  def __init__(self,l,r):
    Assembly.__init__(self,0.,l.xcent,l.ycent,l.zcent,l.condid)
    self.left = l
    self.right = r
  def griddistance(self,ix,iy,iz,xx,yy,zz):
    return (self.left.griddistance(ix,iy,iz,xx,yy,zz) *
            self.right.griddistance(ix,iy,iz,xx,yy,zz))
  def distance(self,xx,yy,zz):
    return (self.left.distance(xx,yy,zz) *
            self.right.distance(xx,yy,zz))
  def isinside(self,xx,yy,zz):
    return (self.left.isinside(xx,yy,zz) *
            self.right.isinside(xx,yy,zz))
  def visualize(self,xmin=None,xmax=None,ymin=None,ymax=None):
    self.left.visualize(xmin,xmax,ymin,ymax)
    self.right.visualize(xmin,xmax,ymin,ymax)


class AssemblyPlus(Assembly):
  """
AssemblyPlus class.  Represents 'or' of assemblies.
  """
  def __init__(self,l,r):
    Assembly.__init__(self,0.,l.xcent,l.ycent,l.zcent,l.condid)
    self.left = l
    self.right = r
  def griddistance(self,ix,iy,iz,xx,yy,zz):
    return (self.left.griddistance(ix,iy,iz,xx,yy,zz) +
            self.right.griddistance(ix,iy,iz,xx,yy,zz))
  def distance(self,xx,yy,zz):
    return (self.left.distance(xx,yy,zz) +
            self.right.distance(xx,yy,zz))
  def isinside(self,xx,yy,zz):
    return (self.left.isinside(xx,yy,zz) +
            self.right.isinside(xx,yy,zz))
  def visualize(self,xmin=None,xmax=None,ymin=None,ymax=None):
    self.left.visualize(xmin,xmax,ymin,ymax)
    self.right.visualize(xmin,xmax,ymin,ymax)


class AssemblyMinus(Assembly):
  """
AssemblyMinus class.
  """
  def __init__(self,l,r):
    Assembly.__init__(self,0.,l.xcent,l.ycent,l.zcent,l.condid)
    self.left = l
    self.right = r
  def griddistance(self,ix,iy,iz,xx,yy,zz):
    return (self.left.griddistance(ix,iy,iz,xx,yy,zz) -
            self.right.griddistance(ix,iy,iz,xx,yy,zz))
  def distance(self,xx,yy,zz):
    return (self.left.distance(xx,yy,zz) -
            self.right.distance(xx,yy,zz))
  def isinside(self,xx,yy,zz):
    return (self.left.isinside(xx,yy,zz) -
            self.right.isinside(xx,yy,zz))
  def visualize(self,xmin=None,xmax=None,ymin=None,ymax=None):
    self.left.visualize(xmin,xmax,ymin,ymax)
    self.right.visualize(xmin,xmax,ymin,ymax)

##############################################################################

class Delta:
  """
Class to hold the set of distances in each of the six directions.
Distances have the sign of the outward normal surface vector, i.e.
distances to outside the surface are positive, inside negative.
  """

  def __init__(self,ix=None,iy=None,iz=None,xx=None,yy=None,zz=None,
                    dels=None,vs=None,ns=None,
                    parity=None,voltage=0.,condid=1,generator=None,kwlist=[]):
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
      self.mglevel = zeros(nn)
    elif generator is not None:
      self.ndata = len(ix)
      self.ix = ix
      self.iy = iy
      self.iz = iz
      self.xx = xx
      self.yy = yy
      self.zz = zz
      self.dels = zeros((6,self.ndata),'d')
      fuzz = 1.e-13
      apply(generator,kwlist + [self.ndata,self.xx,self.yy,self.zz,
                                self.dels[0,:],self.dels[1,:],
                                self.dels[2,:],self.dels[3,:],
                                self.dels[4,:],self.dels[5,:]] + [fuzz])
      self.setvoltages(voltage)
      self.setcondids(condid)
      self.setlevels(0)
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
      self.ns = int(ns)
      self.parity = parity
      self.setlevels(0)
   
  def setvoltages(self,voltage):
    "Routine to set appropriate voltages."
    self.vs = voltage + zeros((6,self.ndata),'d')
   
  def setcondids(self,condid):
    "Routine to setcondid condids."
    self.ns = int(condid) + zeros((6,self.ndata))
   
  def setlevels(self,level):
    self.mglevel = level + zeros(self.ndata)

  def normalize(self,dx,dy,dz):
    """
Normalizes the data with respect to the grid cell sizes.
dx,dy,dz: the grid cell sizes
    """
    self.dels[:,:] = self.dels/array([dx,dx,dy,dy,dz,dz])[:,NewAxis]

  def setparity(self,dfill):
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
                       self.dels,self.parity,fuzz,f3d.fuzzsign,dfill)

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
    self.vs    = take(self.vs,ii,1)
    self.ns    = take(self.ns,ii,1)
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
      mglevel = self.mglevel[:n1]

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
      self.mglevel = zeros(newn)

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
      self.mglevel[:n1] = mglevel

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
    self.mglevel[n1:n1+n2] = d.mglevel[:n2]
    self.ndata = n1 + n2

  def install(self,installrz=1):
    """
Installs the data into the WARP database
    """
    ntot = 0
    nc = f3d.ncond
    nn = sum(where(self.parity[:self.ndata] == -1,1,0))
    ntot = ntot + nn
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
      f3d.icondlevel[nc:nc+nn] = take(self.mglevel,ii)

    ne = f3d.necndbdy
    nn = sum(where(self.parity[:self.ndata] == 0,1,0))
    ntot = ntot + nn
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
      f3d.iecndlevel[ne:ne+nn] = take(self.mglevel,ii)

    no = f3d.nocndbdy
    nn = sum(where(self.parity[:self.ndata] == 1,1,0))
    ntot = ntot + nn
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
      f3d.iocndlevel[no:no+nn] = take(self.mglevel,ii)
    if ntot > 0 and installrz:
      if(w3d.solvergeom == w3d.RZgeom or w3d.solvergeom == w3d.XZgeom):
        frz.install_conductors_rz()

  def __neg__(self):
    "Delta not operator."
    return Delta(self.ix,self.iy,self.iz,self.xx,self.yy,self.zz,
                 -self.dels,self.vs,self.ns)

  def __mul__(self,right):
    "'and' operator, returns maximum of distances to surfaces."
    c = less(self.dels,right.dels)
    return Delta(self.ix,self.iy,self.iz,self.xx,self.yy,self.zz,
                 choose(c,(self.dels,right.dels)),
                 choose(c,(self.vs  ,right.vs)),
                 choose(c,(self.ns  ,right.ns)))

  def __add__(self,right):
    "'or' operator, returns minimum of distances to surfaces."
    c = greater(self.dels,right.dels)
    return Delta(self.ix,self.iy,self.iz,self.xx,self.yy,self.zz,
                 choose(c,(self.dels,right.dels)),
                 choose(c,(self.vs  ,right.vs)),
                 choose(c,(self.ns  ,right.ns)))

  def __sub__(self,right):
    "'or' operator, returns minimum of distances to surfaces."
    rdels = -right.dels
    c = less(self.dels,rdels)
    return Delta(self.ix,self.iy,self.iz,self.xx,self.yy,self.zz,
                 choose(c,(self.dels,rdels)),
                 choose(c,(self.vs  ,right.vs)),
                 choose(c,(self.ns  ,right.ns)))

  def __str__(self):
    "Prints out delta"
    return repr(self.dels)+" "+repr(self.vs)+" "+repr(self.ns)

##############################################################################

class Distance:
  """
Class to hold the distance between points and a conductor
Distances have the sign of the outward normal surface vector, i.e.
distances outside the surface are positive, inside negative.
  """

  def __init__(self,xx=None,yy=None,zz=None,
                    distance=None,generator=None,kwlist=[]):
    if generator is not None:
      self.ndata = len(xx)
      self.xx = xx
      self.yy = yy
      self.zz = zz
      self.distance = zeros(self.ndata,'d')
      fuzz = 1.e-13
      apply(generator,kwlist + [self.ndata,self.xx,self.yy,self.zz,
                                self.distance[:]])
    else:
      self.ndata = len(xx)
      self.xx = xx
      self.yy = yy
      self.zz = zz
      self.distance = distance
   
  def __neg__(self):
    "Delta not operator."
    return Distance(self.xx,self.yy,self.zz, -self.distance)

  def __mul__(self,right):
    "'and' operator, returns maximum of distances to surfaces."
    c = less(abs(self.distance),abs(right.distance))
    return Distance(self.xx,self.yy,self.zz,
                    choose(c,(self.distance,right.distance)))

  def __add__(self,right):
    "'or' operator, returns minimum of distances to surfaces."
    c = greater(abs(self.distance),abs(right.distance))
    dd = Distance(self.xx,self.yy,self.zz,
                  choose(c,(self.distance,right.distance)))
    dd.distance = where((self.distance < 0.) & (right.distance > 0.),
                        self.distance,dd.distance)
    dd.distance = where((self.distance > 0.) & (right.distance < 0.),
                        right.distance,dd.distance)
    return dd

  def __sub__(self,right):
    "'or' operator, returns minimum of distances to surfaces."
    # --- Warning - while this work for many cases, there is no
    # --- gaurantee of robustness! It should only work when the right
    # --- hand side is a cylinder.
    rdistance = -right.distance
    dd = Distance(self.xx,self.yy,self.zz,self.distance)
    dd.distance = where((rdistance >= 0.) & (self.distance >= 0.),
                        sqrt(rdistance**2+self.distance**2),
                        dd.distance)
    dd.distance = where((rdistance >= 0.) & (self.distance <= 0.),
                        rdistance,dd.distance)
    dd.distance = where((rdistance < 0.) & (self.distance <= 0.),
                        maximum(rdistance,self.distance),dd.distance)
    return dd


  def __str__(self):
    "Prints out delta"
    return repr(self.distance)

##############################################################################

class IsInside:
  """
Class to hold flag whether or not a point is inside a conductor.
  """

  def __init__(self,xx=None,yy=None,zz=None,
                    isinside=None,generator=None,condid=1,kwlist=[]):
    self.condid = condid
    if generator is not None:
      self.ndata = len(xx)
      self.xx = xx
      self.yy = yy
      self.zz = zz
      distance = zeros(self.ndata,'d')
      fuzz = 1.e-13
      apply(generator,kwlist + [self.ndata,self.xx,self.yy,self.zz,
                                distance[:]])
      self.isinside = where(distance <= 0.,condid,0.)
    else:
      self.ndata = len(xx)
      self.xx = xx
      self.yy = yy
      self.zz = zz
      self.isinside = isinside*self.condid
   
  def __neg__(self):
    "Delta not operator."
    return IsInside(self.xx,self.yy,self.zz,
                    logical_not(self.isinside),condid=self.condid)

  def __mul__(self,right):
    "'and' operator, returns logical and of isinsides."
    return IsInside(self.xx,self.yy,self.zz,
                    logical_and(self.isinside,right.isinside),
                    condid=self.condid)

  def __add__(self,right):
    "'or' operator, returns logical or of isinsides."
    return IsInside(self.xx,self.yy,self.zz,
                    logical_or(self.isinside,right.isinside),
                    condid=self.condid)

  def __sub__(self,right):
    "'or' operator, returns logical or of isinsides."
    return IsInside(self.xx,self.yy,self.zz,
                    logical_and(self.isinside,logical_not(right.isinside)),
                    condid=self.condid)

  def __str__(self):
    "Prints out delta"
    return repr(self.isinside)

##############################################################################

class Intercept:
  """
Class to hold information about where a trajectory intercepted a conductor.
  """

  def __init__(self,xx=None,yy=None,zz=None,vx=None,vy=None,vz=None,
                    xi=None,yi=None,zi=None,angle=None,distance=None,
                    generator=None,condid=1,kwlist=[]):
    self.condid = condid
    self.ndata = len(xx)
    self.xx = xx
    self.yy = yy
    self.zz = zz
    self.vx = vx
    self.vy = vy
    self.vz = vz
    if generator is not None:
      self.xi = zeros(self.ndata,'d')
      self.yi = zeros(self.ndata,'d')
      self.zi = zeros(self.ndata,'d')
      self.angle = zeros(self.ndata,'d')
      self.distance = zeros(self.ndata,'d')
      apply(generator,kwlist + [self.ndata,self.xx,self.yy,self.zz,
                                self.vx,self.vy,self.vz,
                                self.xi,self.yi,self.zi,self.angle,
                                self.distance])
    else:
      self.xi = xi
      self.yi = yi
      self.zi = zi
      self.angle = angle
      self.distance = distance
   
  def __neg__(self):
    "Delta not operator."
    return Intercept(self.xx,self.yy,self.zz,self.vx,self.vy,self.vz,
                     self.xi,self.yi,self.zi,self.angle,self.distance,
                     condid=self.condid)

  def __mul__(self,right):
    "'and' operator, returns logical and."
    c = less(self.distance,right.distance)
    return Intercept(self.xx,self.yy,self.zz,self.vx,self.vy,self.vz,
                     choose(c,(self.xi,right.xi)),
                     choose(c,(self.yi,right.yi)),
                     choose(c,(self.zi,right.zi)),
                     choose(c,(self.angle,right.angle)),
                     choose(c,(self.distance,right.distance)),
                     condid=self.condid)

  def __add__(self,right):
    "'or' operator, returns logical or."
    c = greater(self.distance,right.distance)
    return Intercept(self.xx,self.yy,self.zz,self.vx,self.vy,self.vz,
                     choose(c,(self.xi,right.xi)),
                     choose(c,(self.yi,right.yi)),
                     choose(c,(self.zi,right.zi)),
                     choose(c,(self.angle,right.angle)),
                     choose(c,(self.distance,right.distance)),
                     condid=self.condid)

  def __sub__(self,right):
    "'or' operator, returns logical or."
    c = greater(self.distance,right.distance)
    return Intercept(self.xx,self.yy,self.zz,self.vx,self.vy,self.vz,
                     choose(c,(self.xi,right.xi)),
                     choose(c,(self.yi,right.yi)),
                     choose(c,(self.zi,right.zi)),
                     choose(c,(self.angle,right.angle)),
                     choose(c,(self.distance,right.distance)),
                     condid=self.condid)

  def __str__(self):
    "Prints out delta"
    return (repr(self.xi)+' '+repr(self.yi)+' '+repr(self.zi)+' '+
            repr(self.angle)+' '+repr(self.distance))

##############################################################################
##############################################################################
##############################################################################

# This is extremely kludgey. The z grid cell size is saved in this variable
# by the getdata function in Grid and then picked up in the Srfrv routines
# to pass into the fortran. There is not simple way of passing this
# information into the conductor object.
_griddzkludge = [0.]

##############################################################################
class Grid:
  """
Class holding the grid info.
Constructor arguments:
  - xmin,xmax,ymin,ymax,zmin,zmax: extent of conductors. Defaults to the
    mesh size. These only need to be set for optimization, to avoid looking
    for conductors where there are none. They can also be used to crop a
    conductor.
  - zbeam=top.zbeam: location of grid frame relative to lab frame
Call getdata(a,dfill) to generate the conductor data. 'a' is a geometry object.
Call installdata(installrz,gridmode) to install the data into the WARP database.
  """

  def __init__(self,xmin=None,xmax=None,ymin=None,ymax=None,
                    zmin=None,zmax=None,zbeam=None):
    """
Creates a grid object which can generate conductor data.
    """
    self.zbeam = where(zbeam==None,top.zbeam,zbeam)[0]
    self.xmin = where(xmin==None,w3d.xmmin,xmin)[0]
    self.xmax = where(xmax==None,w3d.xmmax,xmax)[0]
    self.ymin = where(ymin==None,w3d.ymmin,ymin)[0]
    self.ymax = where(ymax==None,w3d.ymmax,ymax)[0]
    self.zmin = where(zmin==None,w3d.zmmin+self.zbeam,zmin)[0]
    self.zmax = where(zmax==None,w3d.zmmax+self.zbeam,zmax)[0]

    self.nx = w3d.nx
    self.ny = w3d.ny
    self.nz = w3d.nz
    self.nzfull = w3d.nzfull
    self.xmmin = w3d.xmmin
    self.ymmin = w3d.ymmin
    self.zmmin = w3d.zmminglobal

    # --- Calculate dx, dy, and dz in case this is called before
    # --- the generate.
    self.dx = (w3d.xmmax - w3d.xmmin)/w3d.nx
    if(w3d.solvergeom==w3d.RZgeom or w3d.solvergeom==w3d.XZgeom):
      self.dy = self.dx
    else:
      self.dy = (w3d.ymmax - w3d.ymmin)/w3d.ny
    self.dz = (w3d.zmmax - w3d.zmmin)/w3d.nz

    if top.fstype in [7,11,10]:
      if top.fstype in [7,11]:
        setmglevels(self.nx,self.ny,self.nz,self.nzfull,self.dx,self.dy,self.dz)
      if top.fstype == 10:
        setmglevels_rz()
      self.mglevels = f3d.mglevels
      self.mglevelsnx = f3d.mglevelsnx[:f3d.mglevels]
      self.mglevelsny = f3d.mglevelsny[:f3d.mglevels]
      self.mglevelsiz = f3d.mglevelsiz[:f3d.mglevels]
      self.mglevelsnz = f3d.mglevelsnz[:f3d.mglevels]
      self.mglevelslx = f3d.mglevelslx[:f3d.mglevels]
      self.mglevelsly = f3d.mglevelsly[:f3d.mglevels]
      self.mglevelslz = f3d.mglevelslz[:f3d.mglevels]
    else:
      self.mglevels = 1
      self.mglevelsnx = [self.nx]
      self.mglevelsny = [self.ny]
      self.mglevelsiz = [top.izfsslave[me]]
      self.mglevelsnz = [self.nz]
      self.mglevelslx = [1]
      self.mglevelsly = [1]
      self.mglevelslz = [1]

  def getmesh(self,mglevel=0):
    i = mglevel
    dx = self.dx*self.mglevelslx[i]
    dy = self.dy*self.mglevelsly[i]
    dz = self.dz*self.mglevelslz[i]
    _griddzkludge[0] = dz
    nx = self.mglevelsnx[i]
    ny = self.mglevelsny[i]
    iz = self.mglevelsiz[i]
    nz = self.mglevelsnz[i]

    zmmin = self.zmmin + iz*dz

    xmesh = self.xmmin + dx*arange(nx+1)
    ymesh = self.ymmin + dy*arange(ny+1)
    zmesh =      zmmin + dz*arange(nz+1) + self.zbeam
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
    return ix,iy,iz,x,y,z,zmmin,dx,dy,dz,nx,ny,nz,zmesh

  def getdata(self,a,dfill=top.largepos):
    """
Given an Assembly, accumulate the appropriate data to represent that
Assembly on this grid.
 - a: the assembly
 - dfill=top.largepos: points at a depth in the conductor greater than dfill
                       are skipped.
    """
    starttime = wtime()
    tt2 = zeros(8,'d')
    self.dall = Delta()
    for i in range(self.mglevels):
      tt1 = wtime()
      ix,iy,iz,x,y,z,zmmin,dx,dy,dz,nx,ny,nz,zmesh = self.getmesh(i)

      tt2[0] = tt2[0] + wtime() - tt1
      if len(x) == 0: continue
      for zz in zmesh:
        tt1 = wtime()
        z[:] = zz
        iz[:] = nint((zz - zmmin - self.zbeam)/dz)
        tt2[1] = tt2[1] + wtime() - tt1
        tt1 = wtime()
        d = a.griddistance(ix,iy,iz,x,y,z)
        tt2[2] = tt2[2] + wtime() - tt1
        tt1 = wtime()
        d.normalize(dx,dy,dz)
        tt2[3] = tt2[3] + wtime() - tt1
        tt1 = wtime()
        d.setparity(dfill)
        tt2[4] = tt2[4] + wtime() - tt1
        tt1 = wtime()
        d.clean()
        tt2[5] = tt2[5] + wtime() - tt1
        tt1 = wtime()
        d.setlevels(i)
        tt2[6] = tt2[6] + wtime() - tt1
        tt1 = wtime()
        self.dall.append(d)
        tt2[7] = tt2[7] + wtime() - tt1
    endtime = wtime()
    self.generatetime = endtime - starttime
    #print tt2

  def installdata(self,installrz=1,gridmode=1):
    """
Installs the conductor data into the fortran database
    """
    self.dall.install(installrz)
    if gridmode is not None:
      f3d.gridmode = gridmode

  def getdistances(self,a,mglevel=0):
    """
Given an Assembly, accumulate the distances between the assembly and the
grid points.
 - a: the assembly
 - mglevel=0: coarsening level to use
    """
    starttime = wtime()
    tt2 = zeros(4,'d')
    tt1 = wtime()
    ix,iy,iz,x,y,z,zmmin,dx,dy,dz,nx,ny,nz,zmesh = self.getmesh(mglevel)
    try:
      self.distances[0,0,0]
    except:
      self.distances = fzeros((1+nx,1+ny,1+nz),'d')
    ix1 = min(ix)
    ix2 = max(ix)
    iy1 = min(iy)
    iy2 = max(iy)
    tt2[0] = tt2[0] + wtime() - tt1
    if len(x) == 0: return
    for zz in zmesh:
      tt1 = wtime()
      z[:] = zz
      iz[:] = nint((zz - zmmin - self.zbeam)/dz)
      tt2[1] = tt2[1] + wtime() - tt1
      tt1 = wtime()
      d = a.distance(x,y,z)
      tt2[2] = tt2[2] + wtime() - tt1
      tt1 = wtime()
      dd = d.distance
      dd.shape = (ix2-ix1+1,iy2-iy1+1)
      self.distances[ix1:ix2+1,iy1:iy2+1,iz[0]] = dd
      tt2[3] = tt2[3] + wtime() - tt1
    endtime = wtime()
    self.generatetime = endtime - starttime
    #print tt2

  def getisinside(self,a,mglevel=0):
    """
Given an Assembly, set flag for each grid point whether it is inside the
assembly.
 - a: the assembly
 - mglevel=0: coarsening level to use
    """
    starttime = wtime()
    tt2 = zeros(4,'d')
    tt1 = wtime()
    ix,iy,iz,x,y,z,zmmin,dx,dy,dz,nx,ny,nz,zmesh = self.getmesh(mglevel)
    try:
      self.isinside[0,0,0]
    except:
      self.isinside = fzeros((1+nx,1+ny,1+nz),'d')
    ix1 = min(ix)
    ix2 = max(ix)
    iy1 = min(iy)
    iy2 = max(iy)
    tt2[0] = tt2[0] + wtime() - tt1
    if len(x) == 0: return
    for zz in zmesh:
      tt1 = wtime()
      z[:] = zz #####
      iz[:] = nint((zz - zmmin - self.zbeam)/dz)  #####
      tt2[1] = tt2[1] + wtime() - tt1
      tt1 = wtime()
      d = a.isinside(x,y,z)  #####
      tt2[2] = tt2[2] + wtime() - tt1
      tt1 = wtime()
      dd = d.isinside
      dd.shape = (ix2-ix1+1,iy2-iy1+1)
      self.isinside[ix1:ix2+1,iy1:iy2+1,iz[0]] = where(dd>0,dd,
                                   self.isinside[ix1:ix2+1,iy1:iy2+1,iz[0]])
      tt2[3] = tt2[3] + wtime() - tt1
    endtime = wtime()
    self.generatetime = endtime - starttime
    #print tt2

##############################################################################
##############################################################################
##############################################################################

#============================================================================
class Plane(Assembly):
  """
Plane class
  - z0=0: locate of plane relative to zcent
  - zsign=1: when positive, conductor is in the z>0 side
  - theta=0,phi=0: normal of surface defining plane relative to z-axis
    theta is angle in z-x plane
    phi is angle in z-y plane
  - voltage=0: box voltage
  - xcent=0.,ycent=0.,zcent=0.: center of box
  - condid=1: conductor id of box, must be integer
  """
  def __init__(self,z0=0.,zsign=1.,theta=0.,phi=0.,
                    voltage=0.,xcent=0.,ycent=0.,zcent=0.,
                    condid=1):
    kwlist=['z0','zsign','theta','phi']
    Assembly.__init__(self,voltage,xcent,ycent,zcent,condid,kwlist,
                           planeconductorf,planeconductord,planeintercept)
    self.z0 = z0
    self.zsign = zsign
    self.theta = theta
    self.phi = phi

#============================================================================
class Box(Assembly):
  """
Box class
  - xsize,ysize,zsize: box size
  - voltage=0: box voltage
  - xcent=0.,ycent=0.,zcent=0.: center of box
  - condid=1: conductor id of box, must be integer
  """
  def __init__(self,xsize,ysize,zsize,voltage=0.,xcent=0.,ycent=0.,zcent=0.,
                    condid=1):
    kwlist=['xsize','ysize','zsize']
    Assembly.__init__(self,voltage,xcent,ycent,zcent,condid,kwlist,
                           boxconductorf,boxconductord,boxintercept)
    self.xsize = xsize
    self.ysize = ysize
    self.zsize = zsize

#============================================================================
class Cylinder(Assembly):
  """
Cylinder class
  - radius,length: cylinder size
  - theta=0,phi=0: angle of cylinder relative to z-axis
    theta is angle in z-x plane
    phi is angle in z-y plane
  - voltage=0: cylinder voltage
  - xcent=0.,ycent=0.,zcent=0.: center of cylinder
  - condid=1: conductor id of cylinder, must be integer
  """
  def __init__(self,radius,length,theta=0.,phi=0.,
                    voltage=0.,xcent=0.,ycent=0.,zcent=0.,condid=1):
    kwlist = ['radius','length','theta','phi']
    Assembly.__init__(self,voltage,xcent,ycent,zcent,condid,kwlist,
                           cylinderconductorf,cylinderconductord,
                           cylinderintercept)
    self.radius = radius
    self.length = length
    self.theta  = theta
    self.phi    = phi

#============================================================================
class Cylinders(Assembly):
  """
Cylinders class for a list of cylinders
  - radius,length: cylinder size
  - theta=0,phi=0: angle of cylinder relative to z-axis
    theta is angle in z-x plane
    phi is angle in z-y plane
  - voltage=0: cylinder voltage
  - xcent=0.,ycent=0.,zcent=0.: center of cylinder
  - condid=1: conductor id of cylinder, must be integer
  """
  def __init__(self,radius,length,theta=0.,phi=0.,
                    voltage=0.,xcent=0.,ycent=0.,zcent=0.,condid=1):
    kwlist = ['ncylinders','radius','length','theta','phi']
    Assembly.__init__(self,voltage,xcent,ycent,zcent,condid,kwlist,
                           cylindersconductorf,cylindersconductord,
                           cylindersintercept)
    self.ncylinders = 0
    self.radius = radius
    self.length = length
    self.theta  = theta
    self.phi    = phi
    kwlist = self.getkwlist()
    for k in kwlist:
      try:
        self.ncylinders = len(k)
        break
      except:
        pass

    assert self.ncylinders > 0,"At least on of the input arguments must be a list!"
    self.radius = self.radius*ones(self.ncylinders)
    self.length = self.length*ones(self.ncylinders)
    self.theta  = self.theta*ones(self.ncylinders)
    self.phi    = self.phi*ones(self.ncylinders)
    self.xcent  = self.xcent*ones(self.ncylinders)
    self.ycent  = self.ycent*ones(self.ncylinders)
    self.zcent  = self.zcent*ones(self.ncylinders)

#============================================================================
class ZCylinder(Assembly):
  """
Cylinder aligned with z-axis
  - radius,length: cylinder size
  - voltage=0: cylinder voltage
  - xcent=0.,ycent=0.,zcent=0.: center of cylinder
  - condid=1: conductor id of cylinder, must be integer
  """
  def __init__(self,radius,length,voltage=0.,xcent=0.,ycent=0.,zcent=0.,
                    condid=1):
    kwlist = ['radius','length']
    Assembly.__init__(self,voltage,xcent,ycent,zcent,condid,kwlist,
                           zcylinderconductorf,zcylinderconductord,
                           zcylinderintercept)
    self.radius = radius
    self.length = length

#============================================================================
class ZRoundedCylinder(Assembly):
  """
Cylinder with rounded corners aligned with z-axis
  - radius,length: cylinder size
  - radius2: radius of rounded corners
  - voltage=0: cylinder voltage
  - xcent=0.,ycent=0.,zcent=0.: center of cylinder
  - condid=1: conductor id of cylinder, must be integer
  """
  def __init__(self,radius,length,radius2,voltage=0.,xcent=0.,ycent=0.,zcent=0.,
                    condid=1):
    kwlist = ['radius','length','radius2']
    Assembly.__init__(self,voltage,xcent,ycent,zcent,condid,kwlist,
                      zroundedcylinderconductorf,zroundedcylinderconductord,
                      zroundedcylinderintercept)
    self.radius = radius
    self.length = length
    self.radius2 = radius2


#============================================================================
class ZCylinderOut(Assembly):
  """
Outside of a cylinder aligned with z-axis
  - radius,length: cylinder size
  - voltage=0: cylinder voltage
  - xcent=0.,ycent=0.,zcent=0.: center of cylinder
  - condid=1: conductor id of cylinder, must be integer
  """
  def __init__(self,radius,length,voltage=0.,xcent=0.,ycent=0.,zcent=0.,
                    condid=1):
    kwlist = ['radius','length']
    Assembly.__init__(self,voltage,xcent,ycent,zcent,condid,kwlist,
                      zcylinderoutconductorf,zcylinderoutconductord,
                      zcylinderoutintercept)
    self.radius = radius
    self.length = length

#============================================================================
class ZRoundedCylinderOut(Assembly):
  """
Outside of a cylinder with rounded corners aligned with z-axis
  - radius,length: cylinder size
  - radius2: radius of rounded corners
  - voltage=0: cylinder voltage
  - xcent=0.,ycent=0.,zcent=0.: center of cylinder
  - condid=1: conductor id of cylinder, must be integer
  """
  def __init__(self,radius,length,radius2,voltage=0.,xcent=0.,ycent=0.,zcent=0.,
                    condid=1):
    kwlist = ['radius','length','radius2']
    Assembly.__init__(self,voltage,xcent,ycent,zcent,condid,kwlist,
                      zroundedcylinderoutconductorf,
                      zroundedcylinderoutconductord,
                      zroundedcylinderoutintercept)
    self.radius = radius
    self.length = length
    self.radius2 = radius2

#============================================================================
class YCylinder(Assembly):
  """
Cylinder aligned with y-axis
  - radius,length: cylinder size
  - voltage=0: cylinder voltage
  - xcent=0.,ycent=0.,zcent=0.: center of cylinder
  - condid=1: conductor id of cylinder, must be integer
  """
  def __init__(self,radius,length,voltage=0.,xcent=0.,ycent=0.,zcent=0.,
                    condid=1):
    kwlist = ['radius','length']
    Assembly.__init__(self,voltage,xcent,ycent,zcent,condid,kwlist,
                      ycylinderconductorf,ycylinderconductord,
                      ycylinderintercept)
    self.radius = radius
    self.length = length

#============================================================================
class XCylinder(Assembly):
  """
Cylinder aligned with x-axis
  - radius,length: cylinder size
  - voltage=0: cylinder voltage
  - xcent=0.,ycent=0.,zcent=0.: center of cylinder
  - condid=1: conductor id of cylinder, must be integer
  """
  def __init__(self,radius,length,voltage=0.,xcent=0.,ycent=0.,zcent=0.,
                    condid=1):
    kwlist = ['radius','length']
    Assembly.__init__(self,voltage,xcent,ycent,zcent,condid,kwlist,
                      xcylinderconductorf,xcylinderconductord,
                      xcylinderintercept)
    self.radius = radius
    self.length = length

#============================================================================
class Sphere(Assembly):
  """
Sphere
  - radius: radius
  - voltage=0: sphere voltage
  - xcent=0.,ycent=0.,zcent=0.: center of sphere
  - condid=1: conductor id of sphere, must be integer
  """
  def __init__(self,radius,voltage=0.,xcent=0.,ycent=0.,zcent=0.,
                    condid=1):
    kwlist = ['radius']
    Assembly.__init__(self,voltage,xcent,ycent,zcent,condid,kwlist,
                      sphereconductorf,sphereconductord,sphereintercept)
    self.radius = radius

#============================================================================
class Cone(Assembly):
  """
Cone
  - r_zmin: radius at z min
  - r_zmax: radius at z max
  - length: length
  - theta=0,phi=0: angle of cylinder relative to z-axis
    theta is angle in z-x plane
    phi is angle in z-y plane
  - voltage=0: cone voltage
  - xcent=0.,ycent=0.,zcent=0.: center of cone
  - condid=1: conductor id of cone, must be integer
  """
  def __init__(self,r_zmin,r_zmax,length,theta,phi,voltage=0.,
                    xcent=0.,ycent=0.,zcent=0.,condid=1):
    kwlist = ['r_zmin','r_zmax','length','theta','phi']
    Assembly.__init__(self,voltage,xcent,ycent,zcent,condid,kwlist,
                      coneconductorf,coneconductord,coneintercept)
    self.r_zmin = r_zmin
    self.r_zmax = r_zmax
    self.theta = theta
    self.phi = phi
    self.length = length

#============================================================================
class ConeSlope(Assembly):
  """
Cone
  - slope: ratio of radius at zmax minus radius at zmin over length
  - intercept: location where line defining cone crosses the axis, relative
               to zcent
  - r_zmax: radius at z max
  - length: length
  - theta=0,phi=0: angle of cylinder relative to z-axis
    theta is angle in z-x plane
    phi is angle in z-y plane
  - voltage=0: cone voltage
  - xcent=0.,ycent=0.,zcent=0.: center of cone
  - condid=1: conductor id of cone, must be integer
  """
  def __init__(self,slope,length,theta,phi,voltage=0.,
                    xcent=0.,ycent=0.,zcent=0.,condid=1):
    kwlist = ['r_zmin','r_zmax','length','theta','phi']
    Assembly.__init__(self,voltage,xcent,ycent,zcent,condid,kwlist,
                      coneconductorf,coneconductord,coneintercept)
    self.slope = slope
    self.intercept = intercept
    self.theta = theta
    self.phi = phi
    self.length = length
  def getkwlist(self):
    self.r_zmin = self.slope*(-self.length/2. - self.intercept)
    self.r_zmax = self.slope*(+self.length/2. - self.intercept)
    return Assembly.getkwlist(self)

#============================================================================
class Cones(Assembly):
  """
Cones
  - r_zmin: radius at z min
  - r_zmax: radius at z max
  - length: length
  - theta=0,phi=0: angle of cylinder relative to z-axis
    theta is angle in z-x plane
    phi is angle in z-y plane
  - voltage=0: cone voltage
  - xcent=0.,ycent=0.,zcent=0.: center of cone
  - condid=1: conductor id of cone, must be integer
  """
  def __init__(self,r_zmin,r_zmax,length,theta,phi,voltage=0.,
                    xcent=0.,ycent=0.,zcent=0.,condid=1):
    kwlist = ['ncones','r_zmin','r_zmax','length','theta','phi']
    Assembly.__init__(self,voltage,xcent,ycent,zcent,condid,kwlist,
                      conesconductorf,conesconductord,conesintercept)
    self.ncones = 0
    self.r_zmin = r_zmin
    self.r_zmax = r_zmax
    self.length = length
    self.theta = theta
    self.phi = phi
    kwlist = self.getkwlist()
    for k in kwlist:
      try:
        self.ncones = len(k)
        break
      except:
        pass

    assert self.ncones > 0,"At least on of the input arguments must be a list!"
    self.r_zmin = self.r_zmin*ones(self.ncones)
    self.r_zmax = self.r_zmax*ones(self.ncones)
    self.length = self.length*ones(self.ncones)
    self.theta  = self.theta*ones(self.ncones)
    self.phi    = self.phi*ones(self.ncones)
    self.xcent  = self.xcent*ones(self.ncones)
    self.ycent  = self.ycent*ones(self.ncones)
    self.zcent  = self.zcent*ones(self.ncones)

#============================================================================
class ZCone(Assembly):
  """
Cone
  - r_zmin: radius at z min
  - r_zmax: radius at z max
  - length: length
  - voltage=0: cone voltage
  - xcent=0.,ycent=0.,zcent=0.: center of cone
  - condid=1: conductor id of cone, must be integer
  """
  def __init__(self,r_zmin,r_zmax,length,voltage=0.,xcent=0.,ycent=0.,zcent=0.,
                    condid=1):
    kwlist = ['r_zmin','r_zmax','length']
    Assembly.__init__(self,voltage,xcent,ycent,zcent,condid,kwlist,
                      zconeconductorf,zconeconductord,zconeintercept)
    self.r_zmin = r_zmin
    self.r_zmax = r_zmax
    self.length = length

#============================================================================
class ZConeSlope(Assembly):
  """
Cone
  - slope: ratio of radius at zmax minus radius at zmin over length
  - intercept: location where line defining cone crosses the axis, relative
               to zcent
  - length: length
  - voltage=0: cone voltage
  - xcent=0.,ycent=0.,zcent=0.: center of cone
  - condid=1: conductor id of cone, must be integer
  """
  def __init__(self,slope,intercept,length,voltage=0.,
                    xcent=0.,ycent=0.,zcent=0.,condid=1):
    kwlist = ['r_zmin','r_zmax','length']
    Assembly.__init__(self,voltage,xcent,ycent,zcent,condid,kwlist,
                      zconeconductorf,zconeconductord,zconeintercept)
    self.slope = slope
    self.intercept = intercept
    self.length = length
  def getkwlist(self):
    self.r_zmin = self.slope*(-self.length/2. - self.intercept)
    self.r_zmax = self.slope*(+self.length/2. - self.intercept)
    return Assembly.getkwlist(self)

#============================================================================
class ZConeOut(Assembly):
  """
Cone outside
  - r_zmin: radius at z min
  - r_zmax: radius at z max
  - length: length
  - voltage=0: cone voltage
  - xcent=0.,ycent=0.,zcent=0.: center of cone
  - condid=1: conductor id of cone, must be integer
  """
  def __init__(self,r_zmin,r_zmax,length,voltage=0.,xcent=0.,ycent=0.,zcent=0.,
                    condid=1):
    kwlist = ['r_zmin','r_zmax','length']
    Assembly.__init__(self,voltage,xcent,ycent,zcent,condid,kwlist,
                      zconeoutconductorf,zconeoutconductord,zconeoutintercept)
    self.r_zmin = r_zmin
    self.r_zmax = r_zmax
    self.length = length

#============================================================================
class ZConeOutSlope(Assembly):
  """
Cone outside
  - slope: ratio of radius at zmax minus radius at zmin over length
  - intercept: location where line defining cone crosses the axis, relative
               to zcent
  - length: length
  - voltage=0: cone voltage
  - xcent=0.,ycent=0.,zcent=0.: center of cone
  - condid=1: conductor id of cone, must be integer
  """
  def __init__(self,slope,intercept,length,voltage=0.,
                    xcent=0.,ycent=0.,zcent=0.,condid=1):
    kwlist = ['r_zmin','r_zmax','length']
    Assembly.__init__(self,voltage,xcent,ycent,zcent,condid,kwlist,
                      zconeoutconductorf,zconeoutconductord,zconeoutintercept)
    self.slope = slope
    self.intercept = intercept
    self.length = length
  def getkwlist(self):
    self.r_zmin = self.slope*(-self.length/2. - self.intercept)
    self.r_zmax = self.slope*(+self.length/2. - self.intercept)
    return Assembly.getkwlist(self)

#============================================================================
class ZTorus(Assembly):
  """
Torus
  - r1: toroidal radius
  - r2: poloidal radius
  - voltage=0: cone voltage
  - xcent=0.,ycent=0.,zcent=0.: center of cone
  - condid=1: conductor id of cone, must be integer
  """
  def __init__(self,r1,r2,voltage=0.,xcent=0.,ycent=0.,zcent=0.,condid=1):
    kwlist = ['r1','r2']
    Assembly.__init__(self,voltage,xcent,ycent,zcent,condid,kwlist,
                      ztorusconductorf,ztorusconductord,ztorusintercept)
    self.r1 = r1
    self.r2 = r2

#============================================================================
class Beamletplate(Assembly):
  """
Plate from beamlet pre-accelerator
  - za: location of spherical center in the x-plane
  - zb: location of spherical center in the y-plane
  - z0: location of the center of the plate on the z-axis
  - thickness: thickness of the plate
  - voltage=0: cone voltage
  - xcent=0.,ycent=0.,zcent=0.: center of cone
  - condid=1: conductor id of cone, must be integer
  """
  def __init__(self,za,zb,z0,thickness,voltage=0.,
               xcent=0.,ycent=0.,zcent=0.,condid=1):
    kwlist = ['za','zb','z0','thickness']
    Assembly.__init__(self,voltage,xcent,ycent,zcent,condid,kwlist,
                      beamletplateconductorf,beamletplateconductord,
                      beamletplateintercept)
    self.za = za
    self.zb = zb
    self.z0 = z0
    self.thickness = thickness

  def visualize(self,xmin=None,xmax=None,ymin=None,ymax=None):
    xmin = where(xmin==None,w3d.xmmin,xmin)[0]
    xmax = where(xmax==None,w3d.xmmax,xmax)[0]
    ymin = where(ymin==None,w3d.ymmin,ymin)[0]
    ymax = where(ymax==None,w3d.ymmax,ymax)[0]
    nx = min(w3d.nx,20)
    ny = min(w3d.ny,20)
    nz = w3d.nz
    xmmin = w3d.xmmin
    ymmin = w3d.ymmin
    zmmin = w3d.zmmin
    xmmax = w3d.xmmax
    ymmax = w3d.ymmax
    zmmax = w3d.zmmax

    # --- Calculate dx, dy, and dz in case this is called before
    # --- the generate.
    dx = (xmmax - xmmin)/nx
    dy = (ymmax - ymmin)/ny
    dz = (zmmax - zmmin)/nz
    if(w3d.solvergeom==w3d.RZgeom or w3d.solvergeom==w3d.XZgeom):
      dy = dx

    xmesh = xmmin + dx*arange(nx+1)
    ymesh = ymmin + dy*arange(ny+1)
    xmesh = compress(logical_and(xmin-dx <= xmesh,xmesh <= xmax+dx),xmesh)
    ymesh = compress(logical_and(ymin-dy <= ymesh,ymesh <= ymax+dy),ymesh)
    x = ravel(xmesh[:,NewAxis]*ones(len(ymesh)))
    y = ravel(ymesh*ones(len(xmesh))[:,NewAxis])
    ix = nint((x - xmmin)/dx)
    iy = nint((y - ymmin)/dy)
    if len(x) == 0: return

    xx = x
    yy = y
    xx.shape = (len(xmesh),len(ymesh))
    yy.shape = (len(xmesh),len(ymesh))

    # --- Outer face
    z = self.z0*ones(len(xmesh)*len(ymesh),'d') - self.thickness
    iz = nint((z - zmmin)/dz)
    d = self.griddistance(ix,iy,iz,x,y,z)
    zl = z[0] + d.dels[5,:]
    zl.shape = (len(xmesh),len(ymesh))
    ml = VPythonobjects.VisualMesh(xx,yy,zl,twoSided=true)

    # --- Inner face
    z = self.z0*ones(len(xmesh)*len(ymesh),'d') + 0.5*self.za
    iz = nint((z - zmmin)/dz)
    d = self.griddistance(ix,iy,iz,x,y,z)
    zr = z[0] - d.dels[4,:]
    zr.shape = (len(xmesh),len(ymesh))
    mr = VPythonobjects.VisualMesh(xx,yy,zr,twoSided=true)

    # --- Four sides between faces
    xside = xx[:,0]*ones(2)[:,NewAxis]
    yside = yy[:,0]*ones(2)[:,NewAxis]
    zside = array([zl[:,0],zr[:,0]])
    ms1 = VPythonobjects.VisualMesh(xside,yside,zside,twoSided=true)

    xside = xx[:,-1]*ones(2)[:,NewAxis]
    yside = yy[:,-1]*ones(2)[:,NewAxis]
    zside = array([zl[:,-1],zr[:,-1]])
    ms1 = VPythonobjects.VisualMesh(xside,yside,zside,twoSided=true)

    xside = xx[0,:]*ones(2)[:,NewAxis]
    yside = yy[0,:]*ones(2)[:,NewAxis]
    zside = array([zl[0,:],zr[0,:]])
    ms1 = VPythonobjects.VisualMesh(xside,yside,zside,twoSided=true)

    xside = xx[-1,:]*ones(2)[:,NewAxis]
    yside = yy[-1,:]*ones(2)[:,NewAxis]
    zside = array([zl[-1,:],zr[-1,:]])
    ms1 = VPythonobjects.VisualMesh(xside,yside,zside,twoSided=true)

#============================================================================
#============================================================================
class Srfrv(Assembly):
  """
Generic surface of revolution. Contains utility routines to check the
data and make sure it is consistent.
  """
  def checkarcs(self,zz,rr,rad,zc,rc):
    """Utility routine for surface of revoluation routines.
  Checks consistency of input and calculates any parameters not given.
    """
    for i in range(len(zz)-1):
      if ((rad[i] == None or rad[i] == largepos) and
          (zc[i] == None or zc[i] == largepos) and
          (rc[i] == None or rc[i] == largepos)):
        # --- When there is a straight line, then set the radius to a large
        # --- number (used as a flag in the code).
        rad[i] = largepos
        zc[i] = largepos
        rc[i] = largepos
      elif ((zc[i] == None or zc[i] == largepos) or
            (rc[i] == None or rc[i] == largepos)):
        # --- Given a radius and the two endpoints, the center of the
        # --- circle can be found.
        assert 4*rad[i]**2 > ((zz[i] - zz[i+1])**2 + (rr[i] - rr[i+1])**2),\
             "Radius of circle must be larger than the distance between points"
        zm = 0.5*(zz[i] + zz[i+1])
        rm = 0.5*(rr[i] + rr[i+1])
        dbm = sqrt((zm - zz[i+1])**2 + (rm - rr[i+1])**2)
        dcm = sqrt(rad[i]**2 - dbm**2)
        angle1 = arcsin((rm - rr[i+1])/dbm)
        if rad[i] < 0:
          zc[i] = zm + dcm*sin(angle1)
          rc[i] = rm + dcm*cos(angle1)
        else:
          zc[i] = zm - dcm*sin(angle1)
          rc[i] = rm - dcm*cos(angle1)
      elif (rad[i] == None or rad[i] == largepos):
        # --- Given the center, the radius can be found. With the two end
        # --- points of the arc given, the data is redundant, so check
        # --- to be sure it is consistent.
        rad[i] = sqrt((zz[i] - zc[i])**2 + (rr[i] - rc[i])**2)
        rad2 = sqrt((zz[i+1] - zc[i])**2 + (rr[i+1] - rc[i])**2)
        assert (abs(rad[i] - rad2)/rad[i] < 1.e-2),\
           "Points %d and %d are not at the same radius relative to the arc center. The radii are %e and %e"%(i,i+1,rad[i],rad2)
        # --- Make sure the radius has the correct sign.
        if rc[i] > rr[i] or rc[i] > rr[i+1]: rad[i] = -rad[i]

  def setdatadefaults(self,data,ndata,default):
    if data is None:
      data = ndata*[default]
    else:
      assert len(data) == ndata,\
             "Some of the input surface data is not the correct length"
    return data

#============================================================================
class ZSrfrvOut(Srfrv):
  """
Outside of a surface of revolution
  - rofzfunc: name of python function describing surface
  - zmin,zmax: z-extent of the surface
  - rmax=largepos: max radius of the surface
  - voltage=0: cone voltage
  - xcent=0.,ycent=0.,zcent=0.: center of cone
  - condid=1: conductor id of cone, must be integer
  - rofzdata=None: optional tablized data of radius of surface
  - zdata=None: optional tablized data of z locations of rofzdata
      raddata[i] is radius for segment from zdata[i] to zdata[i+1]
  - zcdata=None: z center of circle or curved segment
  - rcdata=None: r center of circle or curved segment
  - raddata=None: optional radius of curvature of segments
      The centers of the circles will be calculated automatically if
      not supplied, or if the centers are supplied, the radii will be
      calculated automatically.
      The length of the radii and centers lists is one less than the length
      of the list of r and z data.
    Note that if tablized data is given, the first argument is ignored.
  """
  def __init__(self,rofzfunc,zmin,zmax,rmax=largepos,
                    voltage=0.,xcent=0.,ycent=0.,zcent=0.,condid=1,
                    rofzdata=None,zdata=None,raddata=None,
                    zcdata=None,rcdata=None):
    kwlist = ['rofzfunc','zmin','zmax','rmax','griddz']
    Assembly.__init__(self,voltage,xcent,ycent,zcent,condid,kwlist,
                      zsrfrvoutconductorf,zsrfrvoutconductord,
                      zsrfrvoutintercept)
    self.rofzfunc = rofzfunc
    self.zmin = zmin
    self.zmax = zmax
    self.rmax = rmax

    # --- Deal with tablized data.
    # --- Make sure the input is consistent
    if operator.isSequenceType(rofzdata):
      self.usedata = true
      self.zdata = zdata
      self.rofzdata = self.setdatadefaults(rofzdata,len(zdata),rmax)
      self.raddata = self.setdatadefaults(raddata,len(zdata)-1,None)
      self.zcdata = self.setdatadefaults(zcdata,len(zdata)-1,None)
      self.rcdata = self.setdatadefaults(rcdata,len(zdata)-1,None)
      self.checkarcs(self.zdata,self.rofzdata,self.raddata,
                     self.zcdata,self.rcdata)
      self.rofzfunc = ' '
    else:
      self.usedata = false

  def getkwlist(self):
    self.griddz = _griddzkludge[0]
    # --- Make sure the rofzfunc is in main.
    # --- Note that this can only really work if a reference to the function
    # --- is passed in (instead of the name).
    if type(self.rofzfunc) == FunctionType:
      import __main__
      __main__.__dict__[self.rofzfunc.__name__] = self.rofzfunc

    # --- Get the name of the input function if a reference to the function
    # --- was passed in.
    if type(self.rofzfunc) == FunctionType:
      self.rofzfunc = self.rofzfunc.__name__

    # --- If data arrays are specified, then put the data in the right place
    if self.usedata:
      f3d.lsrlinr = true
      f3d.npnts_sr = len(self.zdata)
      f3d.forceassign('z_sr',self.zdata)
      f3d.forceassign('r_sr',self.rofzdata)
      f3d.forceassign('rad_sr',self.raddata)
      f3d.forceassign('zc_sr',self.zcdata)
      f3d.forceassign('rc_sr',self.rcdata)
    else:
      f3d.lsrlinr = false

    return Assembly.getkwlist(self)

#============================================================================
class ZSrfrvIn(Srfrv):
  """
Inside of a surface of revolution
  - rofzfunc: name of python function describing surface
  - zmin,zmax: z-extent of the surface
  - rmin=0: min radius of the surface
  - voltage=0: cone voltage
  - xcent=0.,ycent=0.,zcent=0.: center of cone
  - condid=1: conductor id of cone, must be integer
  - rofzdata=None: optional tablized data of radius of surface
  - zdata=None: optional tablized data of z locations of rofzdata
  - raddata=None: optional radius of curvature of segments
      raddata[i] is radius for segment from zdata[i] to zdata[i+1]
  - zcdata=None: z center of circle or curved segment
  - rcdata=None: r center of circle or curved segment
      The centers of the circles will be calculated automatically if
      not supplied, or if the centers are supplied, the radii will be
      calculated automatically.
      The length of the radii and centers lists is one less than the length
      of the list of r and z data.
    Note that if tablized data is given, the first argument is ignored.
  """
  def __init__(self,rofzfunc,zmin,zmax,rmin=0,
                    voltage=0.,xcent=0.,ycent=0.,zcent=0.,condid=1,
                    rofzdata=None,zdata=None,raddata=None,
                    zcdata=None,rcdata=None):
    kwlist = ['rofzfunc','zmin','zmax','rmin','griddz']
    Assembly.__init__(self,voltage,xcent,ycent,zcent,condid,kwlist,
                      zsrfrvinconductorf,zsrfrvinconductord,
                      zsrfrvinintercept)
    self.rofzfunc = rofzfunc
    self.zmin = zmin
    self.zmax = zmax
    self.rmin = rmin

    # --- Deal with tablized data.
    # --- Make sure the input is consistent
    if operator.isSequenceType(rofzdata):
      self.usedata = true
      self.zdata = zdata
      self.rofzdata = self.setdatadefaults(rofzdata,len(zdata),rmin)
      self.raddata = self.setdatadefaults(raddata,len(zdata)-1,None)
      self.zcdata = self.setdatadefaults(zcdata,len(zdata)-1,None)
      self.rcdata = self.setdatadefaults(rcdata,len(zdata)-1,None)
      self.checkarcs(self.zdata,self.rofzdata,self.raddata,
                     self.zcdata,self.rcdata)
      self.rofzfunc = ' '
    else:
      self.usedata = false

  def getkwlist(self):
    self.griddz = _griddzkludge[0]
    # --- Make sure the rofzfunc is in main.
    # --- Note that this can only really work if a reference to the function
    # --- is passed in (instead of the name).
    if type(self.rofzfunc) == FunctionType:
      import __main__
      __main__.__dict__[self.rofzfunc.__name__] = self.rofzfunc

    # --- Get the name of the input function if a reference to the function
    # --- was passed in.
    if type(self.rofzfunc) == FunctionType:
      self.rofzfunc = self.rofzfunc.__name__

    # --- If data arrays are specified, then put the data in the right place
    if self.usedata:
      f3d.lsrlinr = true
      f3d.npnts_sr = len(self.zdata)
      f3d.forceassign('z_sr',self.zdata)
      f3d.forceassign('r_sr',self.rofzdata)
      f3d.forceassign('rad_sr',self.raddata)
      f3d.forceassign('zc_sr',self.zcdata)
      f3d.forceassign('rc_sr',self.rcdata)
    else:
      f3d.lsrlinr = false

    return Assembly.getkwlist(self)

#============================================================================
class ZSrfrvInOut(Srfrv):
  """
Betweem surfaces of revolution
  - rminofz,rmaxofz: names of python functions describing surfaces
  - zmin,zmax: z-extent of the surface
  - voltage=0: cone voltage
  - xcent=0.,ycent=0.,zcent=0.: center of cone
  - condid=1: conductor id of cone, must be integer
  - rminofzdata,rmaxofzdata=None: optional tablized data of radii of surface
  - zmindata,zmaxdata=None: optional tablized data of z locations of r data
  - radmindata,radmaxdata=None: optional radius of curvature of segments
      radmindata[i] is radius for segment from zmindata[i] to zmindata[i+1]
  - zcmindata,zcmaxdata=None: z center of circle or curved segment
  - rcmindata,rcmaxdata=None: r center of circle or curved segment
      The centers of the circles will be calculated automatically if
      not supplied, or if the centers are supplied, the radii will be
      calculated automatically.
      The length of the radii and centers lists is one less than the length
      of the list of r and z data.
    Note that if tablized data is given, the first two arguments are ignored.
  """
  def __init__(self,rminofz,rmaxofz,zmin,zmax,
                    voltage=0.,xcent=0.,ycent=0.,zcent=0.,condid=1,
                    rminofzdata=None,zmindata=None,radmindata=None,
                    rcmindata=None,zcmindata=None,
                    rmaxofzdata=None,zmaxdata=None,radmaxdata=None,
                    rcmaxdata=None,zcmaxdata=None):
    kwlist = ['rminofz','rmaxofz','zmin','zmax','griddz']
    Assembly.__init__(self,voltage,xcent,ycent,zcent,condid,kwlist,
                      zsrfrvinoutconductorf,zsrfrvinoutconductord,
                      zsrfrvinoutintercept)
    self.rminofz = rminofz
    self.rmaxofz = rmaxofz
    self.zmin = zmin
    self.zmax = zmax

    # --- Deal with tablized data.
    # --- Making sure the input is consistent
    if operator.isSequenceType(zmindata):
      self.usemindata = true
      self.zmindata = zmindata
      self.rminofzdata = self.setdatadefaults(rminofzdata,len(zmindata),0.)
      self.radmindata = self.setdatadefaults(radmindata,len(zmindata)-1,None)
      self.rcmindata = self.setdatadefaults(rcmindata,len(zmindata)-1,None)
      self.zcmindata = self.setdatadefaults(zcmindata,len(zmindata)-1,None)
      self.checkarcs(self.zmindata,self.rminofzdata,self.radmindata,
                     self.zcmindata,self.rcmindata)
      self.rminofz = ' '
      self.rmaxofz = ' '
    else:
      self.usemindata = false

    if operator.isSequenceType(zmaxdata):
      self.usemaxdata = true
      self.zmaxdata = zmaxdata
      self.rmaxofzdata = self.setdatadefaults(rmaxofzdata,len(zmaxdata),
                                              largepos)
      self.radmaxdata = self.setdatadefaults(radmaxdata,len(zmaxdata)-1,None)
      self.rcmaxdata = self.setdatadefaults(rcmaxdata,len(zmaxdata)-1,None)
      self.zcmaxdata = self.setdatadefaults(zcmaxdata,len(zmaxdata)-1,None)
      self.checkarcs(self.zmaxdata,self.rmaxofzdata,self.radmaxdata,
                     self.zcmaxdata,self.rcmaxdata)
      self.rminofz = ' '
      self.rmaxofz = ' '
    else:
      self.usemaxdata = false

  def getkwlist(self):
    self.griddz = _griddzkludge[0]
    # --- Make sure the rminofz and rmaxofz are in main.
    # --- Note that this can only really work if a reference to the function
    # --- is passed in (instead of the name).
    if type(self.rminofz) == FunctionType:
      import __main__
      __main__.__dict__[self.rminofz.__name__] = self.rminofz
    if type(self.rmaxofz) == FunctionType:
      import __main__
      __main__.__dict__[self.rmaxofz.__name__] = self.rmaxofz

    # --- Get the name of the input function if a reference to the function
    # --- was passed in.
    if type(self.rminofz) == FunctionType:
      self.rminofz = self.rminofz.__name__
    if type(self.rmaxofz) == FunctionType:
      self.rmaxofz = self.rmaxofz.__name__

    # --- If data arrays are specified, then put the data in the right place
    if self.usemindata:
      f3d.lsrminlinr = true
      f3d.npnts_srmin = len(self.zmindata)
      f3d.forceassign('z_srmin',self.zmindata)
      f3d.forceassign('r_srmin',self.rminofzdata)
      f3d.forceassign('rad_srmin',self.radmindata)
      f3d.forceassign('zc_srmin',self.zcmindata)
      f3d.forceassign('rc_srmin',self.rcmindata)
    else:
      f3d.lsrminlinr = false

    if self.usemaxdata:
      f3d.lsrmaxlinr = true
      f3d.npnts_srmax = len(self.zmaxdata)
      f3d.forceassign('z_srmax',self.zmaxdata)
      f3d.forceassign('r_srmax',self.rmaxofzdata)
      f3d.forceassign('rad_srmax',self.radmaxdata)
      f3d.forceassign('zc_srmax',self.zcmaxdata)
      f3d.forceassign('rc_srmax',self.rcmaxdata)
    else:
      f3d.lsrmaxlinr = false

    return Assembly.getkwlist(self)

#============================================================================
#============================================================================
#============================================================================
class SRFRVLApart:
  """
Class for creating a surface of revolution conductor part using lines and arcs as primitives.  
 - name:   name of conductor part
 - data:   list of lines and arcs (this must be a continuous curve)
            o the first element is the starting point
            o a line is described by the ending point
            o a circle is described by the ending point and the center
            o for both, the starting point is the previous ending point
            o each element is teminated by a letter 't' for top (or out)
              and 'b' for bottom (or in)
  """
  def __init__(self,name,data):
    self.name    = name
    self.data   = data

  def draw(self,ncirc=50,scx=1.,scy=1.,colort='blue',colorb='red',
                 color='none',signx=1.,width=1.):
    """
  Draw lines and arcs of SRFRVLA conductor part. Arcs are decomposed
  into ncirc segments.
    """
    if(color <> 'none'):
      colort=color
      colorb=color
      
    for j,d in enumerate(self.data):
      if(j==0):
        s=[d[0]*scx,d[1]*scy]
      else:
        if(d[-1]=='t'):
          colorn=colort
        else:
          colorn=colorb
        n=[d[0]*scx,d[1]*scy]
        if(len(d)==3):
          self.draw_line([s[0]],[signx*s[1]],[n[0]],[signx*n[1]],colorn,width)
        else:
          c=[d[2]*scx,d[3]*scy]
          self.draw_arc(s[0],signx*s[1],n[0],signx*n[1],c[0],signx*c[1],ncirc,colorn,width)
        s=n

  def shift(self,s,dir):
    """
  Shift SRFRVLA conductor part by s in the direction dir (0=Z, 1=R).
    """
    for d in self.data:
      d[dir]+=s
      if(len(d)>3):
        d[dir+2]+=s

  def draw_line(self,x0,y0,x1,y1,color='black',width=1.):
    """
  Draw a line.
    """
    pldj(x0,y0,x1,y1,color=color,width=width)

  def draw_arc(self,x0,y0,x1,y1,xc,yc,ncirc=50,color='black',width=1.):
    """
  Draw an arc as ncirc segments.
    """
    xy = self.get_xy_arc(x0,y0,x1,y1,xc,yc,ncirc)  
    pldj(xy[0][0:ncirc-1],xy[1][0:ncirc-1],
         xy[0][1:ncirc],  xy[1][1:ncirc],
         color=color,width=width)

  def get_xy_arc(self,x0,y0,x1,y1,xc,yc,ncirc):
    """
  Retuns list of end points of ncirc segments along an arc of
  end points (x0,y0) and (x1,y1), and center (xc,yc).
    """
    xi=x0-xc
    yi=y0-yc
    xf=x1-xc
    yf=y1-yc
    r=sqrt(xi**2+yi**2)
    thetai=arctan2(yi,xi)
    thetaf=arctan2(yf,xf)
    if(abs(thetaf-thetai)>pi):
      if(thetaf-thetai<0.):
        thetaf=thetaf+2.*pi
      else:
        thetaf=thetaf-2.*pi
    dtheta=(thetaf-thetai)/(ncirc-1)
    x=zeros(ncirc+1,'d')
    y=zeros(ncirc+1,'d')
    x[0]=x0
    y[0]=y0
    for i in range(1,ncirc):
      theta=thetai+i*dtheta
      x[i]=xc+r*cos(theta)
      y[i]=yc+r*sin(theta)
    return [x,y]

class SRFRVLAcond:
  """
Class for creating Surface of revolution conductor using lines and arcs as primitives.
A conductor is a structure containing a list of parts. A part is a structure
containing a list of primitives.
 - name: name of conductor
 - parts: list of conductor parts
 - voltage: conductor part voltage
 - condid:  conductor part ID
 - install=1 (optional): flag for installation of conductors immediately after creation.
  """
  def __init__(self,name,parts,voltage,condid,install=1):
    """
  Register SRFRVLA conductor parts using ZSrfrvIn, ZSrfrvOut or ZSrfrvInOut.
    """
    self.name = name
    self.parts = parts
    self.voltage = voltage
    self.condid = condid
    for i,part in enumerate(self.parts):
    # Loop over parts and install them using ZSrfrvIn, ZSrfrvOut or ZSrfrvInOut.
    # The installed parts are stored in the list 'toinstall' and the corresponding
    # IDs and names in the lists 'condids' and 'cnames'.
      data = part.data
      it=0
      ib=0
      nsrmax=0
      nsrmin=0
      l_t=0
      l_b=0
      # Counts number of 'top' and 'bottom' primitives and allocate temporary arrays.
      for d in data[1:]:
        if(d[-1]=='t'):
            nsrmax += 1
        elif(d[-1]=='b'):
            nsrmin += 1
      if(nsrmax<>0):
        l_t = 1
        r_srmax = zeros(nsrmax+1,Float)
        z_srmax = zeros(nsrmax+1,Float)
        rc_srmax = zeros(nsrmax,Float)
        zc_srmax = zeros(nsrmax,Float)
      if(nsrmin<>0):
        l_b = 1
        r_srmin = zeros(nsrmin+1,Float)
        z_srmin = zeros(nsrmin+1,Float)
        rc_srmin = zeros(nsrmin,Float)
        zc_srmin = zeros(nsrmin,Float)

      # Fill arrays with datas from parts.
      do = data[0] 
      for d in data[1:]:
        if(d[-1]=='t'):
          z_srmax[it]=do[0]
          r_srmax[it]=do[1]
          it=it+1
          z_srmax[it]=d[0]
          r_srmax[it]=d[1]
          if(len(d)>3):
            zc_srmax[it-1]=d[2]
            rc_srmax[it-1]=d[3]
          else:
            zc_srmax[it-1]=largepos
            rc_srmax[it-1]=largepos
        elif(d[-1]=='b'):
          z_srmin[ib]=do[0]
          r_srmin[ib]=do[1]
          ib=ib+1
          z_srmin[ib]=d[0]
          r_srmin[ib]=d[1]
          if(len(d)>3):
            zc_srmin[ib-1]=d[2]
            rc_srmin[ib-1]=d[3]
          else:
            zc_srmin[ib-1]=largepos
            rc_srmin[ib-1]=largepos
        do=d

      # Make sure arrays are in Z ascending order.
      if(z_srmin[0]>z_srmin[nsrmin]):
        args =argsort(z_srmin)
        argsc=argsort(z_srmin[1:])
        z_srmin=take(z_srmin,args)
        r_srmin=take(r_srmin,args)
        zc_srmin=take(zc_srmin,argsc)
        rc_srmin=take(rc_srmin,argsc)
      if(z_srmax[0]>z_srmax[nsrmax]):
        args =argsort(z_srmax)
        argsc=argsort(z_srmax[1:])
        z_srmax=take(z_srmax,args)
        r_srmax=take(r_srmax,args)
        zc_srmax=take(zc_srmax,argsc)
        rc_srmax=take(rc_srmax,argsc)

      # Register parts.
      if(l_t==1 and l_b==1):
        part.installed = ZSrfrvInOut('','',
                              min(z_srmax),
                              max(z_srmax),
                              voltage=voltage,
                              condid =condid,
                              rminofzdata=r_srmin,
                              rmaxofzdata=r_srmax,
                              zmindata=z_srmin,
                              zmaxdata=z_srmax,
                              rcmindata=rc_srmin,
                              rcmaxdata=rc_srmax,
                              zcmindata=zc_srmin,
                              zcmaxdata=zc_srmax,
                              )

      elif(l_b==1):
        part.installed = ZSrfrvOut('',
                            min(z_srmin),
                            max(z_srmin),
                            voltage=voltage,
                            condid =condid,
                            rofzdata=r_srmin,
                            zdata=z_srmin)

      elif(l_t==1):
        part.installed = ZSrfrvIn('',
                            min(z_srmax),
                            max(z_srmax),
                            voltage=voltage,
                            condid =condid,
                            rofzdata=r_srmax,
                            zdata=z_srmax)

      # store installed conductor in a list
      if i == 0:
        self.cond = self.parts[i].installed
      else:    
        self.cond += self.parts[i].installed
    if(install):self.install()
    
  def install(self):
    """
  Install SRFRVLA conductors.
    """
    print 'installing',self.name,'( ID=',self.condid,')...'
    for part in self.parts:
      installconductors(part.installed)
    
  def draw(self,ncirc=50,scx=1.,scy=1.,colort='blue',colorb='red',
                 color='none',signx=1.,width=1.):
    """
  Draws a list of conductor parts.
    """
    if(color <> 'none'):
      colort=color
      colorb=color
    for part in self.parts:
      part.draw(ncirc,scx,scy,colort,colorb,color,signx,width)
  
  def shift(self,s,dir):
    """
  Shifts conductor parts using a scalar shift of a list of shifts.
    - s: shift value (scalar or list same length as self.parts)
    - dir: direction of shift (0=Z, 1=R)
    """
    stype = 'list'
    try:
      l=len(s)
    except:
      stype = 'scalar'

    if(stype=='scalar'):
      for part in self.parts:
        part.shift(s,dir)
    else:
      assert len(s)==len(self.parts), 'Error in SRFRVLAcond.shift: shift must have same length as part or be a scalar.'
      for i,part in enumerate(self.parts):
        part.shift(s[i],dir)

class SRFRVLAsystem:
  """
Class for creating a SRFRVLAsystem of conductors (list of conductors).
A SRFRVLA contains two lists:
  - SRFRVLAconds which is the list of SRFRVLA conductors (SRFRVLAcond).
  - conds which is the list of conductors, as defined in generateconductors.py.
  """
  def __init__(self,SRFRVLAconds,install=1):
    self.conds = []
    self.SRFRVLAconds = SRFRVLAconds
    for SRFRVLAcond in SRFRVLAconds:
      self.conds += [SRFRVLAcond.cond]

  def save(self,filename):
    """
  Save conductor into external file:
    - filename: name of external file ('.wob' will be added')
    """
    f=open(filename+'.wob','w')
    f.write("Begin\n")

    for cond in self.SRFRVLAconds:
      f.write("  Conductor "+cond.name+"\n")
      for part in cond.parts:
        f.write("    Part "+part.name+"\n")
        d = part.data[0]
        f.write('     s   %13.6E %13.6E \n'%(d[0],d[1]))
        for d in part.data[1:]:
          if(len(d)==3): # line
            f.write('     l   %13.6E %13.6E                             %c \n'%(d[0],d[1],d[2]))
          else: # arc
            f.write('     a   %13.6E %13.6E %13.6E %13.6E %c \n'%(d[0],d[1],d[2],d[3],d[4]))
        f.write("    Endpart "+part.name+"\n")
      f.write("  End "+cond.name+"\n")

    f.write("End\n")
    f.close()

  def draw(self,ncirc=50,scx=1.,scy=1.,colort='blue',colorb='red',
                 color='none',signx=1.,width=1.):
    """
  Draws a list of conductors.
    """
    if(color <> 'none'):
      colort=color
      colorb=color
    for parts in self.SRFRVLAconds:
        parts.draw(ncirc,scx,scy,colort,colorb,color,signx,width)

class SRFRVLAfromfile(SRFRVLAsystem):
  """
Class for reading SRFRVLA data from file.
  - filename: name of input file
  - voltages: list of voltages of conductor parts 
  - condids: list of IDs of conductor parts 
  - zshifts=None (optional): list of shifts in Z to apply to conductor parts 
  - rshifts=None (optional): list of shifts in R to apply to conductor parts 
  - install=1 (optional): flag for installation of conductors immediately after reading 
  """
  def __init__(self,filename,voltages,condids,zshifts=None,rshifts=None,install=1):
    """
  Reads SRFRVLA conductors from external file. The series of conductors which is stored in a list.
    """
    f=open(filename+'.wob','r')
    cont=1
    ic = -1
    ip = -1
    # create list conductors
    conds=[]
    condnames = []
    while(cont==1):
      line=split(f.readline())
      if(line[0]=='Conductor'):
        ic += 1
        condnames += [line[1]]
        partnames = []
        parts = []
        # create list of data
        data = []
        ocont = 1
        while(ocont==1):
          line=split(f.readline())
          if(line[0]=='End'):
            ocont=0
          elif(line[0]=='Part'):
            ip += 1
            partnames += [line[1]]
          elif(line[0]=='Endpart'):
            parts += [SRFRVLApart(partnames[-1],data)]
            if zshifts is not None:parts[-1].shift(zshifts[ip],0)
            if rshifts is not None:parts[-1].shift(rshifts[ip],1)
            data = []
            partnames = []
          else:
            # append data to list
            if(line[0]=='S' or line[0]=='s'): #start
              data.append([float(line[1]),float(line[2])])
            if(line[0]=='L' or line[0]=='l'): #line
              data.append([float(line[1]),float(line[2]),line[3]])
            if(line[0]=='A' or line[0]=='a'): #arc
              data.append([float(line[1]),float(line[2]),\
                           float(line[3]),float(line[4]),line[5]])
        conds += [parts]
      elif(line[0]=='End'):
        cont=0
    f.close()
    # create conductors
    SRFRVLAconds = []
    for i,parts in enumerate(conds):
      SRFRVLAconds += [SRFRVLAcond(condnames[i],parts,voltages[i],condids[i],install)]
    SRFRVLAsystem.__init__(self,SRFRVLAconds)

class SRFRVLA_circle(SRFRVLApart):
  """
Creates a circle.
  - name: name of conductor part
  - c: center
  - r: radius
  """
  def __init__(self,name,c,r):
    assert r>0., 'In SRFRVLA_circle, r must be >0.' 
    x = c[0]
    y = c[1]
    data = []
    data.append([x-r,y])
    data.append([x+r,y,x,y,'t'])
    data.append([x,y-r,x,y,'b'])
    data.append([x-r,y,x,y,'b'])
    SRFRVLApart.__init__(self,name,data)
  
class SRFRVLA_rectangle(SRFRVLApart):
  """
Creates a rectangle.
  - name: name of conductor part
  - c: center
  - l: length
  - h: height
  """
  def __init__(self,name,c,l,h):
    assert l>0. and h>0., 'In SRFRVLA_rectangle, l, h must be >0.' 
    x1 = c[0]-0.5*l
    x2 = c[0]+0.5*l
    y1 = c[1]-0.5*h
    y2 = c[1]+0.5*h
    data = []
    data.append([x1,y1])
    data.append([x1,y2,'t'])
    data.append([x2,y2,'t'])
    data.append([x2,y1,'t'])
    data.append([x1,y1,'b'])
    SRFRVLApart.__init__(self,name,data)
  
class SRFRVLA_rnd_rectangle(SRFRVLApart):
  """
Creates a rectangle with rounded edge.
  - name: name of conductor part
  - c: center
  - l: length
  - h: height
  - r: radius of rounded edges
  """
  def __init__(self,name,c,l,h,r):
    assert r>=0. and l>=0. and h>=0., 'In SRFRVLA_rnd_rectangle, r, l, h must be >=0.' 
    l = max(l,2.*r)
    h = max(h,2.*r)
    x1 = c[0]-0.5*l
    x2 = c[0]+0.5*l
    y1 = c[1]-0.5*h
    y2 = c[1]+0.5*h
    data = []
    data.append([x1,y1+r])
    data.append([x1,y2-r,'t'])
    if(r>0.):data.append([x1+r,y2,x1+r,y2-r,'t'])
    data.append([x2-r,y2,'t'])
    if(r>0.):data.append([x2,y2-r,x2-r,y2-r,'t'])
    data.append([x2,y1+r,'t'])
    if(r>0.):data.append([x2-r,y1,x2-r,y1+r,'b'])
    data.append([x1+r,y1,'b'])
    if(r>0.):data.append([x1,y1+r,x1+r,y1+r,'b'])
    SRFRVLApart.__init__(self,name,data)
