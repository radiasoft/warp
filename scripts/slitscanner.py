"""
Creates a class which generates slit scanner data.
SlitScan()
"""
from warp import *
slitscanner_version = "$Id: slitscanner.py,v 1.2 2002/09/27 22:13:51 dave Exp $"

def slitscannerdoc():
  import slitscanner
  print slitscanner.__doc__

#---------------------------------------------------------------------------
class SlitData:
  """
Container for data a single slit
  """
  def __init__(self,x,y,z,xp,yp,vz,gi,s,center,width):
    if type(s) == StringType: s = eval(s,locals())
    self.center = center
    self.width = width
    s1 = center - width/2.
    s2 = center + width/2.
    ii = compress(logical_and(s1 < s,s < s2),arange(len(s)))
    self.x = take(x,ii)
    self.y = take(y,ii)
    self.z = take(z,ii)
    self.xp = take(xp,ii)
    self.yp = take(yp,ii)
    self.vz = take(vz,ii)
    self.gi = take(gi,ii)
  def __len__(self):
    return len(self.x)
  def propagate(self,zdist):
    """
Propagate the data some specified distance.
Eventually, this would optionally invoke the full Warp code to do a
self-consistent advance.
    """
    self.xold = self.x + 0.
    self.yold = self.y + 0.
    self.zold = self.z + 0.
    if len(self) > 0:
      dt = zdist/self.vz
      self.x = self.xold + self.xp*self.vz*dt
      self.y = self.yold + self.yp*self.vz*dt
      self.z = self.zold + zdist
  def ppxy(self):
    plp(self.y,self.x)

#---------------------------------------------------------------------------
class Slit1Data:
  """
Container for data from the first slit scan
  """
  def __init__(self,x,y,z,xp,yp,vz,gi,s,center,range,n,width):
    if type(s) == StringType: s = eval(s,locals())
    self.center = center
    self.range = range
    self.n = n
    self.d = range/(n-1.)
    self.width = width
    slist = center + range*(arange(n)/(n-1.) - 0.5)
    self.data = []
    for ss in slist:
      self.data.append(SlitData(x,y,z,xp,yp,vz,gi,s,ss,width))
  def __len__(self):
    return len(self.data)
  def propagate(self,zdist):
    for d in self.data: d.propagate(zdist)
  def setdensity(self,grid):
    grid[:] = map(len,self.data)
  def ppxy(self):
    for d in self.data: d.ppxy()
    
#---------------------------------------------------------------------------
class Slit2Data:
  """
Container for data from the second slit scan
  """
  def __init__(self,slit1data,s,center,range,n,width,slope,lcrossed,zdist):
    self.slit1data = slit1data
    self.center2 = center
    self.range2 = range
    self.n2 = n
    self.width = width
    self.slope = slope
    self.lcrossed = lcrossed
    self.data = []
    for d in slit1data.data:
      sc = center + slope*d.center
      sr = range
      if not lcrossed:
        sc = d.center + sc*zdist
        sr = range*zdist
      self.data.append(Slit1Data(d.x,d.y,d.z,d.xp,d.yp,d.vz,d.gi,
                                 s,sc,sr,n,width))
    self.setdensity()
    self.setmeshes()
  def __len__(self):
    return len(self.data)
  def setdensity(self):
    self.grid = zeros((len(self.slit1data),self.n2),'d')
    for i in range(len(self)):
      self.data[i].setdensity(self.grid[i,:])
  def setmeshes(self):
    xmesh = []
    ymesh = self.center2 + self.range2*(arange(self.n2)/(self.n2-1.) - 0.5)
    ymeshslope = []
    for d1 in self.slit1data.data:
      xmesh.append(d1.center)
      ymeshslope.append(ymesh + self.slope*d1.center)
    self.xmesh = transpose(array(self.n2*[xmesh]))
    self.ymesh = array(len(self.slit1data)*[ymesh])
    self.ymeshslope = array(ymeshslope)
  def ppxy(self):
    for d in self.data: d.ppxy()
  def pplines(self,scale,slope,color=color):
    if scale is None: scale = self.slit1data.d/maxnd(self.grid)*0.9
    x0 = self.xmesh
    x1 = x0 + scale*self.grid
    if slope == 'a': y = self.ymesh
    else:            y = self.ymeshslope
    pldj(x0,y,x1,y,color=color)

#---------------------------------------------------------------------------
class SlitScan:
  """
Generates slit scan data for slice particles.
  center1: center position of first slit
  range1: range of movement of center of first slit
  n1: number of slit locations for the first slit
  width1: width of slit
  s1: plane of first slit, either 'x' or 'y'
  center2: center position of second slit (see info below)
  range2: range of movement of center of second slit (see info below)
  n2: number of slit locations for the second slit
  width2: width of the second slit
  s2: plane of the second slit, either 'x' or 'y'
  slope: slope for second slit (see info below)
  zdist: distance between the two slits
  js=0: particle species to use
The input quantities are best thought of in terms of where the data should be
in the phase-space rather than the actual location of the slits in the experiment.
For the second slit, the center and range are given relative to the quantity
being measured. So for parallel slits (measuring phase space), the center and
range are given in units of radians. Center2 gives the location of the
center of the second slit for the point at the center of the first slit. For
other points, the center is given by product of slope and the first slit
position (plus center2), all multiplied by zdist.
  """

  def __init__(self,center1,range1,n1,width1,s1,
                    center2,range2,n2,width2,s2,slope,
                    zdist,js=0):

    # --- First, save slit parameters
    self.center1 = center1
    self.range1 = range1
    self.n1 = n1
    self.width1 = width1
    self.center2 = center2
    self.range2 = range2
    self.n2 = n2
    self.width2 = width2
    self.slope = slope
    self.zdist = zdist
    self.js = js

    # --- Check if slits are crossed, i.e. s1 != s2
    if s1 != s2: self.lcrossed = true
    else:        self.lcrossed = false

    # --- Get particle data
    x = getx(js=self.js)
    y = gety(js=self.js)
    z = getz(js=self.js)
    xp = getxp(js=self.js)
    yp = getyp(js=self.js)
    vz = getvz(js=self.js)
    gi = getgaminv(js=self.js)

    # --- Gather data from first slit
    self.slit1 = Slit1Data(x,y,z,xp,yp,vz,gi,s1,center1,range1,n1,width1)

    # --- Propagate data from 1st slit to 2nd slit
    self.slit1.propagate(zdist)

    # --- Gather data for second slit
    self.slit2 = Slit2Data(self.slit1,s2,center2,range2,n2,width2,slope,
                           self.lcrossed,zdist)

#----------------------------------------------------------------------------
  def pp(self,slope=None,**kw):
    """
Plots the data. Currently only with the unshearing. Accepts all arguments to
ppgeneric related to plotting gridded data.
    """
    kw["grid"] = self.slit2.grid
    kw["xmin"] = self.center1 - self.range1/2.
    kw["xmax"] = self.center1 + self.range1/2.
    kw["ymin"] = self.center2 - self.range2/2.
    kw["ymax"] = self.center2 + self.range2/2.
    if slope is not None:
      kw['xmesh'] = self.slit2.xmesh
      kw['ymesh'] = self.slit2.ymeshslope
    ppgeneric(kwdict=kw)

#----------------------------------------------------------------------------
  def pplines(self,scale=None,slope='a',color='fg'):
    """
Plots old style line plot of data.
    """
    self.slit2.pplines(scale,slope,color=color)

#----------------------------------------------------------------------------
