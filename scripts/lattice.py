from warp import *
import Ranf
import __main__
lattice_version = "$Id: lattice.py,v 1.1 2000/10/16 18:34:19 dave Exp $"

# Setup classes for MAD style input
# This includes both the elements from hibeam and WARP
# hibeam element names are all lower case, WARP elements are capitalized
# DPG 7/27/99
#import re
#import string
import RandomArray

# --- This function returns a random number from the given distribution.
_errordist_getnextnumber = 0
def errordist(etype):
  global _errordist_getnextnumber,_errordist_next_number
  if etype == 'GAUSSIAN':
    if _errordist_getnextnumber:
      _errordist_getnextnumber = 0
      return _errordist_next_number
    s = RandomArray.random(2)
    phi = 2.*pi*RandomArray.random(2)
    sq = sqrt(-2.*log(s))
    result = tuple(sq*cos(phi))
    _errordist_getnextnumber = 1
    _errordist_next_number = result[1]
    return result[0]
  elif etype == 'UNIFORM':
    return Ranf.ranf()
  elif etype == 'ABSOLUTE':
    return 1.
  else:
    return 1.

#############################################################################
# --- LINE contains a list of lattice elements.
class LINE:
  """
Creates and instance of the LINE lattice type which contains a list of
lattice elements. The argument can either be a single element or a list of
elements.
  """
  def __init__(self,*elems):
    self.type = 'LINE'
    # --- Unravel any imbedded lists.
    i = 0
    elems = list(elems)
    while i < len(elems):
      if type(elems[i]) == type([]):
        elems[i:i+1] = elems[i]
      else:
        i = i + 1
    # --- save list of elements
    self.elems = elems
    # --- initialize expanded list as a blank
    self.elemslist = []
  def expand(self):
    if self.elemslist: return self.elemslist
    for e in self.elems:
      self.elemslist.append(e.expand())
    # --- Unravel any imbedded lists.
    i = 0
    self.elemslist = list(self.elemslist)
    while i < len(self.elemslist):
      if type(self.elemslist[i]) == type([]):
        self.elemslist[i:i+1] = self.elemslist[i]
      else:
        i = i + 1
    return self.elemslist
  def __mul__(self,other):
    # --- This allows multiplication of elements by integers
    return LINE(other*[self])
  def __rmul__(self,other):
    # --- This allows multiplication of elements by integers
    return LINE(other*[self])
  def __add__(self,other):
    # --- This allows addition of elements
    return LINE(self,other)
  def __radd__(self,other):
    # --- This allows addition of elements
    return LINE(other,self)

# --- Create an equivalent class to LINE.
Line = LINE

# --- Child can be any type of lattice element with one or more
# --- parameters changed from the parent.
class child:
  """
Creates an instance of the child lattice type.  A child can be any type of
lattice element with one or more parameters changed from the parent. The first
argument, parent, is the parent element. The remaining elements are a keyword
list of the parameters which are changed.
  """
  def __init__(self,parent,**changes):
    if type(parent)==type(''):
      self.parent = eval(parent,__main__.__dict__,globals())
    else:
      self.parent = parent
    self.__dict__.update(changes)
    self.derivedquantities(self)
  def expand(self):
    return self
  def __getattr__(self,name):
    return eval('self.parent.'+name,locals())
  def __mul__(self,other):
    return LINE(other*[self])
  def __rmul__(self,other):
    return LINE(other*[self])
  def __add__(self,other):
    return LINE(self,other)
  def __radd__(self,other):
    return LINE(other,self)

# --- Create an equivalent class to child
Child = child

# --- Base element class. All elements have the following attributes:
# --- Length, aperture, X offset, Y offset, and error type for the offset.
class Elem:
  """
Base class for the lattice classes. Should never be directly called.
  """
  def __init__(self,l=0,length=0,zshift=0,zs=0,ze=0,aperture=0,
               offset_x=0,offset_y=0,error_type=''):
    self.l = l
    self.length = length
    self.zs = zs
    self.ze = ze
    self.type = ''
    self.zshift = zshift
    self.aperture = aperture
    self.offset_x = offset_x
    self.offset_y = offset_y
    self.error_type = error_type
    Elem.derivedquantities(self,self)
  def derivedquantities(_self,self):
    if self.l or self.length:
      self.length = max(self.l,self.length)
      self.zs = 0.
      self.ze = 0.
    elif self.zs or self.ze:
      self.length = self.ze - self.zs
  def isin(self,z=top.zbeam):
    if self.zs < z and z < self.ze: return 1
    return 0
  def expand(self):
    return self
  def __mul__(self,other):
    return LINE(other*[self])
  def __rmul__(self,other):
    return LINE(other*[self])
  def __add__(self,other):
    return LINE(self,other)
  def __radd__(self,other):
    return LINE(other,self)

#----------------------------------------------------------------------------
# HIBEAM elements
#----------------------------------------------------------------------------
class drift(Elem):
  def __init__(self,l=0,length=0,zshift=0,zs=0,ze=0,aperture=0,error_type='',
               offset_x=0,offset_y=0,i_cap_pointer=0,n_cap_nodes=0):
    Elem.__init__(self,l=l,length=length,zshift=zshift,zs=zs,ze=ze,
                  aperture=aperture,
                  offset_x=offset_x,offset_y=offset_y,error_type=error_type)
    self.type = 'drift'
    self.i_cap_pointer = i_cap_pointer
    self.n_cap_nodes = n_cap_nodes

class box(Elem):
  def __init__(self,l=0,length=0,zshift=0,zs=0,ze=0,aperture=0,
               width_x=0,width_y=0,error_type='',
               offset_x=0,offset_y=0,i_cap_pointer=0,n_cap_nodes=0):
    Elem.__init__(self,l=l,length=length,zshift=zshift,zs=zs,ze=ze,
                  aperture=aperture,
                  offset_x=offset_x,offset_y=offset_y,error_type=error_type)
    self.type = 'box'
    self.aperture = aperture
    self.width_x = width_x
    self.width_y = width_y
    self.i_cap_pointer = i_cap_pointer
    self.n_cap_nodes = n_cap_nodes
    self.derivedquantities(self)
  def derivedquantities(_self,self):
    if aperture:
      self.width_x = aperture
      self.width_y = aperture

class quad(Elem):
  def __init__(self,l=0,length=0,zshift=0,zs=0,ze=0,aperture=0,
               voltage=0,gradient=0,r_elem=0,
               offset_x=0,offset_y=0,error_type='',
               i_cap_pointer=0,n_cap_nodes=0):
    Elem.__init__(self,l=l,length=length,zshift=zshift,zs=zs,ze=ze,
                  aperture=aperture,
                  offset_x=offset_x,offset_y=offset_y,error_type=error_type)
    self.type = 'quad'
    self.voltage = voltage
    self.gradient = gradient
    self.r_elem = r_elem
    self.i_cap_pointer = i_cap_pointer
    self.n_cap_nodes = n_cap_nodes
    self.voltage = voltage
    self.gradient = gradient
    self.derivedquantities(self)
  def derivedquantities(_self,self):
    if self.gradient != 0:
      self.voltage = self.gradient*self.aperture**2
    else:
      self.gradient = self.voltage/self.aperture**2

class hyperb(Elem):
  def __init__(self,l=0,length=0,zshift=0,zs=0,ze=0,aperture=0,voltage=0,
               r_elem=0,offset_x=0,offset_y=0,error_type='',
               i_cap_pointer=0,n_cap_nodes=0):
    Elem.__init__(self,l=l,length=length,zshift=zshift,zs=zs,ze=ze,
                  aperture=aperture,
                  offset_x=offset_x,offset_y=offset_y,error_type=error_type)
    self.type = 'hyperb'
    self.voltage = voltage
    self.r_elem = r_elem
    self.i_cap_pointer = i_cap_pointer
    self.n_cap_nodes = n_cap_nodes

class wire(Elem):
  def __init__(self,l=0,length=0,zshift=0,zs=0,ze=0,aperture=0,
               offset_x=0,offset_y=0,error_type=''):
    Elem.__init__(self,l=l,length=length,zshift=zshift,zs=zs,ze=ze,
                  aperture=aperture,
                  offset_x=offset_x,offset_y=offset_y,error_type=error_type)
    self.type = 'wire'

#----------------------------------------------------------------------------
# WARP elements
#----------------------------------------------------------------------------

class Drft(Elem):
  """
Creates an instance of a Drft lattice element.
  - l=0 drift length
  - length=0 alternate form of the drift length
  - zshift=0 start of element relative to current lattice location
  - ap=0 aperture (can affect location of transverse boundaries)
  - ox=0 offset in x (can affect location of transverse boundaries)
  - oy=0 offset in y (can affect location of transverse boundaries)
  - error_type='' type of error distribution to apply
  """
  def __init__(self,l=0,length=0,zshift=0,zs=0,ze=0,ap=0,ox=0,oy=0,
               error_type=''):
    Elem.__init__(self,l=l,length=length,zshift=zshift,zs=zs,ze=ze,aperture=ap,
                  offset_x=ox,offset_y=oy,error_type=error_type)
    self.type = 'Drft'

class Bend(Elem):
  """
Creates an instance of a Bend lattice element.
  - l=0 drift length
  - length=0 alternate form of the drift length
  - zshift=0 start of element relative to current lattice location
  - ap=0 aperture (can affect location of transverse boundaries)
  - ox=0 offset in x (can affect location of transverse boundaries)
  - oy=0 offset in y (can affect location of transverse boundaries)
  - error_type='' type of error distribution to apply
  - rc=1.e36 radius of curvature of the bend
  """
  def __init__(self,l=0,length=0,zshift=0,zs=0,ze=0,ap=0,ox=0,oy=0,
               error_type='',
               rc=1.e36):
    Elem.__init__(self,l=l,length=length,zshift=zshift,zs=zs,ze=ze,aperture=ap,
                  offset_x=ox,offset_y=oy,error_type=error_type)
    self.type = 'Bend'
    self.rc = rc

class Dipo(Elem):
  """
Creates an instance of a Dipo lattice element.
  - l=0 drift length
  - length=0 alternate form of the drift length
  - zshift=0 start of element relative to current lattice location
  - ap=0 aperture (can affect location of transverse boundaries)
  - ox=0 offset in x (can affect location of transverse boundaries)
  - oy=0 offset in y (can affect location of transverse boundaries)
  - error_type='' type of error distribution to apply
  - ex=0 Ex field (V/m)
  - ey=0 Ey field (V/m)
  - bx=0 Bx field (T/m)
  - by=0 By field (T/m)
  - ta=0 tangent of entrance angle
  - tb=0 tangent of exit angle
  - x1=0 location of first electric plate
  - x2=0 location of second electric plate
  - v1=0 voltage on first electric plate
  - v2=0 voltage on second electric plate
  - l1=0 length of first electric plate
  - l2=0 length of second electric plate
  - w1=0 width of first electric plate
  - w2=0 width of second electric plate
  """
  def __init__(self,l=0,length=0,zshift=0,zs=0,ze=0,ap=0,ox=0,oy=0,
               error_type='',
               ex=0,ey=0,bx=0,by=0,ta=0,tb=0,
               x1=0,x2=0,v1=0,v2=0,l1=0,l2=0,w1=0,w2=0):
    Elem.__init__(self,l=l,length=length,zshift=zshift,zs=zs,ze=ze,aperture=ap,
                  offset_x=ox,offset_y=oy,error_type=error_type)
    self.type = 'Dipo'
    self.ex = ex
    self.ey = ey
    self.bx = bx
    self.by = by
    self.ta = ta
    self.tb = tb
    self.x1 = x1
    self.x2 = x2
    self.v1 = v1
    self.v2 = v2
    self.l1 = l1
    self.l2 = l2
    self.w1 = w1
    self.w2 = w2

class Quad(Elem):
  """
Creates an instance of a Quad lattice element.
  - l=0 drift length
  - length=0 alternate form of the drift length
  - zshift=0 start of element relative to current lattice location
  - ap=0 aperture (can affect location of transverse boundaries)
  - ox=0 offset in x (can affect location of transverse boundaries)
  - oy=0 offset in y (can affect location of transverse boundaries)
  - error_type='' type of error distribution to apply
  - de=0 electric field gradient (V/m**2)
  - db=0 magnetic field gradient (T/m**2)
  - vx=0 voltage on x-rod of an electric quad
  - vy=0 voltage on y-rod of an electric quad
  - rr=0 rod radius of an electric quad
  - rl=0 rod length of an electric quad
  - gl=0 gap length between rod and end plate of an electric quad
  - gp=0 position of the rod to plate gap in the x plane of an electric quad
  - pw=0 plate with of an electric quad
  - pa=0 plate aperture of an electric quad
  - pr=0 plate outer radius of an electric quad
  """
  def __init__(self,l=0,length=0,zshift=0,zs=0,ze=0,ap=0,ox=0,oy=0,
               error_type='',
               de=0,db=0,vx=0,vy=0,rr=0,rl=0,gl=0,gp=0,pw=0,pa=0,pr=0):
    Elem.__init__(self,l=l,length=length,zshift=zshift,zs=zs,ze=ze,aperture=ap,
                  offset_x=ox,offset_y=oy,error_type=error_type)
    self.type = 'Quad'
    self.de = de
    self.db = db
    self.vx = vx
    self.vy = vy
    self.rr = rr
    self.rl = rl
    self.gl = gl
    self.gp = gp
    self.pw = pw
    self.pa = pa
    self.pr = pr

class Sext(Elem):
  """
Creates an instance of a Sext lattice element.
  - l=0 drift length
  - length=0 alternate form of the drift length
  - zshift=0 start of element relative to current lattice location
  - ap=0 aperture (can affect location of transverse boundaries)
  - ox=0 offset in x (can affect location of transverse boundaries)
  - oy=0 offset in y (can affect location of transverse boundaries)
  - error_type='' type of error distribution to apply
  - de=0 electric field gradient
  - db=0 magnetic field gradient
  """
  def __init__(self,l=0,length=0,zshift=0,zs=0,ze=0,ap=0,ox=0,oy=0,
               error_type='',
               de=0,db=0):
    Elem.__init__(self,l=l,length=length,zshift=zshift,zs=zs,ze=ze,aperture=ap,
                  offset_x=ox,offset_y=oy,error_type=error_type)
    self.type = 'Sext'
    self.de = de
    self.db = db

class Hele(Elem):
  """
Creates an instance of a Hele lattice element.
  - l=0 drift length
  - length=0 alternate form of the drift length
  - zshift=0 start of element relative to current lattice location
  - ap=0 aperture (can affect location of transverse boundaries)
  - ox=0 offset in x (can affect location of transverse boundaries)
  - oy=0 offset in y (can affect location of transverse boundaries)
  - error_type='' type of error distribution to apply
  - nn=[] list of n indices of multipole components
  - vv=[] list of v indices of multipole components
  - ae=[] list of magnitudes of electric multipole components
  - am=[] list of magnitudes of magnetic multipole components
  - pe=[] list of phase angles of electric multipole components
  - pm=[] list of phase angles of magnetic multipole components
  """
  def __init__(self,l=0,length=0,zshift=0,zs=0,ze=0,ap=0,ox=0,oy=0,
               error_type='',
               nn=[],vv=[],ae=[],am=[],pe=[],pm=[]):
    Elem.__init__(self,l=l,length=length,zshift=zshift,zs=zs,ze=ze,aperture=ap,
                  offset_x=ox,offset_y=oy,error_type=error_type)
    self.type = 'Hele'
    self.nn = nn
    self.vv = vv
    self.ae = ae
    self.am = am
    self.pe = pe
    self.pm = pm
    self.derivedquantities(self)
  def derivedquantities(_self,self):
    self.nn = array(self.nn)
    self.vv = array(self.vv)
    self.ae = array(self.ae)
    self.am = array(self.am)
    self.pe = array(self.pe)
    self.pm = array(self.pm)

class Accl(Elem):
  """
Creates an instance of a Accl lattice element.
  - l=0 drift length
  - length=0 alternate form of the drift length
  - zshift=0 start of element relative to current lattice location
  - ap=0 aperture (can affect location of transverse boundaries)
  - ox=0 offset in x (can affect location of transverse boundaries)
  - oy=0 offset in y (can affect location of transverse boundaries)
  - error_type='' type of error distribution to apply
  - ez=0 time independent accelerating gradient (V/m)
  - xw=0 x variation in accelerating gradient
  - sw=0 switch sets whether mesh frame is accelerated by a gap
  - et=[] time dependent accelerating gradient (1-D array)
  - ts=0 initial time of time dependent accelerating gradient
  - dt=0 time increment in the time dependent accelerating gradient data (s)
  """
  def __init__(self,l=0,length=0,zshift=0,zs=0,ze=0,ap=0,ox=0,oy=0,
               error_type='',
               ez=0,xw=0,sw=0,et=[],ts=0,dt=0):
    Elem.__init__(self,l=l,length=length,zshift=zshift,zs=zs,ze=ze,aperture=ap,
                  offset_x=ox,offset_y=oy,error_type=error_type)
    self.type = 'Accl'
    self.ez = ez
    self.xw = xw
    self.sw = sw
    self.et = et
    self.ts = ts
    self.dt = dt
    self.derivedquantities(self)
  def derivedquantities(_self,self):
    self.et = array(self.et)

class Emlt(Elem):
  """
Creates an instance of a Emlt lattice element.
  - l=0 drift length
  - length=0 alternate form of the drift length
  - zshift=0 start of element relative to current lattice location
  - ap=0 aperture (can affect location of transverse boundaries)
  - ox=0 offset in x (can affect location of transverse boundaries)
  - oy=0 offset in y (can affect location of transverse boundaries)
  - error_type='' type of error distribution to apply
Either specify the index
  - id=0 index of the data set to use (if one is not supplied)
Or specify the data set
  - dz=0 z increment in the data set (z)
  - e=[] multipole components (1 or 2-D, 1st is data, 2nd is list of components)
  - ep=[] axial derivatives of multipole components
  - eph=[] phase angle
  - nn=[] n indices
  - vv=[] v indices
  - ph=0 overal phase angle
  - sf=0 relative scale factor
  - sc=1 absolute scale factor
  """
  def __init__(self,l=0,length=0,zshift=0,zs=0,ze=0,ap=0,ox=0,oy=0,
               error_type='',
               id=0,dz=0,e=[],ep=[],eph=[],nn=[],vv=[],ph=0,sf=0,sc=1):
    Elem.__init__(self,l=l,length=length,zshift=zshift,zs=zs,ze=ze,aperture=ap,
                  offset_x=ox,offset_y=oy,error_type=error_type)
    self.type = 'Emlt'
    self.id = id
    self.dz = dz
    self.e = e
    self.ep = ep
    self.eph = eph
    self.nn = nn
    self.vv = vv
    self.ph = ph
    self.sf = sf
    self.sc = sc
    self.derivedquantities(self)
  def derivedquantities(_self,self):
    self.e  = array(self.e)
    self.ep = array(self.ep)
    self.eph = array(self.eph)
    if len(shape(self.e)) == 1: self.e = self.e[:,NewAxis]
    if len(shape(self.ep)) == 1: self.ep = self.ep[:,NewAxis]
    if len(shape(self.eph)) == 1: self.eph = self.eph[:,NewAxis]
    self.nn = array(self.nn)
    self.vv = array(self.vv)

class Mmlt(Elem):
  """
Creates an instance of a Mmlt lattice element.
  - l=0 drift length
  - length=0 alternate form of the drift length
  - zshift=0 start of element relative to current lattice location
  - ap=0 aperture (can affect location of transverse boundaries)
  - ox=0 offset in x (can affect location of transverse boundaries)
  - oy=0 offset in y (can affect location of transverse boundaries)
  - error_type='' type of error distribution to apply
Either specify the index
  - id=0 index of the data set to use (if a data set is not supplied)
Or specify the data set
  - dz=0 z increment in the data set (z)
  - m=[] multipole components (1 or 2-D, 1st is data, 2nd is list of components)
  - mp=[] axial derivatives of multipole components
  - mph=[] phase angle
  - nn=[] n indices
  - vv=[] v indices
  - ph=0 overal phase angle
  - sf=0 relative scale factor
  - sc=1 absolute scale factor
  """
  def __init__(self,l=0,length=0,zshift=0,zs=0,ze=0,ap=0,ox=0,oy=0,
               error_type='',
               id=0,dz=0,m=[],mp=[],mph=[],nn=[],vv=[],ph=0,sf=0,sc=1):
    Elem.__init__(self,l=l,length=length,zshift=zshift,zs=zs,ze=ze,aperture=ap,
                  offset_x=ox,offset_y=oy,error_type=error_type)
    self.type = 'Mmlt'
    self.id = id
    self.dz = dz
    self.m = m
    self.mp = mp
    self.mph = mph
    self.nn = nn
    self.vv = vv
    self.ph = ph
    self.sf = sf
    self.sc = sc
    self.derivedquantities(self)
  def derivedquantities(_self,self):
    self.m  = array(self.m)
    self.mp = array(self.mp)
    self.mph = array(self.mph)
    if len(shape(self.m)) == 1: self.m = self.m[:,NewAxis]
    if len(shape(self.mp)) == 1: self.mp = self.mp[:,NewAxis]
    if len(shape(self.mph)) == 1: self.mph = self.mph[:,NewAxis]
    self.nn = array(self.nn)
    self.vv = array(self.vv)

class Bgrd(Elem):
  """
Creates an instance of a Bgrd lattice element.
  - l=0 drift length
  - length=0 alternate form of the drift length
  - zshift=0 start of element relative to current lattice location
  - ap=0 aperture (can affect location of transverse boundaries)
  - ox=0 offset in x (can affect location of transverse boundaries)
  - oy=0 offset in y (can affect location of transverse boundaries)
  - error_type='' type of error distribution to apply
Either specify the index
  - id=0 index of the data set to use (if one is not supplied)
Or specify the data set
  - bx=[] Bx (Tesla)
  - by=[] By (Tesla)
  - bz=[] Bz (Tesla)
  - dx=0 x increment size (m)
  - dy=0 y increment size (m)
  - dz=0 z increment size (m)
  - sf=0 relative scaling factor
  - sc=0 absolute scaling factor
  """
  def __init__(self,l=0,length=0,zshift=0,zs=0,ze=0,ap=0,ox=0,oy=0,
               error_type='',
               id=0,sf=0,sc=1,bx=[],by=[],bz=[],dx=0,dy=0,dz=0):
    Elem.__init__(self,l=l,length=length,zshift=zshift,zs=zs,ze=ze,aperture=ap,
                  offset_x=ox,offset_y=oy,error_type=error_type)
    self.type = 'Bgrd'
    self.id = id
    self.sf = sf
    self.sc = sc
    self.bx = bx
    self.by = by
    self.bz = bz
    self.dx = dx
    self.dy = dy
    self.dz = dz
    self.derivedquantities(self)
  def derivedquantities(_self,self):
    self.bx = array(self.bx)
    self.by = array(self.by)
    self.bz = array(self.bz)

class Pgrd(Elem):
  """
Creates an instance of a Pgrd lattice element.
  - l=0 drift length
  - length=0 alternate form of the drift length
  - zshift=0 start of element relative to current lattice location
  - ap=0 aperture (can affect location of transverse boundaries)
  - ox=0 offset in x (can affect location of transverse boundaries)
  - oy=0 offset in y (can affect location of transverse boundaries)
  - error_type='' type of error distribution to apply
Either specify the index
  - id=0 index of the data set to use (if one is not supplied)
Or specify the data set
  - pp=[] electrostatic potential (Volts)
  - dx=0 x increment size (m)
  - dy=0 y increment size (m)
  - dz=0 z increment size (m)
  - sf=0 relative scaling factor
  - sc=0 absolute scaling factor
  """
  def __init__(self,l=0,length=0,zshift=0,zs=0,ze=0,ap=0,ox=0,oy=0,
               error_type='',
               id=0,sf=0,sc=1,pp=[],dx=0,dy=0,dz=0):
    Elem.__init__(self,l=l,length=length,zshift=zshift,zs=zs,ze=ze,aperture=ap,
                  offset_x=ox,offset_y=oy,error_type=error_type)
    self.type = 'Pgrd'
    self.id = id
    self.sf = sf
    self.sc = sc
    self.pp = pp
    self.dx = dx
    self.dy = dy
    self.dz = dz
    self.derivedquantities(self)
  def derivedquantities(_self,self):
    self.pp = array(self.pp)

############################################################################
############################################################################
# --- Convert and copy the input MAD lattice into a WARP lattice.
def madtowarp(line):
  """
Takes the inputted line obect and sets up the WARP lattice. Any existing
information if the WARP lattice arrays is deleted.
  """
  lattice = line.expand()
  idrft = -1
  iquad = -1
  ibend = -1
  idipo = -1
  isext = -1
  ihele = -1
  ihmlt = -1
  iaccl = -1
  iemlt = -1
  immlt = -1
  ibgrd = -1
  ipgrd = -1
  top.ndrft = 0
  top.nbend = 0
  top.ndipo = 0
  top.nquad = 0
  top.nsext = 0
  top.nhele = 0
  top.nhmlt = 0
  top.nqerr = 0
  top.naccl = 0
  top.ntaccl = 0
  top.nemlt = 0
  top.nmmlt = 0
  top.neerr = 0
  top.nmerr = 0
  top.nbgrd = 0
  top.npgrd = 0
  top.nemltsets = 0
  top.nesmult = 0
  top.nzemltmax = 0
  top.nmmltsets = 0
  top.nmsmult = 0
  top.nzmmltmax = 0
  top.bgrdnx = 0
  top.bgrdny = 0
  top.bgrdnz = 0
  top.bgrdns = 0
  top.pgrdnx = 0
  top.pgrdny = 0
  top.pgrdnz = 0
  top.pgrdns = 0

  # --- Loop through elements and count how many of each type there are.
  # --- This mush be done so that the arrays can be allocated to the correct
  # --- size.
  for e in lattice:
    if (e.type == 'Drft' or e.type == 'drift' or
        e.type == 'box' or e.type == 'wire'):
      top.ndrft = top.ndrft + 1
    elif e.type == 'quad' or e.type == 'hyperb' or e.type == 'Quad':
      top.nquad = top.nquad + 1
      top.nqerr = top.nqerr + 1
    elif e.type == 'Bend':
      top.nbend = top.nbend + 1
    elif e.type == 'Dipo':
      top.ndipo = top.ndipo + 1
    elif e.type == 'Sext':
      top.nsext = top.nsext + 1
    elif e.type == 'Hele':
      top.nhele = top.nhele + 1
    elif e.type == 'Accl':
      top.naccl = top.naccl + 1
    elif e.type == 'Emlt':
      top.nemlt = top.nemlt + 1
      top.neerr = top.neerr + 1
    elif e.type == 'Mmlt':
      top.nmmlt = top.nmmlt + 1
      top.nmerr = top.nmerr + 1
    elif e.type == 'Bgrd':
      top.nbgrd = top.nbgrd + 1
    elif e.type == 'Pgrd':
      top.npgrd = top.npgrd + 1
  
  # --- Allocate the data space needed.
  gchange("Lattice")

  # --- Loop through elements and copy data to WARP variables.
  # --- The zs's and ze's must be calculated on the fly (and not stored
  # --- in the MAD lattice elements.
  # --- Note that for elements with seperate data sets, only the number
  # --- of data points is checked. Arrays are then later allocated and filled.
  zz = 0
  for e in lattice:
    if e.type == 'quad':
      iquad = iquad + 1
      top.quadzs[iquad] = zz + e.zshift
      top.quadze[iquad] = top.quadzs[iquad] + e.length
      zz = top.quadze[iquad]
      top.quadde[iquad] = e.gradient
      top.quadap[iquad] = e.aperture
      top.quadrr[iquad] = e.r_elem
      top.qoffx[iquad] = e.offset_x*errordist(e.error_type)
      top.qoffy[iquad] = e.offset_y*errordist(e.error_type)
    elif e.type == 'hyperb':
      iquad = iquad + 1
      top.quadzs[iquad] = zz + e.zshift
      top.quadze[iquad] = top.quadzs[iquad] + e.length
      zz = top.quadze[iquad]
      top.qoffx[iquad] = e.offset_x*errordist(e.error_type)
      top.qoffy[iquad] = e.offset_y*errordist(e.error_type)
      top.quadde[iquad] = e.gradient
      top.quadap[iquad] = -e.aperture
    elif e.type == 'drift':
      idrft = idrft + 1
      top.drftzs[idrft] = zz + e.zshift
      top.drftze[idrft] = top.drftzs[idrft] + e.length
      zz = top.drftze[idrft]
      top.drftap[idrft] = e.aperture
      top.drftox[idrft] = e.offset_x*errordist(e.error_type)
      top.drftoy[idrft] = e.offset_y*errordist(e.error_type)
    elif e.type == 'box':
      idrft = idrft + 1
      top.drftzs[idrft] = zz + e.zshift
      top.drftze[idrft] = top.drftzs[idrft] + e.length
      zz = top.drftze[idrft]
      top.drftap[idrft] = -e.aperture
      top.drftox[idrft] = e.offset_x*errordist(e.error_type)
      top.drftoy[idrft] = e.offset_y*errordist(e.error_type)
    elif e.type == 'wire':
      idrft = idrft + 1
      top.drftzs[idrft] = zz + e.zshift
      top.drftze[idrft] = top.drftzs[idrft] + e.length
      zz = top.drftze[idrft]
      top.drftap[idrft] = e.aperture
      top.drftox[idrft] = e.offset_x*errordist(e.error_type)
      top.drftoy[idrft] = e.offset_y*errordist(e.error_type)

    elif e.type == 'Quad':
      iquad = iquad + 1
      top.quadzs[iquad] = zz + e.zshift
      top.quadze[iquad] = top.quadzs[iquad] + e.length
      zz = top.quadze[iquad]
      top.quadap[iquad] = e.aperture
      top.qoffx[iquad] = e.offset_x*errordist(e.error_type)
      top.qoffy[iquad] = e.offset_y*errordist(e.error_type)
      top.quadde[iquad] = e.de
      top.quaddb[iquad] = e.db
      top.quadvx[iquad] = e.vx
      top.quadvy[iquad] = e.vy
      top.quadrr[iquad] = e.rr
      top.quadrl[iquad] = e.rl
      top.quadgl[iquad] = e.gl
      top.quadgp[iquad] = e.gp
      top.quadpw[iquad] = e.pw
      top.quadpa[iquad] = e.pa
      top.quadpr[iquad] = e.pr

    elif e.type == 'Drft':
      idrft = idrft + 1
      top.drftzs[idrft] = zz + e.zshift
      top.drftze[idrft] = top.drftzs[idrft] + e.length
      zz = top.drftze[idrft]
      top.drftap[idrft] = e.aperture
      top.drftox[idrft] = e.offset_x*errordist(e.error_type)
      top.drftoy[idrft] = e.offset_y*errordist(e.error_type)

    elif e.type == 'Bend':
      ibend = ibend + 1
      top.bendzs[ibend] = zz + e.zshift
      top.bendze[ibend] = top.bendzs[ibend] + e.length
      zz = top.bendze[ibend]
      top.bendap[ibend] = e.aperture
      #top.bendox[ibend] = e.offset_x*errordist(e.error_type)
      #top.bendoy[ibend] = e.offset_y*errordist(e.error_type)
      top.bendrc[ibend] = e.rc

    elif e.type == 'Dipo':
      idipo = idipo + 1
      top.dipozs[idipo] = zz + e.zshift
      top.dipoze[idipo] = top.dipozs[idipo] + e.length
      zz = top.dipoze[idipo]
      top.dipoap[idipo] = e.aperture
      #top.dipoox[idipo] = e.offset_x*errordist(e.error_type)
      #top.dipooy[idipo] = e.offset_y*errordist(e.error_type)
      top.dipoex[idipo] = e.ex
      top.dipoey[idipo] = e.ey
      top.dipobx[idipo] = e.bx
      top.dipoby[idipo] = e.by
      top.dipota[idipo] = e.ta
      top.dipotb[idipo] = e.tb
      top.dipox1[idipo] = e.x1
      top.dipox2[idipo] = e.x2
      top.dipov1[idipo] = e.v1
      top.dipov2[idipo] = e.v2
      top.dipol1[idipo] = e.l1
      top.dipol2[idipo] = e.l2
      top.dipow1[idipo] = e.w1
      top.dipow2[idipo] = e.w2

    elif e.type == 'Sext':
      isext = isext + 1
      top.sextzs[isext] = zz + e.zshift
      top.sextze[isext] = top.sextzs[isext] + e.length
      zz = top.sextze[isext]
      top.sextap[isext] = e.aperture
      top.sextox[isext] = e.offset_x*errordist(e.error_type)
      top.sextoy[isext] = e.offset_y*errordist(e.error_type)
      top.sextde[isext] = e.de
      top.sextdb[isext] = e.db

    elif e.type == 'Hele':
      ihele = ihele + 1
      top.helezs[ihele] = zz + e.zshift
      top.heleze[ihele] = top.helezs[ihele] + e.length
      zz = top.heleze[ihele]
      top.heleap[ihele] = e.aperture
      top.heleox[ihele] = e.offset_x*errordist(e.error_type)
      top.heleoy[ihele] = e.offset_y*errordist(e.error_type)
      top.nhmlt = max(top.nhmlt,len(e.nn))

    elif e.type == 'Accl':
      iaccl = iaccl + 1
      top.acclzs[iaccl] = zz + e.zshift
      top.acclze[iaccl] = top.acclzs[iaccl] + e.length
      zz = top.acclze[iaccl]
      top.acclap[iaccl] = e.aperture
      top.acclox[iaccl] = e.offset_x*errordist(e.error_type)
      top.accloy[iaccl] = e.offset_y*errordist(e.error_type)
      top.acclez[iaccl] = e.ez
      top.acclxw[iaccl] = e.xw
      top.acclsw[iaccl] = e.sw
      top.acclts[iaccl] = e.ts
      top.accldt[iaccl] = e.dt
      top.ntaccl = max(top.ntaccl,shape(e.et)[0]-1)

    elif e.type == 'Emlt':
      iemlt = iemlt + 1
      top.emltzs[iemlt] = zz + e.zshift
      top.emltze[iemlt] = top.emltzs[iemlt] + e.length
      zz = top.emltze[iemlt]
      top.emltap[iemlt] = e.aperture
      top.emltox[iemlt] = e.offset_x*errordist(e.error_type)
      top.emltoy[iemlt] = e.offset_y*errordist(e.error_type)
      top.emltph[iemlt] = e.ph
      top.emltsf[iemlt] = e.sf
      top.emltsc[iemlt] = e.sc
      if not e.id:
        top.emltid[iemlt] = top.nemltsets + 1
      else:
        top.emltid[iemlt] = e.id
      top.nemltsets = max(top.emltid[iemlt],top.nemltsets)
      if e.e:
        top.nzemltmax = max(shape(e.e)[0]-1,top.nzemltmax)
        if len(shape(e.e)) == 2: top.nesmult = max(shape(e.e)[1],top.nesmult)
      if e.ep:
        top.nzemltmax = max(shape(e.ep)[0]-1,top.nzemltmax)
        if len(shape(e.ep)) == 2: top.nesmult = max(shape(e.ep)[1],top.nesmult)

    elif e.type == 'Mmlt':
      immlt = immlt + 1
      top.mmltzs[immlt] = zz + e.zshift
      top.mmltze[immlt] = top.mmltzs[immlt] + e.length
      zz = top.mmltze[immlt]
      top.mmltap[immlt] = e.aperture
      top.mmltox[immlt] = e.offset_x*errordist(e.error_type)
      top.mmltoy[immlt] = e.offset_y*errordist(e.error_type)
      top.mmltph[immlt] = e.ph
      top.mmltsf[immlt] = e.sf
      top.mmltsc[immlt] = e.sc
      if not e.id:
        top.mmltid[immlt] = top.nmmltsets + 1
      else:
        top.mmltid[immlt] = e.id
      top.nmmltsets = max(top.mmltid[immlt],top.nmmltsets)
      if e.m:
        top.nzmmltmax = max(shape(e.m)[0]-1,top.nzmmltmax)
        if len(shape(e.m)) == 2: top.nmsmult = max(shape(e.m)[1],top.nmsmult)
      if e.mp:
        top.nzmmltmax = max(shape(e.mp)[0]-1,top.nzmmltmax)
        if len(shape(e.mp)) == 2: top.nmsmult = max(shape(e.mp)[1],top.nmsmult)

    elif e.type == 'Bgrd':
      ibgrd = ibgrd + 1
      top.bgrdzs[ibgrd] = zz + e.zshift
      top.bgrdze[ibgrd] = top.bgrdzs[ibgrd] + e.length
      zz = top.bgrdze[ibgrd]
      top.bgrdap[ibgrd] = e.aperture
      top.bgrdox[ibgrd] = e.offset_x*errordist(e.error_type)
      top.bgrdoy[ibgrd] = e.offset_y*errordist(e.error_type)
      top.bgrdsf[ibgrd] = e.sf
      top.bgrdsc[ibgrd] = e.sc
      if not e.id:
        top.bgrdid[immlt] = top.bgrdns + 1
      else:
        top.bgrdid[immlt] = e.id
      top.bgrdns = max(top.bgrdid[immlt],top,bgrdns)
      if e.bx or e.by or e.bz:
        top.bgrdnx=max(shape(e.bx)[0],shape(e.by)[0],shape(e.bz)[0],top.bgrdnx)
        top.bgrdny=max(shape(e.bx)[1],shape(e.by)[1],shape(e.bz)[1],top.bgrdny)
        top.bgrdnz=max(shape(e.bx)[2],shape(e.by)[2],shape(e.bz)[2],top.bgrdnz)

    elif e.type == 'Pgrd':
      ipgrd = ipgrd + 1
      top.pgrdzs[ipgrd] = zz + e.zshift
      top.pgrdze[ipgrd] = top.pgrdzs[ipgrd] + e.length
      zz = top.pgrdze[ipgrd]
      top.pgrdap[ipgrd] = e.aperture
      top.pgrdox[ipgrd] = e.offset_x*errordist(e.error_type)
      top.pgrdoy[ipgrd] = e.offset_y*errordist(e.error_type)
      top.pgrdsf[ipgrd] = e.sf
      top.pgrdsc[ipgrd] = e.sc
      if not e.id:
        top.pgrdid[immlt] = top.pgrdns + 1
      else:
        top.pgrdid[immlt] = e.id
      top.pgrdns = max(top.pgrdid[immlt],top,pgrdns)
      if e.p:
        top.pgrdnx = max(shape(e.p)[0],top.pgrdnx)
        top.pgrdny = max(shape(e.p)[1],top.pgrdny)
        top.pgrdnz = max(shape(e.p)[2],top.pgrdnz)
  
  # --- Allocate the data space needed.
  gchange("Lattice")
  gchange("Mult_data")
  gchange("BGRDdata")
  gchange("PGRDdata")

  # --- Now, go back through the latice and fill in the data that needed
  # --- space allocated.
  iaccl = -1
  ihele = -1
  ihmlt = -1
  iemlt = -1
  immlt = -1
  ibgrd = -1
  ipgrd = -1
  for e in lattice:
    if e.type == 'Accl':
      iaccl = iaccl + 1
      top.acclet[:len(e.et),iaccl] = e.et

    elif e.type == 'Hele':
      ihele = ihele + 1
      top.hele_n[:,ihele] = e.nn
      top.hele_v[:,ihele] = e.vv
      top.heleae[:len(e.ae),ihele] = e.ae
      top.heleam[:len(e.am),ihele] = e.am
      top.helepe[:len(e.pe),ihele] = e.pe
      top.helepm[:len(e.pm),ihele] = e.pm

    elif e.type == 'Emlt':
      iemlt = iemlt + 1
      id = top.emltid[iemlt]
      if (e.e or e.ep) and top.nzemlt[id-1] == 0:
        # --- Only copy data if this is a new data set
        top.nzemlt[id-1] = max(shape(e.e)[0],shape(e.ep)[0]) - 1
        if not e.dz:
          top.dzemlt[id-1] = ((top.emltze[iemlt] - top.emltzs[iemlt])/
                              top.nzemlt[id-1])
        else:
          top.dzemlt[id-1] = e.dz
        top.emlt_n[:len(e.nn)] = e.nn
        top.emlt_v[:len(e.vv)] = e.vv
        top.esemlt[:shape(e.e)[0],:shape(e.e)[1],id-1] = e.e
        top.esemltp[:shape(e.ep)[0],:shape(e.ep)[1],id-1] = e.ep
        top.esemltph[:shape(e.eph)[0],:shape(e.eph)[1],id-1] = e.eph

    elif e.type == 'Mmlt':
      immlt = immlt + 1
      id = top.mmltid[iemlt]
      if (e.m or e.mp) and top.nzmmlt[id-1] == 0:
        # --- Only copy data if this is a new data set
        top.nzmmlt[id-1] = max(shape(e.m)[0],shape(e.mp)[0]) - 1
        if not e.dz:
          top.dzmmlt[id-1] = ((top.mmltze[immlt] - top.mmltzs[immlt])/
                              top.nzmmlt[id-1])
        else:
          top.dzmmlt[id-1] = e.dz
        top.mmlt_n[:len(e.nn)] = e.nn
        top.mmlt_v[:len(e.vv)] = e.vv
        top.msmmlt[:shape(e.m)[0],:shape(e.m)[1],id-1] = e.m
        top.msmmltp[:shape(e.mp)[0],:shape(e.mp)[1],id-1] = e.mp
        top.msmmltph[:shape(e.mph)[0],:shape(e.mph)[1],id-1] = e.mph

    elif e.type == 'Bgrd':
      ibgrd = ibgrd + 1
      id = top.bgrdid[iemlt]
      if (e.bx or e.by or e.bz) and top.bgrddx[id-1] == 0.:
        top.bgrddx[id-1] = e.dx
        top.bgrddy[id-1] = e.dy
        top.bgrddz[id-1] = e.dz
        sbx = shape(bx)
        sby = shape(by)
        sbz = shape(bz)
        top.bgrdbx[:sbx[0],:sbx[1],:sbx[2],id-1] = e.bx
        top.bgrdby[:sby[0],:sby[1],:sby[2],id-1] = e.by
        top.bgrdbz[:sbz[0],:sbz[1],:sbz[2],id-1] = e.bz

    elif e.type == 'Pgrd':
      ipgrd = ipgrd + 1
      id = top.pgrdid[iemlt]
      if e.p and top.pgrddx[id-1] == 0.:
        top.pgrddx[id-1] = e.dx
        top.pgrddy[id-1] = e.dy
        top.pgrddz[id-1] = e.dz
        sp = shape(p)
        top.pgrdp[:sp[0],:sp[1],:sp[2],id-1] = e.p

  # --- Finish by setting some general parameters.
  if iquad >= 2:
    top.tunelen = 0.5*(top.quadzs[2]+top.quadze[2]-top.quadzs[0]-top.quadze[0])
  if ihele >= 2:
    top.tunelen = 0.5*(top.helezs[2]+top.heleze[2]-top.helezs[0]-top.heleze[0])
  if iemlt >= 2:
    top.tunelen = 0.5*(top.emltzs[2]+top.emltze[2]-top.emltzs[0]-top.emltze[0])


###############################################################################
# --- Now, using above classes, read in and parse a MAD lattice file.
def getlattice(file):
  """Reads a MAD style lattice from the file and sets up the WARP lattice."""
  ff = open(file,'r')
  data = ff.readlines()
  ff.close()
  # --- Massage the data removing carriage returns, comments and blank lines
  i = 0
  while i < len(data):
    # --- Remove carriage return at end of line
    data[i] = data[i][:-1]
    # --- Remove all white space
    data[i] = re.sub('\s','',data[i])
    # --- Remove comments
    data[i] = re.sub('!.*','',data[i])
    data[i] = re.sub('#.*','',data[i])
    # --- Delete empty lines
    if data[i]:
      i = i + 1
    else:
      del data[i]
  # --- Massage the data removing line continuations.
  i = 0
  while i < len(data):
    mextend = re.search('&',data[i])
    while mextend:
      data[i] = re.sub('&',data[i+1],data[i])
      del data[i+1]
      mextend = re.search('&',data[i])
    i = i + 1
  # --- Massage the data into python syntax
  for i in xrange(len(data)):
    hibeam = re.search(':',data[i])
    if hibeam:
      # --- Only reformat lines which have colons in them. They are
      # --- hibeam style input. Other lines are WARP style, which
      # --- are direct python commands.
      # --- Replace '=' with '(' in LINE command
      # --- Replace first ',' with '(' in all except LINE commands
      mline = re.search('LINE',data[i])
      if mline:
        data[i] = re.sub('=','(',data[i])
      else:
        data[i] = re.sub(',','(',data[i],1)
      # --- Replace ':' with '=' in all commands
      data[i] = re.sub(':','=',data[i])
      # --- Append ')' to the end of the line
      data[i] = data[i] + ')'
  # --- Now the data is ready to be processed. Simply evaluate each line
  # --- and finish with the line in the 'USE' statement.
  for d in data:
    mfinish = re.search(r'(?P<u>USE\((?P<use>\w*)\))|(?P<end>END)',d)
    if mfinish:
      if mfinish.group('u'):
        lattice = eval(mfinish.group('use'),globals())
    else:
      exec(d,globals())

  # --- Convert MAD lattice to WARP lattice
  madtowarp(lattice)





