from warp import *
import __main__
lattice_version = "$Id: lattice.py,v 1.11 2002/12/02 19:16:25 dave Exp $"

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
    return ranf()
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
  def install(self,zz):
    self.expand()
    for e in self.elemslist: zz = e.install(zz)
    return zz
  def installdata(self):
    self.expand()
    for e in self.elemslist: e.installdata()
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
  def __init__(self,l=0,length=0,zshift=0,zs=0,ze=0,ap=0,
               offset_x=0,offset_y=0,error_type=''):
    self.l = l
    self.length = length
    self.zs = zs
    self.ze = ze
    self.type = ''
    self.zshift = zshift
    self.ap = ap
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
  def install(self):
    pass
  def installdata(self):
    pass
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
                  ap=aperture,
                  offset_x=offset_x,offset_y=offset_y,error_type=error_type)
    self.type = 'drift'
    self.i_cap_pointer = i_cap_pointer
    self.n_cap_nodes = n_cap_nodes
  def install(self,zz):
    top.ndrft = top.ndrft + 1
    self.idrft = top.ndrft
    if top.ndrft > len(top.drftzs)-1:
      idrft = top.ndrft
      top.ndrft = top.ndrft + 100
      gchange("Lattice")
      top.ndrft = idrft
    top.drftzs[top.ndrft] = zz + self.zshift
    top.drftze[top.ndrft] = top.drftzs[top.ndrft] + self.length
    top.drftap[top.ndrft] = self.ap
    top.drftox[top.ndrft] = self.offset_x*errordist(self.error_type)
    top.drftoy[top.ndrft] = self.offset_y*errordist(self.error_type)
    return top.drftze[top.ndrft]

class box(Elem):
  def __init__(self,l=0,length=0,zshift=0,zs=0,ze=0,aperture=0,
               width_x=0,width_y=0,error_type='',
               offset_x=0,offset_y=0,i_cap_pointer=0,n_cap_nodes=0):
    Elem.__init__(self,l=l,length=length,zshift=zshift,zs=zs,ze=ze,
                  ap=aperture,
                  offset_x=offset_x,offset_y=offset_y,error_type=error_type)
    self.type = 'box'
    self.width_x = width_x
    self.width_y = width_y
    self.i_cap_pointer = i_cap_pointer
    self.n_cap_nodes = n_cap_nodes
    self.derivedquantities(self)
  def derivedquantities(_self,self):
    if self.ap:
      self.width_x = self.ap
      self.width_y = self.ap
  def install(self,zz):
    top.ndrft = top.ndrft + 1
    self.idrft = top.ndrft
    if top.ndrft > len(top.drftzs)-1:
      idrft = top.ndrft
      top.ndrft = top.ndrft + 100
      gchange("Lattice")
      top.ndrft = idrft
    top.drftzs[top.ndrft] = zz + self.zshift
    top.drftze[top.ndrft] = top.drftzs[top.ndrft] + self.length
    top.drftap[top.ndrft] = -self.ap
    top.drftox[top.ndrft] = self.offset_x*errordist(self.error_type)
    top.drftoy[top.ndrft] = self.offset_y*errordist(self.error_type)
    return top.drftze[top.ndrft]

class quad(Elem):
  def __init__(self,l=0,length=0,zshift=0,zs=0,ze=0,aperture=0,
               voltage=0,gradient=0,r_elem=0,
               offset_x=0,offset_y=0,error_type='',
               i_cap_pointer=0,n_cap_nodes=0):
    Elem.__init__(self,l=l,length=length,zshift=zshift,zs=zs,ze=ze,
                  ap=aperture,
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
      self.voltage = self.gradient*self.ap**2
    else:
      self.gradient = self.voltage/self.ap**2
  def install(self,zz):
    top.nquad = top.nquad + 1
    self.iquad = top.nquad
    if top.nquad > len(top.quadzs)-1:
      iquad = top.nquad
      top.nquad = top.nquad + 100
      top.nqerr = top.nquad
      gchange("Lattice")
      top.nquad = iquad
    top.quadzs[top.nquad] = zz + self.zshift
    top.quadze[top.nquad] = top.quadzs[top.nquad] + self.length
    top.quadde[top.nquad] = self.gradient
    top.quadap[top.nquad] = self.ap
    top.quadrr[top.nquad] = self.r_elem
    top.qoffx[top.nquad] = self.offset_x*errordist(self.error_type)
    top.qoffy[top.nquad] = self.offset_y*errordist(self.error_type)
    return top.quadze[top.nquad]

class hyperb(Elem):
  def __init__(self,l=0,length=0,zshift=0,zs=0,ze=0,aperture=0,voltage=0,
               r_elem=0,offset_x=0,offset_y=0,error_type='',
               i_cap_pointer=0,n_cap_nodes=0):
    Elem.__init__(self,l=l,length=length,zshift=zshift,zs=zs,ze=ze,
                  ap=aperture,
                  offset_x=offset_x,offset_y=offset_y,error_type=error_type)
    self.type = 'hyperb'
    self.voltage = voltage
    self.r_elem = r_elem
    self.i_cap_pointer = i_cap_pointer
    self.n_cap_nodes = n_cap_nodes
  def install(self,zz):
    top.nquad = top.nquad + 1
    self.iquad = top.nquad
    if top.nquad > len(top.quadzs)-1:
      iquad = top.nquad
      top.nquad = top.nquad + 100
      top.nqerr = top.nquad
      gchange("Lattice")
      top.nquad = iquad
    top.quadzs[top.nquad] = zz + self.zshift
    top.quadze[top.nquad] = top.quadzs[top.nquad] + self.length
    top.qoffx[top.nquad] = self.offset_x*errordist(self.error_type)
    top.qoffy[top.nquad] = self.offset_y*errordist(self.error_type)
    top.quadde[top.nquad] = self.gradient
    top.quadap[top.nquad] = -self.ap
    return top.quadze[top.nquad]

class wire(Elem):
  def __init__(self,l=0,length=0,zshift=0,zs=0,ze=0,aperture=0,
               offset_x=0,offset_y=0,error_type=''):
    Elem.__init__(self,l=l,length=length,zshift=zshift,zs=zs,ze=ze,
                  ap=aperture,
                  offset_x=offset_x,offset_y=offset_y,error_type=error_type)
    self.type = 'wire'
  def install(self,zz):
    top.ndrft = top.ndrft + 1
    self.idrft = top.ndrft
    if top.ndrft > len(top.drftzs)-1:
      idrft = top.ndrft
      top.ndrft = top.ndrft + 100
      gchange("Lattice")
      top.ndrft = idrft
    top.drftzs[top.ndrft] = zz + self.zshift
    top.drftze[top.ndrft] = top.drftzs[top.ndrft] + self.length
    top.drftap[top.ndrft] = self.ap
    top.drftox[top.ndrft] = self.offset_x*errordist(self.error_type)
    top.drftoy[top.ndrft] = self.offset_y*errordist(self.error_type)
    return top.drftze[top.ndrft]

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
    Elem.__init__(self,l=l,length=length,zshift=zshift,zs=zs,ze=ze,ap=ap,
                  offset_x=ox,offset_y=oy,error_type=error_type)
    self.type = 'Drft'
  def install(self,zz):
    top.ndrft = top.ndrft + 1
    self.idrft = top.ndrft
    if top.ndrft > len(top.drftzs)-1:
      top.ndrft = top.ndrft + 100
      gchange("Lattice")
      top.ndrft = self.idrft
    top.drftzs[top.ndrft] = zz + self.zshift
    top.drftze[top.ndrft] = top.drftzs[top.ndrft] + self.length
    top.drftap[top.ndrft] = self.ap
    top.drftox[top.ndrft] = self.offset_x*errordist(self.error_type)
    top.drftoy[top.ndrft] = self.offset_y*errordist(self.error_type)
    return top.drftze[top.ndrft]

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
    Elem.__init__(self,l=l,length=length,zshift=zshift,zs=zs,ze=ze,ap=ap,
                  offset_x=ox,offset_y=oy,error_type=error_type)
    self.type = 'Bend'
    self.rc = rc
  def install(self,zz):
    top.nbend = top.nbend + 1
    self.ibend = top.nbend
    if top.nbend > len(top.bendzs)-1:
      ibend = top.nbend
      top.nbend = top.nbend + 100
      gchange("Lattice")
      top.nbend = ibend
    top.bendzs[top.nbend] = zz + self.zshift
    top.bendze[top.nbend] = top.bendzs[top.nbend] + self.length
    top.bendap[top.nbend] = self.ap
    #top.bendox[top.nbend] = self.offset_x*errordist(self.error_type)
    #top.bendoy[top.nbend] = self.offset_y*errordist(self.error_type)
    top.bendrc[top.nbend] = self.rc
    return top.bendze[ibend]

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
    Elem.__init__(self,l=l,length=length,zshift=zshift,zs=zs,ze=ze,ap=ap,
                  offset_x=ox,offset_y=oy,error_type=error_type)
    self.type = 'Dipo'
    for xx in ['ex','ey','bx','by','ta','tb','x1','x2','v1','v2',
               'l1','l2','w1','w2']:
      self.__dict__[xx] = locals()[xx]
  def install(self,zz):
    top.ndipo = top.ndipo + 1
    self.idipo = top.ndipo
    if top.ndipo > len(top.dipozs)-1:
      idipo = top.ndipo
      top.ndipo = top.ndipo + 100
      gchange("Lattice")
      top.ndipo = idipo
    top.dipozs[top.ndipo] = zz + self.zshift
    top.dipoze[top.ndipo] = top.dipozs[top.ndipo] + self.length
    #top.dipoox[top.ndipo] = self.offset_x*errordist(self.error_type)
    #top.dipooy[top.ndipo] = self.offset_y*errordist(self.error_type)
    for xx in ['ap','ex','ey','bx','by','ta','tb','x1','x2','v1','v2',
               'l1','l2','w1','w2']:
      aa = top.getpyobject('dipo'+xx)
      aa[top.ndipo] = self.__dict__[xx]
    return top.dipoze[top.ndipo]

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
    Elem.__init__(self,l=l,length=length,zshift=zshift,zs=zs,ze=ze,ap=ap,
                  offset_x=ox,offset_y=oy,error_type=error_type)
    self.type = 'Quad'
    for xx in ['de','db','vx','vy','rr','rl','gl','gp','pw','pa','pr']:
      self.__dict__[xx] = locals()[xx]
  def install(self,zz):
    top.nquad = top.nquad + 1
    self.iquad = top.nquad
    if top.nquad > len(top.quadzs)-1:
      top.nquad = top.nquad + 100
      top.nqerr = top.nquad
      gchange("Lattice")
      top.nquad = self.iquad
    top.quadzs[top.nquad] = zz + self.zshift
    top.quadze[top.nquad] = top.quadzs[top.nquad] + self.length
    top.qoffx[top.nquad] = self.offset_x*errordist(self.error_type)
    top.qoffy[top.nquad] = self.offset_y*errordist(self.error_type)
    for xx in ['ap','de','db','vx','vy','rr','rl','gl','gp','pw','pa','pr']:
      aa = top.getpyobject('quad'+xx)
      aa[top.nquad] = self.__dict__[xx]
    return top.quadze[top.nquad]

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
    Elem.__init__(self,l=l,length=length,zshift=zshift,zs=zs,ze=ze,ap=ap,
                  offset_x=ox,offset_y=oy,error_type=error_type)
    self.type = 'Sext'
    self.de = de
    self.db = db
  def install(self,zz):
    top.nsext = top.nsext + 1
    self.isext = top.nsext
    if top.nsext > len(top.sextzs)-1:
      isext = top.nsext
      top.nsext = top.nsext + 100
      gchange("Lattice")
      top.nsext = isext
    top.sextzs[top.nsext] = zz + self.zshift
    top.sextze[top.nsext] = top.sextzs[top.nsext] + self.length
    top.sextap[top.nsext] = self.ap
    top.sextox[top.nsext] = self.offset_x*errordist(self.error_type)
    top.sextoy[top.nsext] = self.offset_y*errordist(self.error_type)
    top.sextde[top.nsext] = self.de
    top.sextdb[top.nsext] = self.db
    return top.sextze[top.nsext]

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
  - ep=[] list of magnitudes of the derivative of electric multipole components
  - mp=[] list of magnitudes of the derivative of magnetic multipole components
  - pe=[] list of phase angles of electric multipole components
  - pm=[] list of phase angles of magnetic multipole components
  - rr=0 rod radius of an electric quad
  - rl=0 rod length of an electric quad
  - gl=0 gap length between rod and end plate of an electric quad
  - gp=0 position of the rod to plate gap in the x plane of an electric quad
  - pw=0 plate with of an electric quad
  - pa=0 plate aperture of an electric quad
  """
  def __init__(self,l=0,length=0,zshift=0,zs=0,ze=0,ap=0,ox=0,oy=0,
               error_type='',
               nn=[],vv=[],ae=[],am=[],ep=[],mp=[],pe=[],pm=[],
               rr=0,rl=0,gl=0,gp=0,pw=0,pa=0):
    Elem.__init__(self,l=l,length=length,zshift=zshift,zs=zs,ze=ze,ap=ap,
                  offset_x=ox,offset_y=oy,error_type=error_type)
    self.type = 'Hele'
    self.nn = nn
    self.vv = vv
    self.ae = ae
    self.am = am
    self.ep = ep
    self.mp = mp
    self.pe = pe
    self.pm = pm
    self.rr = rr
    self.rl = rl
    self.gl = gl
    self.gp = gp
    self.pw = pw
    self.pa = pa
    self.derivedquantities(self)
  def derivedquantities(_self,self):
    self.nn = array(self.nn)
    self.vv = array(self.vv)
    self.ae = array(self.ae)
    self.am = array(self.am)
    self.ep = array(self.ep)
    self.mp = array(self.mp)
    self.pe = array(self.pe)
    self.pm = array(self.pm)
  def install(self,zz):
    top.nhele = top.nhele + 1
    self.ihele = top.nhele
    if top.nhele > len(top.helezs)-1:
      ihele = top.nhele
      top.nhele = top.nhele + 100
      gchange("Lattice")
      top.nhele = ihele
    top.helezs[top.nhele] = zz + self.zshift
    top.heleze[top.nhele] = top.helezs[top.nhele] + self.length
    top.heleap[top.nhele] = self.ap
    top.heleox[top.nhele] = self.offset_x*errordist(self.error_type)
    top.heleoy[top.nhele] = self.offset_y*errordist(self.error_type)
    top.helerr[top.nhele] = self.rr
    top.helerl[top.nhele] = self.rl
    top.helegl[top.nhele] = self.gl
    top.helegp[top.nhele] = self.gp
    top.helepw[top.nhele] = self.pw
    top.helepa[top.nhele] = self.pa
    top.nhmlt = max(top.nhmlt,len(self.nn))
    return top.heleze[top.nhele]
  def installdata(self):
    top.hele_n[:,self.ihele] = self.nn
    top.hele_v[:,self.ihele] = self.vv
    top.heleae[:len(self.ae),self.ihele] = self.ae
    top.heleam[:len(self.am),self.ihele] = self.am
    top.heleep[:len(self.ep),self.ihele] = self.ep
    top.helemp[:len(self.mp),self.ihele] = self.mp
    top.helepe[:len(self.pe),self.ihele] = self.pe
    top.helepm[:len(self.pm),self.ihele] = self.pm

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
    Elem.__init__(self,l=l,length=length,zshift=zshift,zs=zs,ze=ze,ap=ap,
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
  def install(self,zz):
    top.naccl = top.naccl + 1
    self.iaccl = top.naccl
    if top.naccl > len(top.acclzs)-1:
      iaccl = top.naccl
      top.naccl = top.naccl + 100
      gchange("Lattice")
      top.naccl = iaccl
    top.acclzs[top.naccl] = zz + self.zshift
    top.acclze[top.naccl] = top.acclzs[top.naccl] + self.length
    top.acclap[top.naccl] = self.ap
    top.acclox[top.naccl] = self.offset_x*errordist(self.error_type)
    top.accloy[top.naccl] = self.offset_y*errordist(self.error_type)
    top.acclez[top.naccl] = self.ez
    top.acclxw[top.naccl] = self.xw
    top.acclsw[top.naccl] = self.sw
    top.acclts[top.naccl] = self.ts
    top.accldt[top.naccl] = self.dt
    top.ntaccl = max(top.ntaccl,shape(self.et)[0]-1)
    return top.acclze[top.naccl]
  def installdata(self):
    top.acclet[:len(self.et),self.iaccl] = self.et


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
  - rr=0 rod radius of an electric quad
  - rl=0 rod length of an electric quad
  - gl=0 gap length between rod and end plate of an electric quad
  - gp=0 position of the rod to plate gap in the x plane of an electric quad
  - pw=0 plate with of an electric quad
  - pa=0 plate aperture of an electric quad
  """
  def __init__(self,l=0,length=0,zshift=0,zs=0,ze=0,ap=0,ox=0,oy=0,
               error_type='',
               id=0,dz=0,e=[],ep=[],eph=[],nn=[],vv=[],ph=0,sf=0,sc=1,
               rr=0,rl=0,gl=0,gp=0,pw=0,pa=0):
    Elem.__init__(self,l=l,length=length,zshift=zshift,zs=zs,ze=ze,ap=ap,
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
    self.rr = rr
    self.rl = rl
    self.gl = gl
    self.gp = gp
    self.pw = pw
    self.pa = pa
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
  def install(self,zz):
    top.nemlt = top.nemlt + 1
    self.iemlt = top.nemlt
    if top.nemlt > len(top.emltzs)-1:
      iemlt = top.nemlt
      top.nemlt = top.nemlt + 100
      top.neerr = top.nemlt
      gchange("Lattice")
      top.nemlt = iemlt
    top.emltzs[top.nemlt] = zz + self.zshift
    top.emltze[top.nemlt] = top.emltzs[top.nemlt] + self.length
    top.emltap[top.nemlt] = self.ap
    top.emltox[top.nemlt] = self.offset_x*errordist(self.error_type)
    top.emltoy[top.nemlt] = self.offset_y*errordist(self.error_type)
    top.emltph[top.nemlt] = self.ph
    top.emltsf[top.nemlt] = self.sf
    top.emltsc[top.nemlt] = self.sc
    top.emltrr[top.nemlt] = self.rr
    top.emltrl[top.nemlt] = self.rl
    top.emltgl[top.nemlt] = self.gl
    top.emltgp[top.nemlt] = self.gp
    top.emltpw[top.nemlt] = self.pw
    top.emltpa[top.nemlt] = self.pa
    if not self.id:
      top.emltid[top.nemlt] = top.nemltsets + 1
    else:
      top.emltid[top.nemlt] = self.id
    top.nemltsets = max(top.emltid[top.nemlt],top.nemltsets)
    if self.e:
      top.nzemltmax = max(shape(self.e)[0]-1,top.nzemltmax)
      if len(shape(self.e)) == 2:
        top.nesmult = max(shape(self.e)[1],top.nesmult)
    if self.ep:
      top.nzemltmax = max(shape(self.ep)[0]-1,top.nzemltmax)
      if len(shape(self.ep)) == 2:
        top.nesmult = max(shape(self.ep)[1],top.nesmult)
    return top.emltze[top.nemlt]
  def installdata(self):
    id = top.emltid[self.iemlt]
    # --- Only copy data if this is a new data set
    if (self.e or self.ep) and top.nzemlt[id-1] == 0:
      top.nzemlt[id-1] = max(shape(self.e)[0],shape(self.ep)[0]) - 1
      if not self.dz:
        top.dzemlt[id-1] = ((top.emltze[self.iemlt] - top.emltzs[self.iemlt])/
                            top.nzemlt[id-1])
      else:
        top.dzemlt[id-1] = self.dz
      top.emlt_n[:len(self.nn)] = self.nn
      top.emlt_v[:len(self.vv)] = self.vv
      top.esemlt[:shape(self.e)[0],:shape(self.e)[1],id-1] = self.e
      top.esemltp[:shape(self.ep)[0],:shape(self.ep)[1],id-1] = self.ep
      top.esemltph[:shape(self.eph)[0],:shape(self.eph)[1],id-1] = self.eph

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
    Elem.__init__(self,l=l,length=length,zshift=zshift,zs=zs,ze=ze,ap=ap,
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
  def install(self,zz):
    top.nmmlt = top.nmmlt + 1
    self.immlt = top.nmmlt
    if top.nmmlt > len(top.mmltzs)-1:
      immlt = top.nmmlt
      top.nmmlt = top.nmmlt + 100
      top.nmerr = top.nmmlt
      gchange("Lattice")
      top.nmmlt = immlt
    top.mmltzs[top.nmmlt] = zz + self.zshift
    top.mmltze[top.nmmlt] = top.mmltzs[top.nmmlt] + self.length
    top.mmltap[top.nmmlt] = self.ap
    top.mmltox[top.nmmlt] = self.offset_x*errordist(self.error_type)
    top.mmltoy[top.nmmlt] = self.offset_y*errordist(self.error_type)
    top.mmltph[top.nmmlt] = self.ph
    top.mmltsf[top.nmmlt] = self.sf
    top.mmltsc[top.nmmlt] = self.sc
    if not self.id:
      top.mmltid[top.nmmlt] = top.nmmltsets + 1
    else:
      top.mmltid[top.nmmlt] = self.id
    top.nmmltsets = max(top.mmltid[top.nmmlt],top.nmmltsets)
    if self.m:
      top.nzmmltmax = max(shape(self.m)[0]-1,top.nzmmltmax)
      if len(shape(self.m)) == 2:
        top.nmsmult = max(shape(self.m)[1],top.nmsmult)
    if self.mp:
      top.nzmmltmax = max(shape(self.mp)[0]-1,top.nzmmltmax)
      if len(shape(self.mp)) == 2:
        top.nmsmult = max(shape(self.mp)[1],top.nmsmult)
    return top.mmltze[top.nmmlt]
  def installdata(self):
    id = top.mmltid[self.immlt]
    # --- Only copy data if this is a new data set
    if (self.m or self.mp) and top.nzmmlt[id-1] == 0:
      top.nzmmlt[id-1] = max(shape(self.m)[0],shape(self.mp)[0]) - 1
      if not self.dz:
        top.dzmmlt[id-1] = ((top.mmltze[self.immlt] - top.mmltzs[self.immlt])/
                            top.nzmmlt[id-1])
      else:
        top.dzmmlt[id-1] = self.dz
      top.mmlt_n[:len(self.nn)] = self.nn
      top.mmlt_v[:len(self.vv)] = self.vv
      top.msmmlt[:shape(self.m)[0],:shape(self.m)[1],id-1] = self.m
      top.msmmltp[:shape(self.mp)[0],:shape(self.mp)[1],id-1] = self.mp
      top.msmmltph[:shape(self.mph)[0],:shape(self.mph)[1],id-1] = self.mph


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
    Elem.__init__(self,l=l,length=length,zshift=zshift,zs=zs,ze=ze,ap=ap,
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
  def install(self,zz):
    top.nbgrd = top.nbgrd + 1
    self.ibgrd = top.nbgrd
    if top.nbgrd > len(top.bgrdzs)-1:
      ibgrd = top.nbgrd
      top.nbgrd = top.nbgrd + 100
      gchange("Lattice")
      top.nbgrd = ibgrd
    top.bgrdzs[top.nbgrd] = zz + self.zshift
    top.bgrdze[top.nbgrd] = top.bgrdzs[top.nbgrd] + self.length
    top.bgrdap[top.nbgrd] = self.ap
    top.bgrdox[top.nbgrd] = self.offset_x*errordist(self.error_type)
    top.bgrdoy[top.nbgrd] = self.offset_y*errordist(self.error_type)
    top.bgrdsf[top.nbgrd] = self.sf
    top.bgrdsc[top.nbgrd] = self.sc
    if not self.id:
      top.bgrdid[top.nbgrd] = top.bgrdns + 1
    else:
      top.bgrdid[top.nbgrd] = self.id
    top.bgrdns = max(top.bgrdid[top.nbgrd],top.bgrdns)
    if self.bx or self.by or self.bz:
      sbx = shape(self.bx)
      sby = shape(self.by)
      sbz = shape(self.bz)
      top.bgrdnx=max(sbx[0]-1,sby[0]-1,sbz[0]-1,top.bgrdnx)
      top.bgrdny=max(sbx[1]-1,sby[1]-1,sbz[1]-1,top.bgrdny)
      top.bgrdnz=max(sbx[2]-1,sby[2]-1,sbz[2]-1,top.bgrdnz)
    return top.bgrdze[top.nbgrd]
  def installdata(self):
    id = top.bgrdid[self.ibgrd]
    if (self.bx or self.by or self.bz) and top.bgrddx[id-1] == 0.:
      top.bgrddx[id-1] = self.dx
      top.bgrddy[id-1] = self.dy
      top.bgrddz[id-1] = self.dz
      sbx = shape(self.bx)
      sby = shape(self.by)
      sbz = shape(self.bz)
      top.bgrdbx[:sbx[0],:sbx[1],:sbx[2],id-1] = self.bx
      top.bgrdby[:sby[0],:sby[1],:sby[2],id-1] = self.by
      top.bgrdbz[:sbz[0],:sbz[1],:sbz[2],id-1] = self.bz

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
  - xs=0 x start of grid
  - ys=0 y start of grid
  - dx=0 x increment size (m)
  - dy=0 y increment size (m)
  - dz=0 z increment size (m)
  - sf=0 relative scaling factor
  - sc=0 absolute scaling factor
  - rr=0 rod radius of an electric quad
  - rl=0 rod length of an electric quad
  - gl=0 gap length between rod and end plate of an electric quad
  - gp=0 position of the rod to plate gap in the x plane of an electric quad
  - pw=0 plate with of an electric quad
  - pa=0 plate aperture of an electric quad
  """
  def __init__(self,l=0,length=0,zshift=0,zs=0,ze=0,ap=0,ox=0,oy=0,
               error_type='',
               id=0,sf=0,sc=1,pp=[],xs=0,ys=0,dx=0,dy=0,dz=0,
               rr=0,rl=0,gl=0,gp=0,pw=0,pa=0):
    Elem.__init__(self,l=l,length=length,zshift=zshift,zs=zs,ze=ze,ap=ap,
                  offset_x=ox,offset_y=oy,error_type=error_type)
    self.type = 'Pgrd'
    self.id = id
    self.sf = sf
    self.sc = sc
    self.pp = pp
    self.xs = xs
    self.ys = ys
    self.dx = dx
    self.dy = dy
    self.dz = dz
    self.rr = rr
    self.rl = rl
    self.gl = gl
    self.gp = gp
    self.pw = pw
    self.pa = pa
    self.derivedquantities(self)
  def derivedquantities(_self,self):
    self.pp = array(self.pp)
  def install(self,zz):
    top.npgrd = top.npgrd + 1
    self.ipgrd = top.npgrd
    if top.npgrd > len(top.pgrdzs)-1:
      ipgrd = top.npgrd
      top.npgrd = top.npgrd + 100
      gchange("Lattice")
      top.npgrd = ipgrd
    top.pgrdzs[top.npgrd] = zz + self.zshift
    top.pgrdze[top.npgrd] = top.pgrdzs[top.npgrd] + self.length
    top.pgrdap[top.npgrd] = self.ap
    top.pgrdxs[top.npgrd] = self.xs
    top.pgrdys[top.npgrd] = self.ys
    top.pgrdox[top.npgrd] = self.offset_x*errordist(self.error_type)
    top.pgrdoy[top.npgrd] = self.offset_y*errordist(self.error_type)
    top.pgrdsf[top.npgrd] = self.sf
    top.pgrdsc[top.npgrd] = self.sc
    top.pgrdrr[top.npgrd] = self.rr
    top.pgrdrl[top.npgrd] = self.rl
    top.pgrdgl[top.npgrd] = self.gl
    top.pgrdgp[top.npgrd] = self.gp
    top.pgrdpw[top.npgrd] = self.pw
    top.pgrdpa[top.npgrd] = self.pa
    if not self.id:
      top.pgrdid[top.npgrd] = top.pgrdns + 1
    else:
      top.pgrdid[top.npgrd] = self.id
    top.pgrdns = max(top.pgrdid[top.npgrd],top.pgrdns)
    if self.pp:
      top.pgrdnx = max(shape(self.pp)[0]-1,top.pgrdnx)
      top.pgrdny = max(shape(self.pp)[1]-1,top.pgrdny)
      top.pgrdnz = max(shape(self.pp)[2]-1,top.pgrdnz)
    return top.pgrdze[top.npgrd]
  def installdata(self):
    id = top.pgrdid[self.ipgrd]
    if self.pp and top.pgrddx[id-1] == 0.:
      top.pgrddx[id-1] = self.dx
      top.pgrddy[id-1] = self.dy
      top.pgrddz[id-1] = self.dz
      sp = shape(self.pp)
      top.pgrd[:sp[0],:sp[1],1:sp[2]+1,id-1] = self.pp

############################################################################
############################################################################
# --- Convert and copy the input MAD lattice into a WARP lattice.
def madtowarp(line,settunelen=1):
  """
Takes the inputted line obect and sets up the WARP lattice. Any existing
information if the WARP lattice arrays is deleted.
  """
  # --- Clear counts of all lattice elements and data array sizes.
  top.nquad = -1
  top.ndrft = -1
  top.nbend = -1
  top.ndipo = -1
  top.nhele = -1
  top.nhmlt = 0
  top.naccl = -1
  top.ntaccl = 0
  top.nemlt = -1
  top.nemltsets = 0
  top.nesmult = 0
  top.nzemltmax = 0
  top.nmmlt = -1
  top.nmmltsets = 0
  top.nmsmult = 0
  top.nzmmltmax = 0
  top.nbgrd = -1
  top.bgrdnx = 0
  top.bgrdny = 0
  top.bgrdnz = 0
  top.bgrdns = 0
  top.npgrd = -1
  top.pgrdnx = 0
  top.pgrdny = 0
  top.pgrdnz = 0
  top.pgrdns = 0
  top.nsext = -1

  # --- Install the line into fortran.
  # --- For elements with seperate data sets, only the number
  # --- of data points is checked. Arrays are then later allocated and filled.
  zz = line.install(0.)
  
  # --- Make sure element counts are at least zero.
  if top.nquad == -1: top.nquad = 0
  if top.ndrft == -1: top.ndrft = 0
  if top.nbend == -1: top.nbend = 0
  if top.ndipo == -1: top.ndipo = 0
  if top.nhele == -1: top.nhele = 0
  if top.naccl == -1: top.naccl = 0
  if top.nemlt == -1: top.nemlt = 0
  if top.nmmlt == -1: top.nmmlt = 0
  if top.nbgrd == -1: top.nbgrd = 0
  if top.npgrd == -1: top.npgrd = 0
  if top.nsext == -1: top.nsext = 0

  # --- Allocate the data space needed.
  gchange("Lattice")
  gchange("Mult_data")
  gchange("BGRDdata")
  gchange("PGRDdata")

  # --- Now, go back through the latice and fill in the data that needed
  # --- space allocated.
  line.installdata()

  # --- Finish by setting some general parameters.
  if settunelen:
    if top.nquad >= 2:
      top.tunelen=0.5*(top.quadzs[2]+top.quadze[2]-top.quadzs[0]-top.quadze[0])
    elif top.nhele >= 2:
      top.tunelen=0.5*(top.helezs[2]+top.heleze[2]-top.helezs[0]-top.heleze[0])
    elif top.nemlt >= 2:
      top.tunelen=0.5*(top.emltzs[2]+top.emltze[2]-top.emltzs[0]-top.emltze[0])
    elif top.nmmlt >= 2:
      top.tunelen=0.5*(top.mmltzs[2]+top.mmltze[2]-top.mmltzs[0]-top.mmltze[0])
    elif top.nbgrd >= 2:
      top.tunelen=0.5*(top.bgrdzs[2]+top.bgrdze[2]-top.bgrdzs[0]-top.bgrdze[0])
    elif top.npgrd >= 2:
      top.tunelen=0.5*(top.pgrdzs[2]+top.pgrdze[2]-top.pgrdzs[0]-top.pgrdze[0])

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





###############################################################################
###############################################################################
#########################################################################
# --- Utility routines to add in new elements into the middle of a lattice.

# ----------------------------------------------------------------------------
# --- DRFT --- XXX
def addnewdrft(zs,ze,ap=0.,ox=0.,oy=0.):
  """
Adds a new drft element to the lattice. The element will be placed at the
appropriate location.
Required arguments:
  - zs, ze: specify the start and end of the element
The following are all optional and have the same meaning and default as the
drft arrays with the same suffices:
  - ap,ox,oy
  """
  # --- Make sure that at least some of the element is in the proper range,
  # --- z >= 0., and if zlatperi != 0, z <= zlatperi.
  assert (zs < ze),"element start must be less than element end"
  assert (ze > 0.),"element end must be greater than zero"
  assert (top.zlatperi == 0.) or (zs < top.zlatperi),"element start must be less than zlatperi"

  # --- Get a dict of the input arguments and their values.
  ldict = locals()

  # --- Setup the lattice arrays for the insertion of the new element. If
  # --- there are already drfts, then find the place where the new one is to
  # --- be inserted and shift the existing data to open up a space.
  # --- Note that this uses that same check as in resetlat, that zs != ze to
  # --- determine whether or not a drft is defined.
  ie = 0
  # --- Find which element the new one goes before.
  while (ie <= top.ndrft and top.drftzs[ie] <= zs and
         top.drftzs[ie] != top.drftze[ie]):
    ie = ie + 1

  # --- Increase the size of the arrays by one. Except for the case when
  # --- there are no elements yet defined, which is true when the not'ed
  # --- statement is true.
  if ie > top.ndrft or top.drftzs[-1] != top.drftze[-1]:
    top.ndrft = top.ndrft + 100
    gchange("Lattice")

  # --- Setup dictionary relating lattice array with input argument names.
  # --- This is done here so that the references to the lattice arrays
  # --- refer to the updated memory locations after the gchange.
  edict={'zs':top.drftzs,'ze':top.drftze,'ap':top.drftap,
         'ox':top.drftox,'oy':top.drftoy}

  # --- Shift the existing data in the arrays to open up a space for the
  # --- new element.
  if ie <= top.ndrft:
    for e in edict.values():
      e[ie+1:] = e[ie:-1] + 0

  # --- Insert the new element. Note that edict correlates the lattice array
  # --- with the input arguments and ldict correlate the arguements with
  # --- their values.
  for (xx,e) in map(None,edict.keys(),edict.values()):
    e[ie] = ldict[xx]

# ----------------------------------------------------------------------------
# --- BEND --- XXX
def addnewbend(zs,ze,rc,ap=0.,ox=0.,oy=0.):
  """
Adds a new bend element to the lattice. The element will be placed at the
appropriate location.
Required arguments:
  - zs, ze: specify the start and end of the element
  - rc: radius of curvature of the bend
The following are all optional and have the same meaning and default as the
bend arrays with the same suffices:
  - ap,ox,oy
  """
  # --- Make sure that at least some of the element is in the proper range,
  # --- z >= 0., and if zlatperi != 0, z <= zlatperi.
  assert (zs < ze),"element start must be less than element end"
  assert (ze > 0.),"element end must be greater than zero"
  assert (top.zlatperi == 0.) or (zs < top.zlatperi),"element start must be less than zlatperi"

  # --- Get a dict of the input arguments and their values.
  ldict = locals()

  # --- Setup the lattice arrays for the insertion of the new element. If
  # --- there are already bends, then find the place where the new one is to
  # --- be inserted and shift the existing data to open up a space.
  # --- Note that this uses that same check as in resetlat, that zs != ze to
  # --- determine whether or not a bend is defined.
  ie = 0
  # --- Find which element the new one goes before.
  while (ie <= top.nbend and top.bendzs[ie] <= zs and
         top.bendzs[ie] != top.bendze[ie]):
    ie = ie + 1

  # --- Increase the size of the arrays by one. Except for the case when
  # --- there are no elements yet defined, which is true when the not'ed
  # --- statement is true.
  if ie > top.nbend or top.bendzs[-1] != top.bendze[-1]:
    top.nbend = top.nbend + 100
    gchange("Lattice")

  # --- Setup dictionary relating lattice array with input argument names.
  # --- This is done here so that the references to the lattice arrays
  # --- refer to the updated memory locations after the gchange.
  edict={'zs':top.bendzs,'ze':top.bendze,'rc':top.bendrc,'ap':top.bendap,
         'ox':top.bendox,'oy':top.bendoy}

  # --- Shift the existing data in the arrays to open up a space for the
  # --- new element.
  if ie <= top.nbend:
    for e in edict.values():
      e[ie+1:] = e[ie:-1] + 0

  # --- Insert the new element. Note that edict correlates the lattice array
  # --- with the input arguments and ldict correlate the arguements with
  # --- their values.
  for (xx,e) in map(None,edict.keys(),edict.values()):
    e[ie] = ldict[xx]

# ----------------------------------------------------------------------------
# --- DIPO --- XXX
def addnewdipo(zs,ze,by=0.,bx=0.,ta=0.,tb=0.,ex=0.,ey=0.,ap=0.,x1=0.,x2=0.,
               v1=0.,v2=0.,l1=0.,l2=0.,w1=0.,w2=0.):
  """
Adds a new dipo element to the lattice. The element will be placed at the
appropriate location.
Required arguments:
  - zs, ze: specify the start and end of the element
The following are all optional and have the same meaning and default as the
dipo arrays with the same suffices:
  - by,bx,ta,tb,ex,ey,ap,x1,x2,v1,v2,l1,l2,w1,w2
  """
  # --- Make sure that at least some of the element is in the proper range,
  # --- z >= 0., and if zlatperi != 0, z <= zlatperi.
  assert (zs < ze),"element start must be less than element end"
  assert (ze > 0.),"element end must be greater than zero"
  assert (top.zlatperi == 0.) or (zs < top.zlatperi),"element start must be less than zlatperi"

  # --- Get a dict of the input arguments and their values.
  ldict = locals()

  # --- Setup the lattice arrays for the insertion of the new element. If
  # --- there are already dipos, then find the place where the new one is to
  # --- be inserted and shift the existing data to open up a space.
  # --- Note that this uses that same check as in resetlat, that zs != ze to
  # --- determine whether or not a dipo is defined.
  ie = 0
  # --- Find which element the new one goes before.
  while (ie <= top.ndipo and top.dipozs[ie] <= zs and
         top.dipozs[ie] != top.dipoze[ie]):
    ie = ie + 1

  # --- Increase the size of the arrays by one. Except for the case when
  # --- there are no elements yet defined, which is true when the not'ed
  # --- statement is true.
  if ie > top.ndipo or top.dipozs[-1] != top.dipoze[-1]:
    top.ndipo = top.ndipo + 100
    gchange("Lattice")

  # --- Setup dictionary relating lattice array with input argument names.
  # --- This is done here so that the references to the lattice arrays
  # --- refer to the updated memory locations after the gchange.
  edict={'zs':top.dipozs,'ze':top.dipoze,'by':top.dipoby,'bx':top.dipobx,
         'ta':top.dipota,'tb':top.dipotb,'ex':top.dipoex,'ey':top.dipoey,
         'ap':top.dipoap,
         'x1':top.dipox1,'x2':top.dipox2,'v1':top.dipov1,'v2':top.dipov2,
         'l1':top.dipol1,'l2':top.dipol2,'w1':top.dipow1,'w2':top.dipow2}

  # --- Shift the existing data in the arrays to open up a space for the
  # --- new element.
  if ie <= top.ndipo:
    for e in edict.values():
      e[ie+1:] = e[ie:-1] + 0

  # --- Insert the new element. Note that edict correlates the lattice array
  # --- with the input arguments and ldict correlate the arguements with
  # --- their values.
  for (xx,e) in map(None,edict.keys(),edict.values()):
    e[ie] = ldict[xx]

# ----------------------------------------------------------------------------
# --- QUAD --- XXX
def addnewquad(zs,ze,db=0.,de=0.,et=0.,bt=0.,ts=0.,dt=0.,vx=0.,vy=0.,ap=0.,
               rr=0.,rl=0.,gl=0.,gp=0.,pw=0.,pa=0.,pr=0.,sl=0.,ox=0.,oy=0.,
               do=0.,
               glx=0.,gly=0.,axp=0.,axm=0.,ayp=0.,aym=0.,rxp=0.,rxm=0.,
               ryp=0.,rym=0.,vxp=0.,vxm=0.,vyp=0.,vym=0.,oxp=0.,oxm=0.,
               oyp=0.,oym=0.,pwl=0.,pwr=0.,pal=0.,par=0.,prl=0.,prr=0.):
  """
Adds a new quad element to the lattice. The element will be placed at the
appropriate location.
Required arguments:
  - zs, ze: specify the start and end of the element
The following are all optional and have the same meaning and default as the
quad and qdel arrays with the same suffices:
  - db,de,et,bt,ts,dt,vx,vy,ap,rr,rl,gl,gp,pw,pa,pr,sl,ox,oy,do
  - glx,gly,axp,axm,ayp,aym,rxp,rxm,ryp,rym,vxp,vxm,vyp,vym,oxp,oxm,
    oyp,oym,pwl,pwr,pal,par,prl,prr
  """
  # --- Make sure that at least some of the element is in the proper range,
  # --- z >= 0., and if zlatperi != 0, z <= zlatperi.
  assert (zs < ze),"element start must be less than element end"
  assert (ze > 0.),"element end must be greater than zero"
  assert (top.zlatperi == 0.) or (zs < top.zlatperi),"element start must be less than zlatperi"

  # --- Get a dict of the input arguments and their values.
  ldict = locals()

  # --- Setup the lattice arrays for the insertion of the new element. If
  # --- there are already quads, then find the place where the new one is to
  # --- be inserted and shift the existing data to open up a space.
  # --- Note that this uses that same check as in resetlat, that zs != ze to
  # --- determine whether or not a quad is defined.
  ie = 0
  # --- Find which element the new one goes before.
  while (ie <= top.nquad and top.quadzs[ie] <= zs and
         top.quadzs[ie] != top.quadze[ie]):
    ie = ie + 1

  # --- Increase the size of the arrays if it is needed.
  if ie > top.nquad or top.quadzs[-1] != top.quadze[-1]:
    top.nquad = top.nquad + 100
    top.nqerr = top.nqerr + 100
    gchange("Lattice")

  # --- Setup dictionary relating lattice array with input argument names.
  # --- This is done here so that the references to the lattice arrays
  # --- refer to the updated memory locations after the gchange.
  edict={'zs':top.quadzs,'ze':top.quadze,'db':top.quaddb,'de':top.quadde,
       'et':top.quadet,'bt':top.quadbt,'ts':top.quadts,'dt':top.quaddt,
       'vx':top.quadvx,'vy':top.quadvy,'ap':top.quadap,'rr':top.quadrr,
       'rl':top.quadrl,'gl':top.quadgl,'gp':top.quadgp,'pw':top.quadpw,
       'pa':top.quadpa,'pr':top.quadpr,'sl':top.quadsl,'do':top.quaddo,
       'glx':top.qdelglx,'gly':top.qdelgly,'axp':top.qdelaxp,'axm':top.qdelaxm,
       'ayp':top.qdelayp,'aym':top.qdelaym,'rxp':top.qdelrxp,'rxm':top.qdelrxm,
       'ryp':top.qdelryp,'rym':top.qdelrym,'vxp':top.qdelvxp,'vxm':top.qdelvxm,
       'vyp':top.qdelvyp,'vym':top.qdelvym,'oxp':top.qdeloxp,'oxm':top.qdeloxm,
       'oyp':top.qdeloyp,'oym':top.qdeloym,'pwl':top.qdelpwl,'pwr':top.qdelpwr,
       'pal':top.qdelpal,'par':top.qdelpar,'prl':top.qdelprl,'prr':top.qdelprr,
       'ox':top.qoffx,'oy':top.qoffy}

  # --- Shift the existing data in the arrays to open up a space for the
  # --- new element.
  if ie <= top.nquad:
    for e in edict.values():
      if len(shape(e)) == 1:
        e[ie+1:] = e[ie:-1] + 0
      else:
        # --- There are two arrays which are 2-D, quadet, and quadbt.
        e[:,ie+1:] = e[:,ie:-1] + 0

  # --- Insert the new element. Note that edict correlates the lattice array
  # --- with the input arguments and ldict correlate the arguements with
  # --- their values.
  for (xx,e) in map(None,edict.keys(),edict.values()):
    if len(shape(e)) == 1:
      e[ie] = ldict[xx]
    else:
      # --- There are two arrays which are 2-D, quadet, and quadbt.
      e[:,ie] = ldict[xx]

# ----------------------------------------------------------------------------
# --- SEXT --- XXX
def addnewsext(zs,ze,db=0.,de=0.):
  """
Adds a new sext element to the lattice. The element will be placed at the
appropriate location.
Required arguments:
  - zs, ze: specify the start and end of the element
The following are all optional and have the same meaning and default as the
sext arrays with the same suffices:
  - db,de
  """
  # --- Make sure that at least some of the element is in the proper range,
  # --- z >= 0., and if zlatperi != 0, z <= zlatperi.
  assert (zs < ze),"element start must be less than element end"
  assert (ze > 0.),"element end must be greater than zero"
  assert (top.zlatperi == 0.) or (zs < top.zlatperi),"element start must be less than zlatperi"

  # --- Get a dict of the input arguments and their values.
  ldict = locals()

  # --- Setup the lattice arrays for the insertion of the new element. If
  # --- there are already sexts, then find the place where the new one is to
  # --- be inserted and shift the existing data to open up a space.
  # --- Note that this uses that same check as in resetlat, that zs != ze to
  # --- determine whether or not a sext is defined.
  ie = 0
  # --- Find which element the new one goes before.
  while (ie <= top.nsext and top.sextzs[ie] <= zs and
         top.sextzs[ie] != top.sextze[ie]):
    ie = ie + 1

  # --- Increase the size of the arrays by one. Except for the case when
  # --- there are no elements yet defined, which is true when the not'ed
  # --- statement is true.
  if ie > top.nsext or top.sextzs[-1] != top.sextze[-1]:
    top.nsext = top.nsext + 100
    gchange("Lattice")

  # --- Setup dictionary relating lattice array with input argument names.
  # --- This is done here so that the references to the lattice arrays
  # --- refer to the updated memory locations after the gchange.
  edict={'zs':top.sextzs,'ze':top.sextze,'db':top.sextdb,'de':top.sextde}

  # --- Shift the existing data in the arrays to open up a space for the
  # --- new element.
  if ie <= top.nsext:
    for e in edict.values():
      e[ie+1:] = e[ie:-1] + 0

  # --- Insert the new element. Note that edict correlates the lattice array
  # --- with the input arguments and ldict correlate the arguements with
  # --- their values.
  for (xx,e) in map(None,edict.keys(),edict.values()):
    e[ie] = ldict[xx]

# ----------------------------------------------------------------------------
# --- HELE --- XXX
def addnewhele(zs,ze,ap=0.,ox=0.,oy=0.,rr=0.,rl=0.,gl=0.,gp=0.,pw=0.,pa=0.,
               ae=0.,am=0.,ep=0.,mp=0.,nn=0.,vv=0.,pe=0.,pm=0.):
  """
Adds a new hele element to the lattice. The element will be placed at the
appropriate location.
Required arguments:
  - zs, ze: specify the start and end of the element
The following are all optional and have the same meaning and default as the
hele arrays with the same suffices:
  - ap,ox,oy,rr,rl,gl,gp,pw,pa,ae,am,ep,mp,nn,vv,pe,pm
  """
  # --- Make sure that at least some of the element is in the proper range,
  # --- z >= 0., and if zlatperi != 0, z <= zlatperi.
  assert (zs < ze),"element start must be less than element end"
  assert (ze > 0.),"element end must be greater than zero"
  assert (top.zlatperi == 0.) or (zs < top.zlatperi),"element start must be less than zlatperi"

  # --- Get a dict of the input arguments and their values.
  ldict = locals()

  # --- Get the number of multipole components in this element.
  # --- Assume 1 for the case where ae and am are scalars.
  # --- Check top.nhmlt to see if the arrays are big enough for that number.
  nh = 1
  try: nh = len(ae)
  except: pass
  try: nh = len(am)
  except: pass
  if nh > top.nhmlt:
    top.nhmlt = nh
    gchange("Lattice")

  # --- Setup the lattice arrays for the insertion of the new element. If
  # --- there are already heles, then find the place where the new one is to
  # --- be inserted and shift the existing data to open up a space.
  # --- Note that this uses that same check as in resetlat, that zs != ze to
  # --- determine whether or not a hele is defined.
  ie = 0
  # --- Find which element the new one goes before.
  while (ie <= top.nhele and top.helezs[ie] <= zs and
         top.helezs[ie] != top.heleze[ie]):
    ie = ie + 1

  # --- Increase the size of the arrays by one. Except for the case when
  # --- there are no elements yet defined, which is true when the not'ed
  # --- statement is true.
  if ie > top.nhele or top.helezs[-1] != top.heleze[-1]:
    top.nhele = top.nhele + 100
    gchange("Lattice")

  # --- Setup dictionary relating lattice array with input argument names.
  # --- This is done here so that the references to the lattice arrays
  # --- refer to the updated memory locations after the gchange.
  edict={'zs':top.helezs,'ze':top.heleze,'ap':top.heleap,
         'ox':top.heleox,'oy':top.heleoy,'rr':top.helerr,'rl':top.helerl,
         'gl':top.helegl,'gp':top.helegp,'pw':top.helepw,'pa':top.helepa,
         'ae':top.heleae,'am':top.heleam,'ep':top.heleep,'mp':top.helemp,
         'nn':top.hele_n,'vv':top.hele_v,'pe':top.helepe,'pm':top.helepm}

  # --- Shift the existing data in the arrays to open up a space for the
  # --- new element.
  if ie <= top.nhele:
    for e in edict.values():
      if len(shape(e)) == 1:
        e[ie+1:] = e[ie:-1] + 0
      else:
        # --- These are quantities input for each multipole component
        e[:,ie+1:] = e[:,ie:-1] + 0

  # --- Insert the new element. Note that edict correlates the lattice array
  # --- with the input arguments and ldict correlate the arguements with
  # --- their values.
  for (xx,e) in map(None,edict.keys(),edict.values()):
    if len(shape(e)) == 1:
      e[ie] = ldict[xx]
    else:
      # --- These are quantities input for each multipole component
      e[:,ie] = ldict[xx]

# ----------------------------------------------------------------------------
# --- EMLT --- XXX
def addnewemlt(zs,ze,ap=0.,ph=0.,sf=0.,sc=1.,id=None,
               ox=0.,oy=0.,rr=0.,rl=0.,gl=0.,gp=0.,pw=0.,pa=0.,
               es=None,esp=None,phz=None,phpz=None,nn=None,vv=None):
  """
Adds a new emlt element to the lattice. The element will be placed at the
appropriate location.
Required arguments:
  - zs, ze: specify the start and end of the element
One and only one of the following must be supplied (if both are supplied, id
takes precedence):
  - id: data set ID corresponding to already existing emlt multipole data
  - es: 1- or 2-D array containing the multipole data. First dimension is data
    along z, optional second dimension is number of multipole components.
If 'es' is supplied, the following may also be supplied.
  - nn, vv: the multipole indices. If these are not specified, then it is
    assumed that the data in the 'es' array is layed out with the ordering of
    the existing emlt_n and emlt_v. Must have the same len as the second
    dimension of es (can be scalars is es is 1-D).
  - esp: first derivative of es along z. Must have the same shape as es.
  - phz, phpz: phase angle along z and it's derivative. Must have the same
    shape as es.
The following are all optional and have the same meaning and default as the
emlt arrays with the same suffices:
  - ap,ph,sf,sc,ox,oy,rr,rl,gl,gp,pw,pa
  """
  # --- Make sure either an 'id' or a dataset, 'es', was passed in.
  assert (id is not None or es is not None), \
         "either an 'id' or a dataset, 'es', must be passed in"
  # --- If a dataset was passed in, make sure that it has the same
  # --- number of multipole components or less, or that both n and v are
  # --- passed in
  assert (es is None) or \
         ((len(shape(es)) == 1) or (shape(es)[1] <= top.nesmult) or \
          (nn is not None and vv is not None)),\
         "The shape of the dataset must be consistent with the data already created or both n and v must be specified"
  # --- Make sure that at least some of the element is in the proper range,
  # --- z >= 0., and if zlatperi != 0, z <= zlatperi.
  assert (zs < ze),"element start must be less than element end"
  assert (ze > 0.),"element end must be greater than zero"
  assert (top.zlatperi == 0.) or (zs < top.zlatperi),"element start must be less than zlatperi"

  # --- Get a dict of the input arguments and their values.
  ldict = locals()

  # --- Setup the lattice arrays for the insertion of the new element. If
  # --- there are already emlts, then find the place where the new one is to
  # --- be inserted and shift the existing data to open up a space.
  # --- Note that this uses that same check as in resetlat, that emltid > 0,
  # --- to determine whether or not an emlt is defined.
  ie = 0
  # --- Find which element the new one goes before.
  while (ie <= top.nemlt and top.emltzs[ie] <= zs and top.emltid[ie] > 0):
    ie = ie + 1

  # --- Increase the size of the arrays by one. Except for the case when
  # --- there are no elements yet defined, which is true when the not'ed
  # --- statement is true.
  if ie > top.nemlt or top.emltid[-1] != 0:
    top.nemlt = top.nemlt + 100
    top.neerr = top.neerr + 100
    gchange("Lattice")

  # --- Setup dictionary relating lattice array with input argument names.
  # --- This is done here so that the references to the lattice arrays
  # --- refer to the updated memory locations after the gchange.
  edict = {'zs':top.emltzs,'ze':top.emltze,'ap':top.emltap,'ph':top.emltph,
           'sf':top.emltsf,'sc':top.emltsc,
           'ox':top.emltox,'oy':top.emltoy,'rr':top.emltrr,'rl':top.emltrl,
           'gl':top.emltgl,'gp':top.emltgp,'pw':top.emltpw,'pa':top.emltpa}

  # --- Shift the existing data in the arrays to open up a space for the
  # --- new element. The element id must be handled seperately.
  if ie <= top.nemlt:
    top.emltid[ie+1:] = top.emltid[ie:-1] + 0
    for e in edict.values():
      e[ie+1:] = e[ie:-1] + 0

  # --- Insert the new element. Note that edict correlates the lattice array
  # --- with the input arguments and ldict correlate the arguements with
  # --- their values.
  for (xx,e) in map(None,edict.keys(),edict.values()):
    e[ie] = ldict[xx]

  # --- Now setup the multipole component dataset.
  if id is not None:
    # --- If an 'id' was passed in, then just use that.
    top.emltid[ie] = id
  elif es is not None:
    # --- Otherwise, create a new dataset.
    top.nemltsets = top.nemltsets + 1
    top.emltid[ie] = top.nemltsets
    # --- Make sure that es is a 2-D array (first dimension is data versus z,
    # --- second is number of multipole components)
    if len(shape(es)) == 1:
      es = transpose(array([es]))
      if esp is not None: esp = transpose(array([esp]))
      if phz is not None: phz = transpose(array([phz]))
      if phpz is not None: phpz = transpose(array([phpz]))
    # --- Make sure that the first dimension of the arrays is long enough
    if shape(es)[0] > top.nzemltmax+1: top.nzemltmax = shape(es)[0] - 1
    # --- Change the sizes of the arrays
    gchange("Mult_data")
    # --- Set basic parameters
    n0 = shape(es)[0] # --- Number of data points along z
    n1 = shape(es)[1] # --- Number of multipole components
    top.nzemlt[-1] = n0 - 1
    top.dzemlt[-1] = (top.emltze[ie] - top.emltzs[ie])/(n0 - 1.)
    if nn is None and vv is None:
      # --- Assume n and v are ordered correctly and just copy the data in
      top.esemlt[:n0,:n1,-1] = es
      if esp is not None: top.esemltp[:n0,:n1,-1] = esp
      if phz is not None: top.esemltph[:n0,:n1,-1] = phz
      if phpz is not None: top.esemltphp[:n0,:n1,-1] = phpz
    else:
      # --- Make sure that n and v are lists
      if type(nn) in [IntType,FloatType]: nn = list([nn])
      else:                               nn = list(nn)
      if type(vv) in [IntType,FloatType]: vv = list([vv])
      else:                               vv = list(vv)
      # --- Make es a list of arrays
      es = list(transpose(es))
      if esp is not None: esp = list(transpose(esp))
      if phz is not None: phz = list(transpose(phz))
      if phpz is not None: phpz = list(transpose(phpz))
      # --- Loop over existing multipole components
      for i in xrange(top.nesmult):
        # --- Loop over input multipole components checking if any are the same
        for j in xrange(len(nn)):
          if nn[j] == top.emlt_n[i] and vv[j] == top.emlt_v[i]:
            # --- If so, then copy the data to the appropriate place and
            # --- delete the data from the lists.
            top.esemlt[:n0,i,-1] = es[j]
            if esp is not None: top.esemltp[:n0,i,-1] = esp[j]
            if phz is not None: top.esemltph[:n0,i,-1] = phz[j]
            if phpz is not None: top.esemltphp[:n0,i,-1] = phpz[j]
            del nn[j],vv[j],es[j]
            if esp is not None: del esp[j]
            if phz is not None: del phz[j]
            if phpz is not None: del phpz[j]
            break
      # --- Now copy in any left over data, increasing the number of multipole
      # --- components.
      if len(nn) > 0:
        ln = len(nn)
        top.nesmult = top.nesmult + ln
        gchange("Mult_data")
        top.emlt_n[-ln:] = nn
        top.emlt_v[-ln:] = vv
        top.esemlt[:n0,-ln:,-1] = transpose(array(es))
        if esp is not None: top.esemltp[:n0,-ln:,-1] = transpose(array(esp))
        if phz is not None: top.esemltph[:n0,-ln:,-1] = transpose(array(phz))
        if phpz is not None: top.esemltphp[:n0,-ln:,-1] = transpose(array(phpz))

  # --- Return the id of the new dataset. This allows the user to refer to
  # --- this new dataset without having to knowne its actual number.
  return top.emltid[ie]

# ----------------------------------------------------------------------------
# --- MMLT --- XXX
def addnewmmlt(zs,ze,ap=0.,ph=0.,sf=0.,sc=1.,id=None,ox=0.,oy=0.,
               ms=None,msp=None,phz=None,phpz=None,nn=None,vv=None):
  """
Adds a new mmlt element to the lattice. The element will be placed at the
appropriate location.
Required arguments:
  - zs, ze: specify the start and end of the element
One and only one of the following must be supplied (if both are supplied, id
takes precedence):
  - id: data set ID corresponding to already existing mmlt multipole data
  - ms: 1- or 2-D array containing the multipole data. First dimension is data
    along z, optional second dimension is number of multipole components.
If 'ms' is supplied, the following may also be supplied.
  - nn, vv: the multipole indices. If these are not specified, then it is
    assumed that the data in the 'ms' array is layed out with the ordering of
    the existing mmlt_n and mmlt_v. Must have the same len as the second
    dimension of ms (can be scalars is ms is 1-D).
  - msp: first derivative of ms along z. Must have the same shape as ms.
  - phz, phpz: phase angle along z and it's derivative. Must have the same
    shape as ms.
The following are all optional and have the same meaning and default as the
mmlt arrays with the same suffices:
  - ap,ph,sf,sc,ox,oy
  """
  # --- Make sure either an 'id' or a dataset, 'ms', was passed in.
  assert (id is not None or ms is not None), \
         "either an 'id' or a dataset, 'ms', must be passed in"
  # --- If a dataset was passed in, make sure that it has the same
  # --- number of multipole components or less, or that both n and v are
  # --- passed in
  assert (ms is None) or \
         ((len(shape(ms)) == 1) or (shape(ms)[1] <= top.nmsmult) or \
          (nn is not None and vv is not None)),\
         "The shape of the dataset must be consistent with the data already created or both n and v must be specified"
  # --- Make sure that at least some of the element is in the proper range,
  # --- z >= 0., and if zlatperi != 0, z <= zlatperi.
  assert (zs < ze),"element start must be less than element end"
  assert (ze > 0.),"element end must be greater than zero"
  assert (top.zlatperi == 0.) or (zs < top.zlatperi),"element start must be less than zlatperi"

  # --- Get a dict of the input arguments and their values.
  ldict = locals()

  # --- Setup the lattice arrays for the insertion of the new element. If
  # --- there are already mmlts, then find the place where the new one is to
  # --- be inserted and shift the existing data to open up a space.
  # --- Note that this uses that same check as in resetlat, that mmltid > 0,
  # --- to determine whether or not an mmlt is defined.
  ie = 0
  # --- Find which element the new one goes before.
  while (ie <= top.nmmlt and top.mmltzs[ie] <= zs and top.mmltid[ie] > 0):
    ie = ie + 1

  # --- Increase the size of the arrays by one. Except for the case when
  # --- there are no elements yet defined, which is true when the not'ed
  # --- statement is true.
  if ie > top.nmmlt or top.mmltid[-1] != 0:
    top.nmmlt = top.nmmlt + 100
    top.nmerr = top.nmerr + 100
    gchange("Lattice")

  # --- Setup dictionary relating lattice array with input argument names.
  # --- This is done here so that the references to the lattice arrays
  # --- refer to the updated memory locations after the gchange.
  edist = {'zs':top.mmltzs,'ze':top.mmltze,'ap':top.mmltap,'ph':top.mmltph,
           'sf':top.mmltsf,'sc':top.mmltsc,'ox':top.mmltox,'oy':top.mmltoy}

  # --- Shift the existing data in the arrays to open up a space for the
  # --- new element.
  if ie <= top.nmmlt:
    top.mmltid[ie+1:] = top.mmltid[ie:-1] + 0
    for e in edict.values():
      e[ie+1:] = e[ie:-1] + 0

  # --- Insert the new element. Note that edict correlates the lattice array
  # --- with the input arguments and ldict correlate the arguements with
  # --- their values.
  for (xx,e) in map(None,edict.keys(),edict.values()):
    e[ie] = ldict[xx]

  # --- Now setup the multipole component dataset.
  if id is not None:
    # --- If an 'id' was passed in, then just use that.
    top.mmltid[ie] = id
  elif ms is not None:
    # --- Otherwise, create a new dataset.
    top.nmmltsets = top.nmmltsets + 1
    top.mmltid[ie] = top.nmmltsets
    # --- Make sure that ms is a 2-D array (first dimension is data versus z,
    # --- second is number of multipole components)
    if len(shape(ms)) == 1:
      ms = transpose(array([ms]))
      if msp is not None: msp = transpose(array([msp]))
      if phz is not None: phz = transpose(array([phz]))
      if phpz is not None: phpz = transpose(array([phpz]))
    # --- Make sure that the first dimension of the arrays is long enough
    if shape(ms)[0] > top.nzmmltmax+1: top.nzmmltmax = shape(ms)[0] - 1
    # --- Change the sizes of the arrays
    gchange("Mult_data")
    # --- Set basic parameters
    n0 = shape(ms)[0] # --- Number of data points along z
    n1 = shape(ms)[1] # --- Number of multipole components
    top.nzmmlt[-1] = n0 - 1
    top.dzmmlt[-1] = (top.mmltze[ie] - top.mmltzs[ie])/(n0 - 1.)
    if nn is None and vv is None:
      # --- Assume n and v are ordered correctly and just copy the data in
      top.msmmlt[:n0,:n1,-1] = ms
      if msp is not None: top.msmmltp[:n0,:n1,-1] = msp
      if phz is not None: top.msmmltph[:n0,:n1,-1] = phz
      if phpz is not None: top.msmmltphp[:n0,:n1,-1] = phpz
    else:
      # --- Make sure that n and v are lists
      if type(nn) in [IntType,FloatType]: nn = list([nn])
      else:                               nn = list(nn)
      if type(vv) in [IntType,FloatType]: vv = list([vv])
      else:                               vv = list(vv)
      # --- Make ms a list of arrays
      ms = list(transpose(ms))
      if msp is not None: msp = list(transpose(msp))
      if phz is not None: phz = list(transpose(phz))
      if phpz is not None: phpz = list(transpose(phpz))
      # --- Loop over existing multipole components
      for i in xrange(top.nmsmult):
        # --- Loop over input multipole components checking if any are the same
        for j in xrange(len(nn)):
          if nn[j] == top.mmlt_n[i] and vv[j] == top.mmlt_v[i]:
            # --- If so, then copy the data to the appropriate place and
            # --- delete the data from the lists.
            top.msmmlt[:n0,i,-1] = ms[j]
            if msp is not None: top.msmmltp[:n0,i,-1] = msp[j]
            if phz is not None: top.msmmltph[:n0,i,-1] = phz[j]
            if phpz is not None: top.msmmltphp[:n0,i,-1] = phpz[j]
            del nn[j],vv[j],ms[j]
            if msp is not None: del msp[j]
            if phz is not None: del phz[j]
            if phpz is not None: del phpz[j]
            break
      # --- Now copy in any left over data, increasing the number of multipole
      # --- components.
      if len(nn) > 0:
        ln = len(nn)
        top.nmsmult = top.nmsmult + ln
        gchange("Mult_data")
        top.mmlt_n[-ln:] = nn
        top.mmlt_v[-ln:] = vv
        top.msmmlt[:n0,-ln:,-1] = transpose(array(ms))
        if msp is not None: top.msmmltp[:n0,-ln:,-1] = transpose(array(msp))
        if phz is not None: top.msmmltph[:n0,-ln:,-1] = transpose(array(phz))
        if phpz is not None: top.msmmltphp[:n0,-ln:,-1] = transpose(array(phpz))

  # --- Return the id of the new dataset. This allows the user to refer to
  # --- this new dataset without having to knowne its actual number.
  return top.mmltid[ie]

# ----------------------------------------------------------------------------
# --- ACCL --- XXX
def addnewaccl(zs,ze,ez=0.,ap=0.,ox=0.,oy=0.,xw=0.,sw=0.,et=0.,ts=0.,dt=0.):
  """
Adds a new accl element to the lattice. The element will be placed at the
appropriate location.
Required arguments:
  - zs, ze: specify the start and end of the element
The following are all optional and have the same meaning and default as the
accl arrays with the same suffices:
  - ez,ap,ox,oy,xw,sw,et,ts,dt
  """
  # --- Make sure that at least some of the element is in the proper range,
  # --- z >= 0., and if zlatperi != 0, z <= zlatperi.
  assert (zs < ze),"element start must be less than element end"
  assert (ze > 0.),"element end must be greater than zero"
  assert (top.zlatperi == 0.) or (zs < top.zlatperi),"element start must be less than zlatperi"

  # --- Get a dict of the input arguments and their values.
  ldict = locals()

  # --- Setup the lattice arrays for the insertion of the new element. If
  # --- there are already accls, then find the place where the new one is to
  # --- be inserted and shift the existing data to open up a space.
  # --- Note that this uses that same check as in resetlat, that zs != ze to
  # --- determine whether or not a accl is defined.
  ie = 0
  # --- Find which element the new one goes before.
  while (ie <= top.naccl and top.acclzs[ie] <= zs and
         top.acclzs[ie] != top.acclze[ie]):
    ie = ie + 1

  # --- Increase the size of the arrays if it is needed.
  if ie > top.naccl or top.acclzs[-1] != top.acclze[-1]:
    top.naccl = top.naccl + 100
    gchange("Lattice")

  # --- Setup dictionary relating lattice array with input argument names.
  # --- This is done here so that the references to the lattice arrays
  # --- refer to the updated memory locations after the gchange.
  edict={'zs':top.acclzs,'ze':top.acclze,'ez':top.acclez,'ap':top.acclap,
         'ox':top.acclox,'oy':top.accloy,'xw':top.acclxw,'sw':top.acclsw,
         'et':top.acclet,'ts':top.acclts,'dt':top.accldt}

  # --- Shift the existing data in the arrays to open up a space for the
  # --- new element.
  if ie <= top.naccl:
    for e in edict.values():
      if len(shape(e)) == 1:
        e[ie+1:] = e[ie:-1] + 0
      else:
        # --- acclet is 2-D
        e[:,ie+1:] = e[:,ie:-1] + 0

  # --- Insert the new element. Note that edict correlates the lattice array
  # --- with the input arguments and ldict correlate the arguements with
  # --- their values.
  for (xx,e) in map(None,edict.keys(),edict.values()):
    if len(shape(e)) == 1:
      e[ie] = ldict[xx]
    else:
      # --- acclet is 2-D
      e[:,ie] = ldict[xx]

# ----------------------------------------------------------------------------
# --- BGRD --- XXX
def addnewbgrd(zs,ze,id=None,xs=0.,ys=0.,ap=0.,ox=0.,oy=0.,ph=0.,sp=0.,cp=0.,
               sf=0.,sc=1.,sy=0,dx=None,dy=None,bx=None,by=None,bz=None):
  """
Adds a new bgrd element to the lattice. The element will be placed at the
appropriate location.
Required arguments:
  - zs, ze: specify the start and end of the element
Optionally, id may be specified, using a previously defined dataset
takes precedence):
  - id: data set ID corresponding to already existing bgrd multipole data
Or, one or more 3-D field arrays may be specified
  - bx, by, bz
  - dx,dy: transverse grid cell size must also be specified
The following are all optional and have the same meaning and default as the
bgrd arrays with the same suffices:
  - xs,ys,ap,ox,oy,ph,sp,cp,sf,sc,sy
  """
  # --- Make sure either an 'id' or a dataset, 'es', was passed in.
  assert None not in [id,bx,by,bz], \
         "either an 'id' or a dataset, bx, by, or bz, must be passed in"
  # --- Make sure that at least some of the element is in the proper range,
  # --- z >= 0., and if zlatperi != 0, z <= zlatperi.
  assert (zs < ze),"element start must be less than element end"
  assert (ze > 0.),"element end must be greater than zero"
  assert (top.zlatperi == 0.) or (zs < top.zlatperi),"element start must be less than zlatperi"

  # --- Get a dict of the input arguments and their values.
  ldict = locals()

  # --- Setup the lattice arrays for the insertion of the new element. If
  # --- there are already bgrds, then find the place where the new one is to
  # --- be inserted and shift the existing data to open up a space.
  # --- Note that this uses that same check as in resetlat, that bgrdid > 0,
  # --- to determine whether or not an bgrd is defined.
  ie = 0
  # --- Find which element the new one goes before.
  while (ie <= top.nbgrd and top.bgrdzs[ie] <= zs and top.bgrdid[ie] > 0):
    ie = ie + 1

  # --- Increase the size of the arrays by one. Except for the case when
  # --- there are no elements yet defined, which is true when the not'ed
  # --- statement is true.
  if ie > top.nbgrd or top.bgrdid[-1] != 0:
    top.nbgrd = top.nbgrd + 100
    gchange("Lattice")

  # --- Setup dictionary relating lattice array with input argument names.
  # --- This is done here so that the references to the lattice arrays
  # --- refer to the updated memory locations after the gchange.
  edict = {'zs':top.bgrdzs,'ze':top.bgrdze,'xs':top.bgrdxs,'ys':top.bgrdys,
           'ap':top.bgrdap,'ox':top.bgrdox,'oy':top.bgrdoy,'ph':top.bgrdph,
           'sp':top.bgrdsp,'cp':top.bgrdcp,'sf':top.bgrdsf,'sc':top.bgrdsc,
           'sy':top.bgrdsy}

  # --- Shift the existing data in the arrays to open up a space for the
  # --- new element. The element id must be handled seperately.
  if ie <= top.nbgrd:
    top.bgrdid[ie+1:] = top.bgrdid[ie:-1] + 0
    for e in edict.values():
      e[ie+1:] = e[ie:-1] + 0

  # --- Insert the new element. Note that edict correlates the lattice array
  # --- with the input arguments and ldict correlate the arguements with
  # --- their values.
  for (xx,e) in map(None,edict.keys(),edict.values()):
    e[ie] = ldict[xx]

  # --- Now setup the 3-D field grid dataset
  if id is not None:
    # --- If an 'id' was passed in, then just use that.
    top.bgrdid[ie] = id
  else:
    # --- Otherwise, create a new dataset.
    top.bgrdns = top.bgrdns + 1
    # --- Get array size
    if bx is not None: nx,ny,nz = shape(bx)
    if by is not None: nx,ny,nz = shape(by)
    if bz is not None: nx,ny,nz = shape(bz)
    # --- Make sure that the arrays are big enough
    top.bgrdnx = max(nx-1,top.bgrdnx)
    top.bgrdny = max(ny-1,top.bgrdny)
    top.bgrdnz = max(nz-1,top.bgrdnz)
    gchange("Mult_data")
    # --- Copy the data in
    top.bgrddx[-1] = dx
    top.bgrddy[-1] = dy
    top.bgrddz[-1] = (ze - zs)/(nz - 1)
    if bx is not None: top.bgrdbx[:nx,:ny,:nz,-1] = bx
    if by is not None: top.bgrdby[:nx,:ny,:nz,-1] = by
    if bz is not None: top.bgrdbz[:nx,:ny,:nz,-1] = bz

  # --- Return the id of the new dataset. This allows the user to refer to
  # --- this new dataset without having to knowne its actual number.
  return top.bgrdid[ie]

# ----------------------------------------------------------------------------
# --- PGRD --- XXX
def addnewpgrd(zs,ze,id=None,xs=0.,ys=0.,ap=0.,ox=0.,oy=0.,ph=0.,sp=0.,cp=0.,
               sf=0.,sc=1.,rr=0.,rl=0.,gl=0.,gp=0.,pw=0.,pa=0.,
               dx=None,dy=None,phi=None):

  """
Adds a new pgrd element to the lattice. The element will be placed at the
appropriate location.
Required arguments:
  - zs, ze: specify the start and end of the element
Optionally, id may be specified, using a previously defined dataset
takes precedence):
  - id: data set ID corresponding to already existing pgrd multipole data
Or, 3-D phi array may be specified
  - phi
  - dx,dy: transverse grid cell size must also be specified
The following are all optional and have the same meaning and default as the
pgrd arrays with the same suffices:
  - xs,ys,ap,ox,oy,ph,sp,cp,sf,sc,rr,rl,gl,gp,pw,pa
  """
  # --- Make sure either an 'id' or a dataset, 'es', was passed in.
  assert None not in [id,phi], \
         "either an 'id' or the dataset phi must be passed in"
  # --- Make sure that at least some of the element is in the proper range,
  # --- z >= 0., and if zlatperi != 0, z <= zlatperi.
  assert (zs < ze),"element start must be less than element end"
  assert (ze > 0.),"element end must be greater than zero"
  assert (top.zlatperi == 0.) or (zs < top.zlatperi),"element start must be less than zlatperi"

  # --- Get a dict of the input arguments and their values.
  ldict = locals()

  # --- Setup the lattice arrays for the insertion of the new element. If
  # --- there are already pgrds, then find the place where the new one is to
  # --- be inserted and shift the existing data to open up a space.
  # --- Note that this uses that same check as in resetlat, that pgrdid > 0,
  # --- to determine whether or not an pgrd is defined.
  ie = 0
  # --- Find which element the new one goes before.
  while (ie <= top.npgrd and top.pgrdzs[ie] <= zs and top.pgrdid[ie] > 0):
    ie = ie + 1

  # --- Increase the size of the arrays by one. Except for the case when
  # --- there are no elements yet defined, which is true when the not'ed
  # --- statement is true.
  if ie > top.npgrd or top.pgrdid[-1] != 0:
    top.npgrd = top.npgrd + 100
    gchange("Lattice")

  # --- Setup dictionary relating lattice array with input argument names.
  # --- This is done here so that the references to the lattice arrays
  # --- refer to the updated memory locations after the gchange.
  edict = {'zs':top.pgrdzs,'ze':top.pgrdze,'xs':top.pgrdxs,'ys':top.pgrdys,
           'ap':top.pgrdap,'ox':top.pgrdox,'oy':top.pgrdoy,'ph':top.pgrdph,
           'sp':top.pgrdsp,'cp':top.pgrdcp,'id':top.pgrdid,'sf':top.pgrdsf,
           'sc':top.pgrdsc,'rr':top.pgrdrr,'rl':top.pgrdrl,'gl':top.pgrdgl,
           'gp':top.pgrdgp,'pw':top.pgrdpw,'pa':top.pgrdpa}

  # --- Shift the existing data in the arrays to open up a space for the
  # --- new element. The element id must be handled seperately.
  if ie <= top.npgrd:
    top.pgrdid[ie+1:] = top.pgrdid[ie:-1] + 0
    for e in edict.values():
      e[ie+1:] = e[ie:-1] + 0

  # --- Insert the new element. Note that edict correlates the lattice array
  # --- with the input arguments and ldict correlate the arguements with
  # --- their values.
  for (xx,e) in map(None,edict.keys(),edict.values()):
    e[ie] = ldict[xx]

  # --- Now setup the 3-D field grid dataset
  if id is not None:
    # --- If an 'id' was passed in, then just use that.
    top.pgrdid[ie] = id
  else:
    # --- Otherwise, create a new dataset.
    top.pgrdns = top.pgrdns + 1
    # --- Get array size
    nx,ny,nz = shape(phi)
    # --- Make sure that the arrays are big enough
    top.pgrdnx = max(nx-1,top.pgrdnx)
    top.pgrdny = max(ny-1,top.pgrdny)
    top.pgrdnz = max(nz-1,top.pgrdnz)
    gchange("Mult_data")
    # --- Copy the data in
    top.bgrddx[-1] = dx
    top.bgrddy[-1] = dy
    top.bgrddz[-1] = (ze - zs)/(nz - 1)
    top.pgrd[:nx,:ny,:nz,-1] = phi

  # --- Return the id of the new dataset. This allows the user to refer to
  # --- this new dataset without having to knowne its actual number.
  return top.pgrdid[ie]

# ----------------------------------------------------------------------------
# --- Convenient plotting functions
def plotemlt(ie,m=0,p=0,color='fg'):
  """
Plots the field of the emlt element
  - ie: the element to plot
  - m=0: the multipole number to plot
  - color='fg': color of plot
  """
  id = top.emltid[ie] - 1
  dz = top.dzemlt[id]
  nz = top.nzemlt[id]
  zz = top.emltzs[ie] + iota(0,nz)*dz
  if p == 0:
    plg(top.esemlt[:nz+1,m,id],zz,color=color)
  else:
    plg(top.esemltp[:nz+1,m,id],zz,color=color)

def plotmmlt(im,m=0,p=0,color='fg'):
  """
Plots the field of the emlt element
  - im: the element to plot
  - m=0: the multipole number to plot
  - color='fg': color of plot
  """
  id = top.mmltid[im] - 1
  dz = top.dzmmlt[id]
  nz = top.nzmmlt[id]
  zz = top.mmltzs[im] + iota(0,nz)*dz
  if p == 0:
    plg(top.msmmlt[:nz+1,m,id],zz,color=color)
  else:
    plg(top.msmmltp[:nz+1,m,id],zz,color=color)

def plotacclet(ia,color='fg'):
  """
Plots the time dependent field of the accl element
  - ia: the element to plot
  - color='fg': color of plot
  """
  dt = top.accldt[ia]
  nt = top.ntaccl
  tt = top.acclts[ia] + iota(0,nt)*dt
  plg(top.acclet[:,ia],tt,color=color)

