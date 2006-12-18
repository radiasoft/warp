"""Implements adaptive mesh refinement in 3d for the B field solver
"""
from warp import *
from MeshRefinement import *
from magnetostaticMG import MagnetostaticMG

try:
  import psyco
except ImportError:
  pass

#########################################################################
# Note that MRBlockB is psyco.bind at the end of the file
class MRBlockB(MeshRefinement,MagnetostaticMG):
  """
Implements adaptive mesh refinement in 3d for the B field solver
 - parent:
 - refinement=None: amount of refinement along each axis
 - lower,upper: extent of domain in relative to parent, in its own grid
                cell size, and not including any guard cells
 - dims: dimensions of the grid, only used for root block, the one with
         no parents
 - mins,maxs: locations of the grid lower and upper bounds in the beam frame
 - root: coarsest level grid
 - children: list of tuples, each containing three elements,
             (lower,upper,refinement). Children can also be added later
             using addchild.
 - lreducedpickle=true: when true, a small pickle is made by removing all of
                        the big arrays. The information can be regenerated
                        upon restart.
  """
  def __init__(self,parent=None,refinement=None,
                    lower=None,upper=None,
                    fulllower=None,fullupper=None,
                    dims=None,mins=None,maxs=None,
                    nguard=1,
                    children=None,**kw):

    # --- Note that this calls the MultiGrid __init__ as well.
    self.__class__.__bases__[0].__init__(self,
                    parent=parent,refinement=refinement,
                    lower=lower,upper=upper,
                    fulllower=fulllower,fullupper=fullupper,
                    dims=dims,mins=mins,maxs=maxs,
                    nguard=nguard,
                    children=children,**kw)


  def pcaxy(self,kwdict=None,**kw):
    if kwdict is None: kwdict = {}
    kwdict.update(kw)
    self.genericpf('getpotential',kwdict,2,pcaxy)
  def pcazx(self,kwdict=None,**kw):
    if kwdict is None: kwdict = {}
    kwdict.update(kw)
    self.genericpf('getpotential',kwdict,1,pcazx)
  def pcazr(self,kwdict=None,**kw):
    if kwdict is None: kwdict = {}
    kwdict.update(kw)
    kwdict['fullplane'] = 0
    self.genericpf('getpotential',kwdict,1,pcazx)
  def pcazy(self,kwdict=None,**kw):
    if kwdict is None: kwdict = {}
    kwdict.update(kw)
    self.genericpf('getpotential',kwdict,0,pcazy)
  def pcaztheta(self,kwdict=None,**kw):
    if kwdict is None: kwdict = {}
    kwdict.update(kw)
    kwdict['fullplane'] = 0
    self.genericpf('getpotential',kwdict,0,pcazy)
  def pcbxy(self,kwdict=None,**kw):
    if kwdict is None: kwdict = {}
    kwdict.update(kw)
    self.genericpf('getfield',kwdict,2,pcbxy)
  def pcbzx(self,kwdict=None,**kw):
    if kwdict is None: kwdict = {}
    kwdict.update(kw)
    self.genericpf('getfield',kwdict,1,pcbzx)
  def pcbzr(self,kwdict=None,**kw):
    if kwdict is None: kwdict = {}
    kwdict.update(kw)
    kwdict['fullplane'] = 0
    self.genericpf('getfield',kwdict,1,pcbzx)
  def pcbzy(self,kwdict=None,**kw):
    if kwdict is None: kwdict = {}
    kwdict.update(kw)
    self.genericpf('getfield',kwdict,0,pcbzy)
  def pcbztheta(self,kwdict=None,**kw):
    if kwdict is None: kwdict = {}
    kwdict.update(kw)
    kwdict['fullplane'] = 0
    self.genericpf('getfield',kwdict,0,pcbzy)
  def pcjxy(self,kwdict=None,**kw):
    if kwdict is None: kwdict = {}
    kwdict.update(kw)
    self.genericpf('getsource',kwdict,2,pcjxy)
  def pcjzx(self,kwdict=None,**kw):
    if kwdict is None: kwdict = {}
    kwdict.update(kw)
    self.genericpf('getsource',kwdict,1,pcjzx)
  def pcjzr(self,kwdict=None,**kw):
    if kwdict is None: kwdict = {}
    kwdict.update(kw)
    kwdict['fullplane'] = 0
    self.genericpf('getsource',kwdict,1,pcjzx)
  def pcjzy(self,kwdict=None,**kw):
    if kwdict is None: kwdict = {}
    kwdict.update(kw)
    self.genericpf('getsource',kwdict,0,pcjzy)
  def pcjztheta(self,kwdict=None,**kw):
    if kwdict is None: kwdict = {}
    kwdict.update(kw)
    kwdict['fullplane'] = 0
    self.genericpf('getsource',kwdict,0,pcjzy)


  def plaz(self,ix=None,iy=None,comp=2,color='fg',selfonly=0):
    if ix < self.fulllower[0]: return
    if iy < self.fulllower[1]: return
    if ix > self.fullupper[0]: return
    if iy > self.fullupper[1]: return
    if self is self.root: accumulateplotlists()
    try:
      plg(self.potential[comp,ix-self.fulllower[0]+1,iy-self.fulllower[1]+1,1:-1],self.zmesh,
          color=color)
      if not selfonly:
        for child in self.children:
          child.plaz(ix*child.refinement[0],iy*child.refinement[1],comp=comp,
                     color=color)
    finally:
      if self is self.root: plotlistofthings(lturnofflist=1)

  def plax(self,iy=None,iz=None,comp=2,color='fg',selfonly=0):
    if iy < self.fulllower[1]: return
    if iz < self.fulllower[2]: return
    if iy > self.fullupper[1]: return
    if iz > self.fullupper[2]: return
    if self is self.root: accumulateplotlists()
    try:
      plg(self.potential[comp,1:-1,iy-self.fulllower[1]+1,iz-self.fulllower[2]+1],self.xmesh,
          color=color)
      if not selfonly:
        for child in self.children:
          child.plax(iy*child.refinement[1],iz*child.refinement[2],comp=comp,
                     color=color)
    finally:
      if self is self.root: plotlistofthings(lturnofflist=1)

  def play(self,ix=None,iz=None,comp=2,color='fg',selfonly=0):
    if ix < self.fulllower[0]: return
    if iz < self.fulllower[2]: return
    if ix > self.fullupper[0]: return
    if iz > self.fullupper[2]: return
    if self is self.root: accumulateplotlists()
    try:
      plg(self.potential[comp,ix-self.fulllower[0]+1,1:-1,iz-self.fulllower[2]+1],self.ymesh,
          color=color)
      if not selfonly:
        for child in self.children:
          child.play(ix*child.refinement[0],iz*child.refinement[2],comp=comp,
                     color=color)
    finally:
      if self is self.root: plotlistofthings(lturnofflist=1)

  def pljz(self,ix=None,iy=None,comp=2,color='fg',selfonly=0):
    if ix < self.fulllower[0]: return
    if iy < self.fulllower[1]: return
    if ix > self.fullupper[0]: return
    if iy > self.fullupper[1]: return
    if self is self.root: accumulateplotlists()
    try:
      plg(self.source[comp,ix-self.fulllower[0],iy-self.fulllower[1],:],self.zmesh,
          color=color)
      if not selfonly:
        for child in self.children:
          child.pljz(ix*child.refinement[0],iy*child.refinement[1],comp=comp,
                     color=color)
    finally:
      if self is self.root: plotlistofthings(lturnofflist=1)

  def pljx(self,iy=None,iz=None,comp=2,color='fg',selfonly=0):
    if iy < self.fulllower[1]: return
    if iz < self.fulllower[2]: return
    if iy > self.fullupper[1]: return
    if iz > self.fullupper[2]: return
    if self is self.root: accumulateplotlists()
    try:
      plg(self.source[comp,:,iy-self.fulllower[1],iz-self.fulllower[2]],self.xmesh,
          color=color)
      if not selfonly:
        for child in self.children:
          child.pljx(iy*child.refinement[1],iz*child.refinement[2],comp=comp,
                     color=color)
    finally:
      if self is self.root: plotlistofthings(lturnofflist=1)

  def pljy(self,ix=None,iz=None,comp=2,color='fg',selfonly=0):
    if ix < self.fulllower[0]: return
    if iz < self.fulllower[2]: return
    if ix > self.fullupper[0]: return
    if iz > self.fullupper[2]: return
    if self is self.root: accumulateplotlists()
    try:
      plg(self.source[comp,ix-self.fulllower[0],:,iz-self.fulllower[2]],self.ymesh,
          color=color)
      if not selfonly:
        for child in self.children:
          child.pljy(ix*child.refinement[0],iz*child.refinement[2],comp=comp,
                     color=color)
    finally:
      if self is self.root: plotlistofthings(lturnofflist=1)

  def plbz(self,ix=None,iy=None,comp=2,color='fg',selfonly=0,withguard=1):
    if withguard:
      lower,upper = self.fulllower,self.fullupper
      iz = slice(None)
    else:
      lower,upper = self.lower,self.upper
      iz = slice(self.lower[2] - self.fulllower[2],
                 self.upper[2] - self.fulllower[2] + 1)
    if ix < lower[0]: return
    if iy < lower[1]: return
    if ix > upper[0]: return
    if iy > upper[1]: return
    if self is self.root: accumulateplotlists()
    try:
      plg(self.field[comp,ix-self.fulllower[0],iy-self.fulllower[1],iz],
          self.zmesh[iz],color=color)
      if not selfonly:
        for child in self.children:
          child.plbz(comp,ix*child.refinement[0],iy*child.refinement[1],
                         color=color,withguard=withguard)
    finally:
      if self is self.root: plotlistofthings(lturnofflist=1)

  def plbx(self,iy=None,iz=None,comp=2,color='fg',selfonly=0,withguard=1):
    if withguard:
      lower,upper = self.fulllower,self.fullupper
      ix = slice(None)
    else:
      lower,upper = self.lower,self.upper
      ix = slice(self.lower[0] - self.fulllower[0],
                 self.upper[0] - self.fulllower[0] + 1)
    if iy < lower[1]: return
    if iz < lower[2]: return
    if iy > upper[1]: return
    if iz > upper[2]: return
    if self is self.root: accumulateplotlists()
    try:
      plg(self.field[comp,ix,iy-self.fulllower[1],iz-self.fulllower[2]],
                     self.xmesh[ix],color=color)
      if not selfonly:
        for child in self.children:
          child.plbx(comp,iy*child.refinement[1],iz*child.refinement[2],
                         color=color,withguard=withguard)
    finally:
      if self is self.root: plotlistofthings(lturnofflist=1)

  def plby(self,ix=None,iz=None,comp=2,color='fg',selfonly=0,withguard=1):
    if withguard:
      lower,upper = self.fulllower,self.fullupper
      iy = slice(None)
    else:
      lower,upper = self.lower,self.upper
      iy = slice(self.lower[1] - self.fulllower[1],
                 self.upper[1] - self.fulllower[1] + 1)
    if ix < lower[0]: return
    if iz < lower[2]: return
    if ix > upper[0]: return
    if iz > upper[2]: return
    if self is self.root: accumulateplotlists()
    try:
      plg(self.field[comp,ix-self.fulllower[0],iy,iz-self.fulllower[2]],
          self.ymesh[iy],color=color)
      if not selfonly:
        for child in self.children:
          child.plby(comp,ix*child.refinement[0],iz*child.refinement[2],
                         color=color,withguard=withguard)
    finally:
      if self is self.root: plotlistofthings(lturnofflist=1)






