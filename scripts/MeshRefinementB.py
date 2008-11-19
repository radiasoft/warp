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

  def plaz(self,comp=2,ix=None,iy=None,colors=None,selfonly=0,scale=1):
    self.plpotentialz(comp=comp,ix=ix,iy=iy,colors=colors,selfonly=selfonly,
                      scale=scale)

  def plax(self,comp=2,iy=None,iz=None,colors=None,selfonly=0,scale=1):
    self.plpotentialx(comp=comp,iy=iy,iz=iz,colors=colors,selfonly=selfonly,
                      scale=scale)

  def play(self,comp=2,ix=None,iz=None,colors=None,selfonly=0,scale=1):
    self.plpotentialy(comp=comp,ix=ix,iz=iz,colors=colors,selfonly=selfonly,
                      scale=scale)

  def pljz(self,comp=2,ix=None,iy=None,colors=None,selfonly=0,scale=1,
           withboundary=0):
    self.plsourcez(comp=comp,ix=ix,iy=iy,colors=colors,selfonly=selfonly,
                   scale=scale,withboundary=withboundary)

  def pljx(self,comp=2,iy=None,iz=None,colors=None,selfonly=0,scale=1,
           withboundary=0):
    self.plsourcex(comp=comp,iy=iy,iz=iz,colors=colors,selfonly=selfonly,
                   scale=scale,withboundary=withboundary)

  def pljy(self,comp=2,ix=None,iz=None,colors=None,selfonly=0,scale=1,
           withboundary=0):
    self.plsourcey(comp=comp,ix=ix,iz=iz,colors=colors,selfonly=selfonly,
                   scale=scale,withboundary=withboundary)

  def plbz(self,comp=2,ix=None,iy=None,colors=None,selfonly=0,scale=1,
           withguard=1):
    self.plfieldz(comp=comp,ix=ix,iy=iy,colors=colors,selfonly=selfonly,
                  scale=scale,withguard=withguard)

  def plbx(self,comp=2,iy=None,iz=None,colors=None,selfonly=0,scale=1,
           withguard=1):
    self.plfieldx(comp=comp,iy=iy,iz=iz,colors=colors,selfonly=selfonly,
                  scale=scale,withguard=withguard)

  def plby(self,comp=2,ix=None,iz=None,colors=None,selfonly=0,scale=1,
           withguard=1):
    self.plfieldy(comp=comp,ix=ix,iz=iz,colors=colors,selfonly=selfonly,
                  scale=scale,withguard=withguard)



