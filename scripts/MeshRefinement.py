"""Implements adaptive mesh refinement in 3d
"""
from warp import *
from multigrid import MultiGrid
from pyOpenDX import Visualizable,DXCollection,viewboundingbox
import MA
import __main__
import time
try:
  import psyco
except ImportError:
  pass

#########################################################################
# Note that MRBlock is psyco.bind at the end of the file
class MRBlock(MultiGrid,Visualizable):
  """
 - parent:
 - refinement=None: amount of refinement along each axis
 - lower,upper: extent of domain in relative to parent, in its own grid
                cell size, and not including any guard cells
 - ichild: the index number of this child in the given parent
 - dims: dimensions of the grid, only used for root block, the one with
         no parents
 - mins,maxs: locations of the grid lower and upper bounds in the beam frame
 - root: coarsest level grid
 - children: list of tuples, each containing three elements,
             (lower,upper,refinement). Children can also be added later
             using addchild.
  """
  def __init__(self,parent=None,refinement=None,
                    lower=None,upper=None,
                    ichild=None,
                    dims=None,mins=None,maxs=None,
                    nguard=1,
                    root=None,
                    children=None,**kw):

    if parent is None:
      # --- No parents, so just create empty lists
      self.parents = []
      self.ichild = []
      self.root = self
      # --- It is assumed that the root block will be the first one created.
      # --- So clear out the global block list and count.
      self.totalnumberofblocks = 0
      self.listofblocks = []
    else:
      # --- Save the parent and the index number. These are saved in lists
      # --- since a block can have multiple parents.
      self.parents = [parent]
      self.ichild = [ichild]
      self.root = root

    # --- Get the current global block number and increment the counter.
    # --- Also, add self to the global list of blocks.
    self.blocknumber = self.root.totalnumberofblocks
    self.root.totalnumberofblocks += 1
    self.root.listofblocks.append(self)

    self.overlaps = []
    self.nguard = nguard
    self.conductorlist = []

    if parent is None:
      # --- For the root, the dimensions and extent of the grid should
      # --- be specified. If not, they will be taken from w3d.
      self.dims = dims
      self.mins = mins
      self.maxs = maxs
      self.totalrefinement = ones(3)
      self.forcesymmetries = 1

      self.refinement = None

      # --- Make sure the top.nparpgrp is a large number. If it becomes too
      # --- small, fetche becomes inefficient since it is called many times,
      # --- once per each group. The not insignificant function call overhead
      # --- of python begins to use up a sizable chunk of time.
      top.nparpgrp = 100000
      
    else:

      # --- Make sure that refinement is an array of length three. If a scalar
      # --- is input, it is broadcast to all three axis.
      if type(refinement) in [IntType,FloatType]:
        refinement = 3*[refinement]
      self.refinement = array(refinement)

      self.totalrefinement = parent.totalrefinement*self.refinement
      self.deltas = parent.deltas/self.refinement
      self.rootdims = self.root.dims*self.totalrefinement
      self.forcesymmetries = 0

      if lower is None and upper is None:
        # --- The grid mins and maxs are input.
        self.mins = array(mins)
        self.maxs = array(maxs)
        self.lower = nint((self.mins - self.root.mins)/self.deltas)
        self.upper = nint((self.maxs - self.root.mins)/self.deltas)

      else:
        # --- The grid lower and upper bounds are input. The bounds are
        # --- relative to the root grid, but scaled by the total refinement.
        self.lower = array(lower)
        self.upper = array(upper)
        self.mins = self.root.mins + self.lower*self.deltas
        self.maxs = self.root.mins + self.upper*self.deltas

      # --- Now, extend the domain by the given number of guard cells. Checks
      # --- are made so that the domain doesn't extend beyond the original grid.
      self.fulllower = maximum(0,self.lower - nguard*self.refinement)
      self.fullupper = minimum(self.rootdims,
                               self.upper + nguard*self.refinement)

      # --- Recalculate grid quantities, including the guard regions.
      self.dims = self.fullupper - self.fulllower

      # --- Make sure that the number of grid points in each dimension is even
      self.fulllower = where(self.dims%2==1,self.fulllower-1,self.fulllower)
      self.fulllower = maximum(0,self.fulllower)
      self.dims = self.fullupper - self.fulllower
      self.fullupper = where(self.dims%2==1,self.fullupper+1,self.fullupper)
      self.fullupper = minimum(self.rootdims,self.fullupper)
      self.dims = self.fullupper - self.fulllower

      self.mins = self.root.mins + self.fulllower*self.deltas
      self.maxs = self.root.mins + self.fullupper*self.deltas

      # --- First, just use same boundary conditions as root.
      self.bounds = self.root.bounds.copy()
      self.pbounds = self.root.pbounds.copy()

      # --- Check if the mesh doesn't reach the edge of the root grid.
      # --- If not, switch to Dirichlet boundary.
      self.bounds[::2] = where(self.fulllower > 0,0,self.bounds[::2])
      self.bounds[1::2] = where(self.fullupper < self.rootdims,
                                0,self.bounds[1::2])
      self.l2symtry = self.root.l2symtry
      self.l4symtry = self.root.l4symtry
      self.pbounds[::2] = where(self.fulllower > 0,0,self.pbounds[::2])
      self.pbounds[1::2] = where(self.fullupper < self.rootdims,
                                0,self.pbounds[1::2])

      # --- Set so the solve is not parallelized
      self.nslaves = 1

    # --- Set individual quantities based on the values in the arrays,
    # --- if they have been set.
    if self.dims is not None:
      self.nx = self.dims[0]
      self.ny = self.dims[1]
      self.nz = self.dims[2]
    if self.mins is not None:
      self.xmmin = self.mins[0]
      self.ymmin = self.mins[1]
      self.zmmin = self.mins[2]
    if self.maxs is not None:
      self.xmmax = self.maxs[0]
      self.ymmax = self.maxs[1]
      self.zmmax = self.maxs[2]

    # --- Do some further initialization.
    MultiGrid.__init__(self,**kw)

    if parent is None:
      # --- This is only needed by the root grid in cases when the grid
      # --- parameters are obtained from w3d instead of the argument list.
      self.dims = array([self.nx,self.ny,self.nz])
      self.deltas = array([self.dx,self.dy,self.dz])
      self.mins = array([self.xmmin,self.ymmin,self.zmmin])
      self.maxs = array([self.xmmax,self.ymmax,self.zmmax])
      self.lower = zeros(3)
      self.upper = self.dims
      self.fulllower = zeros(3)
      self.fullupper = self.dims
      self.rootdims = self.dims
      # --- This may or may not be a good idea...
      if self.nx == w3d.nx and self.ny == w3d.ny and self.nz == w3d.nz:
        w3d.rho = self.rho
        w3d.phi = self.phi
        w3d.nxp = self.nx
        w3d.nyp = self.ny
        w3d.nzp = self.nz
        w3d.rhop = self.rho
        w3d.phip = self.phi

    # --- childdomains is the node centered grid which keeps track of which
    # --- cells are owned by which children. If there are no children,
    # --- then it is not needed.
    self.childdomains = None

    # --- siblingdomains is the node centered grid which keeps track of which
    # --- cells are owned by this instance and those owned by others at the
    # --- same level. If there are no siblings, then it is not needed.
    self.siblingdomains = None

    # --- Now add any specified children
    self.children = []
    if children is not None:
      for l,u,r in children:
        self.addchild(l,u,refinement=r)

  def __getstate__(self):
    """
Check whether this instance is the registered solver so that upon unpickling
it knows whether to re-register itself.
    """
    dict = self.__dict__.copy()
    if self is getregisteredsolver():
      dict['iamtheregisteredsolver'] = 1
    else:
      dict['iamtheregisteredsolver'] = 0
    return dict

  def __setstate__(self,dict):
    self.__dict__.update(dict)
    if self.iamtheregisteredsolver:
      registersolver(self)
      if self.nx == w3d.nx and self.ny == w3d.ny and self.nz == w3d.nz:
        w3d.rho = self.rho
        w3d.phi = self.phi
        w3d.nxp = self.nx
        w3d.nyp = self.ny
        w3d.nzp = self.nz
        w3d.rhop = self.rho
        w3d.phip = self.phi

  def addchild(self,lower,upper,mins=None,maxs=None,refinement=[2,2,2]):
    """
Add a mesh refined block to this block.
  -lower,upper,mins,maxs,refinement: All have same meanings as for the
                                     constructor.
    """
    child = MRBlock(parent=self,lower=lower,upper=upper,mins=mins,maxs=maxs,
                    refinement=refinement,ichild=len(self.children)+1,
                    nguard=self.nguard,root=self.root)
    self.addblockaschild(child)

  def addblockaschild(self,block):
    """
Given a block instance, installs it as a child.
    """
    self.children.append(block)

  def resetroot(self):
    # --- No parents, so just create empty lists
    self.parents = []
    self.ichild = []
    self.root = self
    self.totalnumberofblocks = 1
    self.listofblocks = [self]
    self.siblingdomains = None
    self.childdomains = None
    self.children = []

  #--------------------------------------------------------------------------
  # --- The next several methods handle conductors
  #--------------------------------------------------------------------------

  def installconductor(self,conductor,dfill=top.largepos):
    if not self.isfirstcall(): return
    MultiGrid.installconductor(self,conductor,dfill=dfill)
    for child in self.children:
      child.installconductor(conductor,dfill=dfill)

  def clearconductors(self):
    if not self.isfirstcall(): return
    MultiGrid.clearconductors(self)
    for child in self.children:
      child.clearconductors()

  def hasconductors(self):
    return self.conductors.interior.n > 0

  def getconductors(self,alllevels=1,result=None):
    if result is None: result = []
    result.append(self.conductors)
    if alllevels:
      for child in self.children:
        child.getconductors(alllevels,result)
    return result

  #--------------------------------------------------------------------------
  # --- The next several methods handle initialization that is done after
  # --- all blocks have been added.
  #--------------------------------------------------------------------------

  def finalize(self):
    # --- This should only be called at the top level.
    assert self.root == self,"finalize must only be called on the root block"
    blocklists = self.generateblocklevellists()
    self.clearparentsandchildren()
    self.findallchildren(blocklists)
    self.initializechilddomains()
    self.findoverlappingsiblings()

  def generateblocklevellists(self,blocklists=None):
    if blocklists is None:
      # --- This will only happen at the top level.
      # --- Create a list of empty lists. Each empty list will get the blocks
      # --- at the appropriate levels appended to it. Note that 100 is
      # --- assumed to be a large enough number - there almost certainly
      # --- will never be 100 levels of refinement.
      blocklists = [[] for i in range(100)]
    # --- Add this instance to the top level of the list and pass the rest
    # --- of it to the children
    if self not in blocklists[0]:
      blocklists[0].append(self)
    for child in self.children:
      b = child.generateblocklevellists(blocklists[1:])
    return blocklists

  def clearparentsandchildren(self):
    self.parents = []
    for child in self.children:
      child.clearparentsandchildren()
    self.children = []

  def findallchildren(self,blocklists):
    for block in blocklists[1]:
      # --- Get extent of possible overlapping domain
      l = maximum(block.fulllower/block.refinement,self.fulllower)
      u = minimum(block.fullupper/block.refinement,self.fullupper)
      if alltrue(u >= l):
        self.addblockaschild(block)
        block.parents.append(self)
        block.ichild.append(len(self.children))
    # --- Only the first block in the list makes the call for the next level.
    # --- This guarantees that this method is called only once for each block.
    if blocklists[0][0] == self:
      for block in blocklists[1]:
        block.findallchildren(blocklists[1:])

  def initializechilddomains(self):
    """
Sets the regions that are covered by the children.
    """
    if not self.isfirstcall(): return
    # --- Loop over the children, first calling each, then setting
    # --- childdomain appropriately.
    for child,ichild in zip(self.children,range(1,1+len(self.children))):
      child.initializechilddomains()

      # --- Set full domain to negative of child number first.
      l = maximum(self.fulllower,child.fulllower/child.refinement)
      u = child.fullupper/child.refinement
      # --- If the child extends to the edge of the parent mesh, it claims the
      # --- grid points on the upper edges.
      u = u + where(u == self.fullupper,1,0)
      # --- The child claims all unclaimed areas.
      ii = self.getchilddomains(l,u)
      #ii[...] = where(ii==0,-ichild,ii)
      ii[...] = where(ii==self.blocknumber,-child.blocknumber,ii)

      # --- Set interior to positive child number.
      l = maximum(self.fulllower,child.lower/child.refinement)
      u = child.upper/child.refinement
      # --- If the child extends to the edge of the parent mesh, it claims the
      # --- grid points on the upper edges.
      u = u + where(u == self.fullupper,1,0)
      # --- The child claims its full interior area
      ii = self.getchilddomains(l,u)
      #ii[...] = +ichild
      ii[...] = +child.blocknumber

  def findoverlappingsiblings(self,parent=None):
    # --- Recursively call the children.
    for child in self.children:
      child.findoverlappingsiblings(self)

    # --- If parent is None, then there are no siblings, so return.
    if parent is None: return

    # --- Note that this routine will be called once from each parent and that
    # --- each parent can have a different ordering of children. To avoid
    # --- problems with areas being claimed multiple times, only the region
    # --- within the calling parent is operated on. When created, the
    # --- siblingdomains is filled with -1 to flag all areas as being unset.
    # --- Get extent of domain within the calling parent
    l = maximum(parent.fulllower*self.refinement,self.fulllower)
    u = minimum(parent.fullupper*self.refinement,self.fullupper)
    # --- Loop over siblings with the parent, checking for overlaps with them.
    for sibling in parent.children:
      if sibling == self: continue

      # --- Find domain which overlaps the sibling. Do nothing if the
      # --- domain is zero size.
      sl = maximum(l,sibling.fulllower)
      su = minimum(u,sibling.fullupper)
      if sometrue(sl > su): continue
      sd = self.getsiblingdomains(sl,su)
      od = sibling.getsiblingdomains(sl,su)

      # --- This instance claims any cells that are unclaimed both here and in
      # --- the sibling
      unclaimed = (od == -1) & (sd == -1)
      sd[...] = where(unclaimed,1,sd)

      # --- Any cells that are marked as claimed here or in the sibling but
      # --- are not set in the other are marked as claimed by others.
      sd[...] = where(sd == -1,0,sd)
      od[...] = where(od == -1,0,od)

      # --- It is possible that some areas will be claimed by multiple
      # --- siblings. Remove claim to those areas from the sibling.
      overclaimed = (od == 1) & (sd == 1)
      od[...] = where(overclaimed,0,od)

      # --- Keep a list of overlapping siblings
      if sibling not in self.overlaps:
        self.overlaps.append(sibling)

    # --- Now, claim any parts of the domain within this parent that have
    # --- not already been set. i.e. those parts not covered by other blocks.
    sd = self.getsiblingdomains(l,u)
    sd[...] = where(sd == -1,1,sd)

  def clearinactiveregions(self,nbcells,parent=None,level=1):
    """
For regions which should not be refined but are included because of the
coalescing of blocks, the childdomains is set so that the E-fields will
not be fetched from there (it is set negative).
 - nbcells: array containing the refinement level of the grid cells.
            The shape will be the same as the intersection of self and the
            calling parent.
 - parent: the calling parent
 - level: the total amount refinement
    """
    if parent is None:
      # --- If there is no parent, it is assumed that this is the root block
      l = self.fulllower
      u = self.fullupper
    else:
      # --- Find intersection of parent and self.
      l = maximum(parent.fulllower*self.refinement,self.fulllower)
      u = minimum(parent.fullupper*self.refinement,self.fullupper)

    for child in self.children:
      # --- Find intersection of parent, self, and child
      cl = maximum(child.fulllower/child.refinement,l)
      cu = minimum(child.fullupper/child.refinement,u)
      if sometrue(cl > cu): continue

      # --- Get childdomains in the intersection region, Wherever the
      # --- the refinement level is lower than the childs, force childdomains
      # --- to be negative.
      ii = self.getchilddomains(cl,cu,1)
      nbc = self.getlocalarray(nbcells,cl,cu,fulllower=l)
      # --- Note that if the upper edge of the child does not extend to the
      # --- upper edge of self, then the childdomains at cu will have the
      # --- blocknumber of self, and so should  not be reset. The check
      # --- for (ii == child.blocknumber) prevents this. The check also
      # --- prevents one child from clearing out the domain of another, though
      # --- that wouldn't break anything since that other child would be
      # --- clearing out the same region itself. Also, with the check,
      # --- the second argument of the where can just be -ii instead
      # --- of -abs(ii) since only positive values of ii will have the sign
      # --- changed.
      r = child.refinement
      ii[...] = where((nbc<max(level*r)) & (ii == child.blocknumber),-ii,ii)

      # --- Stretch out the array so it has the refined cell size of the child
      nbcstretched = zeros(1+(cu-cl)*r)
      for k in range(r[2]):
        if k == 0: ksl = slice(None)
        else:      ksl = slice(-1)
        for j in range(r[1]):
          if j == 0: jsl = slice(None)
          else:      jsl = slice(-1)
          for i in range(r[0]):
            if i == 0: isl = slice(None)
            else:      isl = slice(-1)
            nbcstretched[i::r[0],j::r[1],k::r[2]] = nbc[isl,jsl,ksl]

      child.clearinactiveregions(nbcstretched,self,level*r)

  #--------------------------------------------------------------------------
  # --- The next several methods handle the charge density calculation.
  #--------------------------------------------------------------------------

  def loadrho(self,lzero=true,depositallparticles=0,lrootonly=0):
    """
Loads the charge density from the particles. This should only be called for
the top level grid.
    """
    if lzero: self.zerorho(lrootonly)
    for i,n,q,w in zip(top.ins-1,top.nps,top.sq,top.sw):
      if n == 0: continue
      self.setrho(top.xp[i:i+n],top.yp[i:i+n],top.zp[i:i+n],top.uzp[i:i+n],q,w,
                  depositallparticles,lrootonly)
    if not lrootonly and lzero:
      self.accumulaterhofromsiblings()
      self.getrhofromsiblings()
      if not depositallparticles:
        #self.gatherrhofromchildren_reversed()
        #self.gatherrhofromchildren_python()
        self.gatherrhofromchildren_fortran()
    self.makerhoperiodic()
    self.getrhoforfieldsolve()

  def zerorho(self,lrootonly=0):
    if not self.isfirstcall(): return
    self.rho[...] = 0.
    if not lrootonly:
      for child in self.children:
        child.zerorho()

  def setrho(self,x,y,z,uz,q,w,depositallparticles,lrootonly):
    # --- Choose among the various versions of setrho
    # --- allsort slightly edges out sort. The others are much slower.
    self.setrho_allsort(x,y,z,uz,q,w,lrootonly)
    #self.setrho_sort(x,y,z,uz,q,w,depositallparticles,lrootonly)
    #self.setrho_gather(x,y,z,uz,q,w,depositallparticles,lrootonly)
    #self.setrho_select(x,y,z,uz,q,w,depositallparticles,lrootonly)

  def getichildfromgrid(self,x,y,z,zgrid):
    ichild = zeros(len(x))
    getichild(0,len(x),x,y,z,ichild,
              self.nx,self.ny,self.nz,self.childdomains,
              self.xmmin,self.xmmax,self.ymmin,self.ymmax,
              self.zmmin,self.zmmax,zgrid,
              self.l2symtry,self.l4symtry)
    return ichild

  def sortbyichild(self,ichild,x,y,z,uz):
    if 0:
      # --- Sort the list of which domain the particles are in.
      ichildsorter = argsort(ichild)
      ichildsorted = take(ichild,ichildsorter)

      # --- Given a sorted list of child numbers, get how many of each number
      # --- there are. This would probably be better done in compiled
      # --- code. This should be efficient timewise, but uses two extra
      # --- full-size arrays.  
      ic1 = ichildsorted[1:] - ichildsorted[:-1]
      nn = compress(ic1 > 0,iota(len(ichildsorted-1)))
      nperchild = zeros(len(self.root.listofblocks))
      ss = 0
      for n in nn:
        nperchild[nint(ichildsorted[n-1])] = n - ss
        ss += nperchild[nint(ichildsorted[n-1])]
      nperchild[nint(ichildsorted[-1])] = len(ichildsorted) - sum(nperchild)

      # --- Now use the same sorting on the particle quantities.
      xout = take(x,ichildsorter)
      yout = take(y,ichildsorter)
      zout = take(z,ichildsorter)
      uzout = take(uz,ichildsorter)

    else:
      xout,yout,zout,uzout = zeros((4,len(x)),'d')
      nperchild = zeros(self.root.totalnumberofblocks)
      sortparticlesbyindex(len(x),ichild,x,y,z,uz,self.root.totalnumberofblocks,
                           xout,yout,zout,uzout,nperchild)

    return xout,yout,zout,uzout,nperchild

  def setrho_sort(self,x,y,z,uz,q,w,depositallparticles,lrootonly=0):
    """
Given the list of particles, a charge and a weight, deposits the charge
density of the mesh structure.
This method sorts the particles and then passes each block into the children.

    """
    if len(x) == 0: return
    if len(self.children) > 0 and not lrootonly:

      # --- Find out whether the particles are in the local domain or one of
      # --- the children's.
      ichild = zeros(len(x))
      nperchild = zeros(self.root.totalnumberofblocks)
      getichildandcount(len(x),x,y,z,ichild,
                        len(nperchild),nperchild,
                        self.nx,self.ny,self.nz,self.childdomains,
                        self.xmmin,self.xmmax,self.ymmin,self.ymmax,
                        self.zmmin,self.zmmax,top.zgrid,
                        self.l2symtry,self.l4symtry)

      if max(nperchild) <= len(x):
        xout,yout,zout,uzout = zeros((4,len(x)),'d')
        sortparticlesbyindexwithcounts(len(x),ichild,x,y,z,uz,
                 self.root.totalnumberofblocks,xout,yout,zout,uzout,nperchild)
        x,y,z,uz = xout,yout,zout,uzout

      # --- For each child, pass to it the particles in it's domain.
      ib = nonzero(nperchild)
      cc = take(self.root.listofblocks,ib)
      nn = take(nperchild,ib)
      i = 0
      for child,n in zip(cc,nn):
        i = i + n
        if child == self: continue
        child.setrho_sort(x[i-n:i],y[i-n:i],z[i-n:i],uz[i-n:i],q,w,
                          depositallparticles)

      # --- Now, get remaining particles in this domain
      if not depositallparticles:
        x = x[:nperchild[self.blocknumber]]
        y = y[:nperchild[self.blocknumber]]
        z = z[:nperchild[self.blocknumber]]
        uz = uz[:nperchild[self.blocknumber]]

    # --- Deposit the particles in this domain
    MultiGrid.setrho(self,x,y,z,uz,q,w)

  def setrho_select(self,x,y,z,uz,q,w,depositallparticles,lrootonly=0,
                    ichild=None):
    """
Given the list of particles, a charge and a weight, deposits the charge
density of the mesh structure.
This sends the entire list of particles to setrho3d, but selects out the
ones that deposit by setting the uz flag. Note that this relies on using
the setrho3dselect routine. For the cases tried, it actually isn't all that
slow, and is in fact faster than setrho_sort. That may change if there are a
large number of patches.
    """
    if len(x) == 0: return
    if len(self.children) > 0 and not lrootonly:
      # --- Find out whether the particles are in the local domain or one of
      # --- the children's.
      if ichild is None: ichild = zeros(len(x))
      getichild(self.blocknumber,len(x),x,y,z,ichild,
                self.nx,self.ny,self.nz,self.childdomains,
                self.xmmin,self.xmmax,self.ymmin,self.ymmax,
                self.zmmin,self.zmmax,
                top.zgrid,self.l2symtry,self.l4symtry)

      for block in [self]+self.children:
        # --- Get positions of those particles.
        if depositallparticles and block == self:
          uu = uz
        else:
          uu = ((uz != 0.) & (ichild == block.blocknumber)).astype('d')
        # --- Now deposit rho
        if block == self:
          MultiGrid.setrhoselect(self,x,y,z,uu,q,w)
        else:
          block.setrho_select(x,y,z,uu,q,w,depositallparticles,ichild=ichild)

    else:
      # --- Now deposit rho
      MultiGrid.setrhoselect(self,x,y,z,uz,q,w)

  def setrho_gather(self,x,y,z,uz,q,w,depositallparticles,lrootonly=0):
    """
Given the list of particles, a charge and a weight, deposits the charge
density of the mesh structure.
This loops over the children, gathering the particles that are passed to them.
This so far is the fastest routine, with depositallparticles=1. This is about
twice as fast as setrho_sort, and 50% faster than setrho_select. With
depositallparticles set, it is about 30% faster than with it not.
    """
    if len(x) == 0: return
    if len(self.children) > 0 and not lrootonly:
      # --- Find out whether the particles are in the local domain or one of
      # --- the children's.
      ichild = self.getichildfromgrid(x,y,z,top.zgrid)

      for block in [self]+self.children:
        if depositallparticles and block == self:
          xc = x
          yc = y
          zc = z
          uzc = uz
        else:
          # --- Get particles within the ith domain
          # --- Note that when block==self, the particles selected are the
          # --- ones from this instances domain.
          ii = nonzero(ichild==block.blocknumber)
          if len(ii) == 0: continue
          xc = take(x,ii)
          yc = take(y,ii)
          zc = take(z,ii)
          uzc = take(uz,ii)
        # --- Now deposit rho
        if block == self:
          MultiGrid.setrho(self,xc,yc,zc,uzc,q,w)
        else:
          block.setrho_gather(xc,yc,zc,uzc,q,w,depositallparticles)

    else:
      # --- Now deposit rho
      MultiGrid.setrho(self,x,y,z,uz,q,w)

  def setrho_allsort(self,x,y,z,uz,q,w,lrootonly):
    """
Given the list of particles, a charge and a weight, deposits the charge
density of the mesh structure.
This first gets the blocknumber of the block where each of the particles are
to be deposited. This is then sorted once. The loop is then over the list
of blocks, rather than walking through the tree structure.
This is about as fast as the setrho_select. The sort still takes up about 40%
of the time. Note that this depends on having ichilddomains filled with the
blocknumber rather than the child number relative to the parent.
    """
    if len(self.children) > 0 and not lrootonly:

      ichild = zeros(len(x))
      self.getichild_allsort(x,y,z,ichild)

      x,y,z,uz,nperchild = self.sortbyichild(ichild,x,y,z,uz)
      blocklist = self.root.listofblocks

    else:
      nperchild = [len(x)]
      blocklist = [self]

    # --- For each block, pass to it the particles in it's domain.
    i = 0
    for block,n in zip(blocklist,nperchild):
      MultiGrid.setrho(block,x[i:i+n],y[i:i+n],z[i:i+n],uz[i:i+n],q,w)
      i = i + n

  def getichild_allsort(self,x,y,z,ichild):
    """
Gathers the ichild for the setrho_allsort.
    """
    # --- This must wait until all of the parents have have set ichild
    # --- so that the value in the children takes precedence.
    if not self.islastcall(): return
    if len(x) == 0: return
    if len(self.children) > 0:
      # --- Find out whether the particles are in the local domain or one of
      # --- the children's.
      getichild(self.blocknumber,len(x),x,y,z,ichild,
                self.nx,self.ny,self.nz,self.childdomains,
                self.xmmin,self.xmmax,self.ymmin,self.ymmax,
                self.zmmin,self.zmmax,top.zgrid,
                self.l2symtry,self.l4symtry)
      for child in self.children:
        child.getichild_allsort(x,y,z,ichild)

  def gatherrhofromchildren_python(self):
    """
Python version.
    """
    # --- Do this only the first time this is called. This should only be
    # --- done once and since each parent requires that this be done
    # --- before it can get its rho from here, it must be done on the
    # --- first call.
    if not self.isfirstcall(): return

    # --- Loop over the children
    for child in self.children:

      # --- Make sure that the child has gathered rho from its children.
      child.gatherrhofromchildren_python()

      # --- Get coordinates of child relative to this domain
      l = maximum(child.fulllower/child.refinement,self.fulllower)
      u = minimum(child.fullupper/child.refinement,self.fullupper)

      # --- Check for any Nuemann boundaries
      dopbounds = (sometrue(child.pbounds == 1) and
                  sometrue(l == 0) and
                  sometrue(u == self.rootdims))

      # --- Loop over the three dimensions. The loops loop over all child grid
      # --- cells that contribute to a parent grid cell.
      r = child.refinement
      w = self.getwarrayforrho(r)
      for k in iota(-r[2]+1,r[2]-1):
        for j in iota(-r[1]+1,r[1]-1):
          for i in iota(-r[0]+1,r[0]-1):
            ls = array([i,j,k]) + l*r
            us = array([i,j,k]) + u*r
            pm = (ls < child.fulllower)
            pp = (us > child.fullupper)

            prho = self.getrho(l+pm,u-pp)
            crho = child.getrho(ls+r*pm,us-r*pp,r)
            cowns = child.getsiblingdomains(ls+r*pm,us-r*pp,r)

            # --- Get the rho that will be contributed.
            rh = where(cowns,crho,0.)
            # --- Multiply by the weight in place.
            multiply(rh,w[i+1,j+1,k+1],rh)
            # --- Adjust for symmetries
            if dopbounds:
              if child.pbounds[0] == 1 and i > 0 and l[0] == 0:
                rh[0,:,:] = 2.*rh[0,:,:]
              if child.pbounds[1] == 1 and i < 0 and u[0] == self.rootdims[0]:
                rh[-1,:,:] = 2.*rh[-1,:,:]
              if child.pbounds[2] == 1 and j > 0 and l[1] == 0:
                rh[:,0,:] = 2.*rh[:,0,:]
              if child.pbounds[3] == 1 and j < 0 and u[1] == self.rootdims[1]:
                rh[:,-1,:] = 2.*rh[:,-1,:]
              if child.pbounds[4] == 1 and k > 0 and l[2] == 0:
                rh[:,:,0] = 2.*rh[:,:,0]
              if child.pbounds[5] == 1 and k < 0 and u[2] == self.rootdims[2]:
                rh[:,:,-1] = 2.*rh[:,:,-1]

            # --- Now add in the contribution (in place)
            add(prho,rh,prho)

  def gatherrhofromchildren_fortran(self):
    """
Fortran version
    """
    # --- Do this only the first time this is called. This should only be
    # --- done once and since each parent requires that this be done
    # --- before it can get its rho from here, it must be done on the
    # --- first call.
    if not self.isfirstcall(): return

    # --- Loop over the children
    for child in self.children:

      # --- Make sure that the child has gathered rho from its children.
      child.gatherrhofromchildren_fortran()

      # --- Get coordinates of child relative to this domain
      l = maximum(child.fulllower/child.refinement,self.fulllower)
      u = minimum(child.fullupper/child.refinement,self.fullupper)

      # --- Check for any Nuemann boundaries
      dopbounds = (sometrue(child.pbounds == 1) and
                  (sometrue(l == 0) or
                   sometrue(u == self.rootdims)))

      w = self.getwarrayforrho(child.refinement)
      gatherrhofromchild(self.rho,self.dims[0],self.dims[1],self.dims[2],
                         child.rho,child.dims[0],child.dims[1],child.dims[2],
                         l,u,self.fulllower,child.fulllower,child.fullupper,
                         child.refinement,w,child.siblingdomains,
                         dopbounds,child.pbounds,self.rootdims)

  def gatherrhofromchildren_reversed(self):
    """
Python version with the loop over the child nodes as the inner loop
    """
    # --- Do this only the first time this is called. This should only be
    # --- done once and since each parent requires that this be done
    # --- before it can get its rho from here, it must be done on the
    # --- first call.
    if not self.isfirstcall(): return

    # --- Loop over the children
    for child in self.children:

      # --- Make sure that the child has gathered rho from its children.
      child.gatherrhofromchildren_reversed()

      r = child.refinement
      iwx,iwy,iwz = getmesh3d(-r[0]+1,1,2*r[0]-2,
                              -r[1]+1,1,2*r[1]-2,
                              -r[2]+1,1,2*r[2]-2)

      # --- Get coordinates of child relative to this domain
      l = maximum(child.fulllower/r,self.fulllower)
      u = minimum(child.fullupper/r,self.fullupper)

      # --- Check for any Nuemann boundaries
      dopbounds = (sometrue(child.pbounds == 1) and
                  sometrue(l == 0) and
                  sometrue(u == self.rootdims))

      w = self.getwarrayforrho(r)
      for iz in range(l[2],u[2]+1):
        iz0 = iz - self.fulllower[2]
        iz1 = maximum(iz*r[2] - r[2] + 1,child.fulllower[2])
        iz2 = minimum(iz*r[2] + r[2] - 1,child.fullupper[2])
        iwz1 = iz1 - (iz*r[2] - r[2] + 1)
        iwz2 = iz2 - (iz*r[2] - r[2] + 1) + 1
        for iy in range(l[1],u[1]+1):
          iy0 = iy - self.fulllower[1]
          iy1 = maximum(iy*r[1] - r[1] + 1,child.fulllower[1])
          iy2 = minimum(iy*r[1] + r[1] - 1,child.fullupper[1])
          iwy1 = iy1 - (iy*r[1] - r[1] + 1)
          iwy2 = iy2 - (iy*r[1] - r[1] + 1) + 1
          for ix in range(l[0],u[0]+1):
            ix0 = ix - self.fulllower[0]
            ix1 = maximum(ix*r[0] - r[0] + 1,child.fulllower[0])
            ix2 = minimum(ix*r[0] + r[0] - 1,child.fullupper[0])
            iwx1 = ix1 - (ix*r[0] - r[0] + 1)
            iwx2 = ix2 - (ix*r[0] - r[0] + 1) + 1

            crho = child.getrho([ix1,iy1,iz1],[ix2,iy2,iz2])
            cowns = child.getsiblingdomains([ix1,iy1,iz1],[ix2,iy2,iz2])

            # --- Get the rho that will be contributed.
            rh = where(cowns,crho,0.)

            # --- Multiply by the weight in place.
            wt = w[iwx1:iwx2,iwy1:iwy2,iwz1:iwz2]
            multiply(rh,wt,rh)

            # --- Adjust for symmetries
            if dopbounds:
              if child.pbounds[0] == 1 and ix == 0:
                rh[1:,:,:] = 2.*rh[1:,:,:]
              if child.pbounds[1] == 1 and ix == self.rootdims[0]:
                rh[:-1,:,:] = 2.*rh[:-1,:,:]
              if child.pbounds[2] == 1 and iy == 0:
                rh[:,1:,:] = 2.*rh[:,1:,:]
              if child.pbounds[3] == 1 and iy == self.rootdims[1]:
                rh[:,:-1,:] = 2.*rh[:,:-1,:]
              if child.pbounds[4] == 1 and iz == 0:
                rh[:,:,1:] = 2.*rh[:,:,1:]
              if child.pbounds[5] == 1 and iz == self.rootdims[2]:
                rh[:,:,:-1] = 2.*rh[:,:,:-1]

            # --- Now add in the contribution (in place)
            self.rho[ix0,iy0,iz0] += sum(ravel(rh))

  def accumulaterhofromsiblings(self):
    """
Loop over overlapping siblings and collect the rho from them from regions that
are owned by this instance.
    """
    # --- This should only be done once.
    if not self.isfirstcall(): return
    for child in self.children:
      child.accumulaterhofromsiblings()
    for other in self.overlaps:
      l = maximum(self.fulllower,other.fulllower)
      u = minimum(self.fullupper,other.fullupper)
      srho = self.getrho(l,u)
      orho = other.getrho(l,u)
      sown = self.getsiblingdomains(l,u)
      #srho[...] = srho + where(sown,orho,0.)
      add(srho,where(sown,orho,0.),srho)
      # --- Doesn't seem to work and I didn't take the time to debug it
      #nx,ny,nz = array(shape(srho)) - 1
      #addrhotoowner(nx,ny,nz,srho,sown,orho)

  def getrhofromsiblings(self):
    """
Get rho from overlapping siblings where they own the region.
    """
    if not self.isfirstcall(): return
    for child in self.children:
      child.getrhofromsiblings()
    for other in self.overlaps:
      l = maximum(self.fulllower,other.fulllower)
      u = minimum(self.fullupper,other.fullupper)
      srho = self.getrho(l,u)
      orho = other.getrho(l,u)
      oown = other.getsiblingdomains(l,u)
      srho[...] = where(oown,orho,srho)
      # --- Doesn't seem to work and I didn't take the time to debug it
      #nx,ny,nz = array(shape(srho)) - 1
      #getrhofromowner(nx,ny,nz,srho,oown,orho)

  #--------------------------------------------------------------------------
  # --- Methods to carry out the field solve
  #--------------------------------------------------------------------------

  def solve(self,iwhich=0):
    # --- Wait until all of the parents have called here until actually
    # --- doing the solve. This ensures that the phi in all of the parents
    # --- which is needed on the boundaries will be up to date.
    if not self.islastcall(): return

    self.setphiboundaries()
    MultiGrid.solve(self,iwhich)
    for child in self.children:
      child.solve(iwhich)

  def setphiboundaries(self):
    """
Sets phi on the boundaries, using the values from the parent grid
    """
    if len(self.parents) == 0: return
    for parent in self.parents:
      # --- Coordinates of mesh relative to parent's mesh location
      # --- and refinement. The minimum and maximum are needed in case
      # --- this mesh extends beyond the parent's.
      l = maximum(parent.fulllower*self.refinement,self.fulllower)
      u = minimum(parent.fullupper*self.refinement,self.fullupper)
      pl = l/self.refinement
      pu = u/self.refinement
      sphi = self.getphi(l,u)
      pphi = parent.getphi(pl,pu)
      r = self.refinement
      if l[0] == self.fulllower[0]: sphi[ 0    ,::r[1],::r[2]] = pphi[ 0,:,:]
      if u[0] == self.fullupper[0]: sphi[-1    ,::r[1],::r[2]] = pphi[-1,:,:]
      if l[1] == self.fulllower[1]: sphi[::r[0], 0    ,::r[2]] = pphi[:, 0,:]
      if u[1] == self.fullupper[1]: sphi[::r[0],-1    ,::r[2]] = pphi[:,-1,:]
      if l[2] == self.fulllower[2]: sphi[::r[0],::r[1], 0    ] = pphi[:,:, 0]
      if u[2] == self.fullupper[2]: sphi[::r[0],::r[1],-1    ] = pphi[:,:,-1]

    # --- Do remaining points where interpolation is needed.
    # --- This can be done, now that all of the contributions from the parents
    # --- is done.
    sphi = self.phi[:,:,1:-1]
    self.setslice(sphi[ 0,:,:],self.refinement[1:])
    self.setslice(sphi[-1,:,:],self.refinement[1:])
    self.setslice(sphi[:, 0,:],self.refinement[0::2])
    self.setslice(sphi[:,-1,:],self.refinement[0::2])
    self.setslice(sphi[:,:, 0],self.refinement[:2])
    self.setslice(sphi[:,:,-1],self.refinement[:2])

  def setslice(self,slice,r):
    # --- Does interpolation from cells coincident with the parent to the
    # --- other cells.
    for j in range(r[1]):
      w1 = 1.*j/r[1]
      for i in range(r[0]):
        w0 = 1.*i/r[0]
        if i == 0 and j == 0: continue
        slice[i:-1:r[0],j:-1:r[1]] = (
           slice[    :-1:r[0],    :-1:r[1]]*(1.-w0)*(1.-w1) +
           slice[r[0]:  :r[0],    :-1:r[1]]*    w0 *(1.-w1) +
           slice[    :-1:r[0],r[1]:  :r[1]]*(1.-w0)*    w1  +
           slice[r[0]:  :r[0],r[1]:  :r[1]]*    w0 *    w1)

    w1 = 1.
    for i in range(1,r[0]):
      w0 = 1.*i/r[0]
      slice[i::r[0],-1] = (
         slice[    :-1:r[0],-1]*(1.-w0)*    w1  +
         slice[r[0]:  :r[0],-1]*    w0 *    w1)

    w0 = 1.
    for j in range(1,r[1]):
      w1 = 1.*j/r[1]
      slice[-1,j::r[1]] = (
         slice[-1,    :-1:r[1]]*    w0 *(1.-w1) +
         slice[-1,r[1]:  :r[1]]*    w0 *    w1)

  def optimizeconvergence(self,resetpasses=1):
    if not self.isfirstcall(): return
    MultiGrid.optimizeconvergence(self,resetpasses=resetpasses)
    for child in self.children:
      child.optimizeconvergence(resetpasses=resetpasses)
    
  #--------------------------------------------------------------------------
  # --- Methods to fetch E-fields and potential
  #--------------------------------------------------------------------------

  def fetche(self):
    """
Fetches the E field. This should only be called at the root level grid.
    """
    self.fetchefrompositions_allsort(w3d.xfsapi,w3d.yfsapi,w3d.zfsapi,
                                     w3d.exfsapi,w3d.eyfsapi,w3d.ezfsapi)
    #self.fetchefrompositions_gather(w3d.xfsapi,w3d.yfsapi,w3d.zfsapi,
    #                                w3d.exfsapi,w3d.eyfsapi,w3d.ezfsapi)

  def fetchefrompositions_gather(self,x,y,z,ex,ey,ez):
    if len(x) == 0: return
    if len(self.children) > 0:

      # --- Find out whether the particles are in the local domain or one of
      # --- the children's.
      ichild = zeros(len(x))
      getichildpositiveonly(self.blocknumber,len(x),x,y,z,ichild,
                            self.nx,self.ny,self.nz,self.childdomains,
                            self.xmmin,self.xmmax,self.ymmin,self.ymmax,
                            self.zmmin,self.zmmax,
                            top.zgridprv,self.l2symtry,self.l4symtry)

      for block in [self]+self.children:
        # --- Get list of particles within the ith domain
        # --- Note that when block==self, the particles selected are the
        # --- ones from this instances domain.
        ii = nonzero(ichild==block.blocknumber)
        if len(ii) == 0: continue
        # --- Get positions of those particles.
        xc = take(x,ii)
        yc = take(y,ii)
        zc = take(z,ii)
        # --- Create temporary arrays to hold the E field
        tex,tey,tez = zeros((3,len(xc)),'d')
        # --- Now get the field
        if block == self:
          MultiGrid.fetchefrompositions(self,xc,yc,zc,tex,tey,tez)
        else:
          block.fetchefrompositions_gather(xc,yc,zc,tex,tey,tez)
        # --- Put the E field into the passed in arrays
        #put(ex,ii,tex)
        #put(ey,ii,tey)
        #put(ez,ii,tez)
        putsortedefield(len(tex),isort,tex,tey,tez,ex,ey,ez)

    else:

      # --- Get e-field from this domain
      MultiGrid.fetchefrompositions(self,x,y,z,ex,ey,ez)

  def getichild_positiveonly(self,x,y,z,ichild):
    """
Gathers the ichild for the fetche_allsort.
    """
    # --- This must wait until all of the parents have have set ichild
    # --- so that the value in the children takes precedence.
    if not self.islastcall(): return
    if len(x) == 0: return
    if len(self.children) > 0:

      # --- Find out whether the particles are in the local domain or one of
      # --- the children's.
      getichildpositiveonly(self.blocknumber,len(x),x,y,z,ichild,
                            self.nx,self.ny,self.nz,self.childdomains,
                            self.xmmin,self.xmmax,self.ymmin,self.ymmax,
                            self.zmmin,self.zmmax,top.zgridprv,
                            self.l2symtry,self.l4symtry)
      for child in self.children:
        child.getichild_positiveonly(x,y,z,ichild)

  def sortbyichildgetisort(self,ichild,x,y,z):
    if 0:
      # --- Sort the list of which domain the particles are in.
      ichildsorter = argsort(ichild)
      ichildsorted = take(ichild,ichildsorter)

      # --- Given a sorted list of child numbers, get how many of each number
      # --- there are. This would probably be better done in compiled
      # --- code. This should be efficient timewise, but uses two extra
      # --- full-size arrays.  
      ic1 = ichildsorted[1:] - ichildsorted[:-1]
      nn = compress(ic1 > 0,iota(len(ichildsorted-1)))
      #nperchild = zeros(1+len(self.children))
      nperchild = zeros(len(self.root.listofblocks))
      ss = 0
      for n in nn:
        nperchild[nint(ichildsorted[n-1])] = n - ss
        ss += nperchild[nint(ichildsorted[n-1])]
      nperchild[nint(ichildsorted[-1])] = len(ichildsorted) - sum(nperchild)

      # --- Now use the same sorting on the particle quantities.
      xout = take(x,ichildsorter)
      yout = take(y,ichildsorter)
      zout = take(z,ichildsorter)
      isort = ichildsorter

    else:
      xout,yout,zout = zeros((3,len(x)),'d')
      isort = zeros(len(x))
      nperchild = zeros(self.root.totalnumberofblocks)
      sortparticlesbyindexgetisort(len(x),ichild,x,y,z,
                                   self.root.totalnumberofblocks,
                                   xout,yout,zout,isort,nperchild)

    return xout,yout,zout,isort,nperchild

  def fetchefrompositions_allsort(self,x,y,z,ex,ey,ez):
    """
Given the list of particles, fetch the E fields.
This first gets the blocknumber of the block where each of the particles are
to be deposited. This is then sorted once. The loop is then over the list
of blocks, rather than walking through the tree structure.
This is about as fast as the setrho_select using the python sort. The sort
takes up about 40% of the time. It is significantly faster using the fortran
sort.
Note that this depends on having ichilddomains filled with the
blocknumber rather than the child number relative to the parent.
    """
    if len(self.children) > 0:

      ichild = zeros(len(x))
      # --- This assumes that the root block has blocknumber zero.
      self.getichild_positiveonly(x,y,z,ichild)

      x,y,z,isort,nperchild = self.sortbyichildgetisort(ichild,x,y,z)

      # --- Create temporary arrays to hold the E field
      tex,tey,tez = zeros((3,len(x)),'d')

    else:
      isort = None
      nperchild = [len(x)]
      tex,tey,tez = ex,ey,ez

    # --- For each block, pass to it the particles in it's domain.
    i = 0
    for block,n in zip(self.root.listofblocks,nperchild):
      MultiGrid.fetchefrompositions(block,x[i:i+n],y[i:i+n],z[i:i+n],
                                          tex[i:i+n],tey[i:i+n],tez[i:i+n])
      i = i + n

    # --- Now, put the E fields back into the original arrays, unsorting
    # --- the data
    if isort is not None:
      #put(ex,isort,tex)
      #put(ey,isort,tey)
      #put(ez,isort,tez)
      putsortedefield(len(tex),isort,tex,tey,tez,ex,ey,ez)

  def fetchphi(self):
    """
Fetches the potential. This should only be called at the root level grid.
    """
    self.fetchphifrompositions(w3d.xfsapi,w3d.yfsapi,w3d.zfsapi,w3d.phifsapi)

  def fetchphifrompositions(self,x,y,z,phi):
    """
Fetches the potential, given a list of positions
    """
    if len(x) == 0: return
    if len(self.children) > 0:

      # --- Find out whether the particles are in the local domain or one of
      # --- the children's. It is assumed at first to be from the local
      # --- domain, and is only set to one of the childs domains where
      # --- childdomains is positive (which does not include any guard cells).
      ichild = zeros(len(x))
      add(ichild,self.blocknumber,ichild)
      getichildpositiveonly(self.blocknumber,len(x),x,y,z,ichild,
                            self.nx,self.ny,self.nz,self.childdomains,
                            self.xmmin,self.xmmax,self.ymmin,self.ymmax,
                            self.zmmin,self.zmmax,
                            top.zgridprv,self.l2symtry,self.l4symtry)

      for block in [self]+self.children:
        # --- Get list of particles within the ith domain
        # --- Note that when block==self, the particles selected are the
        # --- ones from this instances domain.
        ii = nonzero(ichild==block.blocknumber)
        if len(ii) == 0: continue
        # --- Get positions of those particles.
        xc = take(x,ii)
        yc = take(y,ii)
        zc = take(z,ii)
        # --- Create temporary arrays to hold the potential
        tphi = zeros(len(xc),'d')
        # --- Now get the field
        if block == self:
          MultiGrid.fetchphifrompositions(self,xc,yc,zc,tphi)
        else:
          block.fetchphifrompositions(xc,yc,zc,tphi)
        # --- Put the potential into the passed in arrays
        put(phi,ii,tphi)

    else:

      # --- Get phi from this domain
      MultiGrid.fetchphifrompositions(self,x,y,z,phi)

  #--------------------------------------------------------------------------
  # --- Utility methods
  #--------------------------------------------------------------------------

  def islastcall(self):
    "Returns true when last parent has called"
    try:                   self.ncallsfromparents
    except AttributeError: self.ncallsfromparents = 0
    self.ncallsfromparents += 1
    if self.ncallsfromparents < len(self.parents): return 0
    self.ncallsfromparents = 0
    return 1

  def isfirstcall(self):
    "Returns true when first parent has called"
    try:                   self.ncallsfromparents
    except AttributeError: self.ncallsfromparents = 0
    self.ncallsfromparents += 1
    if self.ncallsfromparents > 1:
      if self.ncallsfromparents == len(self.parents):
        self.ncallsfromparents = 0
      return 0
    # --- This extra check is needed in case there is one or no parent.
    if self.ncallsfromparents >= len(self.parents):
      self.ncallsfromparents = 0
    return 1

  def setname(self,name='c',ichild=None):
    if not self.isfirstcall(): return
    import __main__
    if ichild is not None:
      name = name + '%d'%ichild
      __main__.__dict__[name] = self
    self.mainname = name
    for child,ichild in zip(self.children,range(1,1+len(self.children))):
      child.setname(name,ichild)

  def getmem(self):
    if not self.isfirstcall(): return
    memtot = product(self.dims + 1)
    for child in self.children:
      memtot = memtot + child.getmem()
    return memtot
      
  def getwarrayforrho(self,r):
    # --- Create weight array needed for rho deposition.
    # --- Is linear falloff in the weights correct for r > 2?
    wi = [0,0,0]
    for i in range(3):
      wi[i] = r[i] - abs(1.*iota(-r[i]+1,+r[i]-1))
      wi[i] = wi[i]/sum(wi[i])
    # --- Expand into 3 dimensions
    w = outerproduct(wi[0],outerproduct(wi[1],wi[2]))
    w.shape = (2*r[0]-1,2*r[1]-1,2*r[2]-1)
    return w

  def getphi(self,lower,upper):
    # --- Note that this takes into account the guard cells in z.
    ix1,iy1,iz1 = lower - self.fulllower
    ix2,iy2,iz2 = upper - self.fulllower + 1 
    iz1 = iz1 + 1
    iz2 = iz2 + 1
    return self.phi[ix1:ix2,iy1:iy2,iz1:iz2]
  def getrho(self,lower,upper,r=[1,1,1]):
    ix1,iy1,iz1 = lower - self.fulllower
    ix2,iy2,iz2 = upper - self.fulllower + 1
    return self.rho[ix1:ix2:r[0],iy1:iy2:r[1],iz1:iz2:r[2]]
  def getselfe(self,lower=None,upper=None,comp=slice(None),r=[1,1,1]):
    if lower is None: lower = self.lower
    if upper is None: upper = self.upper
    if type(comp) == StringType:
      comp = ['x','y','z'].index(comp)
    ix1,iy1,iz1 = lower - self.fulllower
    ix2,iy2,iz2 = upper - self.fulllower + 1
    selfe = MultiGrid.getselfe(self,recalculate=0)
    return selfe[comp,ix1:ix2:r[0],iy1:iy2:r[1],iz1:iz2:r[2]]
  def getchilddomains(self,lower,upper,upperedge=0):
    if self.childdomains is None:
      #self.childdomains = fzeros(1+self.dims)  + self.blocknumber
      self.childdomains = fzeros(1+self.dims)
      add(self.childdomains,self.blocknumber,self.childdomains)
    ix1,iy1,iz1 = lower - self.fulllower
    ix2,iy2,iz2 = upper - self.fulllower + upperedge
    return self.childdomains[ix1:ix2,iy1:iy2,iz1:iz2]
  def getsiblingdomains(self,lower,upper,r=[1,1,1]):
    if self.siblingdomains is None:
      #self.siblingdomains = zeros(1+self.dims) - 1
      self.siblingdomains = zeros(1+self.dims)
      subtract(self.siblingdomains,1,self.siblingdomains)
    ix1,iy1,iz1 = lower - self.fulllower
    ix2,iy2,iz2 = upper - self.fulllower + 1
    return self.siblingdomains[ix1:ix2:r[0],iy1:iy2:r[1],iz1:iz2:r[2]]
  def getlocalarray(self,array,lower,upper,r=[1,1,1],fulllower=None,
                    upperedge=1):
    if fulllower is None: fulllower = self.fulllower
    ix1,iy1,iz1 = lower - fulllower
    ix2,iy2,iz2 = upper - fulllower + upperedge
    return array[ix1:ix2:r[0],iy1:iy2:r[1],iz1:iz2:r[2]]

  def setmgtol(self,mgtol=None):
    """
Sets the convergence tolerance for all blocks. If mgtol is not given, it uses
f3d.mgtol.
    """
    if mgtol is None: mgtol = f3d.mgtol
    self.mgtol = mgtol
    for child in self.children:
      child.setmgtol(mgtol)
  def setmgmaxiters(self,mgmaxiters=None):
    """
Sets the maximum number of iterations for all blocks. If mgmaxiters is
not given, it uses f3d.mgmaxiters.
    """
    if mgmaxiters is None: mgmaxiters = f3d.mgmaxiters
    self.mgmaxiters = mgmaxiters
    for child in self.children:
      child.setmgmaxiters(mgmaxiters)

  def arraysliceoperation(self,ip,idim,arraystring,op,opnd,null,comp=None):
    """
Applies the operator to the array at the specified plane. The blocks only
contribute within their domains of ownership.
    """
    # --- Each block only needs to check once
    # --- XXX This call breaks something
    #if not self.islastcall(): return null
    # --- Don't do anything if the ip is outside the block
    if ip < self.fulllower[idim] or ip > self.fullupper[idim]: return null
    # --- Get the appropriate slice of phi and the childdomains array
    ii = [slice(None),slice(None),slice(None)]
    ii[idim] = ip - self.fulllower[idim]
    ix,iy,iz = ii
    if arraystring == 'phi': getarray = self.getphi
    elif arraystring == 'rho': getarray = self.getrho
    elif arraystring == 'selfe': getarray = self.getselfe
    if comp is None: array = getarray(self.fulllower,self.fullupper)
    else:            array = getarray(self.fulllower,self.fullupper,comp)
    c = self.getchilddomains(self.fulllower,self.fullupper,1)
    # --- Skip points that don't self doesn't own
    if c is not None:
      array = where(c[ix,iy,iz]==self.blocknumber,array[ix,iy,iz],null)
    # --- Find the max of self's and the children's phi
    result = opnd(array)
    for child in self.children:
      ipc = ip*child.refinement[idim]
      cresult = child.arraysliceoperation(ipc,idim,arraystring,op,opnd,null,
                                          comp)
      result = op(result,cresult)
    return result

  def getphislicemin(self,ip,idim):
    """
Finds the minimum value of phi at the specified plane. The blocks only
contribute within their domains of ownership.
    """
    return self.arraysliceoperation(ip,idim,'phi',min,minnd,+largepos)

  def getphislicemax(self,ip,idim):
    """
Finds the maximum value of phi at the specified plane. The blocks only
contribute within their domains of ownership.
    """
    return self.arraysliceoperation(ip,idim,'phi',max,maxnd,-largepos)

  def getrhoslicemin(self,ip,idim):
    """
Finds the minimum value of rho at the specified plane. The blocks only
contribute within their domains of ownership.
    """
    return self.arraysliceoperation(ip,idim,'rho',min,minnd,+largepos)

  def getrhoslicemax(self,ip,idim):
    """
Finds the maximum value of rho at the specified plane. The blocks only
contribute within their domains of ownership.
    """
    return self.arraysliceoperation(ip,idim,'rho',max,maxnd,-largepos)

  def getselfeslicemin(self,ip,idim,comp):
    """
Finds the minimum value of selfe at the specified plane. The blocks only
contribute within their domains of ownership.
    """
    return self.arraysliceoperation(ip,idim,'selfe',min,minnd,
                                    +largepos,comp)

  def getselfeslicemax(self,ip,idim,comp):
    """
Finds the maximum value of selfe at the specified plane. The blocks only
contribute within their domains of ownership.
    """
    return self.arraysliceoperation(ip,idim,'selfe',max,maxnd,
                                    -largepos,comp)

  #--------------------------------------------------------------------------
  # --- The following are used for plotting.
  #--------------------------------------------------------------------------

  def genericpf(self,kw,idim,pffunc,ip=None):
    """
Generic plotting routine. This plots only the local domain. Domains of the
children are also skipped, but the same call is made for them so they will
be plotted.
    """
    # --- Wait until all parents have called so that the child's domain
    # --- if not overlapped by a parent. This only affects the cellaray plots.
    # --- This also avoids the child plotting multiple times.
    if not self.islastcall(): return

    # --- Get the plane to be plotted
    if ip is None:
      ip = kw.get(('ix','iy','iz')[idim],None)
      if ip is None: ip = nint(-self.mins[idim]/self.deltas[idim])
    else:
      ip = ip*self.refinement[idim]
      kw[('ix','iy','iz')[idim]] = ip - self.fulllower[idim]

    # --- Set the values of cmin and cmax for all levels. This must be
    # --- done by the root level.
    if self is self.root:
      cmin = kw.get('cmin',None)
      cmax = kw.get('cmax',None)

      if kw.get('plotselfe',0):
        comp = kw.get('comp',2)
        if cmin is None: cmin = self.getselfeslicemin(ip,idim,comp)
        if cmax is None: cmax = self.getselfeslicemax(ip,idim,comp)
      elif kw.get('plotrho',0):
        if cmin is None: cmin = self.getrhoslicemin(ip,idim)
        if cmax is None: cmax = self.getrhoslicemax(ip,idim)
      else:
        if cmin is None: cmin = self.getphislicemin(ip,idim)
        if cmax is None: cmax = self.getphislicemax(ip,idim)

      kw['cmin'] = cmin
      kw['cmax'] = cmax

    # --- Only make the plot if the plane is included in the domain.
    # --- Even if there is no overlap, the children must be called since
    # --- they may overlap (in the domain of a different parent).
    # --- Note that the full extent is used so that the childdomains slice
    # --- will have the same shape as ireg.
    if self.fulllower[idim] <= ip and ip <= self.fullupper[idim]:
      if self.childdomains is not None:
        # --- Create the ireg array, which will be set to zero in the domain
        # --- of the children.
        ss = list(shape(self.rho))
        del ss[idim]
        ireg = zeros(ss)
        ii = [slice(-1),slice(-1),slice(-1)]
        ii[idim] = ip - self.fulllower[idim]
        ix,iy,iz = ii
        ireg[1:,1:]=equal(self.childdomains[ix,iy,iz],self.blocknumber)
        if idim != 2: ireg = transpose(ireg)
        kw['ireg'] = ireg
      else:
        kw['ireg'] = None
      MultiGrid.genericpf(self,kw,pffunc)
      kw['titles'] = 0
      kw['lcolorbar'] = 0

    for child in self.children:
      child.genericpf(kw,idim,pffunc,ip)

  def pfxy(self,kwdict=None,**kw):
    if kwdict is None: kwdict = {}
    kwdict.update(kw)
    self.genericpf(kwdict,2,pfxy)
  def pfzx(self,kwdict=None,**kw):
    if kwdict is None: kwdict = {}
    kwdict.update(kw)
    self.genericpf(kwdict,1,pfzx)
  def pfzy(self,kwdict=None,**kw):
    if kwdict is None: kwdict = {}
    kwdict.update(kw)
    self.genericpf(kwdict,0,pfzy)
  def pfxyg(self,kwdict=None,**kw):
    if kwdict is None: kwdict = {}
    kwdict.update(kw)
    self.genericpf(kwdict,2,pfxyg)
  def pfzxg(self,kwdict=None,**kw):
    if kwdict is None: kwdict = {}
    kwdict.update(kw)
    self.genericpf(kwdict,1,pfzxg)
  def pfzyg(self,kwdict=None,**kw):
    if kwdict is None: kwdict = {}
    kwdict.update(kw)
    self.genericpf(kwdict,0,pfzyg)


  def plphiz(self,ix=None,iy=None,color='fg',selfonly=0):
    if ix < self.fulllower[0]: return
    if iy < self.fulllower[1]: return
    if ix > self.fullupper[0]: return
    if iy > self.fullupper[1]: return
    plg(self.phi[ix-self.fulllower[0],iy-self.fulllower[1],1:-1],self.zmesh,
        color=color)
    if not selfonly:
      for child in self.children:
        child.plphiz(ix*child.refinement[0],iy*child.refinement[1],color=color)

  def plphix(self,iy=None,iz=None,color='fg',selfonly=0):
    if iy < self.fulllower[1]: return
    if iz < self.fulllower[2]: return
    if iy > self.fullupper[1]: return
    if iz > self.fullupper[2]: return
    plg(self.phi[:,iy-self.fulllower[1],iz-self.fulllower[2]+1],self.xmesh,
        color=color)
    if not selfonly:
      for child in self.children:
        child.plphix(iy*child.refinement[1],iz*child.refinement[2],color=color)

  def plphiy(self,ix=None,iz=None,color='fg',selfonly=0):
    if ix < self.fulllower[0]: return
    if iz < self.fulllower[2]: return
    if ix > self.fullupper[0]: return
    if iz > self.fullupper[2]: return
    plg(self.phi[ix-self.fulllower[0],:,iz-self.fulllower[2]+1],self.ymesh,
        color=color)
    if not selfonly:
      for child in self.children:
        child.plphiy(ix*child.refinement[0],iz*child.refinement[2],color=color)

  def plrhoz(self,ix=None,iy=None,color='fg',selfonly=0):
    if ix < self.fulllower[0]: return
    if iy < self.fulllower[1]: return
    if ix > self.fullupper[0]: return
    if iy > self.fullupper[1]: return
    plg(self.rho[ix-self.fulllower[0],iy-self.fulllower[1],:],self.zmesh,
        color=color)
    if not selfonly:
      for child in self.children:
        child.plrhoz(ix*child.refinement[0],iy*child.refinement[1],color=color)

  def plrhox(self,iy=None,iz=None,color='fg',selfonly=0):
    if iy < self.fulllower[1]: return
    if iz < self.fulllower[2]: return
    if iy > self.fullupper[1]: return
    if iz > self.fullupper[2]: return
    plg(self.rho[:,iy-self.fulllower[1],iz-self.fulllower[2]],self.xmesh,
        color=color)
    if not selfonly:
      for child in self.children:
        child.plrhox(iy*child.refinement[1],iz*child.refinement[2],color=color)

  def plrhoy(self,ix=None,iz=None,color='fg',selfonly=0):
    if ix < self.fulllower[0]: return
    if iz < self.fulllower[2]: return
    if ix > self.fullupper[0]: return
    if iz > self.fullupper[2]: return
    plg(self.rho[ix-self.fulllower[0],:,iz-self.fulllower[2]],self.ymesh,
        color=color)
    if not selfonly:
      for child in self.children:
        child.plrhoy(ix*child.refinement[0],iz*child.refinement[2],color=color)

  def plselfez(self,comp=2,ix=None,iy=None,color='fg',selfonly=0,withguard=1):
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
    plg(self.selfe[comp,ix-self.fulllower[0],iy-self.fulllower[1],iz],
        self.zmesh[iz],color=color)
    if not selfonly:
      for child in self.children:
        child.plselfez(comp,ix*child.refinement[0],iy*child.refinement[1],
                       color=color,withguard=withguard)

  def plselfex(self,comp=2,iy=None,iz=None,color='fg',selfonly=0,withguard=1):
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
    plg(self.selfe[comp,ix,iy-self.fulllower[1],iz-self.fulllower[2]],
                   self.xmesh[ix],color=color)
    if not selfonly:
      for child in self.children:
        child.plselfex(comp,iy*child.refinement[1],iz*child.refinement[2],
                       color=color,withguard=withguard)

  def plselfey(self,comp=2,ix=None,iz=None,color='fg',selfonly=0,withguard=1):
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
    plg(self.selfe[comp,ix-self.fulllower[0],iy,iz-self.fulllower[2]],
        self.ymesh[iy],color=color)
    if not selfonly:
      for child in self.children:
        child.plselfey(comp,ix*child.refinement[0],iz*child.refinement[2],
                       color=color,withguard=withguard)

  def drawbox(self,ip=None,idim=2,withguards=1,color=[],selfonly=0):
    if len(color)==0: color=['red', 'green', 'blue', 'cyan', 'magenta','yellow']
    if ip is None: ip = nint(-self.mins[idim]/self.deltas[idim])
    if ip < self.fulllower[idim] or ip > self.fullupper[idim]: return
    ii = [0,1,2]
    del ii[idim]
    if withguards:
      i01 = self.mins[ii[0]]
      i02 = self.maxs[ii[0]]
      i11 = self.mins[ii[1]]
      i12 = self.maxs[ii[1]]
    else:
      i01 = self.root.mins[ii[0]] + self.lower[ii[0]]*self.deltas[ii[0]]
      i02 = self.root.mins[ii[0]] + self.upper[ii[0]]*self.deltas[ii[0]]
      i11 = self.root.mins[ii[1]] + self.lower[ii[1]]*self.deltas[ii[1]]
      i12 = self.root.mins[ii[1]] + self.upper[ii[1]]*self.deltas[ii[1]]
    yy = [i01,i01,i02,i02,i01]
    xx = [i11,i12,i12,i11,i11]
    if idim==2:
      yy,xx = xx,yy
    else:
      xx = array(xx) + top.zbeam
    plg(yy,xx,color=color[0])
    if not selfonly:
      for child in self.children:
        child.drawbox(ip=ip*child.refinement[idim],idim=idim,
                      withguards=withguards,color=color[1:])

  def drawboxzy(self,ix=None,withguards=1,color=[],selfonly=0):
    self.drawbox(ip=ix,idim=0,withguards=withguards,color=color,
                 selfonly=selfonly)
  def drawboxzx(self,iy=None,withguards=1,color=[],selfonly=0):
    self.drawbox(ip=iy,idim=1,withguards=withguards,color=color,
                 selfonly=selfonly)
  def drawboxxy(self,iz=None,withguards=1,color=[],selfonly=0):
    self.drawbox(ip=iz,idim=2,withguards=withguards,color=color,
                 selfonly=selfonly)

  def drawfilledbox(self,ip=None,idim=2,withguards=1,ibox=None,selfonly=0):
    if ip is None: ip = nint(-self.mins[idim]/self.deltas[idim])
    if ip < self.fulllower[idim] or ip > self.fullupper[idim]: return
    ii = [0,1,2]
    del ii[idim]
    if withguards:
      i01 = self.mins[ii[0]]
      i02 = self.maxs[ii[0]]
      i11 = self.mins[ii[1]]
      i12 = self.maxs[ii[1]]
    else:
      i01 = self.root.mins[ii[0]] + self.lower[ii[0]]*self.deltas[ii[0]]
      i02 = self.root.mins[ii[0]] + self.upper[ii[0]]*self.deltas[ii[0]]
      i11 = self.root.mins[ii[1]] + self.lower[ii[1]]*self.deltas[ii[1]]
      i12 = self.root.mins[ii[1]] + self.upper[ii[1]]*self.deltas[ii[1]]
    xx = [i01,i01,i02,i02,i01]
    yy = [i11,i12,i12,i11,i11]
    if idim==2: xx,yy = yy,xx
    if ibox is None: ibox = ones(1,'b')
    else:            ibox = (ibox+1).astype('b')
    plfp(ibox,yy,xx,[5])
    if not selfonly:
      for child in self.children:
        child.drawfilledbox(ip=ip*child.refinement[idim],idim=idim,
                            withguards=withguards,ibox=ibox)

  def createdxobject(self,kwdict={},**kw):
    """
Create DX object drawing the object.
  - withguards=1: when true, the guard cells are included in the bounding box
    """
    kw.update(kwdict)
    withguards = kw.get('withguards',1)
    xmin,xmax = self.xmmin,self.xmmax
    ymin,ymax = self.ymmin,self.ymmax
    zmin,zmax = self.zmmin,self.zmmax
    if not withguards:
      ng = self.nguard*self.refinement
      xmin,xmax = xmin+ng[0]*self.dx, xmax-ng[0]*self.dx
      ymin,ymax = ymin+ng[1]*self.dy, ymax-ng[1]*self.dy
      zmin,zmax = zmin+ng[2]*self.dz, zmax-ng[2]*self.dz
    dxlist = [viewboundingbox(xmin,xmax,ymin,ymax,zmin,zmax)]
    for child in self.children:
      dxlist.append(child.getdxobject(kwdict=kw))
    self.dxobject = DXCollection(*dxlist)
















  #===========================================================================
  def solve2down(self):
    dxsqi  = 1./self.dx**2
    dysqi  = 1./self.dy**2
    dzsqi  = 1./self.dz**2
    self.phisave[:,:,:] = self.phi
    cond_potmg(self.conductors.interior,
               self.nx,self.ny,self.nz,self.phisave,0,false,
               2,true)
    residual(self.nx,self.ny,self.nz,self.nzfull,dxsqi,dysqi,dzsqi,
             self.phisave,self.rhosave,self.res,
             0,self.bound0,self.boundnz,self.boundxy,
             self.l2symtry,self.l4symtry,
             self.mgparam,2,true,self.lcndbndy,self.icndbndy,self.conductors)
    self.rho[:,:,:] = self.res[:,:,1:-1]
    self.phi[:,:,:] = 0.
    print 1,self.res[10,10,10]

    for child in self.children:
      child.solve2down()

    for parent in self.parents:
      s1 = maximum(self.fulllower,parent.fulllower*self.refinement)
      s2 = minimum(self.fullupper,parent.fullupper*self.refinement)
      sx1,sy1,sz1 = s1 - self.fulllower
      sx2,sy2,sz2 = s2 - self.fulllower
      px1,py1,pz1 = s1/self.refinement - parent.fulllower
      px2,py2,pz2 = s2/self.refinement - parent.fulllower
      restrict3d(sx2-sx1,sy2-sy1,sz2-sz1,pz2-pz1,sz2-sz1,
                 self.res[sx1:sx2+1,sy1:sy2+1,sz1:sz2+3],
                 parent.res[px1:px2+1,py1:py2+1,pz1:pz2+3],
                 self.boundxy,
                 self.bound0,self.boundnz,self.bound0,self.boundnz,
                 0,0,self.l2symtry,self.l4symtry)

  def solve2up(self):
    dxsqi  = 1./self.dx**2
    dysqi  = 1./self.dy**2
    dzsqi  = 1./self.dz**2
    for parent in self.parents:
      s1 = maximum(self.fulllower,parent.fulllower*self.refinement)
      s2 = minimum(self.fullupper,parent.fullupper*self.refinement)
      sx1,sy1,sz1 = s1 - self.fulllower
      sx2,sy2,sz2 = s2 - self.fulllower
      px1,py1,pz1 = s1/self.refinement - parent.fulllower
      px2,py2,pz2 = s2/self.refinement - parent.fulllower
      expand3d(px2-px1,py2-py1,pz2-pz1,sz1-sz2,pz2-pz1,
               parent.phi[px1:px2+1,py1:py2+1,pz1:pz2+3],
               self.phi[sx1:sx2+1,sy1:sy2+1,sz1:sz2+3],
               self.boundxy,self.bound0,self.boundnz,0,0)

    #   for i in range(self.uppasses):
    #     self.sorpass3d(0,self.nx,self.ny,self.nz,self.nzfull,
    #                    self.phi,self.rho,self.rstar,
    #                    dxsqi,dysqi,dzsqi,self.linbend,
    #                    self.l2symtry,self.l4symtry,self.bendx,
    #                    self.bound0,self.boundnz,self.boundxy,self.mgparam,2,
    #                    self.lcndbndy,self.icndbndy,self.conductors)

    childrenserror = 0.
    for child in self.children:
      childerror = child.solve2up()
      childrenserror = max(childrenserror,childerror)

    add(self.phi,self.phisave,self.phi)
    print 2,self.phi[10,10,10]

    # --- When using residual correction form, the other planes do need
    # --- to be set when using other than Dirichlet boundaries since
    # --- those planes are only set with the error of phi.
    if self.bound0  == 1: self.phi[:,:,0] = self.phi[:,:,2]
    if self.boundnz == 1: self.phi[:,:,-1] = self.phi[:,:,-3]
    if self.bound0  == 2: self.phi[:,:,0] = self.phi[:,:,-3]
    if self.boundnz == 2: self.phi[:,:,-1] = self.phi[:,:,2]

    # --- Calculate the change in phi.
    subtract(self.phisave,self.phi,self.phisave)
    absolute(self.phisave,self.phisave)
    self.mgerror[0] = MA.maximum(self.phisave)
    print self.mgerror[0],childrenserror
    print 'err = ',self.mgerror[0]
    return max(childrenserror,self.mgerror[0])

  #===========================================================================
  def solve2init(self):
    # --- Create temp arrays
    self.phisave = fzeros(shape(self.phi),'d')
    self.bendx = fzeros(((self.nx+1)*(self.ny+1)),'d')

    # --- Initialize temporaries
    dxsqi  = 1./self.dx**2
    dysqi  = 1./self.dy**2
    dzsqi  = 1./self.dz**2
    reps0c = self.mgparam/(eps0*2.*(dxsqi+dysqi+dzsqi))
    rdel   = dzsqi/(dxsqi + dysqi + dzsqi)

    checkconductors(self.nx,self.ny,self.nz,self.nzfull,
                    self.dx,self.dy,self.dz,self.conductors,
                    top.my_index,top.nslaves,top.izfsslave,top.nzfsslave)

    # --- Preset rho to increase performance (reducing the number of
    # --- multiplies in the main SOR sweep loop).
    if not self.linbend:
      # --- Do the operation in place (to avoid temp arrays)
      multiply(self.rho,reps0c,self.rho)
    else:
      raise "Bends not yet supported"

    # --- Since using residual correction form, need to save the original rho.
    self.rhosave = self.rho + 0.
    self.res = fzeros(shape(self.phi),'d')

    for child in self.children:
      child.solve2init()

  #===========================================================================
  def solve2(self,iwhich=0):
    # --- No initialization needed
    if iwhich == 1: return

    self.solve2init()

    # --- Initialize temporaries
    dxsqi  = 1./self.dx**2
    dysqi  = 1./self.dy**2
    dzsqi  = 1./self.dz**2
    reps0c = self.mgparam/(eps0*2.*(dxsqi+dysqi+dzsqi))
    rdel   = dzsqi/(dxsqi + dysqi + dzsqi)

    # --- Main multigrid v-cycle loop. Calculate error each iteration since
    # --- very few iterations are done.
    self.mgiters[0] = 0
    self.mgerror[0] = 2.*self.mgtol + 1.
    while (self.mgerror[0] > self.mgtol and self.mgiters[0] < self.mgmaxiters):
      self.mgiters[0] = self.mgiters[0] + 1
 
      self.solve2down()

      # --- Do one vcycle.
      self.vcycle(0,self.nx,self.ny,self.nz,self.nzfull,
                  self.dx,self.dy,self.dz,self.phi,self.rho,
                  self.rstar,self.linbend,self.l2symtry,self.l4symtry,
                  self.bendx,
                  self.boundxy,self.bound0,self.boundnz,
                  self.mgparam,self.mgform,self.mgmaxlevels,
                  self.downpasses,self.uppasses,self.lcndbndy,
                  self.icndbndy,self.conductors)

      self.mgerror[0] = self.solve2up()

      #else
      # mgexchange_phi(nx,ny,nz,nzfull,phi,localb0,localbnz,0,
      #                my_index,nslaves,izfsslave,nzfsslave,
      #                whosendingleft,izsendingleft,
      #                whosendingright,izsendingright)
      # mgexchange_phi(nx,ny,nz,nzfull,phi,localb0,localbnz,-1,
      #                my_index,nslaves,izfsslave,nzfsslave,
      #                whosendingleft,izsendingleft,
      #                whosendingright,izsendingright)
      #endif

    # --- For Dirichlet boundary conditions, copy data into guard planes
    # --- For other boundary conditions, the guard planes are used during
    # --- the solve are so are already set.
    if (self.bound0 == 0): self.phi[:,:,0] = self.phi[:,:,1]
    if (self.boundnz == 0): self.phi[:,:,-1] = self.phi[:,:,-2]

    # --- Make a print out.
    if (self.mgerror[0] > self.mgtol):
      print "Multigrid: Maximum number of iterations reached"
    print ("Multigrid: Error converged to %11.3e in %4d v-cycles"%
           (self.mgerror[0],self.mgiters[0]))

    # --- If using residual correction form, restore saved rho
    self.rho[:,:,:] = self.rhosave

    # --- Restore rho
    if (not self.linbend):
      multiply(self.rho,1./reps0c,self.rho)


  #===========================================================================
  def solve2down1(self):
    dxsqi  = 1./self.dx**2
    dysqi  = 1./self.dy**2
    dzsqi  = 1./self.dz**2
    self.phisave[:,:,:] = self.phi
    cond_potmg(self.conductors.interior,
               self.nx,self.ny,self.nz,self.phisave,0,false,
               2,true)
    residual(self.nx,self.ny,self.nz,self.nzfull,dxsqi,dysqi,dzsqi,
             self.phisave,self.rhosave,self.res,
             0,self.bound0,self.boundnz,self.boundxy,
             self.l2symtry,self.l4symtry,
             self.mgparam,2,true,self.lcndbndy,self.icndbndy,self.conductors)
    self.rho[:,:,:] = self.res[:,:,1:-1]
    self.phi[:,:,:] = 0.

    for i in range(self.downpasses):
      self.sorpass3d(0,self.nx,self.ny,self.nz,self.nzfull,
                     self.phi,self.rho,self.rstar,
                     dxsqi,dysqi,dzsqi,self.linbend,
                     self.l2symtry,self.l4symtry,self.bendx,
                     self.bound0,self.boundnz,self.boundxy,self.mgparam,2,
                     self.lcndbndy,self.icndbndy,self.conductors)

    residual(self.nx,self.ny,self.nz,self.nzfull,dxsqi,dysqi,dzsqi,
             self.phi,self.rho,self.res,
             0,self.bound0,self.boundnz,self.boundxy,
             self.l2symtry,self.l4symtry,
             self.mgparam,2,false,
             self.lcndbndy,self.icndbndy,self.conductors)

    for child in self.children:
      child.solve2down()

    for parent in self.parents:
      s1 = maximum(self.fulllower,parent.fulllower*self.refinement)
      s2 = minimum(self.fullupper,parent.fullupper*self.refinement)
      sx1,sy1,sz1 = s1 - self.fulllower
      sx2,sy2,sz2 = s2 - self.fulllower
      px1,py1,pz1 = s1/self.refinement - parent.fulllower
      px2,py2,pz2 = s2/self.refinement - parent.fulllower
      restrict3d(sx2-sx1,sy2-sy1,sz2-sz1,pz2-pz1,sz2-sz1,
                 self.res[sx1:sx2+1,sy1:sy2+1,sz1:sz2+1],
                 parent.rho[px1:px2+1,py1:py2+1,pz1:pz2+1],
                 self.boundxy,
                 self.bound0,self.boundnz,self.bound0,self.boundnz,
                 0,0,self.l2symtry,self.l4symtry)

  def solve2up1(self):
    dxsqi  = 1./self.dx**2
    dysqi  = 1./self.dy**2
    dzsqi  = 1./self.dz**2
    for parent in self.parents:
      s1 = maximum(self.fulllower,parent.fulllower*self.refinement)
      s2 = minimum(self.fullupper,parent.fullupper*self.refinement)
      sx1,sy1,sz1 = s1 - self.fulllower
      sx2,sy2,sz2 = s2 - self.fulllower
      px1,py1,pz1 = s1/self.refinement - parent.fulllower
      px2,py2,pz2 = s2/self.refinement - parent.fulllower
      expand3d(px2-px1,py2-py1,pz2-pz1,sz2-sz2,pz2-pz1,
               parent.phi[px1:px2+1,py1:py2+1,pz1:pz2+1],
               self.phi[sx1:sx2+1,sy1:sy2+1,sz1:sz2+1],
               self.boundxy,self.bound0,self.boundnz,0,0)

    for i in range(self.uppasses):
      self.sorpass3d(0,self.nx,self.ny,self.nz,self.nzfull,
                     self.phi,self.rho,self.rstar,
                     dxsqi,dysqi,dzsqi,self.linbend,
                     self.l2symtry,self.l4symtry,self.bendx,
                     self.bound0,self.boundnz,self.boundxy,self.mgparam,2,
                     self.lcndbndy,self.icndbndy,self.conductors)

    childrenserror = 0.
    for child in self.children:
      childerror = child.solve2up()
      childrenserror = max(childrenserror,childerror)

    add(self.phi,self.phisave,self.phi)

    # --- When using residual correction form, the other planes do need
    # --- to be set when using other than Dirichlet boundaries since
    # --- those planes are only set with the error of phi.
    if self.bound0  == 1: self.phi[:,:,0] = self.phi[:,:,2]
    if self.boundnz == 1: self.phi[:,:,-1] = self.phi[:,:,-3]
    if self.bound0  == 2: self.phi[:,:,0] = self.phi[:,:,-3]
    if self.boundnz == 2: self.phi[:,:,-1] = self.phi[:,:,2]

    # --- Calculate the change in phi.
    subtract(self.phisave,self.phi,self.phisave)
    absolute(self.phisave,self.phisave)
    self.mgerror[0] = MA.maximum(self.phisave)
    print 'err= ',self.mgerror[0]
    return max(childrenserror,self.mgerror[0])

# --- This can only be done after MRBlock is defined.
try:
  psyco.bind(MRBlock)
except NameError:
  pass

