"""Implements adaptive mesh refinement in 3d
"""
from warp import *
from multigrid import MultiGrid
from pyOpenDX import Visualizable,DXCollection,viewboundingbox
import MA
import __main__

# ---------------------------------------------------------------------------
MRsolver = [None]
def registersolver(solver):
  """
Registers the solver to be used in the particle simulation.
 - solver: is the solver object. It must have the methods loadrho, solve, and
           fetche defined
  """
  MRsolver[0] = solver
def getregisteredsolver():
  return MRsolver[0]
def loadrhoMR():
  assert MRsolver[0] is not None,"No solver has been registered"
  MRsolver[0].loadrho()
def fieldsolMR():
  assert MRsolver[0] is not None,"No solver has been registered"
  MRsolver[0].solve()
def fetcheMR():
  assert MRsolver[0] is not None,"No solver has been registered"
  MRsolver[0].fetche()
def fetchphiMR():
  assert MRsolver[0] is not None,"No solver has been registered"
  MRsolver[0].fetchphi()
def initfieldsolver():
    if w3d.AMRlevels>0:
      import AMR
      AMRtree=AMR.AMRTree()
      __main__.__dict__['AMRtree'] = AMRtree
      gchange('AMR')
      w3d.AMRcoalescing=0.8
      if getcurrpkg()=='w3d' and w3d.solvergeom==w3d.XYZgeomMR:
        registersolver(AMRtree.blocks)
__main__.__dict__['loadrhoMR'] = loadrhoMR
__main__.__dict__['fieldsolMR'] = fieldsolMR
__main__.__dict__['fetcheMR'] = fetcheMR
__main__.__dict__['fetchphiMR'] = fetchphiMR
__main__.__dict__['initfieldsolver'] = initfieldsolver
# ---------------------------------------------------------------------------

#########################################################################
class MRBlock(MultiGrid,Visualizable):
  """
 - parent:
 - refinement=2:
 - lower,upper: extent of domain in relative to parent, in parents grid
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
  def __init__(self,parent=None,refinement=2,
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
    else:
      # --- Save the parent and the index number. These are saved in lists
      # --- since a block can have multiple parents.
      self.parents = [parent]
      self.ichild = [ichild]

    self.overlaps = []
    self.refinement = refinement
    self.nguard = nguard
    self.root = root
    self.conductorlist = []

    if parent is None:
      # --- For the root, the dimensions and extent of the grid should
      # --- be specified. If not, they will be taken from w3d.
      self.dims = dims
      self.mins = mins
      self.maxs = maxs
      self.totalrefinement = 1
      self.root = self
      
    else:

      self.totalrefinement = parent.totalrefinement*self.refinement
      self.deltas = parent.deltas/refinement
      self.rootdims = root.dims*self.totalrefinement

      if lower is None and upper is None:
        # --- The grid mins and maxs are input.
        self.mins = array(mins)
        self.maxs = array(maxs)
        self.lower = nint((self.mins - root.mins)/self.deltas)
        self.upper = nint((self.maxs - root.mins)/self.deltas)

      else:
        # --- The grid lower and upper bounds are input. The bounds are
        # --- relative to the root grid, but scaled by the total refinement.
        self.lower = array(lower)
        self.upper = array(upper)
        self.mins = root.mins + self.lower*self.deltas
        self.maxs = root.mins + self.upper*self.deltas

      # --- Now, extend the domain by the given number of guard cells. Checks
      # --- are made so that the domain doesn't extend beyond the original grid.
      self.fulllower = maximum(0,self.lower - nguard*refinement)
      self.fullupper = minimum(self.rootdims,self.upper + nguard*refinement)

      # --- Recalculate grid quantities, including the guard regions.
      self.dims = self.fullupper - self.fulllower
      self.mins = root.mins + self.fulllower*self.deltas
      self.maxs = root.mins + self.fullupper*self.deltas

      # --- First, just use same boundary conditions as root.
      self.bounds = root.bounds.copy()

      # --- Check if the mesh doesn't reach the edge of the root grid.
      # --- If not, switch to Dirichlet boundary.
      self.bounds[::2] = where(self.fulllower > 0,0,self.bounds[::2])
      self.bounds[1::2] = where(self.fullupper < self.rootdims,
                                0,self.bounds[1::2])
      self.l2symtry = root.l2symtry
      self.l4symtry = root.l4symtry

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
        self.addchild(l,u,r)

    self.createwarrayforrho()

  def createwarrayforrho(self):
    # --- Create weight array needed for rho deposition.
    r = self.refinement
    # --- Is linear falloff in the weights correct for r > 2?
    w0 = r - abs(1.*iota(-r+1,+r-1))
    w0 = w0/sum(w0)
    # --- Expand into 3 dimensions
    self.w = outerproduct(w0,outerproduct(w0,w0))
    self.w.shape = (2*r-1,2*r-1,2*r-1)

  def addchild(self,lower,upper,mins=None,maxs=None,refinement=2):
    """
Add a mesh refined block to this block.
  -lower,upper,mins,maxs,refinement=2: All have same meanings as for the
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

  def installconductor(self,conductor,dfill=top.largepos):
    MultiGrid.installconductor(self,conductor,dfill=dfill)
    for child in self.children:
      child.installconductor(conductor,dfill=dfill)

  def clearconductors(self):
    MultiGrid.clearconductors(self)
    for child in self.children:
      child.clearconductors()

  #--------------------------------------------------------------------------
  # --- The next several methods handle initialization that is done after
  # --- all blocks have been added.
  #--------------------------------------------------------------------------

  def finalize(self):
    # --- This should only be called at the top level.
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
    if not self.islastcall(): return
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
      ii[...] = where(ii==0.,-ichild,ii)

      # --- Set interior to positive child number.
      l = maximum(self.fulllower,child.lower/child.refinement)
      u = child.upper/child.refinement
      # --- If the child extends to the edge of the parent mesh, it claims the
      # --- grid points on the upper edges.
      u = u + where(u == self.fullupper,1,0)
      # --- The child claims its full interior area
      ii = self.getchilddomains(l,u)
      ii[...] = +ichild

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
      od[...] = where(unclaimed,0,od)

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
      r = child.refinement
      cl = maximum(child.fulllower/r,l)
      cu = minimum(child.fullupper/r,u)
      if sometrue(cl > cu): continue

      # --- Get childdomains in the intersection region, Wherever the
      # --- the refinement level is lower than the childs, force childdomains
      # --- to be negative.
      d = self.getchilddomains(cl,cu,1)
      nbc = self.getlocalarray(nbcells,cl,cu,fulllower=l)
      d[...] = where(nbc<level*r,-abs(d),d)

      # --- Stretch out the array so it has the refined cell size of the child
      nbcstretched = zeros(1+(cu-cl)*2)
      nbcstretched[ ::r, ::r, ::r] = nbc[:  ,:  ,:  ]
      nbcstretched[1::r, ::r, ::r] = nbc[:-1,:  ,:  ]
      nbcstretched[ ::r,1::r, ::r] = nbc[:  ,:-1,:  ]
      nbcstretched[1::r,1::r, ::r] = nbc[:-1,:-1,:  ]
      nbcstretched[ ::r, ::r,1::r] = nbc[:  ,:  ,:-1]
      nbcstretched[1::r, ::r,1::r] = nbc[:-1,:  ,:-1]
      nbcstretched[ ::r,1::r,1::r] = nbc[:  ,:-1,:-1]
      nbcstretched[1::r,1::r,1::r] = nbc[:-1,:-1,:-1]
      child.clearinactiveregions(nbcstretched,self,level*r)

  #--------------------------------------------------------------------------
  # --- The next several methods handle the charge density calculation.
  #--------------------------------------------------------------------------

  def loadrho(self,lzero=true):
    """
Loads the charge density from the particles. This should only be called for
the top level grid.
    """
    if lzero: self.zerorho()
    for i,n,q,w in zip(top.ins-1,top.nps,top.sq,top.sw):
      if n == 0: continue
      self.setrho(top.xp[i:i+n],top.yp[i:i+n],top.zp[i:i+n],top.uzp[i:i+n],q,w)
    self.accumulaterhofromsiblings()
    self.getrhofromsiblings()
    self.gatherrhofromchildren()
    self.makerhoperiodic()

  def zerorho(self):
    self.rho[...] = 0.
    for child in self.children:
      child.zerorho()

  def getnperchild(self,ichildsorted):
    """
Given a sorted list of child numbers, returns how many of each number
there are. This would probably be better done in compiled code. This should be
efficient timewise, but uses two extra full-size arrays.
    """
    ic1 = ichildsorted[1:] - ichildsorted[:-1]
    nn = compress(ic1 > 0,iota(len(ichildsorted-1)))
    result = zeros(1+len(self.children))
    ss = 0
    for n in nn:
      result[nint(ichildsorted[n-1])] = n - ss
      ss += result[nint(ichildsorted[n-1])]
    result[nint(ichildsorted[-1])] = len(ichildsorted) - sum(result)
    return result

  def setrho(self,x,y,z,uz,q,w):
    """
Given the list of particles, a charge and a weight, deposits the charge
density of the mesh structure.
    """
    if len(x) == 0: return
    if len(self.children) > 0:

      # --- Find out whether the particles are in the local domain or one of
      # --- the children's.
      ichild = zeros(len(x),'d')
      getgridngp3d(len(x),x,y,z-top.zgrid,ichild,
                   self.nx,self.ny,self.nz,self.childdomains,
                   self.xmmin,self.xmmax,self.ymmin,self.ymmax,
                   self.zmmin,self.zmmax,self.l2symtry,self.l4symtry)

      # --- Take absolute value since the childdomains is signed.
      # --- Can this be done in place to save memory?
      ichild = abs(ichild)

      # --- Sort the list of which domain the particles are in.
      ichildsorter = argsort(ichild)
      ichildsorted = take(ichild,ichildsorter)
      nperchild = self.getnperchild(ichildsorted)

      # --- Now use the same sorting on the particle quantities.
      x = take(x,ichildsorter)
      y = take(y,ichildsorter)
      z = take(z,ichildsorter)
      uz = take(uz,ichildsorter)

      # --- For each child, pass to it the particles in it's domain.
      i = nperchild[0]
      for child,n in zip(self.children,nperchild[1:]):
        child.setrho(x[i:i+n],y[i:i+n],z[i:i+n],uz[i:i+n],q,w)
        i = i + n

      # --- Now, get remaining particles in this domain
      x = x[:nperchild[0]]
      y = y[:nperchild[0]]
      z = z[:nperchild[0]]
      uz = uz[:nperchild[0]]

    # --- Deposit the particles in this domain
    MultiGrid.setrho(self,x,y,z,uz,q,w)

  def gatherrhofromchildren(self):
    # --- Do this only the first time this is called. This should only be
    # --- done once and since each parent requires that this be done
    # --- before it can get its rho from here, it must be done on the
    # --- first call.
    if not self.isfirstcall(): return

    # --- Loop over the children
    for child in self.children:

      # --- Make sure that the child has gathered rho from its children.
      child.gatherrhofromchildren()

      # --- Get coordinates of child relative to this domain
      r = child.refinement
      l = maximum(child.fulllower/r,self.fulllower)
      u = minimum(child.fullupper/r,self.fullupper)

      # --- Grid cells to gather from child along each axis, relative to grid
      # --- cells in the parent mesh.
      dims = iota(-r+1,r-1)

      # --- Loop over the three dimensions. The loops loop over all child grid
      # --- cells that contribute to a parent grid cell.
      for k in dims:
        for j in dims:
          for i in dims:
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
            multiply(rh,self.w[i+1,j+1,k+1],rh)
            # --- Adjust for symmetries
            if child.bounds[0] == 1 and i > 0 and l[0] == 0:
              rh[0,:,:] = 2.*rh[0,:,:]
            if child.bounds[1] == 1 and i < 0 and u[0] == self.rootdims[0]:
              rh[-1,:,:] = 2.*rh[-1,:,:]
            if child.bounds[2] == 1 and j > 0 and l[1] == 0:
              rh[:,0,:] = 2.*rh[:,0,:]
            if child.bounds[3] == 1 and j < 0 and u[1] == self.rootdims[1]:
              rh[:,-1,:] = 2.*rh[:,-1,:]
            if child.bounds[4] == 1 and k > 0 and l[2] == 0:
              rh[:,:,0] = 2.*rh[:,:,0]
            if child.bounds[5] == 1 and k < 0 and u[2] == self.rootdims[2]:
              rh[:,:,-1] = 2.*rh[:,:,-1]
            # --- Now add in the contribution (in place)
            add(prho,rh,prho)

  def accumulaterhofromsiblings(self):
    """
Loop over overlapping siblings and collect the rho from them from regions that
are owned by this instance.
    """
    # --- This should only be done once.
    if not self.islastcall(): return
    for child in self.children:
      child.accumulaterhofromsiblings()
    for other in self.overlaps:
      l = maximum(self.fulllower,other.fulllower)
      u = minimum(self.fullupper,other.fullupper)
      srho = self.getrho(l,u)
      orho = other.getrho(l,u)
      sown = self.getsiblingdomains(l,u)
      srho[...] = srho + where(sown,orho,0.)

  def getrhofromsiblings(self):
    """
Get rho from overlapping siblings where they own the region.
    """
    if not self.islastcall(): return
    for child in self.children:
      child.getrhofromsiblings()
    for other in self.overlaps:
      l = maximum(self.fulllower,other.fulllower)
      u = minimum(self.fullupper,other.fullupper)
      srho = self.getrho(l,u)
      orho = other.getrho(l,u)
      oown = other.getsiblingdomains(l,u)
      srho[...] = where(oown,orho,srho)

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
      r = self.refinement
      # --- Coordinates of mesh relative to parent's mesh location
      # --- and refinement. The minimum and maximum are needed in case
      # --- this mesh extends beyond the parent's.
      l = maximum(parent.fulllower*self.refinement,self.fulllower)
      u = minimum(parent.fullupper*self.refinement,self.fullupper)
      pl = l/self.refinement
      pu = u/self.refinement
      sphi = self.getphi(l,u)
      pphi = parent.getphi(pl,pu)
      if l[0] == self.fulllower[0]: sphi[ 0,::r,::r] = pphi[ 0,:,:]
      if u[0] == self.fullupper[0]: sphi[-1,::r,::r] = pphi[-1,:,:]
      if l[1] == self.fulllower[1]: sphi[::r, 0,::r] = pphi[:, 0,:]
      if u[1] == self.fullupper[1]: sphi[::r,-1,::r] = pphi[:,-1,:]
      if l[2] == self.fulllower[2]: sphi[::r,::r, 0] = pphi[:,:, 0]
      if u[2] == self.fullupper[2]: sphi[::r,::r,-1] = pphi[:,:,-1]

    # --- Do remaining points where interpolation is needed.
    # --- This can be done, now that all of the contributions from the parents
    # --- is done.
    sphi = self.phi[:,:,1:-1]
    self.setslice(sphi[ 0,:,:])
    self.setslice(sphi[-1,:,:])
    self.setslice(sphi[:, 0,:])
    self.setslice(sphi[:,-1,:])
    self.setslice(sphi[:,:, 0])
    self.setslice(sphi[:,:,-1])

  def setslice(self,slice):
    # --- Does interpolation from cells coincident with the parent to the
    # --- other cells.
    r = self.refinement
    slice[1::r,::r] = 0.5*(slice[ :-1:r,::r] + slice[2:  :r,::r])
    slice[::r,1::r] = 0.5*(slice[::r, :-1:r] + slice[::r,2:  :r])
    slice[1::r,1::r] = 0.25*(slice[ :-1:r, :-1:r] + slice[2:  :r, :-1:r] +
                             slice[ :-1:r,2:  :r] + slice[2:  :r,2:  :r])

  #--------------------------------------------------------------------------
  # --- Methods to fetch E-fields and potential
  #--------------------------------------------------------------------------

  def fetche(self):
    """
Fetches the E field. This should only be called at the root level grid.
    """
    self.fetchefrompositions(w3d.xfsapi,w3d.yfsapi,w3d.zfsapi,
                             w3d.exfsapi,w3d.eyfsapi,w3d.ezfsapi)

  def fetchefrompositions(self,x,y,z,ex,ey,ez):
    if len(x) == 0: return
    if len(self.children) > 0:
      # --- Find out whether the particles are in the local domain or one of
      # --- the children's.
      ichild = zeros(len(x),'d')
      getgridngp3d(len(x),x,y,z-top.zgridprv,ichild,
                   self.nx,self.ny,self.nz,self.childdomains,
                   self.xmmin,self.xmmax,self.ymmin,self.ymmax,
                   self.zmmin,self.zmmax,self.l2symtry,self.l4symtry)
      # --- Zero out places where childdomains < 0
      ichild = where(ichild < 0.,0.,ichild)
      
      for i in range(len(self.children)+1):
        # --- Get list of particles within the ith domain
        # --- Note that when i==0, the particles selected are the ones from
        # --- this instances domain.
        ii = compress(ichild == i,arange(len(x)))
        if len(ii) == 0: continue
        # --- Get positions of those particles.
        xc = take(x,ii)
        yc = take(y,ii)
        zc = take(z,ii)
        # --- Create temporary arrays to hold the E field
        n = len(xc)
        tex,tey,tez = zeros(n,'d'),zeros(n,'d'),zeros(n,'d')
        # --- Now get the field
        if i == 0:
          MultiGrid.fetchefrompositions(self,xc,yc,zc,tex,tey,tez)
        else:
          self.children[i-1].fetchefrompositions(xc,yc,zc,tex,tey,tez)
        # --- Put the E field into the passed in arrays
        put(ex,ii,tex)
        put(ey,ii,tey)
        put(ez,ii,tez)

    else:
      # --- Get e-field from this domain
      MultiGrid.fetchefrompositions(self,x,y,z,ex,ey,ez)

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
      # --- the children's.
      ichild = zeros(len(x),'d')
      getgridngp3d(len(x),x,y,z-top.zgridprv,ichild,
                   self.nx,self.ny,self.nz,self.childdomains,
                   self.xmmin,self.xmmax,self.ymmin,self.ymmax,
                   self.zmmin,self.zmmax,self.l2symtry,self.l4symtry)
      # --- Zero out places where childdomains < 0
      ichild = where(ichild < 0.,0.,ichild)
      
      for i in range(len(self.children)+1):
        # --- Get list of particles within the ith domain
        # --- Note that when i==0, the particles selected are the ones from
        # --- this instances domain.
        ii = compress(ichild == i,arange(len(x)))
        if len(ii) == 0: continue
        # --- Get positions of those particles.
        xc = take(x,ii)
        yc = take(y,ii)
        zc = take(z,ii)
        # --- Create temporary arrays to hold the potential
        tphi = zeros(len(xc),'d')
        # --- Now get the field
        if i == 0:
          MultiGrid.fetchphifrompositions(self,xc,yc,zc,tphi)
        else:
          self.children[i-1].fetchphifrompositions(xc,yc,zc,tphi)
        # --- Put the potential into the passed in arrays
        put(phi,ii,tphi)

    else:
      # --- Get e-field from this domain
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
    # --- This extra check is needed in case there is only one parent.
    if self.ncallsfromparents == len(self.parents):
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
    if not self.islastcall(): return
    memtot = product(self.dims + 1)
    for child in self.children:
      memtot = memtot + child.getmem()
    return memtot
      
  def getphi(self,lower,upper):
    # --- Note that this takes into account the guard cells in z.
    ix1,iy1,iz1 = lower - self.fulllower
    ix2,iy2,iz2 = upper - self.fulllower + 1 
    iz1 = iz1 + 1
    iz2 = iz2 + 1
    return self.phi[ix1:ix2,iy1:iy2,iz1:iz2]
  def getrho(self,lower,upper,r=1):
    ix1,iy1,iz1 = lower - self.fulllower
    ix2,iy2,iz2 = upper - self.fulllower + 1
    return self.rho[ix1:ix2:r,iy1:iy2:r,iz1:iz2:r]
  def getchilddomains(self,lower,upper,upperedge=0):
    if self.childdomains is None:
      self.childdomains = fzeros(1+self.dims,'d')
    ix1,iy1,iz1 = lower - self.fulllower
    ix2,iy2,iz2 = upper - self.fulllower + upperedge
    return self.childdomains[ix1:ix2,iy1:iy2,iz1:iz2]
  def getsiblingdomains(self,lower,upper,r=1):
    if self.siblingdomains is None:
      self.siblingdomains = -1*ones(1+self.dims)
    ix1,iy1,iz1 = lower - self.fulllower
    ix2,iy2,iz2 = upper - self.fulllower + 1
    return self.siblingdomains[ix1:ix2:r,iy1:iy2:r,iz1:iz2:r]
  def getlocalarray(self,array,lower,upper,r=1,fulllower=None):
    if fulllower is None: fulllower = self.fulllower
    ix1,iy1,iz1 = lower - fulllower
    ix2,iy2,iz2 = upper - fulllower + 1
    return array[ix1:ix2:r,iy1:iy2:r,iz1:iz2:r]

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
      if idim == 0: ip = kw.get('ix',None)
      if idim == 1: ip = kw.get('iy',None)
      if idim == 2: ip = kw.get('iz',None)
      if ip is None: ip = nint(-self.mins[idim]/self.deltas[idim])
    else:
      ip = ip*self.refinement
      if idim == 0: kw['ix'] = ip - self.fulllower[0]
      if idim == 1: kw['iy'] = ip - self.fulllower[1]
      if idim == 2: kw['iz'] = ip - self.fulllower[2]
      kw['cmin'] = ppgeneric.cmin
      kw['cmax'] = ppgeneric.cmax
    # --- Only make the plot if the plane is included in the domain.
    # --- Even if there is no overlap, the children must be called since
    # --- they may overlap (in the domain of a different parent).
    if self.fulllower[idim] <= ip and ip <= self.fullupper[idim]:
      if self.childdomains is not None:
        # --- Create the ireg array, which will be set to zero in the domain
        # --- of the children.
        ss = list(shape(self.rho))
        del ss[idim]
        ireg = zeros(ss)
        if idim==0:
          ireg[1:,1:]=equal(self.childdomains[ip-self.fulllower[0],:-1,:-1],0)
        if idim==1:
          ireg[1:,1:]=equal(self.childdomains[:-1,ip-self.fulllower[1],:-1],0)
        if idim==2:
          ireg[1:,1:]=equal(self.childdomains[:-1,:-1,ip-self.fulllower[2]],0)
        if idim != 2: ireg = transpose(ireg)
        kw['ireg'] = ireg
      else:
        kw['ireg'] = None
      MultiGrid.genericpf(self,kw,pffunc)
      kw['titles'] = 0
      kw['lcolorbar'] = 0

    for child in self.children:
      child.genericpf(kw,idim,pffunc,ip)

  def pfxy(self,**kw): self.genericpf(kw,2,pfxy)
  def pfzx(self,**kw): self.genericpf(kw,1,pfzx)
  def pfzy(self,**kw): self.genericpf(kw,0,pfzy)
  def pfxyg(self,**kw): self.genericpf(kw,2,pfxyg)
  def pfzxg(self,**kw): self.genericpf(kw,1,pfzxg)
  def pfzyg(self,**kw): self.genericpf(kw,0,pfzyg)


  def plphiz(self,ix=None,iy=None,color='fg'):
    if ix < self.fulllower[0]: return
    if iy < self.fulllower[1]: return
    if ix > self.fullupper[0]: return
    if iy > self.fullupper[1]: return
    plg(self.phi[ix-self.fulllower[0],iy-self.fulllower[1],1:-1],self.zmesh,
        color=color)
    for child in self.children:
      child.plphiz(ix*child.refinement,iy*child.refinement,color=color)

  def plphix(self,iy=None,iz=None,color='fg'):
    if iy < self.fulllower[1]: return
    if iz < self.fulllower[2]: return
    if iy > self.fullupper[1]: return
    if iz > self.fullupper[2]: return
    plg(self.phi[:,iy-self.fulllower[1],iz-self.fulllower[2]+1],self.xmesh,
        color=color)
    for child in self.children:
      child.plphix(iy*child.refinement,iz*child.refinement,color=color)

  def plphiy(self,ix=None,iz=None,color='fg'):
    if ix < self.fulllower[0]: return
    if iz < self.fulllower[2]: return
    if ix > self.fullupper[0]: return
    if iz > self.fullupper[2]: return
    plg(self.phi[ix-self.fulllower[0],:,iz-self.fulllower[2]+1],self.ymesh,
        color=color)
    for child in self.children:
      child.plphiy(ix*child.refinement,iz*child.refinement,color=color)

  def plrhoz(self,ix=None,iy=None,color='fg'):
    if ix < self.fulllower[0]: return
    if iy < self.fulllower[1]: return
    if ix > self.fullupper[0]: return
    if iy > self.fullupper[1]: return
    plg(self.rho[ix-self.fulllower[0],iy-self.fulllower[1],:],self.zmesh,
        color=color)
    for child in self.children:
      child.plrhoz(ix*child.refinement,iy*child.refinement,color=color)

  def plrhox(self,iy=None,iz=None,color='fg'):
    if iy < self.fulllower[1]: return
    if iz < self.fulllower[2]: return
    if iy > self.fullupper[1]: return
    if iz > self.fullupper[2]: return
    plg(self.rho[:,iy-self.fulllower[1],iz-self.fulllower[2]],self.xmesh,
        color=color)
    for child in self.children:
      child.plrhox(iy*child.refinement,iz*child.refinement,color=color)

  def plrhoy(self,ix=None,iz=None,color='fg'):
    if ix < self.fulllower[0]: return
    if iz < self.fulllower[2]: return
    if ix > self.fullupper[0]: return
    if iz > self.fullupper[2]: return
    plg(self.rho[ix-self.fulllower[0],:,iz-self.fulllower[2]],self.ymesh,
        color=color)
    for child in self.children:
      child.plrhoy(ix*child.refinement,iz*child.refinement,color=color)

  def drawbox(self,ip=None,idim=2,withguards=1,color=[]):
    if len(color)==0: color=['red', 'green', 'blue', 'cyan', 'magenta','yellow']
    if ip is None: ip = self.dims[idim]/2
    if ip < self.fulllower[idim] or ip > self.fullupper[idim]: return
    ii = [0,1,2]
    del ii[idim]
    i01 = self.mins[ii[0]]
    i02 = self.maxs[ii[0]]
    i11 = self.mins[ii[1]]
    i12 = self.maxs[ii[1]]
    if not withguards:
      i01 = i01 + self.deltas[ii[0]]*self.nguard*self.refinement
      i02 = i02 - self.deltas[ii[0]]*self.nguard*self.refinement
      i11 = i11 + self.deltas[ii[1]]*self.nguard*self.refinement
      i12 = i12 - self.deltas[ii[1]]*self.nguard*self.refinement
    xx = [i01,i01,i02,i02,i01]
    yy = [i11,i12,i12,i11,i11]
    if idim==2: xx,yy = yy,xx
    plg(xx,yy,color=color[0])
    for child in self.children:
      child.drawbox(ip=ip*child.refinement,idim=idim,withguards=withguards,
                    color=color[1:])

  def drawboxzy(self,ix=None,withguards=1,color=[]):
    self.drawbox(ip=ix,idim=0,withguards=withguards,color=[])
  def drawboxzx(self,iy=None,withguards=1,color=[]):
    self.drawbox(ip=iy,idim=1,withguards=withguards,color=[])
  def drawboxxy(self,iz=None,withguards=1,color=[]):
    self.drawbox(ip=iz,idim=2,withguards=withguards,color=[])

  def drawfilledbox(self,ip=None,idim=2,withguards=1,ibox=None):
    if ip is None: ip = self.dims[idim]/2
    if ip < self.fulllower[idim] or ip > self.fullupper[idim]: return
    ii = [0,1,2]
    del ii[idim]
    i01 = self.mins[ii[0]]
    i02 = self.maxs[ii[0]]
    i11 = self.mins[ii[1]]
    i12 = self.maxs[ii[1]]
    if not withguards:
      i01 = i01 + self.dims[ii[0]]*self.nguard*self.refinement
      i02 = i02 - self.dims[ii[0]]*self.nguard*self.refinement
      i11 = i11 + self.dims[ii[1]]*self.nguard*self.refinement
      i12 = i12 - self.dims[ii[1]]*self.nguard*self.refinement
    xx = [i01,i01,i02,i02,i01]
    yy = [i11,i12,i12,i11,i11]
    if idim==2: xx,yy = yy,xx
    if ibox is None: ibox = ones(1,'b')
    else:            ibox = (ibox+1).astype('b')
    plfp(ibox,yy,xx,[5])
    for child in self.children:
      child.drawfilledbox(ip=ip*child.refinement,idim=idim,
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
      xmin,xmax = xmin+ng*self.dx, xmax-ng*self.dx
      ymin,ymax = ymin+ng*self.dy, ymax-ng*self.dy
      zmin,zmax = zmin+ng*self.dz, zmax-ng*self.dz
    dxlist = [viewboundingbox(xmin,xmax,ymin,ymax,zmin,zmax)]
    for child in self.children:
      dxlist.append(child.getdxobject(kwdict=kw))
    self.dxobject = DXCollection(*dxlist)

