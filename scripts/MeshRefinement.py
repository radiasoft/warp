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
  MRsolver[0] = solver
def loadrhoMR():
  assert MRsolver[0] is not None,"No solver has been registered"
  MRsolver[0].loadrho()
def fieldsolMR():
  assert MRsolver[0] is not None,"No solver has been registered"
  MRsolver[0].solve()
def fetcheMR():
  assert MRsolver[0] is not None,"No solver has been registered"
  MRsolver[0].fetche()
__main__.__dict__['loadrhoMR'] = loadrhoMR
__main__.__dict__['fieldsolMR'] = fieldsolMR
__main__.__dict__['fetcheMR'] = fetcheMR
# ---------------------------------------------------------------------------

class BlockOverlap:
  def __init__(self,block,fulllower,fullupper,otherowns):
    self.block = block
    self.fulllower = fulllower
    self.fullupper = fullupper
    self.otherowns = otherowns

class MRBlock(MultiGrid,Visualizable):
  """
 - parent:
 - refinement=2:
 - lower,upper: extent of domain in absolute, in parents grid cell size
 - children: list of tuples, each containing three elements,
             (lower,upper,refinement).
  """
  def __init__(self,parent=None,refinement=2,
                    lower=None,upper=None,
                    ichild=None,
                    dims=None,mins=None,maxs=None,
                    nguard=1,dimsmax=None,
                    children=None,**kw):

    if parent is None:
      self.parents = []
      self.ichild = []
    else:
      self.parents = [parent]
      self.ichild = [ichild]
    self.overlaps = []
    self.siblings = []
    self.refinement = refinement
    self.nguard = nguard
    self.dimsmax = dimsmax
    self.conductors = []

    if parent is None:
      # --- All grid sizing parameters should be input, otherwise they
      # --- will be taken from w3d.
      self.dims = dims
      self.mins = mins
      self.maxs = maxs
      
    else:

      self.deltas = parent.deltas/refinement

      # --- For now, just use Dirichlet boundaries for all submeshes.
      # --- A future upgrade will allow the proper boundaries for a submesh
      # --- at the edge of the base mesh.
      self.bound0 = dirichlet
      self.boundnz = dirichlet
      self.boundxy = dirichlet

      if lower is None and upper is None:
        # --- The grid mins and maxs are input.
        self.mins = array(mins)
        self.maxs = array(maxs)
        self.dims = (self.maxs - self.mins)/self.deltas
        self.lower = nint((self.mins - parent.mins)/self.deltas) + refinement*parent.lower
        self.upper = nint((self.maxs - parent.mins)/self.deltas) + refinement*parent.lower

      else:
        # --- The grid lower and upper bounds are input. The bounds are
        # --- relative to the parent. They are converted to be relative
        # --- to this instance.
        # --- This is the standard method of specifying a child block.
        self.lower = refinement*(parent.lower + array(lower))
        self.upper = refinement*(parent.lower + array(upper))
        self.dims = self.upper - self.lower
        self.mins = parent.mins+(self.lower-refinement*parent.lower)*self.deltas
        self.maxs = parent.mins+(self.upper-refinement*parent.lower)*self.deltas

      # --- Now, extend the domain by the given number of guard cells. Checks
      # --- are made so that the domain doesn't extend beyond the original grid.
      dimtoolow = minimum(0,self.lower - nguard*refinement)
      dimtoohigh = maximum(dimsmax,self.upper + nguard*refinement) - dimsmax
      self.fulllower = self.lower - nguard*refinement - dimtoolow
      self.fullupper = self.upper + nguard*refinement - dimtoohigh

      # --- Recalculate grid quantities
      self.dims = self.fullupper - self.fulllower
      self.mins = (parent.mins +
                   (self.fulllower-refinement*parent.fulllower)*self.deltas)
      self.maxs = (parent.mins +
                   (self.fullupper-refinement*parent.fulllower)*self.deltas)

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

    MultiGrid.__init__(self,**kw)

    if parent is None:
      # --- This is only needed by the top grid in cases when the grid
      # --- parameters are obtained from w3d instead of the argument list.
      self.dims = array([self.nx,self.ny,self.nz])
      self.deltas = array([self.dx,self.dy,self.dz])
      self.mins = array([self.xmmin,self.ymmin,self.zmmin])
      self.maxs = array([self.xmmax,self.ymmax,self.zmmax])
      self.lower = zeros(3)
      self.upper = self.dims
      self.fulllower = zeros(3)
      self.fullupper = self.dims
      self.dimsmax = self.dims

    # --- idomains is the cell centered grid which keeps track of which
    # --- cells are owned by which children. If there are no children,
    # --- then it is not needed.
    self.idomains = None

    # --- Now add any specified children
    self.children = []
    if children is not None:
      for l,u,r in children:
        self.addchild(l,u,r)

    # --- Create some work data arrays needed.
    # --- For rho deposition
    w0 = array([0.25,0.5,0.25])
    self.w = outerproduct(w0,outerproduct(w0,w0))
    self.w.shape = (3,3,3)
    # --- For rho deposition and field solving.
    self.ncallsfromparents = 0

  def addchild(self,lower,upper,mins=None,maxs=None,refinement=2):
    child = MRBlock(parent=self,lower=lower,upper=upper,mins=mins,maxs=maxs,
                    refinement=refinement,ichild=len(self.children)+1,
                    nguard=self.nguard,dimsmax=self.dimsmax*refinement)
    self.addblockaschild(child)

  def addblockaschild(self,block):
    if self.idomains is None:
      self.idomains = fzeros((self.nx+1,self.ny+1,self.nz+1),'d')
    self.children.append(block)
    # --- Set full domain to negative of child number first.
    ix1,iy1,iz1 = maximum(0,block.fulllower/block.refinement - self.fulllower)
    ix2,iy2,iz2 =           block.fullupper/block.refinement - self.fulllower
    # --- If the child extends to the edge of the parent mesh, it claims the
    # --- grid points on the upper edges.
    if ix2 == self.dims[0]: ix2 += 1
    if iy2 == self.dims[1]: iy2 += 1
    if iz2 == self.dims[2]: iz2 += 1
    ii = self.idomains[ix1:ix2,iy1:iy2,iz1:iz2]
    ii[:,:,:] = where(ii==0.,-len(self.children),ii)
    # --- Set interior to positive child number.
    ix1,iy1,iz1 = maximum(0,block.lower/block.refinement - self.fulllower)
    ix2,iy2,iz2 =           block.upper/block.refinement - self.fulllower
    # --- If the child extends to the edge of the parent mesh, it claims the
    # --- grid points on the upper edges.
    if ix2 == self.dims[0]: ix2 += 1
    if iy2 == self.dims[1]: iy2 += 1
    if iz2 == self.dims[2]: iz2 += 1
    self.idomains[ix1:ix2,iy1:iy2,iz1:iz2] = +len(self.children)

  def setrhoboundaries(self):
    """
Sets rho on the boundaries, using the values from the parent grid
    """
    self.setgridboundaries('rho')
    for child in self.children:
      child.setrhoboundaries()

  def setgridboundaries(self,grid):
    """
Sets grid on the boundaries, using the values from the parent grid
    """
    if len(self.parents) == 0: return
    for parent in self.parents:
      r = self.refinement
      # --- Coordinates of mesh relative to parent's mesh location
      # --- and refinement. The minimum and maximum are needed in case
      # --- this mesh extends beyond the parent's.
      i1 = maximum(0,self.fulllower/self.refinement - parent.fulllower)
      i2 = minimum(parent.dims,
                   self.fullupper/self.refinement - parent.fulllower)
      # --- Coordinates of this mesh that covers the parent's mesh.
      # --- If this mesh is completely within the parents mesh,
      # --- then i3 is zero and i4 is self.dims.
      i3 = (parent.fulllower + i1)*self.refinement - self.fulllower
      i4 = (parent.fulllower + i2)*self.refinement - self.fulllower
      ix1,iy1,iz1 = i1
      ix2,iy2,iz2 = i2
      ix3,iy3,iz3 = i3
      ix4,iy4,iz4 = i4
      # --- First get points directly copied from parent.
      # --- As a convenience, get references to the grid arrays without
      # --- the z guard planes.
      if grid == 'phi':
        sgrid = getattr(self,grid)[:,:,1:-1]
        pgrid = getattr(parent,grid)[:,:,1:-1]
      else:
        sgrid = getattr(self,grid)
        pgrid = getattr(parent,grid)
      # --- Note that if the child extends beyond the parent's domain,
      # --- then certain boundary planes are completely skipped.
      if self.fulllower[0]/self.refinement >= parent.fulllower[0]:
        sgrid[ix3,iy3:iy4+1:r,iz3:iz4+1:r] = pgrid[ix1,iy1:iy2+1,iz1:iz2+1]
      if self.fullupper[0]/self.refinement <= parent.fullupper[0]:
        sgrid[ix4,iy3:iy4+1:r,iz3:iz4+1:r] = pgrid[ix2,iy1:iy2+1,iz1:iz2+1]
      if self.fulllower[1]/self.refinement >= parent.fulllower[1]:
        sgrid[ix3:ix4+1:r,iy3,iz3:iz4+1:r] = pgrid[ix1:ix2+1,iy1,iz1:iz2+1]
      if self.fullupper[1]/self.refinement <= parent.fullupper[1]:
        sgrid[ix3:ix4+1:r,iy4,iz3:iz4+1:r] = pgrid[ix1:ix2+1,iy2,iz1:iz2+1]
      if self.fulllower[2]/self.refinement >= parent.fulllower[2]:
        sgrid[ix3:ix4+1:r,iy3:iy4+1:r,iz3] = pgrid[ix1:ix2+1,iy1:iy2+1,iz1]
      if self.fullupper[2]/self.refinement <= parent.fullupper[2]:
        sgrid[ix3:ix4+1:r,iy3:iy4+1:r,iz4] = pgrid[ix1:ix2+1,iy1:iy2+1,iz2]

    # --- Do remaining points where interpolation is needed.
    # --- This can be done, now that all of the contributions from the parents
    # --- is done.
    if grid == 'phi':
      sgrid = getattr(self,grid)[:,:,1:-1]
    else:
      sgrid = getattr(self,grid)
    self.setslice(sgrid[ 0,:,:])
    self.setslice(sgrid[-1,:,:])
    self.setslice(sgrid[:, 0,:])
    self.setslice(sgrid[:,-1,:])
    self.setslice(sgrid[:,:, 0])
    self.setslice(sgrid[:,:,-1])

  def setslice(self,slice):
    # --- Does interpolation from cells coincident with the parent to the
    # --- other cells.
    r = self.refinement
    slice[1::r,::r] = 0.5*(slice[ :-1:r,::r] + 
                           slice[2:  :r,::r])
    slice[::r,1::r] = 0.5*(slice[::r, :-1:r] + 
                           slice[::r,2:  :r])
    slice[1::r,1::r] = 0.25*(slice[ :-1:r, :-1:r] +
                             slice[2:  :r, :-1:r] +
                             slice[ :-1:r,2:  :r] +
                             slice[2:  :r,2:  :r])

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
    if len(self.children) > 0:
      # --- Find out whether the particles are in the local domain or one of
      # --- the children's.
      ichild = zeros(len(x),'d')
      getgridngp3d(len(x),x,y,z-top.zgrid,ichild,
                   self.nx,self.ny,self.nz,self.idomains,
                   self.xmmin,self.xmmax,self.ymmin,self.ymmax,
                   self.zmmin,self.zmmax,self.l2symtry,self.l4symtry)
      # --- Take absolute value since the idomains is signed.
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
        if n > 0:
          child.setrho(x[i:i+n],y[i:i+n],z[i:i+n],uz[i:i+n],q,w)
        i = i + n
      # --- Now, get remaining particles in this domain
      x = x[:nperchild[0]]
      y = y[:nperchild[0]]
      z = z[:nperchild[0]]
      uz = uz[:nperchild[0]]
    # --- Deposit the particles in this domain
    MultiGrid.setrho(self,x,y,z,uz,q,w)

  def zerorho(self):
    self.rho[...] = 0.
    for child in self.children:
      child.zerorho()

  def gatherrhofromchildren(self):
    # --- Do this only the first time this is called. This should only be
    # --- done once and since each parent requires that this be done
    # --- before it can get its rho from here, it must be done on the
    # --- first call.
    self.ncallsfromparents += 1
    if self.ncallsfromparents > 1:
      if self.ncallsfromparents == len(self.parents): self.ncallsfromparents = 0
      return
    if self.ncallsfromparents == len(self.parents): self.ncallsfromparents = 0
    # --- Loop over the children
    for child,ichild in zip(self.children,range(1,1+len(self.children))):
      # --- Make sure that the child has gathered rho from its children.
      child.gatherrhofromchildren()
      # --- Get coordinates of child relative to this domain
      r = child.refinement
      i1 = maximum(0,child.fulllower/r - self.fulllower)
      i2 = minimum(self.dims,child.fullupper/r - self.fulllower)
      # --- Grid cells to gather from child along each axis, relative to grid
      # --- cells in the parent mesh.
      dims = iota(-r+1,r-1)
      # --- Loop over the three dimensions. The loops loop over all child grid
      # --- cells that contribute to a parent grid cell.
      for k in dims:
        for j in dims:
          for i in dims:
            ii1 = array([i,j,k]) + (i1 + self.fulllower)*r - child.fulllower
            ii2 = array([i,j,k]) + (i2 + self.fulllower)*r - child.fulllower
            # --- Make sure points outside the child's domain are skipped.
            # --- If the ii1 are < 0, then skip the lowest plane since
            # --- it is outside the childs domain. The same for the upper end.
            pm = (ii1 < 0)
            pp = -(ii2 > child.dims)
            # --- Gather pointers to the blocks of data to be operated on
            # --- First, self's rho.
            ix1,iy1,iz1 = i1 + pm
            ix2,iy2,iz2 = i2 + pp + 1
            p = self.rho[ix1:ix2,iy1:iy2,iz1:iz2]
            # --- Then the child's rho
            cx1,cy1,cz1 = ii1 + r*pm
            cx2,cy2,cz2 = ii2 + r*pp + 1
            c = child.rho[cx1:cx2:r,cy1:cy2:r,cz1:cz2:r]
            # --- Note the differing indexing from self.rho. For p(xyz)m==0,
            # --- this has no affect. But or p(xyz)m==1, (i.e. i,j, or k=-1),
            # --- this is needed since, if the child owns the cell that that
            # --- falls in, it will contribute it to the next cell up even
            # --- if its doesn't own it.
            ix1,iy1,iz1 = i1
            ix2,iy2,iz2 = i2 + pp - pm + 1
            d = abs(self.idomains[ix1:ix2,iy1:iy2,iz1:iz2])
            ld = (d == ichild)
            # --- For the upper edges of the childs domain, if no other
            # --- children own the cells, this child must contribute rho to
            # --- there.
            if i <= 0: ld[-1,:,:] = where(d[-1,:,:]==0.,1,ld[-1,:,:])
            if j <= 0: ld[:,-1,:] = where(d[:,-1,:]==0.,1,ld[:,-1,:])
            if k <= 0: ld[:,:,-1] = where(d[:,:,-1]==0.,1,ld[:,:,-1])
            # --- Get the rho that will be contributed.
            rh = where(ld,c,0.)
            # --- Multiply by the weight in place.
            multiply(rh,self.w[i+1,j+1,k+1],rh)
            # --- Now add in the contribution (in place)
            add(p,rh,p)

  def addmyrhotosiblings(self):
    self.ncallsfromparents += 1
    if self.ncallsfromparents < len(self.parents): return
    self.ncallsfromparents = 0
    for sibling in self.overlaps:
      i1 = sibling.fulllower - self.fulllower
      i2 = sibling.fullupper - self.fulllower
      i3 = sibling.fulllower - sibling.block.fulllower
      i4 = sibling.fullupper - sibling.block.fulllower
      myrho = self.rho[i1[0]:i2[0]+1,i1[1]:i2[1]+1,i1[2]:i2[2]+1]
      nrho = sibling.block.rho[i3[0]:i4[0]+1,i3[1]:i4[1]+1,i3[2]:i4[2]+1]
      add(where(sibling.otherowns,myrho,0.),nrho,nrho)
    for child in self.children:
      child.addmyrhotosiblings()

  def getrhofromsiblings(self):
    self.ncallsfromparents += 1
    if self.ncallsfromparents < len(self.parents): return
    self.ncallsfromparents = 0
    for sibling in self.overlaps:
      i1 = sibling.fulllower - self.fulllower
      i2 = sibling.fullupper - self.fulllower
      i3 = sibling.fulllower - sibling.block.fulllower
      i4 = sibling.fullupper - sibling.block.fulllower
      myrho = self.rho[i1[0]:i2[0]+1,i1[1]:i2[1]+1,i1[2]:i2[2]+1]
      nrho = sibling.block.rho[i3[0]:i4[0]+1,i3[1]:i4[1]+1,i3[2]:i4[2]+1]
      myrho[:,:,:] = where(sibling.otherowns,nrho,myrho)
    for child in self.children:
      child.getrhofromsiblings()

  def findoverlappingsiblings(self):
    # --- Loop over parents, checking for overlaps with all of the children
    # --- of each.
    for parent in self.parents:
      for sibling in parent.children:
        if sibling == self: continue
        if sibling in self.siblings: continue
        self.siblings.append(sibling)
        self.checkifsiblingoverlaps(sibling)
    for child in self.children:
      child.findoverlappingsiblings()

  def checkifsiblingoverlaps(self,other):
    # --- Get extent of possible overlapping domain
    l = maximum(other.fulllower,self.fulllower)
    u = minimum(other.fullupper,self.fullupper)
    # --- If any of the lengths are negative, then there is no overlap.
    # --- Don't do anything else.
    doesnotoverlap = sometrue(u < l)
    if doesnotoverlap: return
    # --- Check how much of the overlap region is owned by the other instance
    # --- First, assume that it doesn't own any.
    otherowns = zeros(u-l+1)
    pl = l/self.refinement
    pu = u/self.refinement
    r = self.refinement
    # --- Check the idomains for each parent of the other to find out the
    # --- regions owned by that instance.
    for parent,ichild in zip(other.parents,other.ichild):
      pld = maximum(pl,parent.fulllower)
      pud = minimum(pu,parent.fullupper)
      ix1,iy1,iz1 = pld*other.refinement - l
      ix2,iy2,iz2 = pud*other.refinement - l + 1
      ix3,iy3,iz3 = pld - parent.fulllower
      ix4,iy4,iz4 = pud - parent.fulllower + 1
      otherownsslice = otherowns[ix1:ix2,iy1:iy2,iz1:iz2]
      idomainsslice = abs(parent.idomains[ix3:ix4,iy3:iy4,iz3:iz4])
      # --- For the cells owned by that instance, set otherowns to one.
      iii = where(idomainsslice==ichild,1,otherownsslice[::r,::r,::r])
      # --- That array is relative to the parents mesh size. Copy to each
      # --- of the eight nodes in the lower corner of the parent cell.
      # --- The other nodes are part of different cells in the parent.
      otherownsslice[ ::r, ::r, ::r] = iii
      otherownsslice[1::r, ::r, ::r] = iii[:-1,:  ,:  ]
      otherownsslice[ ::r,1::r, ::r] = iii[:  ,:-1,:  ]
      otherownsslice[1::r,1::r, ::r] = iii[:-1,:-1,:  ]
      otherownsslice[ ::r, ::r,1::r] = iii[:  ,:  ,:-1]
      otherownsslice[1::r, ::r,1::r] = iii[:-1,:  ,:-1]
      otherownsslice[ ::r,1::r,1::r] = iii[:  ,:-1,:-1]
      otherownsslice[1::r,1::r,1::r] = iii[:-1,:-1,:-1]
    # --- Only keep track of overlaps where the other instance owns some of
    # --- the overlap region. If the other instance does not own the any of
    # --- the overlap region, then nothing special needs to be done so there
    # --- is no reason to keep track of the overlap.
    if maxnd(otherowns) == 1:
      self.overlaps.append(BlockOverlap(other,l,u,otherowns))

  def findallparents(self,blocklists):
    for block in blocklists[0]:
      if block in self.parents: continue
      # --- Get extent of possible overlapping domain
      l = maximum(block.fulllower,self.fulllower/self.refinement)
      u = minimum(block.fullupper,self.fullupper/self.refinement)
      if alltrue(u >= l):
        self.parents.append(block)
        block.addblockaschild(self)
        self.ichild.append(len(block.children))
    # --- Now pass rest of lists to children
    for child in self.children:
      child.findallparents(blocklists[1:])

  def generateblocklevellists(self,blocklists=None):
    if blocklists is None:
      # --- This will only happen at the top level.
      # --- Create a list of empty lists. Each empty list will get the blocks
      # --- at the appropriate levels appended to it. Note that 100 is
      # --- assumed to be a large enough number - there almost certainly
      # --- will never be 100 levels of refinement. The topmost list
      # --- will be empty since the top level has not parents.
      blocklists = []
      for i in range(100):
        blocklists.append([])
    # --- Add this instance to the top level of the list and pass the rest
    # --- of it to the children
    if self not in blocklists[1]: blocklists[1].append(self)
    for child in self.children:
      child.generateblocklevellists(blocklists[1:])
    return blocklists

  def finalize(self):
    # --- This should only be called at the top level.
    blocklists = self.generateblocklevellists()
    self.findallparents(blocklists)
    self.findoverlappingsiblings()


  def loadrho(self,lzero=true):
    """
Loads the charge density from the particles. This should only be called for
the top level grid.
    """
    if lzero: self.zerorho()
    for i,n,q,w in zip(top.ins-1,top.nps,top.sq,top.sw):
      if n > 0:
        self.setrho(top.xp[i:i+n],top.yp[i:i+n],top.zp[i:i+n],top.uzp[i:i+n],
                    q,w)
    self.addmyrhotosiblings()
    self.getrhofromsiblings()
    self.gatherrhofromchildren()
    #self.setrhoboundaries()

  def gete(self,x,y,z,ex,ey,ez):
    if len(x) == 0: return
    if len(self.children) > 0:
      # --- Find out whether the particles are in the local domain or one of
      # --- the children's.
      ichild = zeros(len(x),'d')
      getgridngp3d(len(x),x,y,z-top.zgridprv,ichild,
                   self.nx,self.ny,self.nz,self.idomains,
                   self.xmmin,self.xmmax,self.ymmin,self.ymmax,
                   self.zmmin,self.zmmax,self.l2symtry,self.l4symtry)
      # --- Zero out places where idomains < 0
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
          MultiGrid.gete(self,xc,yc,zc,tex,tey,tez)
        else:
          self.children[i-1].gete(xc,yc,zc,tex,tey,tez)
        # --- Put the E field into the passed in arrays
        put(ex,ii,tex)
        put(ey,ii,tey)
        put(ez,ii,tez)

# --- This is the old way with sorting - there's a bug somewhere!!
#     # --- Sort the list of which domain the particles are in.
#     ichildsorter = argsort(ichild)
#     ichildsorted = take(ichild,ichildsorter)
#     nperchild = self.getnperchild(ichildsorted)
#     # --- Now use the same sorting on the particle quantities.
#     x = take(x,ichildsorter)
#     y = take(y,ichildsorter)
#     z = take(z,ichildsorter)
#     # --- Create temporary arrays to hold the E field
#     tex,tey,tez = zeros(len(x),'d'),zeros(len(x),'d'),zeros(len(x),'d')
#     # --- For each child, pass to it the particles in it's domain.
#     i = nperchild[0]
#     for child,n in zip(self.children,nperchild[1:]):
#       if n > 0:
#         child.gete(x[i:i+n],y[i:i+n],z[i:i+n],
#                    tex[i:i+n],tey[i:i+n],tez[i:i+n])
#       i = i + n
#     # --- Get e-field from this domain
#     n = nperchild[0]
#     MultiGrid.gete(self,x[:n],y[:n],z[:n],tex[:n],tey[:n],tez[:n])
#     # --- Put the E field into the passed in arrays
#     put(ex,ichildsorter,tex)
#     put(ey,ichildsorter,tey)
#     put(ez,ichildsorter,tez)

    else:
      # --- Get e-field from this domain
      MultiGrid.gete(self,x,y,z,ex,ey,ez)

  def setphi(self,x,y,z,phi):
    if len(x) == 0: return
    if len(self.children) > 0:
      # --- Find out whether the particles are in the local domain or one of
      # --- the children's.
      ichild = zeros(len(x),'d')
      getgridngp3d(len(x),x,y,z,ichild,
                   self.nx,self.ny,self.nz,self.idomains,
                   self.xmmin,self.xmmax,self.ymmin,self.ymmax,
                   self.zmmin,self.zmmax,self.l2symtry,self.l4symtry)
      # --- Zero out places where idomains < 0
      ichild = where(ichild < 0.,0.,ichild)
      # --- Sort the list of which domain the particles are in.
      ichildsorter = argsort(ichild)
      ichildsorted = take(ichild,ichildsorter)
      nperchild = self.getnperchild(ichildsorted)
      # --- Now use the same sorting on the particle quantities.
      x = take(x,ichildsorter)
      y = take(y,ichildsorter)
      z = take(z,ichildsorter)
      # --- Create temporary arrays to hold the E field
      tphi = zeros(len(x),'d')
      # --- For each child, pass to it the particles in it's domain.
      i = nperchild[0]
      for child,n in zip(self.children,nperchild[1:]):
        if n > 0:
          child.setphi(x[i:i+n],y[i:i+n],z[i:i+n],tphi[i:i+n])
        i = i + n
      # --- Get e-field from this domain
      n = nperchild[0]
      MultiGrid.setphi(self,x[:n],y[:n],z[:n],tphi[:n])
      # --- Put the phi into the passed in arrays
      put(phi,ichildsorter,tphi)
    else:
      # --- Get e-field from this domain
      MultiGrid.setphi(self,x,y,z,phi)

  def fetche(self):
    """
Fetches the E field. This should only be called for
the top level grid.
    """
#   ex,ey,ez = zeros(top.npmax,'d'),zeros(top.npmax,'d'),zeros(top.npmax,'d')
#   for i,n in zip(top.ins-1,top.nps):
#     self.gete(top.xp[i:i+n],top.yp[i:i+n],top.zp[i:i+n],
#               ex[i:i+n],ey[i:i+n],ez[i:i+n])
#   return ex,ey,ez
    self.gete(w3d.xefield,w3d.yefield,w3d.zefield,
              w3d.exefield,w3d.eyefield,w3d.ezefield)

  def solve(self,iwhich=0):
    # --- Wait until all of the parents have called here until actually
    # --- doing the solve. This ensures that the phiboundaries obtained
    # --- will be up to date.
    self.ncallsfromparents += 1
    if self.ncallsfromparents < len(self.parents): return
    self.ncallsfromparents = 0

    self.setgridboundaries('phi')
    MultiGrid.solve(self,iwhich)
    for child in self.children:
      child.solve(iwhich)

  def installconductors(self,conductor):
    self.conductors.append(conductor)
    MultiGrid.installconductors(self,conductor)
    for child in self.children:
      child.installconductors(conductor)

  #############################################################################
  # --- The following are used for plotting.
  def genericpf(self,kw,idim,pffunc,ip=None):
    """
Generic plotting routine. This plots only the local domain. Domains of the
children are also skipped, but the same call is made for them so they will
be plotted.
    """
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
    # --- If that plane is not in the domain, then don't do anything
    if ip < self.fulllower[idim] or self.fullupper[idim] < ip: return
    # --- Create the ireg array, which will be set to zero in the domain
    # --- of the children.
    ss = list(shape(self.rho))
    del ss[idim]
    ireg = zeros(ss)
    if self.idomains is not None:
      if idim==0:
        ireg[1:,1:]=equal(self.idomains[ip-self.fulllower[0],:-1,:-1],0)
      if idim==1:
        ireg[1:,1:]=equal(self.idomains[:-1,ip-self.fulllower[1],:-1],0)
      if idim==2:
        ireg[1:,1:]=equal(self.idomains[:-1,:-1,ip-self.fulllower[2]],0)
    else:
      ireg[1:,1:] = 1
    if idim != 2: ireg = transpose(ireg)
    kw['ireg'] = ireg
    MultiGrid.genericpf(self,kw,pffunc)
    kw['titles'] = 0
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

  def drawbox(self,ip=None,idim=2,withguards=1,color='fg'):
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
    plg(xx,yy,color=color)
    for child in self.children:
      child.drawbox(ip=ip*child.refinement,idim=idim,withguards=withguards,
                    color=color)

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

