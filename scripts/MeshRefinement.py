"""Implements adaptive mesh refinement in 3d
"""
from warp import *
from multigrid import MultiGrid

class BlockOverlap:
  def __init__(self,block,fulllower,fullupper,notmine):
    self.block = block
    self.fulllower = fulllower
    self.fullupper = fullupper
    self.notmine = notmine

class MRBlock(MultiGrid):
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

    if parent is None: self.parents = []
    else:              self.parents = [parent]
    self.overlaps = []
    self.siblings = []
    self.ichild = [ichild]
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
      self.bound0 = dirichlet
      self.boundnz = dirichlet
      self.boundxy = dirichlet

      if lower is None and upper is None:
        # --- The grid mins and maxs are input
        self.mins = mins
        self.maxs = max
        self.dims = (self.maxs - self.mins)/self.deltas
        self.lower = (self.mins - parent.mins)/self.deltas + parent.lower
        self.upper = (self.maxs - parent.mins)/self.deltas + parent.lower

      else:
        # --- The grid lower and upper bounds are input. The bounds are
        # --- relative to the parent. They are converted to be relative
        # --- to this instance.
        self.lower = refinement*(parent.lower + array(lower))
        self.upper = refinement*(parent.lower + array(upper))
        self.dims = self.upper - self.lower
        self.mins = parent.mins+(self.lower-refinement*parent.lower)*self.deltas
        self.maxs = parent.mins+(self.upper-refinement*parent.lower)*self.deltas

      # --- Now, extend the domain by the given number of guard cells. Checks
      # --- are made so that the domain doesn't extend beyon the original grid.
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

    # --- Create some work data arrays needed for rho deposition
    w0 = array([0.25,0.5,0.25])
    self.w = outerproduct(w0,outerproduct(w0,w0))
    self.w.shape = (3,3,3)
    self.ncallsfromparents = 0

  def addchild(self,lower,upper,refinement=2):
    child = MRBlock(parent=self,lower=lower,upper=upper,
                    refinement=refinement,ichild=len(self.children)+1,
                    nguard=self.nguard,dimsmax=self.dimsmax*refinement)
    self.addblockaschild(child)

  def addblockaschild(self,block):
    self.children.append(block)
    ix1,iy1,iz1 = maximum(0,block.fulllower/block.refinement - self.fulllower)
    ix2,iy2,iz2 =           block.fullupper/block.refinement - self.fulllower
    if self.idomains is None:
      self.idomains = fzeros((self.nx,self.ny,self.nz),'d')
    self.idomains[ix1:ix2,iy1:iy2,iz1:iz2] = len(self.children)

  def setphiboundaries(self):
    """
Sets phi on the boundaries, using the values from the parent grid
    """
    for parent in self.parents:
      r = self.refinement
      i1 = maximum(0,self.fulllower/self.refinement - parent.fulllower)
      i2 = minimum(parent.dims,
                   self.fullupper/self.refinement - parent.fulllower)
      i3 = (parent.fulllower + i1)*self.refinement - self.fulllower
      i4 = (parent.fulllower + i2)*self.refinement - self.fulllower
      ix1,iy1,iz1 = i1
      ix2,iy2,iz2 = i2
      ix3,iy3,iz3 = i3
      ix4,iy4,iz4 = i4
      # --- First get points directly copied from parent.
      # --- As a convenience, get references to the phi arrays without
      # --- the z guard planes.
      sphi = self.phi[:,:,1:-1]
      pphi = parent.phi[:,:,1:-1]
      # --- Note that of the child extends beyond the parent's domain,
      # --- then certain boundary planes are completely skipped.
      if self.fulllower[0]/self.refinement >= parent.fulllower[0]:
        sphi[ix3,iy3:iy4+1:r,iz3:iz4+1:r] = pphi[ix1,iy1:iy2+1,iz1:iz2+1]
      if self.fullupper[0]/self.refinement <= parent.fullupper[0]:
        sphi[ix4,iy3:iy4+1:r,iz3:iz4+1:r] = pphi[ix2,iy1:iy2+1,iz1:iz2+1]
      if self.fulllower[1]/self.refinement >= parent.fulllower[1]:
        sphi[ix3:ix4+1:r,iy3,iz3:iz4+1:r] = pphi[ix1:ix2+1,iy1,iz1:iz2+1]
      if self.fullupper[1]/self.refinement <= parent.fullupper[1]:
        sphi[ix3:ix4+1:r,iy4,iz3:iz4+1:r] = pphi[ix1:ix2+1,iy2,iz1:iz2+1]
      if self.fulllower[2]/self.refinement >= parent.fulllower[2]:
        sphi[ix3:ix4+1:r,iy3:iy4+1:r,iz3] = pphi[ix1:ix2+1,iy1:iy2+1,iz1]
      if self.fullupper[2]/self.refinement <= parent.fullupper[2]:
        sphi[ix3:ix4+1:r,iy3:iy4+1:r,iz4] = pphi[ix1:ix2+1,iy1:iy2+1,iz2]

    # --- Do remaining points where interpolation is needed.
    # --- This can be done, now that all of the contributions from the parents
    # --- is done.
    self.setslice(self.phi[ 0,:,1:-1])
    self.setslice(self.phi[-1,:,1:-1])
    self.setslice(self.phi[:, 0,1:-1])
    self.setslice(self.phi[:,-1,1:-1])
    self.setslice(self.phi[:,:, 1])
    self.setslice(self.phi[:,:,-2])

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

# def copyrhofromparent(self):
#   # --- Copies rho from parent, only used for debugging
#   if self.parent is None: return
#   r = self.refinement
#   ix1 = self.lower[0]
#   ix2 = self.upper[0]
#   iy1 = self.lower[1]
#   iy2 = self.upper[1]
#   iz1 = self.lower[2]
#   iz2 = self.upper[2]
#   self.rho[::r,::r,::r] = self.parent.rho[ix1:ix2+1,iy1:iy2+1,iz1:iz2+1]
#   for ix in arange(0,self.dims[0]+1): self.setslice(self.rho[ix,:,:])
#   for iy in arange(0,self.dims[1]+1): self.setslice(self.rho[:,iy,:])
#   for iz in arange(0,self.dims[2]+1): self.setslice(self.rho[:,:,iz])
#   self.rho[1::r,1::r,1::r] = (1./6.)*(self.rho[ :-1:r, :-1:r, :-1:r] +
#                                       self.rho[2:  :r, :-1:r, :-1:r] +
#                                       self.rho[ :-1:r,2:  :r, :-1:r] +
#                                       self.rho[2:  :r,2:  :r, :-1:r] +
#                                       self.rho[ :-1:r, :-1:r,2:  :r] +
#                                       self.rho[2:  :r, :-1:r,2:  :r] +
#                                       self.rho[ :-1:r,2:  :r,2:  :r] +
#                                       self.rho[2:  :r,2:  :r,2:  :r])

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
    result[-1] = len(ichildsorted) - sum(result)
    return result

  def setrho(self,x,y,z,uz,q,w):
    if len(self.children) > 0:
      # --- Find out whether the particles are in the local domain or one of
      # --- the children's.
      ichild = zeros(len(x),'d')
      getgridngp3d(len(x),x,y,z,ichild,
                   self.nx-1,self.ny-1,self.nz-1,self.idomains,
                   self.xmmin,self.xmmax-self.dx,self.ymmin,self.ymmax-self.dy,
                   self.zmmin,self.zmmax-self.dz,self.l2symtry,self.l4symtry)
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
    # --- This needs work - a child can pass the same rho to multiple parents
    # --- who then all give it to a single grandparent, multiplying the
    # --- rho deposited there.
    for child in self.children:
      child.gatherrhofromchildren()
      # --- Get coordinates of child relative to this domain
      i1 = maximum(0,child.fulllower/child.refinement - self.fulllower)
      i2 = minimum(self.dims,
                   child.fullupper/child.refinement - self.fulllower) + 1
      r = child.refinement
      dims = [-1,0,+1]
      # --- Zero out the rho internal to the child's domain.
      # --- Note that the boundary will have charge deposited from particles.
      # --- This is not really needed since zerorho will always be called.
      #self.rho[ix1+1:ix2-1,iy1+1:iy2-1,iz1+1:iz2-1] = 0.
      # --- Loop over the three dimensions
      for k in dims:
        for j in dims:
          for i in dims:
            # --- Make sure points outside the child's domain are skipped.
            pxm,pxp,pym,pyp,pzm,pzp = (i<0),-(i>0),(j<0),-(j>0),(k<0),-(k>0)
            cxm,cxp,cym,cyp,czm,czp = r*array([pxm,pxp,pym,pyp,pzm,pzp])+array([i,i,j,j,k,k])
            cxp,cyp,czp = array([cxp,cyp,czp]) + (i2 - i1)*child.refinement
            ix1,iy1,iz1 = i1
            ix2,iy2,iz2 = i2
            p = self.rho[ix1+pxm:ix2+pxp,iy1+pym:iy2+pyp,iz1+pzm:iz2+pzp]
            c = child.rho[cxm:cxp:r,cym:cyp:r,czm:czp:r]
            # --- Do the addition in place.
            add(p,c*self.w[i+1,j+1,k+1],p)

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
      add(where(sibling.notmine,myrho,0.),nrho,nrho)
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
      myrho[:,:] = where(sibling.notmine,nrho,myrho)
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
    # --- Check how much of the overlap region this instance does not own.
    # --- First, assume that is doesn't own any.
    notmine = ones(u-l+1)
    pl = l/self.refinement
    pu = u/self.refinement
    r = self.refinement
    # --- Check the idomains for each parent to find out the regions
    # --- owned by this instance.
    for parent,ichild in zip(self.parents,self.ichild):
      pld = maximum(pl,parent.fulllower)
      pud = minimum(pu,parent.fullupper)
      ix1,iy1,iz1 = pld*self.refinement - l
      ix2,iy2,iz2 = pud*self.refinement - l
      ix3,iy3,iz3 = pld - parent.fulllower
      ix4,iy4,iz4 = pud - parent.fulllower
      notmineslice = notmine[ix1:ix2,iy1:iy2,iz1:iz2]
      idomainsslice = parent.idomains[ix3:ix4,iy3:iy4,iz3:iz4]
      # --- For the cells owned by this instance, set notmine to zero.
      iii = where(idomainsslice==ichild,0,notmineslice[:-1:r,:-1:r,:-1:r])
      # --- That array is relative to the parents mesh size. Copy to each
      # --- of the eight nodes in the lower corner of the parent cell.
      # --- The other nodes are part of different cells in the parent.
      notmineslice[ :-1:r, :-1:r, :-1:r] = iii
      notmineslice[1:  :r, :-1:r, :-1:r] = iii
      notmineslice[ :-1:r,1:  :r, :-1:r] = iii
      notmineslice[1:  :r,1:  :r, :-1:r] = iii
      notmineslice[ :-1:r, :-1:r,1:  :r] = iii
      notmineslice[1:  :r, :-1:r,1:  :r] = iii
      notmineslice[ :-1:r,1:  :r,1:  :r] = iii
      notmineslice[1:  :r,1:  :r,1:  :r] = iii
    # --- Only keep track of overlaps where this instance does not own the
    # --- whole overlap region. If this instance owns the who overlap
    # --- region, then nothing special needs to be done so there is no
    # --- reason to keep track of the overlap.
    if maxnd(notmine) == 1:
      self.overlaps.append(BlockOverlap(other,l,u,notmine))

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
      self.setrho(top.xp[i:i+n],top.yp[i:i+n],top.zp[i:i+n],top.uzp[i:i+n],q,w)
    self.gatherrhofromchildren()
    self.addmyrhotosiblings()
    self.getrhofromsiblings()

  def solve(self,iwhich=0):
    # --- Wait until all of the parents have called here until actually
    # --- doing the solve. This ensures that the phiboundaries obtained
    # --- will be up to date.
    self.ncallsfromparents += 1
    if self.ncallsfromparents < len(self.parents): return
    self.ncallsfromparents = 0

    self.setphiboundaries()
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
    # --- If that plane is not in the domain, then don't do anything
    if ip < self.fulllower[idim] or self.fullupper[idim] < ip: return
    # --- Create the ireg array, which will be set to zero in the domain
    # --- of the children.
    ss = list(shape(self.rho))
    del ss[idim]
    ireg = zeros(ss)
    if self.idomains is not None:
      if idim==0: ireg[1:,1:]=equal(self.idomains[ip-self.fulllower[0],:,:],0)
      if idim==1: ireg[1:,1:]=equal(self.idomains[:,ip-self.fulllower[1],:],0)
      if idim==2: ireg[1:,1:]=equal(self.idomains[:,:,ip-self.fulllower[2]],0)
    else:
      ireg[1:,1:] = 1
    if idim != 2: ireg = transpose(ireg)
    kw['ireg'] = ireg
    MultiGrid.genericpf(self,kw,pffunc)
    for child in self.children:
      child.genericpf(kw,idim,pffunc,ip)

  def pfxy(self,**kw): self.genericpf(kw,2,pfxy)
  def pfzx(self,**kw): self.genericpf(kw,1,pfzx)
  def pfzy(self,**kw): self.genericpf(kw,0,pfzy)
  def pfxyg(self,**kw): self.genericpf(kw,2,pfxyg)
  def pfzxg(self,**kw): self.genericpf(kw,1,pfzxg)
  def pfzyg(self,**kw): self.genericpf(kw,0,pfzyg)


  def plphiz(self,ix=None,iy=None):
    plg(self.phi[ix-self.fulllower[0],iy-self.fulllower[1],1:-1],
        self.mins[2]+arange(self.dims[2]+1)*self.deltas[2])
    for child in self.children:
      child.plphiz(ix*child.refinement,iy*child.refinement)

  def plrhoz(self,ix=None,iy=None):
    plg(self.rho[ix-self.fulllower[0],iy-self.fulllower[1],:],
        self.mins[2]+arange(self.dims[2]+1)*self.deltas[2])
    for child in self.children:
      child.plrhoz(ix*child.refinement,iy*child.refinement)

