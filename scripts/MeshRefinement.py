"""Implements adaptive mesh refinement in 3d
"""
from warp import *
from multigrid import MultiGrid

class MRBlock(MultiGrid):
  """
 - children: list of tuples, each containing three elements,
             (lower,upper,refinement).
  """
  def __init__(self,parent=None,refinement=2,lower=None,upper=None,
                    dims=None,mins=None,maxs=None,
                    children=None,**kw):

    self.parent = parent
    self.refinement = refinement
    self.conductors = []

    if self.parent is None:
      # --- All grid sizing parameters should be input, otherwise they
      # --- will be taken from w3d.
      self.dims = dims
      self.mins = mins
      self.maxs = maxs
      
    elif lower is None and upper is None:
      # --- The grid mins and maxs are input
      self.mins = mins
      self.maxs = max
      self.dims = self.parent.dims/self.refinement
      self.lower = (self.mins - self.parent.mins)/self.parent.deltas
      self.upper = (self.maxs - self.parent.mins)/self.parent.deltas

    else:
      # --- The grid lower and upper bounds are input. The bounds are
      # --- relative to the parent.
      self.lower = array(lower)
      self.upper = array(upper)
      self.dims = (self.upper - self.lower)*self.refinement
      self.mins = self.parent.mins + self.lower*self.parent.deltas
      self.maxs = self.parent.mins + self.upper*self.parent.deltas

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

    self.dims = array([self.nx,self.ny,self.nz])
    self.deltas = array([self.dx,self.dy,self.dz])
    self.mins = array([self.xmmin,self.ymmin,self.zmmin])
    self.maxs = array([self.xmmax,self.ymmax,self.zmmax])

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
    w0 = array([0.5,1.,0.5])
    w = outerproduct(w0,outerproduct(w0,w0))
    w.shape = (3,3,3)
    self.w = w/sum(sum(sum(w)))

  def addchild(self,lower,upper,refinement=2):
    child = MRBlock(parent=self,lower=lower,upper=upper,
                    refinement=refinement)
    self.children.append(child)
    ix1 = lower[0]
    iy1 = lower[1]
    iz1 = lower[2]
    ix2 = upper[0]
    iy2 = upper[1]
    iz2 = upper[2]
    if self.idomains is None:
      self.idomains = fzeros((self.nx,self.ny,self.nz),'d')
    self.idomains[ix1:ix2,iy1:iy2,iz1:iz2] = len(self.children)

  def setphiboundaries(self):
    """
Sets phi on the boundaries, using the values from the parent grid
    """
    if self.parent is None: return
    r = self.refinement
    ix1 = self.lower[0]
    ix2 = self.upper[0]
    iy1 = self.lower[1]
    iy2 = self.upper[1]
    iz1 = self.lower[2]+1
    iz2 = self.upper[2]+1
    # --- First get points directly copied from parent
    self.phi[ 0,::r,1:-1:r] = self.parent.phi[ix1,iy1:iy2+1,iz1:iz2+1]
    self.phi[-1,::r,1:-1:r] = self.parent.phi[ix2,iy1:iy2+1,iz1:iz2+1]
    self.phi[::r, 0,1:-1:r] = self.parent.phi[ix1:ix2+1,iy1,iz1:iz2+1]
    self.phi[::r,-1,1:-1:r] = self.parent.phi[ix1:ix2+1,iy2,iz1:iz2+1]
    self.phi[::r,::r, 1] =    self.parent.phi[ix1:ix2+1,iy1:iy2+1,iz1]
    self.phi[::r,::r,-2] =    self.parent.phi[ix1:ix2+1,iy1:iy2+1,iz2]

    # --- Do remaining points where interpolation is needed.
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

  def copyrhofromparent(self):
    # --- Copies rho from parent, only used for debugging
    if self.parent is None: return
    r = self.refinement
    ix1 = self.lower[0]
    ix2 = self.upper[0]
    iy1 = self.lower[1]
    iy2 = self.upper[1]
    iz1 = self.lower[2]
    iz2 = self.upper[2]
    self.rho[::r,::r,::r] = self.parent.rho[ix1:ix2+1,iy1:iy2+1,iz1:iz2+1]
    for ix in arange(0,self.dims[0]+1): self.setslice(self.rho[ix,:,:])
    for iy in arange(0,self.dims[1]+1): self.setslice(self.rho[:,iy,:])
    for iz in arange(0,self.dims[2]+1): self.setslice(self.rho[:,:,iz])
    self.rho[1::r,1::r,1::r] = (1./6.)*(self.rho[ :-1:r, :-1:r, :-1:r] +
                                        self.rho[2:  :r, :-1:r, :-1:r] +
                                        self.rho[ :-1:r,2:  :r, :-1:r] +
                                        self.rho[2:  :r,2:  :r, :-1:r] +
                                        self.rho[ :-1:r, :-1:r,2:  :r] +
                                        self.rho[2:  :r, :-1:r,2:  :r] +
                                        self.rho[ :-1:r,2:  :r,2:  :r] +
                                        self.rho[2:  :r,2:  :r,2:  :r])

  def getnperchild(self,ichildsorted):
    """
Given a sorted list of child numbers, returns how many of each number
there are. This would probably be better done in compiled code. This should be
efficient timewise, but uses two extra full-size arrays.
    """
    ic1 = ichildsorted[1:] - ichildsorted[:-1]
    nn = compress(ic1 > 0,iota(len(ichildsorted-1)))
    result = zeros(1+len(self.children))
    for n in nn:
      result[nint(ichildsorted[n-1])] = n
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
    for child in self.children:
      child.gatherrhofromchildren()
      # --- Get coordinates of child relative to this domain
      ix1,iy1,iz1 = child.lower
      ix2,iy2,iz2 = array(child.upper) + 1
      r = child.refinement
      dims = [-1,0,+1]
      # --- Zero out the rho internal to the child's domain.
      # --- Note that the boundary will have charge deposited from particles.
      self.rho[ix1+1:ix2-1,iy1+1:iy2-1,iz1+1:iz2-1] = 0.
      # --- Loop over the three dimensions
      for k in dims:
        for j in dims:
          for i in dims:
            # --- Make sure points outside the child's domain are skipped.
            pxm,pxp,pym,pyp,pzm,pzp = (i<0),-(i>0),(j<0),-(j>0),(k<0),-(k>0)
            cxm,cxp,cym,cyp,czm,czp = r*array([pxm,pxp,pym,pyp,pzm,pzp])+array([i,i,j,j,k,k])
            cxp,cyp,czp = array([cxp,cyp,czp]) + child.dims + 1
            p = self.rho[ix1+pxm:ix2+pxp,iy1+pym:iy2+pyp,iz1+pzm:iz2+pzp]
            c = child.rho[cxm:cxp:r,cym:cyp:r,czm:czp:r]
            # --- Do the addition in place.
            add(p,c*self.w[i+1,j+1,k+1],p)

  def loadrho(self,lzero=true):
    """
Loads the charge density from the particles. This should only be called for
the top level grid.
    """
    if lzero: self.zerorho()
    for i,n,q,w in zip(top.ins-1,top.nps,top.sq,top.sw):
      self.setrho(top.xp[i:i+n],top.yp[i:i+n],top.zp[i:i+n],top.uzp[i:i+n],q,w)
    self.gatherrhofromchildren()

  def solve(self,iwhich=0):
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
  def clearparentsireg(self,ireg,idim):
    # --- Clears out the local mesh domain from the parent's ireg array.
    # --- Plots made by the parent then don't show this instances domain.
    i1 = []
    i2 = []
    for d in range(len(self.dims)):
      if d == idim: continue
      i1.append(self.lower[d]+1)
      i2.append(self.upper[d])
    ireg[i1[0]:i2[0],i1[1]:i2[1]] = 0
    
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
      if ip is None: ip = (-self.mins[idim]/self.deltas[idim])
    else:
      ip = (ip - self.lower[idim])*self.refinement
    # --- If that plane is not in the domain, then don't do anything
    if ip < 0 or self.dims[idim] < ip: return
    # --- Create the ireg array, which will be set to zero in the domain
    # --- of the children.
    ss = list(shape(self.rho))
    del ss[idim]
    ireg = ones(ss)
    ireg[0,:] = 0
    ireg[:,0] = 0
    for child in self.children: child.clearparentsireg(ireg,idim)
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


