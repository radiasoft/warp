"""Implements adaptive mesh refinement in 3d for the B field solver
"""
from warp import *
from magnetostaticMG import MagnetostaticMG
from pyOpenDX import Visualizable,DXCollection,viewboundingbox
import MA
import __main__
import time
import operator
try:
  import psyco
except ImportError:
  pass

#########################################################################
# Note that MRBlockB is psyco.bind at the end of the file
class MRBlockB(MagnetostaticMG,Visualizable):
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
                    dims=None,mins=None,maxs=None,
                    nguard=1,
                    children=None,lreducedpickle=1,**kw):

    # --- Check the dimensionality.
    if w3d.solvergeom in [w3d.XYZgeom]:
      self.nd = 3
    elif w3d.solvergeom in [w3d.RZgeom]:
      self.nd = 2

    # --- Pass the input value of lreducedpickle into the Multigrid class
    kw['lreducedpickle'] = lreducedpickle

    if parent is None:
      # --- No parents, so just create empty lists
      self.parents = []
      self.root = self
      # --- It is assumed that the root block will be the first one created.
      # --- So clear out the global block list and count.
      self.totalnumberofblocks = 0
      self.listofblocks = []
      self.finalized = 0
    else:
      # --- Save the parent and the index number. These are saved in lists
      # --- since a block can have multiple parents.
      self.parents = [parent.blocknumber]
      self.root = parent.root

    # --- Get the current global block number and increment the counter.
    # --- Also, add self to the global list of blocks.
    self.blocknumber = self.root.totalnumberofblocks
    self.root.totalnumberofblocks += 1
    self.root.listofblocks.append(self)

    # --- Note that a dictionary is used for the overlaps so that lookups
    # --- are faster, also, the values in the dictionary are list containing
    # --- the domain of the overlap. The blocks with lower and higher block
    # --- number are treated differently, so create separate lists for each.
    # --- Two separate lists is likely only a small optimization.
    self.overlapslower = {}
    self.overlapshigher = {}
    self.nguard = nguard

    # --- Check the input for errors.
    errorstring = "Input argument %s must be of length %d"
    if refinement is not None and operator.isSequenceType(refinement):
      assert len(refinement)==self.nd,errorstring%('refinement',self.nd)
    if lower is not None:
      assert len(lower)==self.nd,errorstring%('lower',self.nd)
    if upper is not None:
      assert len(upper)==self.nd,errorstring%('upper',self.nd)
    if dims is not None:
      assert len(dims)==self.nd,errorstring%('dims',self.nd)
    if mins is not None:
      assert len(mins)==self.nd,errorstring%('mins',self.nd)
    if maxs is not None:
      assert len(maxs)==self.nd,errorstring%('maxs',self.nd)

    # --- Check the dimensionality. If less than 3D, set the appropriate
    # --- dimension to zero. (The code only works with XYZ and RZ now.)
    if w3d.solvergeom == w3d.RZgeom:
      if refinement is not None:
        if operator.isSequenceType(refinement):
          refinement = [refinement[0],1,refinement[1]]
        else:
          refinement = [refinement,1,refinement]
      if lower is not None:
        lower = [lower[0],0,lower[1]]
      if upper is not None:
        upper = [upper[0],0,upper[1]]
      if dims is not None:
        dims = [dims[0],0,dims[1]]
      if mins is not None:
        mins = [mins[0],0.,mins[1]]
      if maxs is not None:
        maxs = [maxs[0],0.,maxs[1]]

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
        # --- The lower and upper are calculated to be an integer number of
        # --- parent grid cells times the refinement factor. The lower is
        # --- rounded down and upper rounded up to ensure that the patch
        # --- includes the entire extent specified by mins and maxs.
        self.mins = array(mins)
        self.maxs = array(maxs)
        self.lower = (nint(floor((self.mins - self.root.mins)/parent.deltas))*
                     self.refinement)
        self.upper = (nint(ceil((self.maxs - self.root.mins)/parent.deltas))*
                     self.refinement)

      else:
        # --- The grid lower and upper bounds are input. The bounds are
        # --- relative to the root grid, but scaled by the total refinement.
        self.lower = array(lower)
        self.upper = array(upper)

      # --- Now, extend the domain by the given number of guard cells. Checks
      # --- are made so that the domain doesn't extend beyond the original grid.
      self.fulllower = maximum(0,self.lower - nguard*self.refinement)
      self.fullupper = minimum(self.rootdims,
                               self.upper + nguard*self.refinement)

      # --- Get the number of grid points along each dimension
      self.dims = self.fullupper - self.fulllower

      # --- Make sure that the number of grid points is even.
      # --- If it is odd, then enough cells are added to extend to the next
      # --- grid cell of the parent. It is then cutoff at zero.
      self.fulllower = where(self.dims%2==1,self.fulllower-self.refinement,
                                            self.fulllower)
      self.fulllower = maximum(0,self.fulllower)
      self.dims = self.fullupper - self.fulllower

      # --- If it is still odd (which means that the cells added above
      # --- where cutoff at zero) then add some at the top.
      self.fullupper = where(self.dims%2==1,self.fullupper+self.refinement,
                                            self.fullupper)
      self.fullupper = minimum(self.rootdims,self.fullupper)
      self.dims = self.fullupper - self.fulllower

      # --- If it is still odd, then there is some serious problem. The number
      # --- in the base grid may be odd.
      assert alltrue(self.dims%2 == 0),\
             """The number of grid cells in one of the dimensions is odd - they
             all must be even. Check that the number of cells in the base grid
             is even."""

      # --- Now calculate the extent of the grid
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

      # --- Create some temporaries for optimization
      self.fullloweroverrefinement = self.fulllower/self.refinement
      self.fullupperoverrefinement = self.fullupper/self.refinement

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
    MagnetostaticMG.__init__(self,**kw)

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

    # --- childdomains is the node centered grid which keeps track of which
    # --- cells are owned by which children. If there are no children,
    # --- then it is not needed.
    self.childdomains = None

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
    # --- Make sure the the lreducedpickle option gets propagated to all
    # --- of the blocks.
    for child in self.children:
      child.lreducedpickle = self.lreducedpickle
    dict = MagnetostaticMG.__getstate__(self)
    if self.lreducedpickle:
      # --- Remove the big objects from the dictionary. This can be
      # --- regenerated upon the restore.
      dict['childdomains'] = None
    # --- Flag whether this is the registered solver so it know whether
    # --- to reregister itself upon the restore.
    if self is getregisteredsolver():
      dict['iamtheregisteredsolver'] = 1
    else:
      dict['iamtheregisteredsolver'] = 0
    return dict

  def __setstate__(self,dict):
    MagnetostaticMG.__setstate__(self,dict)
   #self.makefortranordered('phi')
   #self.makefortranordered('rho')
   #self.makefortranordered('selfe')
    if self.iamtheregisteredsolver:
      del self.iamtheregisteredsolver
      registersolver(self)
    if self == self.root and self.lreducedpickle:
      # --- It is assumed that at this point, all of the children have been
      # --- restored.
      # --- Regenerate childdomains
      self.initializechilddomains()
      # --- If rho and phi weren't saved, make sure that they are setup.
      # --- Though, this may not always be the right thing to do.
      # --- These can only be done at the end of the restart since only then
      # --- is it gauranteed that the particles are read in.
      installafterrestart(self.loadj)
      installafterrestart(self.solve)

  def makefortranordered(self,vname):
    a = getattr(self,vname)
    if type(a) is ArrayType:
      setattr(self,vname,fzeros(shape(a),a.typecode()))
      getattr(self,vname)[...] = a

  def addchild(self,lower=None,upper=None,mins=None,maxs=None,
                    refinement=2,nslaves=1,my_index=0):
    """
Add a mesh refined block to this block.
  -lower,upper,mins,maxs,refinement: All have same meanings as for the
                                     constructor.
  -nslaves=1: defaults to one so it is not parallelized
    """
    child = MRBlockB(parent=self,lower=lower,upper=upper,mins=mins,maxs=maxs,
                     refinement=refinement,nguard=self.nguard,
                     nslaves=nslaves,my_index=my_index)
    self.children.append(child)

  def addblockaschild(self,block):
    """
Given a block instance, installs it as a child.
    """
    self.children.append(block)

  def resetroot(self):
    # --- No parents, so just create empty lists
    self.parents = []
    self.root = self
    self.totalnumberofblocks = 1
    self.listofblocks = [self]
    self.childdomains = None
    self.children = []

  #--------------------------------------------------------------------------
  # --- The next several methods handle conductors
  #--------------------------------------------------------------------------

  def installconductor(self,conductor,dfill=top.largepos):
    if not self.isfirstcall(): return
    MagnetostaticMG.installconductor(self,conductor,dfill=dfill)
    for child in self.children:
      child.installconductor(conductor,dfill=dfill)

  def clearconductors(self):
    if not self.isfirstcall(): return
    MagnetostaticMG.clearconductors(self)
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

  def setconductorvoltage(self,voltage,condid=0,discrete=false,
                          setvinject=false):
    'Recursively calls setconductorvoltage for base and all children'
    if not self.isfirstcall(): return
    setconductorvoltage(voltage,condid,discrete,setvinject,
                        conductors=self.conductors)
    for child in self.children:
      child.setconductorvoltage(voltage,condid,discrete)

  #--------------------------------------------------------------------------
  # --- The next several methods handle initialization that is done after
  # --- all blocks have been added.
  #--------------------------------------------------------------------------

  def finalize(self):
    # --- This should only be called at the top level.
    if self != self.root or self.finalized: return
    blocklists = self.generateblocklevellists()
    self.clearparentsandchildren()
    self.findallchildren(blocklists)
    self.initializechilddomains()
    self.findoverlappingsiblings(blocklists[1:])
    self.setmgtol(self.bfield.mgtol)
    self.finalized = 1

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
      l = maximum(block.fullloweroverrefinement,self.fulllower)
      u = minimum(block.fullupperoverrefinement,self.fullupper)
      #if alltrue(u >= l):
      if (u[0] >= l[0] and u[1] >= l[1] and u[2] >= l[2]):
        self.children.append(block)
        block.parents.append(self.blocknumber)

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
    for child in self.children:
      child.initializechilddomains()

      # --- Set full domain to negative of child number first.
      l = maximum(self.fulllower,child.fullloweroverrefinement)
      u = child.fullupperoverrefinement
      # --- If the child extends to the edge of the parent mesh, it claims the
      # --- grid points on the upper edges.
      u = u + where(u == self.fullupper,1,0)
      # --- The child claims all unclaimed areas.
      ii = self.getchilddomains(l,u)
      ii[...] = where(ii==self.blocknumber,-child.blocknumber,ii)

      # --- Set interior to positive child number.
      l = maximum(self.fulllower,child.lower/child.refinement)
      u = child.upper/child.refinement
      # --- Check against the case where only the guard cells of the child
      # --- overlap the parent.
      if (u[0] < l[0] or u[1] < l[1] or u[2] < l[2]): continue
      # --- If the child extends to the edge of the parent mesh, it claims the
      # --- grid points on the upper edges.
      u = u + where(u == self.fullupper,1,0)
      # --- The child claims its full interior area
      ii = self.getchilddomains(l,u)
      ii[...] = +child.blocknumber

  def findoverlappingsiblings(self,blocklists):
    """
Recursive routine to find, at each level of refinement, all overlapping
siblings.
This is faster than the original routine since each pair of blocks is checked
only once, rather than twice for each parent as in the original.
    """
    # --- When the list is empty, there are no blocks, so just return.
    if len(blocklists[0]) == 0: return
    # --- Make the call for the next level.
    self.findoverlappingsiblings(blocklists[1:])

    # --- Get a copy of the list (which will be mangled below).
    blocklistscopy = copy.copy(blocklists[0])

    # --- Loop over all blocks.
    for block in blocklists[0]:

      # --- Loop only over blocks that havn't been checked yet.
      del blocklistscopy[0]
      for sibling in blocklistscopy:

        # --- Get the area common to the block and its sibling.
        # --- Don't do anything if there is no overlap.
        sl = maximum(sibling.fulllower,block.fulllower)
        su = minimum(sibling.fullupper,block.fullupper)
        if sl[0] > su[0] or sl[1] > su[1] or sl[2] > su[2]: continue

        if block.blocknumber < sibling.blocknumber:
          block.overlapshigher[sibling.blocknumber] = [sl,su]
          sibling.overlapslower[block.blocknumber] = [sl,su]
        else:
          block.overlapslower[sibling.blocknumber] = [sl,su]
          sibling.overlapshigher[block.blocknumber] = [sl,su]

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
      cl = maximum(child.fullloweroverrefinement,l)
      cu = minimum(child.fullupperoverrefinement,u)
      #if sometrue(cl > cu): continue
      if cl[0] > cu[0] or cl[1] > cu[1] or cl[2] > cu[2]: continue

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

  def loadj(self,lzero=true,depositallparticles=0,lrootonly=0):
    """
Loads the current density from the particles. This should only be called for
the top level grid.
    """
    # --- Make sure that the final setup was done.
    self.finalize()
    # check if subcycling or selfb is turned on for at least one species
    # and create rhospecies if not done already
    if (len(compress(top.pgroup.ndts>1,top.pgroup.ndts))>0):
      if not self.__dict__.has_key('jspecies'):
        self.idts = compress(top.pgroup.ndts>1,arange(top.pgroup.ns))
        self.createjspecies(lrootonly)
        for js in range(top.pgroup.ns):
          if top.pgroup.ndts[js]>1:
            self.initjspecies(1,js,lrootonly)
            i1 = top.pgroup.ins[js]-1
            i2 = top.pgroup.ins[js]+top.pgroup.nps[js]-1
            x=top.pgroup.xp[i1:i2]
            y=top.pgroup.yp[i1:i2]
            z=top.pgroup.zp[i1:i2]
            ux=top.pgroup.uxp[i1:i2]
            uy=top.pgroup.uyp[i1:i2]
            uz=top.pgroup.uzp[i1:i2]
            gaminv=top.pgroup.gaminv[i1:i2]
            self.setj(x,y,z,ux,uy,uz,gaminv,None,
                      top.pgroup.sq[js],top.pgroup.sw[js],lrootonly)
            self.initjspecies(2,js,lrootonly)
    
    # zeros arrays if needed   
    if lzero: self.zeroj(lrootonly)

    # deposit current
    for js,i,n,q,w in zip(arange(top.pgroup.ns),top.pgroup.ins-1,
                          top.pgroup.nps,top.pgroup.sq,top.pgroup.sw):
      if n == 0: continue
      # deposit species according to subcycling rule governed by ndts.
      if top.pgroup.ndts[js]>1:
        if (top.it+1)%top.pgroup.ndts[js]==0: 
          self.pointjtojspecies(js,lrootonly)
          self.zeroj(lrootonly)
        else:
          self.addjspecies(js,lrootonly)
          continue      
      print "calling setj ",self.blocknumber
      self.setj(top.pgroup.xp[i:i+n],top.pgroup.yp[i:i+n],top.pgroup.zp[i:i+n],
                top.pgroup.uxp[i:i+n],top.pgroup.uyp[i:i+n],
                top.pgroup.uzp[i:i+n],top.pgroup.gaminv[i:i+n],None,q,w,
                lrootonly)
      print "setj done ",self.blocknumber
      if (top.pgroup.ndts[js]>1 and top.pgroup.ldts[js]):
        self.pointjtojcopy(lrootonly)
        self.addjspecies(js,lrootonly)
    # distribute current density among blocks
    if not lrootonly and lzero:
      self.propagatejbetweenpatches(depositallparticles)
    self.makejperiodic()
    self.getjforfieldsolve()

  def propagatejbetweenpatches(self,depositallparticles):
      self.getjfromoverlaps()
      if not depositallparticles:
        self.gatherjfromchildren_fortran()
      self.restorejinoverlaps()

  def createjspecies(self,lrootonly=0):
    self.jspecies = {}
    for i in self.idts:
      self.jspecies[i] = fzeros(shape(self.bfield.j),'d')
    if not lrootonly:
      for child in self.children:
        child.idts  =self.idts
        child.createjspecies()

  def initjspecies(self,op,js,lrootonly=0):
    if op==1:
      self.jcopy    = self.bfield.j
      self.bfield.j = self.jspecies[js]
      self.bfield.j[...] = 0.
    if op==2:
      self.bfield.j     = self.jcopy
    if not lrootonly:
      for child in self.children:
        child.initjspecies(op,js)

  def pointjtojspecies(self,js,lselfonly=0):
    # make jcopy point to j, j point to jspecies[js] 
    self.jcopy  = self.bfield.j
    self.bfield.j      = self.jspecies[js]
    if not lselfonly:
      for child in self.children:
        child.pointjtojspecies(js)

  def addjspecies(self,js,lselfonly=0):
    # add jspecies[js] to j
    self.bfield.j += self.jspecies[js]
    if not lselfonly:
      for child in self.children:
        child.addjspecies(js)

  def pointjtojcopy(self,lselfonly=0):
    # make j point to jcopy 
    self.bfield.j     = self.jcopy
    if not lselfonly:
      for child in self.children:
        child.pointjtojcopy()

  def zeroj(self,lselfonly=0):
    if not self.isfirstcall(): return
    self.bfield.j[...] = 0.
    if not lselfonly:
      for child in self.children:
        child.zeroj()

  def sortbyichild(self,ichild,x,y,z,ux,uy,uz,gaminv):
    xout,yout,zout,uxout,uyout,uzout,gout = zeros((7,len(x)),'d')
    nperchild = zeros(self.root.totalnumberofblocks)
    sortposandvelbyindex(len(x),ichild,x,y,z,ux,uy,uz,gaminv,
                         self.root.totalnumberofblocks,
                         xout,yout,zout,uxout,uyout,uzout,gout,nperchild)

    return xout,yout,zout,uxout,uyout,uzout,gout,nperchild

  def setj(self,x,y,z,ux,uy,uz,gaminv,wght,q,w,lrootonly):
    """
Given the list of particles, a charge and a weight, deposits the current
density of the mesh structure.
This first gets the blocknumber of the block where each of the particles are
to be deposited. This is then sorted once. The loop is then over the list
of blocks, rather than walking through the tree structure.
    """
    if len(self.children) > 0 and not lrootonly:

      ichild = zeros(len(x))
      if self.lcylindrical:
        xx = sqrt(x**2 + y**2)
        yy = zeros(len(x),'d')
      else:
        xx = x
        yy = y
      self.getichild_allsort(xx,yy,z,ichild)

      x,y,z,ux,uy,uz,gaminv,nperchild = self.sortbyichild(ichild,x,y,z,
                                                          ux,uy,uz,gaminv)
      blocklist = self.root.listofblocks

    else:
      nperchild = [len(x)]
      blocklist = [self]

    # --- For each block, pass to it the particles in it's domain.
    i = 0
    for block,n in zip(blocklist,nperchild):
      MagnetostaticMG.setj(block,x[i:i+n],y[i:i+n],z[i:i+n],
                           ux[i:i+n],uy[i:i+n],uz[i:i+n],gaminv[i:i+n],
                           None,q,w)
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

  def gatherjfromchildren_fortran(self):
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
      child.gatherjfromchildren_fortran()

      # --- Get coordinates of child relative to this domain
      l = maximum(child.fullloweroverrefinement,self.fulllower)
      u = minimum(child.fullupperoverrefinement,self.fullupper)

      # --- Check for any Nuemann boundaries
      dopbounds = (((child.pbounds[0] == 1 or
                    child.pbounds[1] == 1 or
                    child.pbounds[2] == 1) and
                  (l[0] == 0 or l[1] == 0 or l[2] == 0 or
                   u[0] == self.rootdims[0] or
                   u[1] == self.rootdims[1] or
                   u[2] == self.rootdims[2])) or
                   self.ny == 0)

      w = self.getwarrayforj(child.refinement)
      gatherjfromchild(self.bfield.j,self.dims,child.bfield.j,child.dims,
                         l,u,self.fulllower,child.fulllower,child.fullupper,
                         child.refinement,w,
                         self.xmesh,child.xmesh,self.lcylindrical,
                         dopbounds,child.pbounds,self.rootdims)

    # --- zerorhoinoverlap is call here so that any contribution from
    # --- the children in the overlap regions will get zeroed as necessary.
    self.zerojinoverlap()

  def getjfromoverlaps(self):
    """
Add in the j from overlaping areas. The j is gathered into the block with
the lowerest number. Later on, the j will be copied back to the higher
numbered blocks. This should only ever be called by the root block.
    """
    assert self is self.root,"This should only be called by the root block"
    # --- This loops over the blocks in ascending order to ensure that in any
    # --- area with overlap, the block with the lowest number is the one that
    # --- gets the j. This avoids problems of double counting j. This
    # --- could also be done be zeroing out oj, but that is extra
    # --- (unecessary) computational work, since it already will be done
    # --- at the end of gatherjfromchildren.
    for block in self.listofblocks:
      for othernumber,overlapdomain in block.overlapshigher.items():
        other = block.getblockfromnumber(othernumber)
        l,u = overlapdomain
        sj = block.getj(l,u)
        oj = other.getj(l,u)
        add(sj,oj,sj)

  def restorejinoverlaps(self):
    """
Restore j in overlapping areas for blocks which had the j zeroed out, the
higher numbered blocks. This should only ever be called by the root block.
    """
    assert self is self.root,"This should only be called by the root block"
    # --- The loop does not need to be in ascending order, but this just
    # --- matches the getjfromoverlaps routine.
    for block in self.listofblocks:
      for othernumber,overlapdomain in block.overlapslower.items():
        other = block.getblockfromnumber(othernumber)
        l,u = overlapdomain
        sj = block.getj(l,u)
        oj = other.getj(l,u)
        sj[...] = oj

  def zerojinoverlap(self):
    """
This zeros out j in overlapping regions for higher numbered blocks.  When
j is passed from child to parent, in any overlapping regions only one child
needs to pass the data to the parent.  The choice as to which does the
passing is determined by the blocknumber - the lower gets to do the passing.
For the others, the j in the overlapping region is cleared out. That j
will be restored later by a call to restorejinoverlaps.
Note that this is not recursive, since it is called separately by each block
from gatherjfromchildren.
    """
    for othernumber,overlapdomain in self.overlapslower.items():
      other = self.getblockfromnumber(othernumber)
      l,u = overlapdomain
      sj = self.getj(l,u)
      sj[...] = 0.

  #--------------------------------------------------------------------------
  # --- Methods to carry out the field solve
  #--------------------------------------------------------------------------

  def solve(self,iwhich=0):

    # --- Make sure that the final setup was done.
    self.finalize()

    # --- Wait until all of the parents have called here until actually
    # --- doing the solve. This ensures that the phi in all of the parents
    # --- which is needed on the boundaries will be up to date.
    if not self.islastcall(): return

    # solve on a
    self.setafromparents()
    MagnetostaticMG.solve(self,iwhich)

    # solve for children
    for child in self.children:
      child.solve(iwhich)

  def setafromparents(self):
    """
Sets a, using the values from the parent grid. Setting the full a array
gives a better initial guess for the field solver.
    """
    for parentnumber in self.parents:
      parent = self.getblockfromnumber(parentnumber)
      # --- Coordinates of mesh relative to parent's mesh location
      # --- and refinement. The minimum and maximum are needed in case
      # --- this mesh extends beyond the parent's.
      l = maximum(parent.fulllower*self.refinement,self.fulllower)
      u = minimum(parent.fullupper*self.refinement,self.fullupper)
      # --- The full a arrays are passed in to avoid copying the subsets
      # --- since the fortran needs contiguous arrays.
      gatherafromparents(self.bfield.a,self.dims,l,u,self.fulllower,
                           parent.bfield.a,parent.dims,parent.fulllower,
                           self.refinement)

  def optimizeconvergence(self,resetpasses=1):
    pass
    
  #--------------------------------------------------------------------------
  # --- Methods to fetch B-fields and potential
  #--------------------------------------------------------------------------

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
    xout,yout,zout = zeros((3,len(x)),'d')
    isort = zeros(len(x))
    nperchild = zeros(self.root.totalnumberofblocks)
    sortparticlesbyindexgetisort(len(x),ichild,x,y,z,
                                 self.root.totalnumberofblocks,
                                 xout,yout,zout,isort,nperchild)

    return xout,yout,zout,isort,nperchild

  def fetchb(self):
    """
Fetches the B field. This should only be called at the root level grid.
    """
    if w3d.api_xlf2:
      w3d.xfsapi=top.pgroup.xp[w3d.ipminfsapi-1:w3d.ipminfsapi-1+w3d.ipfsapi]
      w3d.yfsapi=top.pgroup.yp[w3d.ipminfsapi-1:w3d.ipminfsapi-1+w3d.ipfsapi]
      w3d.zfsapi=top.pgroup.zp[w3d.ipminfsapi-1:w3d.ipminfsapi-1+w3d.ipfsapi]
    self.fetchbfrompositions_allsort(w3d.xfsapi,w3d.yfsapi,w3d.zfsapi,
                                     w3d.bxfsapi,w3d.byfsapi,w3d.bzfsapi)

  def fetchbfrompositions_allsort(self,x,y,z,bx,by,bz):
    """
Given the list of particles, fetch the B fields.
This first gets the blocknumber of the block where each of the particles are
to be deposited. This is then sorted once. The loop is then over the list
of blocks, rather than walking through the tree structure.
Note that this depends on having ichilddomains filled with the
blocknumber rather than the child number relative to the parent.
    """
    if len(self.children) > 0:

      ichild = zeros(len(x))
      if self.lcylindrical:
        xx = sqrt(x**2 + y**2)
        yy = zeros(len(x),'d')
      else:
        xx = x
        yy = y
      # --- This assumes that the root block has blocknumber zero.
      self.getichild_positiveonly(xx,yy,z,ichild)

      x,y,z,isort,nperchild = self.sortbyichildgetisort(ichild,x,y,z)

      # --- Create temporary arrays to hold the B field
      tbx,tby,tbz = zeros((3,len(x)),'d')

    else:
      isort = None
      nperchild = [len(x)]
      tbx,tby,tbz = bx,by,bz

    # --- For each block, pass to it the particles in it's domain.
    i = 0
    for block,n in zip(self.root.listofblocks,nperchild):
      MagnetostaticMG.fetchbfrompositions(block,x[i:i+n],y[i:i+n],z[i:i+n],
                                          tbx[i:i+n],tby[i:i+n],tbz[i:i+n])
      i = i + n

    # --- Now, put the B fields back into the original arrays, unsorting
    # --- the data
    if isort is not None:
      n = len(x)
      putsortedefield(len(tbx),isort,tbx,tby,tbz,bx[:n],by[:n],bz[:n])

  def fetcha(self):
    """
Fetches the vector potential. This should only be called at the root level grid.
    """
    if w3d.api_xlf2:
      w3d.xfsapi=top.pgroup.xp[w3d.ipminfsapi-1:w3d.ipminfsapi-1+w3d.ipfsapi]
      w3d.yfsapi=top.pgroup.yp[w3d.ipminfsapi-1:w3d.ipminfsapi-1+w3d.ipfsapi]
      w3d.zfsapi=top.pgroup.zp[w3d.ipminfsapi-1:w3d.ipminfsapi-1+w3d.ipfsapi]
    self.fetchafrompositions_allsort(w3d.xfsapi,w3d.yfsapi,w3d.zfsapi,w3d.afsapi)

  def fetchafrompositions_allsort(self,x,y,z,a):
    """
Given the list of particles, fetch the vector potential.
This first gets the blocknumber of the block where each of the particles are
to be deposited. This is then sorted once. The loop is then over the list
of blocks, rather than walking through the tree structure.
Note that this depends on having ichilddomains filled with the
blocknumber rather than the child number relative to the parent.
    """
    if len(self.children) > 0:

      ichild = zeros(len(x))
      if self.lcylindrical:
        xx = sqrt(x**2 + y**2)
        yy = zeros(len(x),'d')
      else:
        xx = x
        yy = y
      # --- This assumes that the root block has blocknumber zero.
      self.getichild_positiveonly(xx,yy,z,ichild)

      x,y,z,isort,nperchild = self.sortbyichildgetisort(ichild,x,y,z)

      # --- Create temporary arrays to hold the vector potential
      ta = zeros((3,len(x)),'d')

    else:
      isort = None
      nperchild = [len(x)]
      ta = a

    # --- For each block, pass to it the particles in it's domain.
    i = 0
    for block,n in zip(self.root.listofblocks,nperchild):
      MagnetostaticMG.fetchbfrompositions(block,x[i:i+n],y[i:i+n],z[i:i+n],a[:,i:i+n])
      i = i + n

    # --- Now, put the vector potential back into the original arrays, unsorting
    # --- the data
    if isort is not None:
      n = len(x)
      putsortedefield(len(tbx),isort,ta[0,:],ta[1,:],ta[2,:],a[0,:n],a[1,:n],a[2,:n])

  def fetcha_old(self):
    """
Fetches the potential. This should only be called at the root level grid.
    """
    self.fetchafrompositions(w3d.xfsapi,w3d.yfsapi,w3d.zfsapi,w3d.afsapi)

  def fetchafrompositions_old(self,x,y,z,a):
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
      if self.lcylindrical:
        xx = sqrt(x**2 + y**2)
        yy = zeros(len(x),'d')
      else:
        xx = x
        yy = y
      getichildpositiveonly(self.blocknumber,len(x),xx,yy,z,ichild,
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
        ta = zeros(len(3,xc),'d')
        # --- Now get the field
        if block == self:
          MagnetostaticMG.fetchafrompositions(self,xc,yc,zc,ta)
        else:
          block.fetchafrompositions(xc,yc,zc,ta)
        # --- Put the potential into the passed in arrays
        putabyindex(len(ii),ii,len(a),a,ta)
        #put(a,ii,ta)

    else:

      # --- Get a from this domain
      MagnetostaticMG.fetchafrompositions(self,x,y,z,a)

  #--------------------------------------------------------------------------
  # --- Utility methods
  #--------------------------------------------------------------------------

  def getblockfromnumber(self,number):
    return self.root.listofblocks[number]

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
      
  def getwarrayforj(self,r):
    # --- Create weight array needed for j deposition.
    # --- Is linear falloff in the weights correct for r > 2?
    wi = [0,0,0]
    for i in range(3):
      wi[i] = r[i] - abs(1.*iota(-r[i]+1,+r[i]-1))
      wi[i] = wi[i]/sum(wi[i])
    # --- Expand into 3 dimensions
    w = outerproduct(wi[0],outerproduct(wi[1],wi[2]))
    w.shape = (2*r[0]-1,2*r[1]-1,2*r[2]-1)
    result = fzeros(2*r-1,'d')
    result[...] = w
    return result

  def geta(self,lower,upper,comp=slice(None)):
    # --- Note that this takes into account the guard cells
    ix1,iy1,iz1 = lower - self.fulllower + 1
    ix2,iy2,iz2 = upper - self.fulllower + 2 
    if type(comp) == StringType: comp = ['x','y','z'].index(comp)
    return self.bfield.a[comp,ix1:ix2,iy1:iy2,iz1:iz2]
  def getj(self,lower,upper,r=[1,1,1],comp=slice(None)):
    ix1,iy1,iz1 = lower - self.fulllower
    ix2,iy2,iz2 = upper - self.fulllower + 1
    if type(comp) == StringType: comp = ['x','y','z'].index(comp)
    return self.bfield.j[comp,ix1:ix2:r[0],iy1:iy2:r[1],iz1:iz2:r[2]]
  def getb(self,lower=None,upper=None,r=[1,1,1],comp=slice(None)):
    if lower is None: lower = self.lower
    if upper is None: upper = self.upper
    if type(comp) == StringType: comp = ['x','y','z'].index(comp)
    ix1,iy1,iz1 = lower - self.fulllower
    ix2,iy2,iz2 = upper - self.fulllower + 1
    return self.bfield.b[comp,ix1:ix2:r[0],iy1:iy2:r[1],iz1:iz2:r[2]]
  def getchilddomains(self,lower,upper,upperedge=0):
    if self.childdomains is None:
      self.childdomains = fzeros(1+self.dims)
      add(self.childdomains,self.blocknumber,self.childdomains)
    ix1,iy1,iz1 = lower - self.fulllower
    ix2,iy2,iz2 = upper - self.fulllower + upperedge
    return self.childdomains[ix1:ix2,iy1:iy2,iz1:iz2]
  def getlocalarray(self,array,lower,upper,r=[1,1,1],fulllower=None,
                    upperedge=1):
    if fulllower is None: fulllower = self.fulllower
    ix1,iy1,iz1 = lower - fulllower
    ix2,iy2,iz2 = upper - fulllower + upperedge
    return array[ix1:ix2:r[0],iy1:iy2:r[1],iz1:iz2:r[2]]

  def setmgtol(self,mgtol=None):
    """
Sets the convergence tolerance for all blocks. If mgtol is not given, it uses
f3d.bfield.mgtol.
    """
    if mgtol is None: mgtol = f3d.bfield.mgtol
    self.bfield.mgtol = mgtol
    self.mgtol = mgtol
    for child in self.children:
      child.setmgtol(mgtol)

  def setmgmaxiters(self,mgmaxiters=None):
    """
Sets the maximum number of iterations for all blocks. If mgmaxiters is
not given, it uses f3d.bfield.mgmaxiters.
    """
    if mgmaxiters is None: mgmaxiters = f3d.bfield.mgmaxiters
    self.bfield.mgmaxiters = mgmaxiters
    self.mgmaxiters = mgmaxiters
    for child in self.children:
      child.setmgmaxiters(mgmaxiters)

  def arraysliceoperation(self,ip,idim,getdataname,op,opnd,null,comp=None):
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
    getdata = getattr(self,getdataname)
    array = getdata(self.fulllower,self.fullupper,comp=comp)
    if len(self.children) > 0:
      # --- Skip points that don't self doesn't own
      c = self.getchilddomains(self.fulllower,self.fullupper,1)
      array = where(c[ix,iy,iz]==self.blocknumber,array[ix,iy,iz],null)
    else:
      # --- The transpose is taken so that the array is in C ordering so
      # --- the opnd will be faster.
      array = transpose(array[ix,iy,iz])
    # --- Find the max of self's and the children's phi
    result = opnd(array)
    for child in self.children:
      ipc = ip*child.refinement[idim]
      cresult = child.arraysliceoperation(ipc,idim,getdataname,op,opnd,null,
                                          comp)
      result = op(result,cresult)
    return result

  #--------------------------------------------------------------------------
  # --- The following are used for plotting.
  #--------------------------------------------------------------------------

  def genericpf(self,getdataname,kw,idim,pffunc,ip=None):
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

      comp = kw.get('comp',2)
      if type(comp) != IntType:
        try:
          comp = ['x','y','z'].index(comp)
        except ValueError:
          pass
        if self.lcylindrical:
          try:
            comp = ['r','theta','z'].index(comp)
          except ValueError:
            pass
      assert type(comp) == IntType,"Unrecognized component was input"

      if cmin is None:
        cmin = self.arraysliceoperation(ip,idim,getdataname,min,minnd,+largepos,
                                        comp)
      if cmax is None:
        cmax = self.arraysliceoperation(ip,idim,getdataname,max,maxnd,-largepos,
                                        comp)
      cmin = globalmin(cmin)
      cmax = globalmax(cmax)
      kw['cmin'] = cmin
      kw['cmax'] = cmax
      kw['local'] = 1

      accumulateplotlists()

    try:
      # --- Only make the plot if the plane is included in the domain.
      # --- Even if there is no overlap, the children must be called since
      # --- they may overlap (in the domain of a different parent).
      # --- Note that the full extent is used so that the childdomains slice
      # --- will have the same shape as ireg.
      if self.fulllower[idim] <= ip and ip <= self.fullupper[idim]:
        if self.childdomains is not None:
          # --- Create the ireg array, which will be set to zero in the domain
          # --- of the children.
          ss = list(shape(self.bfield.j[0,...]))
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
        MagnetostaticMG.genericpf(self,kw,pffunc)
        kw['titles'] = 0
        kw['lcolorbar'] = 0

      for child in self.children:
        child.genericpf(getdataname,kw,idim,pffunc,ip)

    finally:
      if self is self.root: plotlistofthings(lturnofflist=1)

  def pcaxy(self,kwdict=None,**kw):
    if kwdict is None: kwdict = {}
    kwdict.update(kw)
    self.genericpf('geta',kwdict,2,pcaxy)
  def pcazx(self,kwdict=None,**kw):
    if kwdict is None: kwdict = {}
    kwdict.update(kw)
    self.genericpf('geta',kwdict,1,pcazx)
  def pcazr(self,kwdict=None,**kw):
    if kwdict is None: kwdict = {}
    kwdict.update(kw)
    kwdict['fullplane'] = 0
    self.genericpf('geta',kwdict,1,pcazx)
  def pcazy(self,kwdict=None,**kw):
    if kwdict is None: kwdict = {}
    kwdict.update(kw)
    self.genericpf('geta',kwdict,0,pcazy)
  def pcaztheta(self,kwdict=None,**kw):
    if kwdict is None: kwdict = {}
    kwdict.update(kw)
    kwdict['fullplane'] = 0
    self.genericpf('geta',kwdict,0,pcazy)
  def pcbxy(self,kwdict=None,**kw):
    if kwdict is None: kwdict = {}
    kwdict.update(kw)
    self.genericpf('getb',kwdict,2,pcbxy)
  def pcbzx(self,kwdict=None,**kw):
    if kwdict is None: kwdict = {}
    kwdict.update(kw)
    self.genericpf('getb',kwdict,1,pcbzx)
  def pcbzr(self,kwdict=None,**kw):
    if kwdict is None: kwdict = {}
    kwdict.update(kw)
    kwdict['fullplane'] = 0
    self.genericpf('getb',kwdict,1,pcbzx)
  def pcbzy(self,kwdict=None,**kw):
    if kwdict is None: kwdict = {}
    kwdict.update(kw)
    self.genericpf('getb',kwdict,0,pcbzy)
  def pcbztheta(self,kwdict=None,**kw):
    if kwdict is None: kwdict = {}
    kwdict.update(kw)
    kwdict['fullplane'] = 0
    self.genericpf('getb',kwdict,0,pcbzy)
  def pcjxy(self,kwdict=None,**kw):
    if kwdict is None: kwdict = {}
    kwdict.update(kw)
    self.genericpf('getj',kwdict,2,pcjxy)
  def pcjzx(self,kwdict=None,**kw):
    if kwdict is None: kwdict = {}
    kwdict.update(kw)
    self.genericpf('getj',kwdict,1,pcjzx)
  def pcjzr(self,kwdict=None,**kw):
    if kwdict is None: kwdict = {}
    kwdict.update(kw)
    kwdict['fullplane'] = 0
    self.genericpf('getj',kwdict,1,pcjzx)
  def pcjzy(self,kwdict=None,**kw):
    if kwdict is None: kwdict = {}
    kwdict.update(kw)
    self.genericpf('getj',kwdict,0,pcjzy)
  def pcjztheta(self,kwdict=None,**kw):
    if kwdict is None: kwdict = {}
    kwdict.update(kw)
    kwdict['fullplane'] = 0
    self.genericpf('getj',kwdict,0,pcjzy)




  def plaz(self,ix=None,iy=None,comp=2,color='fg',selfonly=0):
    if ix < self.fulllower[0]: return
    if iy < self.fulllower[1]: return
    if ix > self.fullupper[0]: return
    if iy > self.fullupper[1]: return
    if self is self.root: accumulateplotlists()
    try:
      plg(self.bfield.a[comp,ix-self.fulllower[0]+1,iy-self.fulllower[1]+1,1:-1],self.zmeshlocal,
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
      plg(self.bfield.a[comp,1:-1,iy-self.fulllower[1]+1,iz-self.fulllower[2]+1],self.xmesh,
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
      plg(self.bfield.a[comp,ix-self.fulllower[0]+1,1:-1,iz-self.fulllower[2]+1],self.ymesh,
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
      plg(self.bfield.j[comp,ix-self.fulllower[0],iy-self.fulllower[1],:],self.zmeshlocal,
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
      plg(self.bfield.j[comp,:,iy-self.fulllower[1],iz-self.fulllower[2]],self.xmesh,
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
      plg(self.bfield.j[comp,ix-self.fulllower[0],:,iz-self.fulllower[2]],self.ymesh,
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
      plg(self.bfield.b[comp,ix-self.fulllower[0],iy-self.fulllower[1],iz],
          self.zmeshlocal[iz],color=color)
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
      plg(self.bfield.b[comp,ix,iy-self.fulllower[1],iz-self.fulllower[2]],
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
      plg(self.bfield.b[comp,ix-self.fulllower[0],iy,iz-self.fulllower[2]],
          self.ymesh[iy],color=color)
      if not selfonly:
        for child in self.children:
          child.plby(comp,ix*child.refinement[0],iz*child.refinement[2],
                         color=color,withguard=withguard)
    finally:
      if self is self.root: plotlistofthings(lturnofflist=1)

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
    if self is self.root: accumulateplotlists()
    try:
      plg(yy,xx,color=color[0])
      if not selfonly:
        for child in self.children:
          child.drawbox(ip=ip*child.refinement[idim],idim=idim,
                        withguards=withguards,color=color[1:])
    finally:
      if self is self.root: plotlistofthings(lturnofflist=1)

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
    if self is self.root: accumulateplotlists()
    try:
      plfp(ibox,yy,xx,[5])
      if not selfonly:
        for child in self.children:
          child.drawfilledbox(ip=ip*child.refinement[idim],idim=idim,
                              withguards=withguards,ibox=ibox)
    finally:
      if self is self.root: plotlistofthings(lturnofflist=1)

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





