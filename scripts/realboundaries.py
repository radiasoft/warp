from __future__ import generators
from warp import *
from generateconductors import *
from particlescraper import *
import AMR
import cPickle
realboundaries_version = "$Id: realboundaries.py,v 1.72 2008/11/19 18:30:00 dave Exp $"

##############################################################################
def realboundariesdoc():
  print """
Automatically generates the appropriate boundary conditions for the
slice code using the information supplied in the lattice description. To
use, execute the following command (before the generate is best)

rb_object = RealBoundary(newmesh=1,rodfract=.5)

The function RealBoundary returns an instance of the class that deals
with the transverse boundaries. It has two optional arguments:
  - newmesh When true, a new field mesh is generated which is just big
            enough for the conductor to fit into. Otherwise (the default
            case), the mesh size (i.e. xmmin, xmmax etc.) is not changed.
  - rodfract Only applies when newmesh is true, it is the fraction of
             the quadrupole rods which are included. It defaults to 0.5.

After a restart dump is restarted, the following must be called, passing in
the realboundaries object.
restartRealBoundary

There are two additional functions, saveRealBoundary and
restoreRealBoundary. Use as follows...

saveRealBoundary(rb_object,"rb_datafile")

saves the object in the file rb_datafile, including all of the matrices
which have already be calculated.

rb_object = restoreRealBoundary("rb_datafile")

restores the object in a new run. This saves the time to calculate the
matrices.

The RealBoundary has a publically accessible method, plotcond, which plots
the conductor points in a nice plot. Type doc(rb_object.plotcond) for more
information.

Some additional notes. For the eletrostatic elements, the boundary
produced is always the interdigitated rod structure. For the magnetic
and drifting elements, the boundary produced is always a round pipe. For
the elements which can describe either (quad and hele), a check is made
of whether the electric field component is nonzero, and if so treat as
an electric element, otherwise magnetic. Also important is the order in
which the lattice elements are checked. When the slice location is found
within an element, the routine applies that elements boundary condition
and ignores the rest of the elements. The order in which they are
checked is as follows...
quads, accls, emlts, mmlts, pgrds, egrds, bgrds, heles, bends,
dipos, sexts, drfts
Note that the author can not anticipate all possibilities, so if one of
the assumptions above doesn't work for you, please contact me.

Great pains were taken to minimize the wasted calculation time. For
example, if the boundary conditions of two elements are the same, only
one capacity matrix will be generated which will be shared by the two
elements. Also, matrices are only calculated when needed.
  """
##############################################################################

##############################################################################
# --- Capacity matrix classes
class CapacityMatrix:
  """
CapacityMatrix
Constructor arguments
  - x,y coordinates of conductor points (in meters)
  - v voltage at points, can be a scalar or an array. When it is a scalar, all
      the points have the same voltage (defaults to zero)
      The voltage can also be input in the call to setmatrix so the same
      matrix can be used for conductors with the same geometry but different
      voltages.
  - m optional matrix. When not included, one is automatically calculated.
  - newmesh when true, reset the mesh to be big enough to contain all of the
            conductor points (implies that the matrix will be calculated)
  """
  #----------------------------------------------------------------------------
  def __init__(s,x,y,v=0.,m=None,newmesh=0):
    # --- Check for errors in the input
    if len(x) != len(y):
      raise "Coordinate arrays must be of the same length"
    if type(v) == type(x) and len(v) != len(x):
      raise "Voltage array must be of the same length as the coordinate arrays"
    if m:
      n = shape(m)
      if n[0] != n[1]:
        raise "Input matrix must be square"
      if n[0] != len(x):
        raise "Size of matrix must be same as length of coordinate array"
    if min(s.ingrid(x,y)) == 0:
      raise "All conductor points must be within the grid"
    # --- Save the input into the class
    s.xcond = x
    s.ycond = y
    s.vcond = v*ones(len(x))
    s.newmesh = newmesh
    # --- If a matrix is input, save it, otherwise, calculate one.
    if m:
      s.cmatxy = m
    else:
      if s.newmesh and len(x) > 0:
        s.getmesh(s,min(x),max(x),min(y),max(y),ave(x),ave(y))
      s.getmatrix()
  #----------------------------------------------------------------------------
  def getmatrix(s):
    # --- Get number of conductor points.
    n = len(s.xcond)
    if n > 0:
      fxy.ncxy = n
      fxy.ncxymax = n
      # --- Point fortran arrays to coordinate and voltage data
      fxy.xcond = s.xcond
      fxy.ycond = s.ycond
      fxy.vcond = s.vcond
      # --- Create space for work arrays (in CapMatxy) and point the fortran
      # --- arrays to it.
      s.pcond = fzeros(n,'d')
      s.qcond = fzeros(n,'d')
      s.kpvtxy = fzeros(n,'l')
      fxy.pcond = s.pcond
      fxy.qcond = s.qcond
      fxy.kpvtxy = s.kpvtxy
      # --- Create space for the matrix.
      s.cmatxy = fzeros((n,n),'d')
      fxy.cmatxy = s.cmatxy
      # --- Calculate matrix
      top.fstype = 1
      fieldsol(1)
    else:
      top.fstype = 0
  #----------------------------------------------------------------------------
  def setmatrix(s,vcond=None):
    # --- Reset the mesh if requested
    if s.newmesh: s.setmesh()
    # --- Get number of conductor points.
    n = len(s.xcond)
    if n > 0:
      fxy.ncxy = n
      fxy.ncxymax = n
      # --- Check for new input voltages, otherwise use s.vcond
      if vcond is None: vcond = s.vcond
      # --- Point the fortran arrays to the data arrays. Each instance of this
      # --- object contains its own copies of the arrays in the group CapMatxy.
      # --- The pointers are repointed rather than copying the data.
      fxy.xcond = s.xcond
      fxy.ycond = s.ycond
      fxy.vcond = vcond
      fxy.pcond = s.pcond
      fxy.qcond = s.qcond
      fxy.kpvtxy = s.kpvtxy
      fxy.cmatxy = s.cmatxy
      # --- Turn on capacity matrix field solver
      top.fstype = 1
    else:
      # --- If no points, turn off capacity matrix field solver
      top.fstype = 0
    # --- Recalculate the fields XXX This is now not needed since setmatrix
    #fieldsol(-1)                    is calledimmediately before a field solve
  #----------------------------------------------------------------------------
  def ingrid(s,x,y):
    # --- Determines whether the data points are within the grid.
    xfuzz = 1.e-5*w3d.dx
    yfuzz = 1.e-5*w3d.dy
    if w3d.l4symtry:
      return where(logical_and(logical_and(greater(x,0.),
                                           greater(w3d.xmmax-xfuzz,x)),
                               logical_and(greater(y,0.),
                                           greater(w3d.ymmax-yfuzz,y))),1,0)
    elif w3d.l2symtry:
      return where(logical_and(logical_and(greater(x,w3d.xmmin+xfuzz),
                                           greater(w3d.xmmax-xfuzz,x)),
                               logical_and(greater(y,0.),
                                           greater(w3d.ymmax-yfuzz,y))),1,0)
    else:
      return where(logical_and(logical_and(greater(x,w3d.xmmin+xfuzz),
                                           greater(w3d.xmmax-xfuzz,x)),
                               logical_and(greater(y,w3d.ymmin+yfuzz),
                                           greater(w3d.ymmax-yfuzz,y))),1,0)
  #----------------------------------------------------------------------------
  def insymmetricgrid(s,x,y):
    # --- Determines whether the data points are within the grid.
    xfuzz = 1.e-5*w3d.dx
    yfuzz = 1.e-5*w3d.dy
    if w3d.l4symtry:
      return where(logical_and(logical_and(greater(x,-w3d.xmmax+xfuzz),
                                           greater(w3d.xmmax-xfuzz,x)),
                               logical_and(greater(y,-w3d.ymmax+yfuzz),
                                           greater(w3d.ymmax-yfuzz,y))),1,0)
    elif w3d.l2symtry:
      return where(logical_and(logical_and(greater(x,w3d.xmmin+xfuzz),
                                           greater(w3d.xmmax-xfuzz,x)),
                               logical_and(greater(y,-w3d.ymmax+yfuzz),
                                           greater(w3d.ymmax-yfuzz,y))),1,0)
    else:
      return where(logical_and(logical_and(greater(x,w3d.xmmin+xfuzz),
                                           greater(w3d.xmmax-xfuzz,x)),
                               logical_and(greater(y,w3d.ymmin+yfuzz),
                                           greater(w3d.ymmax-yfuzz,y))),1,0)
  #----------------------------------------------------------------------------
  def getmesh(s,xmin,xmax,ymin,ymax,xcent,ycent):
    # --- This routine redefines the mesh to match the size of the element.
    # --- Extend far enough to cover the conductor plus a little extra space to
    # --- keep the conducting points away from the mesh edge.
    # --- Get number of grid points as if there was no symmetry
    nx = w3d.nx
    ny = w3d.ny
    if w3d.l4symtry: nx = 2*nx
    if w3d.l2symtry or w3d.l4symtry: ny = 2*ny
    # --- Now, calculate the grid cell sizes, assuming that there is an extra
    # --- grid cell beyond the pipe
    w3d.dx = (xmax - xmin)/(nx - 2)
    w3d.dy = (ymax - ymin)/(ny - 2)
    # --- Correct mins and maxes to include the extra grid cell and symmetries.
    w3d.xmmin = xmin - w3d.dx
    w3d.xmmax = xmax + w3d.dx
    w3d.ymmin = ymin - w3d.dy
    w3d.ymmax = ymax + w3d.dy
    if w3d.l4symtry: w3d.xmmin = 0.
    if w3d.l2symtry or w3d.l4symtry: w3d.ymmin = 0.
    # --- Recalculate mesh arrays
    w3d.xmesh = iota(0,w3d.nx)*w3d.dx + w3d.xmmin
    w3d.ymesh = iota(0,w3d.ny)*w3d.dy + w3d.ymmin
    # --- Recalculate ksqx and ksqy (notice that for this part, the capacity
    # --- matrix is turned off since the matrix does not need to be calculated)
    ff = top.fstype
    top.fstype = 0
    fieldsol(1)
    top.fstype = ff
    # --- Now, save all of the pertinent mesh data
    s.nx = w3d.nx
    s.ny = w3d.ny
    s.dx = w3d.dx
    s.dy = w3d.dy
    s.xmmin = w3d.xmmin
    s.xmmax = w3d.xmmax
    s.ymmin = w3d.ymmin
    s.ymmax = w3d.ymmax
    s.xmesh = w3d.xmesh + 0.
    s.ymesh = w3d.ymesh + 0.
    s.kxsq = w3d.kxsq + 0.
    s.kysq = w3d.kysq + 0.
  #----------------------------------------------------------------------------
  def setmesh(s):
    # --- Restore the values to w3d
    # --- Currently, assumes constant nx and ny
    #w3d.nx = s.nx
    #w3d.ny = s.ny
    w3d.dx = s.dx
    w3d.dy = s.dy
    w3d.xmmin = s.xmmin
    w3d.xmmax = s.xmmax
    w3d.ymmin = s.ymmin
    w3d.ymmax = s.ymmax
    w3d.xmesh[:] = s.xmesh
    w3d.ymesh[:] = s.ymesh
    w3d.kxsq[:] = s.kxsq
    w3d.kysq[:] = s.kysq
    # --- Reload the charge density
    w3d.rho = 0.
    loadrho()
  #----------------------------------------------------------------------------
  def issame(s,x,y):
    # --- Checks if the input values would return the same matrix
    if max(abs(s.xcond - x)) > 0.: return 0
    if max(abs(s.ycond - y)) > 0.: return 0
    return 1


##############################################################################
_initialmeshparams = None
class NoPipe(CapacityMatrix):
  """
No boundary condition.
  - newmesh when true, reset the mesh to be the original values
  """
  #----------------------------------------------------------------------------
  def __init__(s,newmesh=0):
    s.xcond = []
    s.ycond = []
    s.vcond = []
    s.newmesh = newmesh
    if s.newmesh:
      # --- Use initial value of the mesh params
      s.getmesh(_initialmeshparams[0],_initialmeshparams[1],
                _initialmeshparams[2],_initialmeshparams[3],
                _initialmeshparams[4],_initialmeshparams[5])
  #----------------------------------------------------------------------------
  def setmatrix(s,v=None):
    CapacityMatrix.setmatrix(s)
  #----------------------------------------------------------------------------
  def issame(s,ap,ax,ay,ox,oy):
    return 1
  #----------------------------------------------------------------------------
  def plotcond(s,color='fg',xscale=1.,yscale=1.):
    pass

##############################################################################
class RoundPipe(CapacityMatrix):
  """
Capacity matrix for a round pipe
Constructor arguments:
  - ap pipe radius
  - v voltage on the pipe
  - ox,oy position of the center of the pipe, defaults to (0,0)
  - newmesh when true, reset the mesh to be big enough to contain all of the
            conductor points (implies that the matrix will be calculated)
  """
  #----------------------------------------------------------------------------
  def __init__(s,ap,ax=0.,ay=0.,v=0.,ox=0.,oy=0.,newmesh=0):
    ap,ax,ay = s.setaperture(ap,ax,ay)
    # --- Check that the radius is positive
    assert ax > 0. and ay > 0.,"Aperture data is inconsistent - either ap must > 0 or both of ax and ay must be > 0"
    # --- Save the input values
    s.ap = ap
    s.ax = ax
    s.ay = ay
    s.v = v
    s.ox = ox
    s.oy = oy
    s.newmesh = newmesh
    # --- Check to make sure that the parameters are OK.
    assert s.paramsok(ap,ax,ay,ox,oy,newmesh),\
           "Error in input - either aperture is zero or too big to fit mesh"
    # --- Reset the mesh if requested.
    if s.newmesh:
      s.getmesh(ox-ax,ox+ax,oy-ay,oy+ay,ox,oy)
    # --- Get the symmetry factor
    symm_fact = 1.
    if (w3d.l2symtry): symm_fact = 0.5
    if (w3d.l4symtry): symm_fact = 0.25
    # --- Now, estimate what a good number of points would be.
    # --- The minimum distance between two points that gaurantess that they
    # --- are not in the same cell is sqrt(dx**2+dy**2). Divide the
    # --- circumference into pieces of that size.
    # --- The circumference of an ellipse is approximated using a formula
    # --- from S. Ramanujan.
    h = (ax - ay)**2/(ax + ay)**2
    circum = pi*(ax + ay)*(1. + 3.*h/(10. + sqrt(4. - 3.*h)))
    n = int(symm_fact*circum/sqrt(w3d.dx**2 + w3d.dy**2))
    # --- Now, get the points on the round pipe
    s.xcond = ax*cos(2.*pi*symm_fact*arange(n)/n) + ox
    s.ycond = ay*sin(2.*pi*symm_fact*arange(n)/n) + oy
    s.vcond = s.v*ones(n,'d')
    # --- and the matrix using those points
    s.getmatrix() 
  #----------------------------------------------------------------------------
  def setaperture(s,ap,ax,ay):
    # --- ax and ay take precedence over ap.
    if ax == 0.: ax = ap
    if ay == 0.: ay = ap
    if ap == 0.: ap = ax
    return ap,ax,ay
  #----------------------------------------------------------------------------
  def setmatrix(s,v=None):
    if v is not None: s.vcond[:] = v
    CapacityMatrix.setmatrix(s,s.vcond)
  #----------------------------------------------------------------------------
  def issame(s,ap,ax,ay,ox,oy):
    ap,ax,ay = s.setaperture(ap,ax,ay)
    # --- Checks if the input values would return the same matrix
    if s.ap == ap and s.ax == ax and s.ay == ay:
      if s.ox == ox and s.oy == oy:
        return 1
    return 0
  #----------------------------------------------------------------------------
  def paramsok(s,ap,ax,ay,ox,oy,newmesh):
    "Check validity of parameter"
    ap,ax,ay = s.setaperture(ap,ax,ay)
    if newmesh:
      # --- If mesh will be resized, only check that aperture is > 0.
      if (ap > 0. or (ax > 0. and ay > 0.)): return 1
      else:                                  return 0
    else:
      # --- If mesh will not be resized, need to check if the aperture will
      # --- fit within the mesh.
      if (ap == 0. and (ax == 0. or ay == 0.)): return 0
      # --- Check if the whole pipe circle is within the grid
      xmax = w3d.xmmax - w3d.dx
      ymax = w3d.ymmax - w3d.dy
      xmin = w3d.xmmin + w3d.dx
      ymin = w3d.ymmin + w3d.dy
      if w3d.l4symtry: xmin = -w3d.xmmax + w3d.dx
      if w3d.l2symtry or w3d.l4symtry: ymin = -w3d.ymmax + w3d.dy
      if ox+ax > xmax or ox-ax < xmin or oy+ay > ymax or oy-ay < ymin:
        return 0
      else:
        return 1
  #----------------------------------------------------------------------------
  def plotcond(s,color='fg',xscale=1.,yscale=1.):
    tt = span(0.,2*pi,101)
    xx = s.ax*cos(tt-pi*0.75) + s.ox
    yy = s.ay*sin(tt-pi*0.75) + s.oy
    ii = s.insymmetricgrid(xx,yy)
    xx = compress(ii,xx)
    yy = compress(ii,yy)
    plg(yy*yscale,xx*xscale,color=color)


##############################################################################
class RoundRods(CapacityMatrix):
  """
Capacity for ESQ rods
Constructor arguments:
  - ap pole tip aperture
  - rr rod radius, defaults to 8/7*a
  - vx, vy voltages on rods in x and y planes respectively, defaults to 0.
  - vxm, vym optional voltages on rods at negative values of x and y,
             default to vx and vy
  - withx, withy flags signalling whether or not to includes the rods in the
                 x and y planes respectively, default to 1
  - ox,oy position of the center of the quadrupole
  - newmesh when true, reset the mesh to be big enough to contain all of the
            conductor points (implies that the matrix will be calculated)
  """
  #----------------------------------------------------------------------------
  def __init__(s,ap,rr=None,vx=0.,vy=0.,vxm=None,vym=None,withx=1,withy=1,
               ox=0.,oy=0.,newmesh=0,rodfract=0.5):
    # --- Save some input values
    s.ap = ap
    s.rr_in = rr
    s.withx = withx
    s.withy = withy
    s.ox = ox
    s.oy = oy
    s.newmesh = newmesh
    s.rodfract = rodfract
    # --- Set default value of the rod radius
    if not rr:
      s.rr = 8./7.*s.ap
    else:
      s.rr = rr
    # --- Only need to make sure that the aperture is greater than zero.
    if s.ap <= 0.:
      raise "Aperture must be greater then zero"
    # --- Reset the mesh if requested.
    if s.newmesh:
      # --- Set gridmax to be big enough to include the specified amount
      # --- of the quadrupole rods
      gridmax = s.ap + s.rr*2*s.rodfract
      s.getmesh(s.ox-gridmax,s.ox+gridmax,s.oy-gridmax,s.oy+gridmax,s.ox,s.oy)
    else:
      # --- The grid cell sizes are calculated here explicitly in case
      # --- they havn't yet been calculated.
      w3d.dx = (w3d.xmmax - w3d.xmmin)/w3d.nx
      w3d.dy = (w3d.ymmax - w3d.ymmin)/w3d.ny
    # --- Check and save voltages
    if vxm is None: vxm = vx
    if vym is None: vym = vy
    s.vx = vx
    s.vy = vy
    s.vxm = vxm
    s.vym = vym
    # --- Set vx = 1, vy = -1. They are later set to the actual values.
    # --- It is done this way so that the voltages can change.
    # --- The integers 0-3 are used to select the appropriate value in
    # --- in the 'choose' statement in setmatrix.
    vx  = 0
    vy  = 1
    vxm = 2
    vym = 3
    # --- Now, estimate what a good number of points would be.
    # --- The minimum distance between two points that gaurantess that they
    # --- are not in the same cell is sqrt(dx**2+dy**2). Take as the minimum
    # --- angle then sqrt(dx**2+dy**2)/r, and divide the circle up evenly.
    n = int(2*pi*s.rr/sqrt(w3d.dx**2 + w3d.dy**2))
    # --- Now calculate the points for each rod (including the full circle).
    # --- Lists are used to allow concatenation.
    x = []
    y = []
    v = []
    # --- First the x rods
    if withx:
      # --- Positive x
      x = x + list(s.rr*cos(2.*pi*arange(n)/n) + s.ox + s.ap + s.rr)
      y = y + list(s.rr*sin(2.*pi*arange(n)/n) + s.oy)
      v = v + n*[vx]
      # --- Negative x
      if not w3d.l4symtry:
        x = x + list(s.rr*cos(2.*pi*arange(n)/n) + s.ox - s.ap - s.rr)
        y = y + list(s.rr*sin(2.*pi*arange(n)/n) + s.oy)
        v = v + n*[vxm]
    if withy:
      # --- Positive y
      x = x + list(s.rr*cos(2.*pi*arange(n)/n) + s.ox)
      y = y + list(s.rr*sin(2.*pi*arange(n)/n) + s.oy + s.ap + s.rr)
      v = v + n*[vy]
      # --- Negative y
      if not (w3d.l2symtry or w3d.l4symtry):
        x = x + list(s.rr*cos(2.*pi*arange(n)/n) + s.ox)
        y = y + list(s.rr*sin(2.*pi*arange(n)/n) + s.oy - s.ap - s.rr)
        v = v + n*[vym]
    # --- Now select only those points within the grid
    ii = s.ingrid(x,y)
    s.xcond = compress(ii,x)
    s.ycond = compress(ii,y)
    s.vcond = compress(ii,v)*1.
    s.vcondunit = int(s.vcond + 0)
    # --- Finally, get the matrix using those points
    s.getmatrix()
  #----------------------------------------------------------------------------
  def setmatrix(s,vlist=None):
    if vlist is not None:
      # --- The exception handling is needed in cases when a RoundRod class
      # --- instance is restore from a run done before vcondunit was added.
      try:
        s.vcond[:] = 1.*choose(s.vcondunit,vlist)
      except AttributeError:
        pass
    CapacityMatrix.setmatrix(s,s.vcond)
  #----------------------------------------------------------------------------
  def issame(s,ap,rr,withx,withy,ox,oy):
    # --- Checks if the input values would return the same matrix
    if s.ap == ap and (s.rr == rr or s.rr_in == rr):
      if s.withx == withx and s.withy == withy and s.ox == ox and s.oy == oy:
        return 1
    return 0
  #----------------------------------------------------------------------------
  def plotcond(s,color='fg',xscale=1.,yscale=1.):
    tt = span(0.,2*pi,101)
    for sx,sy,da,w in [(+1,0,0.,s.withx),(-1,0,pi,s.withx),
                       (0,+1,+pi/2.,s.withy),(0,-1,-pi/2.,s.withy)]:
      if not w: continue
      xx = s.rr*cos(tt+da) + sx*(s.ap + s.rr) + s.ox
      yy = s.rr*sin(tt+da) + sy*(s.ap + s.rr) + s.oy
      ii = s.insymmetricgrid(xx,yy)
      x = compress(ii,xx)
      y = compress(ii,yy)
      plg(y*yscale,x*xscale,color=color)

##############################################################################
##############################################################################
_realboundarycount = 0
class RealBoundary:
  """
Applies appropriate boundaries using capacity matrix solver for the slice code.
Constructor arguments:
  - newmesh=0: when true, the grid is redefined to ensure that the elements
               boundaries are included within the grid
  - rodfract=0.5: when newmesh is true, this gives the fraction of the rods
                  which are included for electric quads.
  - scrapemethod: Method to use to scrape particles. Only applied to 3d.
                  1 uses prwall to scrape on cylindrical boundaries
                  2 (default) uses ParticleScraper module to scrape on
                    actual conductors. The method is precise but can be very
                    inefficient.
                  This input supercedes lscrapeparticles.
  - lscrapeparticles=1: When true, particles are scraped on the conductors.
                        This only applies to the 3d case. (Obsolete)
  - scrapermglevel=1: Coarsening level for index grid used to locate which
                      conductors particles are near. See doc(ParticleScraper)
                      for more info.
  - scraperargs={}: Dictionary of arguments to pass to the scraper creator
                    when scrapemethod==2.
  - dfill=2: parameter passed to installconductors in 3d. sets how much the
             interior of conductors are filled.
  - lclearconductors=1: When true, in the 3d case, all conductor data is cleared
                        out. If false, it is up to the user to clear out the
                        data.
  - pipethickness=largepos: Thickness of surrounding pipes. Only used in 3d.
                            Can be used for efficiency - if the pipe
                            completely surrounds the beam, then it doesn't
                            need much thickness. Setting this can greatly
                            reduce the amount of data generated.
  - conductors=[]: List of additional conductors to include. This list can be
                   further modified using the addconductor and removeconductor
                   methods.
  """
  #----------------------------------------------------------------------------
  def __init__(self,newmesh=0,rodfract=0.5,lscrapeparticles=1,scrapermglevel=1,
                    scrapemethod=2,scraperargs={},
                    dfill=2,lclearconductors=1,pipethickness=largepos,
                    conductors=[]):
    global _realboundarycount
    # --- Only allow one instance of this class.
    if _realboundarycount > 0:
      raise "There can only be one instance of the RealBoundary class"
    _realboundarycount = 1
    self._realboundarycount = _realboundarycount
    # --- Save the input arguments
    self.newmesh = newmesh
    self.rodfract = rodfract
    self.scrapemethod = scrapemethod
    self.lscrapeparticles = lscrapeparticles
    self.scraperargs = scraperargs
    if 'mglevel' not in self.scraperargs: self.scraperargs['mglevel'] = scrapermglevel
    self.dfill = dfill
    self.lclearconductors = lclearconductors
    self.pipethickness = pipethickness
    self.conductors = conductors

    # --- Keep a global lists of all matrices. With a global lists of matrices,
    # --- recalculation of a matrix can be avoided if one with the same values
    # --- have already been calculated.
    self.nomatrixlist = []
    self.pipematrixlist = []
    self.rodmatrixlist = []

    # --- For each element type, create a list of capacity matrices, which are
    # --- all only initially empty. It is done this way so that if the number
    # --- of elements were to change, this coding would not be affected.
    self.drftcm  = []
    self.quadcm  = []
    self.acclcm  = []
    self.emltcm  = []
    self.mmltcm  = []
    self.pgrdcm  = []
    self.egrdcm  = []
    self.bgrdcm  = []
    self.helecm  = []
    self.bendcm  = []
    self.dipocm  = []
    self.sextcm  = []
    self.nonecm  = []

    # --- Keep track of the current matrix. Knowing the current matrix, the
    # --- time wasted reseting the fxy variables to the same values is avoided.
    self.current = None

    # --- Setup so the initialization routine will be called.
    installbeforefs(self.initialsetboundary)

  #----------------------------------------------------------------------------
  def addconductor(self,conductor):
    "Add a conductor to the list of conductors included"
    self.conductors.append(conductor)
  def removeconductor(self,conductor):
    "Remove a conductor from the list of conductors included"
    self.conductors.remove(conductor)

  #----------------------------------------------------------------------------
  def initialsetboundary(self):
    global _initialmeshparams
    # --- This does initialization that depends on which simulation package
    # --- is being used. This happens the first time a field solve is done
    # --- after this instance was created. Normal operation would have the
    # --- instance created before the generate step (and before the package
    # --- was specified).

    # --- This routine removes itself from the beforefs so it is only
    # --- called once.
    uninstallbeforefs(self.initialsetboundary)

    # --- Fetch the initial values of the mesh parameters. This can only be
    # --- done here since they must have been set by the user at this point.
    if _initialmeshparams is None:
      _initialmeshparams = [w3d.xmmin,w3d.xmmax,w3d.ymmin,w3d.ymmax,
                            0.5*(w3d.xmmax+w3d.xmmin),0.5*(w3d.ymmax+w3d.ymmin)]
      self._initialmeshparams = _initialmeshparams

    currpkg = getcurrpkg()
    if currpkg == 'w3d':
      # --- Make sure that an appropriate field solver is turned on.
      if top.fstype not in [3,7,13]: top.fstype = 7
      # --- Make sure the gridmode is set properly
      f3d.gridmode = 1
    elif currpkg == 'wxy':
      if w3d.solvergeom==w3d.XYgeom:
        top.fstype = 10
      else:
        # --- Make sure that fstype = 0 at this point, since fstype
        # --- is set to the proper value by setboundary. This is done so that in
        # --- case this is called before the generate and a user still is setting
        # --- fstype=1 (or 2).
        top.fstype = 0
      # --- Now make the call. This call is needed since it would noramally
      # --- be called at the beginning of the step, which has already passed.

    # --- Make initial call to setboundary
    self.setboundary()

    # --- Turn on the realboundaries for normal operation.
    self.enable()

  #----------------------------------------------------------------------------
  def __del__(self):
    global _realboundarycount
    # --- Class destructor
    _realboundarycount = 0
    self.disable()
  #----------------------------------------------------------------------------
  def disable(self):
    # --- Turns off the realboundaries
    # --- This must be protected since disable may be called by __del__
    # --- on exit from python, at which point 'top' may not exist.
    try:
      # --- Just try uninstalling everything, independent of what package
      # --- is active.
      if isinstalledbeforefs(self.initialsetboundary):
        uninstallbeforefs(self.initialsetboundary)
      if isinstalledbeforefs(self.setboundary):
        uninstallbeforefs(self.setboundary)
      top.fstype = 0
    except:
      pass
  #----------------------------------------------------------------------------
  def enable(self):
    # --- Turns on the realboundaries
    # --- The routine setboundary will be called before every field-solve.
    if not isinstalledbeforefs(self.setboundary):
      installbeforefs(self.setboundary)
  #----------------------------------------------------------------------------
  def __setstate__(self,dict):
    global _initialmeshparams,_realboundarycount
    self.__dict__.update(dict)
    if not isinstalledbeforefs(self.initialsetboundary):
      installbeforefs(self.initialsetboundary)
    # --- Restore these quantities which are otherwise not saved.
    try:
      _initialmeshparams = self._initialmeshparams
      _realboundarycount = self._realboundarycount
    except AttributeError:
      _realboundarycount = 1

  #----------------------------------------------------------------------------
  def setscraper(self,zl,zr,ap,ax,ay,ox,oy):
    if self.scrapemethod != 1: return
    izl = int(ceil((zl - top.zbeam - top.zzmin)/top.dzz))
    izr = int(floor((zr - top.zbeam - top.zzmin)/top.dzz))
    izl = max(0,min(top.nzzarr,izl))
    izr = max(0,min(top.nzzarr,izr)) + 1
    if ax == 0.: ax = ap
    if ay == 0.: ay = ap
    top.prwallz[izl:izr] = ap
    top.prwallxz[izl:izr] = ox
    top.prwallyz[izl:izr] = oy
    top.prwelipz[izl:izr] = ay/ax
    top.prwall = max(top.prwallz)

  #----------------------------------------------------------------------------
  def setmatrix(self,m,v):
    if self.current == m: return
    if m is not None: m.setmatrix(v)
    self.current = m
  #----------------------------------------------------------------------------
  def getnomatrix(self):
    # --- This searches through the list of matrices checking if one with the
    # --- same parameters has already be created. If so, just return that one.
    for m in self.nomatrixlist:
      if m.issame(): return m
    # --- None was found, so create a new one, adding it to the list.
    m = NoPipe(self.newmesh)
    self.nomatrixlist.append(m)
    return m
  #----------------------------------------------------------------------------
  def getpipematrix(self,ap,ax,ay,v,ox,oy):
    # --- This searches through the list of matrices checking if one with the
    # --- same parameters has already be created. If so, just return that one.
    for m in self.pipematrixlist:
      if m.issame(ap,ax,ay,ox,oy): return m
    # --- None was found, so create a new one, adding it to the list.
    m = RoundPipe(ap,ax,ay,v,ox,oy,self.newmesh)
    self.pipematrixlist.append(m)
    return m
  #----------------------------------------------------------------------------
  def getrodmatrix(self,ap,rr,vx,vy,vxm,vym,withx,withy,ox,oy):
    # --- This searches through the list of matrices checking if one with the
    # --- same parameters has already be created. If so, just return that one.
    for m in self.rodmatrixlist:
      if m.issame(ap,rr,withx,withy,ox,oy): return m
    # --- None was found, so create a new one, adding it to the list.
    m = RoundRods(ap,rr,vx,vy,vx,vy,withx,withy,ox,oy,self.newmesh,
                  self.rodfract)
    self.rodmatrixlist.append(m)
    return m
  #----------------------------------------------------------------------------
  def nopipe(self,id,cm):
    currpkg = getcurrpkg()
    if currpkg == 'w3d' or (currpkg == 'wxy' and w3d.solvergeom==w3d.XYgeom):
      return 0
    elif currpkg == 'wxy':
      # --- Check if there is a matrix for this element
      if (len(cm) > id and cm[id] is None) or len(cm) < id+1:
        # --- If the list is too short, add some None's in.
        while len(cm) < id+1: cm.append(None)
        # --- Now, add the capacity matrix.
        cm[id] = self.getnomatrix()
      self.setmatrix(cm[id],0.)
      return 1
  #----------------------------------------------------------------------------
  def roundpipe(self,id,zs,ze,ap,ax,ay,ox,oy,cm):
    currpkg = getcurrpkg()
    if currpkg == 'w3d' or (currpkg == 'wxy' and w3d.solvergeom==w3d.XYgeom):
      return self.roundpipe3d(id,zs,ze,ap,ax,ay,ox,oy,cm)
    elif currpkg == 'wxy':
      return self.roundpipexy(id,zs,ze,ap,ax,ay,ox,oy,cm)
  #----------------------------------------------------------------------------
  def roundpipe3d(self,id,zs,ze,ap,ax,ay,ox,oy,cm):
    if (ze < w3d.zmmin+top.zbeam or
        zs > w3d.zmmax+top.zbeam): return 0
    ax = ax[id]
    ay = ay[id]
    ap = ap[id]
    if ax == 0.: ax = ap
    if ay == 0.: ay = ap
    if ax == 0. or ay == 0.: return 0
    if self.pipethickness < largepos:
      pipe = ZAnnulusElliptic(ay/ax,ax,ax+self.pipethickness,ze-zs,0.,
                              ox[id],oy[id],0.5*(zs+ze))
    else:
      pipe = ZCylinderEllipticOut(ay/ax,ax,ze-zs,0.,ox[id],oy[id],0.5*(zs+ze))
    self.newconductors.append(pipe)
    self.setscraper(zs,ze,ap,ax,ay,ox[id],oy[id])
    return 0
  #----------------------------------------------------------------------------
  def roundpipexy(self,id,zs,ze,ap,ax,ay,ox,oy,cm):
    # --- Only apply matrix is the z location is within the current element
    if zs <= top.zbeam < ze:
      self.setscraper(zs,ze,ap,ax,ay,ox,oy)
      # --- Check if there is a matrix for this element
      if (len(cm) > id and cm[id] is None) or len(cm) < id+1:
        # --- If the list is too short, add some None's in.
        while len(cm) < id+1: cm.append(None)
        # --- Now, add the capacity matrix.
        try:
          cm[id] = self.getpipematrix(ap[id],ax[id],ay[id],0.,ox[id],oy[id])
        except:
          return 0
      self.setmatrix(cm[id],0.)
      return 1
    return None
  #----------------------------------------------------------------------------
  def quadrods(self,elem,elemid,zs,ze,cm):
    currpkg = getcurrpkg()
    if currpkg == 'w3d' or (currpkg == 'wxy' and w3d.solvergeom==w3d.XYgeom):
      return self.quadrods3d(elem,elemid,zs,ze)
    elif currpkg == 'wxy':
      return self.quadrodsxy(elem,elemid,zs,ze,cm)
  #----------------------------------------------------------------------------
  def quadrods3d(self,elem,elemid,zs,ze):
    rl = getattr(top,elem+'rl')[elemid]
    gl = getattr(top,elem+'gl')[elemid]
    pw = getattr(top,elem+'pw')[elemid]
    zc = 0.5*(zs+ze)
    zl = zc-0.5*(rl+gl)-pw
    zr = zc+0.5*(rl+gl)+pw
    if (zr < w3d.zmmin+top.zbeam or
        zl > w3d.zmmax+top.zbeam): return 0
    quad = Quadrupole(zcent=zc,condid=100+elemid,elem=elem,elemid=elemid,
                      splitrodids=1)
    ap = getattr(top,elem+'ap')[elemid]
    ax = getattr(top,elem+'ax')[elemid]
    ay = getattr(top,elem+'ay')[elemid]
    if elem == 'quad':
      ox = top.qoffx[elemid]
      oy = top.qoffy[elemid]
    else:
      ox = getattr(top,elem+'ox')[elemid]
      oy = getattr(top,elem+'oy')[elemid]
    self.setscraper(zl,zr,ap,ax,ay,ox,oy)
    if quad is not None:
      self.newconductors.append(quad)
    return 0
  #----------------------------------------------------------------------------
  def quadrodsxy(self,elem,elemid,zs,ze,cm):
    ap = getattr(top,elem+'ap')[elemid]
    ax = getattr(top,elem+'ax')[elemid]
    ay = getattr(top,elem+'ay')[elemid]
    rl = getattr(top,elem+'rl')[elemid]
    rr = getattr(top,elem+'rr')[elemid]
    gl = getattr(top,elem+'gl')[elemid]
    gp = getattr(top,elem+'gp')[elemid]
    pa = getattr(top,elem+'pa')[elemid]
    pw = getattr(top,elem+'pw')[elemid]
    if elem == 'quad':
      vx = top.quadvx[elemid]
      vy = top.quadvy[elemid]
      ox = top.qoffx[elemid]
      oy = top.qoffy[elemid]
    else:
      vx = 0.
      vy = 0.
      ox = getattr(top,elem+'ox')[elemid]
      oy = getattr(top,elem+'oy')[elemid]

    # --- Find the extent of the element.
    zc = 0.5*(zs + ze)
    if rl > 0.:
      # --- If a rod length is defined, use that.
      zl = zc - 0.5*(rl + gl) - pw
      zr = zc + 0.5*(rl + gl) + pw
    else:
      # --- Otherwise, rely on the zs and ze of the element.
      zl = zs
      zr = ze
    # --- Only apply matrix if the z location is within the current element
    if zl <= top.zbeam < zr and ap > 0.:
      self.setscraper(zl,zr,ap,ax,ay,ox,oy)
      # --- Check if the matrix has already been calculated.
      if (len(cm) > elemid and cm[elemid] is None) or len(cm) < elemid+1:
        # --- If the list is too short, add some None's in.
        while len(cm) < elemid+1: cm.append(None)
        # --- Generate the matrix
        # --- Electric quads: Up to five matrices are generated,
        # --- 2 for the end-plates, 2 for where there are gaps between rods
        # --- and the endplates, and 1 where the rods overlap.
        # --- When one of the matrices is not needed, the element of the
        # --- list is set to None
        cmlist = 5*[None]
        # --- First end-plate
        if pw > 0.:
          if pa == 0.: pa = ap
          if gp > 0.: v1 = vx
          else: v1 = vy
          cmlist[0] = self.getpipematrix(pa,pa,pa,v1,ox,oy)
        # --- First gap location
        if gl > 0.:
          withx = 0
          withy = 0
          if gp > 0.: withx = 1
          else:       withy = 1
          cmlist[1] = self.getrodmatrix(ap,rr,vx,vy,vx,vy,withx,withy,ox,oy)
        # --- Rod overlap
        cmlist[2] = self.getrodmatrix(ap,rr,vx,vy,vx,vy,1,1,ox,oy)
        # --- Second gap location
        if gl > 0.:
          withx = 0
          withy = 0
          if gp > 0.: withy = 1
          else:       withx = 1
          cmlist[3] = self.getrodmatrix(ap,rr,vx,vy,vx,vy,withx,withy,ox,oy)
        # --- Second end-plate
        if pw > 0.:
          if pa == 0.: pa = ap
          if gp > 0.: v2 = vy
          else: v2 = vx
          cmlist[4] = self.getpipematrix(pa,pa,pa,v2,ox,oy)
        # --- Now append the list of 5 to the end of the main list of matrices.
        cm[elemid] = cmlist
      # --- Apply one of the five matrices.
      if zl <= top.zbeam < zl + pw:
        if gp > 0.: v1 = vx
        else:       v1 = vy
        self.setmatrix(cm[elemid][0],v1)
      elif zl+pw <= top.zbeam < zl+pw+gl:
        self.setmatrix(cm[elemid][1],[vx,vy,vx,vy])
      elif zl+pw+gl <= top.zbeam < zr-pw-gl:
        self.setmatrix(cm[elemid][2],[vx,vy,vx,vy])
      elif zr-pw-gl <= top.zbeam < zr-pw:
        self.setmatrix(cm[elemid][3],[vx,vy,vx,vy])
      elif zr-pw <= top.zbeam <= zr:
        if gp > 0.: v2 = vy
        else:       v2 = vx
        self.setmatrix(cm[elemid][4],v2)
      return 1
    return 0
  #----------------------------------------------------------------------------
  def elemlist(self,nol,cid,czs,cze):
    """This little routine returns the info about each element referenced
in the celemid array. It returns each element only once.
    """
    for io in range(nol):
      previd = -1
      prevzs = largepos
      prevze = largepos
      for id,zs,ze in zip(cid[:,io],czs[:,io],cze[:,io]):
        if id == previd or (zs == prevzs and ze == prevze): continue
        previd = id
        prevzs = zs
        prevze = ze
        yield (id,zs,ze)
  #----------------------------------------------------------------------------
  def setboundary(self,lforce=0):
    # --- The lattice element conductors are only setup if the beam frame
    # --- has moved or if a new field solver is being used. This saves the
    # --- previous beam frame location and field solver to compare against
    # --- the current one. Note that only the id of the solver is saved,
    # --- in order to avoid having a full reference to it.
    try:
      self.lastzbeam
    except:
      self.lastzbeam = top.zbeam/2. - 1.
      self.lastsolverid = None
      self.lastsolverparams = None
    solver = getregisteredsolver()
    if solver is None: solver = w3d
    try:
      if solver.__class__ is AMR.AMRTree: solver = solver.blocks
    except AttributeError:
      pass
    solverparams = [solver.nx,solver.ny,solver.nzlocal,
                    solver.dx,solver.dy,solver.dz,
                    solver.xmmin,solver.ymmin,solver.zmminlocal,
                    solver.xmmax,solver.ymmax,solver.zmmaxlocal]
    if (not lforce and
        self.lastzbeam == top.zbeam and
        self.lastsolverid == id(solver) and
        self.lastsolverparams == solverparams): return
    self.lastzbeam = top.zbeam
    self.lastsolverid = id(solver)
    self.lastsolverparams = solverparams

    # --- Reset the list of conductors
    # --- Note that self.conductor is a list of conductors that are always
    # --- installed. Below, new conductors from any lattice elements are
    # --- added to the list.
    self.newconductors = copy.copy(self.conductors)

    if w3d.solvergeom==w3d.XYgeom: del_conductors()

    # --- Clear out the particle scraping radius of using method 1.
    self.setscraper(top.zzmin,top.zzmax,largepos,largepos,largepos,0.,0.)

    # --- Go through each element type and chose the first one covers the
    # --- current location.
    #--------------------------------------------------------------------------
    if top.nquadol > 0:
      for qid,qzs,qze in self.elemlist(top.nquadol,top.cquadid,
                                       top.cquadzs,top.cquadze):
        if (abs(top.quadde[qid]) > 0. or
            abs(top.quadvx[qid]) > 0. or abs(top.quadvy[qid]) > 0. or
            top.quadrr[qid] > 0. or top.quadrl[qid] > 0. or
            top.quadgl[qid] > 0.):
          if self.quadrods('quad',qid,qzs,qze,self.quadcm):
            return
        else:
          if self.roundpipe(qid,qzs,qze,top.quadap,top.quadax,top.quaday,
                            top.qoffx,top.qoffy,self.quadcm):
            return
    #--------------------------------------------------------------------------
    if top.nacclol > 0:
      for aid,azs,aze in self.elemlist(top.nacclol,top.cacclid,
                                       top.cacclzs,top.cacclze):
        if self.roundpipe(aid,azs,aze,top.acclap,top.acclax,top.acclay,
                          top.acclox,top.accloy,self.acclcm):
          return
    #--------------------------------------------------------------------------
    if top.nemltol > 0:
      for eid,ezs,eze in self.elemlist(top.nemltol,top.cemltid,
                                       top.cemltzs,top.cemltze):
        if self.quadrods('emlt',eid,ezs,eze,self.emltcm):
          return
    #--------------------------------------------------------------------------
    if top.nmmltol > 0:
      for mid,mzs,mze in self.elemlist(top.nmmltol,top.cmmltid,
                                       top.cmmltzs,top.cmmltze):
        if top.mmltas[mid] != top.mmltae[mid]:
          # --- Use aperture start and end, adding the same offset that was
          # --- added to mmltzs to get cmmltzs.
          tempmzs = mzs
          mzs = top.mmltas[mid] + (tempmzs - top.mmltzs[mid])
          mze = top.mmltae[mid] + (tempmzs - top.mmltzs[mid])
        if self.roundpipe(mid,mzs,mze,top.mmltap,top.mmltax,top.mmltay,
                          top.mmltox,top.mmltoy,self.mmltcm):
          return
    #--------------------------------------------------------------------------
    if top.npgrdol > 0:
      for pid,pzs,pze in self.elemlist(top.npgrdol,top.cpgrdid,
                                       top.cpgrdzs,top.cpgrdze):
        if self.quadrods('pgrd',pid,pzs,pze,self.pgrdcm):
          return
    #--------------------------------------------------------------------------
    if top.negrdol > 0:
      for eid,ezs,eze in self.elemlist(top.negrdol,top.cegrdid,
                                       top.cegrdzs,top.cegrdze):
        if self.roundpipe(eid,ezs,eze,top.egrdap,top.egrdax,top.egrday,
                          top.egrdox,top.egrdoy,self.egrdcm):
          return
    #--------------------------------------------------------------------------
    if top.nbgrdol > 0:
      for bid,bzs,bze in self.elemlist(top.nbgrdol,top.cbgrdid,
                                       top.cbgrdzs,top.cbgrdze):
        if self.roundpipe(bid,bzs,bze,top.bgrdap,top.bgrdax,top.bgrday,
                          top.bgrdox,top.bgrdoy,self.bgrdcm):
          return
    #--------------------------------------------------------------------------
    if top.nheleol > 0:
      for hid,hzs,hze in self.elemlist(top.nheleol,top.cheleid,
                                       top.chelezs,top.cheleze):
        if (max(abs(top.heleae[:,hid])) > 0. or
            top.helerr[hid] > 0. or top.helerl[hid] > 0. or
            top.helegl[hid] > 0.):
          if self.quadrods('hele',hid,hzs,hze,self.helecm):
            return
        else:
          if self.roundpipe(hid,hzs,hze,top.heleap,top.heleax,top.heleay,
                         top.heleox,top.heleoy,self.helecm):
            return
    #--------------------------------------------------------------------------
    if top.nbendol > 0:
      cbendid = transpose([top.cbendid])
      cbendzs = transpose([top.cbendzs])
      cbendze = transpose([top.cbendze])
      for bid,bzs,bze in self.elemlist(top.nbendol,    cbendid,
                                           cbendzs,    cbendze):
        if self.roundpipe(bid,bzs,bze,top.bendap,top.bendax,top.benday,
                          top.bendox,top.bendoy,self.bendcm):
          return
    #--------------------------------------------------------------------------
    if top.ndipool > 0:
      for did,dzs,dze in self.elemlist(top.ndipool,top.cdipoid,
                                       top.cdipozs,top.cdipoze):
        if self.roundpipe(did,dzs,dze,top.dipoap,top.dipoax,top.dipoay,
                          top.dipoox,top.dipooy,self.dipocm):
          return
    #--------------------------------------------------------------------------
    if top.nsextol > 0:
      for sid,szs,sze in self.elemlist(top.nsextol,top.csextid,
                                       top.csextzs,top.csextze):
        if self.roundpipe(sid,szs,sze,0.,0.,0., #top.sextap,
                          top.sextox,top.sextoy,self.sextcm):
          return
    #--------------------------------------------------------------------------
    # --- Drifts are checked for last in case the other elements might extend
    # --- into the neighboring drifts.
    if top.ndrftol > 0:
      for did,dzs,dze in self.elemlist(top.ndrftol,top.cdrftid,
                                       top.cdrftzs,top.cdrftze):
        if self.roundpipe(did,dzs,dze,top.drftap,top.drftax,top.drftay,
                          top.drftox,top.drftoy,self.drftcm):
          return

    # --- If this part of the code is reached, then there are no applicable
    # --- boundaries, so turn the capacity matrix field solver off.
    self.nopipe(0,self.nonecm)
    # --- This place is always reached in the 3d case.
    if self.lclearconductors:
      if solver is w3d:
        f3d.conductors.interior.n = 0
        f3d.conductors.evensubgrid.n = 0
        f3d.conductors.oddsubgrid.n = 0
      else:
        solver.clearconductors()
    if self.newconductors:
      if solver is w3d:
        installconductors(self.newconductors,dfill=self.dfill,gridmode=1)
      else:
        print "installing conductor"
        solver.installconductor(self.newconductors,dfill=self.dfill)
      if self.scrapemethod == 2 and self.lscrapeparticles:
        try:
          self.scraper.disable()
        except:
          pass
        self.scraper = ParticleScraper(self.newconductors,**self.scraperargs)

  #----------------------------------------------------------------------------
  def plotcond(self,plotphi=1,filled=0,plotedge=1,plotpoints=0,plotsym=1,\
               ccolor='red',ecolor='green',phicolor='fg',\
               xscale=None,yscale=None):
    """
Makes a plot of the conductor.
  - plotphi=1 when true, plots contours of phi
  - filled=0 when true, plot filled contours
  - plotsym=1 when true, plots appropriate reflections when symmeties are used
  - plotedge=1 when true, plots a square around the edge of the mesh (with
               appropriate reflections when plotsym is true)
  - xscale=scale factor for x-coordinates on plot (default = 1) 
  - yscale=scale factor for y-coordinates on plot (default = 1)
    """
    if xscale == None: xscale = 1.
    if yscale == None: yscale = 1.
    xmesh = xscale*w3d.xmesh
    ymesh = yscale*w3d.ymesh
    # --- Plot phi first (so filled contours are on bottom)
    if plotphi:
      plotc(transpose(getphi(iz=0)),ymesh,xmesh,filled=filled,
            color=phicolor)
    if plotsym:
      # --- Plot negative y
      if w3d.l2symtry or w3d.l4symtry:
        if plotphi:
          plotc(transpose(getphi(iz=0)),-ymesh,xmesh,filled=filled,
                color=phicolor)
      # --- Plot negative x
      if w3d.l4symtry:
        if plotphi:
          plotc(transpose(getphi(iz=0)),ymesh,-xmesh,filled=filled,
                color=phicolor)
        if plotphi:
          plotc(transpose(getphi(iz=0)),-ymesh,-xmesh,filled=filled,
                color=phicolor)
    # --- Plot conductor points next.
    if top.fstype == 1:
      self.current.plotcond(ccolor,xscale=xscale,yscale=yscale)
      if plotsym:
        if w3d.l2symtry or w3d.l4symtry:
          self.current.plotcond(color=ccolor,xscale=xscale,yscale=-yscale)
        if w3d.l4symtry:
          self.current.plotcond(color=ccolor,xscale=-xscale,yscale= yscale)
          self.current.plotcond(color=ccolor,xscale=-xscale,yscale=-yscale)

    # --- Plot the edge of the mesh last
    if plotedge:
      xmmax = xscale*w3d.xmmax
      xmmin = xscale*w3d.xmmin
      ymmax = yscale*w3d.ymmax
      ymmin = yscale*w3d.ymmin 
      if w3d.l4symtry and plotsym:
        plg([-ymmax, ymmax,ymmax,-ymmax,-ymmax],
            [-xmmax,-xmmax,xmmax, xmmax,-xmmax],
            color=ecolor)
      elif w3d.l2symtry and plotsym:
        plg([ ymmin, ymmax,ymmax, ymmin, ymmin],
            [-xmmax,-xmmax,xmmax, xmmax,-xmmax],
            color=ecolor)
      else:
        plg([ ymmin, ymmax,ymmax, ymmin, ymmin],
            [ xmmin, xmmin,xmmax, xmmax, xmmin],
            color=ecolor)




##############################################################################
def saveRealBoundary(object,filename):
  """
Save a class to the file filename.
  - object
  - filename
  """
  currpkg = getcurrpkg()
  if currpkg != 'wxy' or (currpkg == 'wxy' and w3d.solvergeom==w3d.XYgeom):return
  # --- This is only do by the first processor on a parallel machine.
  if me == 0:
    ff = open(filename,'w')
    cPickle.dump(object,ff,1)
    ff.close()                                                                   
#-----------------------------------------------------------------------------
def restoreRealBoundary(filename):
  """
Restore a class from the file filename.
  - filename
Returns the object
  """
  currpkg = getcurrpkg()
  # --- XXX Don't know why this check is here - it shouldn't be. This
  # --- function should work before the generate is done, in which case
  # --- the package generally has not been set yet.
  #if currpkg != 'wxy' or (currpkg == 'wxy' and w3d.solvergeom==w3d.XYgeom):return
  ff = open(filename,'r')
  result = cPickle.load(ff)
  ff.close()
  result.enable()
  return result

def restartRealBoundary(object):
  """
When restarting, this function should be called to restart the realboundaries
and to do the correct field-solve.
  - object: is the restored realboundary object
  """
  # --- Disable it in case it was enables by some other means
  object.disable()
  # --- Enable it
  object.enable()
  # --- Set the boundaries for the current location
  object.setboundary()
  # --- Do a field solve with the updated boundaries
  fieldsol(-1)

##############################################################################
