from warp import *
import cPickle
realboundaries_version = "$Id: realboundaries.py,v 1.5 2001/01/16 19:28:02 dave Exp $"

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
quads, accls, emlts, mmlts, mmlt2s, pgrds, bgrds, bgrd2s, heles, bends,
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
      if s.newmesh:
        getmesh(s,min(x),max(x),min(y),max(y),max(x),ave(x),ave(y))
      s.getmatrix()
  #----------------------------------------------------------------------------
  def getmatrix(s):
    # --- Get number of conductor points.
    n = len(s.xcond)
    fxy.ncxy = n
    fxy.ncxymax = n
    # --- Point fortran arrays to coordinate and voltage data
    fxy.forceassign("xcond",s.xcond)
    fxy.forceassign("ycond",s.ycond)
    fxy.forceassign("vcond",s.vcond)
    # --- Create space for work arrays (in CapMatxy) and point the fortran
    # --- arrays to it.
    s.pcond = fzeros(n,'d')
    s.qcond = fzeros(n,'d')
    s.kpvtxy = fzeros(n)
    fxy.forceassign("pcond",s.pcond)
    fxy.forceassign("qcond",s.qcond)
    fxy.forceassign("kpvtxy",s.kpvtxy)
    # --- Create space for the matrix.
    s.cmatxy = fzeros((n,n),'d')
    fxy.forceassign("cmatxy",s.cmatxy)
    # --- Calculate matrix
    fieldsol(1)
  #----------------------------------------------------------------------------
  def setmatrix(s):
    # --- Reset the mesh if requested
    if s.newmesh: s.setmesh()
    # --- Get number of conductor points.
    n = len(s.xcond)
    fxy.ncxy = n
    fxy.ncxymax = n
    # --- Point the fortran arrays to the data arrays. Each instance of this
    # --- object contains its own copies of the arrays in the group CapMatxy.
    # --- The pointers are repointed rather than copying the data.
    fxy.forceassign("xcond",s.xcond)
    fxy.forceassign("ycond",s.ycond)
    fxy.forceassign("vcond",s.vcond)
    fxy.forceassign("pcond",s.pcond)
    fxy.forceassign("qcond",s.qcond)
    fxy.forceassign("kpvtxy",s.kpvtxy)
    fxy.forceassign("cmatxy",s.cmatxy)
    # --- Recalculate the fields
    fieldsol(-1)
  #----------------------------------------------------------------------------
  def ingrid(s,x,y):
    # --- Determines whether the data points are within the grid.
    if w3d.l4symtry:
      return where(logical_and(logical_and(greater(x,0.),
                                           greater(w3d.xmmax,x)),
                               logical_and(greater(y,0.),
                                           greater(w3d.ymmax,y))),1,0)
    elif w3d.l2symtry:
      return where(logical_and(logical_and(greater(x,w3d.xmmin),
                                           greater(w3d.xmmax,x)),
                               logical_and(greater(y,0.),
                                           greater(w3d.ymmax,y))),1,0)
    else:
      return where(logical_and(logical_and(greater(x,w3d.xmmin),
                                           greater(w3d.xmmax,x)),
                               logical_and(greater(y,w3d.ymmin),
                                           greater(w3d.ymmax,y))),1,0)
  #----------------------------------------------------------------------------
  def getmesh(s,xmin,xmax,ymin,ymax,rpart,xcent,ycent):
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
    # --- Reset particle boundaries
    top.prwall = rpart
    top.prwallx = xcent
    top.prwally = ycent
    top.prwallz = top.prwall
    top.prwallxz = top.prwallx
    top.prwallyz = top.prwally
    # --- Recalculate ksqx and ksqy (notice that for this part, the capacity
    # --- matrix is turned off since the matrix does not need to be calculated)
    if top.fstype >= 0:
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
  def issame(s,x,y,v):
    # --- Checks if the input values would return the same matrix
    if max(abs(s.xcond - x)) > 0.: return 0
    if max(abs(s.ycond - y)) > 0.: return 0
    try:
      if max(abs(s.v - v)) > 0.: return 0
    except:
      if abs(s.v - v) > 0.: return 0
    return 1



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
  def __init__(s,ap,v=0.,ox=0.,oy=0.,newmesh=0):
    # --- Check that the radius is positive
    if ap <= 0.: raise "Radius must be greater than zero"
    # --- Save the input values
    s.ap = ap
    s.v = v
    s.ox = ox
    s.oy = oy
    s.newmesh = newmesh
    # --- Reset the mesh if requested. If this is done, there is no need for
    # --- the error checking.
    if s.newmesh:
      s.getmesh(ox-ap,ox+ap,oy-ap,oy+ap,ap,ox,oy)
    else:
      w3d.dx = (w3d.xmmax - w3d.xmmin)/w3d.nx
      w3d.dy = (w3d.ymmax - w3d.ymmin)/w3d.ny
      # --- Check that the whole pipe circle is within the grid
      if (ox+ap > w3d.xmmax - w3d.dx) or (oy+ap > w3d.ymmax - w3d.dy):
        raise "All of round pipe must be within one grid cell of the grid edge"
      if w3d.l2symtry and \
         ((ox-ap < w3d.xmmin + w3d.dx) or (oy-ap < -w3d.ymmax + w3d.dy)):
        raise "All of round pipe must be within one grid cell of the grid edge"
      if w3d.l4symtry and \
         ((ox-ap < -w3d.xmmax + w3d.dx) or (oy-ap < -w3d.ymmax + w3d.dy)):
        raise "All of round pipe must be within one grid cell of the grid edge"
      if (not w3d.l2symtry and not w3d.l4symtry) and \
         ((ox-ap < w3d.xmmin + w3d.dx) or (oy-ap < w3d.ymmin + w3d.dy)):
        raise "All of round pipe must be within one grid cell of the grid edge"
    # --- Get the symmetry factor
    symm_fact = 1.
    if (w3d.l2symtry): symm_fact = 0.5
    if (w3d.l4symtry): symm_fact = 0.25
    # --- Now, estimate what a good number of points would be.
    # --- The minimum distance between two points that gaurantess that they
    # --- are not in the same cell is sqrt(dx**2+dy**2). Take as the minimum
    # --- angle then sqrt(dx**2+dy**2)/r, and divide the circle up evenly.
    n = int(2*pi*symm_fact*ap/sqrt(w3d.dx**2 + w3d.dy**2))
    # --- Now, get the points on the round pipe
    s.xcond = ap*cos(2.*pi*symm_fact*arange(n)/n) + ox
    s.ycond = ap*sin(2.*pi*symm_fact*arange(n)/n) + oy
    s.vcond = s.v*ones(n,'d')
    # --- and the matrix using those points
    s.getmatrix() 
  #----------------------------------------------------------------------------
  def issame(s,ap,v,ox,oy):
    # --- Checks if the input values would return the same matrix
    if s.ap == ap:
      if s.v == v and s.ox == ox and s.oy == oy:
        return 1
    return 0


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
  def __init__(s,ap,rr=None,vx=0.,vy=0.,vxm=None,vym=None,withx=1,withy=1,ox=0.,oy=0.,
               newmesh=0,rodfract=0.5):
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
      s.getmesh(s.ox-gridmax,s.ox+gridmax,s.oy-gridmax,s.oy+gridmax,
                s.ap,s.ox,s.oy)
    else:
      # --- The grid cell sizes are calculated here explicitly in case
      # --- they havn't yet been calculated.
      w3d.dx = (w3d.xmmax - w3d.xmmin)/w3d.nx
      w3d.dy = (w3d.ymmax - w3d.ymmin)/w3d.ny
    # --- Check and save voltages
    if vxm == None: vxm = vx
    if vym == None: vym = vy
    s.vx = vx
    s.vy = vy
    s.vxm = vxm
    s.vym = vym
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
    s.vcond = compress(ii,v)
    # --- Finally, get the matrix using those points
    s.getmatrix()
  #----------------------------------------------------------------------------
  def issame(s,ap,rr,vx,vy,vxm,vym,withx,withy,ox,oy):
    # --- Checks if the input values would return the same matrix
    if s.ap == ap and (s.rr == rr or s.rr_in == rr):
      if s.vx == vx and s.vy == vy and s.vxm == vxm and s.vym == vym and \
         s.withx == withx and s.withy == withy and s.ox == ox and s.oy == oy:
        return 1
    return 0

##############################################################################
##############################################################################
_realboundarycount = 0
class RealBoundary:
  """
Applies appropriate boundaries using capacity matrix solver for the slice code.
Constructor arguments:
  - newmesh=0 when true, the grid is redefined to ensure that the elements
              boundaries are included within the grid
  - rodfract=0.5 when newmesh is true, this gives the fraction of the rods which
                 are included for electric quads.
  """
  #----------------------------------------------------------------------------
  def __init__(s,newmesh=0,rodfract=0.5):
    global _realboundarycount
    # --- Only allow one instance of this class.
    if _realboundarycount > 0:
      raise "There can only be one instance of the RealBoundary class"
    _realboundarycount = 1
    # --- Save the input arguments
    s.newmesh = newmesh
    s.rodfract = rodfract
    # --- Some idiot proofing. Make sure that fstype = 0 at this point. fstype
    # --- is set to the proper value by setboundary. This is done so that in
    # --- case this is called before the generate and a user still is setting
    # --- fstype=1 (or 2). In that case, the capacity matrix solver will be
    # --- initialized, but no conductors have been defined at this point.
    # --- Doing this only avoids a printed notice which may confuse the user.
    top.fstype = 0
    # --- Keep a global lists of all matrices. With a global lists of matrices,
    # --- recalculation of a matrix can be avoided if one with the same values
    # --- have already been calculated.
    s.pipematrixlist = []
    s.rodmatrixlist = []
    # --- For each element type, create a list of capacity matrices, which are
    # --- all only initially empty. It is done this way so that if the number
    # --- of elements were to change, this coding would not be affected.
    s.drftcm  = []
    s.quadcm  = []
    s.acclcm  = []
    s.emltcm  = []
    s.mmltcm  = []
    s.mmlt2cm = []
    s.pgrdcm  = []
    s.bgrdcm  = []
    s.bgrd2cm = []
    s.helecm  = []
    s.bendcm  = []
    s.dipocm  = []
    s.sextcm  = []
    # --- Turn on the realboundaries.
    s.enable()
  #----------------------------------------------------------------------------
  def __del__(s):
    global _realboundarycount
    # --- Class destructor
    _realboundarycount = 0
    s.disable()
  #----------------------------------------------------------------------------
  def disable(s):
    # --- Turns off the realboundaries
    # --- This must be protected since disable may be called by __del__
    # --- on exit from python, at which point 'top' may not exist.
    try:
      top.fstype = 0
    except:
      pass
    try:
      uninstallbeforestep(s.setboundary)
      uninstallbeforefs(s.initialsetboundary)
    except:
      pass
  #----------------------------------------------------------------------------
  def enable(s):
    # --- Turns on the realboundaries
    # --- Keep track of the current matrix. Knowing the current matrix, the
    # --- time wasted reseting the fxy variables to the same values is avoided.
    s.current = None
    # --- The routine setboundary will be called at the start of every
    # --- time step.  It also must be called just be the field solve
    # --- during the generate.  This is done by the initialsetboundary
    # --- routine, which then removes itself from that list.
    installbeforestep(s.setboundary)
    installbeforefs(s.initialsetboundary)
  #----------------------------------------------------------------------------
  def setmatrix(s,m):
    if s.current == m: return
    m.setmatrix()
    s.current = m
  #----------------------------------------------------------------------------
  def getpipematrix(s,ap,v,ox,oy):
    # --- This searches through the list of matrices checking if one with the
    # --- same parameters has already be created. If so, just return that one.
    for m in s.pipematrixlist:
      if m.issame(ap,v,ox,oy): return m
    # --- None was found, so create a new one, adding it to the list.
    m = RoundPipe(ap,v,ox,oy,s.newmesh)
    s.pipematrixlist.append(m)
    return m
  #----------------------------------------------------------------------------
  def getrodmatrix(s,ap,rr,vx,vy,vx,vy,withx,withy,ox,oy):
    # --- This searches through the list of matrices checking if one with the
    # --- same parameters has already be created. If so, just return that one.
    for m in s.rodmatrixlist:
      if m.issame(ap,rr,vx,vy,vx,vy,withx,withy,ox,oy): return m
    # --- None was found, so create a new one, adding it to the list.
    m = RoundRods(ap,rr,vx,vy,vx,vy,withx,withy,ox,oy,s.newmesh,s.rodfract)
    s.rodmatrixlist.append(m)
    return m
  #----------------------------------------------------------------------------
  def roundpipe(s,id,zs,ze,ap,ox,oy,cm):
    # --- Only apply matrix is the z location is within the current element
    if zs <= top.zbeam < ze:
      # --- Check if there is a matrix for this element
      if (len(cm) > id+1 and cm[id] == None) or len(cm) < id+1:
        # --- If the list is too short, add some None's in.
        while len(cm) < id+1: cm.append(None)
        # --- Now, add the capacity matrix.
        if ap[id] > 0.:
          cm[id] = s.getpipematrix(ap[id],0.,ox[id],oy[id])
        else:
          return 0
      s.setmatrix(cm[id])
      return 1
    return None
  #----------------------------------------------------------------------------
  def quadrods(s,id,zs,ze,ap,rr,rl,gl,gp,vx,vy,pa,pw,ox,oy,cm):
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
      # --- Check if the matrix has already been calculated.
      if (len(cm) > id+1 and cm[id] == None) or len(cm) < id+1:
        # --- If the list is too short, add some None's in.
        while len(cm) < id+1: cm.append(None)
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
          cmlist[0] = s.getpipematrix(pa,v1,ox,oy)
        # --- First gap location
        if gl > 0.:
          withx = 0
          withy = 0
          if gp > 0.: withy = 1
          else: withx = 1
          cmlist[1] = s.getrodmatrix(ap,rr,vx,vy,vx,vy,withx,withy,ox,oy)
        # --- Rod overlap
        cmlist[2] = s.getrodmatrix(ap,rr,vx,vy,vx,vy,1,1,ox,oy)
        # --- Second gap location
        if gl > 0.:
          withx = 0
          withy = 0
          if gp > 0.: withx = 1
          else: withx = 1
          cmlist[3] = s.getrodmatrix(ap,rr,vx,vy,vx,vy,withx,withy,ox,oy)
        # --- Second end-plate
        if pw > 0.:
          if pa == 0.: pa = ap
          if gp > 0.: v2 = vy
          else: v2 = vx
          cmlist[4] = s.getpipematrix(pa,v2,ox,oy)
        # --- Now append the list of 5 to the end of the main list of matrices.
        cm[id] = cmlist
      # --- Apply one of the five matrices.
      if zl <= top.zbeam < zl + pw:
        s.setmatrix(cm[id][0])
      elif zl+pw <= top.zbeam < zl+pw+gl:
        s.setmatrix(cm[id][1])
      elif zl+pw+gl <= top.zbeam < zr-pw-gl:
        s.setmatrix(cm[id][2])
      elif zr-pw-gl <= top.zbeam < zr-pw:
        s.setmatrix(cm[id][3])
      elif zr-pw <= top.zbeam <= zr:
        s.setmatrix(cm[id][4])
      return 1
    return 0
  #----------------------------------------------------------------------------
  def initialsetboundary(s):
    # --- This does an initial call to setboundary during the generate. It is
    # --- called just before the field solve.
    # --- This is only needed for the case where an instance of this class is
    # --- created before the generate. During the generate, the capacity
    # --- matrix field solver needs to be initialized but the step command
    # --- has not been called (which would normally make a call to setboundary).
    # --- This routine should only be called during the generate, so steps
    # --- are taken below to ensure that this is the case before the call
    # --- to setboundary is made.
    # --- This routine removes itself from the beforefs list since during
    # --- the run, the capacity matrices are reset at the start of the step.
    uninstallbeforefs(s.initialsetboundary)
    # --- If the generate has already be done (i.e. top.it > 0) then return
    # --- since setboundary doesn't need to be called.
    if top.it > 0: return
    # --- Now make the call.
    s.setboundary()
  #----------------------------------------------------------------------------
  def setboundary(s):
    # --- This routine should only be called by the wxy package
    if package()[0] != 'wxy': return
    # --- Make sure that the capacity matrix field solver is turned on. This is
    # --- here so that in the case where an instance of this class is created
    # --- before the generate, the initial call to the fieldsolver will have
    # --- fstype=0, the default value. Also, if there is a section without
    # --- any defined boundaries, fstype was set to 0.
    top.fstype = 1
    # --- Go through each element type and chose the first one covers the
    # --- current location.
    #--------------------------------------------------------------------------
    if top.quads:
      qid = top.cquadid[0]
      if (abs(top.quadde[qid]) > 0. or
          abs(top.quadvx[qid]) > 0. or abs(top.quadvy[qid]) > 0.):
        if s.quadrods(qid,top.cquadzs[0],top.cquadze[0],top.quadap[qid],
                      top.quadrr[qid],top.quadrl[qid],top.quadgl[qid],
                      top.quadgp[qid],top.quadvx[qid],top.quadvy[qid],
                      top.quadpa[qid],top.quadpw[qid],
                      top.qoffx[qid],top.qoffy[qid],s.quadcm):
          return
      else:
        if s.roundpipe(top.cquadid[0],top.cquadzs[0],top.cquadze[0],top.quadap,
                       top.qoffx,top.qoffy,s.quadcm):
          return
    #--------------------------------------------------------------------------
    if top.accls:
      if s.roundpipe(top.cacclid[0],top.cacclzs[0],top.cacclze[0],top.acclap,
                     top.acclox,top.accloy,s.acclcm):
        return
    #--------------------------------------------------------------------------
    if top.emlts:
      eid = top.cemltid[0]
      if s.quadrods(eid,top.cemltzs[0],top.cemltze[0],top.emltap[eid],
                    top.emltrr[eid],top.emltrl[eid],top.emltgl[eid],
                    top.emltgp[eid],0.,0.,top.emltpa[eid],top.emltpw[eid],
                    top.emltox[eid],top.emltoy[eid],s.emltcm):
        return
    #--------------------------------------------------------------------------
    if top.mmlts:
      if s.roundpipe(top.cmmltid[0],top.cmmltzs[0],top.cmmltze[0],top.mmltap,
                     top.mmltox,top.mmltoy,s.mmltcm):
        return
    #--------------------------------------------------------------------------
    if top.mmlt2s:
      if s.roundpipe(top.cmmlt2id[0],top.cmmlt2zs[0],top.cmmlt2ze[0],
                     top.mmlt2ap,top.mmlt2ox,top.mmlt2oy,s.mmlt2cm):
        return
    #--------------------------------------------------------------------------
    if top.pgrds:
      pid = top.cpgrdid[0]
      if s.quadrods(pid,top.cpgrdzs[0],top.cpgrdze[0],top.pgrdap[pid],
                    top.pgrdrr[pid],top.pgrdrl[pid],top.pgrdgl[pid],
                    top.pgrdgp[pid],0.,0.,top.pgrdpa[pid],top.pgrdpw[pid],
                    top.pgrdox[pid],top.pgrdoy[pid],s.pgrdcm):
        return
    #--------------------------------------------------------------------------
    if top.bgrds:
      if s.roundpipe(top.cbgrdid[0],top.cbgrdzs[0],top.cbgrdze[0],top.bgrdap,
                     top.bgrdox,top.bgrdoy,s.bgrdcm):
        return
    #--------------------------------------------------------------------------
    if top.bgrd2s:
      if s.roundpipe(top.cbgrd2id[0],top.cbgrd2zs[0],top.cbgrd2ze[0],
                     top.bgrd2ap,top.bgrd2ox,top.bgrd2oy,s.bgrd2cm):
        return
    #--------------------------------------------------------------------------
    if top.heles:
      hid = top.cheleid[0]
      if max(abs(top.heleae[:,hid])) > 0.:
        if s.quadrods(hid,top.chelezs[0],top.cheleze[0],top.heleap[hid],
                      top.helerr[hid],top.helerl[hid],top.helegl[hid],
                      top.helegp[hid],0.,0.,top.helepa[hid],top.helepw[hid],
                      top.heleox[hid],top.heleoy[hid],s.helecm):
          return
      else:
        if s.roundpipe(top.cheleid[0],top.chelezs[0],top.cheleze[0],top.heleap,
                       top.heleox,top.heleoy,s.acclcm):
          return
    #--------------------------------------------------------------------------
    if top.bends:
      if s.roundpipe(top.cbendid[0],top.cbendzs[0],top.cbendze[0],top.bendap,
                     top.bendox,top.bendoy,s.bendcm):
        return
    #--------------------------------------------------------------------------
    if top.dipos:
      if s.roundpipe(top.cdipoid[0],top.cdipozs[0],top.cdipoze[0],top.dipoap,
                     top.dipoox,top.dipooy,s.dipocm):
        return
    #--------------------------------------------------------------------------
    if top.sexts:
      if s.roundpipe(top.csextid[0],top.csextzs[0],top.csextze[0],top.sextap,
                     top.sextox,top.sextoy,s.sextcm):
        return
    #--------------------------------------------------------------------------
    # --- Drifts are checked for last in case the other elements might extend
    # --- into the neighboring drifts.
    if top.drfts:
      if s.roundpipe(top.cdrftid[0],top.cdrftzs[0],top.cdrftze[0],top.drftap,
                     top.drftox,top.drftoy,s.drftcm):
        return
    # --- If this part of the code is reached, then there are no applicable
    # --- boundaries, so turn the capacity matrix field solver off.
    top.fstype = 0
  #----------------------------------------------------------------------------
  def plotcond(s,plotphi=1,filled=0,plotedge=1,plotsym=1,ccolor='red',
               ecolor='green'):
    """
Makes a plot of the conductor.
  - plotphi=1 when true, plots contours of phi
  - filled=0 when true, plot filled contours
  - plotsym=1 when true, plots appropriate reflections when symmeties are used
  - plotedge=1 when true, plots a square around the edge of the mesh (with
               appropriate reflections when plotsym is true)
    """
    # --- Plot phi first (so filled contours are on bottom)
    if plotphi:
      plotc(transpose(getphi(iz=0)),w3d.ymesh,w3d.xmesh,filled=filled)
    if plotsym:
      # --- Plot negative y
      if w3d.l2symtry or w3d.l4symtry:
        if plotphi:
          plotc(transpose(getphi(iz=0)),-w3d.ymesh,w3d.xmesh,filled=filled)
      # --- Plot negative x
      if w3d.l4symtry:
        if plotphi:
          plotc(transpose(getphi(iz=0)),w3d.ymesh,-w3d.xmesh,filled=filled)
        if plotphi:
          plotc(transpose(getphi(iz=0)),-w3d.ymesh,-w3d.xmesh,filled=filled)
    # --- Plot conductor points next.
    if fxy.ncxy > 0 and top.fstype == 1:
      plp(fxy.ycond[:fxy.ncxy],fxy.xcond[:fxy.ncxy],color=ccolor,msize=2.)
      if plotsym:
        # --- Plot negative y
        if w3d.l2symtry or w3d.l4symtry:
          plp(-fxy.ycond[:fxy.ncxy],fxy.xcond[:fxy.ncxy],color=ccolor,msize=2.)
        # --- Plot negative x
        if w3d.l4symtry:
          plp(fxy.ycond[:fxy.ncxy],-fxy.xcond[:fxy.ncxy],color=ccolor,msize=2.)
          plp(-fxy.ycond[:fxy.ncxy],-fxy.xcond[:fxy.ncxy],color=ccolor,msize=2.)
    # --- Plot the edge of the mesh last
    if plotedge:
      if w3d.l4symtry and plotsym:
        plg([-w3d.xmmax, w3d.xmmax,w3d.xmmax,-w3d.xmmax,-w3d.xmmax],
            [-w3d.ymmax,-w3d.ymmax,w3d.ymmax, w3d.ymmax,-w3d.ymmax],
            color=ecolor)
      elif w3d.l2symtry and plotsym:
        plg([ w3d.xmmin, w3d.xmmax,w3d.xmmax, w3d.xmmin, w3d.xmmin],
            [-w3d.ymmax,-w3d.ymmax,w3d.ymmax, w3d.ymmax,-w3d.ymmax],
            color=ecolor)
      else:
        plg([ w3d.xmmin, w3d.xmmax,w3d.xmmax, w3d.xmmin, w3d.xmmin],
            [ w3d.ymmin, w3d.ymmin,w3d.ymmax, w3d.ymmax, w3d.ymmin],
            color=ecolor)




##############################################################################
def saveRealBoundary(object,filename):
  """
Save a class to the file filename.
  - object
  - filename
  """
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
  ff = open(filename,'r')
  result = cPickle.load(ff)
  result.enable()
  ff.close()
  return result

##############################################################################
