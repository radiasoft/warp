"""
The module supplies functions for dealing with particles.

These commands returns particle info based on selection criteria.
selectparticles(): return list of indices of particles selected
getn(): get number of particles selected
getx(), gety(), getz(), getr(), gettheta(): get particle position
getvx(), getvy(), getvz(): get particle velocity
getux(), getuy(), getuz(): get particle momentum/mass
getxp(), getyp(), getrp(): get tranverse normalized velocities
getgaminv(): get gamma inverse

getxxpslope(), getyypslope(): get x-x' and y-y' slopes

addparticles(): add particles to the simulation
"""
from warp import *
particles_version = "$Id: particles.py,v 1.2 2002/05/08 14:00:57 dave Exp $"

#-------------------------------------------------------------------------
def particlesdoc():
  import particles
  print particles.__doc__

##########################################################################
#-------------------------------------------------------------------------
# This returns the indices of the particles selected.
def selectparticles(iw=0,kwdict={},**kw):
  """
Selects particles based on either subsets or windows. By default it selects
from window 0, getting all of the live partilces (whose uzp != 0).
  - iw=0: Window to chose from
  - js=0: Species to chose from
  - jslist=None: List of Species to choose from, e.g. [0,3,4]; -1 for all specs
  - win=top.zwindows+top.zbeam: Windows to use (in lab frame)
  - z=top.zp: Coordinate for range selection
  - ix=-1: When 0 <= ix <= nx, picks particles within xmesh[ix]+-wx*dx
  - wx=1.: Width of window around xmesh[ix]
  - iy=-1: When 0 <= iy <= ny, picks particles within ymesh[iy]+-wy*dy
  - wy=1.: Width of window around ymesh[iy]
  - iz=-1: When 0 <= iz <= nz, picks particles within zmesh[iz]+-wz*dz
  - wz=1.: Width of window around zmesh[iz]
  - zl=None: When specified, lower range of selection region
  - zu=None: When specified, upper range of selection region
  """
  # --- Complete dictionary of possible keywords and their default values
  kwdefaults = {"js":0,"jslist":None,"win":None,"z":None,
                "ix":None,"wx":1.,"iy":None,"wy":1.,"iz":None,"wz":1.,
                "zl":None,"zu":None,'checkargs':0,'allowbadargs':0}

  # --- Create dictionary of local values and copy it into local dictionary,
  # --- ignoring keywords not listed in kwdefaults.
  kwvalues = kwdefaults.copy()
  kwvalues.update(kw)
  kwvalues.update(kwdict)
  for arg in kwdefaults.keys(): exec(arg+" = kwvalues['"+arg+"']")

  # --- Check the argument list for bad arguments.
  # --- 'checkargs' allows this routine to be called only to check the
  # --- input for bad arguments.
  # --- 'allowbadargs' allows this routine to be called with bad arguments.
  # --- These are intentionally undocumented features.

  kwvaluescopy = kwvalues.copy()
  for i in kwvaluescopy.keys():
    if i in kwdefaults.keys(): del kwvaluescopy[i]
  badargs = kwvaluescopy
  if checkargs: return badargs
  if badargs and not allowbadargs:
    raise "bad argument ",string.join(badargs.keys())

  # --- If jslist defined, call selectparticles repeatedly for each species
  # --- on the list
  del kwvalues['jslist']    # Remove list so selectparticles is subsequently
                            # called with one species at a time
  if jslist is not None:
    if jslist == -1:    jslist = range(0,top.ns)
    partlist = array([])
    for js in jslist:
        kwvalues['js'] = js
        newparts = selectparticles(iw, kwvalues)
        partlist = array(list(partlist)+list(newparts))
    return partlist

  ir1 = top.ins[js]-1
  ir2 = top.ins[js]+top.nps[js]-1
  if ir2 <= ir1: return array([])
  if zl is not None or zu is not None:
    if z is None: z = top.zp
    if zl is None: zl = -top.largepos
    if zu is None: zu = +top.largepos
    if zl > zu: print "Warning: zl > zu"
    ii=compress(logical_and(less(zl,z[ir1:ir2]),less(z[ir1:ir2],zu)),
                arrayrange(ir1,ir2))
  elif ix is not None:
    xl = w3d.xmmin + ix*w3d.dx - wx*w3d.dx
    xu = w3d.xmmin + ix*w3d.dx + wx*w3d.dx
    ii=compress(logical_and(less(xl,top.xp[ir1:ir2]),less(top.xp[ir1:ir2],xu)),
                arrayrange(ir1,ir2))
  elif iy is not None:
    yl = w3d.ymmin + iy*w3d.dy - wy*w3d.dy
    yu = w3d.ymmin + iy*w3d.dy + wy*w3d.dy
    ii=compress(logical_and(less(yl,top.yp[ir1:ir2]),less(top.yp[ir1:ir2],yu)),
                arrayrange(ir1,ir2))
  elif iz is not None:
    z = top.zp
    if lparallel:
      zl = top.zmslmin[0] + iz*w3d.dz - wz*w3d.dz + top.zbeam
      zu = top.zmslmin[0] + iz*w3d.dz + wz*w3d.dz + top.zbeam
    else:
      zl = w3d.zmmin + iz*w3d.dz - wz*w3d.dz + top.zbeam
      zu = w3d.zmmin + iz*w3d.dz + wz*w3d.dz + top.zbeam
    ii=compress(logical_and(less(zl,z[ir1:ir2]),less(z[ir1:ir2],zu)),
                arrayrange(ir1,ir2))
  elif iw < 0:
    if psubset==[]: setup_subsets()
    if -iw > len(psubset): raise "Bad window number"
    ii = ir1 + compress(less(psubset[-iw-1],top.nps[js]),psubset[-iw-1])
  elif iw == 0:
    ii = xrange(ir1,ir2)
  else:
    if win is None: win = top.zwindows[:,iw] + top.zbeam
    if len(shape(win)) == 2: win = win[:,iw]
    if z is None: z = top.zp
    ii=compress(logical_and(less(win[0],z[ir1:ir2]),less(z[ir1:ir2],win[1])),
                arrayrange(ir1,ir2))
  ii = compress(not_equal(take(top.uzp,ii),0.),ii)
  return ii

#-------------------------------------------------------------------------
# The following return a specific coordinate of the selected particles
# More documetation added after they are declared.
#-------------------------------------------------------------------------
def getn(iw=0,gather=1,**kw):
  "Returns number of particles in selection."
  ii = selectparticles(iw=iw,kwdict=kw)
  if lparallel and gather: return globalsum(len(ii))
  else: return len(ii)
#-------------------------------------------------------------------------
def getx(iw=0,gather=1,**kw):
  "Returns the X positions."
  ii = selectparticles(iw=iw,kwdict=kw)
  result = take(top.xp,ii)
  if lparallel and gather: return gatherarray(result)
  else: return result
#-------------------------------------------------------------------------
def gety(iw=0,gather=1,**kw):
  "Returns the Y positions."
  ii = selectparticles(iw=iw,kwdict=kw)
  result = take(top.yp,ii)
  if lparallel and gather: return gatherarray(result)
  else: return result
#-------------------------------------------------------------------------
def getz(iw=0,gather=1,**kw):
  "Returns the Z positions."
  ii = selectparticles(iw=iw,kwdict=kw)
  result = take(top.zp,ii)
  if lparallel and gather: return gatherarray(result)
  else: return result
#-------------------------------------------------------------------------
def getr(iw=0,gather=1,**kw):
  "Returns the R postions."
  ii = selectparticles(iw=iw,kwdict=kw)
  result = sqrt(take(top.xp,ii)**2 + take(top.yp,ii)**2)
  if lparallel and gather: return gatherarray(result)
  else: return result
#-------------------------------------------------------------------------
def gettheta(iw=0,gather=1,**kw):
  "Returns the theta postions."
  ii = selectparticles(iw=iw,kwdict=kw)
  result = arctan2(take(top.yp,ii),take(top.xp,ii))
  if lparallel and gather: return gatherarray(result)
  else: return result
#-------------------------------------------------------------------------
def getvx(iw=0,gather=1,**kw):
  "Returns the X velocity."
  ii = selectparticles(iw=iw,kwdict=kw)
  result = take(top.uxp*top.gaminv,ii)
  if lparallel and gather: return gatherarray(result)
  else: return result
#-------------------------------------------------------------------------
def getvy(iw=0,gather=1,**kw):
  "Returns the Y velocity."
  ii = selectparticles(iw=iw,kwdict=kw)
  result = take(top.uyp*top.gaminv,ii)
  if lparallel and gather: return gatherarray(result)
  else: return result
#-------------------------------------------------------------------------
def getvz(iw=0,gather=1,**kw):
  "Returns the Z velocity."
  ii = selectparticles(iw=iw,kwdict=kw)
  result = take(top.uzp*top.gaminv,ii)
  if lparallel and gather: return gatherarray(result)
  else: return result
#-------------------------------------------------------------------------
def getux(iw=0,gather=1,**kw):
  "Returns the X momentum over mass."
  ii = selectparticles(iw=iw,kwdict=kw)
  result = take(top.uxp,ii)
  if lparallel and gather: return gatherarray(result)
  else: return result
#-------------------------------------------------------------------------
def getuy(iw=0,gather=1,**kw):
  "Returns the Y momentum over mass."
  ii = selectparticles(iw=iw,kwdict=kw)
  result = take(top.uyp,ii)
  if lparallel and gather: return gatherarray(result)
  else: return result
#-------------------------------------------------------------------------
def getuz(iw=0,gather=1,**kw):
  "Returns the Z momentum over mass."
  ii = selectparticles(iw=iw,kwdict=kw)
  result = take(top.uzp,ii)
  if lparallel and gather: return gatherarray(result)
  else: return result
#-------------------------------------------------------------------------
def getxp(iw=0,gather=1,**kw):
  "Returns the X velocity over the Z velocity (X')."
  ii = selectparticles(iw=iw,kwdict=kw)
  result = take(top.uxp,ii)/take(top.uzp,ii)
  if lparallel and gather: return gatherarray(result)
  else: return result
#-------------------------------------------------------------------------
def getyp(iw=0,gather=1,**kw):
  "Returns the Y velocity over the Z velocity (Y')."
  ii = selectparticles(iw=iw,kwdict=kw)
  result = take(top.uyp,ii)/take(top.uzp,ii)
  if lparallel and gather: return gatherarray(result)
  else: return result
#-------------------------------------------------------------------------
def getrp(iw=0,gather=1,**kw):
  "Returns the radial velocity over the Z velocity (R')."
  ii = selectparticles(iw=iw,kwdict=kw)
  tt = arctan2(take(top.yp,ii),take(top.xp,ii))
  result = (take(top.uxp,ii)*cos(tt)+take(top.uyp,ii)*sin(tt))/take(top.uzp,ii)
  if lparallel and gather: return gatherarray(result)
  else: return result
#-------------------------------------------------------------------------
def getgaminv(iw=0,gather=1,**kw):
  "Returns the gamma inverse."
  ii = selectparticles(iw=iw,kwdict=kw)
  result = take(top.gaminv,ii)
  if lparallel and gather: return gatherarray(result)
  else: return result
#-------------------------------------------------------------------------
# Add the selectparticles documentation to each of the routines.
if sys.version[:5] != "1.5.1":
  if lparallel:
    _gatherdoc = "  gather=1 When 1, all data is gathered to PE0"
  else:
    _gatherdoc = ""
  getn.__doc__ = getn.__doc__ + selectparticles.__doc__ + _gatherdoc
  getx.__doc__ = getx.__doc__ + selectparticles.__doc__ + _gatherdoc
  gety.__doc__ = gety.__doc__ + selectparticles.__doc__ + _gatherdoc
  getz.__doc__ = getz.__doc__ + selectparticles.__doc__ + _gatherdoc
  getr.__doc__ = getr.__doc__ + selectparticles.__doc__ + _gatherdoc
  gettheta.__doc__ = gettheta.__doc__ + selectparticles.__doc__ + _gatherdoc
  getvx.__doc__ = getvx.__doc__ + selectparticles.__doc__ + _gatherdoc
  getvy.__doc__ = getvy.__doc__ + selectparticles.__doc__ + _gatherdoc
  getvz.__doc__ = getvz.__doc__ + selectparticles.__doc__ + _gatherdoc
  getux.__doc__ = getux.__doc__ + selectparticles.__doc__ + _gatherdoc
  getuy.__doc__ = getuy.__doc__ + selectparticles.__doc__ + _gatherdoc
  getuz.__doc__ = getuz.__doc__ + selectparticles.__doc__ + _gatherdoc
  getxp.__doc__ = getxp.__doc__ + selectparticles.__doc__ + _gatherdoc
  getyp.__doc__ = getyp.__doc__ + selectparticles.__doc__ + _gatherdoc
  getrp.__doc__ = getrp.__doc__ + selectparticles.__doc__ + _gatherdoc
#-------------------------------------------------------------------------

##########################################################################
def getxxpslope(iw=0,iz=-1):
  """
Calculates the x-x' slope based on either the window moments in window iw
or the zmoments at iz. This returns a tuple containing (slope,offset,vz).
The product slope*vz gives the slope for x-vx.
  """
  if not lparallel:
    if 0 <= iz <= w3d.nz:
      slope = (top.xxpbarz[iz]-top.xbarz[iz]*top.xpbarz[iz])/top.xrmsz[iz]**2
      offset = top.xpbarz[iz]-slope*top.xbarz[iz]
      vz = top.vzbarz[iz]
    else:
      iiw = max(0,iw)
      slope = (top.xxpbar[iiw]-top.xbar[iiw]*top.xpbar[iiw])/top.xrms[iiw]**2
      offset = top.xpbar[iiw]-slope*top.xbar[iiw]
      vz = top.vzbar[iiw]
  else:
    if 0 <= iz <= w3d.nzfull:
      pe = convertizptope(iz)
      if me == pe:
        iz = iz - top.izpslave[me]
        slope = (top.xxpbarz[iz]-top.xbarz[iz]*top.xpbarz[iz])/top.xrmsz[iz]**2
        offset = top.xpbarz[iz]-slope*top.xbarz[iz]
        vz = top.vzbarz[iz]
      else:
        (slope,offset,vz) = (0.,0.,0.)
      (slope,offset,vz) = tuple(broadcast(array([slope,offset,vz]),pe))
    else:
      iiw = max(0,iw)
      pe = convertiwtope(iiw)
      if me == pe:
        slope = (top.xxpbar[iiw]-top.xbar[iiw]*top.xpbar[iiw])/top.xrms[iiw]**2
        offset = top.xpbar[iiw]-slope*top.xbar[iiw]
        vz = top.vzbar[iiw]
      else:
        (slope,offset,vz) = (0.,0.,0.)
      (slope,offset,vz) = tuple(broadcast(array([slope,offset,vz]),pe))
  return (slope,offset,vz)
#-------------------------------------------------------------------------
def getyypslope(iw=0,iz=-1):
  """
Calculates the y-y' slope based on either the window moments in window iw
or the zmoments at iz. This returns a tuple containing (slope,offset,vz).
The product slope*vz gives the slope for y-vy.
  """
  if not lparallel:
    if 0 <= iz <= w3d.nz:
      slope = (top.yypbarz[iz]-top.ybarz[iz]*top.ypbarz[iz])/top.yrmsz[iz]**2
      offset = top.ypbarz[iz]-slope*top.ybarz[iz]
      vz = top.vzbarz[iz]
    else:
      iiw = max(0,iw)
      slope = (top.yypbar[iiw]-top.ybar[iiw]*top.ypbar[iiw])/top.yrms[iiw]**2
      offset = top.ypbar[iiw]-slope*top.ybar[iiw]
      vz = top.vzbar[iiw]
  else:
    if 0 <= iz <= w3d.nzfull:
      pe = convertizptope(iz)
      if me == pe:
        iz = iz - top.izpslave[me]
        slope = (top.yypbarz[iz]-top.ybarz[iz]*top.ypbarz[iz])/top.yrmsz[iz]**2
        offset = top.ypbarz[iz]-slope*top.ybarz[iz]
        vz = top.vzbarz[iz]
      else:
        (slope,offset,vz) = (0.,0.,0.)
      (slope,offset,vz) = tuple(broadcast(array([slope,offset,vz]),pe))
    else:
      iiw = max(0,iw)
      pe = convertiwtope(iiw)
      if me == pe:
        slope = (top.yypbar[iiw]-top.ybar[iiw]*top.ypbar[iiw])/top.yrms[iiw]**2
        offset = top.ypbar[iiw]-slope*top.ybar[iiw]
        vz = top.vzbar[iiw]
      else:
        (slope,offset,vz) = (0.,0.,0.)
      (slope,offset,vz) = tuple(broadcast(array([slope,offset,vz]),pe))
  return (slope,offset,vz)
#-------------------------------------------------------------------------
def getvzrange():
  "Returns a tuple containg the Vz range for plots"
  if (top.vzrng != 0.):
     vzmax = (1. + top.vtilt)*top.vbeam*(1.+top.vzrng) - top.vzshift
     vzmin = (1. - top.vtilt)*top.vbeam*(1.-top.vzrng) - top.vzshift
  else:
     vzmax = top.vzmaxp + 0.1*(top.vzmaxp-top.vzminp)
     vzmin = top.vzminp - 0.1*(top.vzmaxp-top.vzminp)
  return (vzmin,vzmax)


#-------------------------------------------------------------------------
#-------------------------------------------------------------------------
def addparticles(x=0.,y=0.,z=0.,vx=0.,vy=0.,vz=0.,gi=1.,js=0,
                 lallindomain=false,zmmin=None,zmmax=None,lmomentum=false,
                 resetrho=false,dofieldsol=false,resetmoments=false):
  """
Adds particles to the simulation
  x,y,z,vx,vy,vz,gi: particle coordinates and velocities. Can be a arrays or
                     scalars. Scalars are broadcast to all particles. Any
                     that are unsupplied default to zero, except gi,
                     which defaults to 1.
  js=0: species to which new particles belong
  lallindomain=false: When true, all particles are assumed to be with in the
                      z extent of the domain so particle scraping is not done.
                      Note that no check is done transversely.
  zmmin=w3d.zmmin,zmmax=w3d.zmmax: z extent of the domain
  lmomentum=false: Set to false when velocities are input as velocities, true
                   when input as massless momentum (as WARP stores them).
                   Only used when top.lrelativ is true.
  """

  # --- Get length of arrays, set to one for scalars
  try:              lenx = len(x)
  except TypeError: lenx = 1
  try:              leny = len(y)
  except TypeError: leny = 1
  try:              lenz = len(z)
  except TypeError: lenz = 1
  try:              lenvx = len(vx)
  except TypeError: lenvx = 1
  try:              lenvy = len(vy)
  except TypeError: lenvy = 1
  try:              lenvz = len(vz)
  except TypeError: lenvz = 1
  try:              lengi = len(gi)
  except TypeError: lengi = 1

  # --- Max length of input arrays
  maxlen = max(lenx,leny,lenz,lenvx,lenvy,lenvz,lengi)
  assert lenx==maxlen or lenx==1,"Length of x doesn't match len of others"
  assert leny==maxlen or leny==1,"Length of y doesn't match len of others"
  assert lenz==maxlen or lenz==1,"Length of z doesn't match len of others"
  assert lenvx==maxlen or lenvx==1,"Length of vx doesn't match len of others"
  assert lenvy==maxlen or lenvy==1,"Length of vy doesn't match len of others"
  assert lenvz==maxlen or lenvz==1,"Length of vz doesn't match len of others"
  assert lengi==maxlen or lengi==1,"Length of gi doesn't match len of others"

  # --- Convert all to arrays of length maxlen, broadcasting scalars
  x = array(x)*ones(maxlen,'d')
  y = array(y)*ones(maxlen,'d')
  z = array(z)*ones(maxlen,'d')
  vx = array(vx)*ones(maxlen,'d')
  vy = array(vy)*ones(maxlen,'d')
  vz = array(vz)*ones(maxlen,'d')
  gi = array(gi)*ones(maxlen,'d')

  # --- Set extent of domain
  if not lparallel:
    if zmmin is None: zmmin = w3d.zmmin
    if zmmax is None: zmmax = w3d.zmmax
  else:
    if zmmin is None: zmmin = top.zpslmin[me]
    if zmmax is None: zmmax = top.zpslmax[me]

  # --- Now data can be passed into the fortran addparticles routine.
  addpart(js+1,maxlen,x,y,z,vx,vy,vz,gi,lallindomain,zmmin,zmmax,lmomentum)

  # --- If the slice code is active, then call initdtp
  if package()[0] == 'wxy': initdtp()

  # --- Do followup work if requested
  if resetrho:
    w3d.rho=0.0
    loadrho()
  if dofieldsol:
    fieldsol(-1)
  if resetmoments:
    import getzmom
    getzmom.zmmnt()
    if top.it%top.nhist == 0:
      top.jhist = top.jhist - 1
      savehist(top.time)

































































