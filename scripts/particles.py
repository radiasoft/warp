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
getpid(): get particle identification number

getxxpslope(), getyypslope(): get x-x' and y-y' slopes

addparticles(): add particles to the simulation

setup_subsets(): Create subsets for particle plots (negative window numbers)
clear_subsets(): Clears the subsets for particle plots (negative window
numbers)
"""
from warp import *
import random
particles_version = "$Id: particles.py,v 1.25 2005/04/05 19:00:04 dave Exp $"

#-------------------------------------------------------------------------
def particlesdoc():
  import particles
  print particles.__doc__

##########################################################################
# Setup the random subsets. This is only called when the first plot
# is made so that top.npsplt is known to be set.
# These routines can be called at any time by the user to add more subsets,
# for example subsets with more or fewer particles, or subsets based on
# the number of particles in a species other than 0.
psubset=[]
def setup_subsets(js=0):
  """
Adds plotting subset to the list
  - js=0 is the species to create a subset for
  """
  global psubset
  if lparallel:
    totalnp = parallelsum(top.nps[js])
    if totalnp == 0: totalnp = 1
    fracnp = float(top.nps[js])/float(totalnp)
  else:
    fracnp = 1.
  for i in xrange(0,len(top.npplot)):
    ntopick=min(top.nps[js],int(top.npplot[i]*fracnp+0.5))
    ii = arrayrange(top.nps[0])
    rr = top.nps[0]*RandomArray.random(top.nps[0])
    ii = compress(less(rr,ntopick),ii)
    psubset.append(ii.astype('i'))
#----------------------------------------------------------------------------
def clear_subsets():
  "Clears the particle subsets so that they can be updated."
  global psubset
  psubset = []

#----------------------------------------------------------------------------
# Modified version of the sample routine from the random module of Python2.3
def populationsample(population,k,self=random.Random(0)):
    """Chooses k unique random elements from a population sequence.

    Returns a new list containing elements from the population while
    leaving the original population unchanged.  The resulting list is
    in selection order so that all sub-slices will also be valid random
    samples.  This allows raffle winners (the sample) to be partitioned
    into grand prize and second place winners (the subslices).

    Members of the population need not be hashable or unique.  If the
    population contains repeats, then each occurrence is a possible
    selection in the sample.

    To choose a sample in a range of integers, use xrange as an argument.
    This is especially fast and space efficient for sampling from a
    large population:   sample(xrange(10000000), 60)
    """

    # Sampling without replacement entails tracking either potential
    # selections (the pool) in a list or previous selections in a
    # dictionary.

    # When the number of selections is small compared to the population,
    # then tracking selections is efficient, requiring only a small
    # dictionary and an occasional reselection.  For a larger number of
    # selections, the pool tracking method is preferred since the list takes
    # less space than the dictionary and it doesn't suffer from frequent
    # reselections.

    n = len(population)
    if not 0 <= k <= n:
        raise ValueError, "sample larger than population"
    random = self.random
    _int = int
    # --- Result is an array rather than a list ---
    typecode = array([population[0]]).typecode()
    result = zeros(k,typecode=typecode)
    if n < 6 * k:     # if n len list takes less space than a k len dict
        pool = list(population)
        for i in xrange(k):         # invariant:  non-selected at [0,n-i)
            j = _int(random() * (n-i))
            result[i] = pool[j]
            pool[j] = pool[n-i-1]   # move non-selected item into vacancy
    else:
        try:
            n > 0 and (population[0], population[n//2], population[n-1])
        except (TypeError, KeyError):   # handle sets and dictionaries
            population = tuple(population)
        selected = {}
        for i in xrange(k):
            j = _int(random() * n)
            while j in selected:
                j = _int(random() * n)
            result[i] = selected[j] = population[j]
    return result

# --- Old method which replicates the numbers selected by the fortran
# --- routine psubsets. Note that that method breaks down when
# --- nps/inclump[i] > npplot[i], when nps/inclump[i] particles will be plotted.
# --- This version is maintained in case a user wants a close comparison with
# --- the Basis version.
def setup_subsetsold(js=0):
  """Old subset calculator, do not use"""
  global psubset
  # --- Print warning if npsplt is zero, in which case the subsets won't work
  if (top.npsplt == 0):
    remark("WARNING: npsplt is zero, subsets not calculated")
    return
  d2 = len(top.npplot)
  # --- Create temp ii array and copy top.isubset into it. Then replicate
  # --- the data throughout ii.
  for i in xrange(0,d2):
    n = sum(top.isubset[:,i])
    nsets = int(min(top.np_s[js],top.npmax)/top.npsplt+1)
    ii = zeros(n*nsets,'i') + top.npmax
    ii[0:n] = nonzero(top.isubset[:,i])
    for j in xrange(1,nsets):
      ii[j*n:j*n+n] = ii[0:n] + j*top.npsplt
    ii = compress(less(ii,top.npmax),ii)
    psubset.append(ii)

##########################################################################
#-------------------------------------------------------------------------
def getattrwithsuffix(object=None,name=None,suffix='',checkmain=1,
                      default=None):
  """
Helper function modeled on getattr. Like getattr, given an object and a name,
it will return object.name. Optionally, a suffix can be given, in which case,
it will return the equivalent of getattr(object,name+suffix). The object can
be None, in which case the main dictionary will be searched,
returning __main__.__dict__[name+suffix]. Also, object can be a dictionary
itself. Finally, if the object is specified but name+suffix is not found and
if checkmain is true, then the main dictionary will again be searched. If
nothing if found and default is specified, it is returned, otherwise an
exception is raised.
  """
  assert name is not None,"A name must be supplied"
  if object is not None:
    try:
      result = getattr(object,name+suffix)
      return result
    except AttributeError:
      pass
  if type(object) is DictionaryType:
    try:
      result = object[name+suffix]
      return result
    except KeyError:
      pass
  if object is None or checkmain:
    import __main__
    try:
      result = __main__.__dict__[name+suffix]
      return result
    except KeyError:
      pass

  if default is not None:
    return default
  else:
    raise AttributeError

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
  - zc=None: When specified, picks particles within zc+-wz*dz
  - xc=None: When specified, picks particles within xc+-wx*dx
  - yc=None: When specified, picks particles within yc+-wy*dy
  - ii=None: If ii is supplied, it is just returned.
  - lost=false: When true, returns indices to the lost particles rather than
                the live particles
  """
  # --- Complete dictionary of possible keywords and their default values
  kwdefaults = {"js":0,"jslist":None,"win":None,"z":None,
                "ix":None,"wx":1.,"iy":None,"wy":1.,"iz":None,"wz":1.,
                "zl":None,"zu":None,"zc":None,"xc":None,"yc":None,"ii":None,
                "lost":false,
                'checkargs':0,'allowbadargs':0}

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

  # --- If ii is defined, then just return it. This allows the 'get' routines
  # --- to do the take operation with a pre-defined ii.
  if ii is not None: return ii

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

  # --- If lost is true, append the suffix lost to the variable names
  if lost: suffix = 'lost'
  else:    suffix = ''

  ins = getattrwithsuffix(top,'ins',suffix)
  nps = getattrwithsuffix(top,'nps',suffix)

  ir1 = ins[js] - 1
  ir2 = ins[js] + nps[js] - 1

  if ir2 <= ir1: return array([])

  if zl is not None or zu is not None:
    if z is None: z = getattrwithsuffix(top,'zp',suffix)
    if zl is None: zl = -top.largepos
    if zu is None: zu = +top.largepos
    if zl > zu: print "Warning: zl > zu"
    ii=compress(logical_and(less(zl,z[ir1:ir2]),less(z[ir1:ir2],zu)),
                arrayrange(ir1,ir2))
  elif ix is not None:
    xl = w3d.xmmin + ix*w3d.dx - wx*w3d.dx
    xu = w3d.xmmin + ix*w3d.dx + wx*w3d.dx
    x = getattrwithsuffix(top,'xp',suffix)
    ii=compress(logical_and(less(xl,x[ir1:ir2]),less(x[ir1:ir2],xu)),
                arrayrange(ir1,ir2))
  elif iy is not None:
    yl = w3d.ymmin + iy*w3d.dy - wy*w3d.dy
    yu = w3d.ymmin + iy*w3d.dy + wy*w3d.dy
    y = getattrwithsuffix(top,'yp',suffix)
    ii=compress(logical_and(less(yl,y[ir1:ir2]),less(y[ir1:ir2],yu)),
                arrayrange(ir1,ir2))
  elif iz is not None:
    z = getattrwithsuffix(top,'zp',suffix)
    zl = w3d.zmminglobal + iz*w3d.dz - wz*w3d.dz + top.zbeam
    zu = w3d.zmminglobal + iz*w3d.dz + wz*w3d.dz + top.zbeam
    ii=compress(logical_and(less(zl,z[ir1:ir2]),less(z[ir1:ir2],zu)),
                arrayrange(ir1,ir2))
  elif zc is not None:
    z = getattrwithsuffix(top,'zp',suffix)
    zl = zc - wz*w3d.dz
    zu = zc + wz*w3d.dz
    ii=compress(logical_and(less(zl,z[ir1:ir2]),less(z[ir1:ir2],zu)),
                arrayrange(ir1,ir2))
  elif xc is not None:
    x = getattrwithsuffix(top,'xp',suffix)
    xl = xc - wx*w3d.dx
    xu = xc + wx*w3d.dx
    ii=compress(logical_and(less(xl,x[ir1:ir2]),less(x[ir1:ir2],xu)),
                arrayrange(ir1,ir2))
  elif yc is not None:
    y = getattrwithsuffix(top,'yp',suffix)
    yl = yc - wy*w3d.dy
    yu = yc + wy*w3d.dy
    ii=compress(logical_and(less(yl,y[ir1:ir2]),less(y[ir1:ir2],yu)),
                arrayrange(ir1,ir2))
  elif iw < 0:
    # --- Add some smarts about choosing which method to use.
    # --- If the number of particles is increasing, then use the
    # --- population sampling. Otherwise use the preset subsets.
    # --- Also, if the number of species has changed, then also
    # --- use the population sampling.
    try: selectparticles.npsprev
    except: selectparticles.npsprev = nps + 0
    if len(selectparticles.npsprev) != top.ns:
      selectparticles.npsprev = zeros(top.ns)
    if selectparticles.npsprev[js] >= nps[js]:
      if psubset==[]: setup_subsets()
      if -iw > len(psubset): raise "Bad window number"
      ii = ir1 + compress(less(psubset[-iw-1],nps[js]),psubset[-iw-1])
    else:
      # --- Once this method is used, always use it.
      selectparticles.npsprev = zeros(top.ns)
      if nps[js] <= top.npplot[-iw-1]:
        ii = ir1 + arange(nps[js])
      else:
        ii = ir1 + populationsample(xrange(nps[js]),top.npplot[-iw-1])
  elif iw == 0:
    ii = xrange(ir1,ir2)
  else:
    if win is None: win = top.zwindows[:,iw] + top.zbeam
    if len(shape(win)) == 2: win = win[:,iw]
    if z is None: z = getattrwithsuffix(top,'zp',suffix)
    ii=compress(logical_and(less(win[0],z[ir1:ir2]),less(z[ir1:ir2],win[1])),
                arrayrange(ir1,ir2))
  uz = getattrwithsuffix(top,'uzp',suffix)
  ii = compress(not_equal(take(uz,ii),0.),ii)
  return ii

#-------------------------------------------------------------------------
# The following return a specific coordinate of the selected particles
# More documetation added after they are declared.
# --- The checks of the len of ii are there in case the particle arrays
# --- are unallocated.
#-------------------------------------------------------------------------
def getn(iw=0,gather=1,bcast=0,**kw):
  "Returns number of particles in selection."
  ii = selectparticles(iw=iw,kwdict=kw)
  if lparallel and gather: return globalsum(len(ii))
  else: return len(ii)
#-------------------------------------------------------------------------
def getx(iw=0,gather=1,bcast=0,**kw):
  "Returns the X positions."
  ii = selectparticles(iw=iw,kwdict=kw)
  suffix = (kw.get('lost',0) and 'lost') or ''
  if len(ii) > 0:
    x = getattrwithsuffix(top,'xp',suffix)
    result = take(x,ii)
  else:
    result = array([],'d')
  if lparallel and gather: return gatherarray(result,bcast=bcast)
  else: return result
#-------------------------------------------------------------------------
def gety(iw=0,gather=1,bcast=0,**kw):
  "Returns the Y positions."
  ii = selectparticles(iw=iw,kwdict=kw)
  suffix = (kw.get('lost',0) and 'lost') or ''
  if len(ii) > 0:
    y = getattrwithsuffix(top,'yp',suffix)
    result = take(y,ii)
  else:
    result = array([],'d')
  if lparallel and gather: return gatherarray(result,bcast=bcast)
  else: return result
#-------------------------------------------------------------------------
def getz(iw=0,gather=1,bcast=0,**kw):
  "Returns the Z positions."
  ii = selectparticles(iw=iw,kwdict=kw)
  suffix = (kw.get('lost',0) and 'lost') or ''
  if len(ii) > 0:
    z = getattrwithsuffix(top,'zp',suffix)
    result = take(z,ii)
  else:
    result = array([],'d')
  if lparallel and gather: return gatherarray(result,bcast=bcast)
  else: return result
#-------------------------------------------------------------------------
def getr(iw=0,gather=1,bcast=0,**kw):
  "Returns the R postions."
  ii = selectparticles(iw=iw,kwdict=kw)
  suffix = (kw.get('lost',0) and 'lost') or ''
  if len(ii) > 0:
    x = getattrwithsuffix(top,'xp',suffix)
    y = getattrwithsuffix(top,'yp',suffix)
    result = sqrt(take(x,ii)**2 + take(y,ii)**2)
  else:
    result = array([],'d')
  if lparallel and gather: return gatherarray(result,bcast=bcast)
  else: return result
#-------------------------------------------------------------------------
def gettheta(iw=0,gather=1,bcast=0,**kw):
  "Returns the theta postions."
  ii = selectparticles(iw=iw,kwdict=kw)
  suffix = (kw.get('lost',0) and 'lost') or ''
  if len(ii) > 0:
    x = getattrwithsuffix(top,'xp',suffix)
    y = getattrwithsuffix(top,'yp',suffix)
    result = arctan2(take(y,ii),take(x,ii))
  else:
    result = array([],'d')
  if lparallel and gather: return gatherarray(result,bcast=bcast)
  else: return result
#-------------------------------------------------------------------------
def getvx(iw=0,gather=1,bcast=0,**kw):
  "Returns the X velocity."
  ii = selectparticles(iw=iw,kwdict=kw)
  suffix = (kw.get('lost',0) and 'lost') or ''
  if len(ii) > 0:
    ux = getattrwithsuffix(top,'uxp',suffix)
    gaminv = getattrwithsuffix(top,'gaminv',suffix)
    result = take(ux,ii)*take(gaminv,ii)
  else:
    result = array([],'d')
  if lparallel and gather: return gatherarray(result,bcast=bcast)
  else: return result
#-------------------------------------------------------------------------
def getvy(iw=0,gather=1,bcast=0,**kw):
  "Returns the Y velocity."
  ii = selectparticles(iw=iw,kwdict=kw)
  suffix = (kw.get('lost',0) and 'lost') or ''
  if len(ii) > 0:
    uy = getattrwithsuffix(top,'uyp',suffix)
    gaminv = getattrwithsuffix(top,'gaminv',suffix)
    result = take(uy,ii)*take(gaminv,ii)
  else:
    result = array([],'d')
  if lparallel and gather: return gatherarray(result,bcast=bcast)
  else: return result
#-------------------------------------------------------------------------
def getvz(iw=0,gather=1,bcast=0,**kw):
  "Returns the Z velocity."
  ii = selectparticles(iw=iw,kwdict=kw)
  suffix = (kw.get('lost',0) and 'lost') or ''
  if len(ii) > 0:
    uz = getattrwithsuffix(top,'uzp',suffix)
    gaminv = getattrwithsuffix(top,'gaminv',suffix)
    result = take(uz,ii)*take(gaminv,ii)
  else:
    result = array([],'d')
  if lparallel and gather: return gatherarray(result,bcast=bcast)
  else: return result
#-------------------------------------------------------------------------
def getvr(iw=0,gather=1,bcast=0,**kw):
  "Returns the radial velocity."
  ii = selectparticles(iw=iw,kwdict=kw)
  suffix = (kw.get('lost',0) and 'lost') or ''
  if len(ii) > 0:
    x = getattrwithsuffix(top,'xp',suffix)
    y = getattrwithsuffix(top,'yp',suffix)
    ux = getattrwithsuffix(top,'uxp',suffix)
    uy = getattrwithsuffix(top,'uyp',suffix)
    tt = arctan2(take(y,ii),take(x,ii))
    result = take(ux,ii)*cos(tt) + take(uy,ii)*sin(tt)
  else:
    result = array([],'d')
  if lparallel and gather: return gatherarray(result,bcast=bcast)
  else: return result
#-------------------------------------------------------------------------
def getvtheta(iw=0,gather=1,bcast=0,**kw):
  "Returns the azimuthal velocity."
  ii = selectparticles(iw=iw,kwdict=kw)
  suffix = (kw.get('lost',0) and 'lost') or ''
  if len(ii) > 0:
    x = getattrwithsuffix(top,'xp',suffix)
    y = getattrwithsuffix(top,'yp',suffix)
    ux = getattrwithsuffix(top,'uxp',suffix)
    uy = getattrwithsuffix(top,'uyp',suffix)
    tt = arctan2(take(y,ii),take(x,ii))
    result = -take(ux,ii)*sin(tt) + take(uy,ii)*cos(tt)
  else:
    result = array([],'d')
  if lparallel and gather: return gatherarray(result,bcast=bcast)
  else: return result
#-------------------------------------------------------------------------
def getux(iw=0,gather=1,bcast=0,**kw):
  "Returns the X momentum over mass."
  ii = selectparticles(iw=iw,kwdict=kw)
  suffix = (kw.get('lost',0) and 'lost') or ''
  if len(ii) > 0:
    ux = getattrwithsuffix(top,'uxp',suffix)
    result = take(ux,ii)
  else:
    result = array([],'d')
  if lparallel and gather: return gatherarray(result,bcast=bcast)
  else: return result
#-------------------------------------------------------------------------
def getuy(iw=0,gather=1,bcast=0,**kw):
  "Returns the Y momentum over mass."
  ii = selectparticles(iw=iw,kwdict=kw)
  suffix = (kw.get('lost',0) and 'lost') or ''
  if len(ii) > 0:
    uy = getattrwithsuffix(top,'uyp',suffix)
    result = take(uy,ii)
  else:
    result = array([],'d')
  if lparallel and gather: return gatherarray(result,bcast=bcast)
  else: return result
#-------------------------------------------------------------------------
def getuz(iw=0,gather=1,bcast=0,**kw):
  "Returns the Z momentum over mass."
  ii = selectparticles(iw=iw,kwdict=kw)
  suffix = (kw.get('lost',0) and 'lost') or ''
  if len(ii) > 0:
    uz = getattrwithsuffix(top,'uzp',suffix)
    result = take(uz,ii)
  else:
    result = array([],'d')
  if lparallel and gather: return gatherarray(result,bcast=bcast)
  else: return result
#-------------------------------------------------------------------------
def getxp(iw=0,gather=1,bcast=0,**kw):
  "Returns the X velocity over the Z velocity (X')."
  ii = selectparticles(iw=iw,kwdict=kw)
  suffix = (kw.get('lost',0) and 'lost') or ''
  if len(ii) > 0:
    ux = getattrwithsuffix(top,'uxp',suffix)
    uz = getattrwithsuffix(top,'uzp',suffix)
    result = take(ux,ii)/take(uz,ii)
  else:
    result = array([],'d')
  if lparallel and gather: return gatherarray(result,bcast=bcast)
  else: return result
#-------------------------------------------------------------------------
def getyp(iw=0,gather=1,bcast=0,**kw):
  "Returns the Y velocity over the Z velocity (Y')."
  ii = selectparticles(iw=iw,kwdict=kw)
  suffix = (kw.get('lost',0) and 'lost') or ''
  if len(ii) > 0:
    uy = getattrwithsuffix(top,'uyp',suffix)
    uz = getattrwithsuffix(top,'uzp',suffix)
    result = take(uy,ii)/take(uz,ii)
  else:
    result = array([],'d')
  if lparallel and gather: return gatherarray(result,bcast=bcast)
  else: return result
#-------------------------------------------------------------------------
def getrp(iw=0,gather=1,bcast=0,**kw):
  "Returns the radial velocity over the Z velocity (R')."
  ii = selectparticles(iw=iw,kwdict=kw)
  suffix = (kw.get('lost',0) and 'lost') or ''
  if len(ii) > 0:
    x = getattrwithsuffix(top,'xp',suffix)
    y = getattrwithsuffix(top,'yp',suffix)
    ux = getattrwithsuffix(top,'uxp',suffix)
    uy = getattrwithsuffix(top,'uyp',suffix)
    uz = getattrwithsuffix(top,'uzp',suffix)
    tt = arctan2(take(y,ii),take(x,ii))
    result = ((take(ux,ii)*cos(tt)+take(uy,ii)*sin(tt))/
              take(uz,ii))
  else:
    result = array([],'d')
  if lparallel and gather: return gatherarray(result,bcast=bcast)
  else: return result
#-------------------------------------------------------------------------
def gettp(iw=0,gather=1,bcast=0,**kw):
  "Returns the azimuthal velocity over the Z velocity (R')."
  ii = selectparticles(iw=iw,kwdict=kw)
  suffix = (kw.get('lost',0) and 'lost') or ''
  if len(ii) > 0:
    x = getattrwithsuffix(top,'xp',suffix)
    y = getattrwithsuffix(top,'yp',suffix)
    ux = getattrwithsuffix(top,'uxp',suffix)
    uy = getattrwithsuffix(top,'uyp',suffix)
    uz = getattrwithsuffix(top,'uzp',suffix)
    tt = arctan2(take(y,ii),take(x,ii))
    result = ((-take(ux,ii)*sin(tt)+take(uy,ii)*cos(tt))/
              take(uz,ii))
  else:
    result = array([],'d')
  if lparallel and gather: return gatherarray(result,bcast=bcast)
  else: return result
#-------------------------------------------------------------------------
def getgaminv(iw=0,gather=1,bcast=0,**kw):
  "Returns the gamma inverse."
  ii = selectparticles(iw=iw,kwdict=kw)
  suffix = (kw.get('lost',0) and 'lost') or ''
  if len(ii) > 0:
    gaminv = getattrwithsuffix(top,'gaminv',suffix)
    result = take(gaminv,ii)
  else:
    result = array([],'d')
  if lparallel and gather: return gatherarray(result,bcast=bcast)
  else: return result
#-------------------------------------------------------------------------
def getpid(iw=0,id=0,gather=1,bcast=0,**kw):
  """Returns particle id number.
  -id=0: which pid value to return
  """
  lost = kw.get('lost',0)
  suffix = (kw.get('lost',0) and 'lost') or ''
  if lost: dopid = (top.npidlostmax > 0)
  else:    dopid = (top.npmaxi == top.npmax)
  if dopid:
    ii = selectparticles(iw=iw,kwdict=kw)
    if len(ii) > 0:
      pid = getattrwithsuffix(top,'pid',suffix)
      result = take(pid[:,id],ii)
    else:
      result = array([],'d')
  else:
    result = array([],'d')
  if lparallel and gather: return gatherarray(result,bcast=bcast)
  else: return result
#-------------------------------------------------------------------------
# Add the selectparticles documentation to each of the routines.
if sys.version[:5] != "1.5.1":
  if lparallel:
    _gatherdoc = (
"""  gather=1 When 1, all data is gathered to PE0
  bcast=0: when true, result is broadcast to all processors
""")
  else:
    _gatherdoc = ''
  _gatherdoc = _gatherdoc + ' lost=0: when true, gets lost particles'
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
  gettp.__doc__ = gettp.__doc__ + selectparticles.__doc__ + _gatherdoc
  getgaminv.__doc__ = getgaminv.__doc__ + selectparticles.__doc__ + _gatherdoc
  getpid.__doc__ = getpid.__doc__ + selectparticles.__doc__ + _gatherdoc
#-------------------------------------------------------------------------

##########################################################################
def getxxpslope(iw=0,iz=-1,kwdict={},checkargs=0):
  """
Calculates the x-x' slope from particle moments.
This returns a tuple containing (slope,xoffset,xpoffset,vz).
The product slope*vz gives the slope for x-vx.
  - iw=0: Window number
  - iz=-1: Z grid location
  - zl=None: When specified, lower range of selection region
  - zu=None: When specified, upper range of selection region
  - zc=None: When specified, center of range of selection region
  - slopejs=-1: Species whose moments are used to calculate the slope
                -1 means use data combined from all species.
  """
  if checkargs:
    # --- This is only needed since no other routines take the slopejs
    # --- argument.
    kwdictcopy = kwdict.copy()
    if 'slopejs' in kwdictcopy.keys(): del kwdictcopy['slopejs']
    return kwdictcopy
  zl = kwdict.get('zl')
  zu = kwdict.get('zu')
  zc = kwdict.get('zc')
  js = kwdict.get('js',-1)
  slopejs = kwdict.get('slopejs',js)
  if zl is not None and zu is not None: zc = 0.5*(zl + zu)
  if iw > 0: zc = 0.5*(top.zwindows[0,iw]+top.zwindows[1,iw])
  if zc is not None:
    iz = nint(floor((zc - w3d.zmmin)/w3d.dz))
    wz = (zc - w3d.zmmin)/w3d.dz - iz
  else:
    wz = 0.
  if 0 <= iz <= w3d.nz:
    izp1 = iz + 1
    if iz == w3d.nz: izp1 = iz
    xxpbar = top.xxpbarz[iz,slopejs]*(1.-wz) + top.xxpbarz[izp1,slopejs]*wz
    xbar   = top.xbarz[iz,slopejs]*(1.-wz)   + top.xbarz[izp1,slopejs]*wz
    xpbar  = top.xpbarz[iz,slopejs]*(1.-wz)  + top.xpbarz[izp1,slopejs]*wz
    xrms   = top.xrmsz[iz,slopejs]*(1.-wz)   + top.xrmsz[izp1,slopejs]*wz
    vzbar  = top.vzbarz[iz,slopejs]*(1.-wz)  + top.vzbarz[izp1,slopejs]*wz
    slope = (xxpbar-xbar*xpbar)/xrms**2
    xoffset = xbar
    xpoffset = xpbar
    vz = vzbar
  elif iw <= 0:
    slope = ((top.xxpbar[0,slopejs]-top.xbar[0,slopejs]*top.xpbar[0,slopejs])/
             top.xrms[0,slopejs]**2)
    xoffset = top.xbar[0,slopejs]
    xpoffset = top.xpbar[0,slopejs]
    vz = top.vzbar[0,slopejs]
  else:
    slope = 0.
    xoffset = 0.
    xpoffset = 0.
    vz = 0.
  return (slope,xoffset,xpoffset,vz)
#-------------------------------------------------------------------------
def getyypslope(iw=0,iz=-1,kwdict={},checkargs=0):
  """
Calculates the y-y' slope from particle moments.
This returns a tuple containing (slope,yoffset,ypoffset,vz).
The product slope*vz gives the slope for y-vy.
  - iw=0: Window number
  - iz=-1: Z grid location
  - zl=None: When specified, lower range of selection region
  - zu=None: When specified, upper range of selection region
  - zc=None: When specified, center of range of selection region
  - slopejs=-1: Species whose moments are used to calculate the slope
                -1 means use data combined from all species.
  """
  if checkargs:
    # --- This is only needed since no other routines take the slopejs
    # --- argument.
    kwdictcopy = kwdict.copy()
    if 'slopejs' in kwdictcopy.keys(): del kwdictcopy['slopejs']
    return kwdictcopy
  zl = kwdict.get('zl')
  zu = kwdict.get('zu')
  zc = kwdict.get('zc')
  js = kwdict.get('js',-1)
  slopejs = kwdict.get('slopejs',js)
  if zl is not None and zu is not None: zc = 0.5*(zl + zu)
  if iw > 0: zc = 0.5*(top.zwindows[0,iw]+top.zwindows[1,iw])
  if zc is not None:
    iz = nint(floor((zc - w3d.zmmin)/w3d.dz))
    wz = (zc - w3d.zmmin)/w3d.dz - iz
  else:
    wz = 0.
  if 0 <= iz <= w3d.nz:
    izp1 = iz + 1
    if iz == w3d.nz: izp1 = iz
    yypbar = top.yypbarz[iz,slopejs]*(1.-wz) + top.yypbarz[izp1,slopejs]*wz
    ybar   = top.ybarz[iz,slopejs]*(1.-wz)   + top.ybarz[izp1,slopejs]*wz
    ypbar  = top.ypbarz[iz,slopejs]*(1.-wz)  + top.ypbarz[izp1,slopejs]*wz
    yrms   = top.yrmsz[iz,slopejs]*(1.-wz)   + top.yrmsz[izp1,slopejs]*wz
    vzbar  = top.vzbarz[iz,slopejs]*(1.-wz)  + top.vzbarz[izp1,slopejs]*wz
    slope = (yypbar-ybar*ypbar)/yrms**2
    yoffset = ybar
    ypoffset = ypbar
    vz = vzbar
  elif iw <= 0:
    slope = ((top.yypbar[0,slopejs]-top.ybar[0,slopejs]*top.ypbar[0,slopejs])/
             top.yrms[0,slopejs]**2)
    yoffset = top.ybar[0,slopejs]
    ypoffset = top.ypbar[0,slopejs]
    vz = top.vzbar[0,slopejs]
  else:
    slope = 0.
    yoffset = 0.
    ypoffset = 0.
    vz = 0.
  return (slope,yoffset,ypoffset,vz)
#-------------------------------------------------------------------------
def getvzrange(kwdict={}):
  "Returns a tuple containg the Vz range for plots"
  if (top.vzrng != 0.):
    vzmax = (1. + top.vtilt)*top.vbeam*(1.+top.vzrng) - top.vzshift
    vzmin = (1. - top.vtilt)*top.vbeam*(1.-top.vzrng) - top.vzshift
  else:
    js = kwdict.get('js',-1)
    slopejs = kwdict.get('slopejs',js)
    vzmax = top.vzmaxp[slopejs] + 0.1*(top.vzmaxp[slopejs]-top.vzminp[slopejs])
    vzmin = top.vzminp[slopejs] - 0.1*(top.vzmaxp[slopejs]-top.vzminp[slopejs])
  return (vzmin,vzmax)


#-------------------------------------------------------------------------
#-------------------------------------------------------------------------
def addparticles(x=0.,y=0.,z=0.,vx=0.,vy=0.,vz=0.,gi=1.,pid=1.,js=0,
                 lallindomain=None,zmmin=None,zmmax=None,lmomentum=false,
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
                      This is automatically set to true when the code is in
                      slice mode, i.e. after package('wxy'). Except if the
                      option is explicitly set.
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
  try:              lenpid = len(pid)
  except TypeError: lenpid = 1

  # --- Max length of input arrays
  maxlen = max(lenx,leny,lenz,lenvx,lenvy,lenvz,lengi,lenpid)
  assert lenx==maxlen or lenx==1,"Length of x doesn't match len of others"
  assert leny==maxlen or leny==1,"Length of y doesn't match len of others"
  assert lenz==maxlen or lenz==1,"Length of z doesn't match len of others"
  assert lenvx==maxlen or lenvx==1,"Length of vx doesn't match len of others"
  assert lenvy==maxlen or lenvy==1,"Length of vy doesn't match len of others"
  assert lenvz==maxlen or lenvz==1,"Length of vz doesn't match len of others"
  assert lengi==maxlen or lengi==1,"Length of gi doesn't match len of others"
  assert lenpid==maxlen or lenpid==1,"Length of pid doesn't match len of others"

  # --- Convert all to arrays of length maxlen, broadcasting scalars
  x = array(x)*ones(maxlen,'d')
  y = array(y)*ones(maxlen,'d')
  z = array(z)*ones(maxlen,'d')
  vx = array(vx)*ones(maxlen,'d')
  vy = array(vy)*ones(maxlen,'d')
  vz = array(vz)*ones(maxlen,'d')
  gi = array(gi)*ones(maxlen,'d')
  pid = array(pid)*ones([maxlen,top.npid],'d')

  # --- Set extent of domain
  if not lparallel:
    if zmmin is None: zmmin = w3d.zmmin + top.zbeam
    if zmmax is None: zmmax = w3d.zmmax + top.zbeam
  else:
    if zmmin is None: zmmin = top.zpslmin[me] + top.zbeam
    if zmmax is None: zmmax = top.zpslmax[me] + top.zbeam

  # --- When running in slice mode, automatically set lallindomain to true.
  # --- It is assumed that all particles will be within the specified domain,
  # --- since in the slice mode, the z of the particles is ignored.
  # --- The user can still set lallindomain to false to override this.
  if lallindomain is None:
    if getcurrpkg() == 'wxy': lallindomain = true
    else:                     lallindomain = false

  # --- Now data can be passed into the fortran addparticles routine.
  addpart(maxlen,top.npid,x,y,z,vx,vy,vz,gi,pid,js+1,lallindomain,zmmin,zmmax,lmomentum)
 
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

































































