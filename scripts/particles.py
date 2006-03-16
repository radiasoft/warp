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
particles_version = "$Id: particles.py,v 1.36 2006/03/16 19:02:00 dave Exp $"

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
    #rr = top.nps[0]*RandomArray.random(top.nps[0])
    rr = top.nps[0]*ranf(zeros(top.nps[0]))
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

##########################################################################
#-------------------------------------------------------------------------
def getattrwithsuffix(object=top,name=None,suffix='',checkmain=1,pkg='top',
                      default=None):
  """
Helper function modeled on getattr. Like getattr, given an object and a name,
it will return object.name. Optionally, a suffix can be given, in which case,
it will return the equivalent of getattr(object,name+suffix).  The object can
be None, in which case the main dictionary will be searched, returning
__main__.__dict__[name+suffix]. The object can also be a dictionary or an
opened PDB file. Finally, if the object is specified but name+suffix is not
found and if checkmain is true, then the main dictionary will again be
searched. If nothing is found and default is specified, it is returned,
otherwise an exception is raised.
  """
  assert name is not None,"A name must be supplied"
  if object is not None:
    # --- This will work if object is a Warp package or a generic instance,
    # --- including a pdb file.
    try:
      result = getattr(object,name+suffix)
      return result
    except (AttributeError,KeyError):
      pass
  if type(object) is DictionaryType:
    try:
      result = object[name+suffix]
      return result
    except KeyError:
      pass
  if type(object) is InstanceType and object.__class__ is PR.PR:
    # --- Try to get data from a pdb file that was written as part of the
    # --- specified package.
    try:
      result = object.read(name+suffix+'@'+pkg)
      return result
    except:
      pass
  if object is None or checkmain:
    import __main__
    try:
      result = __main__.__dict__[name+suffix]
      return result
    except KeyError:
      pass

  # --- If all else fails, try finding the name without the suffix.
  if suffix != '':
    return getattrwithsuffix(object=object,name=name,suffix='',
                             checkmain=checkmain,pkg=pkg,default=default)

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
Multiple selection criteria are now supported.
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
  - suffix=None: When specified, variables with the specified suffix will be
                 used rather than the arrays from top.
  - object=top: Object to get particle data from. Besides top, this can be an
                open PDB file, or a dictionary.
  """
  # --- Complete dictionary of possible keywords and their default values
  kwdefaults = {"js":0,"jslist":None,"win":None,"z":None,
                "ix":None,"wx":1.,"iy":None,"wy":1.,"iz":None,"wz":1.,
                "zl":None,"zu":None,"zc":None,"xc":None,"yc":None,"ii":None,
                "lost":false,"suffix":'',"object":top,"w3dobject":None,
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

  # --- If lost is true, append the suffix lost to the variable names
  if lost: suffixparticle = 'lost' + suffix
  else:    suffixparticle = suffix

  ins = getattrwithsuffix(object,'ins',suffixparticle)
  nps = getattrwithsuffix(object,'nps',suffixparticle)
  ns = len(ins)

  # --- If jslist defined, call selectparticles repeatedly for each species
  # --- on the list
  del kwvalues['jslist']    # Remove list so selectparticles is subsequently
                            # called with one species at a time
  if jslist is not None:
    if jslist == -1:    jslist = range(0,ns)
    partlist = array([])
    for js in jslist:
        kwvalues['js'] = js
        newparts = selectparticles(iw, kwvalues)
        partlist = array(list(partlist)+list(newparts))
    return partlist

  # --- If the w3dobject was not passed in, use w3d, or if the object was
  # --- set to a PDB file, use it.
  if w3dobject is None:
    if object is not top: w3dobject = object
    else:                 w3dobject = w3d


  ir1 = ins[js] - 1
  ir2 = ins[js] + nps[js] - 1

  if ir2 <= ir1: return array([])

  indices = arrayrange(ir1,ir2)
  islice = slice(ir1,ir2)

  if zl is not None or zu is not None:
    if z is None: z = getattrwithsuffix(object,'zp',suffixparticle)
    if zl is None: zl = -largepos
    if zu is None: zu = +largepos
    if zl > zu: print "Warning: zl > zu"
    if ii is not None:
      z = take(z,ii)
      islice = slice(len(ii))
      indices = ii
    ii=compress(logical_and(less(zl,z[islice]),less(z[islice],zu)),indices)
  if ix is not None:
    xmmin = getattrwithsuffix(w3dobject,'xmmin',suffix,pkg='w3d')
    dx = getattrwithsuffix(w3dobject,'dx',suffix,pkg='w3d')
    xl = xmmin + ix*dx - wx*dx
    xu = xmmin + ix*dx + wx*dx
    x = getattrwithsuffix(object,'xp',suffixparticle)
    if ii is not None:
      x = take(x,ii)
      islice = slice(len(ii))
      indices = ii
    ii=compress(logical_and(less(xl,x[islice]),less(x[islice],xu)),indices)
  if iy is not None:
    ymmin = getattrwithsuffix(w3dobject,'ymmin',suffix,pkg='w3d')
    dy = getattrwithsuffix(w3dobject,'dy',suffix,pkg='w3d')
    yl = ymmin + iy*dy - wy*dy
    yu = ymmin + iy*dy + wy*dy
    y = getattrwithsuffix(object,'yp',suffixparticle)
    if ii is not None:
      y = take(y,ii)
      islice = slice(len(ii))
      indices = ii
    ii=compress(logical_and(less(yl,y[islice]),less(y[islice],yu)),indices)
  if iz is not None:
    zbeam = getattrwithsuffix(object,'zbeam',suffix)
    zmminglobal = getattrwithsuffix(w3dobject,'zmminglobal',suffix,pkg='w3d')
    dz = getattrwithsuffix(w3dobject,'dz',suffix,pkg='w3d')
    zl = zmminglobal + iz*dz - wz*dz + zbeam
    zu = zmminglobal + iz*dz + wz*dz + zbeam
    z = getattrwithsuffix(object,'zp',suffixparticle)
    if ii is not None:
      z = take(z,ii)
      islice = slice(len(ii))
      indices = ii
    ii=compress(logical_and(less(zl,z[islice]),less(z[islice],zu)),indices)
  if zc is not None:
    z = getattrwithsuffix(object,'zp',suffixparticle)
    dz = getattrwithsuffix(w3dobject,'dz',suffix,pkg='w3d')
    zl = zc - wz*dz
    zu = zc + wz*dz
    if ii is not None:
      z = take(z,ii)
      islice = slice(len(ii))
      indices = ii
    ii=compress(logical_and(less(zl,z[islice]),less(z[islice],zu)),indices)
  if xc is not None:
    x = getattrwithsuffix(object,'xp',suffixparticle)
    dx = getattrwithsuffix(w3dobject,'dx',suffix,pkg='w3d')
    xl = xc - wx*dx
    xu = xc + wx*dx
    if ii is not None:
      x = take(x,ii)
      islice = slice(len(ii))
      indices = ii
    ii=compress(logical_and(less(xl,x[islice]),less(x[islice],xu)),indices)
  if yc is not None:
    y = getattrwithsuffix(object,'yp',suffixparticle)
    dy = getattrwithsuffix(w3dobject,'dy',suffix,pkg='w3d')
    yl = yc - wy*dy
    yu = yc + wy*dy
    if ii is not None:
      y = take(y,ii)
      islice = slice(len(ii))
      indices = ii
    ii=compress(logical_and(less(yl,y[islice]),less(y[islice],yu)),indices)
  if iw < 0:
    # --- Add some smarts about choosing which method to use.
    # --- If the number of particles is increasing, then use the
    # --- population sampling. Otherwise use the preset subsets.
    # --- Also, if the number of species has changed, then also
    # --- use the population sampling.
    try: selectparticles.npsprev
    except: selectparticles.npsprev = nps + 0
    if len(selectparticles.npsprev) != ns:
      selectparticles.npsprev = zeros(ns)
    if selectparticles.npsprev[js] >= nps[js]:
      if psubset==[]: setup_subsets()
      if -iw > len(psubset): raise "Bad window number"
      if ii is None: nn = nps[js]
      else:          nn = len(ii)
      subset = compress(less(psubset[-iw-1],nn),psubset[-iw-1])
      if ii is None: ii = ir1 + subset
      else:          ii = take(ii,subset)
    else:
      # --- Once this method is used, always use it.
      selectparticles.npsprev = zeros(ns)
      npplot = getattrwithsuffix(object,'npplot',suffix)
      if nps[js] <= npplot[-iw-1]:
        ii = indices
      else:
        ii = populationsample(indices,npplot[-iw-1])
  elif iw == 0:
    if ii is None: ii = indices
  else:
    zbeam = getattrwithsuffix(object,'zbeam',suffix)
    zwindows = getattrwithsuffix(object,'zwindows',suffix)
    if win is None: win = zwindows[:,iw] + zbeam
    if len(shape(win)) == 2: win = win[:,iw]
    if z is None: z = getattrwithsuffix(object,'zp',suffixparticle)
    if ii is not None:
      z = take(z,ii)
      islice = slice(len(ii))
      indices = ii
    ii=compress(logical_and(less(win[0],z[islice]),less(z[islice],win[1])),
                indices)
  uz = getattrwithsuffix(object,'uzp',suffixparticle)
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
  suffix = ((kw.get('lost',0) and 'lost') or '') + kw.get('suffix','')
  object = kw.get('object',top)
  if len(ii) > 0:
    x = getattrwithsuffix(object,'xp',suffix)
    result = take(x,ii)
  else:
    result = array([],'d')
  if lparallel and gather: return gatherarray(result,bcast=bcast)
  else: return result
#-------------------------------------------------------------------------
def gety(iw=0,gather=1,bcast=0,**kw):
  "Returns the Y positions."
  ii = selectparticles(iw=iw,kwdict=kw)
  suffix = ((kw.get('lost',0) and 'lost') or '') + kw.get('suffix','')
  object = kw.get('object',top)
  if len(ii) > 0:
    y = getattrwithsuffix(object,'yp',suffix)
    result = take(y,ii)
  else:
    result = array([],'d')
  if lparallel and gather: return gatherarray(result,bcast=bcast)
  else: return result
#-------------------------------------------------------------------------
def getz(iw=0,gather=1,bcast=0,**kw):
  "Returns the Z positions."
  ii = selectparticles(iw=iw,kwdict=kw)
  suffix = ((kw.get('lost',0) and 'lost') or '') + kw.get('suffix','')
  object = kw.get('object',top)
  if len(ii) > 0:
    z = getattrwithsuffix(object,'zp',suffix)
    result = take(z,ii)
  else:
    result = array([],'d')
  if lparallel and gather: return gatherarray(result,bcast=bcast)
  else: return result
#-------------------------------------------------------------------------
def getr(iw=0,gather=1,bcast=0,**kw):
  "Returns the R postions."
  ii = selectparticles(iw=iw,kwdict=kw)
  suffix = ((kw.get('lost',0) and 'lost') or '') + kw.get('suffix','')
  object = kw.get('object',top)
  if len(ii) > 0:
    x = getattrwithsuffix(object,'xp',suffix)
    y = getattrwithsuffix(object,'yp',suffix)
    result = sqrt(take(x,ii)**2 + take(y,ii)**2)
  else:
    result = array([],'d')
  if lparallel and gather: return gatherarray(result,bcast=bcast)
  else: return result
#-------------------------------------------------------------------------
def gettheta(iw=0,gather=1,bcast=0,**kw):
  "Returns the theta postions."
  ii = selectparticles(iw=iw,kwdict=kw)
  suffix = ((kw.get('lost',0) and 'lost') or '') + kw.get('suffix','')
  object = kw.get('object',top)
  if len(ii) > 0:
    x = getattrwithsuffix(object,'xp',suffix)
    y = getattrwithsuffix(object,'yp',suffix)
    result = arctan2(take(y,ii),take(x,ii))
  else:
    result = array([],'d')
  if lparallel and gather: return gatherarray(result,bcast=bcast)
  else: return result
#-------------------------------------------------------------------------
def getvx(iw=0,gather=1,bcast=0,**kw):
  "Returns the X velocity."
  ii = selectparticles(iw=iw,kwdict=kw)
  suffix = ((kw.get('lost',0) and 'lost') or '') + kw.get('suffix','')
  object = kw.get('object',top)
  if len(ii) > 0:
    ux = getattrwithsuffix(object,'uxp',suffix)
    gaminv = getattrwithsuffix(object,'gaminv',suffix)
    result = take(ux,ii)*take(gaminv,ii)
  else:
    result = array([],'d')
  if lparallel and gather: return gatherarray(result,bcast=bcast)
  else: return result
#-------------------------------------------------------------------------
def getvy(iw=0,gather=1,bcast=0,**kw):
  "Returns the Y velocity."
  ii = selectparticles(iw=iw,kwdict=kw)
  suffix = ((kw.get('lost',0) and 'lost') or '') + kw.get('suffix','')
  object = kw.get('object',top)
  if len(ii) > 0:
    uy = getattrwithsuffix(object,'uyp',suffix)
    gaminv = getattrwithsuffix(object,'gaminv',suffix)
    result = take(uy,ii)*take(gaminv,ii)
  else:
    result = array([],'d')
  if lparallel and gather: return gatherarray(result,bcast=bcast)
  else: return result
#-------------------------------------------------------------------------
def getvz(iw=0,gather=1,bcast=0,**kw):
  "Returns the Z velocity."
  ii = selectparticles(iw=iw,kwdict=kw)
  suffix = ((kw.get('lost',0) and 'lost') or '') + kw.get('suffix','')
  object = kw.get('object',top)
  if len(ii) > 0:
    uz = getattrwithsuffix(object,'uzp',suffix)
    gaminv = getattrwithsuffix(object,'gaminv',suffix)
    result = take(uz,ii)*take(gaminv,ii)
  else:
    result = array([],'d')
  if lparallel and gather: return gatherarray(result,bcast=bcast)
  else: return result
#-------------------------------------------------------------------------
def getvr(iw=0,gather=1,bcast=0,**kw):
  "Returns the radial velocity."
  ii = selectparticles(iw=iw,kwdict=kw)
  suffix = ((kw.get('lost',0) and 'lost') or '') + kw.get('suffix','')
  object = kw.get('object',top)
  if len(ii) > 0:
    x = getattrwithsuffix(object,'xp',suffix)
    y = getattrwithsuffix(object,'yp',suffix)
    ux = getattrwithsuffix(object,'uxp',suffix)
    uy = getattrwithsuffix(object,'uyp',suffix)
    gaminv = getattrwithsuffix(object,'gaminv',suffix)
    tt = arctan2(take(y,ii),take(x,ii))
    result = (take(ux,ii)*cos(tt) + take(uy,ii)*sin(tt))*take(gaminv,ii)
  else:
    result = array([],'d')
  if lparallel and gather: return gatherarray(result,bcast=bcast)
  else: return result
#-------------------------------------------------------------------------
def getvtheta(iw=0,gather=1,bcast=0,**kw):
  "Returns the azimuthal velocity."
  ii = selectparticles(iw=iw,kwdict=kw)
  suffix = ((kw.get('lost',0) and 'lost') or '') + kw.get('suffix','')
  object = kw.get('object',top)
  if len(ii) > 0:
    x = getattrwithsuffix(object,'xp',suffix)
    y = getattrwithsuffix(object,'yp',suffix)
    ux = getattrwithsuffix(object,'uxp',suffix)
    uy = getattrwithsuffix(object,'uyp',suffix)
    gaminv = getattrwithsuffix(object,'gaminv',suffix)
    tt = arctan2(take(y,ii),take(x,ii))
    result = (-take(ux,ii)*sin(tt) + take(uy,ii)*cos(tt))*take(gaminv,ii)
  else:
    result = array([],'d')
  if lparallel and gather: return gatherarray(result,bcast=bcast)
  else: return result
#-------------------------------------------------------------------------
def getux(iw=0,gather=1,bcast=0,**kw):
  "Returns the X momentum over mass."
  ii = selectparticles(iw=iw,kwdict=kw)
  suffix = ((kw.get('lost',0) and 'lost') or '') + kw.get('suffix','')
  object = kw.get('object',top)
  if len(ii) > 0:
    ux = getattrwithsuffix(object,'uxp',suffix)
    result = take(ux,ii)
  else:
    result = array([],'d')
  if lparallel and gather: return gatherarray(result,bcast=bcast)
  else: return result
#-------------------------------------------------------------------------
def getuy(iw=0,gather=1,bcast=0,**kw):
  "Returns the Y momentum over mass."
  ii = selectparticles(iw=iw,kwdict=kw)
  suffix = ((kw.get('lost',0) and 'lost') or '') + kw.get('suffix','')
  object = kw.get('object',top)
  if len(ii) > 0:
    uy = getattrwithsuffix(object,'uyp',suffix)
    result = take(uy,ii)
  else:
    result = array([],'d')
  if lparallel and gather: return gatherarray(result,bcast=bcast)
  else: return result
#-------------------------------------------------------------------------
def getuz(iw=0,gather=1,bcast=0,**kw):
  "Returns the Z momentum over mass."
  ii = selectparticles(iw=iw,kwdict=kw)
  suffix = ((kw.get('lost',0) and 'lost') or '') + kw.get('suffix','')
  object = kw.get('object',top)
  if len(ii) > 0:
    uz = getattrwithsuffix(object,'uzp',suffix)
    result = take(uz,ii)
  else:
    result = array([],'d')
  if lparallel and gather: return gatherarray(result,bcast=bcast)
  else: return result
#-------------------------------------------------------------------------
def getxp(iw=0,gather=1,bcast=0,**kw):
  "Returns the X velocity over the Z velocity (X')."
  ii = selectparticles(iw=iw,kwdict=kw)
  suffix = ((kw.get('lost',0) and 'lost') or '') + kw.get('suffix','')
  object = kw.get('object',top)
  if len(ii) > 0:
    ux = getattrwithsuffix(object,'uxp',suffix)
    uz = getattrwithsuffix(object,'uzp',suffix)
    result = take(ux,ii)/take(uz,ii)
  else:
    result = array([],'d')
  if lparallel and gather: return gatherarray(result,bcast=bcast)
  else: return result
#-------------------------------------------------------------------------
def getyp(iw=0,gather=1,bcast=0,**kw):
  "Returns the Y velocity over the Z velocity (Y')."
  ii = selectparticles(iw=iw,kwdict=kw)
  suffix = ((kw.get('lost',0) and 'lost') or '') + kw.get('suffix','')
  object = kw.get('object',top)
  if len(ii) > 0:
    uy = getattrwithsuffix(object,'uyp',suffix)
    uz = getattrwithsuffix(object,'uzp',suffix)
    result = take(uy,ii)/take(uz,ii)
  else:
    result = array([],'d')
  if lparallel and gather: return gatherarray(result,bcast=bcast)
  else: return result
#-------------------------------------------------------------------------
def getrp(iw=0,gather=1,bcast=0,**kw):
  "Returns the radial velocity over the Z velocity (R')."
  ii = selectparticles(iw=iw,kwdict=kw)
  suffix = ((kw.get('lost',0) and 'lost') or '') + kw.get('suffix','')
  object = kw.get('object',top)
  if len(ii) > 0:
    x = getattrwithsuffix(object,'xp',suffix)
    y = getattrwithsuffix(object,'yp',suffix)
    ux = getattrwithsuffix(object,'uxp',suffix)
    uy = getattrwithsuffix(object,'uyp',suffix)
    uz = getattrwithsuffix(object,'uzp',suffix)
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
  suffix = ((kw.get('lost',0) and 'lost') or '') + kw.get('suffix','')
  object = kw.get('object',top)
  if len(ii) > 0:
    x = getattrwithsuffix(object,'xp',suffix)
    y = getattrwithsuffix(object,'yp',suffix)
    ux = getattrwithsuffix(object,'uxp',suffix)
    uy = getattrwithsuffix(object,'uyp',suffix)
    uz = getattrwithsuffix(object,'uzp',suffix)
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
  suffix = ((kw.get('lost',0) and 'lost') or '') + kw.get('suffix','')
  object = kw.get('object',top)
  if len(ii) > 0:
    gaminv = getattrwithsuffix(object,'gaminv',suffix)
    result = take(gaminv,ii)
  else:
    result = array([],'d')
  if lparallel and gather: return gatherarray(result,bcast=bcast)
  else: return result
#-------------------------------------------------------------------------
def getpid(iw=0,id=0,gather=1,bcast=0,**kw):
  """Returns particle id number.
  -id=0: which pid value to return
         if id=-1, returns all pids.
  """
  lost = kw.get('lost',0)
  suffix = ((kw.get('lost',0) and 'lost') or '') + kw.get('suffix','')
  object = kw.get('object',top)
  if lost:
    npidlostmax = getattrwithsuffix(object,'npidlostmax')
    dopid = (npidlostmax > 0)
  else:
    npid = getattrwithsuffix(object,'npid',suffix)
    dopid = (npid > 0)
  if dopid:
    ii = selectparticles(iw=iw,kwdict=kw)
    if len(ii) > 0:
      pid = getattrwithsuffix(object,'pid',suffix)
      if id >= 0: result = take(pid[:,id],ii)
      else:       result = take(pid[:,:],ii)
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
  suffix = kwdict.get('suffix','')
  object = kwdict.get('object',top)
  zwindows = getattrwithsuffix(object,'zwindows',suffix)
  if zl is not None and zu is not None: zc = 0.5*(zl + zu)
  if iw > 0: zc = 0.5*(zwindows[0,iw]+zwindows[1,iw])
  if zc is not None:
    zmmntmin = getattrwithsuffix(object,'zmmntmin',suffix)
    dzm = getattrwithsuffix(object,'dzm',suffix)
    iz = nint(floor((zc - zmmntmin)/dzm))
    wz =            (zc - zmmntmin)/dzm - iz
  else:
    wz = 0.
  nzmmnt = getattrwithsuffix(object,'nzmmnt',suffix)
  if 0 <= iz <= nzmmnt:
    izp1 = iz + 1
    if iz == nzmmnt: izp1 = iz
    xxpbarz = getattrwithsuffix(object,'xxpbarz',suffix)
    xbarz   = getattrwithsuffix(object,'xbarz',suffix)
    xpbarz  = getattrwithsuffix(object,'xpbarz',suffix)
    xrmsz   = getattrwithsuffix(object,'xrmsz',suffix)
    vzbarz  = getattrwithsuffix(object,'vzbarz',suffix)
    xxpbar = xxpbarz[iz,slopejs]*(1.-wz) + xxpbarz[izp1,slopejs]*wz
    xbar   = xbarz[iz,slopejs]*(1.-wz)   + xbarz[izp1,slopejs]*wz
    xpbar  = xpbarz[iz,slopejs]*(1.-wz)  + xpbarz[izp1,slopejs]*wz
    xrms   = xrmsz[iz,slopejs]*(1.-wz)   + xrmsz[izp1,slopejs]*wz
    vzbar  = vzbarz[iz,slopejs]*(1.-wz)  + vzbarz[izp1,slopejs]*wz
    slope = (xxpbar-xbar*xpbar)/dvnz(xrms)**2
    xoffset = xbar
    xpoffset = xpbar
    vz = vzbar
  elif iw <= 0:
    xxpbar = getattrwithsuffix(object,'xxpbar',suffix)
    xbar   = getattrwithsuffix(object,'xbar',suffix)
    xpbar  = getattrwithsuffix(object,'xpbar',suffix)
    xrms   = getattrwithsuffix(object,'xrms',suffix)
    vzbar  = getattrwithsuffix(object,'vzbar',suffix)
    slope = ((xxpbar[0,slopejs]-xbar[0,slopejs]*xpbar[0,slopejs])/
             dvnz(xrms[0,slopejs])**2)
    xoffset = xbar[0,slopejs]
    xpoffset = xpbar[0,slopejs]
    vz = vzbar[0,slopejs]
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
  suffix = kwdict.get('suffix','')
  object = kwdict.get('object',top)
  zwindows = getattrwithsuffix(object,'zwindows',suffix)
  if zl is not None and zu is not None: zc = 0.5*(zl + zu)
  if iw > 0: zc = 0.5*(zwindows[0,iw]+zwindows[1,iw])
  if zc is not None:
    zmmntmin = getattrwithsuffix(object,'zmmntmin',suffix)
    dzm = getattrwithsuffix(object,'dzm',suffix)
    iz = nint(floor((zc - zmmntmin)/dzm))
    wz =            (zc - zmmntmin)/dzm - iz
  else:
    wz = 0.
  nzmmnt = getattrwithsuffix(object,'nzmmnt',suffix)
  if 0 <= iz <= nzmmnt:
    izp1 = iz + 1
    if iz == nzmmnt: izp1 = iz
    yypbarz = getattrwithsuffix(object,'yypbarz',suffix)
    ybarz   = getattrwithsuffix(object,'ybarz',suffix)
    ypbarz  = getattrwithsuffix(object,'ypbarz',suffix)
    yrmsz   = getattrwithsuffix(object,'yrmsz',suffix)
    vzbarz  = getattrwithsuffix(object,'vzbarz',suffix)
    yypbar = yypbarz[iz,slopejs]*(1.-wz) + yypbarz[izp1,slopejs]*wz
    ybar   = ybarz[iz,slopejs]*(1.-wz)   + ybarz[izp1,slopejs]*wz
    ypbar  = ypbarz[iz,slopejs]*(1.-wz)  + ypbarz[izp1,slopejs]*wz
    yrms   = yrmsz[iz,slopejs]*(1.-wz)   + yrmsz[izp1,slopejs]*wz
    vzbar  = vzbarz[iz,slopejs]*(1.-wz)  + vzbarz[izp1,slopejs]*wz
    slope = (yypbar-ybar*ypbar)/dvnz(yrms)**2
    yoffset = ybar
    ypoffset = ypbar
    vz = vzbar
  elif iw <= 0:
    yypbar = getattrwithsuffix(object,'yypbar',suffix)
    ybar   = getattrwithsuffix(object,'ybar',suffix)
    ypbar  = getattrwithsuffix(object,'ypbar',suffix)
    yrms   = getattrwithsuffix(object,'yrms',suffix)
    vzbar  = getattrwithsuffix(object,'vzbar',suffix)
    slope = ((yypbar[0,slopejs]-ybar[0,slopejs]*ypbar[0,slopejs])/
             dvnz(yrms[0,slopejs])**2)
    yoffset = ybar[0,slopejs]
    ypoffset = ypbar[0,slopejs]
    vz = vzbar[0,slopejs]
  else:
    slope = 0.
    yoffset = 0.
    ypoffset = 0.
    vz = 0.
  return (slope,yoffset,ypoffset,vz)
#-------------------------------------------------------------------------
def getvzrange(kwdict={}):
  "Returns a tuple containg the Vz range for plots"
  suffix = kwdict.get('suffix','')
  object = kwdict.get('object',top)
  vzrng  = getattrwithsuffix(object,'vzrng',suffix)
  if (vzrng != 0.):
    vtilt  = getattrwithsuffix(object,'vtilt',suffix)
    vbeam  = getattrwithsuffix(object,'vbeam',suffix)
    vzshift  = getattrwithsuffix(object,'vzshift',suffix)
    vzmax = (1. + vtilt)*vbeam*(1.+vzrng) - vzshift
    vzmin = (1. - vtilt)*vbeam*(1.-vzrng) - vzshift
  else:
    js = kwdict.get('js',-1)
    slopejs = kwdict.get('slopejs',js)
    vzminp  = getattrwithsuffix(object,'vzminp',suffix)
    vzmaxp  = getattrwithsuffix(object,'vzmaxp',suffix)
    vzmax = vzmaxp[slopejs] + 0.1*(vzmaxp[slopejs]-vzminp[slopejs])
    vzmin = vzminp[slopejs] - 0.1*(vzmaxp[slopejs]-vzminp[slopejs])
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
                     which defaults to 1. (gi is 1/gamma, the relatistic
                     paramter)
  pid: additional particle information, such as an ID number or weight.
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

  # --- Set time of creation and ssn
  if top.tpid>0: pid[:,top.tpid-1]=top.time
  if top.spid>0: 
    pid[:,top.spid-1]=top.ssn+arange(maxlen)
    top.ssn += maxlen

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

    
