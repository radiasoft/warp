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
particles_version = "$Id: particles.py,v 1.62 2007/10/29 23:33:34 dave Exp $"

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
def setup_subsets(js=0,pgroup=None):
  """
Adds plotting subset to the list
  - js=0 is the species to create a subset for
  - pgroup=top.pgroup: particle group to base subsets on
  """
  global psubset
  if pgroup is None: pgroup = top.pgroup
  if lparallel:
    totalnp = parallelsum(pgroup.nps[js])
    if totalnp == 0: totalnp = 1
    fracnp = float(pgroup.nps[js])/float(totalnp)
  else:
    fracnp = 1.
  for i in xrange(0,len(top.npplot)):
    ntopick=min(pgroup.nps[js],int(top.npplot[i]*fracnp+0.5))
    ii = arange(pgroup.nps[0])
    #rr = pgroup.nps[0]*RandomArray.random(pgroup.nps[0])
    rr = pgroup.nps[0]*ranf(zeros(pgroup.nps[0]))
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
    typecode = gettypecode(array([population[0]]))
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
def _getobjectpgroup(kw):
  suffix = ((kw.get('lost',0) and 'lost') or '') + kw.get('suffix','')
  object = kw.get('object',top)
  if suffix == '':
    try:
      pgroup = kw.get('pgroup',object.pgroup)
    except AttributeError:
      pgroup = object
  else:
    pgroup = object
  return suffix,object,pgroup

#-------------------------------------------------------------------------
def _setindices(data,ii,islice,indices,ir1,ir2):
  if ii is not None:
    if isinstance(ii,slice): ii = arange(ii.start,ii.stop)
    data = take(data,ii)
    islice = slice(len(ii))
    indices = ii
  if indices is None: indices = arange(ir1,ir2)
  return data,islice,indices

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
  - x=top.xp: Coordinate for x range selection
  - y=top.yp: Coordinate for y range selection
  - z=top.zp: Coordinate for z range selection
  - ix=-1: When 0 <= ix <= nx, picks particles within xmesh[ix]+-wx*dx
  - wx=1.: Width of window around xmesh[ix]
  - iy=-1: When 0 <= iy <= ny, picks particles within ymesh[iy]+-wy*dy
  - wy=1.: Width of window around ymesh[iy]
  - iz=-1: When 0 <= iz <= nz, picks particles within zmesh[iz]+-wz*dz
  - wz=1.: Width of window around zmesh[iz]
  - xl=None: When specified, lower range in x of selection region
  - xu=None: When specified, upper range in x of selection region
  - yl=None: When specified, lower range in y of selection region
  - yu=None: When specified, upper range in y of selection region
  - zl=None: When specified, lower range in z of selection region
  - zu=None: When specified, upper range in z of selection region
  - zc=None: When specified, picks particles within zc+-wz*dz
  - xc=None: When specified, picks particles within xc+-wx*dx
  - yc=None: When specified, picks particles within yc+-wy*dy
  - ssn=None: When specified, returns the particle with the given ssn
              Raises and error if ssn's are not saved.
  - ii=None: If ii is supplied, it is just returned.
  - lost=false: When true, returns indices to the lost particles rather than
                the live particles
  - suffix=None: When specified, variables with the specified suffix will be
                 used rather than the arrays from top.
  - object=top: Object to get particle data from. Besides top, this can be an
                open PDB file, or a dictionary.
  - pgroup=top.pgroup: Particle group to get particles from 
  - usezerovzflag=0: When true, a check is made to skip particles with vz==0.
                     This is also done when top.clearlostpart==0.
  """
  # --- Complete dictionary of possible keywords and their default values
  kwdefaults = {"js":0,"jslist":None,"win":None,
                "x":None,"y":None,"z":None,
                "ix":None,"wx":1.,"iy":None,"wy":1.,"iz":None,"wz":1.,
                "xl":None,"xu":None,"yl":None,"yu":None,"zl":None,"zu":None,
                "zc":None,"xc":None,"yc":None,
                "ssn":None,"ii":None,
                "lost":false,"suffix":'',"object":top,"pgroup":top.pgroup,
                "w3dobject":None,"usezerovzflag":0,
                'checkargs':0,'allowbadargs':0}

  # --- Create dictionary of local values and copy it into local dictionary,
  # --- ignoring keywords not listed in kwdefaults.
  kwvalues = kwdefaults.copy()
  kwvalues.update(kw)
  kwvalues.update(kwdict)
  for arg in kwdefaults.keys(): exec arg+" = kwvalues['"+arg+"']"

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

  suffix,object,pgroup = _getobjectpgroup(kwvalues)

  # --- If lost is true, then strip off the 'lost' part of the suffix
  # --- for nonparticle quantities
  suffixparticle = suffix
  if lost: suffix = suffix[4:]

  ins = getattrwithsuffix(pgroup,'ins',suffixparticle)
  nps = getattrwithsuffix(pgroup,'nps',suffixparticle)
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
        if isinstance(newparts,slice):
          newparts = arange(newparts.start,newparts.stop)
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

  # --- Note that indices is only calculated if needed. If no down selection
  # --- is done, the islice will be returned.
  indices = None
  islice = slice(ir1,ir2)

  if ssn is not None:
    assert top.spid > 0,"ssn's are not used"
    id = getattrwithsuffix(pgroup,'pid',suffixparticle)[:,top.spid-1]
    id,islice,indices = _setindices(id,ii,islice,indices,ir1,ir2)
    ii = compress(nint(id[islice])==ssn,indices)

  if xl is not None or xu is not None:
    if x is None: x = getattrwithsuffix(pgroup,'xp',suffixparticle)
    if xl is None: xl = -largepos
    if xu is None: xu = +largepos
    if xl > xu: print "Warning: xl > xu"
    x,islice,indices = _setindices(x,ii,islice,indices,ir1,ir2)
    ii=compress(logical_and(less(xl,x[islice]),less(x[islice],xu)),indices)
  if yl is not None or yu is not None:
    if y is None: y = getattrwithsuffix(pgroup,'yp',suffixparticle)
    if yl is None: yl = -largepos
    if yu is None: yu = +largepos
    if yl > yu: print "Warning: yl > yu"
    y,islice,indices = _setindices(y,ii,islice,indices,ir1,ir2)
    ii=compress(logical_and(less(yl,y[islice]),less(y[islice],yu)),indices)
  if zl is not None or zu is not None:
    if z is None: z = getattrwithsuffix(pgroup,'zp',suffixparticle)
    if zl is None: zl = -largepos
    if zu is None: zu = +largepos
    if zl > zu: print "Warning: zl > zu"
    z,islice,indices = _setindices(z,ii,islice,indices,ir1,ir2)
    ii=compress(logical_and(less(zl,z[islice]),less(z[islice],zu)),indices)

  if ix is not None:
    xmmin = getattrwithsuffix(w3dobject,'xmmin',suffix,pkg='w3d')
    dx = getattrwithsuffix(w3dobject,'dx',suffix,pkg='w3d')
    xl = xmmin + ix*dx - wx*dx
    xu = xmmin + ix*dx + wx*dx
    if x is None: x = getattrwithsuffix(pgroup,'xp',suffixparticle)
    x,islice,indices = _setindices(x,ii,islice,indices,ir1,ir2)
    ii=compress(logical_and(less(xl,x[islice]),less(x[islice],xu)),indices)
  if iy is not None:
    ymmin = getattrwithsuffix(w3dobject,'ymmin',suffix,pkg='w3d')
    dy = getattrwithsuffix(w3dobject,'dy',suffix,pkg='w3d')
    yl = ymmin + iy*dy - wy*dy
    yu = ymmin + iy*dy + wy*dy
    if y is None: y = getattrwithsuffix(pgroup,'yp',suffixparticle)
    y,islice,indices = _setindices(y,ii,islice,indices,ir1,ir2)
    ii=compress(logical_and(less(yl,y[islice]),less(y[islice],yu)),indices)
  if iz is not None:
    zbeam = getattrwithsuffix(object,'zbeam',suffix)
    zmmin = getattrwithsuffix(w3dobject,'zmmin',suffix,pkg='w3d')
    dz = getattrwithsuffix(w3dobject,'dz',suffix,pkg='w3d')
    zl = zmmin + iz*dz - wz*dz + zbeam
    zu = zmmin + iz*dz + wz*dz + zbeam
    if z is None: z = getattrwithsuffix(pgroup,'zp',suffixparticle)
    z,islice,indices = _setindices(z,ii,islice,indices,ir1,ir2)
    ii=compress(logical_and(less(zl,z[islice]),less(z[islice],zu)),indices)

  if xc is not None:
    if x is None: x = getattrwithsuffix(pgroup,'xp',suffixparticle)
    dx = getattrwithsuffix(w3dobject,'dx',suffix,pkg='w3d')
    xl = xc - wx*dx
    xu = xc + wx*dx
    x,islice,indices = _setindices(x,ii,islice,indices,ir1,ir2)
    ii=compress(logical_and(less(xl,x[islice]),less(x[islice],xu)),indices)
  if yc is not None:
    if y is None: y = getattrwithsuffix(pgroup,'yp',suffixparticle)
    dy = getattrwithsuffix(w3dobject,'dy',suffix,pkg='w3d')
    yl = yc - wy*dy
    yu = yc + wy*dy
    y,islice,indices = _setindices(y,ii,islice,indices,ir1,ir2)
    ii=compress(logical_and(less(yl,y[islice]),less(y[islice],yu)),indices)
  if zc is not None:
    if z is None: z = getattrwithsuffix(pgroup,'zp',suffixparticle)
    dz = getattrwithsuffix(w3dobject,'dz',suffix,pkg='w3d')
    zl = zc - wz*dz
    zu = zc + wz*dz
    z,islice,indices = _setindices(z,ii,islice,indices,ir1,ir2)
    ii=compress(logical_and(less(zl,z[islice]),less(z[islice],zu)),indices)

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
      if indices is None: indices = arange(ir1,ir2)
      if nps[js] <= npplot[-iw-1]:
        ii = indices
      else:
        ii = populationsample(indices,npplot[-iw-1])
  elif iw == 0:
    if ii is None:
      # --- If all particles are being selected, then a slice object is returned
      # --- so that a direct refer to the particle array is returned, avoiding
      # --- memory copying.
      if indices is None: ii = islice
      else:               ii = indices
  else:
    zbeam = getattrwithsuffix(object,'zbeam',suffix)
    zwindows = getattrwithsuffix(object,'zwindows',suffix)
    if win is None: win = zwindows[:,iw] + zbeam
    if len(shape(win)) == 2: win = win[:,iw]
    if z is None: z = getattrwithsuffix(pgroup,'zp',suffixparticle)
    z,islice,indices = _setindices(z,ii,islice,indices,ir1,ir2)
    ii=compress(logical_and(less(win[0],z[islice]),less(z[islice],win[1])),
                indices)

  if usezerovzflag or top.clearlostpart == 0:
    if ii is None and indices is None: ii = arange(ir1,ir2)
    uz = getattrwithsuffix(pgroup,'uzp',suffixparticle)
    ii = compress(not_equal(take(uz,ii),0.),ii)

  return ii

#-------------------------------------------------------------------------
# The following return a specific coordinate of the selected particles
# More documetation added after they are declared.
# --- The checks of the len of ii are there in case the particle arrays
# --- are unallocated.
#-------------------------------------------------------------------------
_particlebcastdefault = [1]
def setgetparticlebcastdefault(defval=1):
  _particlebcastdefault[0] = defval
def getn(iw=0,gather=1,bcast=None,**kw):
  "Returns number of particles in selection."
  if bcast is None: bcast = _particlebcastdefault[0]
  ii = selectparticles(iw=iw,kwdict=kw)
  if isinstance(ii,slice): l = ii.stop - ii.start
  else:                    l = len(ii)
  if lparallel and gather: return globalsum(l)
  else:                    return l
#-------------------------------------------------------------------------
def getx(iw=0,gather=1,bcast=None,**kw):
  "Returns the X positions."
  if bcast is None: bcast = _particlebcastdefault[0]
  ii = selectparticles(iw=iw,kwdict=kw)
  suffix,object,pgroup = _getobjectpgroup(kw)
  if (isinstance(ii,slice) and ii.stop > ii.start) or len(ii) > 0:
    x = getattrwithsuffix(pgroup,'xp',suffix)
  if isinstance(ii,slice) and ii.stop > ii.start:
    result = x[ii]
  elif len(ii) > 0:
    result = take(x,ii)
  else:
    result = array([],'d')
  if lparallel and gather: return gatherarray(result,bcast=bcast)
  else: return result
#-------------------------------------------------------------------------
def gety(iw=0,gather=1,bcast=None,**kw):
  "Returns the Y positions."
  if bcast is None: bcast = _particlebcastdefault[0]
  ii = selectparticles(iw=iw,kwdict=kw)
  suffix,object,pgroup = _getobjectpgroup(kw)
  if (isinstance(ii,slice) and ii.stop > ii.start) or len(ii) > 0:
    y = getattrwithsuffix(pgroup,'yp',suffix)
  if isinstance(ii,slice) and ii.stop > ii.start:
    result = y[ii]
  elif len(ii) > 0:
    result = take(y,ii)
  else:
    result = array([],'d')
  if lparallel and gather: return gatherarray(result,bcast=bcast)
  else: return result
#-------------------------------------------------------------------------
def getz(iw=0,gather=1,bcast=None,**kw):
  "Returns the Z positions."
  if bcast is None: bcast = _particlebcastdefault[0]
  ii = selectparticles(iw=iw,kwdict=kw)
  suffix,object,pgroup = _getobjectpgroup(kw)
  if (isinstance(ii,slice) and ii.stop > ii.start) or len(ii) > 0:
    z = getattrwithsuffix(pgroup,'zp',suffix)
  if isinstance(ii,slice) and ii.stop > ii.start:
    result = z[ii]
  elif len(ii) > 0:
    result = take(z,ii)
  else:
    result = array([],'d')
  if lparallel and gather: return gatherarray(result,bcast=bcast)
  else: return result
#-------------------------------------------------------------------------
def getr(iw=0,gather=1,bcast=None,**kw):
  "Returns the R postions."
  if bcast is None: bcast = _particlebcastdefault[0]
  ii = selectparticles(iw=iw,kwdict=kw)
  suffix,object,pgroup = _getobjectpgroup(kw)
  if (isinstance(ii,slice) and ii.stop > ii.start) or len(ii) > 0:
    x = getattrwithsuffix(pgroup,'xp',suffix)
    y = getattrwithsuffix(pgroup,'yp',suffix)
  if isinstance(ii,slice) and ii.stop > ii.start:
    result = sqrt(x[ii]**2 + y[ii]**2)
  elif len(ii) > 0:
    result = sqrt(take(x,ii)**2 + take(y,ii)**2)
  else:
    result = array([],'d')
  if lparallel and gather: return gatherarray(result,bcast=bcast)
  else: return result
#-------------------------------------------------------------------------
def gettheta(iw=0,gather=1,bcast=None,**kw):
  "Returns the theta postions."
  if bcast is None: bcast = _particlebcastdefault[0]
  ii = selectparticles(iw=iw,kwdict=kw)
  suffix,object,pgroup = _getobjectpgroup(kw)
  if (isinstance(ii,slice) and ii.stop > ii.start) or len(ii) > 0:
    x = getattrwithsuffix(pgroup,'xp',suffix)
    y = getattrwithsuffix(pgroup,'yp',suffix)
  if isinstance(ii,slice) and ii.stop > ii.start:
    result = arctan2(y[ii],x[ii])
  elif len(ii) > 0:
    result = arctan2(take(y,ii),take(x,ii))
  else:
    result = array([],'d')
  if lparallel and gather: return gatherarray(result,bcast=bcast)
  else: return result
#-------------------------------------------------------------------------
def getvx(iw=0,gather=1,bcast=None,**kw):
  "Returns the X velocity."
  if bcast is None: bcast = _particlebcastdefault[0]
  ii = selectparticles(iw=iw,kwdict=kw)
  suffix,object,pgroup = _getobjectpgroup(kw)
  if (isinstance(ii,slice) and ii.stop > ii.start) or len(ii) > 0:
    ux = getattrwithsuffix(pgroup,'uxp',suffix)
    gaminv = getattrwithsuffix(pgroup,'gaminv',suffix)
  if isinstance(ii,slice) and ii.stop > ii.start:
    if top.lrelativ: result = ux[ii]*gaminv[ii]
    else:            result = ux[ii]
  elif len(ii) > 0:
    result = take(ux,ii)*take(gaminv,ii)
  else:
    result = array([],'d')
  if lparallel and gather: return gatherarray(result,bcast=bcast)
  else: return result
#-------------------------------------------------------------------------
def getvy(iw=0,gather=1,bcast=None,**kw):
  "Returns the Y velocity."
  if bcast is None: bcast = _particlebcastdefault[0]
  ii = selectparticles(iw=iw,kwdict=kw)
  suffix,object,pgroup = _getobjectpgroup(kw)
  if (isinstance(ii,slice) and ii.stop > ii.start) or len(ii) > 0:
    uy = getattrwithsuffix(pgroup,'uyp',suffix)
    gaminv = getattrwithsuffix(pgroup,'gaminv',suffix)
  if isinstance(ii,slice) and ii.stop > ii.start:
    if top.lrelativ: result = uy[ii]*gaminv[ii]
    else:            result = uy[ii]
  elif len(ii) > 0:
    result = take(uy,ii)*take(gaminv,ii)
  else:
    result = array([],'d')
  if lparallel and gather: return gatherarray(result,bcast=bcast)
  else: return result
#-------------------------------------------------------------------------
def getvz(iw=0,gather=1,bcast=None,**kw):
  "Returns the Z velocity."
  if bcast is None: bcast = _particlebcastdefault[0]
  ii = selectparticles(iw=iw,kwdict=kw)
  suffix,object,pgroup = _getobjectpgroup(kw)
  if (isinstance(ii,slice) and ii.stop > ii.start) or len(ii) > 0:
    uz = getattrwithsuffix(pgroup,'uzp',suffix)
    gaminv = getattrwithsuffix(pgroup,'gaminv',suffix)
  if isinstance(ii,slice) and ii.stop > ii.start:
    if top.lrelativ: result = uz[ii]*gaminv[ii]
    else:            result = uz[ii]
  elif len(ii) > 0:
    result = take(uz,ii)*take(gaminv,ii)
  else:
    result = array([],'d')
  if lparallel and gather: return gatherarray(result,bcast=bcast)
  else: return result
#-------------------------------------------------------------------------
def getvr(iw=0,gather=1,bcast=None,**kw):
  "Returns the radial velocity."
  if bcast is None: bcast = _particlebcastdefault[0]
  ii = selectparticles(iw=iw,kwdict=kw)
  suffix,object,pgroup = _getobjectpgroup(kw)
  if (isinstance(ii,slice) and ii.stop > ii.start) or len(ii) > 0:
    x = getattrwithsuffix(pgroup,'xp',suffix)
    y = getattrwithsuffix(pgroup,'yp',suffix)
    ux = getattrwithsuffix(pgroup,'uxp',suffix)
    uy = getattrwithsuffix(pgroup,'uyp',suffix)
    gaminv = getattrwithsuffix(pgroup,'gaminv',suffix)
  if isinstance(ii,slice) and ii.stop > ii.start:
    tt = arctan2(y[ii],x[ii])
    result = (ux[ii]*cos(tt) + uy[ii]*sin(tt))*gaminv[ii]
  elif len(ii) > 0:
    tt = arctan2(take(y,ii),take(x,ii))
    result = (take(ux,ii)*cos(tt) + take(uy,ii)*sin(tt))*take(gaminv,ii)
  else:
    result = array([],'d')
  if lparallel and gather: return gatherarray(result,bcast=bcast)
  else: return result
#-------------------------------------------------------------------------
def getvtheta(iw=0,gather=1,bcast=None,**kw):
  "Returns the azimuthal velocity."
  if bcast is None: bcast = _particlebcastdefault[0]
  ii = selectparticles(iw=iw,kwdict=kw)
  suffix,object,pgroup = _getobjectpgroup(kw)
  if (isinstance(ii,slice) and ii.stop > ii.start) or len(ii) > 0:
    x = getattrwithsuffix(pgroup,'xp',suffix)
    y = getattrwithsuffix(pgroup,'yp',suffix)
    ux = getattrwithsuffix(pgroup,'uxp',suffix)
    uy = getattrwithsuffix(pgroup,'uyp',suffix)
    gaminv = getattrwithsuffix(pgroup,'gaminv',suffix)
  if isinstance(ii,slice) and ii.stop > ii.start:
    tt = arctan2(y[ii],x[ii])
    result = (-ux[ii]*sin(tt) + uy[ii]*cos(tt))*gaminv[ii]
  elif len(ii) > 0:
    tt = arctan2(take(y,ii),take(x,ii))
    result = (-take(ux,ii)*sin(tt) + take(uy,ii)*cos(tt))*take(gaminv,ii)
  else:
    result = array([],'d')
  if lparallel and gather: return gatherarray(result,bcast=bcast)
  else: return result
#-------------------------------------------------------------------------
def getux(iw=0,gather=1,bcast=None,**kw):
  "Returns the X momentum over mass."
  if bcast is None: bcast = _particlebcastdefault[0]
  ii = selectparticles(iw=iw,kwdict=kw)
  suffix,object,pgroup = _getobjectpgroup(kw)
  if (isinstance(ii,slice) and ii.stop > ii.start) or len(ii) > 0:
    ux = getattrwithsuffix(pgroup,'uxp',suffix)
  if isinstance(ii,slice) and ii.stop > ii.start:
    result = ux[ii]
  elif len(ii) > 0:
    result = take(ux,ii)
  else:
    result = array([],'d')
  if lparallel and gather: return gatherarray(result,bcast=bcast)
  else: return result
#-------------------------------------------------------------------------
def getuy(iw=0,gather=1,bcast=None,**kw):
  "Returns the Y momentum over mass."
  if bcast is None: bcast = _particlebcastdefault[0]
  ii = selectparticles(iw=iw,kwdict=kw)
  suffix,object,pgroup = _getobjectpgroup(kw)
  if (isinstance(ii,slice) and ii.stop > ii.start) or len(ii) > 0:
    uy = getattrwithsuffix(pgroup,'uyp',suffix)
  if isinstance(ii,slice) and ii.stop > ii.start:
    result = uy[ii]
  elif len(ii) > 0:
    result = take(uy,ii)
  else:
    result = array([],'d')
  if lparallel and gather: return gatherarray(result,bcast=bcast)
  else: return result
#-------------------------------------------------------------------------
def getuz(iw=0,gather=1,bcast=None,**kw):
  "Returns the Z momentum over mass."
  if bcast is None: bcast = _particlebcastdefault[0]
  ii = selectparticles(iw=iw,kwdict=kw)
  suffix,object,pgroup = _getobjectpgroup(kw)
  if (isinstance(ii,slice) and ii.stop > ii.start) or len(ii) > 0:
    uz = getattrwithsuffix(pgroup,'uzp',suffix)
  if isinstance(ii,slice) and ii.stop > ii.start:
    result = uz[ii]
  elif len(ii) > 0:
    result = take(uz,ii)
  else:
    result = array([],'d')
  if lparallel and gather: return gatherarray(result,bcast=bcast)
  else: return result
#-------------------------------------------------------------------------
def getxp(iw=0,gather=1,bcast=None,**kw):
  "Returns the X velocity over the Z velocity (X')."
  if bcast is None: bcast = _particlebcastdefault[0]
  ii = selectparticles(iw=iw,kwdict=kw)
  suffix,object,pgroup = _getobjectpgroup(kw)
  if (isinstance(ii,slice) and ii.stop > ii.start) or len(ii) > 0:
    ux = getattrwithsuffix(pgroup,'uxp',suffix)
    uz = getattrwithsuffix(pgroup,'uzp',suffix)
  if isinstance(ii,slice) and ii.stop > ii.start:
    result = ux[ii]/uz[ii]
  elif len(ii) > 0:
    result = take(ux,ii)/take(uz,ii)
  else:
    result = array([],'d')
  if lparallel and gather: return gatherarray(result,bcast=bcast)
  else: return result
#-------------------------------------------------------------------------
def getyp(iw=0,gather=1,bcast=None,**kw):
  "Returns the Y velocity over the Z velocity (Y')."
  if bcast is None: bcast = _particlebcastdefault[0]
  ii = selectparticles(iw=iw,kwdict=kw)
  suffix,object,pgroup = _getobjectpgroup(kw)
  if (isinstance(ii,slice) and ii.stop > ii.start) or len(ii) > 0:
    uy = getattrwithsuffix(pgroup,'uyp',suffix)
    uz = getattrwithsuffix(pgroup,'uzp',suffix)
  if isinstance(ii,slice) and ii.stop > ii.start:
    result = uy[ii]/uz[ii]
  elif len(ii) > 0:
    result = take(uy,ii)/take(uz,ii)
  else:
    result = array([],'d')
  if lparallel and gather: return gatherarray(result,bcast=bcast)
  else: return result
#-------------------------------------------------------------------------
def getrp(iw=0,gather=1,bcast=None,**kw):
  "Returns the radial velocity over the Z velocity (R')."
  if bcast is None: bcast = _particlebcastdefault[0]
  ii = selectparticles(iw=iw,kwdict=kw)
  suffix,object,pgroup = _getobjectpgroup(kw)
  if (isinstance(ii,slice) and ii.stop > ii.start) or len(ii) > 0:
    x = getattrwithsuffix(pgroup,'xp',suffix)
    y = getattrwithsuffix(pgroup,'yp',suffix)
    ux = getattrwithsuffix(pgroup,'uxp',suffix)
    uy = getattrwithsuffix(pgroup,'uyp',suffix)
    uz = getattrwithsuffix(pgroup,'uzp',suffix)
  if isinstance(ii,slice) and ii.stop > ii.start:
    tt = arctan2(y[ii],x[ii])
    result = ((ux[ii]*cos(tt)+uy[ii]*sin(tt))/
              uz[ii])
  elif len(ii) > 0:
    tt = arctan2(take(y,ii),take(x,ii))
    result = ((take(ux,ii)*cos(tt)+take(uy,ii)*sin(tt))/
              take(uz,ii))
  else:
    result = array([],'d')
  if lparallel and gather: return gatherarray(result,bcast=bcast)
  else: return result
#-------------------------------------------------------------------------
def gettp(iw=0,gather=1,bcast=None,**kw):
  "Returns the azimuthal velocity over the Z velocity (R')."
  if bcast is None: bcast = _particlebcastdefault[0]
  ii = selectparticles(iw=iw,kwdict=kw)
  suffix,object,pgroup = _getobjectpgroup(kw)
  if (isinstance(ii,slice) and ii.stop > ii.start) or len(ii) > 0:
    x = getattrwithsuffix(pgroup,'xp',suffix)
    y = getattrwithsuffix(pgroup,'yp',suffix)
    ux = getattrwithsuffix(pgroup,'uxp',suffix)
    uy = getattrwithsuffix(pgroup,'uyp',suffix)
    uz = getattrwithsuffix(pgroup,'uzp',suffix)
  if isinstance(ii,slice) and ii.stop > ii.start:
    tt = arctan2(y[ii],x[ii])
    result = ((-ux[ii]*sin(tt)+uy[ii]*cos(tt))/
              uz[ii])
  elif len(ii) > 0:
    tt = arctan2(take(y,ii),take(x,ii))
    result = ((-take(ux,ii)*sin(tt)+take(uy,ii)*cos(tt))/
              take(uz,ii))
  else:
    result = array([],'d')
  if lparallel and gather: return gatherarray(result,bcast=bcast)
  else: return result
#-------------------------------------------------------------------------
def getgaminv(iw=0,gather=1,bcast=None,**kw):
  "Returns the gamma inverse."
  if bcast is None: bcast = _particlebcastdefault[0]
  ii = selectparticles(iw=iw,kwdict=kw)
  suffix,object,pgroup = _getobjectpgroup(kw)
  if (isinstance(ii,slice) and ii.stop > ii.start) or len(ii) > 0:
    gaminv = getattrwithsuffix(pgroup,'gaminv',suffix)
  if isinstance(ii,slice) and ii.stop > ii.start:
    result = gaminv[ii]
  elif len(ii) > 0:
    result = take(gaminv,ii)
  else:
    result = array([],'d')
  if lparallel and gather: return gatherarray(result,bcast=bcast)
  else: return result
#-------------------------------------------------------------------------
def getex(iw=0,gather=1,bcast=None,**kw):
  "Returns the Ex field applied to the particles"
  if bcast is None: bcast = _particlebcastdefault[0]
  ii = selectparticles(iw=iw,kwdict=kw)
  suffix,object,pgroup = _getobjectpgroup(kw)
  if (isinstance(ii,slice) and ii.stop > ii.start) or len(ii) > 0:
    ex = getattrwithsuffix(pgroup,'ex',suffix)
  if isinstance(ii,slice) and ii.stop > ii.start:
    result = ex[ii]
  elif len(ii) > 0:
    result = take(ex,ii)
  else:
    result = array([],'d')
  if lparallel and gather: return gatherarray(result,bcast=bcast)
  else: return result
#-------------------------------------------------------------------------
def getey(iw=0,gather=1,bcast=None,**kw):
  "Returns the Ey field applied to the particles"
  if bcast is None: bcast = _particlebcastdefault[0]
  ii = selectparticles(iw=iw,kwdict=kw)
  suffix,object,pgroup = _getobjectpgroup(kw)
  if (isinstance(ii,slice) and ii.stop > ii.start) or len(ii) > 0:
    ey = getattrwithsuffix(pgroup,'ey',suffix)
  if isinstance(ii,slice) and ii.stop > ii.start:
    result = ey[ii]
  elif len(ii) > 0:
    result = take(ey,ii)
  else:
    result = array([],'d')
  if lparallel and gather: return gatherarray(result,bcast=bcast)
  else: return result
#-------------------------------------------------------------------------
def getez(iw=0,gather=1,bcast=None,**kw):
  "Returns the Ez field applied to the particles"
  if bcast is None: bcast = _particlebcastdefault[0]
  ii = selectparticles(iw=iw,kwdict=kw)
  suffix,object,pgroup = _getobjectpgroup(kw)
  if (isinstance(ii,slice) and ii.stop > ii.start) or len(ii) > 0:
    ez = getattrwithsuffix(pgroup,'ez',suffix)
  if isinstance(ii,slice) and ii.stop > ii.start:
    result = ez[ii]
  elif len(ii) > 0:
    result = take(ez,ii)
  else:
    result = array([],'d')
  if lparallel and gather: return gatherarray(result,bcast=bcast)
  else: return result
#-------------------------------------------------------------------------
def geter(iw=0,gather=1,bcast=None,**kw):
  "Returns the Er field applied to the particles"
  if bcast is None: bcast = _particlebcastdefault[0]
  ii = selectparticles(iw=iw,kwdict=kw)
  suffix,object,pgroup = _getobjectpgroup(kw)
  if (isinstance(ii,slice) and ii.stop > ii.start) or len(ii) > 0:
    x = getattrwithsuffix(pgroup,'xp',suffix)
    y = getattrwithsuffix(pgroup,'yp',suffix)
    ex = getattrwithsuffix(pgroup,'ex',suffix)
    ey = getattrwithsuffix(pgroup,'ey',suffix)
  if isinstance(ii,slice) and ii.stop > ii.start:
    theta = arctan2(y[ii],x[ii])
    result = ex[ii]*cos(theta) + ey[ii]*sin(theta)
  elif len(ii) > 0:
    theta = arctan2(take(y,ii),take(x,ii))
    result = take(ex,ii)*cos(theta) + take(ey,ii)*sin(theta)
  else:
    result = array([],'d')
  if lparallel and gather: return gatherarray(result,bcast=bcast)
  else: return result
#-------------------------------------------------------------------------
def getetheta(iw=0,gather=1,bcast=None,**kw):
  "Returns the Etheta field applied to the particles"
  if bcast is None: bcast = _particlebcastdefault[0]
  ii = selectparticles(iw=iw,kwdict=kw)
  suffix,object,pgroup = _getobjectpgroup(kw)
  if (isinstance(ii,slice) and ii.stop > ii.start) or len(ii) > 0:
    x = getattrwithsuffix(pgroup,'xp',suffix)
    y = getattrwithsuffix(pgroup,'yp',suffix)
    ex = getattrwithsuffix(pgroup,'ex',suffix)
    ey = getattrwithsuffix(pgroup,'ey',suffix)
  if isinstance(ii,slice) and ii.stop > ii.start:
    theta = arctan2(y[ii],x[ii])
    result = -ex[ii]*sin(theta) + ey[ii]*cos(theta)
  elif len(ii) > 0:
    theta = arctan2(take(y,ii),take(x,ii))
    result = -take(ex,ii)*sin(theta) + take(ey,ii)*cos(theta)
  else:
    result = array([],'d')
  if lparallel and gather: return gatherarray(result,bcast=bcast)
  else: return result
#-------------------------------------------------------------------------
def getbx(iw=0,gather=1,bcast=None,**kw):
  "Returns the Bx field applied to the particles"
  if bcast is None: bcast = _particlebcastdefault[0]
  ii = selectparticles(iw=iw,kwdict=kw)
  suffix,object,pgroup = _getobjectpgroup(kw)
  if (isinstance(ii,slice) and ii.stop > ii.start) or len(ii) > 0:
    bx = getattrwithsuffix(pgroup,'bx',suffix)
  if isinstance(ii,slice) and ii.stop > ii.start:
    result = bx[ii]
  elif len(ii) > 0:
    result = take(bx,ii)
  else:
    result = array([],'d')
  if lparallel and gather: return gatherarray(result,bcast=bcast)
  else: return result
#-------------------------------------------------------------------------
def getby(iw=0,gather=1,bcast=None,**kw):
  "Returns the By field applied to the particles"
  if bcast is None: bcast = _particlebcastdefault[0]
  ii = selectparticles(iw=iw,kwdict=kw)
  suffix,object,pgroup = _getobjectpgroup(kw)
  if (isinstance(ii,slice) and ii.stop > ii.start) or len(ii) > 0:
    by = getattrwithsuffix(pgroup,'by',suffix)
  if isinstance(ii,slice) and ii.stop > ii.start:
    result = by[ii]
  elif len(ii) > 0:
    result = take(by,ii)
  else:
    result = array([],'d')
  if lparallel and gather: return gatherarray(result,bcast=bcast)
  else: return result
#-------------------------------------------------------------------------
def getbz(iw=0,gather=1,bcast=None,**kw):
  "Returns the Bz field applied to the particles"
  if bcast is None: bcast = _particlebcastdefault[0]
  ii = selectparticles(iw=iw,kwdict=kw)
  suffix,object,pgroup = _getobjectpgroup(kw)
  if (isinstance(ii,slice) and ii.stop > ii.start) or len(ii) > 0:
    bz = getattrwithsuffix(pgroup,'bz',suffix)
  if isinstance(ii,slice) and ii.stop > ii.start:
    result = bz[ii]
  elif len(ii) > 0:
    result = take(bz,ii)
  else:
    result = array([],'d')
  if lparallel and gather: return gatherarray(result,bcast=bcast)
  else: return result
#-------------------------------------------------------------------------
def getpid(iw=0,id=0,gather=1,bcast=None,**kw):
  """Returns particle id number.
  -id=0: which pid value to return
         if id=-1, returns all pids.
  """
  if bcast is None: bcast = _particlebcastdefault[0]
  suffix,object,pgroup = _getobjectpgroup(kw)
  lost = kw.get('lost',0)
  if lost:
    npidlost = getattrwithsuffix(object,'npidlost')
    dopid = (npidlost > 0)
  else:
    npid = getattrwithsuffix(pgroup,'npid',suffix)
    dopid = (npid > 0)
  if dopid:
    ii = selectparticles(iw=iw,kwdict=kw)
    if (isinstance(ii,slice) and ii.stop > ii.start) or len(ii) > 0:
      pid = getattrwithsuffix(pgroup,'pid',suffix)
    if isinstance(ii,slice) and ii.stop > ii.start:
      if id >= 0: result = pid[ii,id]
      else:       result = pid[ii,:]
    elif len(ii) > 0:
      if id >= 0: result = take(pid[:,id],ii)
      else:       result = take(pid[:,:],ii)
    else:
      result = array([],'d')
  else:
    result = array([],'d')
  if lparallel and gather: return gatherarray(result,bcast=bcast)
  else: return result
#-------------------------------------------------------------------------
def getvdrifts(iw=0,js=0,jslist=None,gather=1,bcast=None,edrift=1,bdrift=1,**kw):
  "Returns the velocity drifts."
  if bcast is None: bcast = _particlebcastdefault[0]
  if jslist is None:
    jslist=[js]
  nptot=0
  vxl=[]
  vyl=[]
  vzl=[]
  suffix,object,pgroup = _getobjectpgroup(kw)
  lost = kw.get('lost',0)
  for js in jslist:
    ii = selectparticles(iw=iw,js=js,jslist=None,kwdict=kw)
    if isinstance(ii,slice) and ii.stop <= ii.start: continue
    if not isinstance(ii,slice) and len(ii) == 0: continue
    if isinstance(ii,slice): ii = arange(ii.start,ii.stop)
    np=len(ii)
    nptot+=np
    x  = take(getattrwithsuffix(object,'xp', suffix),ii)
    y  = take(getattrwithsuffix(object,'yp', suffix),ii)
    z  = take(getattrwithsuffix(object,'zp', suffix),ii)
    ux = take(getattrwithsuffix(object,'uxp',suffix),ii)
    uy = take(getattrwithsuffix(object,'uyp',suffix),ii)
    uz = take(getattrwithsuffix(object,'uzp',suffix),ii)
    npinttmp=w3d.npint+0
    npfieldtmp=w3d.npfield+0
    w3d.npint=np
    w3d.npfield=np
    gchange('DKInterptmp')
    setvdrifts(pgroup,np,js+1,x,y,z,ux,uy,uz,"corrector")
    if edrift:
      if bdrift:
        vxl.append(w3d.vdx[:np])
        vyl.append(w3d.vdy[:np])
        vzl.append(w3d.vdz[:np])
      else:
        vxl.append(w3d.vex[:np])
        vyl.append(w3d.vey[:np])
        vzl.append(w3d.vez[:np])
    else:
      if bdrift:
        vxl.append(w3d.vbx[:np])
        vyl.append(w3d.vby[:np])
        vzl.append(w3d.vbz[:np])
    w3d.npint=npinttmp
    w3d.npfield=npfieldtmp
    gchange('DKInterptmp')
  if nptot == 0:
    return array([],'d'),array([],'d'),array([],'d')
  else:
    vx = zeros(nptot,'d')
    vy = zeros(nptot,'d')
    vz = zeros(nptot,'d')
    ip=0
    for il in range(len(vxl)):
      np=shape(vxl[il])
      vx[ip:ip+np]=vxl[il]
      vy[ip:ip+np]=vyl[il]
      vz[ip:ip+np]=vzl[il]
      ip+=np
  if lparallel and gather: 
    return gatherarray(vx,bcast=bcast), \
           gatherarray(vy,bcast=bcast), \
           gatherarray(vz,bcast=bcast)
  else:
    return vx,vy,vz
#-------------------------------------------------------------------------
# Add the selectparticles documentation to each of the routines.
if lparallel:
  _gatherdoc = (
"""  gather=1 When 1, all data is gathered to PE0
  bcast=1: when true, result is broadcast to all processors
           Note that the setgetparticlebcastdefault can be
           used to set the default value to 0, which is
           more efficient if the bcast is not needed.
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
getex.__doc__ = getex.__doc__ + selectparticles.__doc__ + _gatherdoc
getey.__doc__ = getey.__doc__ + selectparticles.__doc__ + _gatherdoc
getez.__doc__ = getez.__doc__ + selectparticles.__doc__ + _gatherdoc
geter.__doc__ = geter.__doc__ + selectparticles.__doc__ + _gatherdoc
getetheta.__doc__ = getetheta.__doc__ + selectparticles.__doc__ + _gatherdoc
getbx.__doc__ = getbx.__doc__ + selectparticles.__doc__ + _gatherdoc
getby.__doc__ = getby.__doc__ + selectparticles.__doc__ + _gatherdoc
getbz.__doc__ = getbz.__doc__ + selectparticles.__doc__ + _gatherdoc
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
def addparticles(x=0.,y=0.,z=0.,vx=0.,vy=0.,vz=0.,gi=1.,pid=0.,w=1.,js=0,sid=None,
                 lallindomain=None,zmmin=None,zmmax=None,lmomentum=false,
                 resetrho=false,dofieldsol=false,resetmoments=false,
                 pgroup=None):
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
  pgroup=top.pgroup: Particle group to add particles too
  """

  # --- Check if this is a new species
  if js+1 > top.ns: setnspecies(js+1,pgroup)

  # --- Set the sid if it hasn't already be done
  if sid is None: sid = js
  if top.pgroup.sid[js] == -1: top.pgroup.sid[js] = sid

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

  # --- If any of the inputs are arrays that are zero length, then return
  if (lenx == 0 or leny == 0 or lenz == 0 or
      lenvx == 0 or lenvy == 0 or lenvz == 0 or
      lengi == 0 or lenpid == 0): return

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
  w = array(w)*ones(maxlen,'d')

  # --- Set time of creation
  if top.tpid>0: pid[:,top.tpid-1]=top.time
  # --- Set weights
  if top.wpid>0: pid[:,top.wpid-1]=w
  # --- Note that ssn is set in addpart

  # --- Set extent of domain
  if not lparallel:
    if zmmin is None: zmmin = w3d.zmmin + top.zbeam
    if zmmax is None: zmmax = w3d.zmmax + top.zbeam
  else:
    if zmmin is None: zmmin = top.zpminlocal + top.zbeam
    if zmmax is None: zmmax = top.zpmaxlocal + top.zbeam

  # --- When running in slice mode, automatically set lallindomain to true.
  # --- It is assumed that all particles will be within the specified domain,
  # --- since in the slice mode, the z of the particles is ignored.
  # --- The user can still set lallindomain to false to override this.
  if lallindomain is None:
    if getcurrpkg() == 'wxy': lallindomain = true
    else:                     lallindomain = false

  # --- Do some error checking
  if not lallindomain and zmmin == zmmax:
    print "=================================================================="
    print "Addparticles: warning - no particles will be loaded - you should"
    print "either set lallindomain=true or set zmmin and zmmax so they are"
    print "different from each other."
    print "=================================================================="

  if pgroup is None: pgroup = top.pgroup

  # --- Now data can be passed into the fortran addparticles routine.
  addpart(pgroup,maxlen,top.npid,x,y,z,vx,vy,vz,gi,pid,js+1,
          lallindomain,zmmin,zmmax,lmomentum)
 
  # --- If the slice code is active, then call initdtp
  if package()[0] == 'wxy': initdtp(top.pgroup)

  # --- Do followup work if requested
  if resetrho:
    loadrho()
  if dofieldsol:
    fieldsol(-1)
  if resetmoments:
    import getzmom
    getzmom.zmmnt()
    if top.it%top.nhist == 0:
      top.jhist = top.jhist - 1
      savehist(top.time)

    
