warp_version = "$Id: warp.py,v 1.59 2003/09/18 18:41:47 dave Exp $"
# import all of the neccesary packages
import __main__
from Numeric import *
import sys
import os.path

# --- Import the RNG module. Older versions have ranf in a seperate module
# --- called Ranf. In newer versions, ranf is part of RNG.
try:
  import Ranf
except ImportError:
  pass
try:
  import RNG
except ImportError:
  pass

# --- Gist needs to be imported before pyBasis since pyBasis calls a function
# --- from gist. Also, since gist is only loaded on PE0 in the parallel
# --- environment, the parallel module must also be loaded in to check
# --- whether this is in fact running in parallel.

# --- This creates a logical lparallel which is true when running in parallel
# --- and false when running serial.
# --- This also creates a number of routines which are needed for parallel
# --- data handling. These routines all have sensible default results when
# --- not running in parallel.
from parallel import *

try:
  if me == 0 and sys.platform != 'mac':
    from gist import *
  else:
    from gistdummy import *
except ImportError:
  pass

# Try importing the warpC shared object which contains all of WARP and pybasisC
try:
  from warpC import *
except ImportError:
  pass

from pyBasis import *

# --- The WARP modules must be imported in the order below because of
# --- linking dependencies.
from toppy import *
from envpy import *
from f3dpy import *
from w3dpy import *
from fxypy import *
from wxypy import *
try:  # RZ code hasn't been installed on all machines yet
  from frzpy import *
  from wrzpy import *
except ImportError:
  pass
try:  # cirpy hasn't been installed on all machines yet
  from cirpy import *
except ImportError:
  pass
try:  # herpy hasn't been installed on all machines yet
  from herpy import *
except ImportError:
  pass
try:  # chopy hasn't been installed on all machines yet
  from chopy import *
except ImportError:
  pass
from ctl import *

# --- Rearrange the list of packages to a more sensible order
package('wxy')
package('w3d')
package('top')

# --- If running in parallel, import the parallel warp module
if lparallel:
  from warpparallel import *

# --- Add stuff to the path
import sys

# --- Set default runid to first filename in the command line, stripping off
# --- the .py suffix.
if sys.argv[0]:
  if sys.argv[0][-3:] == '.py':
    h,t = os.path.split(sys.argv[0][:-3])
    top.runid = t
    del h,t
  else:
    h,t = os.path.split(sys.argv[0])
    top.runid = t
    del h,t
runid = arraytostr(top.runid)

#=============================================================================
# --- Set physical constants which depend on others.
# --- Magnetic constant = 4*pi*1.e-7
top.mu0 = 4*top.pi/10000000
# --- Conversion factor from joules to eV is just echarge
top.jperev = top.echarge
# --- Epsilon_0 calculated from speed of light and mu_0
top.eps0 = 1/(top.mu0*top.clight*top.clight)
# --- Create python versions of the constants
amu       = top.amu
clight    = top.clight
echarge   = top.echarge
emass     = top.emass
eps0      = top.eps0
euler     = top.euler
jperev    = top.jperev
mu0       = top.mu0
boltzmann = top.boltzmann
largepos  = top.largepos
smallpos  = top.smallpos
try:
  dirichlet = top.dirichlet
  neumann   = top.neumann
  periodic  = top.periodic
  absorb    = top.absorb
  reflect   = top.reflect
except AttrbuteError:
  dirichlet = 0
  neumann   = 1
  periodic  = 2
  absorb    = 0
  reflect   = 1

# --- Simple function to calculate Child-Langmuir current density
def childlangmuir(v,d,q=None,m=None):
  if q is None: q = top.sq[0]
  if m is None: m = top.sm[0]
  return 4./9.*eps0*sqrt(2.*q/m)*v**1.5/d**2

# --- Create 


# --- Create python version of dvnz (divisor not zero)
def dvnz(x):
  if x == 0.: return top.smallpos
  else:       return x
  #return sign(abs(x)+top.smallpos,x)

#=============================================================================
# --- Setup and make initial printout of the versions of the packages.
def printversion(v):
  v = arraytostr(v,strip=false)
  lenv = 12
  for i in range(11,19):
    if v[i] == "$": lenv = i
  return v[11:lenv-1]

def versionstext():
  "Returns a string which has the version information of packages loaded."
  r = 'Python WARP\n'
  pkg = package()
  fmt = '******  %s version %s\n'
  if 'fxy' in pkg: r=r+fmt%('Fieldsolver FXY',printversion(fxy.versfxy))
  if 'frz' in pkg: r=r+fmt%('Fieldsolver FRZ',printversion(frz.versfrz))
  if 'f3d' in pkg: r=r+fmt%('Fieldsolver F3D',printversion(f3d.versf3d))
  if 'env' in pkg: r=r+fmt%('Envelope solver ENV',printversion(env.versenv))
  if 'cir' in pkg: r=r+fmt%('Envelope solver CIR',printversion(cir.verscir))
  if 'her' in pkg: r=r+fmt%('Envelope solver HER',printversion(her.versher))
  if 'wxy' in pkg: r=r+fmt%('Particle package WXY',printversion(wxy.verswxy))
  if 'wrz' in pkg: r=r+fmt%('Particle package WRZ',printversion(wrz.verswrz))
  if 'w3d' in pkg: r=r+fmt%('Particle package W3D',printversion(w3d.versw3d))
  if 'top' in pkg: r=r+fmt%('Main package TOP',printversion(top.verstop))
  return r
print versionstext()[:-1] # --- skip last line feed
print 'For more help, type warphelp()'

#=============================================================================
# --- Declare the documentation for the warp module.
def warpdoc():
  print """
Imports the basic modules needed to run WARP, including
Numeric, gist, warpplots

as well as additional modules
histplots, pzplots, plot_conductors, drawlattice

Create python versions of the constants
amu, clight, echarge, emass, eps0, euler, jperev, mu0

Defines following functions...
dump: Creates a dump file containing the current state of the simulation
restart: Retreives the state of a simulation from a dump file
loadrho: Load the particles onto the charge density mesh
fieldsol: Solve for the self-fields
installbeforefs: Install a function which will be called before a field-solve
uninstallbeforefs: Uninstall the function called before a field-solve
isinstalledbeforefs: Checks if a function is installed to be called before a
                     field-solve
installafterfs: Install a function which will be called after a field-solve
isinstalledafterfs: Checks if a function is installed to be called after a
                    field-solve
uninstallafterfs: Uninstall the function called after a field-solve
installbeforestep: Install a function which will be called before a step
uninstallbeforestep: Uninstall the function called before a step
isinstalledbeforefs: Checks if a function is installed to be called before a
                     step
installafterstep: Install a function which will be called after a step
uninstallafterstep: Uninstall the function called after a step
isinstalledafterfs: Checks if a function is installed to be called after a
                    step
installparticlescraper: Installs a function which will be called at the
                        correct time to scrape particles
uninstallparticlescraper: Uninstalls a function which will be called at the
                          correct time to scrape particles
isinstalledparticlescraper: Checks if a function is installed to be called at
                            the correct time to scrape particles
installaddconductor: Installs a function which will be called at the beginning
                    of the field solve so conductors will be added.
uninstalladdconductor: Uninstalls the function which would be called at the
                       beginning of the field solve so conductors will be added.
isinstalledaddconductor: Checks if a function is installed to be called at
                         the beginning of the field solve so conductors will
                         be added.
gethzarrays: Fixes the ordering of hlinechg and hvzofz data from a paralle run
printtimers: Print timers in a nice annotated format
  """

#=============================================================================
# --- Call derivqty to calculate eps0 and jperev
derivqty()

#=============================================================================
# --- Setup mechanism for "before" and "after" python scripts
beforefsfuncs = []
afterfsfuncs = []
callscraperfuncs = []
addconductorfuncs = []
__controlfuncs = {'beforefs':beforefsfuncs,
                  'afterfs':afterfsfuncs,
                  'callscraper':callscraperfuncs,
                  'addconductor':addconductorfuncs,
                  'beforestep':beforestepfuncs,
                  'afterstep':afterstepfuncs}
def beforefs():
  for f in beforefsfuncs: f()
def afterfs():
  for f in afterfsfuncs: f()
def callscraper():
  for f in callscraperfuncs: f()
def calladdconductor():
  for f in addconductorfuncs: f()

def installbeforefs(f):
  "Adds a function to the list of functions called before a field-solve"
  beforefsfuncs.append(f)
  w3d.lbeforefs = true
def uninstallbeforefs(f):
  "Removes the function from the list of functions called before a field-solve"
  if f in beforefsfuncs:
    beforefsfuncs.remove(f)
  else:
    raise 'Warning: uninstallbeforefs: no such function had been installed'
  if len(beforefsfuncs) == 0: w3d.lbeforefs = false
def isinstalledbeforefs(f):
  return f in beforefsfuncs

def installafterfs(f):
  "Adds a function to the list of functions called after a field-solve"
  afterfsfuncs.append(f)
  w3d.lafterfs = true
def uninstallafterfs(f):
  "Removes the function from the list of functions called after a field-solve"
  if f in afterfsfuncs:
    afterfsfuncs.remove(f)
  else:
    raise 'Warning: uninstallafterfs: no such function had been installed'
  if len(afterfsfuncs) == 0: w3d.lafterfs = false
def isinstalledafterfs(f):
  return f in afterfsfuncs

def installparticlescraper(f):
  "Adds a function to the list of functions called to scrape particles"
  callscraperfuncs.append(f)
  w3d.lcallscraper = true
def uninstallparticlescraper(f):
  "Removes the function from the list of functions called to scrape particles"
  if f in callscraperfuncs:
    callscraperfuncs.remove(f)
  else:
    raise 'Warning: uninstallparticlescraper: no such function had been installed'
  if len(callscraperfuncs) == 0: w3d.lcallscraper = false
def isinstalledparticlescraper(f):
  return f in callscraperfuncs

def installaddconductor(f):
  "Adds a function to the list of functions called to add conductors"
  addconductorfuncs.append(f)
  f3d.laddconductor = true
def uninstalladdconductor(f):
  "Removes the function from the list of functions called to add conductors"
  if f in addconductorfuncs:
    addconductorfuncs.remove(f)
  else:
    raise 'Warning: uninstalladdconductor: no such function had been installed'
  if len(addconductorfuncs) == 0: f3d.laddconductor = false
def isinstalledaddconductor(f):
  return f in addconductorfuncs

def installbeforestep(f):
  "Adds a function to the list of functions called before a step"
  beforestepfuncs.append(f)
def uninstallbeforestep(f):
  "Removes the function from the list of functions called before a step"
  if f in beforestepfuncs:
    beforestepfuncs.remove(f)
  else:
    raise 'Warning: uninstallbeforestep: no such function had been installed'
def isinstalledbeforestep(f):
  return f in beforestepfuncs

def installafterstep(f):
  "Adds a function to the list of functions called after a step"
  afterstepfuncs.append(f)
def uninstallafterstep(f):
  "Removes the function from the list of functions called after a step"
  if f in afterstepfuncs:
    afterstepfuncs.remove(f)
  else:
    raise 'Warning: uninstallafterstep: no such function had been installed'
def isinstalledafterstep(f):
  return f in afterstepfuncs

#=============================================================================
# --- Convenience function for random numbers.
def setseed(x=0,y=0):
  RandomArray.seed(x,y)
    
# --- Uniform distribution
def ranf(x=None,i1=None,nbase=None):
  """Returns a pseudo-random number. If the 2nd and 3rd arguments are given,
returns a digit reversed random number.
  x=None: if present, returns and array the same shape fill with random numbers
  i1=None: returns i'th digit reversed number
  nbase=None: base to use for digit reversing
  """
  if not i1:
    if x is None:
      # --- Try various versions of ranf.
      # --- RandomArray only returns single precision
      # --- Ranf is left over from original version of Numeric
      # --- RNG is best and most recent. Returns double precision
      try: return RNG.ranf()
      except: pass
      try: return Ranf.ranf()
      except: pass
      try: return RandomArray.random(shape(array([0])))[0]
      except: raise "No random number modules found"
    else:
      try: return apply(RNG.random_sample,shape(x))
      except: pass
      try: return apply(Ranf.random_sample,shape(x))
      except: pass
      try: return RandomArray.random(shape(x))
      except: raise "No random number modules found"
  else:
    if not nbase: nbase = 2
    n = product(array(shape(array(x))))
    result = zeros(n,'d')
    rnrevarray(n,result,i1,nbase)
    result.shape = shape(x)
    return result

# --- Gaussian distribution
# --- This was in pyBasis but had to be moved here in order to use rnormdig.
# --- First, try and define a normal generator from the RNG module.
try:
  _normaldistribution = RNG.NormalDistribution(0.,1.)
  _normalgenerator = RNG.CreateGenerator(-1,_normaldistribution)
except:
  _normalgenerator = None
def rnormarray(x,i1=None,nbase1=None,nbase2=None):
  if not i1:
    if _normalgenerator is not None:
      return reshape(_normalgenerator.sample(product(array(shape(x)))),shape(x))
    try:
      return RandomArray.standard_normal(x)
    except:
      # --- Use pseudo-random number generator
      s = RandomArray.random(shape(x))
      phi = 2.*pi*RandomArray.random(shape(x))
      sq = sqrt(-2.*log(s))
      return sq*cos(phi)
  else:
    # --- Use digit reversed random number generator
    if not nbase1: nbase1 = 2
    if not nbase2: nbase2 = 3
    n = product(array(shape(x)))
    result = zeros(n,'d')
    rnormdig(i1,n,nbase1,nbase2,0.,result)
    result.shape = shape(x)
    return result

#=============================================================================
def loadrho(ins_i=-1,nps_i=-1,is_i=-1,lzero=true):
  """
loadrho(ins_i=-1,nps_i=-1,is_i=-1,lzero=1)
This routine provides a simple call from the interpreter to load the
rho array.  All of the arguments are optional.
If the species is not specified, all species are loaded, except
when ins or nps are specified, then only species 1 is loaded.
lzero is used to set whether or not rho is zeroed before the load.
The default is to zero out rho.
  """

  # --- if particle location is specified but species is not, set so
  # --- only species number 1 is included
  if (ins_i != -1 and is_i == -1): is_i = 1

  # --- set number of particles
  if (ins_i != -1 and nps_i == -1):
    # --- if particle number is omitted but particle location is specified,
    # --- set nps to get rest of active particles of species
    nps_i = top.nps[is_i] + top.ins[is_i] - ins_i

  # --- if particle number is specified but species is not, set so
  # --- only species number 1 is included
  if (nps_i != -1 and is_i == -1): is_i = 1

  # --- Now call the appropriate compiled interface routine based on the
  # --- current package
  currpkg = package()[0]
  if (currpkg == "w3d"):
    loadrho3d(ins_i,nps_i,is_i,lzero)
  elif (currpkg == "wxy"):
    loadrhoxy(ins_i,nps_i,is_i,lzero)
#=============================================================================
def fieldsol(iwhich=0,lbeforefs=false,lafterfs=false):
  """
This routine provides a simple call from the interpreter to do the fieldsol.
It calls the appropriate compiled interface routine based on the current
package. Only w3d and wxy have field solves defined.
 - iwhich=0: specifies what action to take
 - lbeforefs=false: when true, call functions installed be installbeforefs
 - lafterfs=false:  when true, call functions installed be installafterfs
  """
  if lbeforefs: beforefs()
  currpkg = package()[0]
  if   (currpkg == "w3d"): fieldsol3d(iwhich)
  elif (currpkg == "wxy"): fieldsolxy(iwhich)
  if lafterfs: afterfs()
  # --- Now do extra work, updating arrays which depend directly on phi,
  # --- but only when a complete field solve was done.
  if iwhich == -1 or iwhich == 0:
    if top.efetch == 3:
      getselfe3d(w3d.phi,w3d.nx,w3d.ny,w3d.nz,w3d.selfe,
                 w3d.nx_selfe,w3d.ny_selfe,w3d.nz_selfe,
                 w3d.dx,w3d.dy,w3d.dz,top.pboundxy)
    if top.inject > 0:
      try:
        # --- This routine is not defined in pywarp77
        getinj_phi()
      except NameError:
        pass

##############################################################################
##############################################################################
##############################################################################
##############################################################################
##############################################################################
##############################################################################
##############################################################################
##############################################################################
# --- Dump command
def dump(filename=None,suffix='',attr='dump',serial=0,onefile=1,pyvars=1,
         ff=None,varsuffix=None,histz=2,resizeHist=1,verbose=false,hdf=0):
  """
Creates a dump file
  - filename=(runid+'%06d'%top.it+suffix+'.dump')
  - attr='dump': All variables with the given attribute or group name are
    written to the file. The default attribute makes a restartable dump file.
  - serial=0: When 1, does a dump of only non-parallel data (parallel version
    only).
  - onefile=1: When 1, all processors dump to one file, otherwise each dumps to
    seperate file. The processor number is appended to the dump filename.
  - pyvars=1: When 1, saves user defined python variables to the file.
  - ff=None: Optional file object. When passed in, write to that file instead
             of creating a new one. Note that the inputted file object must be
             closed by the user.
  - varsuffix=None: Suffix to add to the variable names. If none is specified,
                    the suffix '@pkg' is used, where pkg is the package name
                    that the variable is in.
  - resizeHist=1: When true, resize history arrays so that unused locations
                  are removed.
  - hdf=0: when true, dump into an HDF file rather than a PDB.
  """
  timetemp = wtime()
  if not filename:
    # --- Setup default filename based on time step and processor number.
    if serial:
      s = '.sdump'
    else:
      s = '.dump'
    if onefile or not lparallel:
      filename=arraytostr(top.runid)+('%06d'%top.it)+suffix+s
    else:
      filename=arraytostr(top.runid)+('%06d_%05d'%(top.it,me))+suffix+s
  else:
    if not onefile and lparallel:
      # --- Append the processor number to the user inputted filename
      filename = filename + '%05d'%me
  # --- Make list of all of the new python variables.
  interpreter_variables = []
  if pyvars:
    if attr == 'dump' or 'dump' in attr:
      # --- Convert control functions to their names so they can be written
      # --- Note, only functions are converted, not methods of class
      # --- instances. Also, note that functions defined interactively will
      # --- not be restored properly since the source is not available.
      for n,l in __controlfuncs.iteritems():
        __main__.__dict__['control'+n] = []
        for f in l:
          if type(f) is FunctionType:
            __main__.__dict__['control'+n].append(f.__name__)
    # --- Add to the list all variables which are not in the initial list
    for l in __main__.__dict__.keys():
      if l not in initial_global_dict_keys:
        interpreter_variables.append(l)
  # --- Resize history arrays if requested.
  if resizeHist:
    top.lenhist = top.jhist
    gchange("Hist")
  # --- Call routine to make data dump
  if onefile and lparallel:
    paralleldump(filename,attr,interpreter_variables,serial=serial,
                 varsuffix=varsuffix,histz=histz,verbose=verbose)
  else:
    pydump(filename,attr,interpreter_variables,serial=serial,ff=ff,
           varsuffix=varsuffix,verbose=verbose,hdf=hdf)
  # --- Remove control names from main dict
  for n in __controlfuncs.iterkeys():
    if 'control'+n in __main__.__dict__:
      del __main__.__dict__['control'+n]
  # --- Update dump time
  top.dumptime = top.dumptime + (wtime() - timetemp)

# --- Restart command
def restart(filename,onefile=1,verbose=false,dofieldsol=true):
  """
Reads in data from file, redeposits charge density and does field solve
  - filename: restart file name - when restoring parallel run from multiple
              files, filename should only be prefix up to but not including
              the '_' before the processor number.
  - onefile=1: Restores from one file unless 0, then each processor restores
               from seperate file.
  - dofieldsol=true: When true, call fieldsol(0). This allows special cases
                     where just calling fieldsol(0) is not appropriate or
                     optimal
  """
  # --- If each processor is restoring from a seperate file, append
  # --- appropriate suffix, assuming only prefix was passed in
  if not onefile:
    filename = filename + '_%05d.dump'%(me)
  # --- Call different restore routine depending on context
  if onefile and lparallel:
    parallelrestore(filename,verbose=verbose)
  else:
    pyrestore(filename,verbose=verbose)
  # --- Recreate control function lists
  for n,l in __controlfuncs.iteritems():
    if 'control'+n not in __main__.__dict__: continue
    # --- Add items that are defined to the list if not already there.
    for fname in __main__.__dict__['control'+n]:
      if fname in __main__.__dict__:
        if not __main__.__dict__['isinstalled'+n](__main__.__dict__[fname]):
          __main__.__dict__['install'+n](__main__.__dict__[fname])
    # --- Delete the temporary variable
    del __main__.__dict__['control'+n]
  # --- Now that the dump file has been read in, finish up the restart work.
  # --- First set the current packge. Note that currpkg is only ever defined
  # --- in the main dictionary.
  package(__main__.__dict__["currpkg"])
  # --- Allocate all arrays appropriately
  gchange("*")
  # --- Reinitialize some injection stuff if it is needed.
  # --- This is really only needed for the parallel version since some of the
  # --- data saved is only valid for PE0.
  if top.inject > 0: fill_inj(w3d.dx,w3d.dy,w3d.dz,w3d.ix_axis,w3d.iy_axis)
  # --- Set the lattice internal variables. Only needed if reading in a dump
  # --- that was made before the overlapping elements was implemented.
  # --- Otherwise is doesn't hurt anything.
  resetlat()
  setlatt()
  # --- Load the charge density (since it was not saved)
  loadrho()
  # --- Recalculate the fields (since it was not saved)
  if dofieldsol: fieldsol(0)
  # --- Call setup if it is needed.
  if me == 0 and current_window() == -1: setup()

##############################################################################
##############################################################################
##############################################################################
##############################################################################
##############################################################################
##############################################################################
##############################################################################
##############################################################################

# This routine reads in and rearranges the history arrays from
# a parallel run so that they can be interpreted from a serial run.
def gethzarrays(filename,verbose=0):
  # --- Open file and get some numbers
  ff = PR.PR(filename)
  npes = len(ff.read('nps_p@parallel'))
  nz = shape(ff.read('zmesh@w3d'))[0] - 1
  lenhist = ff.read('lenhist@top')
  top.lenhist = lenhist
  top.nzmmnt = nz
  top.nzzarr = nz
  # --- Read in each array
  hlist = ['hlinechg','hvzofz','hepsxz','hepsyz','hepsnxz','hepsnyz','hepsgz',
    'hepshz','hepsngz','hepsnhz','hxbarz','hybarz','hxybarz','hxrmsz','hyrmsz',
    'hxprmsz','hyprmsz','hxsqbarz','hysqbarz','hvxbarz','hvybarz','hxpbarz',
    'hypbarz','hvxrmsz','hvyrmsz','hvzrmsz','hxpsqbarz','hypsqbarz','hxxpbarz',
    'hyypbarz','hxypbarz','hyxpbarz','hxpypbarz','hxvzbarz','hyvzbarz',
    'hvxvzbarz','hvyvzbarz']
  for h in hlist:
    try:
      exec('top.l%s = ff.read("l%s@top")'%(h,h))
      exec('top.i%s = ff.read("i%s@top")'%(h,h))
    except:
      break
  # --- Allocate space
  gchange("Hist")
  for h in hlist:
    # --- If true, then read in data
    if eval('top.l'+h):
      if verbose: print "Reading in ",h
      # --- Read in the whole array
      d = ff.read(h+'@parallel')
      v = eval('top.'+h)
      # --- Copy it to the correct place chunk by chunk
      for i in range(npes):
        v[i*(nz/npes):(i+1)*(nz/npes)+1,:] = d[i,:,:]
  ff.close()

##############################################################################
def printtimers(file=None):
  """Print timers in a nice annotated format
  - file Optional input file. If it is not include, stdout is used. It can
         either be a file name, or a file object. If a file name, a file
         with that name is created. If a file object, the data is written
         to that file (the file remains open).
  """
  if file is None:
    ff = sys.stdout
    closeit = 0
  elif type(file) == type(""):
    ff = open(file,"w")
    closeit = 1
  else:
    ff = file
    closeit = 0
  if not lparallel:
    ff.write('                 Total time')
    if top.it > 0: ff.write('          Time per step')
    ff.write('\n')
    ff.write('                          (s)')
    if top.it > 0: ff.write('                  (s)')
    ff.write('\n')
    ff.write('Generate time      %10.4f'%top.gentime)
    ff.write('\n')
    ff.write('Step time          %10.4f'%top.steptime)
    if top.it > 0: ff.write('           %10.4f'%(top.steptime/top.it))
    ff.write('\n')
    ff.write('Plot time          %10.4f'%top.plottime)
    if top.it > 0: ff.write('           %10.4f'%(top.plottime/(top.it+1)))
    ff.write('\n')
    ff.write('Moments time       %10.4f'%top.momentstime)
    if top.it > 0: ff.write('           %10.4f'%(top.momentstime/(top.it+1)))
    ff.write('\n')
    ff.write('Field Solve time   %10.4f'%top.fstime)
    if top.it > 0: ff.write('           %10.4f'%(top.fstime/(top.it+1)))
    ff.write('\n')
    ff.write('Applied field time %10.4f'%top.latticetime)
    if top.it > 0: ff.write('           %10.4f'%(top.latticetime/(top.it+1)))
    ff.write('\n')
    ff.write('Dump time          %10.4f'%top.dumptime)
    if top.it > 0: ff.write('           %10.4f'%(top.dumptime/top.it))
    ff.write('\n')
  else:
    timers = [['Generate time',      top.gentime              ],
              ['Step time',          top.steptime,    top.it  ],
              ['Plot time',          top.plottime,    top.it+1],
              ['Moments time',       top.momentstime, top.it+1],
              ['Field Solve time',   top.fstime,      top.it+1],
              ['Applied field time', top.latticetime, top.it+1],
              ['Dump time',          top.dumptime             ]]
    timelists = []
    totaltimes = []
    timedevs = []
    for t in timers: timelists.append(array(gather(t[1]))) 
    if me > 0: return
    for t in timelists: totaltimes.append(sum(t))
    for t in timelists: timedevs.append(sqrt(max(0.,ave(t**2) - ave(t)**2)))
    h1a = '                          Total time         Deviation'
    h2a = '                    (all CPUs)   (per CPU)            '
    h3a = '                        (s)         (s)         (s)   '
    h1b = '  Time per step'
    h2b = '    (per CPU)'
    h3b = '       (s)'
    f1a = '%18s  %10.4f  %10.4f  %10.4f'
    f1b = '   %10.4f'
    if top.it > 0:
       h1 = h1a + h1b + '\n'
       h2 = h2a + h2b + '\n'
       h3 = h3a + h3b + '\n'
       f1 = f1a + '\n'
       f2 = f1a + f1b + '\n'
       tt = []
       format = []
       for i in range(len(timelists)):
         if len(timers[i]) == 3:
           tt.append((timers[i][0],totaltimes[i],totaltimes[i]/npes,
                      timedevs[i],totaltimes[i]/npes/timers[i][2]))
           format.append(f2)
         else:
           tt.append((timers[i][0],totaltimes[i],totaltimes[i]/npes,
                      timedevs[i]))
           format.append(f1)
    else:
       h1 = h1a + '\n'
       h2 = h2a + '\n'
       h3 = h3a + '\n'
       f1 = f1a + '\n'
       f2 = f1a + '\n'
       tt = []
       format = []
       for i in range(len(timelists)):
         tt.append((timers[i][0],totaltimes[i],totaltimes[i]/npes,timedevs[i]))
         format.append(f1)
    ff.write(h1)
    ff.write(h2)
    ff.write(h3)
    for i in range(len(timelists)):
      ff.write(format[i]%tt[i])
  if closeit:
    ff.close()


#=============================================================================
# --- Import the convenience routines for plotting various slices and
# --- projections of particles, histories, as well as some line plots.
# --- Import these here near the end so the functions defined above are
# --- included in their dictionaries.
from particles import *
from warpplots import *
from histplots import *
from pzplots import *
from plot_conductor import *
from drawlattice import *

# --- Import some online documentation modules.
from warphelp import *
from warpscripts import *
from warpfortran import *

# --- Import the printparameters modules (which are called from fortran)
from printparameters import *
from printparameters3d import *
from printparametersrz import *

# --- Try to import the GUI
try:
  import os
  import warp
  sys.path=sys.path+[os.path.dirname(warp.__file__)+'/GUI']
  from WarpGUI import *
except:
  pass

##############################################################################
######  Don't put anything below this line!!! ################################
##############################################################################

# --- Save the initial keys in the global dictionary. This allows the pydump
# --- command to save interpreter variables (without saving huge amounts
# --- of stuff that is not needed). Note that initial_global_dict_keys is
# --- declared first as a empty list so that it itself appears in the list
# --- of global keys.
initial_global_dict_keys = []
initial_global_dict_keys = globals().keys()

# --- Get start CPU time
top.starttime = wtime()

# --- Save the versions string here so that it will be dumped into any
# --- dump file.
warpversions = versionstext()
