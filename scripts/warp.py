# import all of the neccesary packages
import __main__
from Numeric import *
import ranlib
warp_version = "$Id: warp.py,v 1.1 2000/10/16 18:34:19 dave Exp $"

# --- Gist needs to be imported before pyBasis since pyBasis calls a function
# --- from gist. Also, since gist is only loaded on PE0 in the parallel
# --- environment, the parallel module must also be loaded in to check
# --- whether this is in fact running in parallel.
try:
  from parallel import *
  lparallel = 1
except ImportError:
  # --- Create place holding variables. These need to be created so that
  # --- scripts can tell whether the code is running serially or in parallel.
  me = 0
  npes = 0
  lparallel = 0

if me == 0:
  from gist import *
else:
  from gistdummy import *

from pyBasis import *

# --- The WARP modules must be imported in the order below because of
# --- linking dependencies.
from toppy import *
from envpy import *
from f3dpy import *
from w3dpy import *
from fxypy import *
from wxypy import *
try:  # cirpy hasn't been installed on all machines yet
  from cirpy import *
except ImportError:
  pass
from ctl import *

# --- Rearrange the list of packages to a more sensible order
package('wxy')
package('w3d')
package('top')

# --- Import some extra packages which are used by various packages.
import PR
import PW
import re
import RandomArray

# --- If running in parallel, import the parallel warp module
if lparallel:
  from warpparallel import *

# --- Add stuff to the path
import sys
sys.path = sys.path + ['/home/ife1/dave/warp/scripts']

# --- Set physical constants which depend on others.
# --- Magnetic constant = 4*pi*1.e-7
top.mu0 = 4*top.pi/10000000
# --- Conversion factor from joules to eV is just echarge
top.jperev = top.echarge
# --- Epsilon_0 calculated from speed of light and mu_0
top.eps0 = 1/(top.mu0*top.clight*top.clight)
# --- Create python versions of the constants
amu     = top.amu
clight  = top.clight
echarge = top.echarge
emass   = top.emass
eps0    = top.eps0
euler   = top.euler
jperev  = top.jperev
mu0     = top.mu0

# --- Import the convenience routines for plotting various slices and
# --- projections of particles as well as some line plots.
from warpplots import *

# --- Import some online documentation modules.
from warphelp import *
from warpscripts import *
from warpfortran import *

# --- Declare the documentation for the warp module.
def warpdoc():
  print """
Imports the basic modules needed to run WARP.

Create python versions of the constants
amu, clight, echarge, emass, eps0, euler, jperev, mu0

Defines following functions...
dump: Creates a dump file containing the current state of the simulation
restart: Retreives the state of a simulation from a dump file
loadrho: Load the particles onto the charge density mesh
fieldsol: Solve for the self-fields
installbeforefs: Install a function which will be called before a field-solve
uninstallbeforefs: Uninstall the function called before a field-solve
installafterfs: Install a function which will be called after a field-solve
uninstallafterfs: Uninstall the function called after a field-solve
installbeforestep: Install a function which will be called before a step
uninstallbeforestep: Uninstall the function called before a step
installafterstep: Install a function which will be called after a step
uninstallafterstep: Uninstall the function called after a step
gethzarrays: Fixes the ordering of hlinechg and hvzofz data from a paralle run
  """

print 'Python WARP'
print '******  Fieldsolver FXY version %s' % arraytostr(fxy.versfxy)
#print '******  Fieldsolver FRZ version %s' % arraytostr(frz.versfrz)
print '******  Fieldsolver F3D version %s' % arraytostr(f3d.versf3d)
print '******  Envelope solver ENV version %s' % arraytostr(env.versenv)
try:
  print '******  Envelope solver CIR version %s' % arraytostr(cir.verscir)
except:
  pass
print '******  Particle package WXY version %s' % arraytostr(wxy.verswxy)
#print '******  Particle package WRZ version %s' % arraytostr(wrz.verswrz)
print '******  Particle package W3D version %s' % arraytostr(w3d.versw3d)
print '******  Main package TOP version %s' % arraytostr(top.verstop)
print 'For more help, type warphelp()'

# --- Call derivqty to calculate eps0 and jperev
derivqty()

# --- Setup mechanism for before and after field-solve python scripts
beforefsfuncs = []
afterfsfuncs = []
def beforefs():
  for f in beforefsfuncs: f()
def afterfs():
  for f in afterfsfuncs: f()
def installbeforefs(f):
  "Adds a function to the list of functions called before a field-solve"
  beforefsfuncs.append(f)
  w3d.lbeforefs = true
def uninstallbeforefs(f):
  "Removes the function from the list of functions called before a field-solve"
  try:
    beforefsfuncs.remove(f)
  except ValueError:
    print 'Warning: uninstallbeforefs: function was not installed'
  if len(beforefsfuncs) == 0: w3d.lbeforefs = false
def installafterfs(f):
  "Adds a function to the list of functions called after a field-solve"
  afterfsfuncs.append(f)
  w3d.lafterfs = true
def uninstallafterfs(f):
  "Removes the function from the list of functions called after a field-solve"
  try:
    afterfsfuncs.remove(f)
  except ValueError:
    print 'Warning: uninstallafterfs: function was not installed'
  if len(afterfsfuncs) == 0: w3d.lafterfs = false
def installbeforestep(f):
  "Adds a function to the list of functions called before a step"
  beforestepfuncs.append(f)
def uninstallbeforestep(f):
  "Removes the function from the list of functions called before a step"
  try:
    beforestepfuncs.remove(f)
  except ValueError:
    print 'Warning: uninstallbeforestep: function was not installed'
def installafterstep(f):
  "Adds a function to the list of functions called after a step"
  afterstepfuncs.append(f)
def uninstallafterstep(f):
  "Removes the function from the list of functions called after a step"
  try:
    afterstepfuncs.remove(f)
  except ValueError:
    print 'Warning: uninstallafterstep: function was not installed'

# --- Convenience function for random numbers.
def setseed(x=0,y=0):
  RandomArray.seed(long(x),long(y))
    
# --- Uniform distribution
# --- Uniform distribution
def ranf(x=None,i1=None,nbase=None):
  if not i1:
    if x == None:
      try:
        return ranlib.sample()
      except:
        return RandomArray.random(shape(array([0])))[0]
    else:
      return RandomArray.random(shape(x))
  else:
    if not nbase: nbase = 2
    n = product(array(shape(array(x))))
    result = zeros(n,'d')
    rnrevarray(n,result,i1,nbase)
    result.shape = shape(x)
    return result

# --- Gaussian distribution
# --- This was in pyBasis but had to be moved here in order to use rnormdig.
def rnormarray(x,i1=None,nbase1=None,nbase2=None):
  if not i1:
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
def loadrho(ins_i=-1,nps_i=-1,is_i=-1,lzero=1):
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
def fieldsol(iwhich=0):
  """
  # --- This routine provides a simple call from the interpreter to do
  # --- the fieldsol.
  # --- Call the appropriate compiled interface routine based on the
  # --- current package
  """
  currpkg = package()[0]
  if (currpkg == "w3d"):
    fieldsol3d(iwhich)
  elif (currpkg == "wxy"):
    fieldsolxy(iwhich)
#=============================================================================

# --- Dump command
def dump(filename=None,suffix='',attr='dump',serial=0,onefile=1,pyvars=1,
         ff=None):
  """
Creates a dump file
  - filename=(runid+'%06d'%top.it+suffix+'.dump')
  - attr='dump': All variables with the given attribute or group name are
    written to the file. The default attribute makes a restartable dump file.
  - serial=0: When 1, does a dump of only non-parallel data (parallel version
    only).
  - onefile=1: When 1, all processors dump to one file, otherwise each dumps to
    seperate file.
  - pyvars=1: When 1, saves user defined python variables to the file.
  - ff=None: Optional file object. When passed in, write to that file instead
             of creating a new one. Note that the inputted file object must be
             closed by the user.
  """
  if not filename:
    # --- Setup default filename based on time step and processor number.
    if serial:
      s = '.sdump'
    else:
      s = '.dump'
    if onefile or npes==0:
      filename=arraytostr(top.runid)+('%06d'%top.it)+suffix+s
    else:
      filename=arraytostr(top.runid)+('%06d_%05d'%(top.it,me))+suffix+s
  # --- Make list of all of the new python variables.
  interpreter_variables = []
  if pyvars:
    for l in __main__.__dict__.keys():
      if l not in initial_global_dict_keys:
        interpreter_variables.append(l)
  # --- Call routine to make data dump
  if onefile and npes > 0:
    paralleldump(filename,attr,interpreter_variables,serial=serial)
  else:
    pydump(filename,attr,interpreter_variables,serial=serial,ff=ff)

# --- Restart command
def restart(filename,onefile=1):
  """
Reads in data from file, redeposits charge density and does field solve
  - filename: restart file name - when restoring parallel run from multiple
    files, filename should only be prefix up to but not including the '_'
    before the processor number.
  - onefile=1: Restores from one file unless 0, then each processor restores
    from seperate file.
  """
  # --- If each processor is restoring from a seperate file, append
  # --- appropriate suffix, assuming only prefix was passed in
  if not onefile:
    filename = filename + '_%05d.dump'%(me)
  # --- Call different restore routine depending on context
  if onefile and npes > 0:
    parallelrestore(filename)
  else:
    pyrestore(filename)
  # --- Finish up the restart work
  loadrho()
  fieldsol(0)
  setlatt()

def restartold(filename):
  restoreold(filename)
  loadrho()
  fieldsol(0)
  setlatt()


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
