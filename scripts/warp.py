# import all of the neccesary packages
import __main__
from Numeric import *
import ranlib
import sys
warp_version = "$Id: warp.py,v 1.17 2001/02/10 01:38:21 dave Exp $"

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

try:
  if me == 0 and sys.platform != 'mac':
    from gist import *
  else:
    from gistdummy import *
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
#from frzpy import *
#from wrzpy import *
try:  # cirpy hasn't been installed on all machines yet
  from cirpy import *
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
sys.path = sys.path + ['/home/ife1/dave/warp/scripts']

#=============================================================================
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

#=============================================================================
# --- Setup and make initial printout of the versions of the packages.
def printversion(v):
  v = arraytostr(v)
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
  if 'wxy' in pkg: r=r+fmt%('Particle package WXY',printversion(wxy.verswxy))
  if 'wrz' in pkg: r=r+fmt%('Particle package WRZ',printversion(wrz.verswrz))
  if 'w3d' in pkg: r=r+fmt%('Particle package W3D',printversion(w3d.versw3d))
  if 'top' in pkg: r=r+fmt%('Main package TOP',printversion(top.verstop))
  return r
print versionstext()[:-1] # --- skip last line feed
print 'For more help, type warphelp()'

#=============================================================================
# --- Import the convenience routines for plotting various slices and
# --- projections of particles, histories, as well as some line plots.
from warpplots import *
from histplots import *
from pzplots import *

# --- Import some online documentation modules.
from warphelp import *
from warpscripts import *
from warpfortran import *

# --- Import the printparameters modules (which are called from fortran)
from printparameters import *
from printparameters3d import *
from printparametersrz import *

#=============================================================================
# --- Declare the documentation for the warp module.
def warpdoc():
  print """
Imports the basic modules needed to run WARP, including
Numeric, gist, warpplots, histplots, pzplots

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
printtimers: Print timers in a nice annotated format
  """

#=============================================================================
# --- Call derivqty to calculate eps0 and jperev
derivqty()

#=============================================================================
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
  if f in beforefsfuncs:
    beforefsfuncs.remove(f)
  else:
    raise 'Warning: uninstallbeforefs: no such function had been installed'
  if len(beforefsfuncs) == 0: w3d.lbeforefs = false
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
def installbeforestep(f):
  "Adds a function to the list of functions called before a step"
  beforestepfuncs.append(f)
def uninstallbeforestep(f):
  "Removes the function from the list of functions called before a step"
  if f in beforestepfuncs:
    beforestepfuncs.remove(f)
  else:
    raise 'Warning: uninstallbeforestep: no such function had been installed'
def installafterstep(f):
  "Adds a function to the list of functions called after a step"
  afterstepfuncs.append(f)
def uninstallafterstep(f):
  "Removes the function from the list of functions called after a step"
  if f in afterstepfuncs:
    afterstepfuncs.remove(f)
  else:
    raise 'Warning: uninstallafterstep: no such function had been installed'

#=============================================================================
# --- Convenience function for random numbers.
def setseed(x=0,y=0):
  RandomArray.seed(long(x),long(y))
    
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
         ff=None,varsuffix=None,histz=0):
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
  - varsuffix=None Suffix to add to the variable names. If none is specified,
                   the suffix '@pkg' is used, where pkg is the package name
                   that the variable is in.
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
  # --- Make list of all of the new python variables.
  interpreter_variables = []
  if pyvars:
    for l in __main__.__dict__.keys():
      if l not in initial_global_dict_keys:
        interpreter_variables.append(l)
  # --- Call routine to make data dump
  if onefile and lparallel:
    paralleldump(filename,attr,interpreter_variables,serial=serial,
                 varsuffix=varsuffix,histz=histz)
  else:
    pydump(filename,attr,interpreter_variables,serial=serial,ff=ff,
           varsuffix=varsuffix)
  top.dumptime = top.dumptime + (wtime() - timetemp)

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
  if onefile and lparallel:
    parallelrestore(filename)
  else:
    pyrestore(filename)
  # --- Now that the dump file has been read in, finish up the restart work.
  # --- First set the current packge. Note that currpkg is only ever defined
  # --- in the main dictionary.
  package(__main__.__dict__["currpkg"])
  # --- Allocate all arrays appropriately
  gchange("*")
  # --- Load the charge density (since it was not saved)
  loadrho()
  # --- Recalculate the fields (since it was not saved)
  fieldsol(0)
  # --- Set the lattice internal variables (probably not really needed)
  setlatt()
  # --- Reinitialize some injection stuff if it is needed.
  # --- This is really only neede for the parallel version since some of the
  # --- data saved is only valid for PE0.
  if top.inject > 0: injctint()
  # --- Call setup if it is needed.
  if me == 0 and current_window() == -1: setup()

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
def printtimers(file=None):
  """Print timers in a nice annotated format
  - file Optional input file. If it is not include, stdout is used. It can
         either be a file name, or a file object. If a file name, a file
         with that name is created. If a file object, the data is written
         to that file (the file remains open).
  """
  if file == None:
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
    for t in timelists: timedevs.append(sqrt(ave(t**2) - ave(t)**2))
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
