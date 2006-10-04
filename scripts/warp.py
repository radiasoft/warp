warp_version = "$Id: warp.py,v 1.126 2006/10/04 01:39:38 dave Exp $"
# import all of the neccesary packages
import __main__
from Numeric import *
import MA
import sys
import os.path
import time

# --- Set this to a more reasonable value
MA.set_print_limit(10000)

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

# --- Since gist is only loaded on PE0 in the parallel
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

# Import the warpC shared object which contains all of WARP
from warpC import *

from Forthon import *
from warputils import *
import RandomArray

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
from em2dpy import *
import controllers
from controllers import *
from ctl import *

# --- Rearrange the list of packages to a more sensible order
package('wxy')
package('w3d')
package('top')

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

# --- Check if the compiler was ifort - if so, set the stacksize unlimited
# --- The fcompname is not yet be available yet if Forthon is not up to date
try:
  if fcompname == 'ifort':
    import resource
    resource.setrlimit(resource.RLIMIT_STACK,(-1,-1))
except:
  pass

#try:
#  from psyco.classes import psyobj
#  import psyco
#  psyco.background(watermark=.3)
#except ImportError:
#  class psyobj:
#    pass

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

# --- Get start time
top.starttime = time.time()
top.starttimedump = top.starttime

# --- Simple function to calculate Child-Langmuir current density
def childlangmuir(v,d,q=None,m=None):
  """
Computes current density from Child-Langmuir formula
 - v: diode voltage
 - d: diode length
 - q: particle charge (defaults to top.sq[0])
 - m: particle mass (defaults to top.sm[0])
  """
  if q is None: q = top.sq[0]
  if m is None: m = top.sm[0]
  return 4./9.*eps0*sqrt(2.*q/m)*v**1.5/d**2

# --- Create 


# --- Create python version of dvnz (divisor not zero)
def dvnz(x):
  if type(x) in [IntType,FloatType]:
    if x == 0.: return top.smallpos
    else:       return x
  else:
    return where(x==0.,top.smallpos,x)
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
print time.ctime(top.starttime)
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
getappliedfields: Gathers the fields from the accelerator lattice at given
                  locations
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
isinstalledbeforestep: Checks if a function is installed to be called before a
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
# --- This had to be moved here in order to use rnormdig.
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
def addspecies(newns=1,pgroup=None,sid=None):
  """Adds one or more new speices. This only allocates the needed arrays.
  - newns=1: the number of new species to add
  - pgroup=top.pgroup: the pgroup to add them to
  - sid=iota(top.ns-newns,top.ns-1): the global id of the species being added
                                     This is only meaningful if pgroup is
                                     specified.
  """
  if newns == 0: return
  top.ns = top.ns + newns
  assert top.ns >= 0,"The total number of species cannot be negative"
  gchange("InPart")
  gchange("InjectVars")
  gchange("LostParticles")
  gchange("ExtPart")
  # --- top.pgroup is special since it always has top.ns species.
  top.pgroup.ns = top.ns
  if sid is None: sid = iota(top.ns-newns,top.ns-1)
  setuppgroup(top.pgroup)
  top.pgroup.sid[-newns:] = sid
  if pgroup is not None:
    pgroup.ns = top.ns
    setuppgroup(pgroup)
    top.pgroup.sid[-newns:] = sid
  if top.lspeciesmoments:
    top.nszarr = top.ns
    gchange("Z_arrays")
    top.nswind = top.ns
    gchange("Win_Moments")
    top.nszmmnt = top.ns
    gchange("Z_Moments")
    top.nslabwn = top.ns
    gchange("Lab_Moments")
    top.nsmmnt = top.ns
    gchange("Moments")
    top.nshist = top.ns
    gchange("Hist")
#=============================================================================
def setnspecies(ns):
  """Set the total number of species."""
  assert ns >= 0,"The total number of species cannot be negative"
  addspecies(ns-top.ns)

##############################################################################
def getappliedfieldsongrid(nx=None,ny=None,nz=None,
                           xmin=None,xmax=None,
                           ymin=None,ymax=None,
                           zmin=None,zmax=None):
  """
Gets the applied fields from the lattice on a grid. 
It returns the tuple (ex,ey,ez,bx,by,bz).
Inputs (optional) are the number of cells in each dimension and the limits of 
the grid: nx, ny, nz, xmin, xmax, ymin, ymax, zmin, zmax. The defaults values 
are respectively: w3d.nx, w3d.ny, w3d.nz, w3d.xmmin, w3d.xmmax,
w3d.ymmin, w3d.ymmax, w3d.zmmin, w3d.zmmax.
  """

  if nx is None: nx=w3d.nx
  if ny is None: ny=w3d.ny
  if nz is None: nz=w3d.nz
  if xmin is None: xmin=w3d.xmmin
  if xmax is None: xmax=w3d.xmmax
  if ymin is None: ymin=w3d.ymmin
  if ymax is None: ymax=w3d.ymmax
  if zmin is None: zmin=w3d.zmmin
  if zmax is None: zmax=w3d.zmmax

  # --- Calculate the grid cell sizes, checking for zero values
  if nx > 0:
    dx = (xmax - xmin)/nx
  else:
    dx = 0.
    xmin = 0.
    xmax = 0.
  if ny > 0:
    dy = (ymax - ymin)/ny
  else:
    dy = 0.
    ymin = 0.
    ymax = 0.
  if nz > 0:
    dz = (zmax - zmin)/nz
  else:
    dz = 0.
    zmin = 0.
    zmax = 0.

  x,y,z = getmesh3d(xmin,dx,nx,ymin,dy,ny,zmin,dz,nz)
  nxyz = product(array(shape(x)))
  x.shape = (nxyz,)
  y.shape = (nxyz,)
  z.shape = (nxyz,)
  ex,ey,ez,bx,by,bz = getappliedfields(x,y,z)
  ex.shape = (1+nx,1+ny,1+nz)
  ey.shape = (1+nx,1+ny,1+nz)
  ez.shape = (1+nx,1+ny,1+nz)
  bx.shape = (1+nx,1+ny,1+nz)
  by.shape = (1+nx,1+ny,1+nz)
  bz.shape = (1+nx,1+ny,1+nz)
  return ex,ey,ez,bx,by,bz

##############################################################################
def getappliedfields(x,y,z,time=0.,js=0):
  """
Gets the applied fields from the lattice at the given locations.
It returns the tuple (ex,ey,ez,bx,by,bz)
 - x,y,z: arrays of position where the field is to be gathered.
          Note that z can be a scalar - all fields are gathered at that
          z position.
 - time=0.: Time to use whan gathering fields - only affects time dependent
            elements.
 - js=0: species to get mass and charge from. Only affects accl elements.
  """
  n = len(x)
  if n == 0: return

  # --- Allow z to be a scalar (as in the slice case)
  if type(z) in [IntType,FloatType]: z = z*ones(n,'d')

  # --- Save existing internal lattice variables so they can be restored
  # --- afterward.
  dzlsave = top.dzl
  dzlisave = top.dzli
  zlframesave = top.zlframe
  zltimesave = top.zltime
  zlminsave = top.zlmin
  zlmaxsave = top.zlmax
  nzlsave = top.nzl
  nzlmaxsave = top.nzlmax

  ex,ey,ez = zeros((3,n),'d')
  bx,by,bz = zeros((3,n),'d')

  zmin = min(z)
  zmax = max(z)

  # --- Set them based on the particle data
  zlen = zmax - zmin
  if zlen == 0.:
    top.zlmin = zmin
    top.zlmax = zmin
    top.nzl = 0
  else:
    top.zlmin = zmin - zlen/10.
    top.zlmax = zmax + zlen/10.
    top.nzl = 100
    top.dzl = (top.zlmax - top.zlmin)/top.nzl
    top.dzli = 1./top.dzl
    top.zlframe = 0.
    top.nzlmax = max(top.nzlmax,top.nzl)
    gchange("LatticeInternal")
  top.zltime = time
    
  # --- Make sure that the lattice is set up. If it is already, this won't
  # --- change anything (maybe).
  setlatt()

  # --- Create other temporaries
  uzp = ones(n,'d')
  gaminv = ones(n,'d')
  dtl = -0.5*top.dt
  dtr = +0.5*top.dt
  m = top.pgroup.sm[js]
  q = top.pgroup.sq[js]
  bendres = ones(n,'d')
  bendradi = ones(n,'d')
  gammabar = 1.
  dt = top.dt

  exteb3d(n,x,y,z,uzp,gaminv,dtl,dtr,
          bx,by,bz,ex,ey,ez,m,q,bendres,bendradi,gammabar,dt)

  # --- Restore these quantities
  top.dzl = dzlsave
  top.dzli = dzlisave
  top.zlframe = zlframesave
  top.zltime = zltimesave
  top.zlmin = zlminsave
  top.zlmax = zlmaxsave
  top.nzl = nzlsave
  top.nzlmax = nzlmaxsave
  gchange("LatticeInternal")
  setlatt()

  return ex,ey,ez,bx,by,bz


##############################################################################
##############################################################################
##############################################################################
def fixrestoresfrombeforeelementoverlaps(filename):
  """This checks if the lattice overlap data is inconsistent with the
lattice input. This should only ever happen if an old dump file is
read in, one created before the element overlaps where implemented.
If this is the case, the lattice is reset and the overlap data generated."""
  doreset = 0
  if top.ndrft >= 0 and sum(top.odrftnn) < top.ndrft+1: doreset = 1
  if top.nbend >= 0 and sum(top.obendnn) < top.nbend+1: doreset = 1
  if top.ndipo >= 0 and sum(top.odiponn) < top.ndipo+1: doreset = 1
  if top.nquad >= 0 and sum(top.oquadnn) < top.nquad+1: doreset = 1
  if top.nsext >= 0 and sum(top.osextnn) < top.nsext+1: doreset = 1
  if top.nhele >= 0 and sum(top.ohelenn) < top.nhele+1: doreset = 1
  if top.nemlt >= 0 and sum(top.oemltnn) < top.nemlt+1: doreset = 1
  if top.nmmlt >= 0 and sum(top.ommltnn) < top.nmmlt+1: doreset = 1
  if top.naccl >= 0 and sum(top.oacclnn) < top.naccl+1: doreset = 1
  if top.nbgrd >= 0 and sum(top.obgrdnn) < top.nbgrd+1: doreset = 1
  if top.npgrd >= 0 and sum(top.opgrdnn) < top.npgrd+1: doreset = 1
  if top.nbsqgrad >= 0 and sum(top.obsqgradnn) < top.nbsqgrad+1: doreset = 1
  if doreset:
    resetlat()
    setlatt()
  
def fixrestoreswithmomentswithoutspecies(filename):
  """This is called automatically by restart to fix the arrays which have
  changed shape. It needs to be called by hand after restore.
  """
  # --- Open the file
  ff = PR.PR(filename,verbose=0)
  # --- First, check if the file has the old moments in it.
  ek = ff.read('ek@top')
  if type(ek) == ArrayType:
    # --- If it has the new ones, do nothing.
    ff.close()
    return
  # --- These variables went from scalars to 1-D arrays
  varlist1 = [ 'ek', 'ekzmbe', 'ekzbeam', 'ekperp', 'pz', 'xmaxp', 'xminp',
               'ymaxp', 'yminp', 'zmaxp', 'zminp', 'vxmaxp', 'vxminp',
               'vymaxp', 'vyminp', 'vzmaxp', 'vzminp']
  # --- These arrays changed shape. There's probably a better way of doing
  # --- this without just explictly listing everything.
  varlist = ['curr','lostpars', 'pnum', 'xbar', 'ybar', 'zbar', 'xpbar',
             'ypbar', 'vxbar', 'vybar', 'vzbar', 'xybar', 'xypbar', 'yxpbar',
             'xpypbar', 'xsqbar', 'ysqbar', 'zsqbar', 'xpsqbar', 'ypsqbar',
             'vxsqbar', 'vysqbar', 'vzsqbar', 'xxpbar', 'yypbar', 'zvzbar',
             'xvzbar', 'yvzbar', 'vxvzbar', 'vyvzbar', 'xrms', 'yrms', 'zrms',
             'rrms', 'xprms', 'yprms', 'epsx', 'epsy', 'epsz', 'epsnx',
             'epsny', 'epsnz', 'epsg', 'epsh', 'epsng', 'epsnh', 'vxrms',
             'vyrms', 'vzrms', 'pnumz', 'xbarz', 'ybarz', 'zbarz', 'xpbarz',
             'ypbarz', 'vxbarz', 'vybarz', 'vzbarz', 'xybarz', 'xypbarz',
             'yxpbarz', 'xpypbarz', 'xsqbarz', 'ysqbarz', 'zsqbarz',
             'xpsqbarz', 'ypsqbarz', 'vxsqbarz', 'vysqbarz', 'vzsqbarz',
             'xxpbarz', 'yypbarz', 'zvzbarz', 'xvzbarz', 'yvzbarz',
             'vxvzbarz', 'vyvzbarz', 'xrmsz', 'yrmsz', 'zrmsz', 'rrmsz',
             'xprmsz', 'yprmsz', 'epsxz', 'epsyz', 'epszz', 'epsnxz',
             'epsnyz', 'epsnzz', 'epsgz', 'epshz', 'epsngz', 'epsnhz',
             'vxrmsz', 'vyrmsz', 'vzrmsz', 'pnumlw', 'xbarlw', 'ybarlw',
             'vzbarlw', 'epsxlw', 'epsylw', 'epszlw', 'vxrmslw', 'vyrmslw',
             'vzrmslw', 'xrmslw', 'yrmslw', 'rrmslw', 'xxpbarlw', 'yypbarlw',
             'currlw', 'lostparslw', 'hbmlen', 'hzbar', 'hekzmbe', 'hekzbeam',
             'hekperp', 'hxmaxp', 'hxminp', 'hymaxp', 'hyminp', 'hzmaxp',
             'hzminp', 'hvxmaxp', 'hvxminp', 'hvymaxp', 'hvyminp', 'hvzmaxp',
             'hvzminp', 'hepsx', 'hepsy', 'hepsz', 'hepsnx', 'hepsny',
             'hepsnz', 'hepsg', 'hepsh', 'hepsng', 'hepsnh', 'hpnum', 'hxbar',
             'hybar', 'hxybar', 'hxrms', 'hyrms', 'hrrms', 'hxprms', 'hyprms',
             'hxsqbar', 'hysqbar', 'hvxbar', 'hvybar', 'hvzbar', 'hxpbar',
             'hypbar', 'hvxrms', 'hvyrms', 'hvzrms', 'hxpsqbar', 'hypsqbar',
             'hxxpbar', 'hyypbar', 'hxypbar', 'hyxpbar', 'hxpypbar',
             'hxvzbar', 'hyvzbar', 'hvxvzbar', 'hvyvzbar', 'hcurrz', 'hepsxz',
             'hepsyz', 'hepsnxz', 'hepsnyz', 'hepsgz', 'hepshz', 'hepsngz',
             'hepsnhz', 'hxbarz', 'hybarz', 'hxybarz', 'hxrmsz', 'hyrmsz',
             'hrrmsz', 'hxprmsz', 'hyprmsz', 'hxsqbarz', 'hysqbarz',
             'hvxbarz', 'hvybarz', 'hvzbarz', 'hxpbarz', 'hypbarz', 'hvxrmsz',
             'hvyrmsz', 'hvzrmsz', 'hxpsqbarz', 'hypsqbarz', 'hxxpbarz',
             'hyypbarz', 'hxypbarz', 'hyxpbarz', 'hxpypbarz', 'hxvzbarz',
             'hyvzbarz', 'hvxvzbarz', 'hvyvzbarz']
  # --- Setup the arrays
  if top.lspeciesmoments and top.ns > 1:
    top.nszmmnt = top.ns
  else:
    top.nszmmnt = 0
  top.nszarr = top.nszmmnt
  top.nswind = top.nszmmnt
  top.nslabwn = top.nszmmnt
  top.nsmmnt = top.nszmmnt
  top.nshist = top.nszmmnt
  gchange('Z_arrays')
  gchange('Win_Moments')
  gchange('Z_Moments')
  gchange('Lab_Moments')
  gchange('Moments')
  gchange('Hist')
  # --- Get list of variables in the file.
  fflist = ff.inquire_names()
  # --- For each one in the file, put the data in the last element.
  for v in varlist1:
    if v+'@top' in fflist:
      d = ff.read(v+'@top')
      a = getattr(top,v)
      a[-1] = d
  for v in varlist:
    if v+'@top' in fflist:
      d = ff.read(v+'@top')
      a = getattr(top,v)
      a[...,-1] = d

def fixrestoreswithoriginalparticlearrays(filename):
  # --- Open the file
  ff = PR.PR(filename,verbose=0)
  # --- Check if it is an old file
  # --- An old file would have top.npmaxb save in it
  if 'xp@pgroup@top' in ff.inquire_names():
    ff.close()
    return
  # --- Setup top.pgroup
  # --- Do this first since ipmax_s needs the correct size
  top.pgroup.ns = top.ns
  # --- Only these needs to be read in.
  top.pgroup.ipmax_s = ff.read('npmax_s@top')
  top.pgroup.npmax = ff.read("npmax@top")
  top.pgroup.npid = ff.read("npid@top")
  top.pgroup.sm = ff.read("sm@top")
  top.pgroup.sq = ff.read("sq@top")
  top.pgroup.sw = ff.read("sw@top")
  top.pgroup.ins = ff.read("ins@top")
  top.pgroup.nps = ff.read("nps@top")
  top.pgroup.ndts = ff.read("ndts@top")
  top.pgroup.ldts = ff.read("ldts@top")
  top.pgroup.dtscale = ff.read("dtscale@top")
  top.pgroup.lselfb = ff.read("lselfb@top")
  top.pgroup.fselfb = ff.read("fselfb@top")

  try:
    top.pgroup.sid = arange(top.ns)
  except:
    top.pgroup.js = arange(top.ns)

  top.pgroup.gaminv = ff.read("gaminv@top")
  top.pgroup.xp = ff.read("xp@top")
  top.pgroup.yp = ff.read("yp@top")
  top.pgroup.zp = ff.read("zp@top")
  top.pgroup.uxp = ff.read("uxp@top")
  top.pgroup.uyp = ff.read("uyp@top")
  top.pgroup.uzp = ff.read("uzp@top")
  top.pgroup.pid = ff.read("pid@top")

  ff.close()

def fixrestoreswitholdparticlearrays(filename):
  # --- Open the file
  ff = PR.PR(filename,verbose=0)
  # --- Check if it is an old file
  if 'xp@top' not in ff.inquire_names():
    ff.close()
    return

  top.pgroup.ns = top.ns
  top.pgroup.npid = ff.read("npid@top")
  top.pgroup.sm = ff.read("sm@top")
  top.pgroup.sq = ff.read("sq@top")
  top.pgroup.sw = ff.read("sw@top")
  top.pgroup.nps = 0
  top.pgroup.ndts = ff.read("ndts@top")
  top.pgroup.ldts = ff.read("ldts@top")
  top.pgroup.dtscale = ff.read("dtscale@top")
  top.pgroup.lselfb = ff.read("lselfb@top")
  top.pgroup.fselfb = ff.read("fselfb@top")

  gaminv = ff.read("gaminv@top")
  xp = ff.read("xp@top")
  yp = ff.read("yp@top")
  zp = ff.read("zp@top")
  uxp = ff.read("uxp@top")
  uyp = ff.read("uyp@top")
  uzp = ff.read("uzp@top")
  if top.pgroup.npid > 0:
    pid = ff.read("pid@top")
  else:
    pid = None

  ins = ff.read("ins@top")
  nps = ff.read("nps@top")
  for js in range(top.pgroup.ns):
    i1 = ins[0] - 1
    i2 = ins[0] + nps[0] - 1
    if pid is None:
      addparticles(x=xp[i1:i2],y=yp[i1:i2],z=zp[i1:i2],
                   vx=uxp[i1:i2],vy=uyp[i1:i2],vz=uzp[i1:i2],
                   gi=gaminv[i1:i2],pid=0.,js=js,
                   lallindomain=false,
                   zmmin=w3d.zmmin,zmmax=w3d.zmmax,lmomentum=true,
                   resetrho=false,dofieldsol=false,resetmoments=false)
    else:
      addparticles(x=xp[i1:i2],y=yp[i1:i2],z=zp[i1:i2],
                   vx=uxp[i1:i2],vy=uyp[i1:i2],vz=uzp[i1:i2],
                   gi=gaminv[i1:i2],pid=pid[:,i1:i2],js=js,
                   lallindomain=false,
                   zmmin=w3d.zmmin,zmmax=w3d.zmmax,lmomentum=true,
                   resetrho=false,dofieldsol=false,resetmoments=false)

def restoreolddump(filename):
  pass
  #fixrestoresfrombeforeelementoverlaps(filename)
  #fixrestoreswithmomentswithoutspecies(filename)
  #fixrestoreswithoriginalparticlearrays(filename)
  #fixrestoreswitholdparticlearrays(filename)

##############################################################################
##############################################################################
##############################################################################
##############################################################################
##############################################################################
##############################################################################
##############################################################################
##############################################################################
# --- Dump command
def dump(filename=None,prefix='',suffix='',attr='dump',serial=0,onefile=0,pyvars=1,
         ff=None,varsuffix=None,histz=2,resizeHist=1,verbose=false,hdf=0):
  """
Creates a dump file
  - filename=(prefix+runid+'%06d'%top.it+suffix+'.dump')
  - attr='dump': All variables with the given attribute or group name are
    written to the file. The default attribute makes a restartable dump file.
  - serial=0: When 1, does a dump of only non-parallel data (parallel version
    only).
  - onefile=0: When 0, all processors dump to one file, otherwise each dumps to
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
      filename=prefix+arraytostr(top.runid)+('%06d'%top.it)+suffix+s
    else:
      filename=prefix+arraytostr(top.runid)+('%06d_%05d_%05d'%(top.it,me,npes))+suffix+s
  else:
    if not onefile and lparallel:
      # --- Append the processor number to the user inputted filename
      filename = filename + '%05d_%05d'%(me,npes)
  print filename
  # --- Make list of all of the new python variables.
  interpreter_variables = []
  if pyvars:
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
  # --- Update dump time
  top.dumptime = top.dumptime + (wtime() - timetemp)

# --- Restart command
def restart(filename,onefile=0,verbose=false,dofieldsol=true):
  """
Reads in data from file, redeposits charge density and does field solve
  - filename: restart file name - when restoring parallel run from multiple
              files, filename should only be prefix up to but not including
              the '_' before the processor number.
  - onefile=0: Restores from one file unless 0, then each processor restores
               from seperate file.
  - dofieldsol=true: When true, call fieldsol(0). This allows special cases
                     where just calling fieldsol(0) is not appropriate or
                     optimal
  """
  # --- If each processor is restoring from a seperate file, append
  # --- appropriate suffix, assuming only prefix was passed in
  if not onefile:
    filename = filename + '_%05d_%05d.dump'%(me,npes)
  # --- Call different restore routine depending on context
  if onefile and lparallel:
    parallelrestore(filename,verbose=verbose)
  else:
    pyrestore(filename,verbose=verbose)
  # --- Fix old dump files.
  restoreolddump(filename)
  # --- Now that the dump file has been read in, finish up the restart work.
  # --- First set the current packge. Note that currpkg is only ever defined
  # --- in the main dictionary.
  package(__main__.__dict__["currpkg"])
  # --- Allocate all arrays appropriately
  gchange("*")
  # --- Reinitialize some injection stuff if it is needed.
  # --- This is really only needed for the parallel version since some of the
  # --- data saved is only valid for PE0.
  if top.inject > 0: fill_inj()
  # --- Do some setup for the RZ solver
  if getcurrpkg() == 'w3d' and w3d.solvergeom==w3d.RZgeom: mk_grids_ptr()
  # --- Load the charge density (since it was not saved)
  if not (w3d.solvergeom in [w3d.RZgeom] or top.fstype == 12):
    loadrho()
    # --- Recalculate the fields (since it was not saved)
    if dofieldsol: fieldsol(0)
  # --- Call setup if it is needed.
  if me == 0 and current_window() == -1: setup()
  # --- Call any functions that had been registered to be called after
  # --- the restart.
  controllers.callafterrestartfuncs()

##############################################################################
##############################################################################
##############################################################################
##############################################################################
##############################################################################
##############################################################################
##############################################################################
##############################################################################

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
from fieldsolver import *
if lparallel: from warpparallel import *
from warpplots import *
from histplots import *
from pzplots import *
from plot_conductor import *
from lattice import *
from drawlattice import *
from MeshRefinement import *
from species import *

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
  warp_path = os.path.dirname(warp.__file__)
  if warp_path <> '':warp_path+='/'
  sys.path=sys.path+[warp_path+'GUI']
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

# --- The controller function container needs to be written out since the
# --- controllers functions may be changed by the user. The container
# --- properly reinstalls any saved controller functions.
# --- The name 'controllerfunctioncontainer' must be the same as what appears
# --- in the controllers module.
initial_global_dict_keys.remove('controllerfunctioncontainer')

# --- Save the versions string here so that it will be dumped into any
# --- dump file.
warpversions = versionstext()
