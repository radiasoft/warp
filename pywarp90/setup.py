#!/usr/bin/env python
# To use:
#       python setup.py install
#
import sys,os,os.path,string
from Forthon.compilers import FCompiler
import getopt

try:
    import distutils
    from distutils.core import setup, Extension
    from distutils.dist import Distribution
    from distutils.command.build import build
except:
    raise SystemExit, "Distutils problem"

optlist,args = getopt.getopt(sys.argv[1:],'gt:F:',['parallel'])
machine = sys.platform
debug   = 0
fcomp   = None
parallel = 0
for o in optlist:
  if   o[0] == '-g': debug = 1
  elif o[0] == '-t': machine = o[1]
  elif o[0] == '-F': fcomp = o[1]
  elif o[0] == '--parallel': parallel = 1

sys.argv = ['setup.py']+args
fcompiler = FCompiler(machine=machine,
                      debug=debug,
                      fcompname=fcomp)

dummydist = Distribution()
dummybuild = build(dummydist)
dummybuild.finalize_options()
builddir = dummybuild.build_temp

warppkgs = ['top','env','w3d','f3d','wxy','fxy','wrz','frz','her','cir','cho']

# --- The behavior of distutils changed from 2.2 to 2.3. In 2.3, the object
# --- files are always put in a build/temp directory relative to where the
# --- source file is, rather than relative to the main build directory.
if sys.hexversion >= 0x020300f0:
  pymodprefix = builddir
else:
  pymodprefix = ''

def makeobjects(pkg):
  return [pkg+'.o',pkg+'_p.o',os.path.join(pymodprefix,pkg+'pymodule.o')]

warpobjects = []
for pkg in warppkgs:
  warpobjects = warpobjects + makeobjects(pkg)

warpobjects = warpobjects + ['dtop.o',
                             'dw3d.o','w3d_injection.o','w3d_interp.o',
                             'w3d_utilities.o','w3d_load.o', 
                             'f3d_mgrid.o','f3d_conductors.o','fft.o','util.o',
                             'fxy_mgrid.o',
                             'dwrz.o',
                             'frz_mgrid.o']
if parallel:
  warpobjects = warpobjects + ['f3dslave.o','frzslave.o','topslave.o',
                               'w3dslave.o']

warpobjects = map(lambda p:os.path.join(builddir,p),warpobjects)

library_dirs = fcompiler.libdirs
libraries = fcompiler.libs
if parallel:
  library_dirs = fcompiler.libdirs + ['/usr/lpp/ppe.poe/lib']
  libraries = fcompiler.libs + ['mpi']
  #warpobjects = warpobjects + ['/usr/local/mpi/ifc_farg.o']

setup (name = "warpC",
       version = '3.0',
       author = 'David P. Grote',
       author_email = "DPGrote@lbl.gov",
       description = "Combines warp's packages into one",
       platforms = "Unix, Windows (cygwin), Mac OSX",
       ext_modules = [Extension('warpC',
                                ['warpC_Forthon.c',
                                 os.path.join(builddir,'Forthon.c'),
                                 'pmath_rng.c','ranf.c','ranffortran.c'],
                                include_dirs=[builddir],
                                library_dirs=library_dirs,
                                libraries=libraries,
                                extra_objects=warpobjects,
                                extra_link_args=['-g']+
                                             fcompiler.extra_link_args)]
       )
