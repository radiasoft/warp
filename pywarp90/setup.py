#!/usr/bin/env python
# To use:
#       python setup.py install
#
import os,os.path
from Forthon.compilers import FCompiler

try:
    import distutils
    from distutils.core import setup, Extension
    from distutils.dist import Distribution
    from distutils.command.build import build
except:
    raise SystemExit, "Distutils problem"

lapacklibdir = ['/home/dave/lapack_ifc/LAPACK']
lapacklibs = ['lapackifc8','blasifc8']

fcomp = FCompiler()

dummydist = Distribution()
dummybuild = build(dummydist)
dummybuild.finalize_options()
builddir = dummybuild.build_temp

warppkgs = ['top','env','w3d','f3d','wxy','fxy','wrz','frz','her','cir','cho']

def makeobjects(pkg):
  return [pkg+'.o',pkg+'_p.o',pkg+'pymodule.o']

warpobjects = []
for pkg in warppkgs:
  warpobjects = warpobjects + makeobjects(pkg)

warpobjects = warpobjects + ['dtop.o',
                             'dw3d.o',
                             'f3d_mgrid.o','f3d_conductors.o','fft.o','util.o',
                             'fxy_mgrid.o',
                             'dwrz.o',
                             'frz_mgrid.o']

warpobjects = map(lambda p:os.path.join(builddir,p),warpobjects)

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
                                library_dirs=fcomp.libdirs+lapacklibdir,
                                libraries=fcomp.libs+lapacklibs,
                                extra_objects=warpobjects,
                                extra_link_args=['-g'])]
       )
