#!/usr/bin/env python
# To use:
#       python setup.py install
#
from distutils.core import setup, Extension

#DX = "/usr/local/dx"
DX = "/usr/local/dx4.3.2/dx"
DXLIBPATH =   ["%s/lib_linux"%DX,"/usr/X11R6/lib","/usr/local/lib"]
DXLIBS = ["DXcallm","DXlite","GL","GLU","X11","Xm","Xt","Magick","dl",
          "netcdf","cdf"]

setup (name = "pyDXObject",
       version = '1',
       description = "Interface to OpenDX",
       extra_path = 'pyDXObject',
       packages = [''],
       include_dirs = ['%s/include'%DX],
       ext_modules = [Extension('pyDXObject',
                                ['pyDXObject.c'],
                                library_dirs=DXLIBPATH,
                                libraries=DXLIBS)
                     ]
       )

