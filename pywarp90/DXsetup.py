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
       version = '1.0',
       author = 'David P. Grote',
       author_email = "DPGrote@lbl.gov",
       description = "Interface to OpenDX",
       platforms = "Unix, Windows (cygwin), Mac OSX",
       ext_modules = [Extension('pyDXObject',
                                ['pyDXObject.c'],
                                include_dirs = ['%s/include'%DX],
                                library_dirs=DXLIBPATH,
                                libraries=DXLIBS)
                     ]
       )

