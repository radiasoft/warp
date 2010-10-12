#!/usr/bin/env python
# To use:
#       python setup.py install
#

import os

setup_version = "$Id: setup.py,v 1.4 2010/10/12 23:31:54 dave Exp $"

try:
    # --- For installing into site-packages, with python setup.py install
    from distutils.core import setup
    # --- For creating an egg file, with python setup.py bdist_egg
    #from setuptools import setup
except:
    raise SystemExit, 'Distutils problem'

# --- Hack warning:
# --- Put warpoptions.py and parallel.py into subdirectories so that they
# --- can be separate packages. This is needed so that the two modules can be
# --- imported independently of (and before) warp.
os.renames('warpoptions.py','warpoptions/warpoptions.py')
ff = open('warpoptions/__init__.py','w')
ff.write('from warpoptions import *\n')
ff.close()

os.renames('parallel.py','parallel/parallel.py')
ff = open('parallel/__init__.py','w')
ff.write('from parallel import *\n')
ff.close()

setup (name = 'warp',
       author = 'David P. Grote, Jean-Luc Vay, et. al.',
       author_email = 'DPGrote@lbl.gov',
       description = 'Warp scripts',
       long_description = """
Warp scripts""",
       platforms = 'Linux, Unix, Windows (cygwin), Mac OSX',
       packages = ['warp','warp.GUI','warpoptions','parallel'],
       package_dir = {'warp': '.'},
       package_data = {'warp':['*.gs','*.gp']}
       )

# --- Undo the hack from above
os.remove('warpoptions/__init__.py')
os.renames('warpoptions/warpoptions.py','warpoptions.py')

os.remove('parallel/__init__.py')
os.renames('parallel/parallel.py','parallel.py')

