#!/usr/bin/env python
# To use:
#       python setup.py install
#

import os

setup_version = "$Id: setup.py,v 1.7 2011/10/27 22:07:22 grote Exp $"

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
os.renames('warpoptions.py','warpoptions/__init__.py')
os.renames('parallel.py','parallel/__init__.py')

try:
    setup (name = 'warp',
           version = '3.0',
           author = 'David P. Grote, Jean-Luc Vay, et. al.',
           author_email = 'dpgrote@lbl.gov',
           description = 'Warp PIC accelerator code',
           long_description = """
Warp is a PIC code designed to model particle accelerators and similar
machines that are space-charge dominated.""",
           url = 'http://warp.lbl.gov',
           platforms = 'Linux, Unix, Windows (cygwin), Mac OSX',
           packages = ['warp','warp.GUI','warpoptions','parallel'],
           package_dir = {'warp': '.'},
           package_data = {'warp':['*.gs','*.gp']}
           )

finally:
    # --- Undo the hack from above. Inside the finally clause, this will
    # --- always happen, even if there is an error during setup.
    os.renames('warpoptions/__init__.py','warpoptions.py')
    os.renames('parallel/__init__.py','parallel.py')

