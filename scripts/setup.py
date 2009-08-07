#!/usr/bin/env python
# To use:
#       python setup.py install
#

try:
    # --- For installing into site-packages, with python setup.py install
    #from distutils.core import setup
    # --- For creating an egg file, with python setup.py bdist_egg
    from setuptools import setup
except:
    raise SystemExit, 'Distutils problem'

setup (name = 'warp',
       author = 'David P. Grote, Jean-Luc Vay, et. al.',
       author_email = 'DPGrote@lbl.gov',
       description = 'Warp scripts',
       long_description = """
Warp scripts""",
       platforms = 'Linux, Unix, Windows (cygwin), Mac OSX',
       packages = ['warp','warp.GUI'], # Note that the GUI info is incomplete
       package_dir = {'warp': '.'},
       package_data = {'warp':['*.gs','*.gp']}
       )
