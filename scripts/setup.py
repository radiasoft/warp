#!/usr/bin/env python
# To use:
#       python setup.py install
#

import os
import sys
import glob

try:
    import distutils
    from distutils.command.install import INSTALL_SCHEMES
    # --- For installing into site-packages, with python setup.py install
    from distutils.core import setup
    # --- For creating an egg file, with python setup.py bdist_egg
    #from setuptools import setup
except:
    raise SystemExit('Distutils problem')

try:
    from distutils.command.build_py import build_py_2to3 as build_py
except ImportError:
    from distutils.command.build_py import build_py

# --- Get around a "bug" in disutils on 64 bit systems. When there is no
# --- extension to be installed, distutils will put the scripts in
# --- /usr/lib/... instead of /usr/lib64. This fix will force the scripts
# --- to be installed in the same place as the .so (if an install was
# --- done in pywarp90 - though that shouldn't be done anymore).
if distutils.sysconfig.get_config_vars()["LIBDEST"].find('lib64') != -1:
    for scheme in INSTALL_SCHEMES.values():
        scheme['purelib'] = scheme['platlib']

# --- With this, the data_files listed in setup, the .so files from pywarp90,
# --- will be installed in the usual place in site-packages along with the scripts.
for scheme in INSTALL_SCHEMES.values():
    scheme['data'] = scheme['platlib']

# --- Get the list of .so files. These should have been copied into scripts
# --- at the end of the make.
data_files = glob.glob('*.so')

# --- Hack warning:
# --- Put warpoptions.py and parallel.py into subdirectories so that they
# --- can be separate packages. This is needed so that the two modules can be
# --- imported independently of (and before) warp.
os.renames('warpoptions.py', 'warpoptions/__init__.py')
os.renames('parallel.py', 'parallel/__init__.py')

try:
    setup(name='warp',
          version='3.0',
          author='David P. Grote, Jean-Luc Vay, et. al.',
          author_email='dpgrote@lbl.gov',
          description='Warp PIC accelerator code',
          long_description="""
Warp is a PIC code designed to model particle accelerators and similar
machines that are space-charge dominated.""",
          url='http://warp.lbl.gov',
          platforms='Linux, Unix, Windows (cygwin), Mac OSX',
          packages=['warp', 'warp.GUI', 'warpoptions', 'parallel'],
          package_dir={'warp': '.'},
          package_data={'warp': ['*.gs', '*.gp', 'aladdin_8.txt']},
          data_files=[('warp', data_files)],
          cmdclass={'build_py': build_py}
          )

finally:
    # --- Undo the hack from above. Inside the finally clause, this will
    # --- always happen, even if there is an error during setup.
    os.renames('warpoptions/__init__.py', 'warpoptions.py')
    os.renames('parallel/__init__.py', 'parallel.py')
