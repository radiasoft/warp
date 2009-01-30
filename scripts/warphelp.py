"""Contains the command warphelp() which prints a list of handy commands.
"""
from warp import *
warphelp_version = "$Id: warphelp.py,v 1.2 2009/01/30 23:45:50 dave Exp $"

def warphelpdoc():
  import warphelp
  print warphelp.__doc__

##########################################################################
def warphelp():
  import warpplots

  print """
WARP help card
For addition documentation on any command or variable, type doc(name) or
doc('name') for fortran variables.
"""

  print warpplots.__doc__

  print """
dump(): creates a dump file
restart(): restart simulation from a dump file
loadrho(): loads charge density array
fieldsol(): perform a Poisson solve
installbeforefs(): installs a function to be called before a field-solve
installafterfs(): installs a function to be called after a field-solve
installbeforestep(): installs a function to be called before a step
installafterstep(): installs a function to be called after a step
"""

  print """
For a description of available scripts, type
warpscripts()

More information about what is in the scripts can be obtained by importing the
script and typing scriptnamedoc(). For example...

from histplot import *
histplotdoc()

For descriptions of useful fortran routines type
warpfortran()
"""
##########################################################################
