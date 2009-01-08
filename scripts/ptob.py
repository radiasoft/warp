from numpy import *
ptob_version = "$Id: ptob.py,v 1.2 2009/01/08 19:27:35 dave Exp $"


def writetobasis(f,xname,x):
  """Writes a variable to the pdb file f with the appropriate mangling
needed to be readable by Basis
  - f is the pdb file which has been opened for writing
  - xname is the string name of the variable to be written
  - x is the variable itself
  """
  if type(x) == type(array([])):
    xtemp = transpose(x) + 0.
    s = list(shape(xtemp))
    s.reverse()
    xtemp.shape = s
    f.write(xname,xtemp)
  else:
    f.write(xname,x)
