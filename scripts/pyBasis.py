# import all of the neccesary packages
from Numeric import *
from pybasis import *
import os
import string
import re
import RandomArray
import PW
import PR
import __main__
import sys
pyBasis_version = "$Id: pyBasis.py,v 1.1.1.1 2000/10/16 18:34:19 dave Exp $"

if sys.platform in ['sn960510','linux-i386']:
  true = -1
  false = 0
else:
  true = 1
  false = 0

# --- Convenience function modeled after the iota of basis
def iota(low,high=None,step=1):
  if high==None:
    if step > 0:
      return arange(1,low+1,step)
    else:
      return arange(low,0,step)
  else:
    if step > 0:
      return arange(low,high+1,step)
    else:
      return arange(low,high-1,step)

# --- Converts an array of characters into a string.
def arraytostr(a):
  if len(shape(a)) == 1:
    result = ''
    for c in a:
      result = result + c
  elif len(shape(a)) == 2:
    result = []
    for i in xrange(shape(a)[1]):
      s = ''
      for c in a[:,i]:
        s = s + c
      result.append(s)
  return result

# --- Convenience function to do printing
def remark(s):
  print s

# --- Allows int operation on arrrays
builtinint = int
def int(x):
  if type(x) == type(array([])):
    return x.astype(Int)
  else:
    return builtinint(x)

# --- Return the nearest integer
def nint(x):
  if type(x) == type(array([])):
    return where(greater(x,0),int(x+0.5),-int(abs(x)+0.5))
  else:
    if x >= 0: return int(x+0.5)
    else: return -int(abs(x)+0.5)

# --- Replicate the sign function
def sign(x,y):
  if type(x) == type(array([])):
    result = where(greater(y,0.),abs(x),-abs(x))
    result = where(equal(y,0.),0.,result)
    return result
  else:
    if y > 0:
      return abs(x)
    elif y < 0:
      return -abs(x)
    else:
      return 0

# --- These are replacements for array creation routines which create
# --- arrays which have the proper ordering for fortran. When arrays created
# --- with these commands are passed to a fortran subroutine, no copies are
# --- needed to get the data into the proper order for fortran.
def fones(shape,typecode=Int):
  try:
    s = list(shape)
  except TypeError:
    s = list([shape])
  s.reverse()
  return transpose(ones(s,typecode))
def fzeros(shape,typecode=Int):
  try:
    s = list(shape)
  except TypeError:
    s = list([shape])
  s.reverse()
  return transpose(zeros(s,typecode))

# --- This function appends a new element to the end of an array.
# --- It is not very efficient since it creates a whole new array each time.
def arrayappend(x,a):
  xshape = list(shape(x))
  if type(a) == type(array([])):
    pass
  elif type(a) == type([]):
    a = array(a)
  else:
    a = array([a])
  ashape = list(shape(a))
  if len(xshape)==1 and len(ashape)==1:
    xshape[0] = xshape[0] + ashape[0]
    y = zeros(xshape,x.typecode())
    y[0:xshape[0]-ashape[0]] = x
    y[xshape[0]-ashape[0]:] = a
  elif len(xshape)>1 and len(ashape)==1 and xshape[-2]==ashape[0]:
    xshape[-1] = xshape[-1] + 1
    y = zeros(xshape,x.typecode())
    y[...,0:-1] = x
    y[...,-1] = a
  return y

# Convenience function which returns true if variable exists
def exists(x):
  try:
    xtemp = eval(x,locals(),globals())
    return true
  except NameError:
    return false

# Returns the average of the input array
def ave(x,index=0):
  if shape(x)[index] > 0:
    return sum(x,index)/shape(x)[index]
  else:
    return 0.

# --- Returns the max of the multiarray
def maxnd(x):
  """Return the max element of an array of any dimension"""
  xtemp = reshape(x,tuple([product(array(x.shape))]))
  return max(xtemp)
# --- Returns the min of the multiarray
def minnd(x):
  """Return the min element of an array of any dimension"""
  xtemp = reshape(x,tuple([product(array(x.shape))]))
  return min(xtemp)

# Gets next available filename with the format 'root.nnn.suffix'.
def getnextfilename(root,suffix):
  dir = string.join(os.listdir('.'))
  i = 0
  name = root+('.%03d.'%i)+suffix
  while re.search(name,dir):
    i = i + 1
    name = root+('.%03d.'%i)+suffix
  return name

# --- Prints out the documentation of the subroutine or variable.
def doc(f):
  if type(f) == type(''):
    try:
      listvar(f)
    except NameError:
      print eval(f+'.__doc__',__main__.__dict__)
  else:
    print f.__doc__

##############################################################################
# Python version of the dump routine. This uses the varlist command to
# list of all of the variables in each package which have the
# attribute attr (and actually attr could be a group name too). It then
# checks on the state of the python object, making sure that unallocated
# arrays are not written out.  Finally, the variable is written out to the
# file with the name in the format vame@pkg.  Additionally, python
# variables can be written to the file by passing in a list of the names
# through vars. The '@' sign is used between the package name and the
# variable name so that no python variable names can be clobbered ('@'
# is not a valid character in python names). The 'ff.write' command is
# used, allowing names with an '@' in them. The writing of python variables
# is put into a 'try' command since some variables cannot be written to
# a pdb file.
# Some fancy foot work is done in the exec and eval commands to get access
# to the global name space. The exec also needed to have the local name
# space explicitly included.
# Note that attr can be a list of attributes and group names.
def pydump(fname,attr=["dump"],vars=[],serial=0,ff=None):
  """
Dump data into a pdb file
  - fname dump file name
  - attr attribute or list of attributes of variables to dump
  - vars list of python variables to dump
  - serial switch between parallel and serial versions
  - ff=None Allows passing in of a file object so that pydump can be called
	    multiple times to pass data into the same file. Note that
	    the file must be explicitly closed by the user.
  """
  # --- Open the file if the file object was not passed in.
  # --- If the file object was passed in, then don't close it.
  if not ff:
    ff = PW.PW(fname)
    closefile = 1
  else:
    closefile = 0
  # --- Convert attr into a list if needed
  if not (type(attr) == type([])): attr = [attr]
  # --- Loop through all of the packages (getting pkg object)
  for pname in package():
    pkg = eval(pname,__main__.__dict__)
    # --- Get variables in this package which have attribute attr.
    vlist = []
    for a in attr: vlist = vlist + pkg.varlist(a)
    # --- Loop over list of variables
    for vname in vlist:
      # --- Check if object is available (i.e. check if dynamic array is
      # --- allocated).
      v = pkg.getpyobject(vname)
      if v!=[]:
        writevar = 1
        # --- If serial flag is set, get attributes and if has the parallel
        # --- attribute, don't write it.
        if serial:
          a = pkg.getvarattr(vname)
          if re.search('parallel',a):
            writevar = 0
        if writevar:
          ff.write(vname+"@"+pname,v)
  # --- Now, write out the python variables (that can be written out).
  for v in vars:
    try:
      exec('ff.'+v+'='+v,__main__.__dict__,locals())
    except:
      pass
  if closefile: ff.close()

# --- Old version which has different naming for variables
def pydumpold(fname,attr="dump",vars=[]):
  ff = PW.PW(fname)
  for p in package():
    vlist = eval(p+'.varlist("'+attr+'")',__main__.__dict__)
    for v in vlist:
      if eval(p+'.getpyobject("'+v+'")',__main__.__dict__)!=[]:
        #exec('ff.'+p+'_'+v+'='+p+'.'+v,__main__.__dict__,locals())
        exec('ff.write("'+p+'@'+v+'",'+p+'.'+v+')',__main__.__dict__,locals())
  for v in vars:
    try:
      exec('ff.'+v+'='+v,__main__.__dict__,locals())
    except:
      pass
  ff.close()


# Python version of the restore routine. It restores all of the variables
# in the pdb file, including any that are not part of a pybasis package.
# An '@' in the name distinguishes between the two. The 'ff.__getattr__' is
# used so that variables with an '@' in the name can be read. The reading
# in of python variables is put in a 'try' command to make it idiot proof.
# More fancy foot work is done to get new variables read in into the
# global dictionary.
def pyrestore(fname,verbose=0):
  """
Restores all of the variables in the specified file.
  - fname file to read in from (assumes PDB format)
  - verbose=0 When true, prints out the names of variables which are read in
  """
  # --- open pdb file
  ff = PR.PR(fname)
  # --- Get a list of all of the variables in the file, loop over that list
  vlist = ff.inquire_ls()
  for v in vlist:
    if verbose: print v
    # --- If the variable has the suffix '@pkg' then it is a warp variable.
    if len(v) > 4 and v[-4]=='@':
      try:
        # --- get the package that the variable is in
        pkg = eval(v[-3:],__main__.__dict__)
        # --- Array assignment is different than scalar assignment.
        if type(ff.__getattr__(v)) == type(array([])):
          # --- forceassign is used, allowing the array read in to have a
          # --- different size than the current size of the warp array.
          pkg.forceassign(v[:-4],ff.__getattr__(v))
        else:
          # --- Simple assignment is done for scalars, using the exec command
          exec(v[-3:]+'.'+v[:-4]+'=ff.__getattr__(v)',
               __main__.__dict__,locals())
      except:
        # --- The catches errors in cases where the variable is not an
        # --- actual warp variable, for example if it had been deleted
        # --- after the dump was originally made.
        print "Warning: There was problem restoring %s"% (v[-3:]+'.'+v[:-4])
    elif v[-7:] == '@global':
      # --- These would be interpreter variables written to the file
      # --- from Basis. A simple assignment is done and the variable
      # --- in put in the main dictionary.
      try:
        #exec('__main__.__dict__[v[:-7]]=ff.read(v)',locals())
        exec('%s=ff.__getattr__("%s");__main__.__dict__["%s"]=%s'%
             (v[:-7],v,v[:-7],v[:-7]))
      except:
        pass
    else:
      # --- These would be interpreter variables written to the file
      # --- from python (or other sources). A simple assignment is done and
      # --- the variable in put in the main dictionary.
      try:
        #exec('__main__.__dict__["%s"]=ff.%s'%(v,v)) # This should work
        exec('%s=ff.%s;__main__.__dict__["%s"]=%s'%(v,v,v,v))
      except:
        pass
  ff.close()

# --- create an alias for pyrestore
restore = pyrestore

def restoreold(fname):
  ff = PR.PR(fname)
  vlist = ff.inquire_ls()
  for v in vlist:
    if len(v) > 4 and v[3]=='@':
      if type(eval('ff.__getattr__("'+v+'")')) == type(array([])):
        exec(v[0:3]+'.forceassign("'+v[4:]+'",ff.__getattr__("'+v+'"))',
             __main__.__dict__,locals())
      else:
        #exec(v[0:3]+'.'+v[4:]+'=ff.'+v,__main__.__dict__,locals())
        exec(v[0:3]+'.'+v[4:]+'=ff.__getattr__("'+v+'")',
             __main__.__dict__,locals())
    else:
      try:
        exec('%s=ff.%s;__main__.__dict__["%s"]=%s'%(v,v,v,v))
      except:
        pass
  ff.close()



def pyBasisdoc():
  print """
iota(): returns a array of sequential integers
arraytostr(): converts an array of chars to a string
remark(): same as print
int(): converts data to integer
nint(): converts data to nearest integer
sign(): emulation of sign function
fones(): returns multi-dimensional array with fortran ordering
fzeros(): returns multi-dimensional array with fortran ordering
arrayappend(): appends to multi-dimensional array
exists(): checks if a variable exists
ave(): averages an array of numbers
maxnd(): finds max of multi-dimensional array
minnd(): finds min of multi-dimensional array
getnextfilename(): finds next available file name in a numeric sequence
doc(): prints info about variables and functions
pydump(): dumps data into pdb format file
pyrestore(): reads data from pdb format file
restore(): equivalent to pyrestore
"""
