# import all of the neccesary packages
from Numeric import *
from types import *
import RNG
import RandomArray
from pybasisC import *
import os
import string
import re
try:
  import PW
  import PR
except ImportError:
  pass
import __main__
import sys
import cPickle
Basis_version = "$Id: pyBasis.py,v 1.21 2002/05/15 00:19:54 dave Exp $"

if sys.platform in ['sn960510','linux-i386','linux2']:
  true = -1
  false = 0
else:
  true = 1
  false = 0

# --- Convenience function modeled after the iota of basis
def iota(low,high=None,step=1):
  if high is None:
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
    # --- Check if it is a WARP variable
    try:
      listvar(f)
      return
    except NameError:
      pass
    # --- Check if it is a module name
    try:
      exec(f+'doc()',__main__.__dict__)
      return
    except:
      pass
    # --- Try to get the actual value of the object
    try:
      f = eval(f,__main__.__dict__)
    except:
      print "Name not found"
      return
  # --- Check if it has a doc string
  try:
    print f.__doc__
    return
  except:
    pass
  print "No documentation found"

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
def pydump(fname=None,attr=["dump"],vars=[],serial=0,ff=None,varsuffix=None,
           verbose=false):
  """
Dump data into a pdb file
  - fname dump file name
  - attr attribute or list of attributes of variables to dump
         Any items that are not strings are skipped. To write no variables,
         use attr=None.
  - vars list of python variables to dump
  - serial switch between parallel and serial versions
  - ff=None Allows passing in of a file object so that pydump can be called
            multiple times to pass data into the same file. Note that
            the file must be explicitly closed by the user.
  - varsuffix=None Suffix to add to the variable names. If none is specified,
                the suffix '@pkg' is used, where pkg is the package name
                that the variable is in.
  - verbose=false When true, prints out the names of the variables as they are
                  written to the dump file
  """
  assert fname is not None or ff is not None,\
         "Either a filename must be specified or a pdb file pointer"
  # --- Open the file if the file object was not passed in.
  # --- If the file object was passed in, then don't close it.
  if not ff:
    ff = PW.PW(fname)
    closefile = 1
  else:
    closefile = 0
  # --- Convert attr into a list if needed
  if not (type(attr) == type([])): attr = [attr]
  # --- Loop through all of the packages in reverse order (getting pkg object)
  # --- Reverse order is used so that if a varsuffix is specified and when
  # --- a variable that have the same names in different packages is written,
  # --- the one in the package with higher precedence is used.
  pkgsuffix = varsuffix
  packagelist = package()
  packagelist.reverse()
  for pname in packagelist:
    pkg = eval(pname,__main__.__dict__)
    if varsuffix is None: pkgsuffix = '@' + pname
    # --- Get variables in this package which have attribute attr.
    vlist = []
    for a in attr:
      if type(a) == StringType: vlist = vlist + pkg.varlist(a)
    # --- Loop over list of variables
    for vname in vlist:
      # --- Check if object is available (i.e. check if dynamic array is
      # --- allocated).
      v = pkg.getpyobject(vname)
      if v is not None:
        writevar = 1
        # --- If serial flag is set, get attributes and if has the parallel
        # --- attribute, don't write it.
        if serial:
          a = pkg.getvarattr(vname)
          if re.search('parallel',a):
            writevar = 0
        # --- Check if variable is a complex array. Currently, these
        # --- can not be written out.
        if type(v) == type(array([0.])) and v.typecode() == Complex:
          writevar = 0
        # --- Write out the variable.
        if writevar:
          if verbose: print "writing "+pname+"."+vname+" as "+vname+pkgsuffix
          ff.write(vname+pkgsuffix,v)
  # --- Now, write out the python variables (that can be written out).
  # --- If supplied, the varsuffix is append to the names here too.
  if varsuffix is None: varsuffix = ''
  for v in vars:
    vval = eval(v,__main__.__dict__,locals())
    # --- Don't try to write out functions or classes. (They don't seem to
    # --- cause problems but this avoids potential problems. The function
    # --- or class body wouldn't be written out anyway.)
    if type(vval) in [FunctionType,ClassType]: continue
    # --- Zero length arrays cannot by written out.
    if type(vval) == type(array([])) and product(array(shape(vval))) == 0:
      continue
    # --- Try writing as normal variable.
    # --- The docontinue temporary is needed since python1.5.2 doesn't
    # --- seem to like continue statements inside of try statements.
    docontinue = 0
    try:
      if verbose: print "writing python variable "+v+" as "+v+varsuffix
      ff.write(v+varsuffix,vval)
      docontinue = 1
    except:
      pass
    if docontinue: continue
    # --- If that didn't work, try writing as a pickled object
    try:
      if verbose:
        print "writing python variable "+v+" as "+v+varsuffix+'@pickle'
      ff.write(v+varsuffix+'@pickle',cPickle.dumps(vval,0))
      docontinue = 1
    except:
      pass
    if docontinue: continue
    # --- All attempts failed so write warning message
    if verbose: print "cannot write python variable "+v
  if closefile: ff.close()

# --- Old version which has different naming for variables
def pydumpold(fname,attr="dump",vars=[]):
  ff = PW.PW(fname)
  for p in package():
    vlist = eval(p+'.varlist("'+attr+'")',__main__.__dict__)
    for v in vlist:
      if eval(p+'.getpyobject("'+v+'")',__main__.__dict__) is not None:
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
def pyrestore(filename=None,fname=None,verbose=0,skip=[]):
  """
Restores all of the variables in the specified file.
  - filename: file to read in from (assumes PDB format)
  - verbose=0: When true, prints out the names of variables which are read in
  - skip=[]: list of variables to skip
  """
  # --- The original had fname, but changed to filename to be consistent
  # --- with restart and dump.
  if filename is None: filename = fname
  # --- Make sure a filename was inputted.
  assert filename is not None,"A filename must be specified"
  # --- open pdb file
  ff = PR.PR(filename)
  # --- Get a list of all of the variables in the file, loop over that list
  vlist = ff.inquire_ls()
  # --- vlist is looped over twice. The first time reads in all of the scalar
  # --- variables and the python variables. The second reads in the arrays.
  # --- This is done so that the integers which specify the dimensions of
  # --- of arrays are read in before the array itself is. In some cases,
  # --- not having the integers read in first would cause problems
  # --- (so far only in the f90 version).
  for v in vlist:
    # --- If v in the skip list, then continue
    if v in skip:
      if verbose: print "skipping "+v
      continue
    # --- If the variable has the suffix '@pkg' then it is a warp variable.
    if len(v) > 4 and v[-4]=='@':
      # --- If v in the skip list, then continue
      if v[:-4] in skip or v[-3:]+'.'+v[:-4] in skip:
        if verbose: print "skipping "+v
        continue
      try:
        # --- get the package that the variable is in
        # --- Array assignment is different than scalar assignment.
        if type(ff.__getattr__(v)) != type(array([])):
          # --- Simple assignment is done for scalars, using the exec command
          if verbose: print "reading in "+v[-3:]+"."+v[:-4]
          #pkg.__setattr__(v[-3:],ff.__getattr__(v))
          exec(v[-3:]+'.'+v[:-4]+'=ff.__getattr__(v)',
               __main__.__dict__,locals())
      except:
        # --- The catches errors in cases where the variable is not an
        # --- actual warp variable, for example if it had been deleted
        # --- after the dump was originally made.
        print "Warning: There was problem restoring %s"% (v[-3:]+'.'+v[:-4])
    elif v[-7:] == '@pickle':
      # --- Thses would be interpreter variables written to the file
      # --- as pickled objects. The data is unpickled and the variable
      # --- in put in the main dictionary.
      try:
        if verbose: print "reading in python variable "+v[:-7]
        __main__.__dict__[v[:-7]] = cPickle.loads(ff.__getattr__(v))
      except:
        if verbose: print "error with variable "+v[:-7]
    elif v[-7:] == '@global':
      # --- These would be interpreter variables written to the file
      # --- from Basis. A simple assignment is done and the variable
      # --- in put in the main dictionary.
      try:
        if verbose: print "reading in python variable "+v[:-7]
        __main__.__dict__[v[:-7]] = ff.__getattr__(v)
      except:
        if verbose: print "error with variable "+v[:-7]
    elif v[-9:] == '@parallel':
      # --- Don't read in variables with suffix @parallel
      pass
    else:
      # --- These would be interpreter variables written to the file
      # --- from python (or other sources). A simple assignment is done and
      # --- the variable in put in the main dictionary.
      try:
        if verbose: print "reading in python variable "+v
        __main__.__dict__[v] = ff.__getattr__(v)
      except:
        if verbose: print "error with variable "+v[:-7]
  # --- Now loop again to read in the arrays.
  for v in vlist:
    # --- If v in the skip list, then continue
    if v in skip:
      if verbose: print "skipping "+v
      continue
    if len(v) > 4 and v[-4]=='@':
      # --- If v in the skip list, then continue
      if v[:-4] in skip or v[-3:]+'.'+v[:-4] in skip:
        if verbose: print "skipping "+v
        continue
      try:
        pkg = eval(v[-3:],__main__.__dict__)
        if type(ff.__getattr__(v)) == type(array([])):
          # --- forceassign is used, allowing the array read in to have a
          # --- different size than the current size of the warp array.
          if verbose: print "reading in "+v[-3:]+"."+v[:-4]
          pkg.forceassign(v[:-4],ff.__getattr__(v))
      except:
        print "Warning: There was problem restoring %s"% (v[-3:]+'.'+v[:-4])
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



def Basisdoc():
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
